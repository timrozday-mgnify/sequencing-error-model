#!/usr/bin/env nextflow
// Long runs of plan phase 6b, on Slurm (README.md beside this file).
//   real:  SRR24523812's three-mode comparison rerun at k >= 17 (item 1), plus the Latent(S) reference specs
//          the synthetic sweep uses as truths (item 2).
//   sweep: one task per row of --sweep (id,module,args,truth): `python -m sequencing_error_model.<module>
//          --skiver ... <args> [--spec <truth>] --output <id>.json`. truth names a reference spec from the
//          real arm (e.g. latent2), or is empty for the built-in example spec.
nextflow.enable.dsl = 2

// Every Nth pair, never a contiguous block: at equal read counts a head slice of this run reads up to
// 2.5x the error rate of a tail slice (plan phase 6b item 2), so a block would make every cross-mode ratio
// below a measurement of flowcell position. The same stride on both mates keeps them aligned, and all three
// modes are fed from this one output, so they see the same reads.
process SUBSET {
    tag "$name"
    label 'small'
    input:
    tuple val(name), val(every), path(r1), path(r2)
    output:
    tuple val(name), path("${name}_1.fastq.gz"), path("${name}_2.fastq.gz")
    script:
    """
    python -m sequencing_error_model.sources.fastq_quality $r1 --every-n $every --output s1
    python -m sequencing_error_model.sources.fastq_quality $r2 --every-n $every --output s2
    mv s1/* ${name}_1.fastq.gz
    mv s2/* ${name}_2.fastq.gz
    """
}

process ALIGN {
    tag "$name"
    label 'medium'
    input:
    tuple val(name), path(r1), path(r2)
    path assembly
    output:
    tuple val(name), path("${name}.sam")
    script:
    """
    minibwa index $assembly
    minibwa map -t $task.cpus $assembly $r1 $r2 > ${name}.sam
    """
}

process FIT_REFERENCE {
    tag "latent$latent"
    label 'fit'
    publishDir "${params.outdir}/specs", mode: 'copy'
    input:
    val latent
    tuple val(name), path(sam)
    path assembly
    output:
    tuple val("latent$latent"), path("reference_latent$latent")
    script:
    // Latent(S) with the phase 5 tokens (Position(8)); S = 0 is the phase 6 spec (the CLI's defaults).
    def tokens = latent ? "--latent $latent --error-tokens 'QualityWindow(1)' 'Context(1,1)' 'Position(8)' Homopolymer Mate" : ''
    """
    python -m sequencing_error_model.sources.bam $sam $assembly --unclip 2 8 12 2 \\
        --max-reads ${params.train_pairs * 2} --every-n ${params.train_every_n} \\
        $tokens --output reference_latent$latent
    """
}

process FIT_OVERLAP {
    label 'fit'
    publishDir "${params.outdir}/specs", mode: 'copy'
    input:
    tuple val(name), path(r1), path(r2)
    output:
    path 'pe_overlap'
    script:
    """
    python -m sequencing_error_model.sources.pe_overlap $r1 $r2 --max-pairs ${params.train_pairs} \\
        --every-n ${params.train_every_n} --output pe_overlap
    """
}

process KMER {
    tag "${name} k$k"
    label 'medium'
    publishDir "${params.outdir}/kmer", mode: 'copy'
    input:
    tuple val(name), path(r1), path(r2), val(k)
    output:
    tuple val(name), val(k), path("${name}_k$k")
    script:
    // skiver takes one file: both mates pooled, as in phase 6.
    def out = "${name}_k$k"
    """
    mkdir $out
    cat $r1 $r2 > reads.fastq.gz
    ${params.skiver} analyze reads.fastq.gz -k $k -v ${params.v} -c ${params.c} -t $task.cpus -o $out/analyze
    rm reads.fastq.gz
    python -m sequencing_error_model.sources.kmer $out/analyze $r1 $r2 --output $out/kmer_spec
    """
}

process COMPARE {
    tag "${kname} k$k"
    label 'fit'
    publishDir "${params.outdir}/compare", mode: 'copy'
    input:
    tuple val(kname), val(k), path(kmer)
    tuple val(name), path(r1), path(r2), path(sam)
    path assembly
    tuple val(rname), path(reference_spec)
    path overlap_spec
    output:
    path "${kname}_k${k}.json"
    script:
    """
    python -m sequencing_error_model.compare $r1 $r2 $sam $assembly --overlap-spec $overlap_spec \\
        --reference-spec $reference_spec --skiver $kmer/analyze --kmer-spec $kmer/kmer_spec \\
        --max-pairs ${params.train_pairs} --max-reads ${params.train_pairs * 2} \\
        --every-n ${params.train_every_n} --output ${kname}_k${k}.json
    """
}

process SWEEP {
    tag "$id"
    label 'sweep'
    publishDir "${params.outdir}/sweep", mode: 'copy'
    input:
    tuple val(id), val(module), val(args), path(truth)
    output:
    path "${id}.json"
    script:
    def spec = truth.name == 'NO_TRUTH' ? '' : "--spec $truth"
    // Exit 1 is a report whose tolerances failed: that is a result, so keep it.
    """
    python -m sequencing_error_model.$module --skiver ${params.skiver} --skiver-arg=-t --skiver-arg=$task.cpus \\
        $spec $args --output ${id}.json || { rc=\$?; [ \$rc -eq 1 ] && [ -s ${id}.json ] || exit \$rc; }
    """
}

workflow {
    def ks = params.k.toString().tokenize(',')*.toInteger()
    def latents = params.latent.toString().tokenize(',')*.toInteger()
    def rows = params.sweep ? channel.fromPath(params.sweep).splitCsv(header: true) : channel.empty()
    def specs = channel.empty()

    if (params.real) {
        def assembly = file(params.assembly, checkIfExists: true)
        def slices = [['subset', params.every_n]]
        if (params.depth_every_n) slices << ['depth', params.depth_every_n]
        def reads = SUBSET(channel.fromList(slices).map { it + [file(params.r1), file(params.r2)] })
        def subset = reads.filter { it[0] == 'subset' }
        def aligned = ALIGN(subset, assembly)
        specs = FIT_REFERENCE(channel.fromList(latents), aligned.first(), assembly)
        def overlap = FIT_OVERLAP(subset)
        def kmer = KMER(reads.combine(channel.fromList(ks)))
        // .first(): one-item channels become values, so every k reuses them.
        COMPARE(kmer, subset.join(aligned).first(), assembly, specs.filter { it[0] == 'latent0' }.first(), overlap.first())
    }

    // Rows without a truth run at once; rows naming one wait for that spec.
    def none = file("${projectDir}/NO_TRUTH")
    def bare = rows.filter { !it.truth }.map { [it.id, it.module, it.args, none] }
    def truthed = rows.filter { it.truth }.map { [it.truth, it] }
        .combine(specs, by: 0).map { _t, r, spec -> [r.id, r.module, r.args, spec] }
    SWEEP(bare.mix(truthed))
}
