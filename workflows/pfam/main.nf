#!/usr/bin/env nextflow
// P2 of docs/pfam_mode_plan.md: the pfam mode against libraries whose genomes are known (README.md beside this file).
// Per library: strided pairs at each dose -> six-frame ORFs -> kfpfam + HMMER -> pfam fits; the undosed pairs are also
// aligned to the genomes, and every fit is scored against them. Last: the go / no-go over all libraries.
nextflow.enable.dsl = 2

// Every Nth pair over the whole run, never a head slice (observations.subset). Named by number, so the alignment
// and the search speak of the same reads. Standard additions go in here, before any step can select on them.
process READS {
    tag "$run $dose"
    label 'small'
    input:
    tuple val(run), val(every), path(r1), path(r2), val(dose)
    output:
    tuple val(run), val(dose), path("${run}.kept_1.fastq.gz"), path("${run}.kept_2.fastq.gz")
    script:
    // .kept: an output named as the staged input would be written through its link, onto the source
    def limit = params.max_pairs ? "--max-pairs ${params.max_pairs}" : ''
    """
    python -m sequencing_error_model.sources.pfam reads $r1 $r2 --every-n $every $limit --dose $dose --output ${run}.kept
    """
}

process ORFS {
    tag "$run $dose"
    label 'small'
    input:
    tuple val(run), val(dose), path(r1), path(r2)
    output:
    tuple val(run), val(dose), path('orfs.faa')
    script:
    """
    python -m sequencing_error_model.sources.pfam orfs $r1 $r2 --output orfs.faa
    """
}

// kfpfam's Mode A, HMMER at E <= 1000: hits far below any anchor, which the exposure needs. Never a GA-cut table.
// The candidates are several times the ORFs' size, so they stay inside the task.
process SEARCH {
    tag "$run $dose"
    label 'search'
    input:
    tuple val(run), val(dose), path(orfs)
    path prep
    path index
    path hmm
    path pairs_py
    output:
    tuple val(run), val(dose), path('pairs.tsv')
    script:
    """
    export RAYON_NUM_THREADS=$task.cpus FILTER_MIN_GAPPED=${params.gapped}
    ${params.filter_bin} join $prep $index ${params.bar} ${params.gate} 400000 $orfs candidates.tsv
    ${params.pyhmmer_python} $pairs_py $hmm $prep/models.tsv $orfs candidates.tsv ${params.gapped} $task.cpus pairs.tsv
    rm candidates.tsv
    """
}

process FIT {
    tag "$run $name"
    label 'fit'
    publishDir { "${params.outdir}/fits/$run" }, mode: 'copy'
    input:
    tuple val(run), val(dose), path(r1), path(r2), path(pairs), val(name), val(anchors)
    path prep
    output:
    tuple val(run), path(name), emit: spec
    path "${name}.json"
    script:
    """
    python -m sequencing_error_model.sources.pfam fit $r1 $r2 --pairs $pairs --prep $prep --anchors $anchors \\
        --quality-every ${params.quality_every} --t1 ${params.t1} --t2 ${params.t2} --output $name > ${name}.json
    """
}

process ALIGN {
    tag "$run"
    label 'align'
    input:
    tuple val(run), val(dose), path(r1), path(r2), path(genomes)
    output:
    tuple val(run), path('aligned.sam')
    script:
    """
    minimap2 -ax sr -t $task.cpus --secondary=no $genomes $r1 $r2 > aligned.sam
    """
}

process SCORE {
    tag "$run"
    label 'score'
    publishDir "${params.outdir}/score", mode: 'copy'
    input:
    tuple val(run), path(r1), path(r2), path(pairs), path(sam), path(genomes), path(specs)
    path prep
    output:
    path "${run}.json"
    script:
    def named = specs.collect { "${it.name}=${it}" }.join(' ')
    """
    python -m sequencing_error_model.pfam_truth score $r1 $r2 --pairs $pairs --prep $prep --sam $sam \\
        --reference $genomes --t1 ${params.t1} --t2 ${params.t2} --spec $named --output ${run}.json
    """
}

// Reported beside the verdict, never judged (plan section 12): the truth, where there is a SAM, is for the undosed reads.
process LOCUS {
    tag "$run $dose"
    label 'locus'
    publishDir { "${params.outdir}/locus" }, mode: 'copy'
    input:
    tuple val(run), val(dose), path(r1), path(r2), path(pairs), path(truth)
    path prep
    output:
    path "${run}.${dose}.json"
    script:
    def known = truth ? "--sam ${truth[0]} --reference ${truth[1]}" : ''
    """
    python -m sequencing_error_model.pfam_locus $r1 $r2 --pairs $pairs --prep $prep $known \\
        --t1 ${params.t1} --t2 ${params.t2} --output ${run}.${dose}.json
    """
}

process VERDICT {
    label 'small'
    publishDir params.outdir, mode: 'copy'
    input:
    path reports
    path loci
    output:
    path 'verdict.*'
    script:
    """
    python -m sequencing_error_model.pfam_truth verdict $reports --loci $loci --output verdict
    """
}

def fastq(url) {
    return params.fastq_dir ? file("${params.fastq_dir}/${url.tokenize('/').last()}", checkIfExists: true) : file(url)
}

workflow {
    def need = ['prep', 'index', 'hmm', 'pairs_py'].findAll { !params[it] }
    if (need) error "set params.${need.join(', params.')} (site.config)"
    def prep = file(params.prep, type: 'dir', checkIfExists: true)
    def index = file(params.index, checkIfExists: true)
    def hmm = file(params.hmm, checkIfExists: true)
    def pairs_py = file(params.pairs_py, checkIfExists: true)
    def keep = params.only ? params.only.toString().tokenize(',') : null
    def rows = channel.fromPath(params.runs, checkIfExists: true).splitCsv(header: true)
        .filter { keep == null || it.run in keep }
    def doses = params.doses.toString().tokenize(',')
    // The spec's name says what it is: `both`, `stops` and `overlap` undosed, `doseX` the standard additions.
    def fits = [['0', 'both', 'both'], ['0', 'stops', 'stops'], ['0', 'overlap', 'overlap']] +
        doses.collect { [it, "dose$it".toString(), 'both'] }

    def reads = READS(rows.combine(channel.fromList(['0'] + doses))
        .map { r, dose -> [r.run, r.every_n, fastq(r.r1), fastq(r.r2), dose] })
    def pairs = SEARCH(ORFS(reads), prep, index, hmm, pairs_py)
    def ready = reads.join(pairs, by: [0, 1])
    def specs = FIT(ready.combine(channel.fromList(fits)).filter { it[1] == it[5] }
        .map { run, dose, r1, r2, p, _d, name, anchors -> [run, dose, r1, r2, p, name, anchors] }, prep).spec

    def genomes = rows.map { r ->
        if (!params.genomes[r.genomes]) error "${r.run}: no params.genomes.${r.genomes} (site.config)"
        [r.run, file(params.genomes[r.genomes], checkIfExists: true)]
    }
    def undosed = ready.filter { it[1] == '0' }
    def sam = ALIGN(reads.filter { it[1] == '0' }.combine(genomes, by: 0))
    SCORE(undosed.map { run, _d, r1, r2, p -> [run, r1, r2, p] }
        .join(sam).join(genomes).join(specs.groupTuple(size: fits.size(), remainder: true)), prep)
    def with_truth = undosed.combine(sam.join(genomes), by: 0)
        .map { run, dose, r1, r2, p, s, g -> [run, dose, r1, r2, p, [s, g]] }
    def without = ready.filter { it[1] != '0' }.map { it + [[]] }
    LOCUS(with_truth.mix(without), prep)
    VERDICT(SCORE.out.collect(), LOCUS.out.collect())
}
