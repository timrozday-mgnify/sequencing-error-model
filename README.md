# sequencing-error-model

> Working name. Early development: the native generator and the `pe-overlap` mode work; the `reference`
> and `kmer` modes and the exporters are planned.

Train sequencing error models from read evidence. Models are
quality-aware: they estimate the error profile of a base from its reported quality and the
surrounding bases and qualities (quality is a feature, never treated as the true error rate),
and generated reads carry both bases and qualities. The trained models
can then be used to simulate reads, either with the built-in generator or by
exporting them to widely used read simulators (ART/art_modern, InSilicoSeq, Badread,
PBSIM3, NEAT, ...).

Three evidence modes train the same kind of model, so they can validate each other:

- **pe-overlap**: disagreements between R1 and R2 where paired reads overlap (Illumina;
  no reference).
- **reference**: alignment to a reference genome, spike-in, or a self-assembly of the sample.
- **kmer**: [skiver](https://github.com/GZHoffie/skiver) (k,v)-mer consensus (no reference),
  either from an unmodified, pinned skiver release (`skiver analyze`, the **default** build)
  or from a modified skiver (`skiver dump`, the **enhanced** build).

The pe-overlap and reference modes are built first and used to check the kmer mode.
Further sources (ONT duplex reads, GATK/DADA2 error tables) come later.

See [docs/IMPLEMENTATION_PLAN.md](docs/IMPLEMENTATION_PLAN.md) for scope, feasibility,
limitations and the phased plan.

## Usage

Fit a model from paired Illumina reads whose mates overlap (amplicon or short-insert libraries). See
[docs/pe_overlap.md](docs/pe_overlap.md) for what it identifies and its limits:

```bash
uv run python -m sequencing_error_model.sources.pe_overlap reads_R1.fastq.gz reads_R2.fastq.gz --output spec/
```

Simulate reads (bases, qualities and CIGAR) from a spec, as single templates or as standalone pairs:

```bash
uv run sem-generate --model spec/ --input genome.fasta --output reads.fastq --pairs 10000 --read-length 125 --insert-mean 300 --insert-sd 50
```

## Development

```bash
uv sync
task lint
task test
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the PR policy.
