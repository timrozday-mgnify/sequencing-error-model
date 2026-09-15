# sequencing-error-model

> Working name. Early planning stage: there is no usable functionality yet.

Train sequencing error models from [skiver](https://github.com/GZHoffie/skiver)
(k,v)-mer outputs, which need no reference genome or alignment. Optional extra evidence
sources — raw FASTQ quality profiles, paired-end overlaps, self-assembled references,
ONT duplex reads, GATK/DADA2 error tables — fill the gaps skiver can't cover. The trained models
can then be used to simulate reads, either with the built-in generator or by
exporting them to widely used read simulators (ART/art_modern, InSilicoSeq, Badread,
PBSIM3, NEAT, ...).

Two input modes:

- **default**: the CSV outputs of an unmodified, pinned skiver release (`skiver analyze`).
- **enhanced**: per-observation outputs from a modified skiver (`skiver dump`), which
  enable richer model components such as Phred, read-position, strand, homopolymer and
  read-level latent-state effects.

See [docs/IMPLEMENTATION_PLAN.md](docs/IMPLEMENTATION_PLAN.md) for scope, feasibility,
limitations and the phased plan.

## Development

```bash
uv sync
task lint
task test
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the PR policy.
