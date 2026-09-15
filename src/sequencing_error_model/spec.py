"""ErrorModelSpec: the simulator-neutral model the generator and exporters read (plan §7).

A spec is a directory holding `spec.json` (schema version, quality alphabet, provenance,
component tokens, flags and meta) and `arrays.npz` (every parameter and cached marginal).
Loading never unpickles.

Head Q models P(q_t | ...) and head E models P(op_t | ...) (§5.2); components are named by
spec tokens (§5.4), e.g. `Context(2,2)`. `identified=False` marks parameters the data could
not identify and that were fixed instead (e.g. neighbouring-Q weights in default mode, §5.6);
`generative=False` marks training-only components (e.g. `FragmentOverdispersion`).
"""

import json
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

SCHEMA_VERSION = 1
MODES = ("default", "enhanced")
MAX_Q = 93
# Component name -> heads it may appear in (E: error, Q: quality). The shared `Latent(S)`
# layer sits outside both heads.
# ponytail: names and heads only; argument arity is checked by each component's fitter.
COMPONENTS = {
    **dict.fromkeys(("Context", "ContextTable", "Position", "Mate", "Strand", "Homopolymer"), "EQ"),
    **dict.fromkeys(
        ("QualityWindow", "QualityxContext", "GC", "Weibull", "IndelLength", "FragmentOverdispersion"), "E"
    ),
    **dict.fromkeys(("QualityMarkov", "InsertionQuality"), "Q"),
}
_TOKEN = re.compile(r"(\w+)(?:\((\d+(?:,\d+)*)\))?")

Array = npt.NDArray[Any]


class SpecError(ValueError):
    """A spec is invalid, malformed on disk, or from an unsupported schema version."""


@dataclass(frozen=True)
class Component:
    token: str
    params: dict[str, Array] = field(default_factory=dict)
    generative: bool = True
    identified: bool = True
    meta: dict[str, Any] = field(default_factory=dict)  # JSON values only

    @property
    def name(self) -> str:
        return self.token.partition("(")[0]

    @property
    def args(self) -> tuple[int, ...]:
        inner = self.token.partition("(")[2].rstrip(")")
        return tuple(int(a) for a in inner.split(",")) if inner else ()


@dataclass(frozen=True)
class ErrorModelSpec:
    quality_alphabet: tuple[int, ...]  # the only Q values generation may emit
    provenance: dict[str, Any]  # needs "mode"; sources, skiver version, k, v, c, input hashes, base profile
    quality_head: tuple[Component, ...]
    error_head: tuple[Component, ...]
    latent: Component | None = None
    marginals: dict[str, Array] = field(default_factory=dict)  # cached for exporters and reports

    def __post_init__(self) -> None:
        a = self.quality_alphabet
        if not a or list(a) != sorted(set(a)) or not all(type(q) is int and 0 <= q <= MAX_Q for q in a):
            raise SpecError(f"quality alphabet must be increasing integers in [0, {MAX_Q}], got {a}")
        if self.provenance.get("mode") not in MODES:
            raise SpecError(f"provenance mode must be one of {MODES}, got {self.provenance.get('mode')!r}")
        for head, components in (("Q", self.quality_head), ("E", self.error_head)):
            names = [c.name for c in components]
            for c in components:
                m = _TOKEN.fullmatch(c.token)
                if not m or head not in COMPONENTS.get(m[1], ""):
                    raise SpecError(f"{c.token!r} is not a head {head} component")
            if len(set(names)) != len(names):
                raise SpecError(f"head {head} repeats a component: {names}")
        if self.latent is not None and not re.fullmatch(r"Latent\(\d+\)", self.latent.token):
            raise SpecError(f"latent layer must be Latent(S), got {self.latent.token!r}")
        for key, arr in self._arrays().items():
            if not isinstance(arr, np.ndarray) or arr.dtype.hasobject:
                raise SpecError(f"{key}: parameters must be non-object numpy arrays")

    def _components(self) -> Iterator[tuple[str, Component]]:
        yield from ((f"quality/{i}", c) for i, c in enumerate(self.quality_head))
        yield from ((f"error/{i}", c) for i, c in enumerate(self.error_head))
        if self.latent is not None:
            yield "latent", self.latent

    def _arrays(self) -> dict[str, Array]:
        out = {f"marginals/{k}": v for k, v in self.marginals.items()}
        for prefix, c in self._components():
            out |= {f"{prefix}/{k}": v for k, v in c.params.items()}
        return out

    def save(self, path: Path) -> None:
        def entry(c: Component) -> dict[str, Any]:
            return {
                "token": c.token,
                "generative": c.generative,
                "identified": c.identified,
                "meta": c.meta,
                "params": list(c.params),
            }

        doc = {
            "schema_version": SCHEMA_VERSION,
            "quality_alphabet": list(self.quality_alphabet),
            "provenance": self.provenance,
            "quality_head": [entry(c) for c in self.quality_head],
            "error_head": [entry(c) for c in self.error_head],
            "latent": None if self.latent is None else entry(self.latent),
            "marginals": list(self.marginals),
        }
        path.mkdir(parents=True, exist_ok=True)
        (path / "spec.json").write_text(json.dumps(doc, indent=2) + "\n")
        np.savez_compressed(path / "arrays.npz", **self._arrays())  # type: ignore[arg-type]


def load(path: Path) -> ErrorModelSpec:
    try:
        doc = json.loads((path / "spec.json").read_text())
        if doc.get("schema_version") != SCHEMA_VERSION:
            raise SpecError(f"{path}: schema version {doc.get('schema_version')!r}, supported: {SCHEMA_VERSION}")
        with np.load(path / "arrays.npz", allow_pickle=False) as npz:
            arrays = {k: npz[k] for k in npz.files}

        def component(prefix: str, d: dict[str, Any]) -> Component:
            params = {k: arrays.pop(f"{prefix}/{k}") for k in d["params"]}
            return Component(d["token"], params, d["generative"], d["identified"], d["meta"])

        spec = ErrorModelSpec(
            tuple(doc["quality_alphabet"]),
            doc["provenance"],
            tuple(component(f"quality/{i}", d) for i, d in enumerate(doc["quality_head"])),
            tuple(component(f"error/{i}", d) for i, d in enumerate(doc["error_head"])),
            None if doc["latent"] is None else component("latent", doc["latent"]),
            {k: arrays.pop(f"marginals/{k}") for k in doc["marginals"]},
        )
    except (KeyError, TypeError) as e:
        raise SpecError(f"{path}: malformed spec: {e!r}") from e
    if arrays:
        raise SpecError(f"{path}: arrays not referenced by spec.json: {sorted(arrays)}")
    return spec
