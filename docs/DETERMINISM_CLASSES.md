# Determinism classes (bead huq.5)

Every crate's `CONTRACT.md` declares one determinism class in its
"Determinism class" section (lint-enforced by `xtask check-contracts`).
The class is a CLAIM SURFACE: it decides which code rules bind
(`docs/CONVENTIONS.md`, "Determinism tiers and the libm doctrine"), which
Gauntlet tier verifies it, and how it interacts with `fs_exec::ExecMode`.
This document is the taxonomy; CONVENTIONS carries the code rules.

## The classes

| Class | Guarantee | Verified by | ExecMode interaction |
| --- | --- | --- | --- |
| `Deterministic` | Bit-identical results across runs, worker counts, build modes, AND ISAs (the four-quadrant bar: aarch64/x86-64 × debug/release) | G5 determinism audits: replay + golden hashes reproduced on both ISAs; goldens registered in `golden-couplings.json` per `docs/GOLDEN_POLICY.md` | Valid under `ExecMode::Deterministic` only; a `Fast` event stream cannot discharge a `Deterministic` claim |
| `DeterministicPerIsa` | Bit-identical across runs, worker counts, and build modes ON ONE ISA; last-ULP drift across ISAs/libm versions is admitted and scoped | G5 audits on a single host class: same-ISA replay + worker-count invariance; cross-ISA comparison is explicitly NOT evidence for or against | Valid under `ExecMode::Deterministic`; the CONTRACT must say "same-ISA" — an unqualified determinism claim at this class is a doctrine violation |
| `Fast` | Statistical/tolerance envelopes only; no bit-stability claim. Every result records how it was made | G0 property laws and G3 metamorphic/tolerance gates on the envelopes; G5 applies only to the RECORDING (the mode tag itself must be deterministic) | `ExecMode::Fast`; the mode is stamped into every event so downstream consumers can refuse to launder a Fast result into a bit-stable claim |

## Rules of use

- The class binds the WHOLE public claim surface of the crate. Mixed
  crates scope the exception explicitly in the CONTRACT (e.g. one Fast
  lane inside a Deterministic crate) and the exception is part of the
  no-claim boundary.
- `Deterministic` crates route every transcendental through
  `fs_math::det` and register in `LIBM_DOCTRINE_CRATES`
  (`xtask check-libm`). Platform libm is admissible only at
  `DeterministicPerIsa` and below (CONVENTIONS, bead lyms).
- Class PROMOTIONS (`DeterministicPerIsa` → `Deterministic`) shift
  last-ULP outputs: every golden in the crate re-freezes under the
  golden-bump protocol in the same change.
- Conformance suites (`fs-casebook`) record the evidence per case; the
  class decides which tolerance specs are admissible: `Exact`/`Ulps(0)`
  claims belong to the deterministic classes, envelope specs
  (`RelativeLe`/`AbsoluteLe`) to any class, and a `Fast` crate claiming
  `Exact` cross-run stability is lying about its class.
- Gauntlet tier names in tests/docs are the "definition of done"
  (AGENTS.md): a determinism claim without its G5 lane is documented as
  targeted, not achieved.

## Tooling

The G5 audit harness is `fs-detaudit`: worker-matrix bit-identity audits
with first-divergence localization, the cross-ISA divergence classifier
(report of record: `docs/G5_CROSS_ISA_REPORT.md`), and measured ExecMode
fast/deterministic deltas.
