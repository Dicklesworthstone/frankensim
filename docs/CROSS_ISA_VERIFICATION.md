# Cross-ISA Verification Playbook (x86-64 ⇄ aarch64)

FrankenSim promises **deterministic bit-stable behavior** and **roofline-honest
performance on both reference ISA families** — Apple aarch64 (16 KiB pages, always
has native FMA) and x86-64 Threadripper/EPYC (4 KiB pages, FMA only with the target
feature). The dev + DSR-CI fleet is aarch64-heavy, so x86-specific issues are
*structurally invisible* until someone builds or runs on x86. This is how to verify
cross-ISA correctly — and, just as important, how to avoid the traps that make you
diagnose the wrong thing.

Related: `docs/CI_GATES.md` (the automated gates), `docs/GOLDEN_POLICY.md`,
`memory/linux-perf-machines` (SSH/setup for the Threadrippers).

---

## 1. The two-layer gate: what catches what

| Layer | What it catches | Needs x86 hardware? |
| --- | --- | --- |
| **Pre-commit x86 cross-COMPILE** (`scripts/ci/x86_cross_check.sh`, bead `ebro`) | x86 compile + **test-compile** breakage (missing `#[cfg]` arms, renamed symbols not followed, `cargo check` lib-only misses test targets) | **No** — cross-compiles from the M4 |
| **Runtime x86** (the Threadrippers) | page-size/alignment asserts, SIMD-capsule bit-divergence, golden reproduction, perf attainment, fsqlite/storage behavior | **Yes** |

`cargo check --workspace --all-targets --target x86_64-unknown-linux-gnu` from an
aarch64 box is the cheap first line — `--all-targets` is load-bearing (a plain
`cargo check` never compiles test targets, so it misses breaks like a test using a
type that doesn't `impl Debug`). Only *runtime* classes need a real Threadripper.

---

## 2. Running on the Threadrippers (trj / ts1 / ts2)

SSH aliases + specs live in `memory/linux-perf-machines`. Essentials:

- **Check load first.** `ssh <host> uptime` — ts1/ts2 are usually quiet; trj is
  often saturated (load 200+). A roofline number from a loaded box is noise.
- **Get the code via `git clone`, never rsync/tar the working trees.** frankensim +
  its **8 sibling path-deps** must be siblings for cargo to resolve:
  `frankensim, franken_networkx, franken_numpy, frankenscipy, frankensqlite,
  frankentorch, asupersync, franken_engine`. `git clone --depth 1
  https://github.com/Dicklesworthstone/<repo>.git` each into one dir (~1.4 GB
  tracked). **frankensqlite's working tree is ~15 GB of untracked test corpora** —
  a naive tar/rsync times out just walking it; the clone leaves that behind.
- **`asupersync` is the build long-pole** (~3 min release, `codegen-units=1`).
- **ts2 has a shared `CARGO_TARGET_DIR=/data/tmp/cargo-target`** — concurrent agents
  contend on its lock. Use a private `CARGO_TARGET_DIR` to isolate.

---

## 3. Verification lanes

```bash
# Compile (all crates, all targets) — the ebro gate, runnable anywhere via cross:
cargo check --workspace --all-targets --target x86_64-unknown-linux-gnu

# SIMD capsules are bitwise-correct vs their scalar twin (AVX2/BMI2/etc.):
cargo test -p fs-simd --release           # tier_equivalence_battery
cargo test -p fs-substrate --release      # BMI2 Morton, os_affinity, prefetch capsules

# Determinism: aarch64-frozen golden hashes reproduce on x86 (run per-crate; a
# full-workspace `cargo test` times out under swarm build-lock contention):
cargo test -p fs-fft -p fs-sparse -p fs-la -p fs-topo -p fs-evidence --release

# Perf attainment (release, #[ignore]'d; authority configuration is in §5):
cargo test -p fs-fft --release --test perf_lane -- --ignored --nocapture
```

---

## 4. THE TRAPS (each one cost a real misdiagnosis)

1. **Version skew across machines masquerades as an ISA bug.** The constellation
   path-deps float on HEAD, and different machines' clones drift. A day-stale
   `frankensqlite` clone on the M4 turned an fsqlite HEAD regression into a
   convincing "passes on aarch64, fails on x86" — it was really *old dep vs new
   dep*. **Before believing any ISA split, `git log -1` the relevant path-dep clone
   on BOTH machines.** (This is how bead `u8og` got misfiled as x86-specific.)

2. **Stale goldens masquerade as a flag/change breaking determinism.** A golden
   test failing under a new build flag is often a *pre-existing stale golden*, not
   the flag's doing. Confirm the frozen constant is current (git-blame it against
   the code that feeds it) before inferring a mechanism. (A stale `fs-cheb` golden,
   already out of date after an fs-fft commit, got wrongly blamed on `+fma`
   contracting `a*b+c`.)

3. **The libm-fma trap.** On baseline `x86_64-*` (no `+fma`), `f64::mul_add` lowers
   to a per-element **libm `fma()` CALL** (~1 GFLOP/s), NOT a `vfmadd` — Rust can't
   emit hardware FMA without the target feature, and won't contract `a*b+c` on its
   own. Invisible on aarch64 (native `fmadd`). It caps any mul_add-heavy x86 kernel
   (bead `cwjn`: fs-feec apply 2.6% on x86 vs 44% on M4). **Fix per-kernel** with a
   `#[target_feature(enable="fma")]` capsule (runtime `is_x86_feature_detected` +
   scalar fallback, registered in `unsafe-capsules.json`), NOT a global
   `+fma` — global breaks the runtime ISA-admission/dispatch contract and baseline
   portability. Routing through a dispatched `ops().axpy` does NOT help tiny inner
   loops (indirect-call + per-call feature-detect + scalar-tail-libm swamps the one
   vfmadd). See `.cargo/config.toml` (documents the rejection) and bead `xlvx`
   (fs-roofline's FMA probe uses the capsule pattern).

4. **Page size / allocation alignment differs.** aarch64 (Apple) = 16 KiB pages,
   x86-64 = 4 KiB. A "1 MiB / 4 KiB = 256 pages" assertion is 257 on x86 when the
   heap allocation isn't 4 KiB-aligned. Compute expected page counts from the actual
   start offset — `(start % page + len).div_ceil(page)` — never hardcode. (bead
   `9e6d1eb` fixed exactly this in `fs-substrate`.)

5. **Contention on the shared M4 poisons roofline numbers.** The dev M4 runs the
   whole agent fleet; single-core attainment swings run-to-run (one gate measured
   0.69 / 0.92 / 0.93 across three back-to-back runs) and best-of-N does not rescue
   it. **Measure perf on a quiet Threadripper**, and check `uptime` right before,
   not after.

---

## 5. Perf-baseline and authority inputs

Perf lanes (`fs-fft`, `fs-feec`, …) use `FRANKENSIM_BASELINE_STORE` plus
`FRANKENSIM_FIRMWARE_ID` for historical-axis comparison. The committed files in
`perf-baselines/` are plain promoted candidates: by themselves they are always
report-only and cannot authorize a positive gate.

An authority-admitted run also requires an attested store envelope, a configured
promotion-authority policy, and the retained source-receipt set named by that
baseline. Missing, partial, malformed, denied, revoked, tampered, or cross-machine
authority inputs downgrade the run to measured/report-only; they never silently
fall back to a positive gate. The retained JSON output binds the one frozen
authority decision used for that run. The FEEC and FFT lanes also require
`FRANKENSIM_ROOFLINE_LEDGER=<db>` before they can emit a positive gate. Their
shared recorder commits the exact admission receipt and exact final-gate JSON
atomically; a missing ledger path is report-only and a failed write fails
closed.

The retained-source input is a protected operator inventory declaration, not a
byte-availability oracle: these lanes check exact hash membership but do not
fetch or independently rehash the source artifacts.

For the custom perf lanes, configure
`FRANKENSIM_PROMOTION_AUTHORITY_POLICY=<tsv-file>` and
`FRANKENSIM_RETAINED_SOURCE_RECEIPTS=<strictly-sorted-lowerhex-lines-file>`, plus
`FRANKENSIM_ROOFLINE_LEDGER=<db>` for a citable FEEC or FFT gate.
The `roofline` CLI takes the corresponding `--authority-policy <tsv-file>` and
`--retained-receipts <strictly-sorted-lowerhex-lines-file>` flags and its own
`--ledger <db>` output argument. Authority configuration and durable output are
separate requirements; stdout alone is never citable retained evidence.

Promotion needs at least three spaced probes that mutually agree within a drift
band — **promote only on a genuinely quiet host**; a contention-deflated baseline
poisons future "quiet" runs by flagging them Suspect. `roofline promote` creates a
plain candidate. Attestation and authority admission are separate operations; do
not manufacture an in-tree tag and call it independent authority.

---

## 6. Escalate, don't whack-a-mole

If you find yourself catching one transient x86 break after another, the leverage is
at the *gate*, not the instances. The pre-commit cross-check (`ebro`) closed the
compile-class gap. Runtime-class gaps (page-size, storage, perf) still need a
Threadripper lane — if that recurs, propose wiring it into the nightly/DSR flow
rather than re-catching it by hand.
