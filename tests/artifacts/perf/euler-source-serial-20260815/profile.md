# Euler source-generation live profile — 2026-08-15

This profile ranks the exact serial mechanics/audio source path before the
32-pass optimization campaign.  It is evidence for optimization selection,
not a product artifact or a broad machine benchmark.

## Scenario fingerprint

- Host: `yto` (`ns3365308`), Intel Xeon E3-1245 v2, 4 cores / 8 threads,
  Linux 6.17.0-8-generic.
- Process: PID 3277189, one thread pinned to CPU 7, approximately 99.9% CPU,
  no swap and no physical reads or writes during the sample.
- Binary SHA-256:
  `a27cfa401ec0d3616da4c546c238dc97667911acc74d16a101f3f2886cbd21fe`.
- GNU build ID: `90b64d1faab70bfb62a0ab4052b5eff672df2918`.
- Toolchain: `rustc 1.99.0-nightly (3659db0d3 2026-07-05)`.
- Workload: 6.144 MHz mechanics, 2 s physical preroll, 8 s published source,
  48 kHz controls, 0.12 rad initial inclination, 10 nm / 2 nm self-affine
  disc/support surfaces, H=0.75, cycles 4..256, one 16x9 SPP1 proof frame.
- Preroll baseline: 12,288,000 accepted steps in 10,643.9 s, or about
  1,154.5 steps/s and 1.478 wall-hours per simulated second.

## Sampling method

Linux `perf` was unavailable because `perf_event_paranoid=4`.  Read-only
sampling used `bpftrace` at 199 Hz for 12 seconds against the existing process.
The usable resolved-leaf population was 2,388 samples; an independent 20-second
run observed the expected 3,980 total ticks.  Percentages below are lower-bound
leaf shares because related callees can appear under several symbols.

```text
sudo env BPFTRACE_MAX_MAP_KEYS=16384 bpftrace -q -e '
profile:hz:199 /pid == 3277189/ {
  @leaf[usym(reg("ip"))] = count();
  @samples = count();
}
interval:s:12 { print(@samples); print(@leaf); exit(); }'
```

## Ranked hotspots

| Rank | Hot path | Samples | Leaf share | Direct interpretation |
| ---: | --- | ---: | ---: | --- |
| 1 | Deterministic/libm math, chiefly FMA fallback | 626 | 26.2% | Recovered caller chain reaches `axisymmetric_mass::arc_integral`; immutable profile mass is recomputed during every bind. |
| 2 | Float/`Debug`/String formatting | 432 | 18.1% | Contact request/receipt identities and full production checkpoint fingerprints format structured state as text. |
| 3 | Allocation, free, and `memcpy` | 400 | 16.8% | Per-step template/nested-state cloning and formatted identity preimages allocate and copy. |
| 4 | BLAKE3 compression/update | 378 | 15.8% | Contact and checkpoint identities are recomputed inside the accepted-step loop. |
| 5 | Axisymmetric validation excluding callee math | 197 | 8.2% | Global chart construction validation, feature intersections, profile identity, and mass admission repeat during binding. |
| 6 | Surface trace filtering | 87 | 3.6% | Roughness excitation is a smaller, genuine physics cost and is not the first target. |

Static call-path inspection sharpens the first row: one accepted compliant
midpoint step performs ten complete `validate_profile_with_cx` passes and two
axisymmetric mass integrations.  The same step performs two redundant full
checkpoint validations plus the one successor fingerprint needed for the
accepted state.

## Optimization order and falsifiable gates

1. Seal the immutable profile once per trajectory run and reuse the exact
   admitted mass properties and chart query authority.  Preserve actual
   support/curvature/contact arithmetic, ordering, cancellation checkpoints,
   identities, and the checked public API.  This attacks the measured 8.2%
   direct validation plus much of the 26.2% math envelope.
2. Add a private trusted-trajectory checkpoint path that validates the entry
   checkpoint once while retaining canonical successor fingerprints and all
   checked public entry points.
3. Remove formatted identity preimages and hot-loop cloning only after their
   individual output equivalence and timing are measured.

Every accepted change must pass focused correctness checks and an alternating
same-host short-run A/B.  Physics fields, branch topology, final checkpoint,
modal audio, and retained artifact bytes must match exactly.  Neutral or slower
changes do not land.

## Pass 1 result — sealed profile admission

Commit `73a888942ebe3f30faaafcbf9c4f86f17c7b639a` admits the immutable
axisymmetric profile, density, mass properties, and model mass agreement once
per cinematic trajectory and reuses the admitted mass properties in the
start/midpoint bindings. Public checked entry points remain unchanged.

The exact focused test ran remotely on allowed worker `vmi1293453`:

```text
cargo test -q -p fs-euler-disc-e2e --features cinematic-render --lib \
  cinematic_fixture::tests::g0_parameterized_fixture_advances_real_coupling_into_render_audio_controls \
  -- --exact --nocapture --test-threads=1

test result: ok. 1 passed; 0 failed; 162 filtered out
```

With `FS_CINEMATIC_PASS1_AB=1`, the same test used the production 6.144 MHz
rate, warmed both paths for 32 accepted steps, and alternated six checked and
admitted runs of 256 accepted steps with 128:1 control reduction. Every paired
trajectory compared exactly equal.

| Round | Order | Checked (s) | Admitted (s) | Ratio |
| ---: | --- | ---: | ---: | ---: |
| 0 | checked first | 3.209064 | 2.888960 | 1.1108x |
| 1 | admitted first | 2.982242 | 3.008055 | 0.9914x |
| 2 | checked first | 3.420870 | 3.091011 | 1.1067x |
| 3 | admitted first | 2.903053 | 3.294085 | 0.8813x |
| 4 | checked first | 3.007315 | 2.987226 | 1.0067x |
| 5 | admitted first | 3.180566 | 3.014425 | 1.0551x |

Totals were 18.703110 s checked and 18.283762 s admitted, an aggregate
**1.0229x** speedup; the median of paired ratios is approximately **1.0309x**.
The focused debug timing is noisy and establishes only a modest retained win.
It is not a release/full-source throughput claim; the final integrated release
pass must remeasure the complete source workload.

## Pass 2 result — trusted internal successors

The admitted serial runner now validates the caller-provided start checkpoint
once, then trusts only successors it constructed and sealed itself. Public
checked entry points remain checked. Version, input, gas, surface, time, scalar,
and transaction checks are unchanged, as is the one required successor
fingerprint per accepted step.

The same exact focused remote test passed on allowed worker `vmi1293453`. With
`FS_CINEMATIC_PASS2_AB=1`, it warmed both paths and alternated six checked and
admitted runs of 256 accepted steps at 6.144 MHz. Every paired trajectory was
exactly equal.

| Round | Order | Checked (s) | Admitted (s) | Ratio |
| ---: | --- | ---: | ---: | ---: |
| 0 | checked first | 3.006497 | 2.274077 | 1.3221x |
| 1 | admitted first | 2.998740 | 2.204006 | 1.3606x |
| 2 | checked first | 2.822939 | 2.000070 | 1.4114x |
| 3 | admitted first | 2.900341 | 1.979522 | 1.4652x |
| 4 | checked first | 2.795126 | 2.014507 | 1.3875x |
| 5 | admitted first | 2.804001 | 2.001701 | 1.4008x |

Totals were 17.327644 s checked and 12.473884 s admitted, aggregate **1.3891x**.
Normalizing by Pass 1's 1.0229x checked/admitted ratio gives about **1.3580x**
for the new trusted-successor path. This is still a focused debug measurement,
not the final release/full-source throughput claim.
