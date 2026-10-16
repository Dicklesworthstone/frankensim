# Hereditary drumhead bending, not an audio fade

`--head-relaxation material.fshr` supplies a generalized-Maxwell bending law
for each drumhead. The existing equilibrium elastic FEM, installed tension,
head basis, geometry, and contact projections are retained. The material adds
internal strain history to the SAME mechanical time solve. It changes actual
restoring forces and removes physical energy; no pressure filter or per-mode
output envelope stands in for that loss.

```sh
# Explicitly illustrative inputs, NOT identified polyester measurements.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 \
  --drum-spec crates/fs-couple/examples/percussion/estimated-relaxing-drum.fsd \
  --head-relaxation crates/fs-couple/examples/percussion/estimated-head-relaxation.fshr \
  --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 2 > material-head.csv)

# Same material mechanics, complete snare bank and two exterior receivers.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 \
  --drum-spec crates/fs-couple/examples/percussion/estimated-relaxing-drum.fsd \
  --head-relaxation crates/fs-couple/examples/percussion/estimated-head-relaxation.fshr \
  --head-stretching --analytic-newton --impact-substeps 8 511 --cavity-modes \
  --microphone-right -0.12,0.05,0.4 > material-snare.wav)
```

These are supported invocations, not claims of completed or calibrated renders.
The supplied drum file MUST have zero damping ratios on both heads, so material
loss is not counted again as authored modal damping. Both Young moduli must
match their Maxwell equilibrium moduli exactly. No parameter is silently reset.
Unselected commands preserve the previous model. The `drum-modal` image does
not have material memory and rejects this option; use `drum` instead.

## Complete SI input

The bounded UTF-8 file starts with `frankensim-head-relaxation-v1`. Blank lines
and `#` comments are allowed. Exactly one initial-state row, frequency-band
row and material row for each head are required. Branches follow their head.
Unknown, duplicate, missing and nonfinite parameters reject without defaults.

```text
frankensim-head-relaxation-v1
initial,relaxed
band_hz,20,2000
head,batter,4000000000
branch,batter,1000000000,0.001
head,resonant,4000000000
branch,resonant,600000000,0.002
```

`head,name,E_equilibrium_Pa` and `branch,name,E_branch_Pa,tau_seconds` use the
existing `fs-material::visco::GeneralizedMaxwell` admission. Equilibrium modulus
and relaxation times are positive; branch moduli are nonnegative. A material
with no branches is explicitly elastic. There are at most eight branches per
head and 256 added scalar memory states overall; every retained bending basis
row needs one state per positive material branch. Overflow refuses rather than
truncating either spectrum or mechanics. A zero branch adds no state.

`initial,relaxed` equilibrates each added branch at the declared initial
geometry. `initial,unrelaxed` means zero prior viscous strain, so an initially
displaced head stores additional branch energy. Both are physical initial
conditions, not decay presets. With initially flat resting heads they coincide.
This is cold construction, not an unaccounted in-flight material replacement.

`band_hz,lower,upper` declares the material-use window. The requested head-mode
window and actual equilibrium frequencies must fit. A conservative upper bound
on instantaneous head stiffness must fit both that window and the original
mechanical Nyquist guard. Nonlinear harmonics and impact transients may still
leave the retained band; this is not an out-of-band accuracy certificate.

## Energy-consistent spatial coupling

For the actual modal basis Phi, the adapter assembles the existing DKT bending
matrix with ZERO prestress and projects `K_b = Phi^T K_DKT Phi`. It does not
subtract nearly equal tension-dominated matrices or apply material loss to
installed membrane tension. A Cholesky factor `K_b = L L^T` supplies the complete
strain-energy map `L^T q`; all cross-modal terms are retained.

Each supplied modulus ratio E_j/E_equilibrium scales that energy map and is
passed to the EXISTING `fs-phs::with_relaxation_branches`. The owner stores
`H_j=(sqrt(E_j/E_equilibrium)*row.q-z_j)^2/2` and evolves internal strain with
resistance `1/tau_j`. The branch stiffness and loss are thus the same Maxwell
complex modulus, not a loss factor frozen at one modal frequency. The analytic
Jacobian differentiates that exact storage; ordinary finite differences and
bounded substeps retain the same equation. No second integrator is introduced.

Material states follow mechanical and existing felt/Kelvin states, with no new
external force or direct gas/radiation participation. Sticks, carried wires,
compliant pads on supported commands, cavity losses, and prescribed-vent/mono/
stereo observation keep their existing physical addresses and clocks. Refused
steps roll back the complete material state together with the other mechanics.

CSV adds `head_memory_energy_j` (ALREADY included in total storage) and
`head_relaxation_power_w` (an endpoint diagnostic, NOT interval loss). The
ordinary accepted loss/work ledger includes material dissipation exactly once.

## Deliberate limits and tests

This is proportional isotropic LINEAR BENDING relaxation at fixed conditions.
Installed tension is held, and optional geometric membrane stretching remains
elastic with the supplied equilibrium modulus. Nonlinear membrane viscosity,
stress relaxation of the tuning preload, separate shear/bulk spectra, changing
temperature/moisture and a fully polymeric large-strain law are not implemented.
The example coefficients do not acquire measured/coupon authority through import.
Radiation feedback, acoustic bandwidth and real-time qualification are unchanged.

Six library regressions cover initial history/energy, analytic owner parity,
continuous relaxation refinement, joint felt/contact work, substep rollback and
empty selection. Five example regressions cover strict material inputs, spatial
bending work and tension separation, changed real contact-driven head motion,
complete snare/air/acoustic addresses, and elastic/numerical-image boundaries.

```sh
cargo test --release -p fs-couple --lib render::plate::impact::relaxation::tests
cargo test --release -p fs-couple --example percussion head_relaxation::tests -- --test-threads=1
```

Native execution is required to establish these tests pass. Independent matrix
arithmetic and source checks are not Rust compilation or audio-render evidence.
