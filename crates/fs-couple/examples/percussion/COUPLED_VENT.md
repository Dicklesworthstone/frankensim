# Exterior pressure feedback on a drum vent

`--coupled-vent-radiation` connects the existing head/mouth BEM boundary to the
same joint mechanical time owner that advances the heads, wire contacts,
interior acoustic modes and neck slug. It implies `--radiation-feedback`.
This is a stationary, linear exterior impedance coupled to the retained
nonlinear mechanics, not a resolved exterior fluid simulation or a calibrated
instrument.

The accompanying `--cavity-neck-internal` declaration is mandatory:

```text
--cavity-neck-internal radius_m length_internal_m resistance_internal_Pa_s_m3 azimuth_rad z_m
```

The length includes the duct and any declared **interior** end correction; the
resistance includes only retained internal losses. Both **exclude exterior
radiation**, which the BEM-derived load supplies. No correction is guessed from
the radius, subtracted from an effective length, or added a second time. The
numbers below are explicit illustrative inputs, not measurements of a specimen.
The axial coordinate is measured from the top head, as for `--cavity-neck`.

From the workspace root, this selects both stretching heads and a bounded
40..500 Hz acoustic demonstration:

```sh
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 192 20 --cavity-modes --coupled-vent-radiation \
  --cavity-neck-internal 0.003 0.008 1000 0.35 0.08 \
  --radiation-spec crates/fs-couple/examples/percussion/coupled-vent.radiation \
  --analytic-newton --microphone-right -0.12,0.05,0.4 > coupled-vent.wav)
```

The run still has to pass its native modal, contact, BEM resolution, multiport
passivity and independent fit-audit gates. Native run 38391471556 passed the
cavity/neck mechanics but refused the collocated head-mouth power matrix; it did
not validate coupled-vent audio. The variational source image below targets that
failure without loosening the power check. The native reciprocal-core test step
passed on commit `0c8af3ac`; this playback integration still needs its own native
result. A declared band is not a convergence result; failed fitting does not
switch to prescribed flow or invent a replacement impedance.

## Reciprocal source discretization

Coupled vents now use `fs_bem::reciprocal_radiation`: continuous linear pressure
on the retained triangle vertices and constant auxiliary single-layer density
on each panel. The mixed symmetric Galerkin operator uses weak surface curls,
analytic static triangle moments and regular Helmholtz quadrature. Its outgoing
part is integrated with positive antipodal plane-wave quadrature. Reciprocity
and discrete outgoing-flux balance come from this operator, not a symmetrized
impedance, clipped power eigenvalues or added damping.

Both source work and receiver pressure use that same P1/P0 solution. Panel mean
pressure is valid for work against the prescribed constant normal velocity; it
is **not** passed to a different constant-pressure field representation. The
auxiliary density is not labelled as physical mouth velocity. Every source also
checks the mixed solve residual and surface work against outgoing sound flux.
The existing full Hermitian multiport check and all load/observer fit audits are
unchanged. Unvented feedback and prescribed-flow one-way rendering retain their
existing collocation image.

This first variational image admits a band below the enclosing slab's interior
Dirichlet bound; it does not replace Burton-Miller at resonant frequencies. The
example caps pressure-plus-density unknowns at 4096 and retained geometry at
256 MiB, and charges the whole frequency sweep against `max_dense_work`. Angular
refinement 0..2 is selected from the declared band and boundary extent, never
from the sign of a failed power eigenvalue. Mesh and angular convergence remain
required. Independent Python sphere checks reduced impedance error from
8.3..10.7% on 20 panels to 2.3..3.0% on 80; those are equation checks, not native
drum validation or error bounds for an arbitrary instrument.

Finite-point, far-field and adaptive/directional `--microphone-spec` receivers
use that same P1/P0 source. The adaptive path integrates each original vertex
shape and its analytic spatial derivative; it never substitutes panel means or
treats auxiliary density as normal velocity. Original barycentric coordinates
survive every integration subdivision. Pressure and derivative components have
separate discrepancy checks under the existing per-receiver work/depth limits.
Prepared rows are reused for all modal sources at a frequency. Geometry and
source identity, cancellation, exterior winding and actual clearance still
refuse invalid requests before playing.

For close miking, replace `--microphone-right` in the command above with:

```text
--microphone-spec crates/fs-couple/examples/percussion/coupled-vent-microphones.frm
```

That file declares a cardioid 7.45 mm above the reference batter head and an
omni 7.45 mm below its resonant head. Their positions/axes are inputs, not
calibrated microphone capsules. The cardioid uses the calculated fluid particle
velocity, including its reactive near-field component, with no assumed source
direction or distance-based EQ. Both channels observe one loaded trajectory.
No artificial two-sample latency is inserted at these short clearances; the
physical propagation stays in the fitted transfer. The existing band, passive
load and independent observer-fit checks still apply. The new core/routing
regressions and extended close-stereo render test require native execution;
this declaration is not evidence of a completed render or instrument fidelity.

## One physical source map and one playing clock

Each head's original modal normal-velocity row is retained. The additional mouth
row is normalized by its actual polygon area so its integrated outward flow is
`Q = volume_weight * p_neck`, with the same neck coordinate used by the cavity.
The original closed computational boundary is preserved; a prescribed Neumann
mouth closes it acoustically but is not a rigid cap. The geometry adds no neck
inertance or damping.

For every source, BEM pressure is integrated against every source row to obtain
the **full** impedance. Head-mouth cross-reactions are retained, including their
signs. The full Hermitian power form must pass the existing check; positive
self-resistances alone are not sufficient. The passive realization then attaches
through `CavityCoupling::attach_exterior_radiation`, not a parallel solver or a
one-step-delayed pressure force.

The resisting mouth pressure obeys
`F_neck = -volume_weight * p_exterior`, hence
`F_neck * p_neck = -p_exterior * Q`. Its acoustic storage and dissipation are
already part of the joint energy ledger. They must not be added again as an
external force or a second loss. Internal cavity pressure and the mouth's load
pressure remain distinct from pressure observed at a microphone.

The same BEM scene supplies the mono or stereo pressure filters. Both receiver
channels observe one loaded trajectory. Preparing filters neither strikes the
drum nor advances its clock; failed/cancelled mechanical trials retain the
existing rollback contract. Analytic Newton, optional impact subdivision and
player-force programs remain execution controls on that same state.

## Preserved options and refusal boundaries

`--cavity-neck` keeps its original effective-impedance, zero-exterior-pressure
semantics. `--prescribed-vent-radiation` with that declaration keeps the existing
one-way observation model. It cannot be combined with feedback. A coupled
request cannot reuse the old effective declaration, omit a mouth source, expose
an interior standing-wave coordinate as an exterior source, or use a source
column twice.

Coupled vent audio is available on the nonlinear-capable `drum`, `drum-stretch`,
`snare` and `snare-off` microphone/WAV families. It does not silently convert the
linear-only `drum-modal` image. Snare wire count, contacts, second-stick addresses
and any selected head-stretching law remain explicit. The compact-neck chart is
checked against the requested acoustic upper frequency, not only the cavity
basis or the default observer window.

The current mouth must lie inside one fixed planar sidewall strip. Moving
apertures and elastic barrels with necks remain refused. The modal, source-count,
fit-order, pole and time-step limits are unchanged; the variational preparation
has the additional explicit bounds above. Turbulent
jets, finite-amplitude exterior flow, measured material/neck identification,
full-band convergence and real-time qualification remain outside this increment.

## Focused native checks

```sh
cargo test --release -p fs-couple --example percussion coupled_selection -- --test-threads=1
cargo test --release -p fs-couple --example percussion feedback:: -- --test-threads=1
cargo test --release -p fs-bem --lib reciprocal_radiation
```

The new checks cover explicit declaration/admission, omitted and duplicate mouth
columns, accidental interior sources, out-of-chart frequency bands, the actual
two-stick stretching-head/snare/cavity/neck composition, work-conjugate exterior
pressure, total energy and refusal/retry. Two checks additionally execute actual
BEM projection and loaded stereo PCM preparation. The declared passive-load
fixture isolates mechanical coupling; it is not evidence that the separate BEM
fit succeeds. The new playback integration still requires native execution.
