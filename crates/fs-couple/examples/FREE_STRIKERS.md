# Untethered strikers and free translations

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/free-striker.performance \
  /tmp/free-striker.wav --block 37
```

The example launches a 40-gram mass at 1 m/s across a 0.5 mm gap toward a
supplied two-mode resonator. The striker has **zero stiffness**, not a low
artificial natural frequency. Contact decelerates it and excites the initially
resting receiver. It rebounds and coasts; the later 1 N force at sample 300,
released at sample 340, acts on the striker only. Only the receiver contributes
pressure. Disabling contact, or moving its gap out of reach, produces silence.
The model and contact coefficients are illustrative authored inputs, not a
measured hammer, mallet, material pair, or calibrated instrument.

## Explicit physical mass input

All existing modal performance versions accept this additive voice form:

```text
voice free-mass 1 PORT_COUNT
mass MASS_KG DISPLACEMENT_M VELOCITY_M_PER_S
port INITIAL_FORCE_N SHAPE
```

Repeat `port` for the declared count. The coordinate count must be exactly one.
The mass record uses physical kilograms, metres and metres per second; its
runtime coordinate is `q=sqrt(m)*x` and velocity is `sqrt(m)*v`. Existing `limits`
remain mass-normalized displacement/velocity ceilings, plus joules and pascals.
There is no pressure-transfer field for a free mass. Ordinary `mode` records
still require strictly positive frequency; a zero in an old row is not retyped.

**Port and attachment shapes keep their existing mass-normalized units.** For
a unit physical translation in the declared direction, use `1/sqrt(m)`:
40 grams gives shape 5; 160 grams gives shape 2.5. Changing mass requires changing
its corresponding shapes too. Signed shapes may describe a chosen direction or
leverage, but the reader does not infer geometry or silently renormalize them.
This same map converts applied newtons to generalized force and generalized
motion back to physical attachment displacement, preserving work conjugacy.

The existing force events, bilateral connections, normal contact and supported
friction paths operate on the same accepted states. `free-mass` retains supplied
motion. Assemblies containing free coordinates cannot use the existing static
preload inverse, even if a more general constrained equilibrium might exist.
Such requests refuse rather than inventing a support or selecting a static pose.

## Native API and scope

`ModalAcousticTimeModel::try_free_mass` constructs this one-coordinate physical
image. `try_new_with_free_coordinates` admits an explicitly mass-normalized bank
mixing zero-stiffness and elastic coordinates; zero stiffness requires zero
modal damping ratio and zero acoustic transfer. `try_new` keeps its original
positive-frequency admission. Both paths use the existing exact held-force
transition and transactional state, work/energy and pressure checks.

This supplies fixed-axis translation, not six-degree-of-freedom rotation, moving
attachment directions, collision discovery, rigid impulses, or inferred gravity.
The contact remains compliant and time-discretized. Refine the timestep for
contact accuracy; the linear modal bandwidth guard does not eliminate nonlinear
aliasing. Free-body acoustic radiation and measured real-time performance are
not claimed. Existing input/output hashes and the WAV encoder are unchanged.

Focused native checks (not executed in the authoring environment):

```bash
cargo test -p fs-couple --test modal_free_mass --test music_render_free_mass
cargo test -p fs-couple --test modal_contact --test modal_multi_contact
```
