# Supported masses: preload, release and resonator contact

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/supported-mass-preload.performance \
  /tmp/supported-mass.wav --block 37
```

A 40-gram mass is held against a two-mode resonator by an authored 2 N load.
Its restoring force comes from a **declared** 400 N/m spring to a fixed support,
not a fabricated modal frequency. Only the resonator emits observer pressure.
The complete spring/contact equilibrium is established before sample zero;
the load releases before sample 37, returns at 1 N before sample 1201, then
releases again before sample 1800. The 4801-sample window retains its short final
callback. Removing the contact stiffness makes the unforced receiver silent.
All parameters are authored examples, not measured specimen data.

## Explicit physical mass initialization

The existing reader adds one voice form without changing its schema versions:

```text
voice free-mass-preload 1 PORT_COUNT
mass MASS_KG 0 0
port INITIAL_FORCE_N SHAPE
```

Repeat the port row as declared. The mass record retains the physical units of
`free-mass`: kilograms, metres, metres per second. Both initial position and
velocity must be zero for preload; nonzero motion is rejected, never discarded.
Every other component must request `static-preload` (elastic) or
`free-mass-preload` (translational). Ordinary `free-mass` still retains its
supplied motion and does not request equilibrium.

The modal coordinate remains `q = sqrt(mass) * x`. Physical actuator and
attachment shapes are still in `1/sqrt(kg)`; unit translation for a 0.04 kg mass
uses shape 5. Changing the mass to 0.16 kg requires changing its unit-translation
shapes to 2.5 as well. No shape is silently normalized or inferred. A free
coordinate still has zero intrinsic stiffness, zero intrinsic damping and no
direct acoustic transfer. All actual spring/dashpot properties live in the
connection records.

A fixed support is represented by a zero attachment shape. In the example,
`left 0 5` and `right 0 0` bind the spring to the moving mass and an immobile
reference, respectively. The nonzero stress-free extension is retained.
Other legal supports include chains of free masses grounded by springs, and
free masses connected to elastic components with positive intrinsic stiffness.
A dashpot alone cannot provide static support.

## One equilibrium owner, bounded extra work

`CoupledModalSystem::initialize_static_equilibrium` and
`initialize_contact_equilibrium` share the same static response. Elastic
coordinates keep the existing diagonal-compliance/connection-factor path.
Zero-stiffness coordinates are resolved through the declared spring support
matrix, using the existing `fs-la` Cholesky implementation. The additional dense
factor has at most as many rows as bilateral connections (currently at most 64),
not one row per retained elastic mode. No full modal stiffness is assembled.
The original all-elastic arithmetic and the dynamic integrator are unchanged.

Let `n` be total modes, `k` bilateral connections and `z` free coordinates.
For `z > 0`, the static setup screen charged to `coupling_limits` is

```text
n * (k + 1)^2 + z * k^2 + z^2 * k + z^3
```

It must fit the existing `MAX_SETUP_TERMS`; the allowance is not enlarged.
The normal-contact preload retains its separate contact-set work limits.
The equilibrated free-support factor also refuses squared pivots no larger than
`256 * f64::EPSILON * z`: this avoids using rounded singularity as a hidden
pose constraint. It is a conservative floating-point degeneracy screen, not
an interval rank certificate, and can reject extremely ill-conditioned supports.

Force residuals are checked both in the condensed equations and on the actual
rounded component states. All original component state limits, bilateral force
limits, contact laws, penetration limits and total-energy limits remain active.
Cancellation or refusal publishes no initial state. Static preload does not
reconstruct actuator work before the declared window or simulate a warm-up.

## Boundaries that still refuse

Every free coordinate must be constrained by the bilateral spring network,
possibly through elastic components. A free assembly with a remaining rigid
translation is refused even under zero net load; no arbitrary pose is selected.
Redundant links cannot stand in for an independent missing support. A body held
up **only** by unilateral contact remains outside this inverse-stiffness path;
contact may alter the preload but is not used to make a singular bilateral
network invertible. Friction-loaded static equilibrium remains unsupported.
The independent v1 path cannot preload an unsupported free mass.

The command uses its existing input hashing, physical-force scheduling and PCM
writer. No new output gain, contact law, integrator or format version is added.
Ordinary retained free-mass inputs keep their original behavior. The separate
high-rate `--decimate` path remains available under its existing clock rules.
This feature is fixed-axis translation and authored reduced acoustics, not
rigid-body rotation, collision discovery, experimental validation or a measured
real-time guarantee.

Focused native checks:

```bash
cargo test -p fs-couple --test supported_free_preload --test music_render_supported_mass
cargo test -p fs-couple --test modal_contact_preload --test modal_free_mass
cargo test -p fs-couple --test coupled_render --test music_render_preload
```

These Rust tests were added but not executed in the authoring environment,
which lacks Cargo/rustc. Independent Python full-stiffness and dynamic numerical
comparisons do not establish that the Rust implementation compiles.
