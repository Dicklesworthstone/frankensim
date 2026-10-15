# Independent material blades → coupled tube → physical WAV

`paired_reed` composes two independently reduced plate specimens, one shared
slit/contact, the regional waveguide owner, a gas/geometry-derived radiative
termination, and the existing stationary baffled pressure observer and PCM writer.
There is no prescribed note, sample bank, output equalizer or peak normalization.

From the repository root, use two **new** destination paths and the same pressure
full-scale for both runs:

```sh
cargo run --release -p fs-couple --example paired_reed -- \
  crates/fs-couple/examples/paired_reed/equal-blades.paired equal-blades.wav 0.005
cargo run --release -p fs-couple --example paired_reed -- \
  crates/fs-couple/examples/paired_reed/heavier-second-blade.paired heavier-blade.wav 0.005
```

Only the second blade's density differs. Both files retain the same elasticity,
thickness, loss ratio, atmosphere, geometry and blowing-pressure program. Mass
and natural frequencies are derived by the original plate eigensolver; the
shared slit, both pressure faces and bore react to that change during stepping.
The explicit loss ratio is unchanged, so its dimensional viscous coefficient
also changes consistently with mass. No stiffness or pressure is retuned to
hide the substitution. The second file is a **synthetic property experiment**,
not a claim that a real species/alloy has those independent properties.

Change a `blade` row to supply Young's modulus [Pa], Poisson ratio, density
[kg/m³], thickness [m] and viscous damping ratio. Both rows are independent.
Change `gas` to set moist-air Kelvin temperature, static pressure [Pa], and
relative humidity. Those values rebuild a coherent gas used by the nonlinear
source, tube transit/impedance, radiative load and receiver. The library also
accepts other complete `GasState` values and prescribed section-by-section gases;
this small input example intentionally names moist air rather than inventing
its own species parser or thermodynamic model.

`audio` gives mechanical rate and total sample count. `pressure` gives peak
blowing pressure [Pa], attack samples, release-start sample and release-duration
samples. It is a physical piecewise-linear source evaluated at each interval's
midpoint. The trailing window keeps both blade motions, traveling waves and
radiation-load history; source release does not reset them. Increasing pressure
is not guaranteed to produce a stable sustained musical note.

Each blade is a 25 mm × 10 mm clamped-root plate on the existing 6 × 2 rectangle
mesh, with a 0.2 mm total mean rest slit. The existing DKT/modal owners select
the first mode in a fixed search window up to 450 Hz; an absent or unresolved
mode refuses rather than prescribing a replacement frequency. The source
retains individual linear-slope and slit-uniformity limits. The shared contact
is explicitly authored, and is counted once, not once per blade.

The tube is 0.25 m long with a 7 mm radius and an explicit 5 mm integer-transit
length-error allowance. Its termination uses `BaffledRadiationLoad`, not a fixed
reflection or added end-length adjustment. The receiver is on axis, 0.25 m from
that outlet, using the outlet's actual gas and volume flow. The bounded receiver
and load band is at most 1 kHz and 0.05 of the mechanical sample rate. This band
is an approximation allowance, not proof that the nonlinear source has no
higher-frequency content. Distributed wall losses are **not selected** in this
example; the library's `PairedApertureNetwork::with_viscothermal` selects them
through the same owner used by the single-blade and flue systems.

The WAV carries exterior pressure with the requested fixed full-scale, and
stderr reports clip count, derived blade parameters, represented tube length,
remaining stored energy and maximum balance residual. Existing destinations
are refused. No file is created until the whole physical run, receiver and
encoding have succeeded; an I/O failure while writing may leave a partial file.

## Library use and scope

`bernoulli_aperture::paired::PairedAperture::from_plates` accepts two original
`PlateApertureReduction` values, including ones built with `from_material_chart`.
Their separate material receipts, geometry, states and nodal observations remain
attached. `PairedApertureNetwork` composes that source with existing regional
tubes, branches and reactive/viscothermal loads. Independent blade forces have
separate physical work; failed trials publish neither blade nor any wave state.

This is two one-mode linear blades sharing one **uniform effective** slit. It
retains unequal and antisymmetric motion, unlike the existing symmetric-pair
image, which remains available. It does not yet supply independently evolving
hereditary material memory, curved-reed shells, opposing-mesh contact, a conical
staple's pressure recovery, real tone-hole geometry or an unbaffled bell/room.
It is **not a calibrated oboe or bassoon model**, and makes no measured timbre,
hard-real-time, full-band or piano-level-accuracy claim.

## Focused checks

```sh
cargo test -p fs-couple bernoulli_aperture::paired::tests
cargo test -p fs-couple bernoulli_aperture::paired_network::tests
cargo test -p fs-couple --example paired_reed
```

The implementation environment has no `cargo` or `rustc`: these native tests,
the example executable and its WAV output have not been executed there.
Independent Python/SciPy equations checked both blade momenta, contact/flow
closure, regional wave feedback and material-only substitution. An additional
12,000-step reference uses a cubic-Hermite cantilever and compact observer,
**not this Rust DKT reduction or its disk quadrature**; it tests the physical
scale and coupled energy accounting, not implementation or experimental fidelity.
