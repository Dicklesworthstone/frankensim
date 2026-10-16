# Distributed wire–head contact, not a noise layer

The `snare` family retains the original double-head drum geometry, film
constants/tensions, stick, volume coupling and prepared modal mechanics. It
adds 20 **independent** tensioned strands below the resonant head. Each has
eight mass-normalized transverse modes and 12 length-weighted contact stations:
160 wire coordinates and 240 wire–head reactions, plus the original stick
contact. All reactions see the same candidate body states in one joint solve.
No wire is driven by a prescribed force, sampled buzz, noise source or EQ.

```
cargo run --release -p fs-couple --example percussion -- snare 4096 > snare.csv
cargo run --release -p fs-couple --example percussion -- snare-off 4096 > snare-off.csv
```

The microphone uses the existing finite-point BEM evaluation and causal filter
pipeline, unchanged. Wire reactions change actual head motion, which changes
its radiated pressure. Direct radiation from the thin wires is **omitted**, not
assigned an arbitrary gain. No radiation force is fed back into the mechanics.

```
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 48000 20 0.08 0.05 0.35 > snare-mic.wav)
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-off-mic 48000 20 0.08 0.05 0.35 > snare-off-mic.wav)
```

The arguments are output frames, full-scale pressure [Pa] and receiver XY Z
coordinates [m], exactly as in `AUDIO.md`. `snare-wav` and `snare-off-wav` use
the separate far-field approximation. The original `drum`, `drum-modal` and
`splash` families are unchanged. No count is silently increased or truncated.
`snare-off` starts with 3 mm clearance instead of 20 micrometres; it is an
explicit comparison of initial installation geometry, **not** a live throw-off
mechanism or a guarantee of separation under every possible excitation.

## What comes from a real object

D'Addario's official PureSound Custom Pro steel product page identifies a
nominal 14-inch, 20-strand steel-coil assembly (CPS1420), also offered in other
strand counts. It describes end plates lifting coils from the rim. It does
not supply the coil geometry, tension, bending law or loss parameters used here:

https://www.daddario.com/en-nl/products/puresound-custom-pro-snare-wire-steel

Only strand count, steel-coil construction and nominal drum size are anchors.
`SnareSet::reference` labels every remaining parameter as an editable estimate:

| Parameter | Explicit estimate |
|---|---:|
| Speaking length / bank width | 300 mm / 40 mm |
| Metal wire radius | 0.15 mm |
| Coil centreline radius / axial pitch | 0.55 mm / 0.85 mm |
| Metal density | 7,800 kg/m³ |
| Installed effective tension per strand | 0.7 N |
| Effective flexural stiffness per strand | 1e-6 N m² |
| Modal viscous coefficient | 4 s⁻¹ |
| Distributed contact K / exponent | 5e8 N/m^(alpha+1) / 1.5 |
| Hunt–Crossley loss coefficient | 0.05 s/m |

Helix arc length determines mass per **axial** metre. It does not identify the
coil's effective bending or tension law; using the straight metal wire's EI as
that law would be a separate, unjustified assumption. These estimates imply
about 2.308 g/m per strand and 13.85 g total wire mass, excluding end plates.
A surveyed specimen can supply different values without changing the solver.

## Generic building blocks

`impact::linear::wire::{HelicalWire, WireSpan, LineContact, film_shapes}` builds
ordinary `ImpactBody` and `fs_dcontact::Obstacle` values. Frequency and sine-mode
normalization come from the existing prestressed-beam and contact owners.
The receiver can be any supplied finite point-major modal shape table; the
film helper evaluates the actual XY location barycentrically instead of snapping
to a nearby mesh node. Receiver and filament have the same positive transverse
axis; closure is receiver displacement minus filament displacement minus gap.
The exact same signed row distributes force, preserving reciprocal work.

Contact K is **per unit line length**. Positive quadrature weights are metres,
so contact refinement does not accidentally multiply the total stiffness.
Wire modes, receiving-surface modes and contact stations must all be refined
independently; 12 stations and eight modes are a bounded example, not a
convergence result. Pinning a wire's ends is an explicit mechanical restriction.
No end-plate flexibility, coil torsion, inter-wire friction, coil-scale contacts,
snare-bed shape or nonuniform preload is inferred.

The shared normal-contact ceiling is now 512, still subject to the original
quadratic setup-work and bounded root/sweep limits. Existing caller limits stay
unchanged; this example explicitly requests 256 modes and 50 million contact
setup terms. A request that does not fit refuses rather than dropping strands.
Scalar contact-law scratch is retained between calls rather than allocated on
every step. This is not a whole-runtime allocation-free claim; other parts of
the original contact path still allocate. No real-time timing was established.

The energy-consistent wire/membrane collision formulation has primary-literature
precedent in Bilbao, Torin and Chatziioannou, *Numerical Modeling of Collisions
in Musical Instruments* (2014): https://arxiv.org/abs/1405.2589 . This code reuses
FrankenSim's existing modal/contact integrators, not that paper's finite-
difference implementation, and does not inherit its verification by citation.

## Executed evidence versus outstanding native validation

Seven library Rust regressions cover mass/dispersion, reciprocal contact work,
contact refinement, actual film interpolation, excitation/backreaction, silence,
resume/refusal and >32-point admission/work limits. Two example regressions
check the explicit configuration and original head assembly with all 241
reactions. **Native compilation and these tests have not been executed here:**
Cargo/Rust are unavailable. Run the focused native tests before trusting audio:

```
cargo test -p fs-couple --test wire_contact
cargo test -p fs-couple --example percussion snare
```

An independent Python/SciPy reference ran 1,024 samples of a 20-strand,
240-point bundle against an analytical single receiving mode. It uses the same
estimated helix mass but a separate fixture damping of 2 s⁻¹. All 240 contacts
entered and released; wires moved and reacted back on the receiver. Disabled
contact left the wire states exactly zero. Maximum per-step energy residual was
2.722e-15 J. This does not execute the Rust root solver, full film eigensolve,
BEM, vector fitter or native audio. A smaller 2,048-step reference also verified
replay and line-strength refinement. No listening or specimen calibration is
claimed. The existing limited head/radiation bands are unchanged: these commands
are not yet a full-band realistic snare or a real-time instrument.
