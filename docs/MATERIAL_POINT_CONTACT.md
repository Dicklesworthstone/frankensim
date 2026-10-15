# Source-bound contact in coupled impact dynamics

`fs_contact::normal_patch::source::SourcePointContact` binds the existing
`BoundNormalContactModel` to an explicit local relative-gap geometry and a
bounded indentation/rate/step envelope. Density, elasticity, yield data,
ordered surface identity, adhesion and the selected constitutive model remain
in the original immutable material/interface bundles. A material name does not
select a force law. The model card selects elastic Hertz or Hunt--Crossley;
the geometry selects sphere/plane or the actual elliptic paraboloid.

The Hertz coefficient is obtained from the existing finite-patch evaluator,
not from a second elliptic solver. Compilation checks the maximum indentation
at zero and both extreme rates against the original applicability rules. In
particular the separating Hunt--Crossley force factor must remain positive;
invalid contact states are not rescued by clamping tensile force to zero.
The retained closure-corner receipt records the original pressure footprint,
uncertainty and authority. These are model checks, not measured qualification.

## One mechanical equation

Use `fs_couple::render::plate::impact::material::MaterialContact::new` to bind
that source to a unique contact key, a complete mass-normalized contact map
and a reference gap. For two translating bodies of masses ma and mb, the map
is `[1/sqrt(ma), -1/sqrt(mb)]`; with additional modal bodies, include their
actual contact traces and zeros for unrelated modes. Source-derived body mass
can come from the existing `fs-mbd-material::resolve_body` owner using the
same material cards. The geometry owner is responsible for signs and work
conjugacy; a list of finite numbers does not certify a physical contact map.

`ImpactSystem::new_with_material_contacts(bodies, contacts, pads, volumes,
config)` lowers every contact into the existing `fs-dcontact::Obstacle`:
K comes from the source, exponent is 3/2, and the interface card supplies chi.
The existing contact potential and passive dissipative port enter the original
`fs-phs` solve together with all body, felt and cavity storage. There is no
lagged force injection, independent collision time integrator, or caller
stiffness/restitution override. Duplicate keys and incomplete maps refuse.

The original source contacts stay attached to the runtime and survive
`prepare`, `prepare_analytic`, and `into_reference`. Initial and candidate
endpoint states are checked, and the nonlinear callback checks the actual
Gonzalez discrete-gradient rate used by its contact-loss port. Even elastic
source contacts use that same controlled callback with zero extra resistance.
This preserves one equation while preventing a derived velocity outside the
source box from being silently accepted. The fixed rate envelope is applied
to dry states too; it is a conservative operating box, not contact detection.

Accepted-state publication uses the existing common energy/history/clock gate.
Source refusal or cancellation leaves the accepted state unchanged. Prepared
success-path checks use bounded existing contact maps without allocating.
Newton trials outside the admitted domain can refuse; reducing the step or
choosing a physically supported model is explicit caller work.

`material_contacts()` exposes the immutable sources.
`material_contact_observation(i)` evaluates endpoint pressure, footprint,
force and instantaneous loss through the original finite-patch law. It is a
cold allocating observation. Its power-times-step is NOT the integrated
step loss: use `ImpactFrame::dissipated_energy_j` and its uncorrected balance
for the numerical owner's accepted work accounting.

## Focused verification

```
cargo test -p fs-contact --test source_point
cargo test -p fs-couple --test material_impact
cargo test -p fs-couple --lib render::plate::impact::material::tests
```

The fixtures are explicitly synthetic. They exercise original-law parity,
elliptic geometry, source-only stiffness/density/damping substitutions,
reference/prepared/analytic execution, rebound, energy accounting, immutable
source retention and refused transactions. They are not measurements of real
metals or a claim of an experimentally realistic complete instrument.

This path assumes frozen isotropic material states, nonadhesive local
half-space contact and fixed contact maps. It does not implement collision
detection/CCD, moving contact geometry, finite-layer corrections, plasticity,
friction, wear, temperature evolution or melting. Contact energy is shared
with the existing mechanics; no new thermal destination for dissipated work
is inferred.
