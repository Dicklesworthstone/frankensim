# Electrical interface heating and spatial beam response

`fs_couple::strip_enthalpy::bridge::electrical::step_with_contacts_candidate`
connects the existing source-pinned electrical-contact network to the original
finite-storage circuit, spatial Maxwell beam and enthalpy solve. Interface
resistance changes current and heat; the declared A/B partition changes each
body's temperature, relaxation, compatible curvature and sag.

## Existing owners, one accepted endpoint

Prepare an `ElectricalContactNetwork` containing the existing bulk material
network and `MaterialElectricalContact` records. The accepted circuit must
already retain those contact IDs and material sources. The prepared network
retains both bulk and interface source claims across trials and accepted steps.

```rust,ignore
use fs_couple::strip_enthalpy::bridge::electrical;
let prepared = contact_network.prepare(&thermal, cancelled)?;
let candidate = electrical::step_with_contacts_candidate(
    &thermal, &beam, &bridge_bindings, &circuit, &prepared, &held_ports,
    target_s, gravity_m_s2, &additional_downward_n_m,
    electrical_controls, &energy_guesses_j, cancelled,
)?;
let (thermal_next, beam_next, circuit_next, report) = candidate.into_parts();
```

Both bulk-only and contact-capable entry points share beam admission, physical
Maxwell heating and combined energy checks. They delegate to the existing
contact-capable circuit driver, whose empty-contact path uses the original bulk
circuit update. No second nonlinear, circuit, constitutive or thermal solver is
introduced. Every trial begins with the same accepted charge/flux, beam history
and caloric state. The final trial carries all three endpoints into acceptance.

The thermal links supply fixed area and source A/B order. The contact's explicit
heat fraction remains a model input, not a guessed 50/50 or effusivity split.
Bulk and contact destinations come from the same prepared source binding;
callers do not maintain a second independent heat map. Interface source records
and losses are retained in `report.electrical.contacts`; exact receiving cells
and fractions are in `report.contact_heat_destinations`.

`deposited_joule_j` includes bulk and interface heat once. Physical Maxwell heat
is separate. Numerical mechanical damping is not physical heating. Independent
boundary heat and mechanical load work retain their original calculation. The
final gate checks electrical, elastic and thermal storage plus numerical loss
against electrical supply, load work and external heat. Missing or changed
contact sources, work exhaustion and late cancellation return no endpoints.

## Verification

Five new G0/G1/G3 regressions use the original component implementations:
analytic finite-capacitor heating for both charge signs and three partitions;
interface-only resistance changes affecting sag; partition-only changes with
unchanged charge; endpoint resistance and owner replay through unloading; and
source-switch/budget/final-cancellation refusal. Existing bulk tests are unchanged.

```sh
rch exec -- cargo test -p fs-couple --test bridge_electrical_contacts \
  --test bridge_electrothermal --test circuit_electrical_contacts
```

Native execution was unavailable in the editing environment: DSR, its documented
fallback, RCH and Cargo could not start. Independent Python/SciPy equations
exercised 36 parameter cases, six analytic partition/sign cases and four accepted
loading/unloading steps. These are reference equations, not Rust execution or
experimental material validation. The source module loader also uses ordinary
path modules rather than including files with inner module documentation as
macro item fragments; public exports and state ownership are unchanged.

## Boundaries

The existing fixed-geometry, quasistatic, small-strain, single-solid-phase beam
and caloric limits remain. Contact resistance may depend on admitted thermal
coordinates; area, heat fraction, circuit topology, C/L and conductor geometry
stay fixed in an interval. No pressure-dependent area, adhesion, opening, arc,
contact energy storage or field force is inferred. The bulk-only entry point
continues to reject contact-bearing circuits without their heat-routing owner.
