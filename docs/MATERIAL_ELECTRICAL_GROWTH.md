# Finite electrical supply and atomic spatial deposition

`GrowingBridge::step_electrical_candidate` joins the existing electrical/beam
feedback interval and the existing spatial birth transaction. A rejected final
birth cannot consume capacitor charge, age a Maxwell history, or accept a heat
field while leaving the beam unchanged.

```rust,ignore
let (next_specimen, next_circuit, report) = specimen.step_electrical_candidate(
    &circuit, &prepared_material_network, &held_ports,
    target_s, gravity_m_s2, &old_cell_line_loads_n_m,
    electrical_controls, &old_cell_energy_guesses_j, cancelled,
)?;
```

All inputs describe the pre-event topology. The original electrical/beam driver
computes one trial-consistent circuit, beam and caloric endpoint. The original
thermal scheduler then applies endpoint events; the original bridge appends
stress-free tip material in schedule order. Neither electrical nor mechanical
history is aged during birth. All endpoints are returned only after the complete
birth batch, material admission, cancellation and combined energy checks succeed.

The report retains the original electrical/beam interval and a separate birth
audit. The final check reconstructs storage across the entire transaction:

```
delta electrical + delta elastic + delta thermal + numerical mechanical loss
- electrical supply - mechanical load work - independent boundary heat
- incoming enthalpy - zero-duration event external heat = 0.
```

Joule and physical Maxwell losses are internal transfers; neither is added again
to this total. Passive thermal bodies retain incoming enthalpy without becoming
beam segments or electrical stores. `birth_equilibrium_pending` indicates that
new tip material has not equilibrated under the enlarged mechanical load yet.

Circuit topology, charges and fluxes do not change during zero-duration births.
New material is not automatically a resistor. The next interval requires a
complete material-network binding for the expanded thermal topology and the
returned circuit, not a reconstructed initial capacitor. Ordinary preparation
is an explicit source-selection epoch; carry explicit property pins when source
continuity is required. `extend_thermal` is not an electrical topology migration
API. Existing bulk-only electrical/contact restrictions remain unchanged.

Four regressions compare original-owner replay followed by another coupled
interval, a final invalid birth after a valid append, final cancellation/replay,
and passive incoming enthalpy without an electrical debit. They use original
material, circuit, enthalpy, nonlinear and bridge components, not mock steppers.

```sh
rch exec -- cargo test -p fs-couple --test bridge_electrical_growth \
  --test bridge_electrothermal --test bridge_enthalpy_growth
```

Native execution remains unverified: DSR, its documented fallback, RCH and Cargo
are absent in the editing environment. Independent numerical component checks do
not execute or validate the Rust transaction and cancellation paths.

This inherits the single-layer clamped seed, quasistatic/small-strain beam,
fixed thermal geometry and single-solid caloric restrictions. No new constitutive
or nonlinear solver, contact detection, circuit-topology transfer, melting,
impact or instantaneous birth-equilibrium model is introduced. Split intervals
at scheduled events and electrical or mechanical source discontinuities.

## Contact-heated spatial growth

`GrowingBridge::step_electrical_with_contacts_candidate` takes a retained
`PreparedElectricalContactNetwork` in place of `PreparedMaterialNetwork`.
The contact-capable interval uses the original source-pinned interface law and
explicit A/B heat partition. Both electrical growth entry points then call one
private completion path: scheduled births, shared material admission and the
combined electrical/elastic/thermal energy check are identical.

```rust,ignore
let (next_specimen, next_circuit, report) = specimen.step_electrical_with_contacts_candidate(
    &circuit, &prepared_contact_network, &held_ports,
    target_s, gravity_m_s2, &old_cell_line_loads_n_m,
    electrical_controls, &old_cell_energy_guesses_j, cancelled,
)?;
```

A passive contact receiver is heated by the original contact partition but does
not become a beam segment. A failure in the last scheduled birth discards that
receiver's proposed temperature change as well as the circuit debit, beam update
and earlier appends. The returned interval report retains the actual interface
heat destinations. No extra circuit step occurs during zero-duration births.

The caller still supplies a complete successor network after topology changes.
This API neither creates new circuit branches nor migrates contact-source pins
across topology epochs. Interface area, partition and electrical topology remain
fixed during each interval; this is not a moving-contact or impact model.

Three further Rust regressions exercise original-owner replay, a late invalid
birth after a valid append, and final cancellation/replay with real contact heat
in a passive receiving body. They are authored but not executed in this session;
Cargo, rustc and RCH are unavailable. Patch checks are not native validation.
