# fs-couple: legacy tests

The test targets in this folder test the modules in `src/legacy/`, which no live code uses. Among
them are the regression references for the composition kernel. For example, `tribothermal` and
`rigid_damping_heat` hold the graph tests that check a kernel graph against the hand-coupled
module.

Each target keeps its name. Cargo.toml gives its path in a `[[test]]` entry, so
`cargo test -p fs-couple --test tribothermal` runs it as before.

Fixtures and helpers that other tests share stay in `tests/`, and these targets reach them with
`../`. They include `support/`, `material_exchange/`, `electrical_contacts/`, `circuit_enthalpy/` and
`receipts/`.

Two targets whose subject is live stay in `tests/`, although each uses a module from `src/legacy/`
as a helper: spatial_conduction_enthalpy and spatial_induction_rotor.
