# Shared material snapshots in coupled networks

Gas-wall conduction consumes the complete material snapshot already resolved by
`MaterialNetwork`, not a separate conductivity-only query. The same snapshot
therefore carries electrical transport, thermal conductivity, physical condition
coordinates and original source receipts into the gas contact. Single-surface,
multi-surface and existing multi-chamber users of `resolve_contact` share this
handoff. The material-layer owner still checks the exact card, phase conditions,
quantity descriptors, conductivity key, source pin and isolation behavior.

## Reuse source selections, not frozen coefficients

`MaterialNetwork::prepare` selects the legacy resistivity basis.
`MaterialNetwork::prepare_with_bases` accepts one explicit conductivity/resistivity
choice per original electrical branch. Both return a borrowed
`PreparedMaterialNetwork`. It retains only the original network and its selected
transport claim IDs. It does not retain temperatures, physical histories or
cached coefficient values.

```rust,ignore
let selected = material_network.prepare_with_bases(
    &accepted_thermal, &electrical_bases, || cancelled(),
)?;

// In an existing coupled solver's material callback:
let mut coefficients = selected.resolve(trial_thermal, || cancelled())?;
// Compose original radiation/contact/source adapters here, preserving their
// separate source results. Do not add circuit heat twice.
```

Every resolve uses the same existing network evaluator on the actual trial
enthalpies. A selected source that leaves its supported domain refuses; a second
in-domain claim is not silently substituted. Prepared bindings may be retained
across callbacks or accepted steps when that source continuity is required.
Preparing again explicitly permits selection at a new starting state. This is
not a new material registry, property law, numerical solver or thermal history.

The existing standalone circuit/enthalpy step now uses this same preparation
path. Single- and multi-surface gas-network conveniences also prepare once at
the beginning of each proposed interval. Their legacy entry points retain the
resistivity-only choice; new entry points expose explicit electrical bases:

- `GasPiston::material_network_candidate_with_bases`
- `GasPiston::material_network_surfaces_candidate_with_bases`

Each takes `&[(original_branch_index, ElectricalPropertyBasis)]` immediately
after the external electrical ports. All remaining input geometry, nonlinear
controls, cell heat destinations and candidate-publication semantics are the
same as the corresponding legacy path. No fallback from missing resistivity to
conductivity (or conversely) is inferred.

These changes do not alter heat-transfer equations, circuit integration, gas
thermodynamics, piston dynamics, or their fixed-geometry/uniform-state limits.
The contact-aware electrical driver and custom multi-chamber callbacks are not
automatically replaced by these conveniences; they can explicitly reuse the
prepared source evaluator through their original composition APIs.

Focused checks:

```sh
cargo test -p fs-couple --test circuit_material_network
cargo test -p fs-couple --lib cavity_phs::piston::thermal::boundary::
```

New regressions cover fresh values with retained claims, refusal to switch either
transport source, preparation cancellation/retry, mixed electrical bases, and
complete snapshots through real single-/multi-surface gas trajectories. The
fixtures are synthetic. Native Rust compilation and tests were not available
during authoring; source review is not a native test pass.
