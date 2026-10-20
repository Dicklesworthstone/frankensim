# Layer-resolved self-heating

`fs_couple::strip_enthalpy::step_candidate` couples the existing `LayeredStrip`
material histories to the shared `fs_thermal` enthalpy network. Each mechanical
layer has one explicit thermal destination; different layers can have different
materials, heat capacities and temperatures. Passive thermal cells can exchange
heat with those layers without acquiring invented mechanical behavior.

## Shared physics and accounting

Every nonlinear trial starts from the same accepted material-point histories and
cell enthalpies. The original compatible section solve supplies endpoint strain,
curvature and each layer's actual physical dashpot dissipation. That dissipation
enters the corresponding thermal cell exactly once. The existing backward-Euler
network residual and Newton-Krylov solver determine the coupled endpoint; no
second heat-capacity law, Maxwell update or compatible section solver is used.

Numerical dissipation is NOT heat. The final fresh evaluation checks local cell
residuals, total thermal balance, structural equilibrium and the combined account:

```
elastic storage change + thermal storage change + numerical loss
    - imposed mechanical work - independently evaluated other heat
```

Prescribed force/moment controls supply imposed work. Prescribed strain/curvature
uses actual reactions. Other heat is measured with the original thermal rate
owner before adding material dissipation; it is not inferred by forcing a balance.
The returned thermal audit already includes physical dissipation in its external
heat field. Do not add that heat a second time.

Bindings require the actual retained mechanical and caloric objects, compatible
source cards and conditions, matching reference volume and mass, and matching
accepted temperature. A material name or matching scalar coefficient is not a
substitute for these bindings. Endpoint callbacks can supply independent heat
sources and temperature-dependent replacements for existing link conductances;
they must not mutate external physical history.

## Timeline

Process start-time events before a coupled interval. Energy-increment guesses are
in sorted thermal CellId order, not layer order. An interior event requires an
interval cut; target-time events remain pending until the coupled endpoint is
accepted. `into_parts()` returns both histories together. Use the existing thermal
birth operation and `LayeredStrip::deposit` to add stress-free material at that
accepted time. Rejected trials never age or add material.

## Focused verification

Seven integration regressions exercise exact heat partition, the actual two-cell
transport solve, a heat-capacity-only substitution that changes compatible
curvature, numerical-loss separation, retained load/unload history, mapping and
budget refusals, cancellation, and a subsequent material birth. They reuse the
original material-card, structural, thermal and nonlinear implementations.

```
rch exec -- cargo test -p fs-solid --test deposition_layer_energy
rch exec -- cargo test -p fs-couple --test strip_enthalpy
```

Native Rust execution is not established by this change. RCH, Cargo and rustc
were absent in the editing environment. The independent Python/SciPy reference
script was rerun successfully for 324 parameter cases, with maximum combined
energy discrepancy 1.48e-15 J. Those are equation checks, not Rust execution or
experimental material validation. The synthetic low heat capacities in the
fixtures are deliberate sensitivity inputs, not properties of a named material.

## Physical scope

Small-strain, axially uniform, perfectly bonded Euler-Bernoulli sections; one
uniform temperature per layer; fixed thermal geometry; single solid phase;
constant density and temperature-independent spring energy and free thermal
strain. Relaxation times and heat capacity may vary. Temperature-dependent
elastic free energy, expansion work, melting, yielding, inertia, moving contact
and full three-dimensional printing are not supplied by this coupling. Source
agreement does not establish empirical validity, and nonlinear convergence does
not establish time accuracy. Refine timestep and spatial representation separately.
