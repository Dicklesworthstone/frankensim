# Growing layered self-heating on one accepted history

`fs_couple::strip_enthalpy::growth::GrowingStrip` connects the existing coupled
section/heat solve to actual scheduled thermal births. This closes the gap where
an application could accept new thermal material and then fail to add its
mechanical layer, leaving the two histories inconsistent.

## Use

Construct the normal `DepositionSimulation` with its physical schedule and
material-backed caloric laws. Supply a `LayeredStrip`, the selected physical body
ID, complete bindings for active layers, and `LayerBirth` descriptors for future
cells. Each descriptor holds an existing `LayerThermalBinding` and reference
`LayerGeometry`. `GrowingStrip::try_new` accepts a maximum registered layer count
and an explicit event-energy tolerance. It returns the initialized owner and the
actual current-time event audit. Empty initial states are allowed when scheduled
current-time births create a nonempty selected-body seed.

Call `state.step_candidate(target_s, control, config, energy_guesses, resolve,
cancelled)` for each physical interval. The result contains a new complete owner
and a `GrowthReport`. Replace the accepted owner only after success. Both the
initialization and step APIs retain the existing source, caloric, reference-mass,
volume, material-history and temperature checks. There is one shared admission
implementation, not a second density resolver or a relaxed birth-only check.

The caller must split at scheduled event times and load/source discontinuities.
The original thermal solver rejects an event strictly inside an interval; this
owner does not silently subcycle or skip that event. Energy guesses remain in
sorted thermal CellId order, including passive thermal cells.

## Ordering and rollback

1. The original section/enthalpy coupling solves the interval from immutable
   accepted histories. No Newton trial deposits material or advances a schedule.
2. The thermal timeline applies all events at the accepted target in their
   original order, with zero thermal substeps. Actual incoming enthalpy determines
   each newborn temperature; the existing deposit operation initializes its
   stress-free reference on the current compatible strain and curvature.
3. The same material admission checks validate the complete resulting body. The
   full storage/work account includes incoming event enthalpy exactly once.
   Cancellation or failure in the last birth discards the entire candidate,
   including the solved interval and every earlier birth in that batch.

Every birth in the selected body requires an explicit mechanical binding. A
binding cannot switch to another physical body. Other thermal bodies may be
born and retain their incoming mass/enthalpy without acquiring invented
mechanical layers. This distinction supports finite thermal supports or other
passive cells on the same network.

`report.interval` describes the pre-event coupled endpoint. `report.events`
contains the actual event mass/energy account. Neither should be added twice.
`birth_equilibrium_pending` marks a newly born or initially supplied specimen
that this owner has not yet equilibrated in a positive-time coupled step. Birth
itself is not an instantaneous equilibrium, fracture, welding or impact solve.
The next interval uses the expanded layer set and its original retained history.

## Verification and physical boundary

Six focused regressions compare the composed operation with the original manual
owners, then evolve the new material in another coupled interval. They also
exercise a final invalid birth after an earlier valid birth, final cancellation,
checkpoint replay, initial scheduled activation, passive-body enthalpy, and
complete/bounded registrations.

```
rch exec -- cargo test -p fs-couple --test strip_enthalpy --test strip_enthalpy_growth
```

Native Rust tests were not executed in the editing environment: RCH, Cargo and
rustc are unavailable. The 324-case independent thermal/Maxwell reference script
was rerun for the underlying unchanged interval equations; it does not execute
this Rust event composition. Source review verified that the nonlinear solve
after admission is unchanged. No native or empirical validation is claimed.

The original small-strain, axially uniform, perfectly bonded section and
constant-density/spring-energy/free-strain caloric restrictions still apply.
This is not the spatial growing-road/bridge solver, a new event scheduler,
automatic timestep selection, moving-contact feedback, phase change or a complete
printing model. Geometry, material applicability, thermal uniformity and the
perfect-bond assumption require independent justification. Callbacks must not
mutate external physical history; arbitrary external side effects cannot be
rolled back by cloning this owner.
