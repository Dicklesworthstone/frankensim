# Material contact refresh and between-impact thermal feedback

The source-bound contact runtime can now refresh an unloaded contact without
rebuilding the simulation or resetting motion, mechanical memory, heat or time.
The optional thermal-feedback path joins the existing owners in this order:

```
accepted reservoir enthalpy/temperature
  -> original material and interface queries
  -> unloaded contact source refresh
  -> existing coupled mechanical step
  -> actual contact work and existing thermal transport
  -> joint acceptance
```

## Generic source and mechanical APIs

`fs_contact::normal_patch::source::state::PointContactMaterialCards` retains the
two original bulk cards, the ordered interface card and explicit property/model
selections. `resolve(point, geometry, envelope)` calls the existing bulk,
interface, normal-model and Hertz-geometry owners. It does not infer properties
from material names, replace missing data with a typical material, or duplicate
a constitutive formula. Every returned `SourcePointContact` owns its original
source records and admitted operating envelope.

`ImpactSystem::refresh_unloaded_contacts` and the prepared equivalent accept
`ContactSourceUpdate` values by existing contact key. Every changed contact must
have **exactly zero indentation**. Geometry, gap, coordinate map, body inertia,
state, sample counter and histories remain unchanged. A batch is admitted before
anything changes. Duplicate or unknown keys, changed geometry, a source-rate
excursion, compression or cancellation refuse the complete batch.

A shared live slot holds the original `fs-dcontact` owner underneath the existing
storage wrappers. This matters when acoustic or relaxation storage has already
been attached: replacing a separate contact field would leave its Hamiltonian
using old stiffness while forces and tangents used new stiffness. The slot keeps
all consumers on the same owner without nesting another wrapper on each refresh.

Generic refresh changes contact constitutive data only. It is not whole-body
material replacement: changing density, geometry or inertia requires those
owners to supply a consistent mechanical reconfiguration.

## Thermal feedback with joint rollback

`HeatedImpactSystem::step_with_material_feedback` takes a material-resolution
callback and the same thermal-network callback as `step_with_coupled_transport`.
The first callback sees an unloaded contact and all accepted finite thermal
reservoirs. It returns a new source, or `None` when unchanged. Loaded contacts
hold their existing source. The result records `refreshed_contacts` and
`held_loaded_contacts`, alongside the unchanged thermal/mechanical reports.

A later thermal, source, solver or cancellation refusal also restores the
previous contact source and loss-path selection. It does not leave refreshed
constitutive data behind a rolled-back motion/heat step. Callbacks must propose
rather than mutate outside physical state. Arbitrary callback side effects or
panics are not transactional.

`resolve_isothermal_contact` supplies a stricter convenience binding. It requires
the exact original ordered bulk/interface cards, fully solid thermal states,
the matching caloric-card identities, unchanged nonthermal coordinates and axis
schemas, and consistent unchanged bulk/caloric densities. Both thermal stores
must have the same temperature, because the existing isotropic contact binder
requires a common physical state. Unequal temperatures are rejected, never
averaged into an invented interface state. General callers may supply another
physically justified source resolver under its own declared assumptions.

```rust,ignore
use fs_couple::render::plate::impact::material::heating::feedback::resolve_isothermal_contact;

// The same cards and source selection used for this contact's initial binding.
// thermal_network returns the canonical joint thermal proposals, not live mutations.
let frame = heated.step_with_material_feedback(
    &forces,
    &gate,
    |contact, stores| {
        resolve_isothermal_contact(
            contact, &cards, &reference_point, [&stores[0], &stores[1]],
        ).map(Some)
    },
    thermal_network,
)?;
```

A repeated-impact regression changes only synthetic heat-capacity data. The
first impact has the same mechanical response. Its heat changes the accepted
material temperature; the refreshed source then changes the next impact. The
same original mechanical and thermal solvers run in both cases.

## Exact boundary of this model

This is **sample-and-hold feedback between contact episodes**, not continuous
loaded thermoelastic coupling. A compressed contact stores energy; changing its
stiffness without constitutive work accounting would manufacture or remove
energy. The implementation refuses direct loaded refresh and holds that source
through a compression episode. The caller must establish that within-episode
thermal changes and the frozen material envelope are acceptable for its problem.
There is no automatic accuracy estimate or guarantee for this approximation.

Thermal expansion, density-changing geometry, contact plasticity, unequal-body
temperature laws, evolving adhesion and moving melt are not supplied here. The
existing reference inertia and thermal geometry are retained. Heat partition,
thermal lumping and contact geometry remain explicit physical inputs. No
experimental accuracy, hard-real-time performance or full-project material
substitutability is claimed.

Focused native checks, including the existing contact and heat regressions:

```
cargo test -p fs-contact --test source_point -- --test-threads=1
cargo test -p fs-couple --test material_contact_heat -- --test-threads=1
cargo test -p fs-couple --test material_impact -- --test-threads=1
```
