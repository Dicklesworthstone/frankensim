## Distributed dielectric loss and finite thermal state

The default-off `resolved-dielectric-thermal` feature exposes
`ConductionThermalNetwork::step_dielectric_field_candidate` and the same method
on `PreparedHeatNetwork`. Both consume the physical interval of a borrowed
native `DielectricFieldCandidate<F,D>`, then jointly return its field endpoint
and the advanced `DepositionSimulation`. Planar, axisymmetric and tetrahedral
fields use one thermal residual and their original physical cell volumes.
The planar `step_debye_candidate` convenience methods also prepare the field
candidate from its original field and reusable operator. Every electrical cell
has one explicitly mapped thermal receiver; the mapped electrical cells together
must cover each receiver's complete physical volume. `ConductionRegion::triangle` indexes the original
electrical cell, including tetrahedra on the 3-D path.
This opt-in also enables `fs-em/resolved-leaky-dielectric` for explicitly
source-backed Ohmic leakage; the baseline Debye field remains independently usable.

The native field step supplies local irreversible loss. The coupling deposits
that loss once through the existing backward-Euler enthalpy residual; neither
instantaneous electric storage nor polarization storage is heat. The retained
field receipt preserves separate polarization and Ohmic losses. The joint
physical balance is `delta U_field + delta H_thermal = W_electrical + Q_external`.
The final balance sums the native displacement and conduction work separately
with the existing compensated thermal-account arithmetic.
Closed-enclosure radiation remains internal exchange, with its computed net
reported separately and removed from external supply. Fourier links, sourced
thermal interfaces, integrated conduction layers and radiation use their
existing material and thermal owners.

The bridge admits fully solid fixed-volume charts and validity-wide constant
electrical storage, relaxation and DC conductivity coefficients. An explicitly
sourced zero conductivity retains the same source-validation obligations.
Conductivity lookup uses its original property name and conductivity dimensions.
All original source claims are
revalidated at every trial temperature, and the accepted endpoint snapshots are
retained alongside the unchanged reference snapshots in the field. Invalid material ownership, changed
external nonthermal context, invalid volume mapping, mismatched clocks,
expired source validity or unowned graph electrical work refuses. No phase polarization, dielectric
parameter work, electrocaloric derivative or moving geometry is inferred.

One already prepared native field candidate serves the entire thermal solve,
because these electrical coefficients remain fixed. Native midpoint dielectric
work/heat and backward-Euler thermal transport retain their separate
discretizations; the coupled transport is first order in general. Failed energy
checks, source queries or cancellation publish neither field nor thermal state.
The shared report's legacy `deposited_joule_j` name means all deposited physical
dielectric loss on this lane.

Focused regressions cover analytical polarization heating over repeated steps,
heterogeneous series dielectric layers with independent finite-cell Fourier
transport, separate DC leakage and polarization heat, revolved and tetrahedral
physical volumes, closed radiation exchange, the lossless storage limit,
source-domain refusal and deterministic retry after late cancellation. Synthetic source
fixtures establish numerical identities, not measured material accuracy.

Targeted validation: `cargo test -p fs-couple --features resolved-dielectric-thermal --test dielectric_field_thermal`.
