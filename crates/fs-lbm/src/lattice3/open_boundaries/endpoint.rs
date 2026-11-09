//! Density- and velocity-dependent sources in the original open stream pass.

use super::*;
use crate::d3q19::DarcyForchheimer3;
use crate::lattice3::sources::resistance::admit_drag;

/// The held field is immutable while neighboring boundary donors still read it.
#[derive(Clone, Copy)]
pub(super) enum EndpointSource3<'a> {
    Acceleration(&'a [[f64; 3]]),
    Resistance(&'a [DarcyForchheimer3]),
}

/// Each worker resolves its own force epoch. Constitutive drag is retained
/// separately so a much larger independent body force cannot erase its receipt.
pub(super) struct EndpointWorkspace3<'a> {
    pub(super) forces: &'a mut [[f64; 3]],
    pub(super) drag: Option<&'a mut [[f64; 3]]>,
}

impl EndpointWorkspace3<'_> {
    pub(super) fn part(&mut self, range: core::ops::Range<usize>) -> EndpointWorkspace3<'_> {
        EndpointWorkspace3 {
            forces: &mut self.forces[range.clone()],
            drag: self.drag.as_deref_mut().map(|field| &mut field[range]),
        }
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    pub(super) fn admit_endpoint<'a, 'b>(
        &self,
        field: Option<(EndpointSource3<'a>, EndpointWorkspace3<'b>)>,
        collisions: Option<&[L::Collision]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(Option<EndpointSource3<'a>>, Option<EndpointWorkspace3<'b>>), GridError3<L::Error>>
    {
        let Some((source, resolved)) = field else {
            return Ok((None, None));
        };
        let cells = self.grid.f.len();
        if resolved.forces.len() != cells
            || resolved.drag.as_ref().is_some_and(|drag| drag.len() != cells)
        {
            return Err(GridError3::Input("endpoint force workspace length"));
        }
        match source {
            EndpointSource3::Acceleration(acceleration) => {
                crate::lattice3::sources::admit_forces(
                    &self.grid, acceleration, collisions, poll,
                )?;
            }
            EndpointSource3::Resistance(laws) => {
                // The resistance owner admits accepted or next-endpoint laws
                // and seam topology before staging; no callback runs here.
                if laws.len() != cells || resolved.drag.is_none() {
                    return Err(GridError3::Input("resistance endpoint workspace length"));
                }
            }
        }
        Ok((Some(source), Some(resolved)))
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenStep<'_, L, Q> {
    // Donor and final observations use their OWN density, raw momentum and law.
    // Only the original native implicit constitutive solve resolves resistance.
    pub(super) fn endpoint_state(
        &self,
        populations: &[f64; Q],
        cell: usize,
    ) -> Result<(f64, [f64; 3], [f64; 3], [f64; 3]), GridError3<L::Error>> {
        match self.endpoint {
            Some(EndpointSource3::Acceleration(field)) => {
                let (rho, momentum) = raw_moments_cell::<L, Q>(populations, cell)?;
                let force = field[cell].map(|a| rho * a);
                let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
                Ok((rho, velocity, force, [0.0; 3]))
            }
            Some(EndpointSource3::Resistance(laws)) => {
                let (rho, momentum) = raw_moments_cell::<L, Q>(populations, cell)?;
                let (force, drag) = laws[cell]
                    .endpoint_force(rho, momentum)
                    .ok_or(GridError3::NonFiniteCell { cell })?;
                let collision = self.data.collisions.map_or(self.data.collision, |field| field[cell]);
                if force != [0.0; 3] && !L::supports_body_force(collision) {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "collision model does not admit body force",
                    });
                }
                let (_, velocity) =
                    macro_cell::<L, Q>(populations, force, cell, self.data.budget.maximum_mach)?;
                admit_drag::<L, Q>(laws[cell], velocity, drag, cell)?;
                Ok((rho, velocity, force, drag))
            }
            None => {
                let force = self.data.forces.map_or(self.data.force, |field| field[cell]);
                let (rho, velocity) =
                    macro_cell::<L, Q>(populations, force, cell, self.data.budget.maximum_mach)?;
                Ok((rho, velocity, force, [0.0; 3]))
            }
        }
    }
}
