//! Held original-source bath traces for the native molecular face operator.

use super::*;

/// One explicit bank selection shared by velocity gradients and face physics.
#[derive(Clone, Copy)]
pub(super) struct View<'a> {
    pub(super) reservoirs: &'a Reservoirs,
    pub(super) trial: bool,
}

impl<'a> From<&'a Reservoirs> for View<'a> {
    fn from(reservoirs: &'a Reservoirs) -> Self {
        Self {
            reservoirs,
            trial: false,
        }
    }
}

impl<'a> View<'a> {
    pub(super) fn primitive_face(&self, side: BoundarySide, face: usize) -> Option<Primitive> {
        self.reservoirs.primitive_face_at(side, face, self.trial)
    }

    pub(super) fn normal_face(
        self,
        axis: usize,
        sign: f64,
        face: usize,
        species: usize,
    ) -> Option<SpatialPrimitive<'a>> {
        self.reservoirs
            .normal_face_at(axis, sign, face, species, self.trial)
    }
}

pub(super) fn admit(
    boundaries: Boundaries3d,
    reservoirs: &Reservoirs,
    pressure_outlets: PressureOutletPolicy,
) -> Result<(), DissipationError> {
    if reservoirs.has_stagnation() {
        return Err(DissipationError::Invalid(
            "compact stagnation reservoirs do not specify a molecular face bath",
        ));
    }
    for side in BoundarySide::ALL {
        match physical_condition(boundaries, side) {
            None | Some(BoundaryCondition::Reflecting) => {}
            Some(BoundaryCondition::Transmissive)
                if reservoirs.primitive_face(side, 0).is_some() => {}
            Some(
                BoundaryCondition::PressureOutlet { .. }
                | BoundaryCondition::CharacteristicOutlet { .. },
            )
                if pressure_outlets == PressureOutletPolicy::ZeroDiffusiveFlux => {}
            _ => {
                return Err(DissipationError::Invalid(
                    "open molecular sides require complete source-owned reservoir traces",
                ));
            }
        }
    }
    Ok(())
}

/// Same tangential ordering as CartesianGrid::face_center and boundary loads.
pub(super) fn face_index(grid: &CartesianGrid, cell: usize, axis: usize) -> usize {
    let c = grid.coordinates(cell);
    let [a, b] = tangents(axis);
    c[a] + grid.shape[a] * c[b]
}

#[allow(clippy::too_many_arguments)]
pub(super) fn velocity_neighbour(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: View<'_>,
    primitive: &[Primitive],
    cell: usize,
    axis: usize,
    upper: bool,
    walls: &[MolecularWall3d; 6],
) -> Result<[f64; 3], DissipationError> {
    let coordinate = grid.coordinates(cell)[axis];
    let crossing = if upper {
        coordinate + 1 == grid.shape[axis]
    } else {
        coordinate == 0
    };
    let side = side_for(axis, upper);
    if crossing && physical_condition(boundaries, side) == Some(BoundaryCondition::Transmissive) {
        let exterior = reservoirs
            .primitive_face(side, face_index(grid, cell, axis))
            .ok_or(DissipationError::Invalid(
                "unbound molecular reservoir gradient",
            ))?;
        // The trace lives on the physical face, half a cell away. Mirroring
        // about it preserves the existing centered derivative on cell centers.
        let ghost =
            std::array::from_fn(|d| 2.0 * exterior.velocity[d] - primitive[cell].velocity[d]);
        if !ghost.iter().all(|x| x.is_finite()) {
            return Err(DissipationError::Invalid(
                "molecular reservoir ghost velocity overflow",
            ));
        }
        return Ok(ghost);
    }
    let (neighbour, reflected) = cartesian::neighbour(grid, boundaries, cell, axis, upper);
    wall_velocity_sample(
        primitive,
        neighbour,
        cell,
        axis,
        upper,
        reflected,
        walls,
        0.5 * grid.spacing[axis],
    )
}

/// Midpoint tangential derivative between the actual cell and exterior traces.
/// A uniform bath has zero tangential derivatives. Face profiles use centered
/// differences, periodic wrapping, or quadratic one-sided end differences;
/// two faces use a linear difference and one face resolves no variation.
pub(super) fn tangential_gradient(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: View<'_>,
    side: BoundarySide,
    cell: usize,
    interior_gradient: &[[f64; 3]; 3],
) -> Result<[[f64; 3]; 2], DissipationError> {
    let axis = side.axis();
    let [a, b] = tangents(axis);
    let c = grid.coordinates(cell);
    let mut result = [[0.0; 3]; 2];
    for (row, derivative) in [a, b].into_iter().enumerate() {
        let n = grid.shape[derivative];
        let i = c[derivative];
        let value = |index: usize| -> Result<[f64; 3], DissipationError> {
            let mut p = c;
            p[derivative] = index;
            reservoirs
                .primitive_face(side, p[a] + grid.shape[a] * p[b])
                .map(|q| q.velocity)
                .ok_or(DissipationError::Invalid(
                    "unbound tangential reservoir trace",
                ))
        };
        let boundary_gradient = if n == 1 {
            [0.0; 3]
        } else if boundaries.axis(derivative) == AxisBoundaries::Periodic {
            let low = value(if i == 0 { n - 1 } else { i - 1 })?;
            let high = value(if i + 1 == n { 0 } else { i + 1 })?;
            std::array::from_fn(|d| (high[d] - low[d]) / (2.0 * grid.spacing[derivative]))
        } else if n == 2 {
            let low = value(0)?;
            let high = value(1)?;
            std::array::from_fn(|d| (high[d] - low[d]) / grid.spacing[derivative])
        } else if i == 0 {
            let at = value(0)?;
            let next = value(1)?;
            let far = value(2)?;
            std::array::from_fn(|d| {
                (4.0 * (next[d] - at[d]) - (far[d] - at[d])) / (2.0 * grid.spacing[derivative])
            })
        } else if i + 1 == n {
            let at = value(i)?;
            let next = value(i - 1)?;
            let far = value(i - 2)?;
            std::array::from_fn(|d| {
                (4.0 * (at[d] - next[d]) - (at[d] - far[d])) / (2.0 * grid.spacing[derivative])
            })
        } else {
            let low = value(i - 1)?;
            let high = value(i + 1)?;
            std::array::from_fn(|d| (high[d] - low[d]) / (2.0 * grid.spacing[derivative]))
        };
        for (column, component) in [axis, a, b].into_iter().enumerate() {
            result[row][column] =
                0.5 * interior_gradient[derivative][component] + 0.5 * boundary_gradient[component];
        }
    }
    if !result.iter().flatten().all(|x| x.is_finite()) {
        return Err(DissipationError::Invalid(
            "tangential reservoir gradient overflow",
        ));
    }
    Ok(result)
}
