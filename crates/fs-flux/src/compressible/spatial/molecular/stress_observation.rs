//! Borrow the original accepted-state stress for acoustic source extraction.
use super::*;

enum FaceStorage<'a> {
    Frozen { flux: &'a [Vec<SpatialFaceFlux>; 3], states: &'a [State] },
    Mixture { flux: &'a [Vec<f64>; 3], states: &'a [f64], species: usize },
}

/// Accepted-state primitives and the ORIGINAL conservative viscous face stress.
///
/// Positive traction is `tau[:,axis]`, while the original momentum flux is
/// `-tau[:,axis]`. This includes explicit bulk viscosity, actual material laws,
/// configured walls/reservoirs and any selected SGS closure. It contains no
/// inviscid pressure or advective momentum flux. Heat/species diffusion stays
/// in its native owner and is not relabeled as an acoustic force.
///
/// This view borrows reusable operator storage after a complete successful
/// evaluation. It allocates nothing and cannot outlive or advance its owner.
/// It describes the held accepted endpoint, not a time-integrated step load.
pub struct StressField3d<'a> {
    grid: &'a CartesianGrid,
    boundaries: Boundaries3d,
    time_s: f64,
    primitive: &'a [Primitive],
    faces: FaceStorage<'a>,
}

impl<'a> StressField3d<'a> {
    pub(in crate::compressible::spatial) fn frozen(
        grid: &'a CartesianGrid,
        boundaries: Boundaries3d,
        time_s: f64,
        primitive: &'a [Primitive],
        flux: &'a [Vec<SpatialFaceFlux>; 3],
        states: &'a [State],
    ) -> Self {
        Self { grid, boundaries, time_s, primitive, faces: FaceStorage::Frozen { flux, states } }
    }

    pub(in crate::compressible::spatial) fn mixture(
        grid: &'a CartesianGrid,
        boundaries: Boundaries3d,
        time_s: f64,
        primitive: &'a [Primitive],
        flux: &'a [Vec<f64>; 3],
        states: &'a [f64],
        species: usize,
    ) -> Self {
        Self { grid, boundaries, time_s, primitive, faces: FaceStorage::Mixture { flux, states, species } }
    }

    pub const fn grid(&self) -> &CartesianGrid { self.grid }
    pub const fn boundaries(&self) -> Boundaries3d { self.boundaries }
    pub const fn time_s(&self) -> f64 { self.time_s }

    /// Original source-decoded primitive in native x-fast cell order.
    pub fn primitive(&self, cell: usize) -> Option<Primitive> {
        self.primitive.get(cell).copied()
    }

    /// Original conserved momentum density [kg/(m² s)], never reconstructed
    /// by multiplying the rounded primitive velocity and density.
    pub fn momentum(&self, cell: usize) -> Option<[f64; 3]> {
        self.primitive.get(cell)?;
        Some(match self.faces {
            FaceStorage::Frozen { states, .. } => states[cell].momentum,
            FaceStorage::Mixture { states, species, .. } => {
                let start = cell * (species + 4) + species;
                [states[start], states[start + 1], states[start + 2]]
            }
        })
    }

    /// Actual face traction [Pa], in the POSITIVE native coordinate direction.
    ///
    /// `coordinates[axis]` is a face coordinate in `0..=shape[axis]`; both
    /// other coordinates are cell coordinates in `0..shape[tangent]`.
    /// Multiply by the outward normal sign before a closed-surface integral.
    /// This reads the face directly; no averaging, gradient or constitutive
    /// evaluation is performed. Invalid axes or coordinates return `None`.
    pub fn face_traction_pa(&self, axis: usize, coordinates: [usize; 3]) -> Option<[f64; 3]> {
        if axis >= 3 || (0..3).any(|a| coordinates[a] >= self.grid.shape[a] + usize::from(a == axis)) {
            return None;
        }
        let face = self.grid.cell_face(axis, coordinates);
        Some(match self.faces {
            FaceStorage::Frozen { flux, .. } => {
                let f = super::super::evolve::flux_array(flux[axis][face], axis);
                [-f[1], -f[2], -f[3]]
            }
            FaceStorage::Mixture { flux, species, .. } => {
                let start = face * (species + 4) + species;
                [-flux[axis][start], -flux[axis][start + 1], -flux[axis][start + 2]]
            }
        })
    }

    /// Cell-center symmetric stress [Pa], ordered xx, yy, zz, xy, xz, yz.
    ///
    /// Opposing original faces are averaged to the cell center; off-diagonal
    /// pairs are then symmetrized. This is a spatial reconstruction, exact for
    /// affine stress and second order for smooth resolved interior fields.
    /// Walls, coefficient jumps and modeled SGS retain their original limits.
    /// Surface sources should read `face_traction_pa` directly instead.
    pub fn cell_stress_pa(&self, cell: usize) -> Option<[f64; 6]> {
        self.primitive.get(cell)?;
        let coordinates = self.grid.coordinates(cell);
        let mut tensor = [[0.0; 3]; 3];
        for axis in 0..3 {
            let mut upper = coordinates;
            upper[axis] += 1;
            let low = self.face_traction_pa(axis, coordinates)?;
            let high = self.face_traction_pa(axis, upper)?;
            for component in 0..3 {
                tensor[component][axis] = 0.5 * low[component] + 0.5 * high[component];
            }
        }
        Some([
            tensor[0][0], tensor[1][1], tensor[2][2],
            0.5 * tensor[0][1] + 0.5 * tensor[1][0],
            0.5 * tensor[0][2] + 0.5 * tensor[2][0],
            0.5 * tensor[1][2] + 0.5 * tensor[2][1],
        ])
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    pub(super) fn prepare_stress_field(
        &mut self, flow: &Euler3d, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        poll(cancelled)?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid("stress observer and accepted flow differ"));
        }
        spatial(
            &self.gas, &self.grid, self.boundaries, &self.walls,
            self.reservoir_exchange.then(|| flow.reservoirs.accepted()),
            &flow.states, &mut self.law, &mut self.work.primitive,
            &mut self.work.gradient, &mut self.work.viscosity,
            &mut self.work.flux, &mut self.work.binary,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            cancelled,
        )?;
        poll(cancelled)
    }
}
