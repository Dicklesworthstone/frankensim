//! Native two-dimensional BDM1 flow kinematics.
//!
//! H(div) conformity makes the normal velocity continuous, but allows jumps
//! in the tangential velocity. The distributional curl therefore contains
//! both the constant curl inside each triangle and an affine vortex sheet on
//! every interior edge. Omitting those sheets violates Stokes' theorem for a
//! general discrete solution. No finite-difference stencil or nearest-cell
//! sampling is used here; gradients and traces come from the original basis.
//!
//! Material-flow callers use the same native `FluxSystem` they bound at
//! construction, passing their accepted solution's `velocity_dofs()` directly.
//! Those coefficients are already in m/s: no pressure/density conversion is
//! applied to vorticity, and SI mesh geometry gives rates in 1/s.

use super::FluxSystem;
use crate::bdm::eval_basis;

/// Constant velocity gradient inside one triangle, in the mesh's coordinates.
/// With mesh coordinates in metres and velocity in m/s, all rates are in 1/s.
/// This is the broken (cell-interior) gradient; it excludes interface sheets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellKinematics2 {
    gradient: [[f64; 2]; 2],
}

impl CellKinematics2 {
    /// `gradient[i][j] = partial_j u_i` for the affine BDM1 velocity.
    #[must_use]
    pub const fn gradient(self) -> [[f64; 2]; 2] {
        self.gradient
    }

    /// Signed out-of-plane curl `partial_x u_y - partial_y u_x`.
    /// Positive values correspond to counterclockwise fluid rotation.
    /// The angular velocity of a rigidly rotating parcel is half this value.
    #[must_use]
    pub fn vorticity(self) -> f64 {
        self.gradient[1][0] - self.gradient[0][1]
    }

    /// Native cell divergence `partial_x u_x + partial_y u_y`.
    #[must_use]
    pub fn divergence(self) -> f64 {
        self.gradient[0][0] + self.gradient[1][1]
    }

    /// Objective symmetric strain invariant `sqrt(2 D:D)`.
    /// Rigid rotation has vorticity but contributes no viscous shear rate.
    #[must_use]
    pub fn shear_rate(self) -> f64 {
        (core::f64::consts::SQRT_2 * self.gradient[0][0])
            .hypot(core::f64::consts::SQRT_2 * self.gradient[1][1])
            .hypot(self.gradient[0][1] + self.gradient[1][0])
    }
}

/// Affine distributional vorticity on one interior edge.
///
/// Let `t = [-normal_y, normal_x]`, with the native normal pointing out of
/// `edge.tris.0` into `edge.tris.1`. The sheet strength is `(u_1 - u_0).t`,
/// evaluated at the endpoints in `edge.verts` order. Multiplying it by the
/// line delta distribution gives curl. Its SI unit is m/s, and its integral
/// is circulation [m²/s]; it is not a bulk vorticity in 1/s. In particular,
/// squaring a sheet strength is not a mesh-independent enstrophy density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VorticitySheet2 {
    /// Native interior-edge index.
    pub edge: usize,
    /// Signed strengths at the two globally ordered edge vertices.
    pub endpoint_strengths: [f64; 2],
    /// Native edge length.
    pub length: f64,
}

impl VorticitySheet2 {
    /// Exact signed integral of the affine sheet strength.
    #[must_use]
    pub fn circulation(self) -> f64 {
        self.length * (0.5 * self.endpoint_strengths[0] + 0.5 * self.endpoint_strengths[1])
    }
}

/// Stokes accounting for the represented BDM field, in velocity × length.
/// Outer boundaries are counterclockwise and hole boundaries clockwise,
/// following the fluid cells' outward normals. Disconnected domains add.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CirculationBalance2 {
    /// Integral of the constant cell-interior curl over all triangles.
    pub cell_circulation: f64,
    /// Integral of all interior tangential-jump sheets.
    pub sheet_circulation: f64,
    /// Integral of the native one-sided velocity along the fluid boundary.
    /// Weakly imposed wall velocities need not equal this solved trace.
    pub boundary_circulation: f64,
}

impl CirculationBalance2 {
    /// Integrated distributional curl, including interior sheets.
    #[must_use]
    pub fn total_circulation(self) -> f64 {
        self.cell_circulation + self.sheet_circulation
    }

    /// Signed Stokes defect; zero up to floating-point reconstruction error.
    /// This identity checks accounting, not spatial convergence to a flow.
    #[must_use]
    pub fn stokes_defect(self) -> f64 {
        self.total_circulation() - self.boundary_circulation
    }
}

/// An incomplete, nonfinite or out-of-range kinematics query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FluxKinematicsError {
    /// Input length, edge index or scalar test-function values are invalid.
    Invalid(&'static str),
    /// The local gradient, curl or divergence cannot be represented finitely.
    NonFiniteCell(usize),
    /// A local trace or sheet strength cannot be represented finitely.
    NonFiniteEdge(usize),
    /// The integrated result cannot be represented finitely.
    NonFiniteIntegral,
}

impl core::fmt::Display for FluxKinematicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "flow kinematics: {reason}"),
            Self::NonFiniteCell(cell) => write!(f, "flow kinematics: nonfinite cell {cell}"),
            Self::NonFiniteEdge(edge) => write!(f, "flow kinematics: nonfinite edge {edge}"),
            Self::NonFiniteIntegral => write!(f, "flow kinematics: nonfinite integral"),
        }
    }
}

impl std::error::Error for FluxKinematicsError {}

impl FluxSystem<'_> {
    /// Allocation-free gradient/curl/divergence in one native triangle.
    /// Only its six velocity coefficients are read; trailing pressure entries
    /// are ignored. Returns `None` for a missing cell, an incomplete velocity
    /// vector or a nonfinite local rate. No mesh search or stencil is needed.
    #[must_use]
    pub fn try_cell_kinematics(&self, velocity: &[f64], cell: usize) -> Option<CellKinematics2> {
        if velocity.len() < self.n_u {
            return None;
        }
        let basis = self.bases.get(cell)?;
        let mut gradient = [[0.0; 2]; 2];
        for i in 0..6 {
            let coefficient = velocity[self.dof(cell, i)];
            for (row, basis_row) in gradient.iter_mut().zip(basis.grad[i]) {
                for (value, derivative) in row.iter_mut().zip(basis_row) {
                    *value += coefficient * derivative;
                }
            }
        }
        let result = CellKinematics2 { gradient };
        (gradient.iter().flatten().all(|value| value.is_finite())
            && result.vorticity().is_finite()
            && result.divergence().is_finite()
            && result.shear_rate().is_finite())
        .then_some(result)
    }

    /// Gradients and rates in native triangle order, using cached BDM bases.
    ///
    /// # Errors
    /// Refuses an incomplete velocity vector or any nonfinite cell rate.
    pub fn cell_kinematics(
        &self,
        velocity: &[f64],
    ) -> Result<Vec<CellKinematics2>, FluxKinematicsError> {
        self.check_kinematics_length(velocity)?;
        (0..self.bases.len())
            .map(|cell| {
                self.try_cell_kinematics(velocity, cell)
                    .ok_or(FluxKinematicsError::NonFiniteCell(cell))
            })
            .collect()
    }

    /// Exact affine vortex sheet on a native interior edge; no allocation.
    /// Boundary edges return `Ok(None)` because this diagnostic does not
    /// extend the velocity by zero outside the fluid or invent a wall sheet.
    ///
    /// # Errors
    /// Refuses incomplete velocity data, an absent edge or nonfinite traces.
    pub fn vorticity_sheet(
        &self,
        velocity: &[f64],
        edge: usize,
    ) -> Result<Option<VorticitySheet2>, FluxKinematicsError> {
        self.check_kinematics_length(velocity)?;
        let face = self
            .mesh
            .edges
            .get(edge)
            .ok_or(FluxKinematicsError::Invalid("edge index out of range"))?;
        if face.tris.1 == usize::MAX {
            return Ok(None);
        }
        let mut endpoint_strengths = [0.0; 2];
        for (strength, vertex) in endpoint_strengths
            .iter_mut()
            .zip([face.verts.0, face.verts.1])
        {
            let point = self.mesh.verts[vertex];
            *strength = self.kinematics_sheet_strength(velocity, edge, point);
        }
        if endpoint_strengths.iter().any(|value| !value.is_finite()) {
            return Err(FluxKinematicsError::NonFiniteEdge(edge));
        }
        Ok(Some(VorticitySheet2 {
            edge,
            endpoint_strengths,
            length: face.len,
        }))
    }

    /// Exact Stokes balance for the discrete field, including tangential jumps.
    /// Work is linear in cells and edges, deterministic and allocation-free.
    /// Midpoint integration is exact for affine BDM boundary and jump traces;
    /// only one point per edge is needed for total circulation.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite velocity data or nonfinite integrals.
    pub fn vorticity_balance(
        &self,
        velocity: &[f64],
    ) -> Result<CirculationBalance2, FluxKinematicsError> {
        self.check_kinematics_length(velocity)?;
        let mut cells = Sum::default();
        let mut sheets = Sum::default();
        let mut boundary = Sum::default();
        for (cell, &area) in self.mesh.areas.iter().enumerate() {
            let flow = self
                .try_cell_kinematics(velocity, cell)
                .ok_or(FluxKinematicsError::NonFiniteCell(cell))?;
            cells.add(area * flow.vorticity());
        }
        for (edge, face) in self.mesh.edges.iter().enumerate() {
            if face.tris.1 != usize::MAX {
                let circulation =
                    face.len * self.kinematics_sheet_strength(velocity, edge, face.mid);
                if !circulation.is_finite() {
                    return Err(FluxKinematicsError::NonFiniteEdge(edge));
                }
                sheets.add(circulation);
            } else {
                let trace = self.kinematics_velocity(velocity, face.tris.0, face.mid);
                let circulation =
                    face.len * (trace[1] * face.normal[0] - trace[0] * face.normal[1]);
                if !circulation.is_finite() {
                    return Err(FluxKinematicsError::NonFiniteEdge(edge));
                }
                boundary.add(circulation);
            }
        }
        let balance = CirculationBalance2 {
            cell_circulation: cells.value(),
            sheet_circulation: sheets.value(),
            boundary_circulation: boundary.value(),
        };
        if !balance.cell_circulation.is_finite()
            || !balance.sheet_circulation.is_finite()
            || !balance.boundary_circulation.is_finite()
            || !balance.stokes_defect().is_finite()
        {
            return Err(FluxKinematicsError::NonFiniteIntegral);
        }
        Ok(balance)
    }

    /// Integrate distributional curl against a continuous P1 scalar test.
    /// Supply one finite scalar value per mesh vertex. Cell and edge products
    /// are integrated analytically, including sheets with zero mean but
    /// nonzero first moment. A constant-one test gives total circulation.
    /// This returns a weak integral, not a smoothed bulk vorticity field.
    /// Work is linear and allocation-free; trailing pressure entries are ignored.
    ///
    /// # Errors
    /// Refuses invalid vector lengths, nonfinite scalar values/rates/integrals.
    pub fn weak_vorticity(
        &self,
        velocity: &[f64],
        scalar_vertices: &[f64],
    ) -> Result<f64, FluxKinematicsError> {
        self.check_kinematics_length(velocity)?;
        if scalar_vertices.len() != self.mesh.verts.len()
            || scalar_vertices.iter().any(|value| !value.is_finite())
        {
            return Err(FluxKinematicsError::Invalid(
                "one finite scalar test value per vertex required",
            ));
        }
        let mut integral = Sum::default();
        for (cell, vertices) in self.mesh.tris.iter().enumerate() {
            let flow = self
                .try_cell_kinematics(velocity, cell)
                .ok_or(FluxKinematicsError::NonFiniteCell(cell))?;
            let mean = vertices
                .iter()
                .map(|&v| scalar_vertices[v] / 3.0)
                .sum::<f64>();
            integral.add(self.mesh.areas[cell] * flow.vorticity() * mean);
        }
        for (edge, face) in self.mesh.edges.iter().enumerate() {
            if let Some(sheet) = self.vorticity_sheet(velocity, edge)? {
                let [a, b] = sheet.endpoint_strengths;
                let (p, q) = (scalar_vertices[face.verts.0], scalar_vertices[face.verts.1]);
                integral.add(face.len * (a * (p / 3.0 + q / 6.0) + b * (p / 6.0 + q / 3.0)));
            }
        }
        let result = integral.value();
        if !result.is_finite() {
            return Err(FluxKinematicsError::NonFiniteIntegral);
        }
        Ok(result)
    }

    fn check_kinematics_length(&self, velocity: &[f64]) -> Result<(), FluxKinematicsError> {
        if velocity.len() < self.n_u {
            return Err(FluxKinematicsError::Invalid("velocity dof count mismatch"));
        }
        Ok(())
    }

    fn kinematics_velocity(&self, velocity: &[f64], cell: usize, point: [f64; 2]) -> [f64; 2] {
        let mut value = [0.0; 2];
        for i in 0..6 {
            let basis = eval_basis(&self.bases[cell], i, point);
            let coefficient = velocity[self.dof(cell, i)];
            value[0] += coefficient * basis[0];
            value[1] += coefficient * basis[1];
        }
        value
    }

    fn kinematics_sheet_strength(&self, velocity: &[f64], edge: usize, point: [f64; 2]) -> f64 {
        let face = &self.mesh.edges[edge];
        let owner = self.kinematics_velocity(velocity, face.tris.0, point);
        let neighbour = self.kinematics_velocity(velocity, face.tris.1, point);
        (neighbour[1] - owner[1]) * face.normal[0] - (neighbour[0] - owner[0]) * face.normal[1]
    }
}

// Neumaier accumulation keeps cancellation of opposite-signed cell/sheet
// circulations from introducing an unnecessary mesh-size-dependent drift.
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}
