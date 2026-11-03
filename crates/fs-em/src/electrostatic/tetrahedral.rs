//! Linear dielectric weak forms on the shared three-dimensional P1 complex.
//!
//! `D = epsilon E + P`, `E = -grad(phi)`, and `div(D) = rho_free` give
//! `K phi = q_free + integral grad(N) . P`, with
//! `K = d0^T M_epsilon d0`. The same constant barycentric gradients supply
//! the polarization load and field reconstruction. This keeps the reciprocal
//! piezoelectric block above the electrical and mechanical field owners.
//!
//! Coefficients are explicit absolute SI tensors in the mesh coordinate frame.
//! They carry no measured-material authority here: the constitutive caller
//! retains the source, frame transform, validity domain and charge convention.
//! This numerical primitive does not choose electrode constraints, add Ohmic
//! conduction, or solve a coupled material law independently.

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
    whitney::{ElementGeometryError, try_element_geometry},
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

/// Explicit mesh/work envelope checked before geometry and matrix allocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetDielectricBudget {
    /// Maximum number of potential nodes.
    pub maximum_nodes: usize,
    /// Maximum number of tetrahedra; at most 36 edge-mass triplets per cell.
    pub maximum_tetrahedra: usize,
    /// Minimum `abs(det(J))/longest_edge^3` on every cell.
    pub minimum_scaled_jacobian: f64,
}

/// Mesh, material-weight, representability or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum TetDielectricError {
    /// Invalid bounded numerical input.
    Invalid(&'static str),
    /// Original checked Whitney geometry refusal.
    Geometry(ElementGeometryError),
    /// Original weighted-Hodge coefficient or assembly refusal.
    Assembly(WeightedError),
    /// Cancellation observed before publication.
    Cancelled,
}

impl core::fmt::Display for TetDielectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral dielectric: {self:?}")
    }
}
impl std::error::Error for TetDielectricError {}

fn invalid(reason: &'static str) -> TetDielectricError {
    TetDielectricError::Invalid(reason)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TetDielectricError> {
    if cancelled() {
        Err(TetDielectricError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, TetDielectricError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("dielectric arithmetic is not representable"))
    }
}
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

/// Prepared immutable dielectric operator and its original element metrics.
#[derive(Clone, Debug)]
pub struct TetDielectricOperator {
    node_count: usize,
    tetrahedra: Vec<[u32; 4]>,
    permittivity_f_m: Vec<[[f64; 3]; 3]>,
    gradients: Vec<[[f64; 3]; 4]>,
    volumes_m3: Vec<f64>,
    stiffness: Csr,
}

impl TetDielectricOperator {
    /// Assemble the native FEEC dielectric matrix once for repeated coupling.
    ///
    /// Tensor SPD admission uses the original interval-pivot weighted-Hodge
    /// check. Negative stored cell orientation is allowed; physical volumes
    /// are positive and gradients preserve the original vertex order.
    /// Cancellation is polled during admission and weighted assembly, around
    /// the bounded native batched geometry/sparse products, and at publication.
    pub fn assemble(
        mesh: &TetComplex,
        nodes_m: &[[f64; 3]],
        permittivity_f_m: &[[[f64; 3]; 3]],
        budget: TetDielectricBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, TetDielectricError> {
        poll(cancelled)?;
        let n = mesh.vertex_count;
        let cells = mesh.tets.len();
        let triplets = cells
            .checked_mul(36)
            .ok_or_else(|| invalid("cell extent overflow"))?;
        let maximum_edges = cells
            .checked_mul(6)
            .ok_or_else(|| invalid("edge extent overflow"))?;
        if n < 4
            || cells == 0
            || n > budget.maximum_nodes
            || cells > budget.maximum_tetrahedra
            || nodes_m.len() != n
            || permittivity_f_m.len() != cells
            || mesh.edges.len() > maximum_edges
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
        {
            return Err(invalid(
                "nonempty bounded mesh, one absolute tensor per cell and positive quality threshold required",
            ));
        }
        // d0 only consumes this public edge table. Reject malformed edges
        // before passing them to the native CSR constructors.
        for (i, &[a, b]) in mesh.edges.iter().enumerate() {
            poll(cancelled)?;
            if a >= b || b as usize >= n || i > 0 && mesh.edges[i - 1] >= [a, b] {
                return Err(invalid("canonical valid dielectric edge table required"));
            }
        }
        let mut incident = vec![false; n];
        for tet in &mesh.tets {
            poll(cancelled)?;
            for (local, &node) in tet.iter().enumerate() {
                if node as usize >= n || tet[..local].contains(&node) {
                    return Err(invalid("valid distinct tetrahedral vertices required"));
                }
                incident[node as usize] = true;
            }
        }
        if incident.iter().any(|seen| !seen) {
            return Err(invalid("isolated dielectric nodes are not admitted"));
        }
        let geo = try_element_geometry(mesh, nodes_m).map_err(TetDielectricError::Geometry)?;
        poll(cancelled)?;
        let mut volumes_m3 = Vec::with_capacity(cells);
        for (cell, tet) in mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            let mut longest_squared = 0.0_f64;
            for a in 0..4 {
                for b in a + 1..4 {
                    let edge = core::array::from_fn(|j| {
                        nodes_m[tet[a] as usize][j] - nodes_m[tet[b] as usize][j]
                    });
                    longest_squared = longest_squared.max(dot(&edge, &edge));
                }
            }
            let longest = fs_math::det::sqrt(longest_squared);
            let volume = geo.vol_signed[cell].abs();
            let quality = (6.0 * volume / longest) / longest_squared;
            if !quality.is_finite() || quality < budget.minimum_scaled_jacobian {
                return Err(invalid(
                    "dielectric tetrahedron is outside its geometric quality envelope",
                ));
            }
            volumes_m3.push(volume);
        }
        let weights: Vec<_> = permittivity_f_m
            .iter()
            .copied()
            .map(CellWeight::Tensor)
            .collect();
        let mass = weighted_mass_matrix(
            mesh,
            &geo,
            1,
            &weights,
            ELECTRIC_PERMITTIVITY_DIMS,
            WeightedAssemblyLimits {
                max_cells: budget.maximum_tetrahedra,
                max_dofs: maximum_edges,
                max_triplets: triplets,
            },
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => TetDielectricError::Cancelled,
            other => TetDielectricError::Assembly(other),
        })?;
        poll(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&mesh.d0()), mass.matrix());
        for row in 0..n {
            poll(cancelled)?;
            if stiffness.row(row).1.iter().any(|v| !v.is_finite()) {
                return Err(invalid("assembled dielectric stiffness is not finite"));
            }
        }
        poll(cancelled)?;
        Ok(Self {
            node_count: n,
            tetrahedra: mesh.tets.clone(),
            permittivity_f_m: permittivity_f_m.to_vec(),
            gradients: geo.grads,
            volumes_m3,
            stiffness,
        })
    }

    /// Original vertex count, including constrained electrode nodes.
    pub const fn node_count(&self) -> usize {
        self.node_count
    }
    /// Original connectivity and local vertex ordering.
    pub fn tetrahedra(&self) -> &[[u32; 4]] {
        &self.tetrahedra
    }
    /// Global-frame constant barycentric gradients [1/m].
    pub fn gradients(&self) -> &[[[f64; 3]; 4]] {
        &self.gradients
    }
    /// Positive physical reference-cell volumes [m^3].
    pub fn volumes_m3(&self) -> &[f64] {
        &self.volumes_m3
    }
    /// Original absolute dielectric tensors, one per cell [F/m].
    pub fn permittivity_f_m(&self) -> &[[[f64; 3]; 3]] {
        &self.permittivity_f_m
    }
    /// Full unconstrained nodal matrix `d0^T M_epsilon d0` [F].
    pub const fn stiffness(&self) -> &Csr {
        &self.stiffness
    }

    /// Work-conjugate nodal polarization load `sum V grad(N) . P` [C].
    ///
    /// Add this to the free-charge RHS of `K phi = q_free + load`.
    /// Thus `K phi - load` is the nodal free-charge reaction; for a surface
    /// electrode it equals minus the outward dielectric displacement flux.
    pub fn polarization_load_c(
        &self,
        polarization_c_m2: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, TetDielectricError> {
        poll(cancelled)?;
        if polarization_c_m2.len() != self.tetrahedra.len()
            || polarization_c_m2.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite polarization vector per dielectric cell required",
            ));
        }
        let mut load = vec![0.0; self.node_count];
        for (cell, tet) in self.tetrahedra.iter().enumerate() {
            poll(cancelled)?;
            for (a, &node) in tet.iter().enumerate() {
                let value = finite(
                    self.volumes_m3[cell] * dot(&self.gradients[cell][a], &polarization_c_m2[cell]),
                )?;
                load[node as usize] = finite(load[node as usize] + value)?;
            }
        }
        poll(cancelled)?;
        Ok(load)
    }

    /// Recover `E = -grad(phi)` using voltage differences within each cell.
    /// A constant reference voltage therefore cannot create a spurious field.
    pub fn electric_field_v_m(
        &self,
        voltage_v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, TetDielectricError> {
        poll(cancelled)?;
        if voltage_v.len() != self.node_count || voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid("one finite electric potential per node required"));
        }
        let mut fields = Vec::with_capacity(self.tetrahedra.len());
        for (cell, tet) in self.tetrahedra.iter().enumerate() {
            poll(cancelled)?;
            let mut field = [0.0; 3];
            for a in 1..4 {
                let difference = finite(voltage_v[tet[a] as usize] - voltage_v[tet[0] as usize])?;
                for (j, value) in field.iter_mut().enumerate() {
                    *value = finite((-self.gradients[cell][a][j]).mul_add(difference, *value))?;
                }
            }
            fields.push(field);
        }
        poll(cancelled)?;
        Ok(fields)
    }

    /// Recover the original `D = epsilon E + P` separately in every cell.
    /// No averaging across material interfaces is performed.
    pub fn displacement_c_m2(
        &self,
        voltage_v: &[f64],
        polarization_c_m2: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, TetDielectricError> {
        if polarization_c_m2.len() != self.tetrahedra.len()
            || polarization_c_m2.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite polarization vector per dielectric cell required",
            ));
        }
        let mut displacement = self.electric_field_v_m(voltage_v, cancelled)?;
        for (cell, value) in displacement.iter_mut().enumerate() {
            poll(cancelled)?;
            let electric = *value;
            for (j, component) in value.iter_mut().enumerate() {
                *component = finite(
                    dot(&self.permittivity_f_m[cell][j], &electric) + polarization_c_m2[cell][j],
                )?;
            }
        }
        poll(cancelled)?;
        Ok(displacement)
    }
}
