//! Distributed stranded windings on a fixed planar magnetic problem.
//!
//! Each profile is signed J_z / terminal_current [1/m2] in triangle order.
//! It already includes turns and the outgoing/return orientation. Area times
//! this density integrates to a signed turn count; volume times A_z times the
//! density is flux LINKAGE [Wb], including the problem's actual axial length.
//! Each profile must balance outgoing and return current on EVERY connected
//! cross-section component (each has its own constant-potential gauge).
//! End closure is outside this 2-D model: callers explicitly supply nonnegative
//! series end/leakage self-inductances [H], including zero when neglected.
//!
//! Unit-current solves reuse the existing FEEC operator. The symmetric magnetic
//! energy Gram matrix supplies active-section self/mutual inductance. Every
//! entry is compared with the independently evaluated A.J_unit pairing before
//! publication; the residual is an estimate, not a discretization certificate.
//! Geometry/material state stays fixed. No moving winding, nonlinear magnetic
//! law, skin/proximity loss, end mutual coupling or temperature-work law is
//! inferred. The optional `ohmic` owner derives intrinsic series resistance
//! from the same spatial material/profile data and retains explicit graph losses.

pub mod circuit;
/// Spatial material-derived winding resistance and native-circuit Joule heat.
pub mod ohmic;

use crate::{
    planar::{MagneticError, PlanarField, PlanarMagnetostatics, poll},
    transfer::{CellSample, PlanarFrame},
};
use fs_la::factor::cholesky;
use fs_solver::CheckedCgConfig;

/// Bound the small dense coupled magnetic block independently of mesh size.
pub const MAX_WINDINGS: usize = 16;

/// Allocation and measured consistency tolerances for winding extraction.
#[derive(Clone, Copy, Debug)]
pub struct WindingConfig {
    /// Maximum profiles, in 1..=MAX_WINDINGS.
    pub max_windings: usize,
    /// Ceiling on retained scalar basis values: windings*(2*vertices+7*cells),
    /// including each solved field's original impressed-current snapshot.
    pub max_basis_values: usize,
    /// Absolute source/energy flux-pairing agreement tolerance [H].
    pub pairing_tolerance_h: f64,
    /// Absolute outgoing-minus-return current tolerance [turns].
    pub current_balance_tolerance_turns: f64,
    /// Maximum row-sum residual of L*L_inverse-I; finite, in (0,1).
    pub inverse_residual_tolerance: f64,
}

/// Source, physical-model or bounded numerical refusal.
#[derive(Debug)]
pub enum WindingError {
    /// Existing material/FEEC/field-solve refusal.
    Magnetic(MagneticError),
    /// Invalid winding, correspondence, limits or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Unit-field energy and original source pairing do not agree.
    Pairing {
        /// Source-current basis index.
        source: usize,
        /// Receiving winding index.
        receiver: usize,
        /// Direct integral of A_source.J_receiver [H].
        flux_h: f64,
        /// Magnetic energy Gram entry [H].
        energy_h: f64,
    },
    /// The magnetic terminal block is singular or its inverse was not admitted.
    SingularInductance,
    /// Existing Kirchhoff/consistent-initial-state refusal.
    Circuit(fs_phs::circuit::CircuitError),
    /// Existing descriptor/quadratic-storage admission refusal.
    Phs(fs_phs::PhsError),
}
impl core::fmt::Display for WindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar winding: {self:?}")
    }
}
impl std::error::Error for WindingError {}
impl From<MagneticError> for WindingError {
    fn from(error: MagneticError) -> Self {
        Self::Magnetic(error)
    }
}
impl From<fs_phs::circuit::CircuitError> for WindingError {
    fn from(error: fs_phs::circuit::CircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<fs_phs::PhsError> for WindingError {
    fn from(error: fs_phs::PhsError) -> Self {
        Self::Phs(error)
    }
}

fn finite(value: f64) -> Result<f64, WindingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WindingError::Invalid(
            "winding arithmetic is not representable",
        ))
    }
}

/// Fixed-mesh/material unit-current fields and their reciprocal magnetic block.
/// Original problem, complete material snapshots, signed winding data and field
/// residual reports remain reachable. The physical Hamiltonian is
/// H(lambda)=lambda^T L^-1 lambda/2; terminal currents are L^-1 lambda.
pub struct PlanarWindingBasis<'a> {
    problem: &'a PlanarMagnetostatics,
    prescribed_vertices: Vec<usize>,
    profiles_per_m2: Vec<Vec<f64>>,
    fields: Vec<PlanarField<'a>>,
    active_inductance_h: Vec<f64>,
    inductance_h: Vec<f64>,
    inverse_inductance_per_h: Vec<f64>,
    end_inductance_h: Vec<f64>,
    net_turns_by_component: Vec<Vec<f64>>,
    maximum_pairing_defect_h: f64,
    inverse_residual: f64,
}
impl<'a> PlanarWindingBasis<'a> {
    /// Solve one 1-A excitation per profile with HOMOGENEOUS prescribed A_z.
    /// Boundary vertex selection uses the planar owner's complete-edge and
    /// connected-component admission. Background fields/permanent magnets and
    /// nonzero boundary bias are separate source-energy problems, not an L.
    ///
    /// End/leakage entries are caller-declared series self-inductances; zero
    /// explicitly neglects unrepresented ends. No extra leakage is invented to
    /// regularize dependent winding modes. A singular total block refuses.
    pub fn new(
        problem: &'a PlanarMagnetostatics,
        profiles_per_m2: Vec<Vec<f64>>,
        prescribed_vertices: &[usize],
        end_inductance_h: Vec<f64>,
        config: WindingConfig,
        linear: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindingError> {
        poll(cancelled)?;
        if problem.has_permanent_magnets() {
            return Err(WindingError::Invalid(
                "permanent-magnet bias cannot enter a homogeneous unit-current inductance basis",
            ));
        }
        let count = profiles_per_m2.len();
        let cells = problem.mesh().faces().len();
        let vertices = problem.mesh().vertices().len();
        let thickness_m = match problem.mesh().metric() {
            fs_rep_mesh::Metric2::Planar { thickness } => thickness,
            fs_rep_mesh::Metric2::Axisymmetric { .. } => unreachable!("admitted planar problem"),
        };
        let values = cells
            .checked_mul(7)
            .and_then(|n| vertices.checked_mul(2).and_then(|v| n.checked_add(v)))
            .and_then(|n| n.checked_mul(count));
        if count == 0
            || count > config.max_windings
            || config.max_windings == 0
            || config.max_windings > MAX_WINDINGS
            || values.is_none_or(|n| n > config.max_basis_values)
            || end_inductance_h.len() != count
            || prescribed_vertices.len() > vertices
            || !config.pairing_tolerance_h.is_finite()
            || config.pairing_tolerance_h < 0.0
            || !config.current_balance_tolerance_turns.is_finite()
            || config.current_balance_tolerance_turns < 0.0
            || !(config.inverse_residual_tolerance.is_finite()
                && config.inverse_residual_tolerance > 0.0
                && config.inverse_residual_tolerance < 1.0)
            || end_inductance_h.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(WindingError::Invalid(
                "winding count, storage, end model or tolerances are invalid",
            ));
        }
        let components = component_labels(problem, cancelled)?;
        let component_count = components.iter().copied().max().unwrap_or(0) + 1;
        let mut net_turns_by_component = Vec::with_capacity(count);
        for profile in &profiles_per_m2 {
            poll(cancelled)?;
            if profile.len() != cells || profile.iter().any(|v| !v.is_finite()) {
                return Err(WindingError::Invalid(
                    "one finite current density per triangle is required",
                ));
            }
            let mut net = vec![0.0; component_count];
            let mut positive = vec![0.0; component_count];
            let mut negative = vec![0.0; component_count];
            for (cell, &density) in profile.iter().enumerate() {
                poll(cancelled)?;
                let component = components[problem.mesh().faces()[cell][0] as usize];
                let turns =
                    finite(density * (problem.mesh().face_measure(cell).unwrap() / thickness_m))?;
                net[component] = finite(net[component] + turns)?;
                if turns > 0.0 {
                    positive[component] = finite(positive[component] + turns)?;
                }
                if turns < 0.0 {
                    negative[component] = finite(negative[component] - turns)?;
                }
            }
            if positive.iter().all(|&v| v == 0.0) {
                return Err(WindingError::Invalid(
                    "winding needs nonzero outgoing and return current",
                ));
            }
            for component in 0..component_count {
                if (positive[component] > 0.0) != (negative[component] > 0.0)
                    || net[component].abs() > config.current_balance_tolerance_turns
                {
                    return Err(WindingError::Invalid(
                        "unit current must close separately on every connected mesh component",
                    ));
                }
            }
            net_turns_by_component.push(net);
        }
        let boundary: Vec<_> = prescribed_vertices.iter().map(|&v| (v, 0.0)).collect();
        let mut fields = Vec::with_capacity(count);
        for profile in &profiles_per_m2 {
            poll(cancelled)?;
            fields.push(problem.solve(profile, &boundary, linear, cancelled)?);
        }
        let mut active = vec![0.0; count * count];
        for i in 0..count {
            for j in i..count {
                poll(cancelled)?;
                let mut value = 0.0;
                for cell in 0..cells {
                    poll(cancelled)?;
                    let a = fields[i].magnetic_flux_density_t()[cell];
                    let b = fields[j].magnetic_flux_density_t()[cell];
                    let volume = problem.mesh().face_measure(cell).unwrap();
                    let density =
                        finite(a[0].mul_add(b[0], a[1] * b[1]) / problem.permeability_h_m()[cell])?;
                    value = finite(volume.mul_add(density, value))?;
                }
                active[i * count + j] = value;
                active[j * count + i] = value;
            }
            if active[i * count + i] <= 0.0 {
                return Err(WindingError::Invalid(
                    "winding has no representable resolved magnetic energy",
                ));
            }
        }
        let mut maximum_pairing_defect_h = 0.0_f64;
        for i in 0..count {
            for j in 0..count {
                poll(cancelled)?;
                let flux = pairing(
                    problem,
                    fields[i].potential_wb_m(),
                    &profiles_per_m2[j],
                    cancelled,
                )?;
                let defect = finite(flux - active[i * count + j])?.abs();
                if defect > config.pairing_tolerance_h {
                    return Err(WindingError::Pairing {
                        source: i,
                        receiver: j,
                        flux_h: flux,
                        energy_h: active[i * count + j],
                    });
                }
                maximum_pairing_defect_h = maximum_pairing_defect_h.max(defect);
            }
        }
        let mut inductance_h = active.clone();
        for i in 0..count {
            inductance_h[i * count + i] =
                finite(inductance_h[i * count + i] + end_inductance_h[i])?;
        }
        let (inverse_inductance_per_h, inverse_residual) = inverse(
            &inductance_h,
            count,
            config.inverse_residual_tolerance,
            cancelled,
        )?;
        let mut prescribed_vertices = prescribed_vertices.to_vec();
        prescribed_vertices.sort_unstable();
        poll(cancelled)?;
        Ok(Self {
            problem,
            prescribed_vertices,
            profiles_per_m2,
            fields,
            active_inductance_h: active,
            inductance_h,
            inverse_inductance_per_h,
            end_inductance_h,
            net_turns_by_component,
            maximum_pairing_defect_h,
            inverse_residual,
        })
    }

    /// Exact original spatial material/geometry owner.
    pub fn problem(&self) -> &PlanarMagnetostatics {
        self.problem
    }
    /// Canonical vertices carrying the original homogeneous A_z boundary.
    /// Material replacement retains this physical boundary, not just the mesh.
    pub fn prescribed_vertices(&self) -> &[usize] {
        &self.prescribed_vertices
    }
    /// Number of terminal-current coordinates.
    pub fn winding_count(&self) -> usize {
        self.fields.len()
    }
    /// Signed J_z/terminal-current profiles [1/m2], in winding and cell order.
    pub fn profiles_per_m2(&self) -> &[Vec<f64>] {
        &self.profiles_per_m2
    }
    /// Original unit-current solutions, retaining their actual solver reports.
    pub fn unit_fields(&self) -> &[PlanarField<'a>] {
        &self.fields
    }
    /// Symmetric row-major active-section inductance [H].
    pub fn active_inductance_h(&self) -> &[f64] {
        &self.active_inductance_h
    }
    /// Total row-major inductance [H], including explicit end/leakage diagonals.
    pub fn inductance_h(&self) -> &[f64] {
        &self.inductance_h
    }
    /// Symmetric inverse used by the original circuit's Hamiltonian [1/H].
    pub fn inverse_inductance_per_h(&self) -> &[f64] {
        &self.inverse_inductance_per_h
    }
    /// Declared series end/leakage self-inductance [H], not a fitted field result.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.end_inductance_h
    }
    /// Measured signed closure residual [turns] for each winding and connected
    /// component. Components are ordered by their first mesh vertex.
    pub fn net_turns_by_component(&self) -> &[Vec<f64>] {
        &self.net_turns_by_component
    }
    /// Largest original source-pairing versus energy-Gram discrepancy [H].
    pub fn maximum_pairing_defect_h(&self) -> f64 {
        self.maximum_pairing_defect_h
    }
    /// Measured infinity-norm residual of the symmetric inverse.
    pub fn inverse_residual(&self) -> f64 {
        self.inverse_residual
    }

    /// Terminal flux linkage lambda=L*i [Wb], including declared end storage.
    pub fn flux_linkages_wb(&self, currents_a: &[f64]) -> Result<Vec<f64>, WindingError> {
        multiply(&self.inductance_h, currents_a, self.winding_count())
    }
    /// Terminal current i=L^-1*lambda [A].
    pub fn currents_a(&self, flux_linkages_wb: &[f64]) -> Result<Vec<f64>, WindingError> {
        multiply(
            &self.inverse_inductance_per_h,
            flux_linkages_wb,
            self.winding_count(),
        )
    }
    /// Evaluate the ORIGINAL distributed source pairing on any finite nodal A.
    /// This exposes constant-gauge and independently solved-field cross-checks.
    /// It includes only the active section; declared end storage is separate.
    pub fn pair_flux_linkages_wb(
        &self,
        potential_wb_m: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, WindingError> {
        self.profiles_per_m2
            .iter()
            .map(|profile| pairing(self.problem, potential_wb_m, profile, cancelled))
            .collect()
    }
    /// Reconstruct the same fixed linear field from terminal currents, then
    /// transfer it under the existing explicit triangle/position/frame contract.
    /// End/leakage energy has no field in this cross-section. This is still
    /// point sampling, with no inferred conservative voxel projection.
    pub fn sample_for_currents(
        &self,
        currents_a: &[f64],
        frame: &PlanarFrame,
        samples: &[CellSample],
        position_tolerance_m: f64,
        max_samples: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, WindingError> {
        if currents_a.len() != self.winding_count() || currents_a.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid(
                "one finite terminal current per winding is required",
            ));
        }
        let mut result = Vec::new();
        for (index, field) in self.fields.iter().enumerate() {
            poll(cancelled)?;
            let contribution =
                field.sample_cells(frame, samples, position_tolerance_m, max_samples, cancelled)?;
            if index == 0 {
                result = vec![[0.0; 3]; contribution.len()];
            }
            for (out, value) in result.iter_mut().zip(contribution) {
                for axis in 0..3 {
                    out[axis] = finite(currents_a[index].mul_add(value[axis], out[axis]))?;
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }
}

fn multiply(matrix: &[f64], vector: &[f64], n: usize) -> Result<Vec<f64>, WindingError> {
    if vector.len() != n || vector.iter().any(|v| !v.is_finite()) {
        return Err(WindingError::Invalid(
            "complete finite winding coordinates are required",
        ));
    }
    let mut result = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            result[i] = finite(matrix[i * n + j].mul_add(vector[j], result[i]))?;
        }
    }
    Ok(result)
}
fn pairing(
    problem: &PlanarMagnetostatics,
    potential: &[f64],
    profile: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, WindingError> {
    if potential.len() != problem.mesh().vertices().len()
        || potential.iter().any(|v| !v.is_finite())
    {
        return Err(WindingError::Invalid(
            "flux pairing needs a complete finite nodal potential",
        ));
    }
    let mut flux = 0.0;
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        poll(cancelled)?;
        // P1 exact cell average; divide before summing to limit intermediate scale.
        let mean = finite(face.iter().map(|&v| potential[v as usize] / 3.0).sum())?;
        let load = finite(profile[cell] * problem.mesh().face_measure(cell).unwrap())?;
        flux = finite(load.mul_add(mean, flux))?;
    }
    Ok(flux)
}
fn inverse(
    matrix: &[f64],
    n: usize,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64), WindingError> {
    poll(cancelled)?;
    let factor = cholesky(matrix, n).map_err(|_| WindingError::SingularInductance)?;
    let mut inverse = vec![0.0; n * n];
    for column in 0..n {
        poll(cancelled)?;
        let mut rhs = vec![0.0; n];
        rhs[column] = 1.0;
        factor.solve(&mut rhs);
        for row in 0..n {
            inverse[row * n + column] = finite(rhs[row])?;
        }
    }
    // Only floating-point inverse roundoff is averaged. The physical L was
    // already an exactly symmetric energy Gram matrix, checked against J.A.
    for i in 0..n {
        for j in 0..i {
            let value = 0.5 * inverse[i * n + j] + 0.5 * inverse[j * n + i];
            inverse[i * n + j] = value;
            inverse[j * n + i] = value;
        }
    }
    cholesky(&inverse, n).map_err(|_| WindingError::SingularInductance)?;
    let mut residual = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        let mut row = 0.0;
        for j in 0..n {
            let mut entry = if i == j { -1.0 } else { 0.0 };
            for k in 0..n {
                entry = finite(matrix[i * n + k].mul_add(inverse[k * n + j], entry))?;
            }
            row = finite(row + entry.abs())?;
        }
        residual = residual.max(row);
    }
    if residual > tolerance {
        return Err(WindingError::SingularInductance);
    }
    Ok((inverse, residual))
}

fn component_labels(
    problem: &PlanarMagnetostatics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, WindingError> {
    let n = problem.mesh().vertices().len();
    let mut neighbors = vec![Vec::new(); n];
    for &[a, b] in problem.mesh().edges() {
        poll(cancelled)?;
        neighbors[a as usize].push(b as usize);
        neighbors[b as usize].push(a as usize);
    }
    let mut labels = vec![usize::MAX; n];
    let mut stack = Vec::new();
    let mut component = 0;
    for start in 0..n {
        if labels[start] != usize::MAX {
            continue;
        }
        labels[start] = component;
        stack.push(start);
        while let Some(v) = stack.pop() {
            poll(cancelled)?;
            for &w in &neighbors[v] {
                if labels[w] == usize::MAX {
                    labels[w] = component;
                    stack.push(w);
                }
            }
        }
        component += 1;
    }
    Ok(labels)
}
