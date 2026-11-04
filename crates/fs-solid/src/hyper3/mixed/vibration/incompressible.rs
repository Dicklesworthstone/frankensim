//! Physical modes of the original exactly constrained Taylor–Hood solid.
//!
//! The accepted equilibrium supplies [Kuu, B^T; B, 0]. Native Householder QR
//! supplies the complete displacement nullspace Z, and fs-modal solves
//! (Z^T Kuu Z) q = lambda (Z^T M Z) q. Pressure has ZERO inertia. Its shape is
//! recovered from B^T psi = lambda M phi - Kuu phi, with phi=Z q.
//!
//! Only the pressure row already removed by the original, admitted gauge is
//! omitted from QR. There is no local rank truncation, bulk penalty or pressure
//! mass. Every original weak pressure equation, including the gauge row, is
//! checked after lifting; the full free-displacement momentum equation is also
//! checked. A gauge fixes the pressure increment there to zero, not the original
//! equilibrium pressure to a new value.
//!
//! Exact refers to the original finite-element weak constraint at the held
//! constitutive/source state, not pointwise continuum incompressibility. A
//! thermal material can therefore retain its own reference volume. Original
//! geometry, quadrature, dead loads, absolute supports and stable-KKT admission
//! remain in force. This bounded dense route admits at most 128 FREE
//! displacement DOFs BEFORE constraint reduction. Signed eigenvalues and their
//! projected-pencil intervals are retained; unstable/singular equilibria still
//! refuse through the original static owner.

use super::super::pressure::UniformBoundaryPressure;
use super::super::{
    MixedUpError, MixedUpProblem, MixedUpSolution, PreparedMixedUpProblem, checkpoint,
    inertia_error, invalid,
};
use super::{MixedVibrationError, MixedVibrationSettings, symmetric_dense, tolerance};
use crate::hyper3::vibration::MAX_DENSE_FREE_DOFS;
use crate::inertia::ReferenceMass;
use fs_exec::Cx;
use fs_la::factor::qr;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_modal::{ModePair, eigh_gen_dense};
use fs_sparse::Csr;

/// One constrained mode, with separate native projected-pencil evidence.
#[derive(Clone, Debug)]
pub struct IncompressibleVibrationMode {
    /// Native mass-normalized mode in Z coordinates. Its residual and interval
    /// refer to the retained projected pencil, not the unreduced Kuu alone.
    pub coordinates: ModePair,
    /// Lifted displacement shape, in free_dofs order, normalized with the
    /// ORIGINAL reference mass. Multiply by the modal coordinate to get metres.
    pub displacement_shape: Vec<f64>,
    /// Original pressure-node order, with zero increment at any declared gauge.
    /// Multiply by the SAME modal coordinate as displacement_shape to get Pa.
    /// There is no separate pressure mass or normalization.
    pub pressure_shape: Vec<f64>,
    /// ||B_all phi||_infinity / max_i sum_j |B_all_ij phi_j|.
    /// All original pressure rows participate, including the gauge row.
    pub constraint_relative_residual: f64,
    /// Free momentum backward error, scaled by uncancelled terms from Kuu phi,
    /// Kup psi and lambda M phi rather than by their nearly cancelling actions.
    pub momentum_relative_residual: f64,
}

/// Immutable physical problem and pencil, the complete constrained subspace,
/// and every mode of its bounded dense projected pencil.
#[derive(Debug)]
pub struct IncompressibleMixedVibration<
    'a,
    L: HyperelasticMaterial + MixedFiniteStrainPotential = Hyperelastic,
> {
    problem: MixedUpProblem<'a, L>,
    boundary_pressure: Option<UniformBoundaryPressure>,
    settings: MixedVibrationSettings,
    reference_density_kg_m3: Vec<f64>,
    equilibrium: MixedUpSolution,
    free_dofs: Vec<usize>,
    independent_pressure_indices: Vec<usize>,
    mixed_tangent: Csr,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    basis: Vec<f64>,
    projected_stiffness_n_m: Vec<f64>,
    projected_mass_kg: Vec<f64>,
    modes: Vec<IncompressibleVibrationMode>,
}

fn constraint_residual(
    full: &Csr,
    nu: usize,
    free: &[usize],
    phi: &[f64],
    cx: &Cx<'_>,
) -> Result<f64, MixedUpError> {
    let (mut residual, mut scale) = (0.0_f64, 0.0_f64);
    for row in nu..full.nrows() {
        checkpoint(cx)?;
        let (mut action, mut magnitude) = (0.0_f64, 0.0_f64);
        for (&dof, &value) in free.iter().zip(phi) {
            let term = full.get(row, dof) * value;
            action += term;
            magnitude += term.abs();
        }
        if !action.is_finite() || !magnitude.is_finite() {
            return Err(invalid("nonfinite incompressible constraint residual"));
        }
        residual = residual.max(action.abs());
        scale = scale.max(magnitude);
    }
    Ok(if scale == 0.0 { 0.0 } else { residual / scale })
}

/// Small dense Galerkin products only; QR and the eigensolve remain native.
fn project(
    matrix: &[f64],
    basis: &[f64],
    n: usize,
    r: usize,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, MixedUpError> {
    let mut action = vec![0.0; n * r];
    for i in 0..n {
        checkpoint(cx)?;
        for j in 0..n {
            for k in 0..r {
                action[i * r + k] = matrix[i * n + j].mul_add(basis[j * r + k], action[i * r + k]);
            }
        }
    }
    let mut reduced = vec![0.0; r * r];
    for i in 0..r {
        checkpoint(cx)?;
        for j in 0..r {
            for k in 0..n {
                reduced[i * r + j] =
                    basis[k * r + i].mul_add(action[k * r + j], reduced[i * r + j]);
            }
        }
    }
    symmetric_dense(&reduced, r)?;
    Ok(reduced)
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    /// Equilibrate the original exact-constraint solid and return its complete
    /// constrained spectrum within the 128-free-displacement dense cap.
    ///
    /// The original pressure gauge, if any, must be admitted by the native
    /// static owner. No constraint is dropped to rescue a rank-deficient KKT.
    pub fn solve_incompressible_vibration(
        &self,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepare(cx)?
            .solve_incompressible_vibration_impl(None, density, settings, cx)
    }

    /// Re-equilibrate a starting iterate under the original absolute supports,
    /// held dead loads and exact constraint before constructing its modes.
    pub fn solve_incompressible_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepare(cx)?.solve_incompressible_vibration_impl(
            Some((displacement_m, pressure_pa)),
            density,
            settings,
            cx,
        )
    }
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> PreparedMixedUpProblem<'a, L> {
    pub(in crate::hyper3::mixed) fn solve_incompressible_vibration_impl(
        &self,
        initial: Option<(&[[f64; 3]], &[f64])>,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        checkpoint(cx)?;
        let prepared = self;
        let problem = &self.problem;
        if problem.volumetric != VolumetricTreatment::ExactIncompressible {
            return Err(
                invalid("nullspace vibration requires the original exact constraint").into(),
            );
        }
        if !tolerance(settings.spectrum.relative_residual_tolerance)
            || !tolerance(settings.pressure_residual_tolerance)
        {
            return Err(
                invalid("modal and pressure tolerances must be finite and in (0,1)").into(),
            );
        }
        let n = prepared.free_displacements;
        let nu = 3 * prepared.mesh.nodes_count;
        let np = prepared.mesh.pressure_nodes.len();
        let m = prepared.free.len() - n;
        if n == 0 || n > settings.spectrum.maximum_free_dofs || n > MAX_DENSE_FREE_DOFS {
            return Err(invalid(
                "incompressible dense vibration requires 1..=128 free displacement DOFs",
            )
            .into());
        }
        if m == 0 || m > n {
            return Err(invalid(
                "independent pressure constraints must fit the free displacement space",
            )
            .into());
        }
        // Before extra dense allocations: QR, original/projected pencils,
        // basis, conservative eigensolver scratch/output and lifted pressures.
        // Original assembly, geometry and KKT factor budgets remain separate.
        let bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(48))
            .and_then(|v| {
                np.checked_mul(n)
                    .and_then(|q| q.checked_mul(4))
                    .and_then(|q| v.checked_add(q))
            })
            .and_then(|v| np.checked_mul(16).and_then(|q| v.checked_add(q)))
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .ok_or_else(|| invalid("incompressible modal storage overflow"))?;
        if bytes > settings.maximum_condensation_bytes {
            return Err(invalid(format!(
                "incompressible modal reduction needs at most {bytes} payload bytes; cap {}",
                settings.maximum_condensation_bytes
            ))
            .into());
        }
        let free_dofs = prepared.free[..n].to_vec();
        let independent_pressure_indices: Vec<_> =
            prepared.free[n..].iter().map(|&d| d - nu).collect();
        let mass = ReferenceMass::from_tet10_quadrature(
            prepared.mesh.elements.iter().map(|element| {
                (
                    element.connectivity,
                    element.points.iter().map(|q| (q.shapes, q.weight_m3)),
                )
            }),
            density,
            nu,
            settings.maximum_mass_assembly_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        let mut mass_kg = Vec::with_capacity(n * n);
        for &row in &free_dofs {
            checkpoint(cx)?;
            for &col in &free_dofs {
                mass_kg.push(mass.matrix().get(row, col));
            }
        }
        symmetric_dense(&mass_kg, n)?;
        drop(mass);
        // This solve retains final KKT inertia admission, even for a starting
        // iterate that already satisfies equilibrium exactly.
        let equilibrium = match initial {
            None => prepared.solve(cx)?,
            Some((u, p)) => prepared.solve_from(u, p, cx)?,
        };
        let full = prepared
            .evaluate(
                &equilibrium.displacement_m,
                &equilibrium.pressure_pa,
                true,
                cx,
            )?
            .tangent
            .ok_or_else(|| invalid("native mixed tangent is missing"))?;
        let mut stiffness_n_m = Vec::with_capacity(n * n);
        for &row in &free_dofs {
            checkpoint(cx)?;
            for &col in &free_dofs {
                stiffness_n_m.push(full.get(row, col));
            }
        }
        symmetric_dense(&stiffness_n_m, n)?;
        for p in 0..np {
            checkpoint(cx)?;
            let (columns, values) = full.row(nu + p);
            for (&col, &value) in columns.iter().zip(values) {
                if col >= nu && value != 0.0 {
                    return Err(invalid("exact constraint has nonzero pressure compliance").into());
                }
            }
            for &dof in &free_dofs {
                let b = full.get(nu + p, dof);
                if !b.is_finite() || b != full.get(dof, nu + p) {
                    return Err(invalid(
                        "exact mixed pressure/displacement blocks must be finite transposes",
                    )
                    .into());
                }
            }
        }
        let mut transpose = Vec::with_capacity(n * m);
        for &dof in &free_dofs {
            checkpoint(cx)?;
            for &p in &independent_pressure_indices {
                transpose.push(full.get(nu + p, dof));
            }
        }
        // One global scaling protects native QR's sums of squares. It changes
        // no row, pressure coordinate or nullspace. Native KKT already admits
        // every independent pressure row; no numerical rank threshold is used.
        let constraint_scale = transpose.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if constraint_scale == 0.0 || !constraint_scale.is_finite() {
            return Err(invalid("unrepresentable exact pressure constraint").into());
        }
        for value in &mut transpose {
            *value /= constraint_scale;
        }
        checkpoint(cx)?;
        let factor = qr(&transpose, n, m);
        drop(transpose);
        checkpoint(cx)?;
        for i in 0..m {
            if !factor.r(i, i).is_finite() || factor.r(i, i) == 0.0 {
                return Err(
                    invalid("native constraint QR cannot represent all admitted rows").into(),
                );
            }
        }
        let r = n - m;
        let mut basis = vec![0.0; n * r];
        for column in 0..r {
            checkpoint(cx)?;
            let mut vector = vec![0.0; n];
            vector[m + column] = 1.0;
            factor.apply_q(&mut vector);
            if vector.iter().any(|v| !v.is_finite()) {
                return Err(invalid("nonfinite exact-constraint nullspace").into());
            }
            let relative = constraint_residual(&full, nu, &free_dofs, &vector, cx)?;
            if relative > settings.pressure_residual_tolerance {
                return Err(MixedVibrationError::PressureResidual {
                    stage: "nullspace constraints",
                    index: column,
                    relative,
                });
            }
            for row in 0..n {
                basis[row * r + column] = vector[row];
            }
        }
        let projected_stiffness_n_m = project(&stiffness_n_m, &basis, n, r, cx)?;
        let projected_mass_kg = project(&mass_kg, &basis, n, r, cx)?;
        checkpoint(cx)?;
        let projected_modes = eigh_gen_dense(&projected_stiffness_n_m, &projected_mass_kg, r)
            .map_err(MixedVibrationError::Modal)?;
        checkpoint(cx)?;
        if projected_modes.len() != r {
            return Err(invalid("incomplete constrained modal basis").into());
        }
        let eigenvalue_scale = projected_modes
            .iter()
            .map(|mode| mode.lambda.abs())
            .fold(0.0, f64::max);
        let mut modes = Vec::with_capacity(r);
        for (index, coordinates) in projected_modes.into_iter().enumerate() {
            checkpoint(cx)?;
            if coordinates.phi.len() != r
                || coordinates.phi.iter().any(|v| !v.is_finite())
                || !coordinates.lambda.is_finite()
                || !coordinates.residual.is_finite()
                || coordinates.residual < 0.0
                || !coordinates.interval.0.is_finite()
                || !coordinates.interval.1.is_finite()
                || coordinates.residual
                    > settings.spectrum.relative_residual_tolerance * eigenvalue_scale
            {
                return Err(MixedVibrationError::ModeResidual {
                    index,
                    eigenvalue_per_s2: coordinates.lambda,
                    residual_per_s2: coordinates.residual,
                });
            }
            let displacement_shape: Vec<_> = (0..n)
                .map(|i| {
                    (0..r)
                        .map(|j| basis[i * r + j] * coordinates.phi[j])
                        .sum::<f64>()
                })
                .collect();
            let constraint_relative_residual =
                constraint_residual(&full, nu, &free_dofs, &displacement_shape, cx)?;
            if constraint_relative_residual > settings.pressure_residual_tolerance {
                return Err(MixedVibrationError::PressureResidual {
                    stage: "mode constraints",
                    index,
                    relative: constraint_relative_residual,
                });
            }
            let mut rhs = vec![0.0; n];
            let mut magnitudes = vec![0.0; n];
            for i in 0..n {
                checkpoint(cx)?;
                let (mut stiffness_action, mut mass_action) = (0.0_f64, 0.0_f64);
                for j in 0..n {
                    let k = stiffness_n_m[i * n + j] * displacement_shape[j];
                    let mass = mass_kg[i * n + j] * displacement_shape[j];
                    stiffness_action += k;
                    mass_action += mass;
                    magnitudes[i] += k.abs() + (coordinates.lambda * mass).abs();
                }
                rhs[i] = coordinates.lambda.mul_add(mass_action, -stiffness_action);
            }
            // (B^T/s) z = rhs, hence the ORIGINAL pressure shape is psi=z/s.
            // Q also removes the small projected eigen-residual from this LS
            // solve; the full original momentum equation is checked below.
            let recovered = factor.solve_ls(&rhs);
            let mut pressure_shape = vec![0.0; np];
            for (&p, value) in independent_pressure_indices.iter().zip(recovered) {
                pressure_shape[p] = value / constraint_scale;
            }
            if pressure_shape.iter().any(|v| !v.is_finite()) {
                return Err(invalid("nonfinite recovered incompressible pressure mode").into());
            }
            let (mut residual, mut scale) = (0.0_f64, 0.0_f64);
            for (i, &dof) in free_dofs.iter().enumerate() {
                checkpoint(cx)?;
                let mut action = 0.0_f64;
                for (p, &value) in pressure_shape.iter().enumerate() {
                    let term = full.get(dof, nu + p) * value;
                    action += term;
                    magnitudes[i] += term.abs();
                }
                if !action.is_finite() || !rhs[i].is_finite() || !magnitudes[i].is_finite() {
                    return Err(invalid("nonfinite constrained momentum residual").into());
                }
                residual = residual.max((action - rhs[i]).abs());
                scale = scale.max(magnitudes[i]);
            }
            let momentum_relative_residual = if scale == 0.0 { 0.0 } else { residual / scale };
            if momentum_relative_residual > settings.spectrum.relative_residual_tolerance {
                return Err(MixedVibrationError::PressureResidual {
                    stage: "mode momentum",
                    index,
                    relative: momentum_relative_residual,
                });
            }
            modes.push(IncompressibleVibrationMode {
                coordinates,
                displacement_shape,
                pressure_shape,
                constraint_relative_residual,
                momentum_relative_residual,
            });
        }
        let result = IncompressibleMixedVibration {
            problem: *problem,
            boundary_pressure: prepared.boundary_pressure,
            settings,
            reference_density_kg_m3: density.to_vec(),
            equilibrium,
            free_dofs,
            independent_pressure_indices,
            mixed_tangent: full,
            stiffness_n_m,
            mass_kg,
            basis,
            projected_stiffness_n_m,
            projected_mass_kg,
            modes,
        };
        checkpoint(cx)?;
        Ok(result)
    }
}

impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> IncompressibleMixedVibration<'_, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'_, L> {
        &self.problem
    }
    /// Held load on the complete material-domain boundary, if selected. Its
    /// q Hessian(V) is retained in Kuu, while equilibrium.potential_j remains
    /// the original material mixed potential. The constraint and pressure
    /// convention are unchanged.
    pub const fn boundary_pressure(&self) -> Option<UniformBoundaryPressure> {
        self.boundary_pressure
    }
    pub const fn settings(&self) -> MixedVibrationSettings {
        self.settings
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    pub const fn equilibrium(&self) -> &MixedUpSolution {
        &self.equilibrium
    }
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Indices into equilibrium.pressure_nodes retained by native gauge handling.
    pub fn independent_pressure_indices(&self) -> &[usize] {
        &self.independent_pressure_indices
    }
    /// Original full physical mixed Hessian: all displacement rows, then all
    /// original pressure nodes. No gauge row or constrained support is erased.
    pub const fn mixed_tangent(&self) -> &Csr {
        &self.mixed_tangent
    }
    /// Row-major free Kuu [N/m], including the original equilibrium pressure's
    /// geometric tangent. This matrix alone is NOT the constrained mode pencil.
    pub fn stiffness_n_m(&self) -> &[f64] {
        &self.stiffness_n_m
    }
    /// Row-major original consistent P2 reference displacement mass [kg].
    pub fn mass_kg(&self) -> &[f64] {
        &self.mass_kg
    }
    /// Number of physical modes n_free_displacements-n_independent_pressures.
    /// It can be zero for a fully constrained displacement space.
    pub fn dimension(&self) -> usize {
        self.free_dofs.len() - self.independent_pressure_indices.len()
    }
    /// Row-major Z, with free_dofs.len() rows and dimension() columns.
    /// Its Euclidean-orthonormal columns span the admitted constraint nullspace.
    pub fn basis(&self) -> &[f64] {
        &self.basis
    }
    pub fn projected_stiffness_n_m(&self) -> &[f64] {
        &self.projected_stiffness_n_m
    }
    pub fn projected_mass_kg(&self) -> &[f64] {
        &self.projected_mass_kg
    }
    pub fn modes(&self) -> &[IncompressibleVibrationMode] {
        &self.modes
    }
}
