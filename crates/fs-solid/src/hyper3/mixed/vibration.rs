//! Physical vibration of the original finite-compliance Taylor–Hood solid.
//!
//! The accepted mixed equilibrium owns raw physical tangent blocks. Pressure
//! has no inertia: H=-Kpp, R=H^-1 Kpu, dp=R du and K_eff=Kuu+Kup R. A single
//! native sparse factor supplies all pressure solves; fs-modal owns the bounded
//! complete dense displacement spectrum. No pressure penalty, bulk replacement,
//! lumped quadratic mass or additional initial-stress matrix is introduced.
//!
//! Only OriginalMaterial on its admitted positive bulk-curvature branch is
//! supported. Exact incompressibility requires a constraint/nullspace treatment
//! and refuses here. The original static owner also refuses unstable/singular
//! KKT equilibria; this adapter does not bypass that gate. Returned eigenvalues
//! remain signed, with their residual intervals, without frequency clipping.
//!
//! Constitutive/source conditions and dead loads are held fixed. The original
//! P2/P1 weak pressure space, geometry and quadrature remain authoritative;
//! curved geometry and nonlinear coefficients still require quadrature and
//! mesh convergence. This route is not an inf-sup or continuum-stability proof.

use super::{MixedUpError, MixedUpProblem, MixedUpSolution, checkpoint, inertia_error, invalid};
use crate::hyper3::vibration::{MAX_DENSE_FREE_DOFS, VibrationSettings};
use crate::inertia::ReferenceMass;
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::{Coo, Csr, DirectOrdering, Inertia, LdltOptions, SymbolicLdlt};

/// Explicit modal accuracy and extra storage admission. Original mixed mesh,
/// quadrature, Newton, assembly and sparse-factor budgets remain in force.
#[derive(Clone, Copy, Debug)]
pub struct MixedVibrationSettings {
    pub spectrum: VibrationSettings,
    /// Native full signed P2 reference-mass assembly payload, separate from KKT.
    pub maximum_mass_assembly_bytes: usize,
    /// Conservative pressure-block staging, dense work and retained modal
    /// payload. Excludes the original mixed assembly and native sparse factor,
    /// which retain their own budgets; this is not an OS memory/RSS limit.
    pub maximum_condensation_bytes: usize,
    /// Relative backward-error gate for pressure elimination and recovered
    /// pressure equations, finite and strictly between zero and one.
    pub pressure_residual_tolerance: f64,
}

#[derive(Debug)]
pub enum MixedVibrationError {
    Mixed(MixedUpError),
    Modal(ModalError),
    PressureBlockInertia {
        got: Inertia,
    },
    PressureResidual {
        stage: &'static str,
        index: usize,
        relative: f64,
    },
    ModeResidual {
        index: usize,
        eigenvalue_per_s2: f64,
        residual_per_s2: f64,
    },
}
impl From<MixedUpError> for MixedVibrationError {
    fn from(error: MixedUpError) -> Self {
        Self::Mixed(error)
    }
}
impl core::fmt::Display for MixedVibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mixed solid vibration: {self:?}")
    }
}
impl core::error::Error for MixedVibrationError {}

/// One mass-normalized displacement mode and its algebraic pressure response.
#[derive(Debug, Clone)]
pub struct MixedVibrationMode {
    pub displacement: ModePair,
    /// Pressure-node order is equilibrium.pressure_nodes. Multiply this and
    /// displacement.phi by the SAME modal coordinate to obtain Pa and metres.
    /// This pressure shape has no independent kinetic energy or normalization.
    pub pressure_shape: Vec<f64>,
    /// ||H dp-Kpu phi||_infinity divided by the maximum row sum of
    /// |Kpu_ij phi_j| + |H_ip dp_p|, before cancellation within either action.
    pub pressure_equation_relative_residual: f64,
}

/// Immutable original problem, accepted equilibrium and full physical pencil.
#[derive(Debug)]
pub struct MixedUpVibration<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential = Hyperelastic>
{
    problem: MixedUpProblem<'a, L>,
    settings: MixedVibrationSettings,
    reference_density_kg_m3: Vec<f64>,
    equilibrium: MixedUpSolution,
    free_dofs: Vec<usize>,
    mixed_tangent: Csr,
    pressure_response_pa_per_m: Vec<f64>,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    modes: Vec<MixedVibrationMode>,
}

fn tolerance(value: f64) -> bool {
    value.is_finite() && value > 0.0 && value < 1.0
}

fn symmetric_dense(matrix: &[f64], n: usize) -> Result<(), MixedUpError> {
    if matrix.iter().any(|v| !v.is_finite()) {
        return Err(invalid("nonfinite mixed vibration pencil"));
    }
    let scale = matrix.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let tol = 64.0 * f64::EPSILON * n as f64 * scale;
    for i in 0..n {
        for j in 0..i {
            if (matrix[i * n + j] - matrix[j * n + i]).abs() > tol {
                return Err(invalid(
                    "mixed vibration pencil is not symmetric to roundoff",
                ));
            }
        }
    }
    Ok(())
}

fn pressure_residual(
    h: &Csr,
    x: &[f64],
    rhs: &[f64],
    rhs_magnitudes: Option<&[f64]>,
    cx: &Cx<'_>,
) -> Result<f64, MixedUpError> {
    let (mut residual, mut scale) = (0.0_f64, 0.0_f64);
    for row in 0..h.nrows() {
        checkpoint(cx)?;
        let (columns, values) = h.row(row);
        let (mut action, mut magnitude) = (
            0.0_f64,
            rhs_magnitudes.map_or_else(|| rhs[row].abs(), |magnitudes| magnitudes[row]),
        );
        for (&col, &value) in columns.iter().zip(values) {
            action = value.mul_add(x[col], action);
            magnitude += (value * x[col]).abs();
        }
        if !action.is_finite() || !magnitude.is_finite() || !rhs[row].is_finite() {
            return Err(invalid("nonfinite pressure-elimination residual"));
        }
        residual = residual.max((action - rhs[row]).abs());
        scale = scale.max(magnitude);
    }
    Ok(if scale == 0.0 { 0.0 } else { residual / scale })
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    /// Equilibrate the original mixed solid, eliminate its pressure increments,
    /// and return every physical displacement mode within the dense DOF cap.
    pub fn solve_vibration(
        &self,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.solve_vibration_impl(None, density, settings, cx)
    }

    /// Re-equilibrate an initial displacement/pressure iterate before forming
    /// the pencil. Original absolute support targets and load stepping apply;
    /// no foreign accepted state or unstable-KKT bypass is admitted.
    pub fn solve_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.solve_vibration_impl(Some((displacement_m, pressure_pa)), density, settings, cx)
    }

    fn solve_vibration_impl(
        &self,
        initial: Option<(&[[f64; 3]], &[f64])>,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        checkpoint(cx)?;
        if self.volumetric != VolumetricTreatment::OriginalMaterial || self.pressure_gauge.is_some()
        {
            return Err(invalid("mixed modal condensation requires original finite pressure compliance without a gauge").into());
        }
        if !tolerance(settings.spectrum.relative_residual_tolerance)
            || !tolerance(settings.pressure_residual_tolerance)
        {
            return Err(
                invalid("modal and pressure tolerances must be finite and in (0,1)").into(),
            );
        }
        let prepared = self.prepare(cx)?;
        let n = prepared.free_displacements;
        let np = prepared.mesh.pressure_nodes.len();
        let nu = 3 * prepared.mesh.nodes_count;
        if n == 0 || n > settings.spectrum.maximum_free_dofs || n > MAX_DENSE_FREE_DOFS {
            return Err(invalid(
                "mixed dense vibration requires 1..=128 bounded free displacement DOFs",
            )
            .into());
        }
        // Before any extra dense buffer: R and recovered pressures, conservative
        // dense eigensolver scratch/output, pressure vectors and sparse staging.
        let bytes = np
            .checked_mul(n)
            .and_then(|v| v.checked_mul(2))
            .and_then(|v| {
                n.checked_mul(n)
                    .and_then(|q| q.checked_mul(32))
                    .and_then(|q| v.checked_add(q))
            })
            .and_then(|v| np.checked_mul(16).and_then(|q| v.checked_add(q)))
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|v| {
                self.elements
                    .len()
                    .checked_mul(16 * 128)
                    .and_then(|q| v.checked_add(q))
            })
            .ok_or_else(|| invalid("mixed modal condensation storage overflow"))?;
        if bytes > settings.maximum_condensation_bytes {
            return Err(invalid(format!(
                "mixed modal condensation needs at most {bytes} payload bytes; cap {}",
                settings.maximum_condensation_bytes
            ))
            .into());
        }
        let free_dofs = prepared.free[..n].to_vec();
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
        let mut coo = Coo::new(np, np);
        for row in 0..np {
            checkpoint(cx)?;
            let (columns, values) = full.row(nu + row);
            for (&col, &value) in columns.iter().zip(values) {
                if col >= nu {
                    let transpose = full.get(col, nu + row);
                    if !value.is_finite()
                        || !transpose.is_finite()
                        || (value - transpose).abs()
                            > 64.0 * f64::EPSILON * np as f64 * value.abs().max(transpose.abs())
                    {
                        return Err(invalid(
                            "pressure-compliance block must be finite and symmetric",
                        )
                        .into());
                    }
                    coo.push(row, col - nu, -value);
                }
            }
        }
        let h = coo.assemble();
        drop(coo);
        checkpoint(cx)?;
        let symbolic =
            SymbolicLdlt::analyze(&h, DirectOrdering::Amd).map_err(MixedUpError::Direct)?;
        checkpoint(cx)?;
        let factor = symbolic
            .factor(
                &h,
                &LdltOptions {
                    max_factor_bytes: Some(self.settings.maximum_factor_bytes),
                    ..LdltOptions::default()
                },
            )
            .map_err(MixedUpError::Direct)?;
        if factor.inertia()
            != (Inertia {
                positive: np,
                negative: 0,
            })
        {
            return Err(MixedVibrationError::PressureBlockInertia {
                got: factor.inertia(),
            });
        }
        let mut response = vec![0.0; np * n];
        for (column, &dof) in free_dofs.iter().enumerate() {
            checkpoint(cx)?;
            let mut rhs = Vec::with_capacity(np);
            for p in 0..np {
                let value = full.get(nu + p, dof);
                if !value.is_finite() || value != full.get(dof, nu + p) {
                    return Err(invalid(
                        "mixed pressure/displacement blocks must be finite transposes",
                    )
                    .into());
                }
                rhs.push(value);
            }
            let x = factor.solve(&rhs);
            let relative = pressure_residual(&h, &x, &rhs, None, cx)?;
            if !relative.is_finite() || relative > settings.pressure_residual_tolerance {
                return Err(MixedVibrationError::PressureResidual {
                    stage: "condensation",
                    index: column,
                    relative,
                });
            }
            for (p, value) in x.into_iter().enumerate() {
                response[p * n + column] = value;
            }
        }
        drop(factor);
        drop(symbolic);
        for (row, &dof) in free_dofs.iter().enumerate() {
            checkpoint(cx)?;
            let (columns, values) = full.row(dof);
            for (&col, &value) in columns.iter().zip(values) {
                if col >= nu {
                    for j in 0..n {
                        let entry = &mut stiffness_n_m[row * n + j];
                        *entry = value.mul_add(response[(col - nu) * n + j], *entry);
                    }
                }
            }
        }
        symmetric_dense(&stiffness_n_m, n)?;
        checkpoint(cx)?;
        let displacement_modes =
            eigh_gen_dense(&stiffness_n_m, &mass_kg, n).map_err(MixedVibrationError::Modal)?;
        checkpoint(cx)?;
        if displacement_modes.len() != n {
            return Err(invalid("incomplete mixed modal basis").into());
        }
        let scale = displacement_modes
            .iter()
            .map(|m| m.lambda.abs())
            .fold(0.0, f64::max);
        let mut modes = Vec::with_capacity(n);
        for (index, displacement) in displacement_modes.into_iter().enumerate() {
            checkpoint(cx)?;
            if displacement.phi.len() != n
                || displacement.phi.iter().any(|v| !v.is_finite())
                || !displacement.lambda.is_finite()
                || !displacement.residual.is_finite()
                || displacement.residual < 0.0
                || !displacement.interval.0.is_finite()
                || !displacement.interval.1.is_finite()
                || displacement.residual > settings.spectrum.relative_residual_tolerance * scale
            {
                return Err(MixedVibrationError::ModeResidual {
                    index,
                    eigenvalue_per_s2: displacement.lambda,
                    residual_per_s2: displacement.residual,
                });
            }
            let mut pressure_shape = Vec::with_capacity(np);
            let mut rhs = Vec::with_capacity(np);
            let mut rhs_magnitudes = Vec::with_capacity(np);
            for p in 0..np {
                checkpoint(cx)?;
                pressure_shape.push(
                    (0..n)
                        .map(|j| response[p * n + j] * displacement.phi[j])
                        .sum(),
                );
                rhs.push(
                    (0..n)
                        .map(|j| full.get(nu + p, free_dofs[j]) * displacement.phi[j])
                        .sum(),
                );
                // Isochoric modes can cancel Kpu*phi almost exactly. Retain
                // the uncancelled block work in their backward-error scale.
                rhs_magnitudes.push(
                    (0..n)
                        .map(|j| (full.get(nu + p, free_dofs[j]) * displacement.phi[j]).abs())
                        .sum(),
                );
            }
            let relative = pressure_residual(&h, &pressure_shape, &rhs, Some(&rhs_magnitudes), cx)?;
            if !relative.is_finite() || relative > settings.pressure_residual_tolerance {
                return Err(MixedVibrationError::PressureResidual {
                    stage: "mode",
                    index,
                    relative,
                });
            }
            modes.push(MixedVibrationMode {
                displacement,
                pressure_shape,
                pressure_equation_relative_residual: relative,
            });
        }
        let result = MixedUpVibration {
            problem: *self,
            settings,
            reference_density_kg_m3: density.to_vec(),
            equilibrium,
            free_dofs,
            mixed_tangent: full,
            pressure_response_pa_per_m: response,
            stiffness_n_m,
            mass_kg,
            modes,
        };
        checkpoint(cx)?;
        Ok(result)
    }
}

impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedUpVibration<'_, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'_, L> {
        &self.problem
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
    /// Original physical mixed Hessian, including constrained displacement rows.
    /// All original displacement DOFs precede pressure_nodes in the equilibrium.
    pub const fn mixed_tangent(&self) -> &Csr {
        &self.mixed_tangent
    }
    /// Row-major pressure-node by free-displacement response dp=R du [Pa/m].
    pub fn pressure_response_pa_per_m(&self) -> &[f64] {
        &self.pressure_response_pa_per_m
    }
    /// Complete row-major condensed free-displacement stiffness [N/m].
    pub fn stiffness_n_m(&self) -> &[f64] {
        &self.stiffness_n_m
    }
    /// Complete row-major consistent reference displacement mass [kg].
    pub fn mass_kg(&self) -> &[f64] {
        &self.mass_kg
    }
    pub fn modes(&self) -> &[MixedVibrationMode] {
        &self.modes
    }
}
