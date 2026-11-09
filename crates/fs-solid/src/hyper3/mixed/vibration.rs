//! Physical vibration of the original finite-compliance Taylor–Hood solid.
//!
//! The accepted mixed equilibrium owns raw physical tangent blocks. Pressure
//! has no inertia: H=-Kpp, R=H^-1 Kpu, dp=R du and K_eff=Kuu+Kup R. A single
//! native sparse factor supplies all pressure solves; fs-modal owns the bounded
//! complete dense displacement spectrum. No pressure penalty, bulk replacement,
//! lumped quadratic mass or additional initial-stress matrix is introduced.
//!
//! The condensation methods support OriginalMaterial on its admitted positive
//! bulk-curvature branch. The separate [incompressible] route supplies physical
//! modes of the original exact constraint through its displacement nullspace.
//! The original static owner refuses unstable/singular KKT equilibria; neither
//! adapter bypasses that gate. Returned eigenvalues remain signed, with their
//! residual intervals, without frequency clipping.
//!
//! Constitutive/source conditions and dead loads are held fixed. The pressure
//! adapter also retains its held boundary load and exact q Hessian(V). The original
//! P2/P1 weak pressure space, geometry and quadrature remain authoritative;
//! curved geometry and nonlinear coefficients still require quadrature and
//! mesh convergence. This route is not an inf-sup or continuum-stability proof.

/// Full forced response with pressure recovery and original support reactions.
pub mod harmonic;
pub mod incompressible;

use super::gravity::HeldMixedGravity;
use super::pressure::UniformBoundaryPressure;
use super::{
    MixedUpError, MixedUpProblem, MixedUpSolution, PreparedMixedUpProblem, checkpoint,
    inertia_error, invalid,
};
use crate::hyper3::vibration::{MAX_DENSE_FREE_DOFS, VibrationSettings};
use crate::inertia::ReferenceMass;
use crate::linear3::tet10::TET10_EDGES;
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::{Coo, Csr, DirectOrdering, Inertia, LdltFactor, LdltOptions, SymbolicLdlt};

/// Explicit modal accuracy and extra storage admission. Original mixed mesh,
/// quadrature, Newton, assembly and sparse-factor budgets remain in force.
#[derive(Clone, Copy, Debug)]
pub struct MixedVibrationSettings {
    pub spectrum: VibrationSettings,
    /// Native full signed P2 reference-mass assembly and retained operator
    /// payload, separate from KKT and dense modal storage.
    pub maximum_mass_assembly_bytes: usize,
    /// Conservative pressure reduction (finite-compliance condensation or
    /// exact-constraint nullspace), dense work and retained modal payload.
    /// Excludes the original mixed assembly and retained native sparse factor,
    /// which retain their own budgets; this is not an OS memory/RSS limit.
    pub maximum_condensation_bytes: usize,
    /// Relative backward-error gate for pressure elimination and recovered
    /// pressure equations, or all original exact-constraint rows. Finite and
    /// strictly between zero and one. Exact-constrained momentum uses the
    /// spectrum relative tolerance as its separate backward-error gate.
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
    boundary_pressure: Option<UniformBoundaryPressure>,
    gravity: Option<HeldMixedGravity<'a>>,
    settings: MixedVibrationSettings,
    reference_mass: ModalReferenceMass,
    equilibrium: MixedUpSolution,
    free_dofs: Vec<usize>,
    mixed_tangent: Csr,
    // Retain the admitted H=-Kpp factor for inhomogeneous pressure forcing.
    // Its storage remains under problem.settings.maximum_factor_bytes.
    pressure_factor: LdltFactor,
    pressure_response_pa_per_m: Vec<f64>,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    modes: Vec<MixedVibrationMode>,
}

#[derive(Clone, Copy)]
enum ReferenceDensity<'a> {
    PerElement(&'a [f64]),
    Corners(&'a [[f64; 4]]),
}

/// The complete original operator remains available after support reduction.
/// Cell means describe inventory; they never replace a graded mass or load.
#[derive(Debug)]
struct ModalReferenceMass {
    operator: ReferenceMass,
    means_kg_m3: Vec<f64>,
    corners_kg_m3: Option<Vec<[f64; 4]>>,
}

impl ModalReferenceMass {
    fn assemble<L: HyperelasticMaterial + MixedFiniteStrainPotential>(
        prepared: &PreparedMixedUpProblem<'_, L>,
        density: ReferenceDensity<'_>,
        maximum_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedUpError> {
        checkpoint(cx)?;
        let elements = &prepared.mesh.elements;
        let ndofs = 3 * prepared.mesh.nodes_count;
        let ReferenceDensity::Corners(corners) = density else {
            let ReferenceDensity::PerElement(values) = density else {
                unreachable!()
            };
            let operator = ReferenceMass::from_tet10_quadrature(
                elements.iter().map(|element| {
                    (
                        element.connectivity,
                        element.points.iter().map(|q| (q.shapes, q.weight_m3)),
                    )
                }),
                values,
                ndofs,
                maximum_bytes,
                cx,
            )
            .map_err(inertia_error)?;
            return Ok(Self {
                operator,
                means_kg_m3: values.to_vec(),
                corners_kg_m3: None,
            });
        };
        if corners.len() != elements.len() {
            return Err(invalid(
                "one reference-density corner profile per Tet10 required",
            ));
        }
        for (element, values) in elements.iter().zip(corners) {
            checkpoint(cx)?;
            if values.iter().any(|rho| !rho.is_finite() || *rho <= 0.0) {
                return Err(invalid(
                    "reference corner densities must be finite and positive",
                ));
            }
            if values.iter().all(|rho| *rho == values[0]) {
                continue;
            }
            let nodes = element.connectivity;
            let affine = TET10_EDGES.iter().enumerate().all(|(edge, &[a, b])| {
                (0..3).all(|c| {
                    prepared.problem.nodes_m[nodes[edge + 4]][c]
                        == f64::midpoint(
                            prepared.problem.nodes_m[nodes[a]][c],
                            prepared.problem.nodes_m[nodes[b]][c],
                        )
                })
            });
            if !affine && prepared.problem.quadrature_points_per_dir < 6 {
                return Err(invalid(
                    "curved Tet10 graded reference density requires at least six Duffy points",
                ));
            }
        }
        let sample = |values: [f64; 4], lambda: [f64; 4]| {
            if values.iter().all(|rho| *rho == values[0]) {
                values[0]
            } else {
                (0..4).map(|a| values[a] * lambda[a]).sum()
            }
        };
        let operator = ReferenceMass::from_tet10_density_quadrature(
            elements.iter().zip(corners).map(|(element, &values)| {
                (
                    element.connectivity,
                    element
                        .points
                        .iter()
                        .map(move |q| (q.shapes, q.weight_m3, sample(values, q.pressure_shapes))),
                )
            }),
            ndofs,
            maximum_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        let mut means = Vec::with_capacity(corners.len());
        for (element, &values) in elements.iter().zip(corners) {
            checkpoint(cx)?;
            let mean = if values.iter().all(|rho| *rho == values[0]) {
                values[0]
            } else {
                element
                    .points
                    .iter()
                    .map(|q| {
                        sample(values, q.pressure_shapes)
                            * (q.weight_m3 / element.reference_volume_m3)
                    })
                    .sum::<f64>()
            };
            if !mean.is_finite() || mean <= 0.0 {
                return Err(invalid("reference-volume mean density is unrepresentable"));
            }
            means.push(mean);
        }
        checkpoint(cx)?;
        Ok(Self {
            operator,
            means_kg_m3: means,
            corners_kg_m3: Some(corners.to_vec()),
        })
    }

    fn acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, MixedUpError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.operator.ndofs() / 3
            || acceleration_m_s2.iter().flatten().any(|a| !a.is_finite())
        {
            return Err(invalid(
                "one finite acceleration vector per original node required",
            ));
        }
        let forces = self
            .operator
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)?;
        let mut loads = Vec::with_capacity(forces.len());
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if force != 0.0 {
                loads.push((dof, force));
            }
        }
        checkpoint(cx)?;
        Ok(loads)
    }
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

// Conservative source qualification precedes canonical triangle storage.
// Schur products can differ across triangles by accumulation roundoff; a
// materially asymmetric source must still refuse, never be silently averaged.
fn canonical_symmetric_dense(matrix: &mut [f64], n: usize) -> Result<(), MixedUpError> {
    symmetric_dense(matrix, n)?;
    for i in 0..n {
        for j in 0..i {
            matrix[i * n + j] = matrix[j * n + i];
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
        self.prepare(cx)?
            .solve_vibration_impl(None, density, settings, cx)
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
        self.prepare(cx)?.solve_vibration_impl(
            Some((displacement_m, pressure_pa)),
            density,
            settings,
            cx,
        )
    }

    /// Equilibrate with positive P1 reference density at each Tet10's four
    /// corners. The signed P2 inertia retains the complete original profile;
    /// density is not averaged or updated on the deformed configuration.
    /// Graded curved maps require at least six Duffy points per direction.
    pub fn solve_vibration_with_reference_density(
        &self,
        corner_density_kg_m3: &[[f64; 4]],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepare(cx)?.solve_vibration_density_impl(
            None,
            ReferenceDensity::Corners(corner_density_kg_m3),
            settings,
            cx,
        )
    }

    /// Re-equilibrate the original iterate with the same graded reference
    /// inventory used for the full mass and acceleration loading.
    pub fn solve_vibration_from_with_reference_density(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        corner_density_kg_m3: &[[f64; 4]],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepare(cx)?.solve_vibration_density_impl(
            Some((displacement_m, pressure_pa)),
            ReferenceDensity::Corners(corner_density_kg_m3),
            settings,
            cx,
        )
    }
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> PreparedMixedUpProblem<'a, L> {
    pub(super) fn solve_vibration_impl(
        &self,
        initial: Option<(&[[f64; 3]], &[f64])>,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.solve_vibration_density_impl(
            initial,
            ReferenceDensity::PerElement(density),
            settings,
            cx,
        )
    }

    fn solve_vibration_density_impl(
        &self,
        initial: Option<(&[[f64; 3]], &[f64])>,
        density: ReferenceDensity<'_>,
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        checkpoint(cx)?;
        let prepared = self;
        let problem = &self.problem;
        if problem.volumetric != VolumetricTreatment::OriginalMaterial
            || problem.pressure_gauge.is_some()
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
                problem
                    .elements
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
        let reference_mass = ModalReferenceMass::assemble(
            prepared,
            density,
            settings.maximum_mass_assembly_bytes,
            cx,
        )?;
        let mut mass_kg = Vec::with_capacity(n * n);
        for &row in &free_dofs {
            checkpoint(cx)?;
            for &col in &free_dofs {
                mass_kg.push(reference_mass.operator.matrix().get(row, col));
            }
        }
        canonical_symmetric_dense(&mut mass_kg, n)?;
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
                    max_factor_bytes: Some(problem.settings.maximum_factor_bytes),
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
        canonical_symmetric_dense(&mut stiffness_n_m, n)?;
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
            problem: *problem,
            boundary_pressure: prepared.boundary_pressure,
            gravity: prepared.gravity,
            settings,
            reference_mass,
            equilibrium,
            free_dofs,
            mixed_tangent: full,
            pressure_factor: factor,
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
    /// A held uniform load on the entire material-domain boundary. Its q V
    /// loading potential is separate from equilibrium.potential_j, which keeps
    /// the original material mixed potential. The retained tangent includes
    /// its follower Hessian exactly once; no pressure inertia is introduced.
    pub const fn boundary_pressure(&self) -> Option<UniformBoundaryPressure> {
        self.boundary_pressure
    }
    /// Original stationary field and reference density held during these
    /// modes. The retained tangent includes minus its external force gradient.
    pub const fn gravity(&self) -> Option<HeldMixedGravity<'_>> {
        self.gravity
    }
    pub const fn settings(&self) -> MixedVibrationSettings {
        self.settings
    }
    /// Reference-volume mean densities for reporting. The original full field
    /// supplies graded inertia and acceleration loading.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_mass.means_kg_m3
    }
    pub fn reference_corner_density_kg_m3(&self) -> Option<&[[f64; 4]]> {
        self.reference_mass.corners_kg_m3.as_deref()
    }
    /// Total original mass before support or pressure reduction [kg].
    pub fn reference_mass_kg(&self) -> f64 {
        self.reference_mass.operator.reference_mass_kg()
    }
    /// Project P2 nodal acceleration through the SAME full reference mass.
    /// Returned original displacement forces [N] include support coordinates
    /// and signed corner weights. They have no pressure inertia. Select the
    /// free entries for modal forcing; do not add this load a second time as
    /// body-force density. Evolving force laws keep their own work convention.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, MixedVibrationError> {
        Ok(self
            .reference_mass
            .acceleration_load_n(acceleration_m_s2, cx)?)
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
