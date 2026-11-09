//! Work-conjugate local mechanics on original affine lock/gear manifolds.
//!
//! For the admitted physical lift dq=N*z, the same native RNEA/load derivatives
//! supply M_r=N^T*M*N, C_r=N^T*C*N and K_r=N^T*K*N. Original constant ideal
//! reactions do no tangent work. Neither multiplier coordinates, penalty masses
//! nor fitted damping enter the existing fs-modal pencil. Affine rows have no
//! curvature term; configuration-dependent constraints need their own owner.

/// Physical forcing and original ideal reactions on this admitted local pencil.
pub mod forced;

use super::ConstrainedArticulatedEquilibrium;
use super::super::{EquilibriumBushings, SpatialMotionLoads};
use super::super::linearized::{
    self as native, ArticulatedLinearization, CentralDifferenceControls,
    JointTangentAction, LinearizationControls, LinearizationError,
};
use crate::articulated::ArticulatedModel;
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, PreparedJointConstraints,
    admit_constraints, apply_transpose_rows,
};
use fs_la::subspace::{LinearSubspace, SubspaceControls, SubspaceError};
use fs_math::{c64::C64, dd::Dd};

/// Explicit bounds for original row admission and the represented physical lift.
#[derive(Clone, Copy, Debug)]
pub struct ConstrainedLinearizationControls {
    /// Bounded dense QR and original row mobility, in 1..=64.
    pub maximum_physical_dofs: usize,
    /// Original normalized row-mobility floor, in [4096*EPSILON,1].
    pub minimum_reciprocal_condition: f64,
    /// Absolute diagonal floor for original QR of row-normalized (A*S)^T.
    /// This is an additional numerical gate, not an alternative rank decision.
    pub minimum_scaled_qr_diagonal: f64,
    /// Maximum normalized original A*N entry.
    pub maximum_constraint_residual: f64,
    /// Maximum absolute entry of (S^-1*N)^T*(S^-1*N)-I.
    pub maximum_orthogonality_error: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ConstrainedLinearizationError {
    Original(LinearizationError),
    Constraint(ConstraintError),
    Subspace(SubspaceError),
    Invalid(&'static str),
}
impl From<LinearizationError> for ConstrainedLinearizationError {
    fn from(value: LinearizationError) -> Self { Self::Original(value) }
}
impl From<ConstraintError> for ConstrainedLinearizationError {
    fn from(value: ConstraintError) -> Self { Self::Constraint(value) }
}
impl From<SubspaceError> for ConstrainedLinearizationError {
    fn from(value: SubspaceError) -> Self { Self::Subspace(value) }
}
impl core::fmt::Display for ConstrainedLinearizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constrained articulated response: {self:?}")
    }
}
impl std::error::Error for ConstrainedLinearizationError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConstrainedLinearizationError> {
    if cancelled() { Err(ConstraintError::Cancelled.into()) } else { Ok(()) }
}

fn represented(value: Dd) -> Result<f64, ConstrainedLinearizationError> {
    let result = value.to_f64();
    if !result.is_finite() || (value != Dd::ZERO && result == 0.0) {
        Err(ConstrainedLinearizationError::Invalid("physical subspace arithmetic is not representable"))
    } else {
        Ok(result)
    }
}

fn product(a: f64, b: f64) -> Result<Dd, ConstrainedLinearizationError> {
    let value = Dd::from_f64(a)*Dd::from_f64(b);
    represented(value)?;
    if a != 0.0 && b != 0.0 && value == Dd::ZERO {
        return Err(ConstrainedLinearizationError::Invalid("physical subspace product underflow"));
    }
    Ok(value)
}

/// Rebind a public equilibrium result to its original rows, loads and materials.
///
/// Cached report fields cannot admit the root. Original row position/rate
/// tolerances, native independent-row mobility and RNEA force balance against
/// the actual original A^T*lambda all run again. The unconstrained public
/// response entry retains its original no-reaction balance requirement.
///
/// Rows are fixed affine relations. Mean reactions are constant in their
/// original joint coordinates; subtracting that constant from native residual
/// probes reduces cancellation without changing their derivatives. Perturbed
/// ideal reactions vanish under N^T. Only allowed directions are differentiated.
///
/// The original Householder QR chooses n-m columns after original admission
/// proves m independent rows. The shared subspace owner independently gates
/// the represented physical A*N and scaled column orthogonality.
#[allow(clippy::too_many_arguments)]
pub fn prepare<'a>(
    model: &'a ArticulatedModel,
    equilibrium: &ConstrainedArticulatedEquilibrium,
    loads: SpatialMotionLoads<'a>,
    bushings: EquilibriumBushings<'a>,
    rows: &[JointConstraint],
    controls: LinearizationControls<'_>,
    constrained: ConstrainedLinearizationControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedLinearization<'a>, ConstrainedLinearizationError> {
    poll(&mut cancelled)?;
    let n = model.dof_count();
    let m = rows.len();
    if !(1..=64).contains(&constrained.maximum_physical_dofs)
        || n > constrained.maximum_physical_dofs || m > n
        || !constrained.minimum_scaled_qr_diagonal.is_finite()
        || !(4096.0*f64::EPSILON..=1.0).contains(&constrained.minimum_scaled_qr_diagonal)
    {
        return Err(ConstrainedLinearizationError::Invalid("explicit bounded row/nullspace controls required"));
    }
    let admitted_rows = admit_constraints(rows, &equilibrium.state.joint_position,
        &equilibrium.state.joint_velocity, &mut cancelled)?;
    let mobility = PreparedJointConstraints::prescribed_trial(model,
        &equilibrium.state.joint_position, rows,
        constrained.minimum_reciprocal_condition, &mut cancelled)?;
    let reciprocal_condition = mobility.reciprocal_condition;
    let reaction = apply_transpose_rows(rows, &equilibrium.multipliers, n, &mut cancelled)?;
    let original = native::prepare_with_reaction(model, &equilibrium.state,
        loads, bushings, controls, Some(&reaction), &mut cancelled)?;
    let reduced = n-m;
    let mut basis = vec![0.0; n*reduced];
    if m == 0 {
        for i in 0..n { basis[i*n+i] = controls.coordinate_scales[i]; }
    } else if reduced != 0 {
        // Normalize each row before the unchanged native Householder owner.
        // The actual original row is retained for the final A*N admission.
        let mut transpose = vec![0.0; n*m];
        for (j, row) in rows.iter().enumerate() {
            poll(&mut cancelled)?;
            let mut scale = 0.0_f64;
            for (&a, &s) in row.coefficients.iter().zip(controls.coordinate_scales) {
                scale = scale.max(represented(product(a, s)?)?.abs());
            }
            if scale == 0.0 {
                return Err(ConstrainedLinearizationError::Invalid("scaled original row is zero"));
            }
            for i in 0..n {
                transpose[i*m+j] = represented(
                    product(row.coefficients[i], controls.coordinate_scales[i])?
                    / Dd::from_f64(scale))?;
                if row.coefficients[i] != 0.0 && transpose[i*m+j] == 0.0 {
                    return Err(ConstrainedLinearizationError::Invalid("normalized row entry underflow"));
                }
            }
        }
        poll(&mut cancelled)?;
        let qr = fs_la::factor::qr(&transpose, n, m);
        poll(&mut cancelled)?;
        for i in 0..m {
            let diagonal = qr.r(i, i).abs();
            if !diagonal.is_finite() || diagonal < constrained.minimum_scaled_qr_diagonal {
                return Err(ConstrainedLinearizationError::Invalid("original scaled row QR is unresolved"));
            }
        }
        let mut column = vec![0.0; n];
        for j in 0..reduced {
            poll(&mut cancelled)?;
            column.fill(0.0);
            column[m+j] = 1.0;
            qr.apply_q(&mut column);
            for i in 0..n {
                basis[i*reduced+j] = represented(
                    product(controls.coordinate_scales[i], column[i])?)?;
            }
        }
    }
    let original_rows: Vec<_> = rows.iter().flat_map(|row| row.coefficients.iter().copied()).collect();
    let subspace = LinearSubspace::try_new(controls.coordinate_scales, reduced,
        basis, &original_rows, SubspaceControls {
            maximum_dimension: constrained.maximum_physical_dofs,
            maximum_constraint_residual: constrained.maximum_constraint_residual,
            maximum_orthogonality_error: constrained.maximum_orthogonality_error,
        }, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(ConstrainedArticulatedLinearization {
        original, subspace, mean_reaction: reaction, reciprocal_condition,
        difference: controls.difference,
        rows: rows.to_vec(), coordinate_scales: controls.coordinate_scales.to_vec(),
        mean_position_residuals: admitted_rows.position,
        minimum_scaled_qr_diagonal: constrained.minimum_scaled_qr_diagonal,
    })
}

/// Reduced effort and the unchanged native physical action used to compute it.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstrainedTangentAction {
    pub generalized_effort: Vec<f64>,
    /// |N|^T times original successive-stencil estimates; not error enclosures.
    pub disagreement: Vec<f64>,
    pub original: JointTangentAction,
}

/// Immutable original source binding plus its admitted physical tangent lift.
pub struct ConstrainedArticulatedLinearization<'a> {
    original: ArticulatedLinearization<'a>,
    subspace: LinearSubspace,
    mean_reaction: Vec<f64>,
    reciprocal_condition: f64,
    difference: CentralDifferenceControls,
    rows: Vec<JointConstraint>,
    coordinate_scales: Vec<f64>,
    mean_position_residuals: Vec<f64>,
    minimum_scaled_qr_diagonal: f64,
}

impl ConstrainedArticulatedLinearization<'_> {
    pub fn state(&self) -> &PrescribedBaseMotionState { self.original.state() }
    pub fn subspace(&self) -> &LinearSubspace { &self.subspace }
    pub fn mean_reaction_effort(&self) -> &[f64] { &self.mean_reaction }
    pub fn reciprocal_condition(&self) -> f64 { self.reciprocal_condition }
    pub fn scaled_equilibrium_residual(&self) -> f64 {
        self.original.scaled_equilibrium_residual()
    }

    /// Physical q/velocity/acceleration perturbation for a reduced direction.
    pub fn lift_perturbation(&self, reduced: &[f64], cancelled: impl FnMut()->bool)
        -> Result<Vec<f64>, ConstrainedLinearizationError> {
        Ok(self.subspace.lift(reduced, cancelled)?)
    }

    /// Work-conjugate drive projection. Ideal original reactions annihilate
    /// to the independently admitted original A*N residual, not by assumption.
    pub fn project_effort(&self, physical: &[f64], cancelled: impl FnMut()->bool)
        -> Result<Vec<f64>, ConstrainedLinearizationError> {
        Ok(self.subspace.project(physical, cancelled)?)
    }

    pub fn inertia_action(&self, acceleration: &[f64], mut cancelled: impl FnMut()->bool)
        -> Result<Vec<f64>, ConstrainedLinearizationError> {
        let physical = self.subspace.lift(acceleration, &mut cancelled)?;
        if self.subspace.dimension() == 0 {
            poll(&mut cancelled)?;
            return Ok(Vec::new());
        }
        let effort = self.original.inertia_action(&physical, &mut cancelled)?;
        Ok(self.subspace.project(&effort, cancelled)?)
    }

    pub fn stiffness_action(&self, displacement: &[f64], cancelled: impl FnMut()->bool)
        -> Result<ConstrainedTangentAction, ConstrainedLinearizationError> {
        self.action(displacement, false, self.difference, cancelled)
    }

    pub fn velocity_action(&self, velocity: &[f64], cancelled: impl FnMut()->bool)
        -> Result<ConstrainedTangentAction, ConstrainedLinearizationError> {
        self.action(velocity, true, self.difference, cancelled)
    }

    fn action(&self, direction: &[f64], rate: bool, options: CentralDifferenceControls,
        mut cancelled: impl FnMut()->bool)
        -> Result<ConstrainedTangentAction, ConstrainedLinearizationError> {
        let physical = self.subspace.lift(direction, &mut cancelled)?;
        if self.subspace.dimension() == 0 {
            poll(&mut cancelled)?;
            return Ok(ConstrainedTangentAction {
                generalized_effort: Vec::new(), disagreement: Vec::new(),
                original: JointTangentAction {
                    effort: vec![0.0; self.subspace.physical_dimension()],
                    disagreement: vec![0.0; self.subspace.physical_dimension()],
                    residual_evaluations: 0, refinements: 0, final_scaled_step: 0.0,
                    maximum_stencil_distortion: 0.0,
                },
            });
        }
        let original = self.original.action(&physical, rate, options, &mut cancelled)?;
        let generalized_effort = self.subspace.project(&original.effort, &mut cancelled)?;
        let n = self.subspace.physical_dimension();
        let r = self.subspace.dimension();
        let mut disagreement = vec![0.0; r];
        for j in 0..r {
            poll(&mut cancelled)?;
            let mut sum = Dd::ZERO;
            for i in 0..n {
                sum = sum + product(self.subspace.basis()[i*r+j].abs(),
                    original.disagreement[i])?;
            }
            disagreement[j] = represented(sum)?;
        }
        poll(&mut cancelled)?;
        Ok(ConstrainedTangentAction { generalized_effort, disagreement, original })
    }

    /// Export the same allowed-direction operators to the existing modal owner.
    /// One native mass extraction and every original residual query are counted.
    /// K/C retain their complete nonsymmetric entries and signed stiffness.
    /// A fully locked model is a legitimate empty pencil and uses no probes.
    pub fn quadratic_model(&self, maximum_reduced_dofs: usize, maximum_evaluations: usize,
        mut cancelled: impl FnMut()->bool)
        -> Result<ConstrainedQuadraticModel, ConstrainedLinearizationError> {
        poll(&mut cancelled)?;
        let r = self.subspace.dimension();
        if maximum_reduced_dofs > 64 || r > maximum_reduced_dofs
            || maximum_evaluations > 16384 || (r != 0 && maximum_evaluations == 0)
        {
            return Err(ConstrainedLinearizationError::Invalid("explicit reduced dimension and source-work bounds required"));
        }
        let mut result = ConstrainedQuadraticModel {
            subspace: self.subspace.clone(), mass: Vec::new(),
            stiffness: vec![0.0; r*r], velocity: vec![0.0; r*r],
            stiffness_disagreement: vec![0.0; r*r],
            velocity_disagreement: vec![0.0; r*r],
            source_evaluations: 0, maximum_stencil_distortion: 0.0,
            original_mass: Vec::new(),
            original_stiffness_lift: vec![0.0; self.subspace.physical_dimension()*r],
            original_velocity_lift: vec![0.0; self.subspace.physical_dimension()*r],
            rows: self.rows.clone(), coordinate_scales: self.coordinate_scales.clone(),
            mean_position_residuals: self.mean_position_residuals.clone(),
            minimum_scaled_qr_diagonal: self.minimum_scaled_qr_diagonal,
        };
        if r == 0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        let mass = self.original.mass_matrix(&mut cancelled)?;
        result.mass = self.subspace.symmetric_congruence(&mass, &mut cancelled)?;
        result.original_mass = mass;
        result.source_evaluations = 1;
        let mut direction = vec![0.0; r];
        for column in 0..r {
            direction.fill(0.0);
            direction[column] = 1.0;
            for rate in [false, true] {
                poll(&mut cancelled)?;
                let remaining = maximum_evaluations.saturating_sub(result.source_evaluations);
                if remaining < 4 {
                    return Err(LinearizationError::Budget {
                        used: result.source_evaluations, maximum: maximum_evaluations,
                    }.into());
                }
                let mut options = self.difference;
                options.maximum_evaluations = options.maximum_evaluations.min(remaining);
                let value = self.action(&direction, rate, options, &mut cancelled)?;
                result.source_evaluations += value.original.residual_evaluations;
                result.maximum_stencil_distortion = result.maximum_stencil_distortion
                    .max(value.original.maximum_stencil_distortion);
                // Retain the same native allowed-direction force columns for
                // complete physical reaction recovery; never probe locked DOFs.
                let physical = if rate { &mut result.original_velocity_lift }
                    else { &mut result.original_stiffness_lift };
                for (row, &effort) in value.original.effort.iter().enumerate() {
                    physical[row*r+column] = effort;
                }
                let (matrix, errors) = if rate {
                    (&mut result.velocity, &mut result.velocity_disagreement)
                } else {
                    (&mut result.stiffness, &mut result.stiffness_disagreement)
                };
                for row in 0..r {
                    matrix[row*r+column] = value.generalized_effort[row];
                    errors[row*r+column] = value.disagreement[row];
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

/// Immutable allowed-coordinate pencil. Reduced coordinates z are dimensionless,
/// physical joint perturbations are N*z, and generalized efforts are N^T*f.
/// C is not assumed reciprocal/dissipative, K is not assumed positive, and the
/// reported finite-difference disagreements are not stability certificates.
pub struct ConstrainedQuadraticModel {
    subspace: LinearSubspace,
    mass: Vec<f64>,
    velocity: Vec<f64>,
    stiffness: Vec<f64>,
    stiffness_disagreement: Vec<f64>,
    velocity_disagreement: Vec<f64>,
    source_evaluations: usize,
    maximum_stencil_distortion: f64,
    original_mass: Vec<f64>,
    original_stiffness_lift: Vec<f64>,
    original_velocity_lift: Vec<f64>,
    rows: Vec<JointConstraint>,
    coordinate_scales: Vec<f64>,
    mean_position_residuals: Vec<f64>,
    minimum_scaled_qr_diagonal: f64,
}
impl ConstrainedQuadraticModel {
    pub fn dimension(&self) -> usize { self.subspace.dimension() }
    pub fn subspace(&self) -> &LinearSubspace { &self.subspace }
    pub fn mass(&self) -> &[f64] { &self.mass }
    pub fn velocity_operator(&self) -> &[f64] { &self.velocity }
    pub fn stiffness(&self) -> &[f64] { &self.stiffness }
    pub fn stiffness_disagreement(&self) -> &[f64] { &self.stiffness_disagreement }
    pub fn velocity_disagreement(&self) -> &[f64] { &self.velocity_disagreement }
    pub fn source_evaluations(&self) -> usize { self.source_evaluations }
    pub fn maximum_stencil_distortion(&self) -> f64 { self.maximum_stencil_distortion }

    /// Existing complete quadratic poles, including unstable signed growth.
    /// The original eigensolver admits a zero-dimensional empty pole set.
    pub fn poles(&self, maximum_sweeps_per_window: usize, cancelled: impl FnMut()->bool)
        -> Result<Vec<C64>, fs_modal::ModalError> {
        fs_modal::quadratic_eigenvalues_with_cancel(&self.mass, &self.velocity,
            &self.stiffness, self.dimension(), maximum_sweeps_per_window, cancelled)
    }

    /// Existing local modes/residues on this exact M/C/K. Fully locked systems
    /// have no local mode and return None rather than a fictitious zero mass.
    pub fn prepare_modes(&self, options: fs_modal::quadratic_modes::PoleOptions,
        mut cancelled: impl FnMut()->bool)
        -> Result<Option<fs_modal::quadratic_modes::PreparedQuadraticModes<'_>>,
            fs_modal::quadratic_modes::QuadraticModeError> {
        if cancelled() {
            return Err(fs_modal::quadratic_modes::QuadraticModeError::Cancelled);
        }
        if self.dimension() == 0 {
            return Ok(None);
        }
        Ok(Some(fs_modal::quadratic_modes::PreparedQuadraticModes::new(
            &self.mass, &self.velocity, &self.stiffness, self.dimension(), options, cancelled)?))
    }

    /// Expand a reduced complex displacement or velocity shape by the same
    /// original physical lift. No finite-amplitude trajectory is inferred.
    pub fn expand_mode(&self, reduced: &[C64], mut cancelled: impl FnMut()->bool)
        -> Result<Vec<C64>, ConstrainedLinearizationError> {
        let real: Vec<_> = reduced.iter().map(|v| v.re).collect();
        let imaginary: Vec<_> = reduced.iter().map(|v| v.im).collect();
        let real = self.subspace.lift(&real, &mut cancelled)?;
        let imaginary = self.subspace.lift(&imaginary, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(real.into_iter().zip(imaginary).map(|(re, im)| C64::new(re, im)).collect())
    }

    pub fn project_effort(&self, physical: &[f64], cancelled: impl FnMut()->bool)
        -> Result<Vec<f64>, ConstrainedLinearizationError> {
        Ok(self.subspace.project(physical, cancelled)?)
    }
}
