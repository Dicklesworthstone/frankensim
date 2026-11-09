//! Original small-signal mechanics about loaded or uniformly rotating equilibria.
//!
//! M*dqdd + C*dqd + K*dq = d_effort, in the original joint coordinates.
//! RNEA supplies M and the SAME equilibrium residual supplies numerical
//! position/rate actions. C retains material damping and gyroscopic terms;
//! K retains preload geometry, follower loads and centrifugal softening.
//! No symmetry, positive stiffness, stability or derivative certificate is inferred.

use super::*;
use crate::articulated::{BaseState, FreeFloatingBaseState};
use crate::articulated::generalized::{GeneralizedError, free_floating_terms};
use fs_math::{c64::C64, dd::Dd};
pub use fs_solver::finite_difference::central::CentralDifferenceControls;
use fs_solver::finite_difference::central::{self, CentralDifference, CentralDifferenceError};

/// Explicit physical chart scales and original equilibrium admission.
#[derive(Clone, Copy, Debug)]
pub struct LinearizationControls<'a> {
    /// Positive joint position scales, metres or radians by joint kind.
    pub coordinate_scales: &'a [f64],
    /// Positive joint rate scales, metres/second or radians/second.
    pub velocity_scales: &'a [f64],
    /// Positive conjugate effort scales, newtons or newton-metres.
    pub effort_scales: &'a [f64],
    /// Maximum 2-norm of independently re-evaluated scaled equilibrium rows.
    pub equilibrium_tolerance: f64,
    /// Maximum relative difference between requested and actually represented
    /// physical coordinate offsets, in (0,1). Large absolute coordinates must
    /// not silently change a numerical stencil by a significant fraction.
    pub maximum_stencil_distortion: f64,
    /// Central actions operate on scaled coordinates and scaled effort rows.
    pub difference: CentralDifferenceControls,
}

/// Original physical refusal or unresolved local numerical differentiation.
#[derive(Clone, Debug, PartialEq)]
pub enum LinearizationError {
    Invalid(&'static str),
    Equilibrium(EquilibriumError),
    Difference(CentralDifferenceError<EquilibriumError>),
    Generalized(GeneralizedError),
    Unbalanced { scaled_residual: f64, tolerance: f64 },
    Budget { used: usize, maximum: usize },
}
impl From<EquilibriumError> for LinearizationError {
    fn from(e: EquilibriumError) -> Self { Self::Equilibrium(e) }
}
impl From<PrescribedBaseMotionError> for LinearizationError {
    fn from(e: PrescribedBaseMotionError) -> Self { EquilibriumError::from(e).into() }
}
impl From<ArticulatedError> for LinearizationError {
    fn from(e: ArticulatedError) -> Self { EquilibriumError::from(e).into() }
}
impl From<FreeFloatingMotionError> for LinearizationError {
    fn from(e: FreeFloatingMotionError) -> Self { EquilibriumError::from(e).into() }
}
impl From<GeneralizedError> for LinearizationError {
    fn from(e: GeneralizedError) -> Self { Self::Generalized(e) }
}
impl core::fmt::Display for LinearizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated equilibrium response: {self:?}")
    }
}
impl std::error::Error for LinearizationError {}

/// One original-unit K*dq or C*dqd action, with two-stencil diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct JointTangentAction {
    pub effort: Vec<f64>,
    /// Absolute successive-stencil differences, in physical effort units.
    /// These are estimates, not derivative-error enclosures.
    pub disagreement: Vec<f64>,
    pub residual_evaluations: usize,
    pub refinements: usize,
    pub final_scaled_step: f64,
    /// Largest observed relative physical-offset representation discrepancy.
    pub maximum_stencil_distortion: f64,
}

/// Immutable original model/load binding and independently admitted root.
pub struct ArticulatedLinearization<'a> {
    model: &'a ArticulatedModel,
    state: PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'a>,
    bushings: EquilibriumBushings<'a>,
    position_scales: Vec<f64>,
    velocity_scales: Vec<f64>,
    effort_scales: Vec<f64>,
    difference: CentralDifferenceControls,
    maximum_stencil_distortion: f64,
    scaled_equilibrium_residual: f64,
    // Only the independently admitted affine-constraint owner supplies this.
    // Subtract its constant mean reaction before stencils to reduce cancellation.
    reaction_effort: Option<Vec<f64>>,
}

/// Bind a stationary or uniformly rotating root to its original loads/materials.
///
/// An equilibrium's public fields cannot replace independent balance admission.
/// The support has zero linear velocity/acceleration and all joint rates are
/// zero. World loads pass the original rotating-symmetry admission. Constant
/// held actuator efforts are supported; arbitrary callbacks, contact switches
/// and closed-loop constraints require their own complete tangents.
pub fn prepare<'a>(
    model: &'a ArticulatedModel,
    state: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'a>,
    bushings: EquilibriumBushings<'a>,
    controls: LinearizationControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedLinearization<'a>, LinearizationError> {
    prepare_with_reaction(model, state, loads, bushings, controls, None, &mut cancelled)
}

/// Internal common preparation. The constrained caller admits original rows and
/// reconstructs their original transpose before supplying the mean reaction.
/// The public unconstrained entry always passes None and retains its own gate.
pub(super) fn prepare_with_reaction<'a>(
    model: &'a ArticulatedModel,
    state: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'a>,
    bushings: EquilibriumBushings<'a>,
    controls: LinearizationControls<'_>,
    reaction: Option<&[f64]>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedLinearization<'a>, LinearizationError> {
    poll(&mut cancelled)?;
    let n = model.dof_count();
    if (n == 0 && reaction.is_none()) || model.link_count() > MAX_FREE_FLOATING_STEP_LINKS
        || !state.time_s.is_finite()
        || state.base.twist_body.linear != Twist::zero().linear
        || state.base.acceleration_body != Twist::zero()
        || state.joint_velocity.iter().any(|v| *v != 0.0)
        || !controls.equilibrium_tolerance.is_finite() || controls.equilibrium_tolerance <= 0.0
        || !controls.maximum_stencil_distortion.is_finite()
        || controls.maximum_stencil_distortion <= 0.0 || controls.maximum_stencil_distortion >= 1.0
        || [controls.coordinate_scales, controls.velocity_scales, controls.effort_scales]
            .iter().any(|v| v.len()!=n || v.iter().any(|x| !x.is_finite() || *x<=0.0))
    {
        return Err(LinearizationError::Invalid("bounded equilibrium and complete positive physical scales required"));
    }
    if reaction.is_some_and(|r| r.len()!=n || r.iter().any(|v| !v.is_finite())) {
        return Err(LinearizationError::Invalid("complete finite original reaction required"));
    }
    rotating::admit_frame(state, loads, &mut cancelled)?;
    validate_forward_dynamics_header(model, state.base.twist_body,
        &state.joint_position, &state.joint_velocity, loads.held.generalized_force,
        loads.held.gravity_world, loads.held.external_body_wrench)?;
    match bushings {
        EquilibriumBushings::Diagonal(b) => bushing::admit(model, loads, b, &mut cancelled)?,
        EquilibriumBushings::Coupled(b) => bushing::admit(model, loads, b, &mut cancelled)?,
    }
    let mut result = ArticulatedLinearization {
        model, state: state.clone(), loads, bushings,
        position_scales: controls.coordinate_scales.to_vec(),
        velocity_scales: controls.velocity_scales.to_vec(),
        effort_scales: controls.effort_scales.to_vec(), difference: controls.difference,
        maximum_stencil_distortion: controls.maximum_stencil_distortion,
        scaled_equilibrium_residual: 0.0,
        reaction_effort: reaction.map(|r| r.to_vec()),
    };
    // A genuinely empty physical model has no numerical derivative to probe.
    if n != 0 {
        central::differentiate(&vec![0.0;n], n, controls.difference,
            |_| Ok::<_, EquilibriumError>(vec![0.0;n]), || {
                poll(&mut cancelled).map_err(EquilibriumError::from)
            }).map_err(LinearizationError::Difference)?;
    }
    let value = result.evaluate(state.clone(), &mut cancelled)?;
    let scaled = scaled_rows(&value.physical_residual, &result.effort_scales)?;
    let residual = fs_solver::norm2(&scaled);
    if !residual.is_finite() || residual > controls.equilibrium_tolerance {
        return Err(LinearizationError::Unbalanced {
            scaled_residual: residual, tolerance: controls.equilibrium_tolerance,
        });
    }
    result.scaled_equilibrium_residual = residual;
    poll(&mut cancelled)?;
    Ok(result)
}

impl ArticulatedLinearization<'_> {
    pub fn state(&self) -> &PrescribedBaseMotionState { &self.state }
    pub fn scaled_equilibrium_residual(&self) -> f64 { self.scaled_equilibrium_residual }

    fn evaluate(&self, state: PrescribedBaseMotionState, cancelled: &mut impl FnMut() -> bool)
        -> Result<Evaluation, EquilibriumError> {
        let mut value = match self.bushings {
            EquilibriumBushings::Diagonal(b) =>
                evaluate_native(self.model, state, self.loads, b, None, cancelled),
            EquilibriumBushings::Coupled(b) =>
                evaluate_native(self.model, state, self.loads, b, None, cancelled),
        }?;
        if let Some(reaction) = &self.reaction_effort {
            for (i, &reaction) in reaction.iter().enumerate() {
                if i % 64 == 0 { poll(cancelled)?; }
                let residual = (Dd::from_f64(value.required[i])
                    - Dd::from_f64(value.applied[i]) - Dd::from_f64(reaction)).to_f64();
                if !residual.is_finite() {
                    return Err(EquilibriumError::Invalid("unrepresentable reaction-adjusted balance"));
                }
                value.physical_residual[i] = residual;
            }
        }
        Ok(value)
    }

    /// Original recursive inertia action. No mass matrix or unit velocity is
    /// constructed; the direction occupies the native acceleration argument.
    pub fn inertia_action(&self, acceleration: &[f64], mut cancelled: impl FnMut() -> bool)
        -> Result<Vec<f64>, LinearizationError> {
        poll(&mut cancelled)?;
        let n = self.model.dof_count();
        if acceleration.len()!=n || acceleration.iter().any(|v| !v.is_finite()) {
            return Err(LinearizationError::Invalid("complete finite acceleration direction required"));
        }
        let zero = vec![0.0;n];
        let value = inverse_dynamics(self.model,
            BaseState::stationary(self.state.base.world_from_base), &self.state.joint_position,
            &zero, acceleration, fs_ga::Vec3::default(),
            &vec![Wrench::default();self.model.link_count()])?;
        poll(&mut cancelled)?;
        Ok(value.generalized_force)
    }

    /// K*dq, retaining actual preload geometry and rotating-frame softening.
    pub fn stiffness_action(&self, displacement: &[f64], cancelled: impl FnMut() -> bool)
        -> Result<JointTangentAction, LinearizationError> {
        self.action(displacement, false, self.difference, cancelled)
    }

    /// C*dqd, retaining material damping AND the original velocity bias.
    /// C is not assumed symmetric, dissipative, or proportional to M or K.
    pub fn velocity_action(&self, velocity: &[f64], cancelled: impl FnMut() -> bool)
        -> Result<JointTangentAction, LinearizationError> {
        self.action(velocity, true, self.difference, cancelled)
    }

    pub(super) fn action(&self, direction: &[f64], rate: bool, controls: CentralDifferenceControls,
        cancelled: impl FnMut() -> bool) -> Result<JointTangentAction, LinearizationError> {
        let n = self.model.dof_count();
        if direction.len()!=n || direction.iter().any(|v| !v.is_finite()) {
            return Err(LinearizationError::Invalid("complete finite tangent direction required"));
        }
        let scales = if rate { &self.velocity_scales } else { &self.position_scales };
        let scaled: Vec<_> = direction.iter().zip(scales).map(|(v,s)| v/s).collect();
        if scaled.iter().zip(direction).any(|(scaled, physical)|
            !scaled.is_finite() || (*physical != 0.0 && *scaled == 0.0)) {
            return Err(LinearizationError::Invalid("scaled tangent direction is not representable"));
        }
        let cancel = RefCell::new(cancelled);
        let distortion = Cell::new(0.0_f64);
        let derivative = central::differentiate(&scaled, n, controls, |offset| {
            let mut state = self.state.clone();
            let target = if rate { &mut state.joint_velocity } else { &mut state.joint_position };
            for ((value, &increment), &scale) in target.iter_mut().zip(offset).zip(scales) {
                let next = increment.mul_add(scale, *value);
                if !next.is_finite() || (increment!=0.0 && next==*value) {
                    return Err(EquilibriumError::Invalid("physical tangent stencil is not representable"));
                }
                if increment != 0.0 {
                    let requested = Dd::from_f64(increment) * Dd::from_f64(scale);
                    let represented = Dd::from_f64(next) - Dd::from_f64(*value);
                    let error = ((represented-requested)/requested).abs().to_f64();
                    if !error.is_finite() || error > self.maximum_stencil_distortion {
                        return Err(EquilibriumError::Invalid("physical tangent stencil exceeds representation allowance"));
                    }
                    distortion.set(distortion.get().max(error));
                }
                *value = next;
            }
            let evaluated = self.evaluate(state, &mut *cancel.borrow_mut())?;
            scaled_rows(&evaluated.physical_residual, &self.effort_scales)
        }, || poll(&mut *cancel.borrow_mut()).map_err(EquilibriumError::from))
            .map_err(LinearizationError::Difference)?;
        physical_action(derivative, &self.effort_scales, distortion.get())
    }

    /// The original composite-rigid-body joint mass block, shared by both
    /// full and constrained response exports. Counts as one native mass query.
    pub(super) fn mass_matrix(&self, mut cancelled: impl FnMut() -> bool)
        -> Result<Vec<f64>, LinearizationError> {
        poll(&mut cancelled)?;
        let n = self.model.dof_count();
        let native = free_floating_terms(self.model,
            FreeFloatingBaseState {
                world_from_base:self.state.base.world_from_base, twist_body:Twist::zero(),
            }, &self.state.joint_position, &vec![0.0;n], &mut cancelled)?;
        let mut mass = vec![0.0;n*n];
        for i in 0..n { for j in 0..n { mass[i*n+j] = native.mass[6+i][6+j]; } }
        poll(&mut cancelled)?;
        Ok(mass)
    }

    /// Export a bounded full pencil to the EXISTING fs-modal spectral owners.
    /// Matrix-free actions above remain available for larger trees. The original
    /// composite-rigid-body owner supplies its symmetric joint/joint mass block;
    /// no independently differenced mass columns are projected into symmetry.
    /// K and C retain every numerical nonsymmetric entry. Work caps include
    /// every original residual query plus the one native mass extraction.
    pub fn quadratic_model(&self, maximum_dofs: usize, maximum_evaluations: usize,
        mut cancelled: impl FnMut() -> bool) -> Result<EquilibriumQuadraticModel, LinearizationError> {
        poll(&mut cancelled)?;
        let n = self.model.dof_count();
        if !(1..=64).contains(&maximum_dofs) || n>maximum_dofs
            || maximum_evaluations==0 || maximum_evaluations>16384 {
            return Err(LinearizationError::Invalid("explicit dense dimension and source-work bounds required"));
        }
        let mass = self.mass_matrix(&mut cancelled)?;
        let mut used = 1;
        let mut stiffness = vec![0.0;n*n];
        let mut velocity = vec![0.0;n*n];
        let mut stiffness_disagreement = vec![0.0;n*n];
        let mut velocity_disagreement = vec![0.0;n*n];
        let mut direction = vec![0.0;n];
        let mut maximum_stencil_distortion = 0.0_f64;
        for column in 0..n {
            for rate in [false,true] {
                poll(&mut cancelled)?;
                let remaining = maximum_evaluations.saturating_sub(used);
                if remaining<4 {
                    return Err(LinearizationError::Budget { used, maximum:maximum_evaluations });
                }
                let scales = if rate { &self.velocity_scales } else { &self.position_scales };
                direction.fill(0.0); direction[column] = scales[column];
                let mut options = self.difference;
                options.maximum_evaluations = options.maximum_evaluations.min(remaining);
                let action = self.action(&direction, rate, options, &mut cancelled)?;
                used += action.residual_evaluations;
                maximum_stencil_distortion = maximum_stencil_distortion.max(action.maximum_stencil_distortion);
                let (matrix, errors) = if rate {
                    (&mut velocity,&mut velocity_disagreement)
                } else { (&mut stiffness,&mut stiffness_disagreement) };
                for row in 0..n {
                    matrix[row*n+column] = action.effort[row]/scales[column];
                    errors[row*n+column] = action.disagreement[row]/scales[column];
                    if !matrix[row*n+column].is_finite() || !errors[row*n+column].is_finite()
                        || (action.effort[row] != 0.0 && matrix[row*n+column] == 0.0)
                        || (action.disagreement[row] != 0.0 && errors[row*n+column] == 0.0) {
                        return Err(LinearizationError::Invalid("physical tangent matrix is not representable"));
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(EquilibriumQuadraticModel {
            n,mass,velocity,stiffness,stiffness_disagreement,velocity_disagreement,
            source_evaluations:used,maximum_stencil_distortion,
        })
    }
}

fn scaled_rows(physical: &[f64], scales: &[f64]) -> Result<Vec<f64>, EquilibriumError> {
    let scaled: Vec<_> = physical.iter().zip(scales).map(|(v,s)| v/s).collect();
    if scaled.iter().zip(physical).any(|(scaled, physical)|
        !scaled.is_finite() || (*physical != 0.0 && *scaled == 0.0)) {
        return Err(EquilibriumError::Invalid("scaled physical effort is not representable"));
    }
    Ok(scaled)
}

fn physical_action(value: CentralDifference, scales: &[f64], distortion: f64)
    -> Result<JointTangentAction,LinearizationError> {
    let effort: Vec<_> = value.action.iter().zip(scales).map(|(v,s)|
        (Dd::from_f64(*v)*Dd::from_f64(*s)).to_f64()).collect();
    let disagreement: Vec<_> = value.disagreement.iter().zip(scales).map(|(v,s)| v*s).collect();
    if effort.iter().chain(&disagreement).any(|v| !v.is_finite())
        || effort.iter().zip(&value.action).any(|(physical, scaled)| *scaled != 0.0 && *physical == 0.0)
        || disagreement.iter().zip(&value.disagreement).any(|(physical, scaled)| *scaled != 0.0 && *physical == 0.0) {
        return Err(LinearizationError::Invalid("physical tangent action is not representable"));
    }
    Ok(JointTangentAction {
        effort,disagreement,residual_evaluations:value.evaluations,
        refinements:value.refinements,final_scaled_step:value.final_step,
        maximum_stencil_distortion:distortion,
    })
}

/// Complete local mechanical pencil in original physical joint coordinates.
/// Matrices are immutable; estimated stencil disagreement remains visible.
pub struct EquilibriumQuadraticModel {
    n:usize, mass:Vec<f64>, velocity:Vec<f64>, stiffness:Vec<f64>,
    stiffness_disagreement:Vec<f64>, velocity_disagreement:Vec<f64>,
    source_evaluations:usize, maximum_stencil_distortion:f64,
}
impl EquilibriumQuadraticModel {
    pub fn dimension(&self)->usize { self.n }
    pub fn mass(&self)->&[f64] { &self.mass }
    pub fn velocity_operator(&self)->&[f64] { &self.velocity }
    pub fn stiffness(&self)->&[f64] { &self.stiffness }
    pub fn stiffness_disagreement(&self)->&[f64] { &self.stiffness_disagreement }
    pub fn velocity_disagreement(&self)->&[f64] { &self.velocity_disagreement }
    pub fn source_evaluations(&self)->usize { self.source_evaluations }
    pub fn maximum_stencil_distortion(&self)->f64 { self.maximum_stencil_distortion }

    /// Native complete quadratic poles, including unstable growth rates.
    /// Neither the numerical stencil nor QR certifies spectral/stability error.
    pub fn poles(&self, maximum_sweeps_per_window:usize, cancelled:impl FnMut()->bool)
        ->Result<Vec<C64>,fs_modal::ModalError> {
        fs_modal::quadratic_eigenvalues_with_cancel(&self.mass,&self.velocity,&self.stiffness,
            self.n,maximum_sweeps_per_window,cancelled)
    }

    /// Native local mode shapes, residues and sensitivity preparation on the
    /// same M/C/K. No new eigensolver or proportional-damping fit exists.
    pub fn prepare_modes(&self,options:fs_modal::quadratic_modes::PoleOptions,
        cancelled:impl FnMut()->bool)
        ->Result<fs_modal::quadratic_modes::PreparedQuadraticModes<'_>,fs_modal::quadratic_modes::QuadraticModeError> {
        fs_modal::quadratic_modes::PreparedQuadraticModes::new(&self.mass,&self.velocity,
            &self.stiffness,self.n,options,cancelled)
    }
}
