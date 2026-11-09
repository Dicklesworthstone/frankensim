//! Equilibrium temperature response through the original native mechanical solve.
//!
//! Pair the complete mechanical mixed block J_T with the PHYSICAL bulk entropy
//! block S_u at one admitted displacement. One original K solve gives
//! du = -K^-1 J_T dT and dS = (C_epsilon/T) dT + S_u du. The transpose uses
//! S_u^T and J_T^T separately; stabilization is never reclassified as entropy.

use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::{
    PreparedThermalResidualJacobian3, PreparedThermoelasticEntropy3,
    ThermalJacobianLimits3, ThermoelasticCellPotential3, ThermoelasticEntropyLimits3,
    ThermoelasticMotion3,
};
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::{CutElasticity3, ElasticityError3, ElasticitySolution3};

/// Explicit preparation bounds; the two native blocks retain their own point
/// and stencil allowances. No geometry or displacement-sized cache is copied.
#[derive(Debug, Clone, Copy)]
pub struct ThermoelasticLinearResponseLimits3 {
    pub mechanical: ThermalJacobianLimits3,
    pub entropy: ThermoelasticEntropyLimits3,
    /// Bound BOTH original physical and independent displacement dimensions.
    pub maximum_native_dofs: usize,
    /// Aggregate preparation peak, including both cached blocks and their
    /// scratch. Scalar/scale slots are charged as eight bytes; fixed stack
    /// scratch, Vec/Arc metadata and allocator overhead are excluded.
    pub maximum_workspace_components: usize,
}

/// Bounds for one original controlled K solve and the surrounding actions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticLinearSolveSettings3 {
    pub relative_residual_tolerance: f64,
    /// Total native CG iterations, including any true-residual corrections.
    /// Zero permits only a right-hand side already solved by the zero iterate.
    pub maximum_iterations: usize,
    pub poll_iterations: usize,
    /// Conservative scalar storage allowance for native CG, its bounded
    /// history/report, adaptive physical expansion, the RHS and returned
    /// vectors. Excludes the prepared cache, borrowed inputs, constant stack
    /// scratch, Vec metadata and allocator overhead. No matrix is assembled.
    pub maximum_workspace_components: usize,
}

/// Coupled displacement/physical-entropy response to a cell temperature vector.
///
/// The displacement is an ALGEBRAIC direction, not a new admitted physical
/// state. Its original solver reports a recomputed true-Euclidean residual.
/// For a direction in kelvin, entropy_direction_j_k is in joules per kelvin.
/// Cell rows retain the original active-cell order.
#[derive(Debug, Clone)]
pub struct ThermoelasticTemperatureResponse3 {
    pub displacement: ElasticitySolution3,
    pub entropy_direction_j_k: Vec<f64>,
}

/// Exact implicit transpose of the paired response, to native solve tolerance.
///
/// temperature_derivative pairs with the original cell temperature direction;
/// its units are those of the authored scalar objective per kelvin. The
/// algebraic adjoint retains original true-residual evidence and iteration use.
#[derive(Debug, Clone)]
pub struct ThermoelasticTemperaturePullback3 {
    pub temperature_derivative: Vec<f64>,
    pub adjoint: ElasticitySolution3,
}

/// Coherent cached J_T, S_u and signed elastic fixed-total-strain capacity.
///
/// Both blocks are prepared internally at the SAME physical displacement and
/// original material state. Geometry, source selection, scales, reference
/// density and independent external forces are frozen. Optional embedded
/// motion includes its declared local temperature slope; strong box supports
/// remain homogeneous. Material/geometry owner, space and exact scale checks
/// are inherited from each original block before any solve can run.
///
/// Preparation admits physical strain but does not certify that the supplied
/// displacement solves an equilibrium. For equilibrium derivatives, supply
/// the accepted native equilibrium coefficients and the same support law.
/// Reprepare after changing displacement or support law, even without a material
/// refresh. No reference to caller displacement/motion is retained.
///
/// This is an elastic contribution only. The original caloric owner must add
/// its own entropy/capacity exactly once, and transport remains external.
/// Multiplying dS by each original cell temperature gives its elastic heat-
/// capacity action. The condensed stabilized block need not be symmetric or
/// positive, and this API makes no complete heat-step or passivity claim.
#[derive(Debug)]
pub struct PreparedThermoelasticLinearResponse3 {
    mechanical: PreparedThermalResidualJacobian3,
    entropy: PreparedThermoelasticEntropy3,
    workspace_components: usize,
    points: usize,
    action_terms: usize,
}

#[derive(Clone, Copy)]
enum Owner<'a> {
    Cartesian(&'a CutElasticity3),
    Adaptive(&'a AdaptiveElasticity3),
}

impl<'a> Owner<'a> {
    fn raw(self) -> &'a CutElasticity3 {
        match self {
            Self::Cartesian(operator) => operator,
            Self::Adaptive(operator) => &operator.raw,
        }
    }
}

macro_rules! native {
    ($owner:expr, $method:ident $(, $argument:expr)* $(,)?) => {
        match $owner {
            Owner::Cartesian(operator) => operator.$method($($argument),*),
            Owner::Adaptive(operator) => operator.$method($($argument),*),
        }
    };
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

fn poll(
    completed: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    if checkpoint(completed).is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(
    values: &[f64],
    completed: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    for chunk in values.chunks(192) {
        poll(completed, checkpoint)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite thermoelastic linear response field"));
        }
    }
    Ok(())
}

impl PreparedThermoelasticLinearResponse3 {
    #[must_use]
    pub fn n(&self) -> usize {
        self.mechanical.n()
    }

    #[must_use]
    pub fn cells(&self) -> usize {
        self.entropy.cells()
    }

    /// Physical source potential at the same admitted displacement.
    /// The signed capacity holds TOTAL strain fixed, before equilibration.
    #[must_use]
    pub fn cell_potential(&self, cell: usize) -> Option<&ThermoelasticCellPotential3> {
        self.entropy.cell_potential(cell)
    }

    /// Sum of the two original point passes, including mechanical supports.
    #[must_use]
    pub const fn points(&self) -> usize {
        self.points
    }

    /// Conservative peak scalar slots for the paired preparation.
    #[must_use]
    pub const fn workspace_components(&self) -> usize {
        self.workspace_components
    }

    /// Combined cached stencil and vector-combination work per action.
    /// Original K applications and linear finite checks are additional.
    #[must_use]
    pub const fn action_terms(&self) -> usize {
        self.action_terms
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        owner: Owner<'_>,
        u: &[f64],
        motion: Option<ThermoelasticMotion3<'_>>,
        maximum_kinematic_strain: f64,
        limits: ThermoelasticLinearResponseLimits3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(0, &mut checkpoint)?;
        let raw = owner.raw();
        let n = native!(owner, n);
        if n == 0 || n > limits.maximum_native_dofs
            || raw.n() > limits.maximum_native_dofs
        {
            return Err(invalid("thermoelastic linear response native dimension allowance exhausted"));
        }
        // The native mixed block charges 27 cell slots + 5 ghost slots.
        // The physical entropy block charges 30 cell slots. Charge both peaks
        // together before allocation, source evaluation or support callbacks.
        let workspace_components = raw.cells.len().checked_mul(57)
            .and_then(|cells| raw.ghosts.len().checked_mul(5)
                .and_then(|ghosts| cells.checked_add(ghosts)))
            .ok_or_else(|| invalid("thermoelastic linear response cache count overflow"))?;
        if workspace_components > limits.maximum_workspace_components {
            return Err(invalid("thermoelastic linear response cache allowance exhausted"));
        }
        let mechanical = match motion {
            Some(motion) => native!(owner, prepare_thermal_residual_jacobian_with_motion,
                u, motion.displacement, motion.temperature_derivative,
                maximum_kinematic_strain, limits.mechanical, || checkpoint(0))?,
            None => native!(owner, prepare_thermal_residual_jacobian,
                u, maximum_kinematic_strain, limits.mechanical, || checkpoint(0))?,
        };
        let entropy = native!(owner, prepare_thermoelastic_entropy,
            u, maximum_kinematic_strain, limits.entropy, || checkpoint(0))?;
        for cell in 0..entropy.cells() {
            poll(0, &mut checkpoint)?;
            let potential = entropy.cell_potential(cell).expect("prepared entropy cell");
            if potential.temperature_k <= 0.0
                || !(potential.heat_capacity_j_k / potential.temperature_k).is_finite()
            {
                return Err(invalid("nonfinite native elastic entropy-temperature partial"));
            }
        }
        let points = mechanical.points().checked_add(entropy.points())
            .ok_or_else(|| invalid("thermoelastic linear response point count overflow"))?;
        let action_terms = mechanical.action_terms().checked_add(entropy.action_terms())
            .and_then(|terms| terms.checked_add(entropy.cells()))
            .and_then(|terms| terms.checked_add(n))
            .ok_or_else(|| invalid("thermoelastic linear response action count overflow"))?;
        poll(0, &mut checkpoint)?;
        Ok(Self { mechanical, entropy, workspace_components, points, action_terms })
    }

    fn admit_solve(
        &self,
        owner: Owner<'_>,
        settings: ThermoelasticLinearSolveSettings3,
        checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(0, checkpoint)?;
        if native!(owner, n) != self.n() || native!(owner, cells) != self.cells()
            || !settings.relative_residual_tolerance.is_finite()
            || settings.relative_residual_tolerance <= 0.0
            || settings.relative_residual_tolerance >= 1.0
            || settings.poll_iterations == 0
        {
            return Err(invalid("invalid thermoelastic linear response solve controls or space"));
        }
        // Original Cartesian CG uses at most ten independent vector slots
        // including our RHS. Adaptive applies add two physical vectors.
        // Reserve one cell output and three bounded-history slots, including
        // the native report clone and Vec growth (minimum capacity four).
        let history = settings.maximum_iterations.min(settings.poll_iterations).max(4);
        let workspace = self.n().checked_mul(10)
            .and_then(|n| owner.raw().n().checked_mul(2).and_then(|raw| n.checked_add(raw)))
            .and_then(|n| n.checked_add(self.cells()))
            .and_then(|n| history.checked_mul(3).and_then(|history| n.checked_add(history)))
            .ok_or_else(|| invalid("thermoelastic linear response solve workspace overflow"))?;
        if workspace > settings.maximum_workspace_components {
            return Err(invalid("thermoelastic linear response solve workspace exhausted"));
        }
        Ok(())
    }

    fn forward(
        &self,
        owner: Owner<'_>,
        temperature_direction: &[f64],
        settings: ThermoelasticLinearSolveSettings3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticTemperatureResponse3, ElasticityError3> {
        self.admit_solve(owner, settings, &mut checkpoint)?;
        if temperature_direction.len() != self.cells() {
            return Err(invalid("thermoelastic temperature direction shape mismatch"));
        }
        finite(temperature_direction, 0, &mut checkpoint)?;
        let mut rhs = vec![0.0; self.n()];
        // Native binding/scale checks precede the original solve. The immutable
        // owner cannot change between the two cached block applications.
        native!(owner, thermal_residual_temperature_apply_into,
            &self.mechanical, temperature_direction, &mut rhs, || checkpoint(0))?;
        for chunk in rhs.chunks_mut(192) {
            poll(0, &mut checkpoint)?;
            for value in chunk {
                *value = -*value;
            }
        }
        let displacement = native!(owner, solve_controlled, &rhs,
            settings.relative_residual_tolerance, settings.maximum_iterations,
            settings.poll_iterations, &mut checkpoint)?;
        let completed = displacement.iterations();
        drop(rhs);
        let mut entropy_direction_j_k = vec![0.0; self.cells()];
        native!(owner, thermoelastic_entropy_displacement_apply_into,
            &self.entropy, displacement.coefficients(), &mut entropy_direction_j_k,
            || checkpoint(completed))?;
        for (cell, value) in entropy_direction_j_k.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll(completed, &mut checkpoint)?;
            }
            let potential = self.entropy.cell_potential(cell).expect("prepared entropy cell");
            *value += (potential.heat_capacity_j_k / potential.temperature_k)
                * temperature_direction[cell];
        }
        finite(&entropy_direction_j_k, completed, &mut checkpoint)?;
        poll(completed, &mut checkpoint)?;
        Ok(ThermoelasticTemperatureResponse3 { displacement, entropy_direction_j_k })
    }

    fn transpose(
        &self,
        owner: Owner<'_>,
        displacement_cotangent: Option<&[f64]>,
        entropy_cotangent: &[f64],
        settings: ThermoelasticLinearSolveSettings3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticTemperaturePullback3, ElasticityError3> {
        self.admit_solve(owner, settings, &mut checkpoint)?;
        if entropy_cotangent.len() != self.cells()
            || displacement_cotangent.is_some_and(|dual| dual.len() != self.n())
        {
            return Err(invalid("thermoelastic response cotangent shape mismatch"));
        }
        finite(entropy_cotangent, 0, &mut checkpoint)?;
        if let Some(dual) = displacement_cotangent {
            finite(dual, 0, &mut checkpoint)?;
        }
        let mut rhs = vec![0.0; self.n()];
        native!(owner, thermoelastic_entropy_displacement_transpose_into,
            &self.entropy, entropy_cotangent, &mut rhs, || checkpoint(0))?;
        if let Some(dual) = displacement_cotangent {
            for (i, (value, dual)) in rhs.iter_mut().zip(dual).enumerate() {
                if i % 192 == 0 {
                    poll(0, &mut checkpoint)?;
                }
                *value += dual;
            }
            finite(&rhs, 0, &mut checkpoint)?;
        }
        // The ORIGINAL displacement operator is symmetric, including its
        // Nitsche form and T^T K T map. Hence the same native controlled solve
        // is K^T z = u_bar + S_u^T S_bar. The condensed block need not be symmetric.
        let adjoint = native!(owner, solve_controlled, &rhs,
            settings.relative_residual_tolerance, settings.maximum_iterations,
            settings.poll_iterations, &mut checkpoint)?;
        let completed = adjoint.iterations();
        drop(rhs);
        let mut temperature_derivative = vec![0.0; self.cells()];
        native!(owner, thermal_residual_temperature_transpose_into,
            &self.mechanical, adjoint.coefficients(), &mut temperature_derivative,
            || checkpoint(completed))?;
        for (cell, value) in temperature_derivative.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll(completed, &mut checkpoint)?;
            }
            let potential = self.entropy.cell_potential(cell).expect("prepared entropy cell");
            *value = (potential.heat_capacity_j_k / potential.temperature_k)
                * entropy_cotangent[cell] - *value;
        }
        finite(&temperature_derivative, completed, &mut checkpoint)?;
        poll(completed, &mut checkpoint)?;
        Ok(ThermoelasticTemperaturePullback3 { temperature_derivative, adjoint })
    }
}

macro_rules! linear_response_api {
    ($operator:ty, $variant:ident) => {
        impl $operator {
            /// Prepare coherent complete-mechanical and physical-entropy blocks.
            ///
            /// The base displacement is physically admitted; algebraic response
            /// directions are not restricted by its strain-amplitude budget.
            /// Supply accepted equilibrium coefficients and the same optional
            /// embedded support motion for equilibrium derivatives. This does
            /// not solve or independently certify the base equilibrium.
            pub fn prepare_thermoelastic_linear_response(
                &self,
                u: &[f64],
                motion: Option<ThermoelasticMotion3<'_>>,
                maximum_kinematic_strain: f64,
                limits: ThermoelasticLinearResponseLimits3,
                mut checkpoint: impl FnMut() -> ControlFlow<()>,
            ) -> Result<PreparedThermoelasticLinearResponse3, ElasticityError3> {
                PreparedThermoelasticLinearResponse3::prepare(
                    Owner::$variant(self), u, motion, maximum_kinematic_strain,
                    limits, |_| checkpoint(),
                )
            }

            /// Solve K du = -J_T dT, then return physical
            /// dS = (C_epsilon/T) dT + S_u du in original cell order.
            ///
            /// Uses the original controlled native solver and true-residual
            /// gate. No source, geometry or support callbacks run during reuse.
            /// Every error/cancellation discards all tentative output; neither
            /// the operator nor the prepared state is changed.
            pub fn solve_thermoelastic_temperature_response(
                &self,
                prepared: &PreparedThermoelasticLinearResponse3,
                temperature_direction: &[f64],
                settings: ThermoelasticLinearSolveSettings3,
                checkpoint: impl FnMut(usize) -> ControlFlow<()>,
            ) -> Result<ThermoelasticTemperatureResponse3, ElasticityError3> {
                prepared.forward(Owner::$variant(self), temperature_direction, settings, checkpoint)
            }

            /// Transpose of the SAME coupled (du,dS) response in one K solve.
            ///
            /// None means a zero displacement cotangent. Entropy cotangents
            /// pair with cell-integrated physical entropy, without extra volume
            /// or temperature factors. Strong-support displacement cotangents
            /// do no work. For the elastic heat-capacity action T*dS, supply
            /// T times its output cotangent. Caloric/transport terms stay with
            /// their original owners.
            pub fn solve_thermoelastic_temperature_pullback(
                &self,
                prepared: &PreparedThermoelasticLinearResponse3,
                displacement_cotangent: Option<&[f64]>,
                entropy_cotangent: &[f64],
                settings: ThermoelasticLinearSolveSettings3,
                checkpoint: impl FnMut(usize) -> ControlFlow<()>,
            ) -> Result<ThermoelasticTemperaturePullback3, ElasticityError3> {
                prepared.transpose(Owner::$variant(self), displacement_cotangent,
                    entropy_cotangent, settings, checkpoint)
            }
        }
    };
}

linear_response_api!(CutElasticity3, Cartesian);
linear_response_api!(AdaptiveElasticity3, Adaptive);
