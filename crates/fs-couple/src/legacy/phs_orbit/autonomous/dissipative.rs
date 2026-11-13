//! Self-sustained original physical feedback with internal nonlinear loss.

use super::AutonomousPortFeedback;
use crate::phs_orbit::{
    PhsFeedbackStateJacobian, PhsStorageHessianAction,
    dissipative::{DissipativeError, DissipativePower, NonlinearDissipation},
    port_flows_into,
    tangent::StorageTangent,
};
use core::cell::RefCell;
use fs_math::dd::Dd;
use fs_orbit::{
    AdaptiveHbOrbit, HbBudget, HbRefinementBudget, HbScales, OrbitProblem, ScaledHbOrbit,
    ShootBudget, VariationalShootControls, VariationalShootError, VariationalShootOrbit,
    solve_hb_autonomous_adaptive_scaled_seeded, solve_hb_autonomous_scaled_seeded,
    solve_shooting_variational,
};
use fs_phs::{PortHamiltonian, harmonic::equilibrium::LinearizationError};

type Result<T> = core::result::Result<T, DissipativeError>;

mod sensitivity;
pub use sensitivity::{
    DissipativeAutonomousParameterDirection, DissipativeAutonomousPowerObjective,
    DissipativeAutonomousSensitivityError, PreparedDissipativeAutonomousAdjoint,
    PreparedDissipativeAutonomousPowerAdjoint, PreparedDissipativeAutonomousSensitivity,
};

struct Scratch {
    effort: Vec<f64>,
    resistance: Vec<f64>,
    state_direction: Vec<f64>,
    effort_direction: Vec<f64>,
    resistance_direction: Vec<f64>,
    field: Vec<f64>,
    input: Vec<f64>,
    flows: Vec<f64>,
    input_direction: Vec<f64>,
    feedback_derivative: Vec<f64>,
}

/// Complete self-excited source `xdot=(J-R)e-D(x,e)+G*u(x,G^T e)`.
///
/// Original storage, all material-memory states, physical ports and the SAME
/// nonlinear resistance used by fs-phs time stepping remain present. The
/// feedback supplies physical port work; D absorbs independent internal work,
/// even outside the range of G. Neither is replaced by an effective damping.
/// The feedback Jacobian is TOTAL, as in [`PhsFeedbackStateJacobian`]; D receives
/// separate state and effort directions, whose chain is added exactly once.
///
/// Preparation requires complete storage and source tangents and admits at most
/// 256 states/ports. It retains `6*n+3*m+n*m` scalar scratch entries plus the
/// existing storage-tangent owner. Warm field/tangent/power evaluation allocates
/// nothing in this adapter; source callbacks retain their own resource behavior.
/// Native solver controls bound additional numerical work and memory.
///
/// Its source-bound sensitivity owner includes original storage, J/R/G,
/// feedback and independent D parameter laws in waveform, natural-frequency
/// and physical power/work gradients. Smooth, read-only sources are required;
/// hidden history, corners and hybrid events are unsupported.
pub struct DissipativeAutonomousHamiltonian<'a> {
    system: &'a PortHamiltonian,
    feedback: &'a dyn AutonomousPortFeedback,
    feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
    dissipation: NonlinearDissipation<'a>,
    tangent: StorageTangent<'a>,
    scratch: RefCell<Scratch>,
    failure: RefCell<Option<DissipativeError>>,
}

impl<'a> DissipativeAutonomousHamiltonian<'a> {
    /// Use the original storage's declared COMPLETE continuous tangent.
    pub fn new(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
        dissipation: NonlinearDissipation<'a>,
    ) -> Result<Self> {
        Self::prepare(system, feedback, feedback_jacobian, dissipation, None)
    }

    /// Supply the COMPLETE original Hessian, including every memory cross term.
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
        dissipation: NonlinearDissipation<'a>,
        hessian: &'a PhsStorageHessianAction<'a>,
    ) -> Result<Self> {
        Self::prepare(
            system,
            feedback,
            feedback_jacobian,
            dissipation,
            Some(hessian),
        )
    }

    fn prepare(
        system: &'a PortHamiltonian,
        feedback: &'a dyn AutonomousPortFeedback,
        feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
        dissipation: NonlinearDissipation<'a>,
        hessian: Option<&'a PhsStorageHessianAction<'a>>,
    ) -> Result<Self> {
        let (n, m) = (system.state_dim(), system.port_dim());
        if !(1..=256).contains(&n)
            || m > 256
            || feedback.port_count() != m
            || dissipation.max_tangent_actions < n
        {
            return Err(DissipativeError::Invalid(
                "bounded original states, complete feedback ports and all loss tangent actions required",
            ));
        }
        let tangent = match hessian {
            Some(action) => StorageTangent::hessian(n, action),
            None => StorageTangent::from_storage(system),
        };
        if !tangent.available() {
            return Err(DissipativeError::Invalid(
                "complete original storage tangent required",
            ));
        }
        Ok(Self {
            system,
            feedback,
            feedback_jacobian,
            dissipation,
            tangent,
            scratch: RefCell::new(Scratch {
                effort: vec![0.0; n],
                resistance: vec![0.0; n],
                state_direction: vec![0.0; n],
                effort_direction: vec![0.0; n],
                resistance_direction: vec![0.0; n],
                field: vec![0.0; n],
                input: vec![0.0; m],
                flows: vec![0.0; m],
                input_direction: vec![0.0; m],
                feedback_derivative: vec![0.0; n * m],
            }),
            failure: RefCell::new(None),
        })
    }

    fn evaluate(&self, state: &[f64], scratch: &mut Scratch) -> Result<f64> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(DissipativeError::Invalid(
                "complete finite original state required",
            ));
        }
        scratch.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut scratch.effort);
        if scratch.effort.iter().any(|v| !v.is_finite())
            || !port_flows_into(self.system, &scratch.effort, &mut scratch.flows)
        {
            return Err(DissipativeError::SourceRefused);
        }
        scratch.input.fill(f64::NAN);
        if !self
            .feedback
            .input(state, &scratch.flows, &mut scratch.input)
            || scratch.input.iter().any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused);
        }
        self.dissipation
            .evaluate_into(state, &scratch.effort, &mut scratch.resistance, || false)
            .map_err(DissipativeError::Dissipation)
    }

    /// Complete original autonomous field, with one final rounding per row.
    pub fn vector_field_into(&self, state: &[f64], out: &mut [f64]) -> Result<()> {
        let mut scratch = self.scratch.borrow_mut();
        self.evaluate(state, &mut scratch)?;
        self.system
            .vector_field_with_dissipation_from_effort_into(
                &scratch.effort,
                &scratch.input,
                &scratch.resistance,
                out,
            )
            .map_err(|_| DissipativeError::SourceRefused)
    }

    /// `(J-R)Hxx-D_x-D_e Hxx+G du_total/dx` in original state coordinates.
    /// Passive pointwise loss need not have a positive differential slope; the
    /// actual source tangent is retained without secants or symmetrization.
    pub fn state_jacobian_into(&self, state: &[f64], out: &mut [f64]) -> Result<()> {
        let (n, m) = (self.system.state_dim(), self.system.port_dim());
        if out.len() != n * n {
            return Err(DissipativeError::Invalid(
                "complete row-major state Jacobian required",
            ));
        }
        let mut scratch = self.scratch.borrow_mut();
        self.evaluate(state, &mut scratch)?;
        let Scratch {
            effort,
            flows,
            feedback_derivative,
            state_direction,
            effort_direction,
            resistance_direction,
            input_direction,
            field,
            ..
        } = &mut *scratch;
        feedback_derivative.fill(f64::NAN);
        if !(self.feedback_jacobian)(state, flows, feedback_derivative)
            || feedback_derivative.iter().any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused);
        }
        state_direction.fill(0.0);
        for column in 0..n {
            state_direction[column] = 1.0;
            if !self
                .tangent
                .apply_hessian(state, state_direction, effort_direction)
            {
                return Err(DissipativeError::SourceRefused);
            }
            resistance_direction.fill(f64::NAN);
            if !(self.dissipation.tangent)(
                state,
                effort,
                state_direction,
                effort_direction,
                resistance_direction,
            ) || resistance_direction.iter().any(|v| !v.is_finite())
            {
                return Err(DissipativeError::Dissipation(
                    LinearizationError::DissipationRefused {
                        column: Some(column),
                    },
                ));
            }
            for port in 0..m {
                input_direction[port] = feedback_derivative[port * n + column];
            }
            self.system
                .vector_field_with_dissipation_from_effort_into(
                    effort_direction,
                    input_direction,
                    resistance_direction,
                    field,
                )
                .map_err(|_| DissipativeError::SourceRefused)?;
            for row in 0..n {
                out[row * n + column] = field[row];
            }
            state_direction[column] = 0.0;
        }
        Ok(())
    }

    /// Original energy and separate feedback, R and D powers. Energy rate is
    /// independently contracted from the complete field. This instantaneous
    /// identity does not certify a cycle, passivity or temporal convergence.
    pub fn power(&self, state: &[f64]) -> Result<DissipativePower> {
        let mut scratch = self.scratch.borrow_mut();
        let nonlinear_power = self.evaluate(state, &mut scratch)?;
        let Scratch {
            effort,
            input,
            resistance,
            field,
            ..
        } = &mut *scratch;
        self.system
            .vector_field_with_dissipation_from_effort_into(effort, input, resistance, field)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let energy = self.system.hamiltonian(state);
        let port_power = self
            .system
            .supplied_power_from_effort(effort, input)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let resistance_power = self
            .system
            .resistance_power_from_effort(effort)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let energy_rate = effort
            .iter()
            .zip(field)
            .fold(Dd::ZERO, |s, (&e, f)| {
                s + Dd::from_f64(e) * Dd::from_f64(*f)
            })
            .to_f64();
        if [
            energy,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused);
        }
        Ok(DissipativePower {
            energy,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        })
    }

    fn finish<T>(&self, result: Result<T>) -> Result<T> {
        match self.failure.borrow_mut().take() {
            Some(failure) => Err(failure),
            None => result,
        }
    }

    /// Resolve amplitude, frequency and every original harmonic together. The
    /// seed packs original DC/complex coefficients followed by omega in rad/s;
    /// its full waveform defines the native Parseval phase condition. No first
    /// coordinate is assumed to oscillate and no amplitude is artificially fixed.
    pub fn solve_harmonic_balance_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<ScaledHbOrbit> {
        *self.failure.borrow_mut() = None;
        self.finish(
            solve_hb_autonomous_scaled_seeded(&Adapter(self), seed, budget, scales)
                .map_err(DissipativeError::Orbit),
        )
    }

    /// Native retained/omitted-equation and waveform/frequency refinement of
    /// this SAME source. Finite-band acceptance does not prove an infinite tail,
    /// primitive period, stability, attraction or physical excitation validity.
    pub fn solve_harmonic_balance_adaptive_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
        scales: &HbScales,
    ) -> Result<AdaptiveHbOrbit> {
        *self.failure.borrow_mut() = None;
        self.finish(
            solve_hb_autonomous_adaptive_scaled_seeded(
                &Adapter(self),
                seed,
                budget,
                refinement,
                scales,
            )
            .map_err(DissipativeError::Orbit),
        )
    }

    /// Independent original-ODE shooting with solved period and the exact
    /// derivative of the SAME discrete RK4 map. The full-state phase plane and
    /// excursion test remain native; no neutral multiplier or memory is removed.
    pub fn solve_shooting_variational(
        &self,
        seed: &[f64],
        period_guess: f64,
        budget: &ShootBudget,
        controls: &VariationalShootControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VariationalShootOrbit> {
        if cancelled() {
            return Err(DissipativeError::Shooting(VariationalShootError::Cancelled));
        }
        *self.failure.borrow_mut() = None;
        let result = solve_shooting_variational(
            &Adapter(self),
            seed,
            period_guess,
            budget,
            controls,
            &mut cancelled,
        );
        if matches!(result, Err(VariationalShootError::Cancelled)) {
            return Err(DissipativeError::Shooting(VariationalShootError::Cancelled));
        }
        self.finish(result.map_err(DissipativeError::Shooting))
    }
}

struct Adapter<'s, 'a>(&'s DissipativeAutonomousHamiltonian<'a>);
impl Adapter<'_, '_> {
    fn record(&self, result: Result<()>) -> bool {
        if let Err(error) = result {
            let mut failure = self.0.failure.borrow_mut();
            if failure.is_none() {
                *failure = Some(error);
            }
            false
        } else {
            true
        }
    }
    fn jacobian(&self, _time: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.record(self.0.state_jacobian_into(state, out))
    }
}
impl OrbitProblem for Adapter<'_, '_> {
    fn dim(&self) -> usize {
        self.0.system.state_dim()
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _time: f64, state: &[f64], out: &mut [f64]) {
        if !self.record(self.0.vector_field_into(state, out)) {
            out.fill(f64::NAN);
        }
    }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        Some(Self::jacobian)
    }
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<fs_orbit::TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}

#[cfg(test)]
mod tests;
