//! Full finite-deformation motion on the original sealed-gas volume chart.
//!
//! `fs-time` owns the matrix-free implicit midpoint trajectory (generalized
//! alpha with rho_inf=1). The membrane and original ideal/PR parcel own every
//! trial force, tangent, geometry, pressure and storage evaluation. Its signed
//! discrete energy defect is measured and gated, never reported as heat.

mod operator;

use super::*;
use fs_math::dd::Dd;
use fs_time::galpha::{
    ImplicitSolveConfig, ImplicitStepTelemetry, OperatorGeneralizedAlpha, SecondOrderState,
    TimeSolveError,
};
use operator::Operator;

/// Explicit source, nonlinear-work and publication bounds for native motion.
#[derive(Clone, Copy, Debug)]
pub struct GasMembraneMotionControls {
    /// Original isentropic source inversion; no caloric law is substituted.
    pub gas: GasEntropySolve,
    /// Original Newton--Krylov controls for each implicit time step.
    pub time: ImplicitSolveConfig,
    /// State-admission, force, tangent and final-source callbacks per step.
    /// Repeated callback coordinates reuse their native gas trial.
    pub maximum_evaluations: usize,
    /// Bound on retained accepted steps and their original native telemetry.
    pub maximum_steps: usize,
    /// Caller Krylov storage cap, additionally limited to 256 directions.
    pub maximum_krylov_restart: usize,
    /// Maximum absolute enclosing discrete energy defect [J] per interval.
    /// This is a numerical acceptance gate, not a time-accuracy certificate.
    pub energy_tolerance_j: f64,
}

impl GasMembraneMotionControls {
    fn validate(self) -> Result<(), GasMembraneMotionError> {
        self.time
            .newton
            .validate()
            .map_err(GasMembraneError::Newton)?;
        if self.time.max_newton_iterations == 0
            || self.maximum_evaluations == 0
            || self.maximum_steps == 0
            || self.maximum_krylov_restart == 0
            || self.maximum_krylov_restart > MAX_VOLUME_COORDINATES
            || self.time.newton.linear_restart > self.maximum_krylov_restart
            || !self.energy_tolerance_j.is_finite()
            || self.energy_tolerance_j <= 0.0
        {
            return Err(GasMembraneMotionError::Invalid(
                "positive motion budgets, admitted Krylov restart and finite energy allowance required",
            ));
        }
        Ok(())
    }
}

/// No refusal changes the accepted mechanics, gas inventory or source clock.
#[derive(Debug)]
pub enum GasMembraneMotionError {
    Invalid(&'static str),
    /// Original material, geometry, pressure or gas-state refusal.
    Physical(GasMembraneError),
    /// Original structural integrator/nonlinear solve refusal.
    Time(TimeSolveError),
    EvaluationBudget,
    StepBudget,
    /// Actual signed storage-minus-work defect exceeded the caller allowance.
    Energy {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for GasMembraneMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sealed gas/membrane motion: {self:?}")
    }
}
impl std::error::Error for GasMembraneMotionError {}
impl From<GasMembraneError> for GasMembraneMotionError {
    fn from(error: GasMembraneError) -> Self {
        match error {
            GasMembraneError::Cancelled => Self::Cancelled,
            other => Self::Physical(other),
        }
    }
}
impl From<MembraneError> for GasMembraneMotionError {
    fn from(error: MembraneError) -> Self {
        GasMembraneError::from(error).into()
    }
}
impl From<TimeSolveError> for GasMembraneMotionError {
    fn from(error: TimeSolveError) -> Self {
        match error {
            TimeSolveError::Cancelled => Self::Cancelled,
            other => Self::Time(other),
        }
    }
}

/// Measured endpoint changes and work on one accepted physical interval.
#[derive(Clone, Debug)]
pub struct GasMembraneMotionStep {
    /// Original accepted integrator clock and nonlinear convergence report.
    pub time: ImplicitStepTelemetry,
    pub kinetic_energy_change_j: f64,
    pub membrane_energy_change_j: f64,
    /// Difference of the native reference-relative gas energy changes.
    pub gas_energy_change_j: f64,
    /// Reversible work by the gas, equal to minus its internal-energy change.
    pub work_by_gas_j: f64,
    /// Native fixed-perimeter swept volume between the accepted shapes [m3].
    pub swept_volume_change_m3: f64,
    /// Held absolute exterior pressure times the actual native swept volume.
    pub exterior_pressure_work_j: f64,
    /// Held nodal dead forces dotted with the accepted displacement increment.
    pub dead_load_work_j: f64,
    /// Delta kinetic + Delta membrane + Delta gas + exterior work - dead work.
    /// This numerical defect is not material damping or a heat source.
    pub integration_defect_j: f64,
    pub evaluations: usize,
    pub gas_trials: usize,
}

/// Continuing full-coordinate motion on one immutable cavity/source binding.
///
/// The same prescribed perimeter, original reference masses, source entropy,
/// composition and constant absolute exterior pressure apply throughout. Every
/// accepted shape carries the actual native gas volume, pressure and temperature.
/// Dead forces may change between held intervals. No source is advanced until
/// the native solve, fresh endpoint source check and energy gate all succeed.
///
/// This is a spatially uniform reversible gas coupled to nonlinear membrane
/// mechanics. It adds no gas inertia/waves, heat or species exchange, material
/// damping, enclosure contact, remeshing or moving supports. Time refinement is
/// required; the energy gate alone does not resolve an oscillation's frequency.
#[derive(Clone)]
pub struct GasMembraneMotion<'m, 'p> {
    model: &'m GasMembrane<'p>,
    state: SecondOrderState,
    accepted: GasMembraneEvaluation,
    masses_kg: Vec<f64>,
}

impl<'m, 'p> GasMembraneMotion<'m, 'p> {
    /// Start at time zero from an admitted full shape and velocity. Prescribed
    /// coordinates retain their bound values and have exactly zero velocity.
    /// Initial acceleration is the original force balance divided by the
    /// membrane owner's positive diagonal reference mass.
    pub fn new(
        model: &'m GasMembrane<'p>,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        dead_load_n: &[[f64; 3]],
        controls: GasMembraneMotionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasMembraneMotionError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !finite_vectors(velocity_m_s, model.reference_m.len())
            || model
                .fixed
                .iter()
                .enumerate()
                .any(|(i, fixed)| fixed.is_some() && velocity_m_s[i / 3][i % 3] != 0.0)
        {
            return Err(GasMembraneMotionError::Invalid(
                "finite nodal velocities with stationary prescribed coordinates required",
            ));
        }
        let accepted = model.evaluate(displacement_m, dead_load_n, controls.gas, &mut cancelled)?;
        let masses_kg: Vec<_> = model
            .free_dofs
            .iter()
            .map(|&i| model.patch.nodal_mass_kg()[i / 3])
            .collect();
        if masses_kg.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(GasMembraneMotionError::Invalid(
                "positive native reference mass required",
            ));
        }
        let q: Vec<_> = model
            .free_dofs
            .iter()
            .map(|&i| displacement_m[i / 3][i % 3])
            .collect();
        let v: Vec<_> = model
            .free_dofs
            .iter()
            .map(|&i| velocity_m_s[i / 3][i % 3])
            .collect();
        let a: Vec<_> = model
            .free_dofs
            .iter()
            .zip(&masses_kg)
            .map(|(&i, &m)| -accepted.residual_n[i / 3][i % 3] / m)
            .collect();
        if a.iter().any(|a| !a.is_finite()) {
            return Err(GasMembraneMotionError::Invalid(
                "initial acceleration is unrepresentable",
            ));
        }
        kinetic(&masses_kg, &v)?;
        poll(&mut cancelled)?;
        Ok(Self {
            model,
            state: SecondOrderState::new(0.0, &q, &v, &a),
            accepted,
            masses_kg,
        })
    }

    /// Immutable original geometry, supports and material/gas source binding.
    pub fn model(&self) -> &GasMembrane<'p> {
        self.model
    }

    /// Native free-coordinate checkpoint and accepted convergence history.
    /// q is absolute displacement [m], v is velocity [m/s]. After a nonlinear
    /// step, a is the native algorithmic generalized-alpha acceleration; it is
    /// not an independently force-balanced instantaneous endpoint acceleration.
    pub fn state(&self) -> &SecondOrderState {
        &self.state
    }

    /// Fresh native gas and membrane evaluation at the accepted shape and the
    /// most recent held dead load. The endpoint instantaneous force residual
    /// does not assert static equilibrium or balance with the algorithmic
    /// endpoint acceleration; fixed rows retain support reactions.
    pub fn accepted(&self) -> &GasMembraneEvaluation {
        &self.accepted
    }

    /// Full nodal velocity in the original ordering, with prescribed rows zero.
    pub fn velocity_m_s(&self) -> Vec<[f64; 3]> {
        self.model.expand(&self.state.v, true)
    }

    pub fn kinetic_energy_j(&self) -> Result<f64, GasMembraneMotionError> {
        kinetic(&self.masses_kg, &self.state.v)
    }

    /// Advance through the original matrix-free midpoint owner. The supplied
    /// dead forces are held for this interval; pressure follows the actual
    /// geometry and the source isentrope at every admitted stage and endpoint.
    /// All physical participants, clock and telemetry publish together.
    pub fn step(
        &mut self,
        duration_s: f64,
        dead_load_n: &[[f64; 3]],
        controls: GasMembraneMotionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneMotionStep, GasMembraneMotionError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.state.steps >= controls.maximum_steps {
            return Err(GasMembraneMotionError::StepBudget);
        }
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !finite_vectors(dead_load_n, self.model.reference_m.len())
        {
            return Err(GasMembraneMotionError::Invalid(
                "positive interval and finite nodal dead loads required",
            ));
        }
        let operator = Operator::new(self.model, &self.masses_kg, controls, &mut cancelled);
        let forcing: Vec<_> = self
            .model
            .free_dofs
            .iter()
            .map(|&i| dead_load_n[i / 3][i % 3])
            .collect();
        // Copy only the current mechanical vectors. Previously accepted native
        // telemetry is moved, without a trajectory-sized clone, at publication.
        let mut candidate =
            SecondOrderState::new(self.state.t, &self.state.q, &self.state.v, &self.state.a);
        candidate.steps = self.state.steps;
        let integrator =
            OperatorGeneralizedAlpha::new(forcing.len(), duration_s, 1.0, controls.time);
        let result = integrator.step_controlled(&mut candidate, &operator, &forcing, &mut || {
            operator.cancelled()
        });
        if let Some(error) = operator.take_stop() {
            return Err(error);
        }
        let time = result.map_err(|error| {
            operator.take_failure().map_or_else(
                || GasMembraneMotionError::from(error),
                GasMembraneMotionError::from,
            )
        })?;
        let accepted = operator.final_evaluation(&candidate.q, dead_load_n)?;
        kinetic(&self.masses_kg, &candidate.v)?;
        let mut kinetic_change = Dd::from_f64(0.0);
        let mut dead_work = Dd::from_f64(0.0);
        for (i, &dof) in self.model.free_dofs.iter().enumerate() {
            kinetic_change = kinetic_change
                + Dd::from_f64(
                    (0.5 * candidate.v[i] + 0.5 * self.state.v[i])
                        * self.masses_kg[i]
                        * (candidate.v[i] - self.state.v[i]),
                );
            dead_work = dead_work
                + Dd::from_f64(dead_load_n[dof / 3][dof % 3] * (candidate.q[i] - self.state.q[i]));
        }
        let membrane_change =
            Dd::from_f64(accepted.strain_energy_j) - Dd::from_f64(self.accepted.strain_energy_j);
        let gas_change = Dd::from_f64(accepted.gas_volume.stored_energy_change_j)
            - Dd::from_f64(self.accepted.gas_volume.stored_energy_change_j);
        let swept = self.model.patch.pressure_work_j(
            &self.accepted.displacement_m,
            &accepted.displacement_m,
            1.0,
            &mut || operator.cancelled(),
        )?;
        let exterior_work = self.model.exterior_pressure_pa * swept;
        let defect = (kinetic_change + membrane_change + gas_change + Dd::from_f64(exterior_work)
            - dead_work)
            .to_f64();
        let frame = GasMembraneMotionStep {
            time,
            kinetic_energy_change_j: kinetic_change.to_f64(),
            membrane_energy_change_j: membrane_change.to_f64(),
            gas_energy_change_j: gas_change.to_f64(),
            work_by_gas_j: -gas_change.to_f64(),
            swept_volume_change_m3: swept,
            exterior_pressure_work_j: exterior_work,
            dead_load_work_j: dead_work.to_f64(),
            integration_defect_j: defect,
            evaluations: operator.evaluations(),
            gas_trials: operator.gas_trials(),
        };
        if ![
            frame.kinetic_energy_change_j,
            frame.membrane_energy_change_j,
            frame.gas_energy_change_j,
            swept,
            exterior_work,
            frame.dead_load_work_j,
            defect,
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(GasMembraneMotionError::Invalid(
                "motion energy/work account is unrepresentable",
            ));
        }
        if defect.abs() > controls.energy_tolerance_j {
            return Err(GasMembraneMotionError::Energy {
                defect_j: defect,
                allowance_j: controls.energy_tolerance_j,
            });
        }
        if operator.cancelled() {
            return Err(operator.take_stop().expect("stopped operator"));
        }
        let entry = candidate
            .history
            .pop()
            .expect("native accepted-step telemetry");
        candidate.history = std::mem::take(&mut self.state.history);
        candidate.history.push(entry);
        self.state = candidate;
        self.accepted = accepted;
        Ok(frame)
    }
}

fn kinetic(masses: &[f64], velocity: &[f64]) -> Result<f64, GasMembraneMotionError> {
    let energy: f64 = masses
        .iter()
        .zip(velocity)
        .map(|(&m, &v)| 0.5 * m * v * v)
        .sum();
    if !energy.is_finite() {
        return Err(GasMembraneMotionError::Invalid(
            "kinetic energy is unrepresentable",
        ));
    }
    Ok(energy)
}
