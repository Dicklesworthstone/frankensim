//! Continuing motion with source-owned simultaneous heating and volume work.
mod gas;
mod operator;
#[cfg(test)]
mod tests;

use super::*;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::midpoint::GasHeatVolumeStep;
use operator::ThermalOperator;

/// Independent source controls and bounds on the coupled native proposal.
#[derive(Clone, Copy, Debug)]
pub struct GasMembraneThermalControls {
    pub gas: GasEnergySolve,
    /// Required exactly when the opposing chamber is sealed.
    pub exterior_gas: Option<GasEnergySolve>,
    pub time: ImplicitSolveConfig,
    pub maximum_evaluations: usize,
    pub maximum_steps: usize,
    pub maximum_krylov_restart: usize,
    /// Finite nonnegative enclosing storage-minus-work-minus-heat allowance.
    pub energy_tolerance_j: f64,
    /// Separate finite nonnegative midpoint-traction/volume-work allowance.
    pub pressure_work_tolerance_j: f64,
}
impl GasMembraneThermalControls {
    fn validate(self) -> Result<(), GasMembraneMotionError> {
        self.gas.validate()?;
        if let Some(gas) = self.exterior_gas {
            gas.validate()?;
        }
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
            || ![self.energy_tolerance_j, self.pressure_work_tolerance_j]
                .iter()
                .all(|x| x.is_finite() && *x >= 0.0)
        {
            return Err(GasMembraneMotionError::Invalid(
                "bounded thermal motion controls and nonnegative work/energy allowances required",
            ));
        }
        Ok(())
    }
}

/// Actual accepted thermodynamic states and instantaneous endpoint mechanics.
/// A heated state is not represented by an isentropic work receipt or potential.
#[derive(Clone, Debug, PartialEq)]
pub struct GasMembraneThermalEvaluation {
    pub displacement_m: Vec<[f64; 3]>,
    pub gas: MembraneGas,
    pub exterior_gas: Option<MembraneGas>,
    pub pressure: GasPressureResponse,
    pub exterior_pressure_pa: f64,
    /// Internal minus dead and endpoint pressure loads; fixed rows are reactions.
    /// This is not balance against the native algorithmic endpoint acceleration.
    pub residual_n: Vec<[f64; 3]>,
    pub strain_energy_j: f64,
}

/// Accepted interval accounting; supplied heat is counted once per chamber.
#[derive(Clone, Debug)]
pub struct GasMembraneThermalStep {
    pub time: ImplicitStepTelemetry,
    /// Original heat, work, storage, entropy and discrete derivative receipt.
    pub gas: GasHeatVolumeStep,
    pub exterior_gas: Option<GasHeatVolumeStep>,
    pub midpoint_pressure_pa: f64,
    pub exterior_midpoint_pressure_pa: f64,
    pub kinetic_energy_change_j: f64,
    pub membrane_energy_change_j: f64,
    pub swept_volume_change_m3: f64,
    pub exterior_pressure_work_j: f64,
    pub dead_load_work_j: f64,
    /// Stage-traction work minus native gas work plus held-reservoir work.
    pub pressure_work_defect_j: f64,
    /// Delta K + Delta E_membrane + both Delta U + reservoir work - dead work - heat.
    pub integration_defect_j: f64,
    pub evaluations: usize,
    /// Coupled proposals, each containing one or two native midpoint/endpoint inversions.
    pub gas_trials: usize,
}

/// Existing structural trajectory continued with prescribed signed interval heat.
///
/// The native gas owner uses midpoint in (internal energy, volume), while native
/// generalized-alpha at rho_inf=1 uses midpoint mechanical coordinates. For a
/// curved chart, V(q_mid) need not equal (V_old+V_new)/2 and midpoint traction
/// need not integrate to exact swept-volume work. Both the resulting pressure
/// work mismatch and total numerical energy defect are measured and gated;
/// neither is heat or physical damping. Time refinement remains necessary.
///
/// No isentropic tangent is substituted under heating. Heat is prescribed
/// energy from/to an external owner, not a reservoir-temperature or transport
/// law. Composition, gas mass, geometry/support binding and reference masses
/// are retained. Gas waves, phase change, heat-transfer/source dynamics and
/// moving supports are outside this path.
#[derive(Clone)]
pub struct GasMembraneThermalMotion<'m, 'p> {
    model: &'m GasMembrane<'p>,
    state: SecondOrderState,
    accepted: GasMembraneThermalEvaluation,
    masses_kg: Vec<f64>,
}
impl<'m, 'p> GasMembraneThermalMotion<'m, 'p> {
    pub(super) fn from_motion(motion: GasMembraneMotion<'m, 'p>) -> Self {
        let old = motion.accepted;
        Self {
            model: motion.model,
            state: motion.state,
            masses_kg: motion.masses_kg,
            accepted: GasMembraneThermalEvaluation {
                displacement_m: old.displacement_m,
                gas: old.gas,
                exterior_gas: old.exterior_gas.map(|gas| gas.gas),
                pressure: old.pressure,
                exterior_pressure_pa: old.exterior_pressure_pa,
                residual_n: old.residual_n,
                strain_energy_j: old.strain_energy_j,
            },
        }
    }
    /// Original geometry/support/reference binding. Current gases are in `accepted()`.
    pub fn model(&self) -> &GasMembrane<'p> {
        self.model
    }
    /// Preserved native state and full history; a remains algorithmic after a step.
    pub fn state(&self) -> &SecondOrderState {
        &self.state
    }
    pub fn accepted(&self) -> &GasMembraneThermalEvaluation {
        &self.accepted
    }
    pub fn velocity_m_s(&self) -> Vec<[f64; 3]> {
        self.model.expand(&self.state.v, true)
    }
    pub fn kinetic_energy_j(&self) -> Result<f64, GasMembraneMotionError> {
        kinetic(&self.masses_kg, &self.state.v)
    }

    /// Solve simultaneous native gas heat/work and membrane motion. Positive heat
    /// enters each parcel; heat and dead forces are prescribed for this interval.
    /// `exterior_heat_j` and its independent controls are required for sealed
    /// opposing gas and must both be absent for a fixed-pressure reservoir.
    /// Native midpoint/endpoint source admission, work/energy gates and final
    /// cancellation precede publishing gases, mechanics, clock and history.
    pub fn step(
        &mut self,
        duration_s: f64,
        dead_load_n: &[[f64; 3]],
        heat_j: f64,
        exterior_heat_j: Option<f64>,
        controls: GasMembraneThermalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneThermalStep, GasMembraneMotionError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.state.steps >= controls.maximum_steps {
            return Err(GasMembraneMotionError::StepBudget);
        }
        let sealed = self.accepted.exterior_gas.is_some();
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !heat_j.is_finite()
            || exterior_heat_j.is_some_and(|q| !q.is_finite())
            || exterior_heat_j.is_some() != sealed
            || controls.exterior_gas.is_some() != sealed
            || !finite_vectors(dead_load_n, self.model.reference_m.len())
        {
            return Err(GasMembraneMotionError::Invalid(
                "positive interval, finite loads/heat and matching exterior heat controls required",
            ));
        }
        let operator = ThermalOperator::new(
            self.model,
            &self.masses_kg,
            &self.accepted,
            &self.state.q,
            heat_j,
            exterior_heat_j,
            controls,
            &mut cancelled,
        );
        let forcing: Vec<_> = self
            .model
            .free_dofs
            .iter()
            .map(|&i| dead_load_n[i / 3][i % 3])
            .collect();
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
        let time =
            result.map_err(|error| operator.take_failure().unwrap_or_else(|| error.into()))?;
        let trial = operator.final_trial(&candidate.q)?;
        kinetic(&self.masses_kg, &candidate.v)?;
        let mut kinetic_change = Dd::from_f64(0.0);
        let mut dead_work = Dd::from_f64(0.0);
        let mut stage_volume_work = Dd::from_f64(0.0);
        for (i, &dof) in self.model.free_dofs.iter().enumerate() {
            let dq = candidate.q[i] - self.state.q[i];
            kinetic_change = kinetic_change
                + Dd::from_f64(
                    (0.5 * candidate.v[i] + 0.5 * self.state.v[i])
                        * self.masses_kg[i]
                        * (candidate.v[i] - self.state.v[i]),
                );
            dead_work = dead_work + Dd::from_f64(dead_load_n[dof / 3][dof % 3] * dq);
            stage_volume_work =
                stage_volume_work + Dd::from_f64(trial.stage.gradient[dof / 3][dof % 3] * dq);
        }
        let membrane_change = Dd::from_f64(trial.endpoint.strain_energy_j)
            - Dd::from_f64(self.accepted.strain_energy_j);
        let swept = self.model.patch.pressure_work_j(
            &self.accepted.displacement_m,
            &trial.endpoint.displacement,
            1.0,
            &mut || operator.cancelled(),
        )?;
        let reservoir_work = self.model.exterior.reservoir_work_j(swept);
        let opposing = trial.exterior.as_ref().map(|gas| gas.report);
        let extra_storage = opposing.map_or(0.0, |gas| gas.stored_energy_change_j);
        let extra_work = opposing.map_or(0.0, |gas| gas.work_by_gas_j);
        let pressure_defect = (Dd::from_f64(trial.net_pressure_pa * stage_volume_work.to_f64())
            - Dd::from_f64(trial.gas.report.work_by_gas_j)
            - Dd::from_f64(extra_work)
            + Dd::from_f64(reservoir_work))
        .to_f64();
        let defect = (kinetic_change
            + membrane_change
            + Dd::from_f64(trial.gas.report.stored_energy_change_j)
            + Dd::from_f64(extra_storage)
            + Dd::from_f64(reservoir_work)
            - dead_work
            - Dd::from_f64(heat_j)
            - Dd::from_f64(exterior_heat_j.unwrap_or(0.0)))
        .to_f64();
        let frame = GasMembraneThermalStep {
            time,
            gas: trial.gas.report,
            exterior_gas: opposing,
            midpoint_pressure_pa: trial.gas.midpoint_pressure_pa,
            exterior_midpoint_pressure_pa: trial.exterior_midpoint_pa,
            kinetic_energy_change_j: kinetic_change.to_f64(),
            membrane_energy_change_j: membrane_change.to_f64(),
            swept_volume_change_m3: swept,
            exterior_pressure_work_j: reservoir_work,
            dead_load_work_j: dead_work.to_f64(),
            pressure_work_defect_j: pressure_defect,
            integration_defect_j: defect,
            evaluations: operator.evaluations(),
            gas_trials: operator.gas_trials(),
        };
        if ![
            frame.kinetic_energy_change_j,
            frame.membrane_energy_change_j,
            swept,
            reservoir_work,
            frame.dead_load_work_j,
            pressure_defect,
            defect,
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(GasMembraneMotionError::Invalid(
                "thermal motion work/energy overflow",
            ));
        }
        if pressure_defect.abs() > controls.pressure_work_tolerance_j {
            return Err(GasMembraneMotionError::PressureWork {
                defect_j: pressure_defect,
                allowance_j: controls.pressure_work_tolerance_j,
            });
        }
        if defect.abs() > controls.energy_tolerance_j {
            return Err(GasMembraneMotionError::Energy {
                defect_j: defect,
                allowance_j: controls.energy_tolerance_j,
            });
        }
        let accepted = trial.into_evaluation(dead_load_n)?;
        if operator.cancelled() {
            return Err(operator.take_stop().expect("stopped operator"));
        }
        drop(operator);
        let entry = candidate.history.pop().expect("native accepted telemetry");
        candidate.history = std::mem::take(&mut self.state.history);
        candidate.history.push(entry);
        self.state = candidate;
        self.accepted = accepted;
        Ok(frame)
    }
}
