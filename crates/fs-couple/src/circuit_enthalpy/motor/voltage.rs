//! Finite-inductance motor, continuing rotor motion and actual thermal feedback.
//!
//! Every endpoint trial re-evaluates the original winding R(T,phase), rebinds
//! only that resistance on the cached native voltage motor and restarts the
//! SAME accepted current and rotor coordinates. Current/voltage/source limits,
//! winding magnetic storage, reciprocal back-EMF/torque and source regeneration
//! therefore remain in the coupled solve. No mechanical solve or factorization
//! is needed to rebind R; thermal callbacks never advance accepted history.
//!
//! The closed boundary contains kinetic, magnetic and caloric storage. Signed
//! native C/K/G work remains an explicit external operator exchange; shaft work
//! cancels internally. Both winding and bridge heat enter their authored cells
//! exactly once. L, reciprocal Kt and all mechanical coefficients remain fixed,
//! so no unaccounted magnetic or mechanical parameter work is introduced.
//!
//! Winding heat is the native MIDPOINT quadrature R I_mid² dt, not exact RMS
//! heat. Thermal rates and material coefficients use the endpoint. Discrete
//! energy closure does not establish temporal accuracy or measured material
//! validity; time refinement and original source domains remain mandatory.

use super::*;
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorVoltageError, RotorVoltageMotor, RotorVoltageStep, RotorWindingMotor,
    RotorWindingState,
};

/// Accepted current/rotor state with its actual source-bound winding resistance.
/// The owned native cache retains the original stepper and model lifetimes.
#[derive(Clone, Debug)]
pub struct MaterialRotorVoltageDrive<'step, 'model> {
    native: RotorVoltageMotor<'step, 'model>,
    time_s: f64,
    state: RotorWindingState,
    winding: Arc<ResolvedConductor>,
    winding_cell: CellId,
    bridge_cell: CellId,
}
impl<'step, 'model> MaterialRotorVoltageDrive<'step, 'model> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn state(&self) -> &RotorWindingState {
        &self.state
    }
    pub fn parameters(&self) -> RotorWindingMotor {
        self.native.motor()
    }
    pub fn winding(&self) -> &ResolvedConductor {
        &self.winding
    }
    pub const fn winding_cell(&self) -> CellId {
        self.winding_cell
    }
    pub const fn bridge_cell(&self) -> CellId {
        self.bridge_cell
    }
}

/// Native winding/motion accounts with the complete kinetic/magnetic boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorVoltageThermalStep {
    pub drive: RotorVoltageStep,
    /// Signed native C/K/G work exported from this boundary [J], never heat.
    pub operator_work_j: f64,
    /// Electrical net work + external loads - operator work [J], evaluated
    /// independently of kinetic, magnetic and caloric storage changes.
    pub supplied_to_stored_and_heat_j: f64,
}

pub struct RotorVoltageEnthalpyCandidate<'t, 'step, 'model> {
    thermal: DepositionSimulation<'t>,
    motor: MaterialRotorVoltageDrive<'step, 'model>,
    report: CircuitEnthalpyReport<MotorThermalSample, RotorVoltageThermalStep>,
}
impl<'t, 'step, 'model> RotorVoltageEnthalpyCandidate<'t, 'step, 'model> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn motor(&self) -> &MaterialRotorVoltageDrive<'step, 'model> {
        &self.motor
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<MotorThermalSample, RotorVoltageThermalStep> {
        &self.report
    }
    /// Publish current, rotor and heat histories together. Both reported loss
    /// channels have already entered their receiving cells.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialRotorVoltageDrive<'step, 'model>,
        CircuitEnthalpyReport<MotorThermalSample, RotorVoltageThermalStep>,
    ) {
        (self.thermal, self.motor, self.report)
    }
}

fn voltage_error(error: RotorVoltageError) -> Error {
    match error {
        RotorVoltageError::Drive(RotorDriveError::Power(RotorPowerError::Dynamics(
            RotorError::Cancelled,
        ))) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MotorThermalError::VoltageDrive(other)),
    }
}

struct VoltageDriver<'step, 'model> {
    motor: MaterialRotorVoltageDrive<'step, 'model>,
    step: RotorVoltageThermalStep,
    heat_j: f64,
    stored_change_j: Dd,
    supplied_j: Dd,
}
impl<'step, 'model> ThermalDriver for VoltageDriver<'step, 'model> {
    type Model = MaterialRotorVoltageDrive<'step, 'model>;
    type Step = RotorVoltageThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.stored_change_j.to_f64(),
            supplied_j: self.supplied_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        self.stored_change_j - self.supplied_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.motor, self.step)
    }
}

fn physical_account(
    drive: RotorVoltageStep,
    tolerance_j: f64,
) -> Result<(RotorVoltageThermalStep, f64, Dd, Dd), Error> {
    for heat in [
        drive.winding_heat_j,
        drive.bridge_heat_j,
        drive.electrical_draw_j,
        drive.electrical_return_j,
    ] {
        if !heat.is_finite() || heat < 0.0 {
            return Err(invalid(
                "inductive motor heat and directional work must be finite and nonnegative",
            ));
        }
    }
    let heat = Dd::from_f64(drive.winding_heat_j) + Dd::from_f64(drive.bridge_heat_j);
    let electrical =
        Dd::from_f64(drive.electrical_draw_j) - Dd::from_f64(drive.electrical_return_j);
    let external = Dd::from_f64(finite(drive.external_work_j)?);
    let shaft = Dd::from_f64(finite(drive.shaft_work_j)?);
    let motion = &drive.motion;
    let operator = Dd::from_f64(finite(motion.velocity_force_work_j)?)
        + Dd::from_f64(finite(motion.stiffness_force_work_j)?)
        + Dd::from_f64(finite(motion.gyroscopic_work_j)?);
    let kinetic = Dd::from_f64(finite(motion.kinetic_energy_change_j)?);
    let magnetic = Dd::from_f64(finite(drive.magnetic_energy_change_j)?);
    let load = Dd::from_f64(finite(motion.input_work_j)?);
    let stored = kinetic + magnetic;
    let supplied = electrical + external - operator;
    scalar(stored)?;
    for residual in [
        kinetic + operator - load,
        load - external - shaft,
        electrical - magnetic - shaft - heat,
        stored + heat - supplied,
    ] {
        if scalar(residual)?.abs() > tolerance_j {
            return Err(invalid(
                "inductive motor magnetic, shaft and thermal work must close independently",
            ));
        }
    }
    let step = RotorVoltageThermalStep {
        drive,
        operator_work_j: scalar(operator)?,
        supplied_to_stored_and_heat_j: scalar(supplied)?,
    };
    Ok((step, scalar(heat)?, stored, supplied))
}

impl MotorThermalNetwork<'_> {
    /// Bind the actual winding source without resetting stored current or
    /// rotor motion. Only resistance changes; original L/Kt and all mechanical
    /// cache data remain intact. Complete state/current/speed admission stays
    /// with the native candidate, as for the current-driven binding.
    pub fn bind_voltage<'step, 'model>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &RotorVoltageMotor<'step, 'model>,
        state: RotorWindingState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialRotorVoltageDrive<'step, 'model>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let parameters = native.motor();
        if !state.current_a.is_finite()
            || state.rotor.velocity.is_empty()
            || state.rotor.velocity.len() > MAX_ROTOR_DOFS
            || state.rotor.displacement.len() != state.rotor.velocity.len()
            || parameters.drive.shaft_coordinate >= state.rotor.velocity.len()
            || state
                .rotor
                .displacement
                .iter()
                .chain(&state.rotor.velocity)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "inductive motor needs finite continuing current and original rotor coordinates",
            ));
        }
        scalar(
            Dd::from_f64(parameters.inductance_h)
                * Dd::from_f64(0.5 * state.current_a)
                * Dd::from_f64(state.current_a),
        )?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let sample = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let winding = Arc::new(self.resolve_winding(thermal, &sample.source, None, &cancel)?);
        let native = native
            .with_winding_resistance(winding.resistance_ohm(), || (*cancel.borrow_mut())())
            .map_err(voltage_error)?;
        check(&cancel)?;
        Ok(MaterialRotorVoltageDrive {
            native,
            time_s: thermal.time_s(),
            state,
            winding,
            winding_cell: self.winding.cell,
            bridge_cell: self.bridge_cell,
        })
    }

    /// Close finite current dynamics, magnetic/kinetic storage and both loss
    /// receivers against the original finite enthalpy residual. Each trial
    /// rebinds only R(T,phase), then uses native voltage/current/source control.
    /// Circuit-trial units count actual native motor solves. The prepared
    /// physical duration must span the same thermal interval; source, event,
    /// budget and cancellation refusals leave all accepted states unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_candidate<'t, 'step, 'model>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialRotorVoltageDrive<'step, 'model>,
        external_force: &[f64],
        commanded_voltage_v: f64,
        power: RotorVoltageControl,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorVoltageEnthalpyCandidate<'t, 'step, 'model>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        let h = motor.native.duration_s();
        let represented_end = finite(thermal.time_s() + h)?;
        let clock_roundoff = 0.5 * (clock_spacing(thermal.time_s()) + clock_spacing(target_s));
        if motor.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || represented_end <= thermal.time_s()
            || (h - dt).abs() > clock_roundoff + 64.0 * f64::EPSILON * h.max(dt)
            || motor.winding_cell != self.winding.cell
            || motor.bridge_cell != self.bridge_cell
        {
            return Err(invalid(
                "inductive motor and heat require unchanged receivers and the same native interval",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let original = self.resolve_winding(thermal, &accepted.source, None, &cancel)?;
        if !same_conductor(&original, &motor.winding) {
            return Err(invalid(
                "inductive winding source and resistance must match accepted enthalpy",
            ));
        }
        let (thermal, motor, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let winding = Arc::new(self.resolve_winding(
                    trial,
                    &coefficients.source,
                    Some(&motor.winding),
                    &cancel,
                )?);
                let native = motor
                    .native
                    .with_winding_resistance(winding.resistance_ohm(), || (*cancel.borrow_mut())())
                    .map_err(voltage_error)?;
                let drive = native
                    .step(
                        &motor.state,
                        external_force,
                        commanded_voltage_v,
                        power,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(voltage_error)?;
                let (step, heat_j, stored_change_j, supplied_j) =
                    physical_account(drive, controls.total_energy_tolerance_j)?;
                let state = RotorWindingState {
                    rotor: step.drive.motion.state_after.clone(),
                    current_a: step.drive.winding_current_after_a,
                };
                check(&cancel)?;
                Ok(DrivenCoefficients {
                    circuit: VoltageDriver {
                        motor: MaterialRotorVoltageDrive {
                            native,
                            time_s: target_s,
                            state,
                            winding: Arc::clone(&winding),
                            winding_cell: motor.winding_cell,
                            bridge_cell: motor.bridge_cell,
                        },
                        step,
                        heat_j,
                        stored_change_j,
                        supplied_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: MotorThermalSample {
                        network: coefficients.source,
                        winding,
                        winding_cell: self.winding.cell,
                        bridge_cell: self.bridge_cell,
                    },
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                let mut heat_by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (cell, heat) in [
                    (
                        candidate.motor.winding_cell,
                        candidate.step.drive.winding_heat_j,
                    ),
                    (
                        candidate.motor.bridge_cell,
                        candidate.step.drive.bridge_heat_j,
                    ),
                ] {
                    checkpoint()?;
                    let total = heat_by_cell.entry(cell).or_insert(Dd::ZERO);
                    *total = *total + Dd::from_f64(heat);
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in heat_by_cell {
                    checkpoint()?;
                    let heat = scalar(heat)?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "inductive motor heat is unrepresentable as a thermal rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(RotorVoltageEnthalpyCandidate {
            thermal,
            motor,
            report,
        })
    }
}
