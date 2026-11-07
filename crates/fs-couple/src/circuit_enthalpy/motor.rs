//! A current-regulated motor and finite thermal bodies on one accepted timeline.
//!
//! Winding resistance comes from the original material at each proposed thermal
//! endpoint. The native reciprocal rotor drive then enforces current, voltage
//! and source draw/return limits with that resistance. Its independently paired
//! winding and bridge heats enter their explicit cells exactly once. Every
//! thermal trial restarts the same accepted motion; the original cached rotor
//! factorization is reused, including under thermal feedback and curtailment.
//!
//! The combined account encloses rotor KINETIC energy and these caloric bodies.
//! Original stiffness, damping and gyroscopic work are signed exchanges across
//! that boundary. In particular, stiffness work is not presumed to be elastic
//! storage and active-bearing work is not called heat. Shaft work is internal
//! motor/rotor exchange. Electrical input and external loads remain independent.
//!
//! This is the native zero-inductance, held-current motor with fixed reciprocal
//! torque constant, bridge efficiencies, geometry and prepared M/C/K/G/spin.
//! Thermal rates and material coefficients use the endpoint, hence first-order
//! coupled accuracy in general. It does not infer magnetic saturation, winding
//! inductance, switching loss, temperature-dependent magnets, bearing heat
//! routing or a temperature protection law. Material-domain and native drive
//! refusals publish neither history; refinement is still needed for time error.

/// Finite winding-current and magnetic storage on the same material/enthalpy residual.
pub mod voltage;
/// Complete source-owned contact, radiation, layer and coolant heat transfer.
pub mod heat;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate, finite,
    magnetic::same_material,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use crate::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductor},
    state_point::MaterialPropertySelection,
    thermal::source::ThermalMaterialError,
};
use fs_math::dd::Dd;
use fs_mbd::rotordynamics::{
    MAX_ROTOR_DOFS, RotorError,
    transient::{
        RotorMidpoint, RotorPowerError, RotorTransientState,
        drive::{
            RotorCurrentDrive, RotorCurrentDriveStep, RotorDriveError, RotorDrivePowerControl,
        },
    },
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// Whole-cell uniform winding conductor, R=rho(T,phase)*length/area.
/// Length times area must equal the receiving cell's physical volume.
#[derive(Clone, Debug)]
pub struct MotorWinding {
    pub cell: CellId,
    pub basis: ElectricalPropertyBasis,
    pub selection: MaterialPropertySelection,
    pub length_m: f64,
    pub area_m2: f64,
}

/// Original caloric/transport inventory and explicit winding/bridge receivers.
/// Both channels may share one cell only when that homogenization is authored.
/// Passive thermal cells, conduction and reservoirs retain their native owners.
pub struct MotorThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub winding: MotorWinding,
    pub bridge_cell: CellId,
}

/// Accepted motion and actual winding source, borrowing one prepared rotor.
/// No electrical state is invented for the native held-current approximation.
#[derive(Clone, Debug)]
pub struct MaterialRotorCurrentDrive<'r> {
    native: &'r RotorMidpoint<'r>,
    time_s: f64,
    state: RotorTransientState,
    parameters: RotorCurrentDrive,
    winding: Arc<ResolvedConductor>,
    winding_cell: CellId,
    bridge_cell: CellId,
}
impl MaterialRotorCurrentDrive<'_> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn state(&self) -> &RotorTransientState {
        &self.state
    }
    /// Original authored drive with its actual accepted material resistance.
    pub const fn parameters(&self) -> RotorCurrentDrive {
        self.parameters
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

#[derive(Clone, Debug)]
pub struct MotorThermalSample<S = MaterialNetworkSample> {
    /// Complete original transport results at the accepted thermal endpoint.
    pub network: S,
    pub winding: Arc<ResolvedConductor>,
    pub winding_cell: CellId,
    pub bridge_cell: CellId,
}

/// Original drive and the explicit kinetic/thermal accounting boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct MotorThermalStep {
    pub drive: RotorCurrentDriveStep,
    /// Signed native C, K and gyroscopic work, exported from this boundary [J].
    /// Can be negative for active/nonsymmetric operators; never deposited as heat.
    pub operator_work_j: f64,
    /// Electrical draw minus return plus external load work minus operator work
    /// [J], evaluated from those independent native pairings, not from storage.
    pub supplied_to_kinetic_and_heat_j: f64,
}

pub struct MotorEnthalpyCandidate<'t, 'r, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    motor: MaterialRotorCurrentDrive<'r>,
    report: CircuitEnthalpyReport<MotorThermalSample<S>, MotorThermalStep>,
}
impl<'t, 'r, S> MotorEnthalpyCandidate<'t, 'r, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn motor(&self) -> &MaterialRotorCurrentDrive<'r> {
        &self.motor
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<MotorThermalSample<S>, MotorThermalStep> {
        &self.report
    }
    /// Publish motion and thermal history together. Winding and bridge heat in
    /// the native drive report are already deposited and must not be added again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialRotorCurrentDrive<'r>,
        CircuitEnthalpyReport<MotorThermalSample<S>, MotorThermalStep>,
    ) {
        (self.thermal, self.motor, self.report)
    }
}

#[derive(Debug)]
pub enum MotorThermalError {
    Material(MaterialNetworkError),
    Source(ThermalMaterialError),
    Conductor(JouleHeatingError),
    Drive(RotorDriveError),
    VoltageDrive(fs_mbd::rotordynamics::transient::drive::voltage::RotorVoltageError),
    Invalid(&'static str),
}
impl core::fmt::Display for MotorThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rotor motor material/enthalpy: {self:?}")
    }
}
impl std::error::Error for MotorThermalError {}
type Error = CircuitEnthalpyError<MotorThermalError>;
fn motor_error<E: From<MotorThermalError>>(error: Error) -> CircuitEnthalpyError<E> {
    error.map_source(E::from)
}
fn invalid(reason: &'static str) -> Error {
    CircuitEnthalpyError::Source(MotorThermalError::Invalid(reason))
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancel.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(MotorThermalError::Material(other)),
    }
}
fn drive_error(error: RotorDriveError) -> Error {
    match error {
        RotorDriveError::Power(RotorPowerError::Dynamics(RotorError::Cancelled)) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(MotorThermalError::Drive(other)),
    }
}
fn source(error: ThermalMaterialError) -> Error {
    CircuitEnthalpyError::Source(MotorThermalError::Source(error))
}
fn scalar(value: Dd) -> Result<f64, Error> {
    finite(value.hi)?;
    finite(value.lo)?;
    Ok(finite(value.to_f64())?)
}
fn same_conductor(a: &ResolvedConductor, b: &ResolvedConductor) -> bool {
    a.source_basis() == b.source_basis()
        && a.source_key() == b.source_key()
        && a.length_m().to_bits() == b.length_m().to_bits()
        && a.area_m2().to_bits() == b.area_m2().to_bits()
        && a.resistance_ohm().to_bits() == b.resistance_ohm().to_bits()
        && same_material(a.material(), b.material())
}
fn clock_spacing(time: f64) -> f64 {
    let magnitude = time.abs();
    let next = magnitude.next_up();
    if next.is_finite() {
        next - magnitude
    } else {
        magnitude - magnitude.next_down()
    }
}

struct MotorDriver<'r> {
    motor: MaterialRotorCurrentDrive<'r>,
    step: MotorThermalStep,
    heat_j: f64,
    supplied_j: Dd,
}
impl<'r> ThermalDriver for MotorDriver<'r> {
    type Model = MaterialRotorCurrentDrive<'r>;
    type Step = MotorThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.drive.motion.kinetic_energy_change_j,
            supplied_j: self.supplied_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.drive.motion.kinetic_energy_change_j) - self.supplied_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.motor, self.step)
    }
}

// Each physical balance is reconstructed from original work/heat contractions.
// The heat deposition is never an electrical or mechanical residual remainder.
fn physical_account(
    drive: RotorCurrentDriveStep,
    tolerance_j: f64,
) -> Result<(MotorThermalStep, f64, Dd), Error> {
    for heat in [
        drive.winding_heat_j,
        drive.bridge_heat_j,
        drive.electrical_draw_j,
        drive.electrical_return_j,
    ] {
        if !heat.is_finite() || heat < 0.0 {
            return Err(invalid(
                "motor heat and directional source energy must be finite and nonnegative",
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
    let load = Dd::from_f64(finite(motion.input_work_j)?);
    let supplied = electrical + external - operator;
    for residual in [
        kinetic + operator - load,
        load - external - shaft,
        electrical - shaft - heat,
        kinetic + heat - supplied,
    ] {
        if scalar(residual)?.abs() > tolerance_j {
            return Err(invalid(
                "native motor, load and kinetic/thermal work accounts must close independently",
            ));
        }
    }
    let step = MotorThermalStep {
        drive,
        operator_work_j: scalar(operator)?,
        supplied_to_kinetic_and_heat_j: scalar(supplied)?,
    };
    Ok((step, scalar(heat)?, supplied))
}

impl MotorThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        }
    }

    fn resolve_winding(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        original: Option<&ResolvedConductor>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ResolvedConductor, Error> {
        check(cancel)?;
        // A bridge remains a physical heat receiver even at ideal efficiency.
        for id in [self.winding.cell, self.bridge_cell] {
            if !thermal.cells().contains_key(&id) || !self.cells.iter().any(|c| c.cell == id) {
                return Err(invalid(
                    "winding and bridge heat require active original material cells",
                ));
            }
        }
        let body = self
            .cells
            .iter()
            .find(|c| c.cell == self.winding.cell)
            .ok_or_else(|| invalid("winding lacks its original material binding"))?;
        let phase = sample
            .cells
            .iter()
            .find(|c| c.cell == self.winding.cell)
            .ok_or_else(|| invalid("winding lacks its actual caloric phase"))?
            .phase;
        let volume = thermal.cells()[&self.winding.cell].geometry.volume_m3;
        let ohmic = MaterialOhmicBody {
            material: body.point,
            selection: self.winding.selection.clone(),
            length_m: self.winding.length_m,
            area_m2: self.winding.area_m2,
        };
        let conductor = if let Some(original) = original {
            if original.source_basis() != self.winding.basis
                || original.length_m().to_bits() != self.winding.length_m.to_bits()
                || original.area_m2().to_bits() != self.winding.area_m2.to_bits()
            {
                return Err(invalid(
                    "accepted winding geometry and electrical source basis cannot change",
                ));
            }
            let point = body.point.at(phase).map_err(source)?;
            let snapshot = original
                .material()
                .at_query_point(body.point.material, &point)
                .map_err(|e| source(ThermalMaterialError::Material(e)))?;
            ohmic.conductor_from_material(phase, volume, snapshot, self.winding.basis)
        } else {
            ohmic.resolve_conductor_with_basis(phase, volume, self.winding.basis)
        }
        .map_err(|e| CircuitEnthalpyError::Source(MotorThermalError::Conductor(e)))?;
        check(cancel)?;
        Ok(conductor)
    }

    /// Attach actual winding material to finite motion without advancing time.
    /// The material replaces the drive's numeric winding-resistance placeholder.
    /// Finite original coordinates are retained; complete rotor/drive admission
    /// remains with the native stepper on every candidate. The borrowed prepared
    /// rotor fixes model/spin and duration for this history, avoiding silent
    /// mechanical coefficient changes or a factorization inside thermal trials.
    pub fn bind<'r>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &'r RotorMidpoint<'r>,
        state: RotorTransientState,
        parameters: RotorCurrentDrive,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialRotorCurrentDrive<'r>, Error> {
        let cancel = RefCell::new(cancelled);
        self.bind_resolved(thermal, native, state, parameters, || {
            let network = self.network();
            let prepared = network
                .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
                .map_err(material)?;
            let sample = prepared
                .resolve(thermal, || (*cancel.borrow_mut())())
                .map_err(material)?;
            self.resolve_winding(thermal, &sample.source, None, &cancel)
        }, &cancel)
    }

    // Both transport frontends bind the same native motion and material winding.
    fn bind_resolved<'r, E: From<MotorThermalError>>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &'r RotorMidpoint<'r>,
        state: RotorTransientState,
        mut parameters: RotorCurrentDrive,
        resolve: impl FnOnce() -> Result<ResolvedConductor, CircuitEnthalpyError<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MaterialRotorCurrentDrive<'r>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(motor_error::<E>)?;
        if state.velocity.is_empty()
            || state.velocity.len() > MAX_ROTOR_DOFS
            || state.displacement.len() != state.velocity.len()
            || parameters.shaft_coordinate >= state.velocity.len()
            || state
                .displacement
                .iter()
                .chain(&state.velocity)
                .any(|x| !x.is_finite())
        {
            return Err(motor_error(invalid(
                "motor must retain finite original rotor coordinates and a valid shaft index",
            )));
        }
        let winding = Arc::new(resolve()?);
        parameters.winding_resistance_ohm = winding.resistance_ohm();
        check(cancel).map_err(motor_error::<E>)?;
        Ok(MaterialRotorCurrentDrive {
            native,
            time_s: thermal.time_s(),
            state,
            parameters,
            winding,
            winding_cell: self.winding.cell,
            bridge_cell: self.bridge_cell,
        })
    }

    /// Solve native current/voltage/source-limited actuation and finite heat
    /// together. R(T,phase), torque curtailment, winding/bridge losses, thermal
    /// conduction and reservoirs participate in every endpoint trial. No frozen
    /// resistance, frozen selected current or post-acceptance heat update occurs.
    ///
    /// The interval must match the borrowed native duration, within represented
    /// clock roundoff. Each native motor solve consumes one circuit-trial unit.
    /// Original source claims remain pinned; events and source discontinuities
    /// require interval cuts as for the common enthalpy driver. No small-update,
    /// source-domain or nonlinear-failure fallback can publish partial history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, 'r>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialRotorCurrentDrive<'r>,
        external_force: &[f64],
        commanded_current_a: f64,
        power: RotorDrivePowerControl,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MotorEnthalpyCandidate<'t, 'r>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
            thermal, target_s, motor, external_force, commanded_current_a, power,
            controls, initial_energy_increments_j,
            |state, original| {
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let winding = self.resolve_winding(state, &coefficients.source, original, &cancel)?;
                Ok((coefficients, winding))
            },
            &cancel,
        )
    }

    // A source callback changes only transport composition. Every trial still
    // restarts this original native motor and uses one physical heat deposition.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, 'r, S, E: From<MotorThermalError>>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialRotorCurrentDrive<'r>,
        external_force: &[f64],
        commanded_current_a: f64,
        power: RotorDrivePowerControl,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>, Option<&ResolvedConductor>)
            -> Result<(CircuitThermalCoefficients<S>, ResolvedConductor), CircuitEnthalpyError<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MotorEnthalpyCandidate<'t, 'r, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(motor_error::<E>)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        let h = motor.native.duration_seconds();
        let represented_end = finite(thermal.time_s() + h)?;
        let clock_roundoff = 0.5 * (clock_spacing(thermal.time_s()) + clock_spacing(target_s));
        if motor.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || represented_end <= thermal.time_s()
            || (h - dt).abs() > clock_roundoff + 64.0 * f64::EPSILON * h.max(dt)
            || motor.winding_cell != self.winding.cell
            || motor.bridge_cell != self.bridge_cell
        {
            return Err(motor_error(invalid(
                "motor and heat require unchanged receivers and the same positive native interval",
            )));
        }
        let (_, original) = resolve(thermal, None)?;
        if !same_conductor(&original, &motor.winding) {
            return Err(motor_error(invalid(
                "winding source, geometry and resistance must match the accepted thermal state",
            )));
        }
        let (thermal, motor, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _remaining| {
                check(cancel).map_err(motor_error::<E>)?;
                let (coefficients, winding) = resolve(trial, Some(&motor.winding))?;
                let winding = Arc::new(winding);
                let parameters = RotorCurrentDrive {
                    winding_resistance_ohm: winding.resistance_ohm(),
                    ..motor.parameters
                };
                let drive = motor
                    .native
                    .step_current_drive(
                        &motor.state,
                        external_force,
                        commanded_current_a,
                        parameters,
                        power,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(drive_error).map_err(motor_error::<E>)?;
                let (step, heat_j, supplied_j) =
                    physical_account(drive, controls.total_energy_tolerance_j)
                        .map_err(motor_error::<E>)?;
                let next = MaterialRotorCurrentDrive {
                    native: motor.native,
                    time_s: target_s,
                    state: step.drive.motion.state_after.clone(),
                    parameters,
                    winding: Arc::clone(&winding),
                    winding_cell: motor.winding_cell,
                    bridge_cell: motor.bridge_cell,
                };
                check(cancel).map_err(motor_error::<E>)?;
                Ok(DrivenCoefficients {
                    circuit: MotorDriver {
                        motor: next,
                        step,
                        heat_j,
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
            |candidate, duration, checkpoint| receiver_powers(
                [
                    (candidate.motor.winding_cell, candidate.step.drive.winding_heat_j),
                    (candidate.motor.bridge_cell, candidate.step.drive.bridge_heat_j),
                ],
                duration, checkpoint,
            ),
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(motor_error::<E>)?;
        Ok(MotorEnthalpyCandidate {
            thermal,
            motor,
            report,
        })
    }
}


// The two drive approximations share exactly-once receiver aggregation, including
// an explicitly homogenized winding/bridge cell. No residual supplies this heat.
fn receiver_powers<E: From<MotorThermalError>>(
    receivers: [(CellId, f64); 2],
    duration: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<E>> {
    let mut heat_by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
    for (cell, heat) in receivers {
        checkpoint()?;
        let total = heat_by_cell.entry(cell).or_insert(Dd::ZERO);
        *total = *total + Dd::from_f64(heat);
    }
    let mut powers = BTreeMap::new();
    for (cell, heat) in heat_by_cell {
        checkpoint()?;
        let heat = scalar(heat).map_err(motor_error::<E>)?;
        let power = finite(heat / duration)?;
        if heat > 0.0 && power == 0.0 {
            return Err(motor_error(invalid("motor heat is unrepresentable as a thermal rate")));
        }
        powers.insert(cell, power);
    }
    Ok(powers)
}
