//! Original squirrel-cage induction and finite, source-bound winding heat.
//!
//! Each thermal endpoint resolves the original stator and rotor conductor
//! sources. Their actual resistances enter the native seven-state machine,
//! which restarts the same accepted fluxes, angle and momenta on every trial.
//! The native stator/rotor copper and bearing heat enter their declared cells
//! once through the existing circuit/enthalpy residual. No torque or electrical
//! loss is fitted from the thermal balance, and no history advances in a trial.
//!
//! Magnetic inductances, shaft inertias, pole pairs, bearing drag and geometry
//! remain fixed. Only resistance changes, introducing no magnetic parameter
//! work. Material and thermal rates use the endpoint, so the coupled method is
//! generally first order even though the native frozen-coefficient machine
//! step is symmetric. Original source domains, thermal events and time-error
//! refinement remain the caller's responsibility. This adds no skin effect,
//! saturation, cage geometry, insulation failure or experimental validation.

use super::super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate, finite,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use super::{
    MotorThermalError, MotorThermalNetwork, MotorWinding, receiver_powers, same_conductor,
};
use fs_material::conductor::ResolvedConductor;
use fs_math::dd::Dd;
use fs_phs::circuit::induction::{
    InductionMotor, InductionMotorParameters, InductionOperatingPoint, InductionStep,
};
use fs_phs::{PhsError, PreparedStepError, SolveBudget};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, sync::Arc};

/// One whole-cell conductor source and its explicitly authored winding basis.
#[derive(Clone, Debug)]
pub struct InductionWinding {
    /// Actual material, path geometry and finite receiving volume.
    pub conductor: MotorWinding,
    /// Positive fixed conversion from this whole-conductor R to the native
    /// power-invariant winding R. For three equal stator phases represented by
    /// their total series wire length and physical wire area, use 1/3. A cage
    /// referral/connection factor must come from the authored machine model;
    /// no turns ratio, bar network or phase factor is inferred from volume.
    /// The native winding loss already sums the two power-invariant axes and
    /// must NOT be multiplied by this factor again during thermal deposition.
    pub resistance_multiplier: f64,
}

/// Distinct finite stator/rotor bodies with an explicit bearing-heat receiver.
/// The latter may be either winding cell or another bound thermal cell; all
/// bearing heat goes there. This is an authored partition, not an inferred
/// contact or coolant model. Passive cells, links and reservoirs retain their
/// existing material/enthalpy owners.
pub struct InductionThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub stator: InductionWinding,
    pub rotor: InductionWinding,
    pub bearing_cell: CellId,
}

/// Original immutable machine history and its actual accepted material sources.
#[derive(Clone, Debug)]
pub struct MaterialInductionMotor {
    time_s: f64,
    state: [f64; 7],
    parameters: InductionMotorParameters,
    operating_point: InductionOperatingPoint,
    windings: [Arc<ResolvedConductor>; 2],
    winding_cells: [CellId; 2],
    resistance_multipliers: [f64; 2],
    bearing_cell: CellId,
}
impl MaterialInductionMotor {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn state(&self) -> &[f64; 7] {
        &self.state
    }
    pub const fn parameters(&self) -> InductionMotorParameters {
        self.parameters
    }
    pub const fn operating_point(&self) -> &InductionOperatingPoint {
        &self.operating_point
    }
    pub fn stator_winding(&self) -> &ResolvedConductor {
        &self.windings[0]
    }
    pub fn rotor_winding(&self) -> &ResolvedConductor {
        &self.windings[1]
    }
    /// Stator then rotor, independent of the thermal cell's sorted ID order.
    pub const fn winding_cells(&self) -> [CellId; 2] {
        self.winding_cells
    }
    pub const fn resistance_multipliers(&self) -> [f64; 2] {
        self.resistance_multipliers
    }
    pub const fn bearing_cell(&self) -> CellId {
        self.bearing_cell
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InductionEnthalpySolve {
    /// Shared finite-enthalpy Newton, trial, per-cell and total-work budgets.
    pub thermal: CircuitEnthalpySolve,
    /// Native seven-state induction-machine residual/update admission.
    pub machine: SolveBudget,
    pub backtracking_halvings: u8,
}

#[derive(Clone, Debug)]
pub struct InductionThermalSample {
    pub network: MaterialNetworkSample,
    pub stator: Arc<ResolvedConductor>,
    pub rotor: Arc<ResolvedConductor>,
    /// Actual native stator and referred rotor coefficients after conversion.
    pub resistance_ohm: [f64; 2],
    pub winding_cells: [CellId; 2],
    pub bearing_cell: CellId,
}

pub struct InductionEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    motor: MaterialInductionMotor,
    report: CircuitEnthalpyReport<InductionThermalSample, InductionStep>,
}
impl<'t> InductionEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn motor(&self) -> &MaterialInductionMotor {
        &self.motor
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<InductionThermalSample, InductionStep> {
        &self.report
    }
    /// Publish the original machine and caloric histories together. Every
    /// native loss in the report is already in the receiving thermal cells.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialInductionMotor,
        CircuitEnthalpyReport<InductionThermalSample, InductionStep>,
    ) {
        (self.thermal, self.motor, self.report)
    }
}

#[derive(Debug)]
pub enum InductionThermalError {
    Material(MotorThermalError),
    Native(PreparedStepError),
    Invalid(&'static str),
}
impl From<MotorThermalError> for InductionThermalError {
    fn from(error: MotorThermalError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for InductionThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "induction machine material/enthalpy: {self:?}")
    }
}
impl std::error::Error for InductionThermalError {}
type Error = CircuitEnthalpyError<InductionThermalError>;

fn invalid(reason: &'static str) -> Error {
    CircuitEnthalpyError::Source(InductionThermalError::Invalid(reason))
}
fn native(error: PreparedStepError) -> Error {
    match error {
        PreparedStepError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(InductionThermalError::Native(other)),
    }
}
fn phs(error: PhsError) -> Error {
    native(error.into())
}
fn material(error: MaterialNetworkError) -> Error {
    super::material(error).map_source(InductionThermalError::Material)
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancel.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn scalar(value: Dd) -> Result<f64, Error> {
    finite(value.hi)?;
    finite(value.lo)?;
    Ok(finite(value.to_f64())?)
}

struct InductionDriver {
    motor: MaterialInductionMotor,
    step: InductionStep,
    heat_j: f64,
    supplied_j: Dd,
}
impl ThermalDriver for InductionDriver {
    type Model = MaterialInductionMotor;
    type Step = InductionStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.delta_h,
            supplied_j: self.supplied_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.native.delta_h) - self.supplied_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.motor, self.step)
    }
}

// Paired source/shaft work and original copper/bearing losses are independent
// native outputs. None is recovered from the combined thermal residual.
fn physical_account(step: &InductionStep, tolerance_j: f64) -> Result<(f64, Dd), Error> {
    let mut heat = Dd::ZERO;
    for value in [
        step.stator_copper_heat_j,
        step.rotor_copper_heat_j,
        step.bearing_heat_j,
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid(
                "induction winding and bearing heat must be finite nonnegative",
            ));
        }
        heat = heat + Dd::from_f64(value);
    }
    let mut supplied = Dd::ZERO;
    for &value in step.electrical_work_j.iter().chain(&step.shaft_work_j) {
        supplied = supplied + Dd::from_f64(finite(value)?);
    }
    for residual in [
        Dd::from_f64(step.native.dissipated) - heat,
        Dd::from_f64(step.native.supplied) - supplied,
        Dd::from_f64(step.native.delta_h) + heat - supplied,
    ] {
        if scalar(residual)?.abs() > tolerance_j {
            return Err(invalid(
                "induction field/kinetic, source/shaft and physical heat accounts disagree",
            ));
        }
    }
    Ok((scalar(heat)?, supplied))
}

impl InductionThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        }
    }
    fn multipliers(&self) -> [f64; 2] {
        [
            self.stator.resistance_multiplier,
            self.rotor.resistance_multiplier,
        ]
    }
    fn winding_cells(&self) -> [CellId; 2] {
        [self.stator.conductor.cell, self.rotor.conductor.cell]
    }
    fn admit(&self, thermal: &DepositionSimulation<'_>) -> Result<(), Error> {
        let cells = self.winding_cells();
        if cells[0] == cells[1]
            || self
                .multipliers()
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
            || cells.iter().chain(&[self.bearing_cell]).any(|id| {
                !thermal.cells().contains_key(id) || !self.cells.iter().any(|cell| cell.cell == *id)
            })
        {
            return Err(invalid(
                "induction needs distinct active winding cells, a bound bearing receiver and positive fixed resistance conversions",
            ));
        }
        Ok(())
    }
    fn resolve_windings(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        original: Option<&[Arc<ResolvedConductor>; 2]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<[ResolvedConductor; 2], Error> {
        self.admit(thermal)?;
        let resolve = |binding: &InductionWinding, index: usize| {
            // Reuse the original whole-cell geometry, source-selection and
            // pinned-material resolver. The extra receiver is only its existing
            // active-cell admission; no bridge physics or bridge heat is added.
            let network = MotorThermalNetwork {
                cells: self.cells,
                links: self.links,
                winding: binding.conductor.clone(),
                bridge_cell: self.bearing_cell,
            };
            network
                .resolve_winding(thermal, sample, original.map(|w| w[index].as_ref()), cancel)
                .map_err(|error| error.map_source(InductionThermalError::Material))
        };
        Ok([resolve(&self.stator, 0)?, resolve(&self.rotor, 1)?])
    }
    fn resistances(&self, windings: &[ResolvedConductor; 2]) -> Result<[f64; 2], Error> {
        let factors = self.multipliers();
        let values = std::array::from_fn(|i| factors[i] * windings[i].resistance_ohm());
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid(
                "source-bound induction resistance must remain positive and representable",
            ));
        }
        Ok(values)
    }

    /// Bind two original conductor sources without changing any flux, angle,
    /// momentum, enthalpy or clock. Numeric winding-resistance placeholders are
    /// replaced by the actual material values in the declared winding basis.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        state: [f64; 7],
        mut parameters: InductionMotorParameters,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialInductionMotor, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit(thermal)?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let coefficients = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let windings = self.resolve_windings(thermal, &coefficients.source, None, &cancel)?;
        let resistance = self.resistances(&windings)?;
        parameters.stator_resistance_ohm = resistance[0];
        parameters.rotor_resistance_ohm = resistance[1];
        let machine = InductionMotor::new(parameters).map_err(phs)?;
        let operating_point = machine.operating_point(&state).map_err(phs)?;
        finite(machine.system().hamiltonian(&state))?;
        check(&cancel)?;
        Ok(MaterialInductionMotor {
            time_s: thermal.time_s(),
            state,
            parameters,
            operating_point,
            windings: windings.map(Arc::new),
            winding_cells: self.winding_cells(),
            resistance_multipliers: self.multipliers(),
            bearing_cell: self.bearing_cell,
        })
    }

    /// Jointly solve native induction and finite material heating for held
    /// inputs [v_alpha,v_beta,tau_rotor,tau_stator]. Every endpoint trial uses
    /// both original R(T,phase) sources and the same accepted machine state.
    /// Stator and rotor remain distinct finite bodies; bearing heat may join
    /// either winding or another explicit receiver. Thermal links and authored
    /// boundaries are evaluated by the original shared enthalpy owner.
    ///
    /// One native workspace is prepared per coupled interval and reused across
    /// every thermal/Newton probe by rebinding only its four winding loss rows.
    /// Source selections, path geometry and resistance conversions stay pinned.
    /// Source/domain, event, energy, budget and cancellation refusals publish
    /// neither machine nor caloric history. L, inertias, pole pairs and bearing
    /// drag remain fixed; source-temperature dependence is not inferred for them.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialInductionMotor,
        inputs: &[f64; 4],
        controls: InductionEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<InductionEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit(thermal)?;
        admit_driver_unknowns(thermal, controls.thermal, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if dt <= 0.0
            || motor.time_s.to_bits() != thermal.time_s().to_bits()
            || motor.winding_cells != self.winding_cells()
            || motor.bearing_cell != self.bearing_cell
            || motor.resistance_multipliers.map(f64::to_bits)
                != self.multipliers().map(f64::to_bits)
        {
            return Err(invalid(
                "induction and heat require one positive interval and unchanged winding/receiver correspondence",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let original = self.resolve_windings(thermal, &accepted.source, None, &cancel)?;
        let resistance = self.resistances(&original)?;
        if original
            .iter()
            .zip(&motor.windings)
            .any(|(a, b)| !same_conductor(a, b))
            || resistance
                != [
                    motor.parameters.stator_resistance_ohm,
                    motor.parameters.rotor_resistance_ohm,
                ]
        {
            return Err(invalid(
                "both induction conductor sources must match the accepted thermal state",
            ));
        }
        let machine = RefCell::new(
            InductionMotor::new(motor.parameters)
                .map_err(phs)?
                .prepare_controlled(controls.machine, controls.backtracking_halvings, || {
                    (*cancel.borrow_mut())()
                })
                .map_err(native)?,
        );
        let (thermal, motor, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls.thermal,
            initial_energy_increments_j,
            |trial, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let windings = self.resolve_windings(
                    trial,
                    &coefficients.source,
                    Some(&motor.windings),
                    &cancel,
                )?;
                let resistance = self.resistances(&windings)?;
                let windings = windings.map(Arc::new);
                let mut native_machine = machine.borrow_mut();
                native_machine
                    .set_winding_resistances(resistance)
                    .map_err(phs)?;
                let (mut state, mut ports) = ([0.0; 7], [0.0; 4]);
                let step = native_machine
                    .step_controlled(&motor.state, inputs, dt, &mut state, &mut ports, || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(native)?;
                let (heat_j, supplied_j) =
                    physical_account(&step, controls.thermal.total_energy_tolerance_j)?;
                let operating_point = native_machine.operating_point(&state).map_err(phs)?;
                let parameters = native_machine.parameters();
                let next = MaterialInductionMotor {
                    time_s: target_s,
                    state,
                    parameters,
                    operating_point,
                    windings: windings.clone(),
                    winding_cells: motor.winding_cells,
                    resistance_multipliers: motor.resistance_multipliers,
                    bearing_cell: motor.bearing_cell,
                };
                check(&cancel)?;
                Ok(DrivenCoefficients {
                    circuit: InductionDriver {
                        motor: next,
                        step,
                        heat_j,
                        supplied_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: InductionThermalSample {
                        network: coefficients.source,
                        stator: Arc::clone(&windings[0]),
                        rotor: Arc::clone(&windings[1]),
                        resistance_ohm: resistance,
                        winding_cells: motor.winding_cells,
                        bearing_cell: motor.bearing_cell,
                    },
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                receiver_powers(
                    [
                        (
                            candidate.motor.winding_cells[0],
                            candidate.step.stator_copper_heat_j,
                        ),
                        (
                            candidate.motor.winding_cells[1],
                            candidate.step.rotor_copper_heat_j,
                        ),
                        (candidate.motor.bearing_cell, candidate.step.bearing_heat_j),
                    ],
                    duration,
                    checkpoint,
                )
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(InductionEnthalpyCandidate {
            thermal,
            motor,
            report,
        })
    }
}
