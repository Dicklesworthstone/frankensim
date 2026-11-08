//! Finite battery work and two physical loss channels on native enthalpy cells.
//!
//! The original battery owns chemical charge, polarization, bus storage and
//! discrete-gradient electrical work. The original thermal network owns heat
//! transport and calorics. Only the two resistances are refreshed at each
//! proposed thermal endpoint; OCV/SOC and both capacitances remain unchanged.
//! No midpoint OCV approximation, second electrical integrator, or additional
//! heat source is inferred from the energy-balance residual.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate, finite,
    magnetic::same_material,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use fs_matdb::PropertyKey;
use fs_material::{
    state_point::{
        MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    },
    thermal::source::ThermalMaterialError,
};
use fs_math::dd::Dd;
use fs_phs::{
    PhsError, PreparedStepError, SolveBudget,
    circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource, BatteryStep},
};
use fs_qty::{Dims, QuantitySpec, Resistance};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// A declared equivalent-circuit resistance, never inferred from conductivity.
#[derive(Clone, Debug)]
pub enum BatteryResistance<'a> {
    /// Retain the positive resistance explicitly supplied to the native battery.
    Fixed,
    /// Positive scalar resistance in ohms from the receiver's original card.
    /// Its selected claim and complete T/phase/fixed-condition domain survive
    /// every refresh. Current-, SOC- or frequency-dependent laws are not frozen
    /// into this temperature-only equivalent circuit.
    Source {
        key: &'a PropertyKey,
        selection: MaterialPropertySelection,
    },
}

/// One original loss channel and its sole finite receiving body.
pub struct BatteryThermalChannel<'a> {
    pub cell: CellId,
    pub resistance: BatteryResistance<'a>,
}

/// Existing caloric/conduction owners and the series/polarization heat routes.
///
/// Channel order is `[series, polarization]`. Both may explicitly name the same
/// cell when the battery is thermally homogenized; neither heat path is dropped.
pub struct BatteryThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub channels: &'a [BatteryThermalChannel<'a>; 2],
}

/// Actual resistance and the original selected source at one accepted state.
#[derive(Clone, Debug, PartialEq)]
pub struct BatteryChannelMaterial {
    cell: CellId,
    resistance_ohm: f64,
    source: Option<(PropertyKey, ResolvedMaterialStatePoint)>,
}
impl BatteryChannelMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }
    pub const fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }
    /// None denotes an explicitly fixed parameter, without material authority.
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source.as_ref().map(|(_, material)| material)
    }
    pub fn source_key(&self) -> Option<&PropertyKey> {
        self.source.as_ref().map(|(key, _)| key)
    }
}

/// One accepted chemical/polarization/bus charge history and its physical clock.
#[derive(Clone, Debug)]
pub struct MaterialBattery {
    time_s: f64,
    curve: BatteryOcvCurve,
    parameters: BatteryParameters,
    state: [f64; 3],
    channels: Arc<[BatteryChannelMaterial; 2]>,
}
impl MaterialBattery {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn curve(&self) -> &BatteryOcvCurve {
        &self.curve
    }
    pub const fn parameters(&self) -> BatteryParameters {
        self.parameters
    }
    pub const fn state(&self) -> &[f64; 3] {
        &self.state
    }
    pub fn channels(&self) -> &[BatteryChannelMaterial; 2] {
        &self.channels
    }
}

#[derive(Clone, Debug)]
pub struct BatteryThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[BatteryChannelMaterial; 2]>,
}

/// Both returned histories cover the same accepted interval. Heat is included.
pub struct BatteryEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    battery: MaterialBattery,
    report: CircuitEnthalpyReport<BatteryThermalSample, BatteryStep>,
}
impl<'t> BatteryEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn battery(&self) -> &MaterialBattery {
        &self.battery
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<BatteryThermalSample, BatteryStep> {
        &self.report
    }
    /// Publish both histories together; do not deposit their electrical heat again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialBattery,
        CircuitEnthalpyReport<BatteryThermalSample, BatteryStep>,
    ) {
        (self.thermal, self.battery, self.report)
    }
}

#[derive(Debug)]
pub enum BatteryThermalError {
    Material(MaterialNetworkError),
    Source(ThermalMaterialError),
    Model(PhsError),
    Electrical(PreparedStepError),
    Invalid(&'static str),
}
impl core::fmt::Display for BatteryThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "battery material/enthalpy: {self:?}")
    }
}
impl std::error::Error for BatteryThermalError {}
type Error = CircuitEnthalpyError<BatteryThermalError>;
fn invalid(what: &'static str) -> Error {
    CircuitEnthalpyError::Source(BatteryThermalError::Invalid(what))
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
        other => CircuitEnthalpyError::Source(BatteryThermalError::Material(other)),
    }
}
fn point(error: ThermalMaterialError) -> Error {
    CircuitEnthalpyError::Source(BatteryThermalError::Source(error))
}
fn source(error: MaterialStatePointError) -> Error {
    point(ThermalMaterialError::Material(error))
}
fn electrical(error: PreparedStepError) -> Error {
    match error {
        PreparedStepError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(BatteryThermalError::Electrical(other)),
    }
}
fn model(error: PhsError) -> Error {
    CircuitEnthalpyError::Source(BatteryThermalError::Model(error))
}

fn with_resistances(
    mut p: BatteryParameters,
    channels: &[BatteryChannelMaterial; 2],
) -> BatteryParameters {
    p.series_resistance_ohm = channels[0].resistance_ohm;
    p.polarization_resistance_ohm = channels[1].resistance_ohm;
    p
}
fn same_channel(a: &BatteryChannelMaterial, b: &BatteryChannelMaterial) -> bool {
    a.cell == b.cell
        && a.resistance_ohm.to_bits() == b.resistance_ohm.to_bits()
        && match (&a.source, &b.source) {
            (None, None) => true,
            (Some((ka, ma)), Some((kb, mb))) => ka == kb && same_material(ma, mb),
            _ => false,
        }
}

struct Driver {
    battery: MaterialBattery,
    step: BatteryStep,
}
impl ThermalDriver for Driver {
    type Model = MaterialBattery;
    type Step = BatteryStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.battery.heat_j(),
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.battery.storage_change_j(),
            supplied_j: self.step.battery.terminal_work_j,
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        let work = &self.step.battery;
        Dd::from_f64(work.chemical_energy_change_j)
            + Dd::from_f64(work.polarization_energy_change_j)
            + Dd::from_f64(work.bus_energy_change_j)
            - Dd::from_f64(work.terminal_work_j)
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.battery, self.step)
    }
}

impl BatteryThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        }
    }

    fn resolve_channels(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        parameters: BatteryParameters,
        retained: Option<&[BatteryChannelMaterial; 2]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<[BatteryChannelMaterial; 2], Error> {
        let fixed = [
            parameters.series_resistance_ohm,
            parameters.polarization_resistance_ohm,
        ];
        let mut resolved = Vec::with_capacity(2);
        for (i, binding) in self.channels.iter().enumerate() {
            check(cancel)?;
            let body = self
                .cells
                .iter()
                .find(|cell| cell.cell == binding.cell)
                .ok_or_else(|| {
                    invalid("battery heat receiver lacks its original material binding")
                })?;
            let phase = sample
                .cells
                .iter()
                .find(|cell| cell.cell == binding.cell)
                .ok_or_else(|| invalid("battery heat receiver is inactive"))?
                .phase;
            if !thermal.cells().contains_key(&binding.cell) {
                return Err(invalid("battery heat receiver is inactive"));
            }
            let old = retained.map(|channels| &channels[i]);
            let (resistance_ohm, selected) = match &binding.resistance {
                BatteryResistance::Fixed => {
                    if old.is_some_and(|channel| channel.source.is_some()) {
                        return Err(invalid(
                            "battery resistance cannot discard its selected source",
                        ));
                    }
                    (fixed[i], None)
                }
                BatteryResistance::Source { key, selection } => {
                    if key.quantity() != QuantitySpec::dimensional(Resistance::DIMS) {
                        return Err(invalid(
                            "battery impedance requires a scalar resistance in ohms",
                        ));
                    }
                    let query = body.point.at(phase).map_err(point)?;
                    let snapshot = if let Some(old) = old {
                        let Some((original_key, original)) = &old.source else {
                            return Err(invalid(
                                "battery resistance cannot change its original source law",
                            ));
                        };
                        if original_key != *key {
                            return Err(invalid(
                                "battery resistance cannot change its original property key",
                            ));
                        }
                        original
                            .at_query_point(body.point.material, &query)
                            .map_err(source)?
                    } else {
                        let requirement = ScalarPropertyRequirement::try_with_key(
                            key,
                            ScalarAdmissibility::StrictlyPositive,
                        )
                        .map_err(source)?;
                        resolve_material_state_point(
                            body.point.material,
                            &query,
                            &[requirement],
                            selection.clone(),
                        )
                        .map_err(source)?
                    };
                    let property = snapshot
                        .property_by_key(key)
                        .ok_or_else(|| invalid("battery resistance lost its original property"))?;
                    body.point.admit_static(property).map_err(point)?;
                    let static_source = match body.point.liquid_fraction_axis {
                        Some(axis) => property.has_only_owned_and_fixed_external_conditions(
                            axis,
                            QuantitySpec::dimensional(Dims::NONE),
                        ),
                        None => property.is_static_at_fixed_external_conditions(),
                    };
                    if !static_source {
                        return Err(invalid(
                            "battery resistance cannot freeze unowned electrical conditions",
                        ));
                    }
                    (property.value_si(), Some(((*key).clone(), snapshot)))
                }
            };
            if !resistance_ohm.is_finite() || resistance_ohm <= 0.0 {
                return Err(invalid("battery resistances must be finite and positive"));
            }
            resolved.push(BatteryChannelMaterial {
                cell: binding.cell,
                resistance_ohm,
                source: selected,
            });
        }
        resolved
            .try_into()
            .map_err(|_| invalid("both battery loss channels are required"))
    }

    /// Bind original charges and OCV/C without advancing or resetting history.
    /// Selected resistances replace the native fixed numbers at the actual
    /// accepted temperature. Fixed parameters gain no material-source authority.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatterySource,
        state: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialBattery, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|value| !value.is_finite())
            || !native.system().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "battery must retain finite charges inside its original SOC support",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let coefficients = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let channels = Arc::new(self.resolve_channels(
            thermal,
            &coefficients.source,
            native.parameters(),
            None,
            &cancel,
        )?);
        let parameters = with_resistances(native.parameters(), &channels);
        BatterySource::new(native.curve().clone(), parameters).map_err(model)?;
        check(&cancel)?;
        Ok(MaterialBattery {
            time_s: thermal.time_s(),
            curve: native.curve().clone(),
            parameters,
            state,
            channels,
        })
    }

    /// Solve finite charge and material self-heating at one shared endpoint.
    ///
    /// Positive terminal current enters the bus (charging); negative current
    /// draws power. The canonical electrical step integrates OCV across every
    /// crossed SOC knot and owns both nonnegative heat channels. Every thermal
    /// trial restarts from the accepted charge state, refreshes R(T,phase), and
    /// deposits only those native joules. Source exhaustion, unsupported data,
    /// events, work budgets or cancellation publish neither history.
    ///
    /// This combines native discrete-gradient electrical work with backward-
    /// Euler thermal rates and endpoint material coefficients: general coupled
    /// accuracy is first order. OCV(T), reversible entropic heat, capacitance
    /// drift, electrochemical aging and protection controls are not supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        battery: &MaterialBattery,
        terminal_current_a: f64,
        electrical_solve: SolveBudget,
        backtracking_halvings: u8,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if battery.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || !terminal_current_a.is_finite()
        {
            return Err(invalid(
                "battery and thermal histories need the same positive interval and finite current",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let initial =
            self.resolve_channels(thermal, &accepted.source, battery.parameters, None, &cancel)?;
        if initial
            .iter()
            .zip(battery.channels.iter())
            .any(|(a, b)| !same_channel(a, b))
        {
            return Err(invalid(
                "battery sources and resistances must match the accepted thermal state",
            ));
        }
        let (thermal, battery, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let channels = Arc::new(self.resolve_channels(
                    trial,
                    &coefficients.source,
                    battery.parameters,
                    Some(&battery.channels),
                    &cancel,
                )?);
                let parameters = with_resistances(battery.parameters, &channels);
                let native =
                    BatterySource::new(battery.curve.clone(), parameters).map_err(model)?;
                let mut solver = native
                    .prepare_controlled(electrical_solve, backtracking_halvings, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(electrical)?;
                let mut next = [0.0; 3];
                let step = solver
                    .step_controlled(
                        &battery.state,
                        terminal_current_a,
                        dt,
                        &mut next,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(electrical)?;
                let account = &step.battery;
                let tolerance = controls.total_energy_tolerance_j;
                for residual in [
                    account.balance_residual_j(),
                    account.storage_change_j() - step.native.delta_h,
                    account.heat_j() - step.native.dissipated,
                    account.terminal_work_j - step.native.supplied,
                ] {
                    if !residual.is_finite() || residual.abs() > tolerance {
                        return Err(invalid(
                            "native battery storage, terminal work and both heat paths must close",
                        ));
                    }
                }
                Ok(DrivenCoefficients {
                    circuit: Driver {
                        battery: MaterialBattery {
                            time_s: target_s,
                            curve: battery.curve.clone(),
                            parameters,
                            state: next,
                            channels: channels.clone(),
                        },
                        step,
                    },
                    circuit_trials: 1,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: BatteryThermalSample {
                        network: coefficients.source,
                        channels,
                    },
                })
            },
            |driver, duration, checkpoint| {
                let mut powers = BTreeMap::new();
                let account = &driver.step.battery;
                for (channel, heat) in driver
                    .battery
                    .channels
                    .iter()
                    .zip([account.series_heat_j, account.polarization_heat_j])
                {
                    checkpoint()?;
                    if !heat.is_finite() || heat < 0.0 {
                        return Err(invalid("battery heat must be finite and nonnegative"));
                    }
                    let power = finite(heat / duration)?;
                    if heat != 0.0 && power == 0.0 {
                        return Err(invalid("battery heat disappeared during power conversion"));
                    }
                    let total = powers.entry(channel.cell).or_insert(0.0);
                    *total = finite((Dd::from_f64(*total) + Dd::from_f64(power)).to_f64())?;
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(BatteryEnthalpyCandidate {
            thermal,
            battery,
            report,
        })
    }
}
