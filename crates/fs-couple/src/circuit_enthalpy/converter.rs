//! The native buck/boost trajectory and five physical heat paths on one
//! material/enthalpy endpoint. Electrical integration stays with fs-phs;
//! source evaluation and finite heat stay with their original owners.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, cell_driver_candidate, finite,
    magnetic::same_material,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use crate::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use fs_matdb::PropertyKey;
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductor},
    state_point::{
        MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    },
    thermal::source::ThermalMaterialError,
};
use fs_phs::circuit::converter::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters, CONDUCTION_CHANNELS, ConverterStep,
    switching::{PreparedSwitchCycle, SwitchCycleControls, SwitchInterval},
};
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PhsError, PreparedStepError};
use fs_qty::{Dims, QuantitySpec, Resistance};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

/// Explicit interpretation of one native channel's linear on-state resistance.
/// No option supplies semiconductor switching, recovery, gate or diode physics.
#[derive(Clone, Debug)]
pub enum ConverterResistance<'a> {
    /// Retain the caller-admitted native resistance; its heat still has an owner.
    Fixed,
    /// Original whole-cell uniform conductor law, R=rho(T,f)*length/area.
    Bulk {
        basis: ElectricalPropertyBasis,
        selection: MaterialPropertySelection,
        length_m: f64,
        area_m2: f64,
    },
    /// An explicitly supplied static scalar resistance source in ohms. Its
    /// original card, claim and complete applicability remain authoritative;
    /// a resistivity, efficiency or switching-loss property cannot substitute.
    Resistance {
        key: &'a PropertyKey,
        selection: MaterialPropertySelection,
    },
}

/// One original channel and its sole finite receiving body.
pub struct ConverterThermalChannel<'a> {
    pub cell: CellId,
    pub resistance: ConverterResistance<'a>,
}

/// Original caloric/transport inventory and five distinct heat receivers.
/// Channel order is winding, input high, input low, output high, output low,
/// exactly as in the native `ConductionHeat`. Off channels still require valid
/// material sources; zero duty does not permit an out-of-domain device claim.
pub struct ConverterThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub channels: &'a [ConverterThermalChannel<'a>; CONDUCTION_CHANNELS],
}

#[derive(Clone, Debug, PartialEq)]
enum ChannelSource {
    Fixed,
    Bulk(ResolvedConductor),
    Resistance {
        key: PropertyKey,
        material: ResolvedMaterialStatePoint,
    },
}

/// Original source snapshot retained with the exact receiving channel.
#[derive(Clone, Debug, PartialEq)]
pub struct ConverterChannelMaterial {
    cell: CellId,
    resistance_ohm: f64,
    source: ChannelSource,
}
impl ConverterChannelMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }
    pub const fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }
    pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
        match &self.source {
            ChannelSource::Fixed => None,
            ChannelSource::Bulk(c) => Some(c.material()),
            ChannelSource::Resistance { material, .. } => Some(material),
        }
    }
    pub fn conductor(&self) -> Option<&ResolvedConductor> {
        match &self.source {
            ChannelSource::Bulk(c) => Some(c),
            _ => None,
        }
    }
}

/// Accepted original flux/charges and their actual material resistances. L/C
/// remain unchanged; no separate current or voltage history is integrated.
#[derive(Clone, Debug)]
pub struct MaterialBuckBoost {
    time_s: f64,
    parameters: BuckBoostParameters,
    state: [f64; 3],
    channels: Arc<[ConverterChannelMaterial; CONDUCTION_CHANNELS]>,
}
impl MaterialBuckBoost {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn state(&self) -> &[f64; 3] {
        &self.state
    }
    pub fn channels(&self) -> &[ConverterChannelMaterial; CONDUCTION_CHANNELS] {
        &self.channels
    }
}

/// Either native electrical lane feeds the SAME material/thermal residual.
#[derive(Clone, Copy, Debug)]
pub enum ConverterElectricalDrive<'a> {
    Averaged(BridgeModulation),
    /// The complete native complementary-gate cycle must span this thermal
    /// interval, within floating-point time roundoff. Its current ripple and
    /// channel heat are resolved by the original exact segment propagators.
    Switching {
        intervals: &'a [SwitchInterval],
        controls: SwitchCycleControls,
    },
}

#[derive(Debug)]
pub enum ConverterThermalError {
    Material(MaterialNetworkError),
    Source(ThermalMaterialError),
    Conductor(JouleHeatingError),
    Model(PhsError),
    Electrical(PreparedStepError),
    Invalid(&'static str),
}
impl core::fmt::Display for ConverterThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "converter material/enthalpy: {self:?}")
    }
}
impl std::error::Error for ConverterThermalError {}
type Error = CircuitEnthalpyError<ConverterThermalError>;
fn invalid(what: &'static str) -> Error {
    CircuitEnthalpyError::Source(ConverterThermalError::Invalid(what))
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
        other => CircuitEnthalpyError::Source(ConverterThermalError::Material(other)),
    }
}
fn point(error: ThermalMaterialError) -> Error {
    CircuitEnthalpyError::Source(ConverterThermalError::Source(error))
}
fn source(error: MaterialStatePointError) -> Error {
    point(ThermalMaterialError::Material(error))
}
fn electrical(error: PreparedStepError) -> Error {
    match error {
        PreparedStepError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConverterThermalError::Electrical(other)),
    }
}
fn model(error: PhsError) -> Error {
    CircuitEnthalpyError::Source(ConverterThermalError::Model(error))
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
fn resistances(p: BuckBoostParameters) -> [f64; CONDUCTION_CHANNELS] {
    [
        p.winding_resistance_ohm,
        p.input_on_resistance_ohm[0],
        p.input_on_resistance_ohm[1],
        p.output_on_resistance_ohm[0],
        p.output_on_resistance_ohm[1],
    ]
}
fn with_resistances(
    mut p: BuckBoostParameters,
    channels: &[ConverterChannelMaterial; CONDUCTION_CHANNELS],
) -> BuckBoostParameters {
    p.winding_resistance_ohm = channels[0].resistance_ohm;
    p.input_on_resistance_ohm = [channels[1].resistance_ohm, channels[2].resistance_ohm];
    p.output_on_resistance_ohm = [channels[3].resistance_ohm, channels[4].resistance_ohm];
    p
}
fn same_channel(a: &ConverterChannelMaterial, b: &ConverterChannelMaterial) -> bool {
    a.cell == b.cell
        && a.resistance_ohm.to_bits() == b.resistance_ohm.to_bits()
        && match (&a.source, &b.source) {
            (ChannelSource::Fixed, ChannelSource::Fixed) => true,
            (ChannelSource::Bulk(x), ChannelSource::Bulk(y)) => {
                x.source_basis() == y.source_basis()
                    && x.source_key() == y.source_key()
                    && x.length_m().to_bits() == y.length_m().to_bits()
                    && x.area_m2().to_bits() == y.area_m2().to_bits()
                    && same_material(x.material(), y.material())
            }
            (
                ChannelSource::Resistance {
                    key: x,
                    material: a,
                },
                ChannelSource::Resistance {
                    key: y,
                    material: b,
                },
            ) => x == y && same_material(a, b),
            _ => false,
        }
}

/// Actual endpoint material sources and original finite transport observations.
#[derive(Clone, Debug)]
pub struct ConverterThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[ConverterChannelMaterial; CONDUCTION_CHANNELS]>,
}

pub struct ConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialBuckBoost,
    report: CircuitEnthalpyReport<ConverterThermalSample, ConverterStep>,
}
impl<'t> ConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialBuckBoost {
        &self.converter
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<ConverterThermalSample, ConverterStep> {
        &self.report
    }
    /// Publish the electrical and thermal histories together. Native conduction
    /// heat already entered these cells; it must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialBuckBoost,
        CircuitEnthalpyReport<ConverterThermalSample, ConverterStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}
struct Driver {
    converter: MaterialBuckBoost,
    step: ConverterStep,
}
impl ThermalDriver for Driver {
    type Model = MaterialBuckBoost;
    type Step = ConverterStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.native.dissipated,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.delta_h,
            supplied_j: self.step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

impl ConverterThermalNetwork<'_> {
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
        parameters: BuckBoostParameters,
        retained: Option<&[ConverterChannelMaterial; CONDUCTION_CHANNELS]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<[ConverterChannelMaterial; CONDUCTION_CHANNELS], Error> {
        let fixed = resistances(parameters);
        let mut seen = BTreeSet::new();
        let mut channels = Vec::with_capacity(CONDUCTION_CHANNELS);
        for (index, binding) in self.channels.iter().enumerate() {
            check(cancel)?;
            if !seen.insert(binding.cell) {
                return Err(invalid(
                    "converter channels need distinct original thermal receivers",
                ));
            }
            let body = self
                .cells
                .iter()
                .find(|c| c.cell == binding.cell)
                .ok_or_else(|| {
                    invalid("converter heat receiver lacks its original material binding")
                })?;
            let cell = thermal
                .cells()
                .get(&binding.cell)
                .ok_or_else(|| invalid("converter heat receiver is inactive"))?;
            let phase = sample
                .cells
                .iter()
                .find(|c| c.cell == binding.cell)
                .ok_or_else(|| invalid("converter receiver lacks its actual caloric phase"))?
                .phase;
            let query = body.point.at(phase).map_err(point)?;
            let old = retained.map(|rows| &rows[index]);
            let (resistance_ohm, resolved) = match &binding.resistance {
                ConverterResistance::Fixed => (fixed[index], ChannelSource::Fixed),
                ConverterResistance::Bulk {
                    basis,
                    selection,
                    length_m,
                    area_m2,
                } => {
                    let ohmic = MaterialOhmicBody {
                        material: body.point,
                        selection: selection.clone(),
                        length_m: *length_m,
                        area_m2: *area_m2,
                    };
                    let conductor = if let Some(old) = old {
                        let ChannelSource::Bulk(original) = &old.source else {
                            return Err(invalid(
                                "converter channel cannot change its original source law",
                            ));
                        };
                        let snapshot = original
                            .material()
                            .at_query_point(body.point.material, &query)
                            .map_err(source)?;
                        ohmic.conductor_from_material(
                            phase,
                            cell.geometry.volume_m3,
                            snapshot,
                            *basis,
                        )
                    } else {
                        ohmic.resolve_conductor_with_basis(phase, cell.geometry.volume_m3, *basis)
                    }
                    .map_err(|e| {
                        CircuitEnthalpyError::Source(ConverterThermalError::Conductor(e))
                    })?;
                    (conductor.resistance_ohm(), ChannelSource::Bulk(conductor))
                }
                ConverterResistance::Resistance { key, selection } => {
                    if key.quantity() != QuantitySpec::dimensional(Resistance::DIMS) {
                        return Err(invalid(
                            "converter on-state source must be a static scalar resistance in ohms",
                        ));
                    }
                    let snapshot = if let Some(old) = old {
                        let ChannelSource::Resistance {
                            key: original_key,
                            material,
                        } = &old.source
                        else {
                            return Err(invalid(
                                "converter channel cannot change its original source law",
                            ));
                        };
                        if original_key != *key {
                            return Err(invalid(
                                "converter channel cannot change its original resistance key",
                            ));
                        }
                        material
                            .at_query_point(body.point.material, &query)
                            .map_err(source)?
                    } else {
                        let requirement = ScalarPropertyRequirement::try_with_key(
                            key,
                            ScalarAdmissibility::NonNegative,
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
                    let property = snapshot.property_by_key(key).ok_or_else(|| {
                        invalid("converter channel lost its original resistance property")
                    })?;
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
                            "converter resistance cannot freeze an unowned field, current, gate or frequency condition",
                        ));
                    }
                    (
                        property.value_si(),
                        ChannelSource::Resistance {
                            key: (*key).clone(),
                            material: snapshot,
                        },
                    )
                }
            };
            if !resistance_ohm.is_finite() || resistance_ohm < 0.0 {
                return Err(invalid(
                    "converter channel resistance is not finite and nonnegative",
                ));
            }
            channels.push(ConverterChannelMaterial {
                cell: binding.cell,
                resistance_ohm,
                source: resolved,
            });
        }
        channels
            .try_into()
            .map_err(|_| invalid("all five converter heat channels are required"))
    }

    /// Attach the existing native L/C and flux/charge state to actual finite
    /// material owners. Source bindings replace their corresponding numeric
    /// resistances; Fixed bindings retain them explicitly. Resistance changes
    /// no stored energy. No initial current, bus charge or temperature is reset.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &AveragedBuckBoost,
        state: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialBuckBoost, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|v| !v.is_finite()) || !native.system().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "converter must retain a finite original flux/charge state",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let sample = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let channels = Arc::new(self.resolve_channels(
            thermal,
            &sample.source,
            native.parameters(),
            None,
            &cancel,
        )?);
        let parameters = with_resistances(native.parameters(), &channels);
        AveragedBuckBoost::new(parameters, native.modulation()).map_err(model)?;
        check(&cancel)?;
        Ok(MaterialBuckBoost {
            time_s: thermal.time_s(),
            parameters,
            state,
            channels,
        })
    }

    /// Advance all three native electrical coordinates and five finite heat
    /// receivers inside the original enthalpy/Newton owner. Each trial starts
    /// from the SAME accepted flux/charges and evaluates original pinned R(T,f)
    /// at its proposed thermal endpoint. Averaged or fully switched conduction
    /// is integrated by the native exact electrical trajectory, never by an
    /// endpoint I²R estimate or an efficiency multiplier.
    ///
    /// Native signed terminal work, storage change and five channel heats close
    /// independently before entering the shared local/global thermal gates.
    /// Each exact electrical segment counts against maximum_circuit_trials.
    /// Source, budget, event and cancellation refusals publish neither history.
    ///
    /// L/C, geometry and nonthermal source coordinates stay fixed. No reversible
    /// storage-parameter heat, semiconductor gate/recovery/switching loss, dead
    /// time or phase mixing is inferred. The averaged lane omits PWM RMS loss;
    /// the switched lane resolves its actual conduction ripple. Source and
    /// first-order thermal-endpoint accuracy require their own valid regimes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialBuckBoost,
        drive: ConverterElectricalDrive<'_>,
        terminal_current_a: &[f64; 2],
        electrical_budget: LinearSplitBudget,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConverterEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits() || dt <= 0.0 {
            return Err(invalid(
                "converter and thermal owners need the same positive physical interval",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let initial = self.resolve_channels(
            thermal,
            &accepted.source,
            converter.parameters,
            None,
            &cancel,
        )?;
        if initial
            .iter()
            .zip(converter.channels.iter())
            .any(|(a, b)| !same_channel(a, b))
        {
            return Err(invalid(
                "converter source identities and resistances must match the accepted thermal state",
            ));
        }
        let segments = match drive {
            ConverterElectricalDrive::Averaged(_) => 1,
            ConverterElectricalDrive::Switching { intervals, .. } => intervals.len(),
        };
        if segments == 0 {
            return Err(invalid("converter requires a nonempty electrical interval"));
        }
        let (thermal, converter, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel)?;
                if segments > remaining {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let channels = Arc::new(self.resolve_channels(
                    trial,
                    &coefficients.source,
                    converter.parameters,
                    Some(&converter.channels),
                    &cancel,
                )?);
                let parameters = with_resistances(converter.parameters, &channels);
                let (mut next, mut mean_voltage) = ([0.0; 3], [0.0; 2]);
                let step = match drive {
                    ConverterElectricalDrive::Averaged(modulation) => {
                        let mut native = AveragedBuckBoost::new(parameters, modulation)
                            .map_err(model)?
                            .prepare_controlled(dt, electrical_budget, || (*cancel.borrow_mut())())
                            .map_err(electrical)?;
                        native
                            .step_controlled(
                                &converter.state,
                                terminal_current_a,
                                &mut next,
                                &mut mean_voltage,
                                || (*cancel.borrow_mut())(),
                            )
                            .map_err(electrical)?
                    }
                    ConverterElectricalDrive::Switching {
                        intervals,
                        controls: switching,
                    } => {
                        let mut native = PreparedSwitchCycle::new_controlled(
                            parameters,
                            intervals,
                            electrical_budget,
                            switching,
                            || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?;
                        let period = native.period_s();
                        // An absolute endpoint clock can lose low bits of a
                        // short physical period. Admit that representational
                        // error, but never let tolerance invent a clock tick
                        // for a cycle that cannot advance the accepted time.
                        let represented_end = finite(thermal.time_s() + period)?;
                        let clock_roundoff =
                            0.5 * (clock_spacing(thermal.time_s()) + clock_spacing(target_s));
                        let duration_roundoff = 64.0 * f64::EPSILON * period.max(dt);
                        if represented_end <= thermal.time_s()
                            || (period - dt).abs() > clock_roundoff + duration_roundoff
                        {
                            return Err(invalid(
                                "complete switch cycle must span the same thermal interval",
                            ));
                        }
                        native
                            .step_controlled(
                                &converter.state,
                                terminal_current_a,
                                &mut next,
                                &mut mean_voltage,
                                || (*cancel.borrow_mut())(),
                            )
                            .map_err(electrical)?
                    }
                };
                let mut heat = 0.0;
                for &value in &step.conduction.joules {
                    if finite(value)? < 0.0 {
                        return Err(invalid("negative native converter conduction heat"));
                    }
                    heat = finite(heat + value)?;
                }
                let work = finite(step.terminal_work_j[0] + step.terminal_work_j[1])?;
                finite(step.native.delta_h)?;
                finite(step.native.supplied)?;
                if finite(step.native.dissipated)? < 0.0
                    || (heat - step.native.dissipated).abs() > controls.total_energy_tolerance_j
                    || (work - step.native.supplied).abs() > controls.total_energy_tolerance_j
                {
                    return Err(invalid(
                        "native converter work and five-channel heat do not close",
                    ));
                }
                check(&cancel)?;
                Ok(DrivenCoefficients {
                    circuit: Driver {
                        converter: MaterialBuckBoost {
                            time_s: target_s,
                            parameters,
                            state: next,
                            channels: Arc::clone(&channels),
                        },
                        step,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: ConverterThermalSample {
                        network: coefficients.source,
                        channels,
                    },
                    circuit_trials: segments,
                })
            },
            |candidate, duration, checkpoint| {
                let mut powers = BTreeMap::new();
                for (source, &heat) in candidate
                    .converter
                    .channels
                    .iter()
                    .zip(&candidate.step.conduction.joules)
                {
                    checkpoint()?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "converter channel heat is unrepresentable as a thermal rate",
                        ));
                    }
                    powers.insert(source.cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(ConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
