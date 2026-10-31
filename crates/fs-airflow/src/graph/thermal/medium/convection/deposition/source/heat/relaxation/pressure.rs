//! Quasi-steady pressure, coolant enthalpy and wall heat at ONE trial endpoint.
//! The original series resistance, convection and advection owners supply all
//! physical equations. There is no lagged flow or second temperature timeline.
use super::*;
use crate::graph::thermal::medium::series::hydraulic::{MaterialFlowResistance,
    pressure::{PressureSection, SeriesPressureBudget, SeriesPressureError, SeriesPressureFlow}};
use crate::graph::thermal::medium::convection::ForcedConvectionSettings;
use crate::graph::thermal::medium::convection::deposition::source::{
    ConvectionPropertyTemperature, SourceCellConvectionLaw,
};
use fs_conduction::transient::deposition::{CellId, DepositionError, advection::MassFlow};
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::transport::RelaxationFlowCoefficients;
use fs_convection::CorrelationId;
use fs_material::MaterialError;
use fs_material::liquid::ResolvedLiquidState;
use crate::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpOperatingControls, PumpOperatingError, PumpOperatingPoint, PumpSupplyState,
};
use fs_qty::Length;
use std::collections::{BTreeMap, BTreeSet};

/// One original finite coolant cell and the resistance of its loop segment.
/// The section's bulk temperature resolves hydraulic rho/mu; it is NOT a film
/// temperature. Section order defines the positive-flow circulation direction.
#[derive(Clone, Copy, Debug)]
pub struct PressureCoolantSection {
    pub cell: CellId,
    pub resistance: MaterialFlowResistance,
}

/// A wall patch cooled by one original hydraulic section. There is no separate
/// speed or flow area: velocity comes from the solved local Q and that section's
/// actual flow area. Heat geometry/correlation and property convention stay explicit.
#[derive(Clone, Copy, Debug)]
pub struct PressureConvectiveLink {
    pub link: u64,
    pub section: usize,
    pub characteristic_length: Length,
    pub property_temperature: ConvectionPropertyTemperature,
    pub correlation: CorrelationId,
    pub settings: ForcedConvectionSettings,
}

/// Held geometry/drive for one interval. Closed circulation needs at least two
/// distinct fixed-inventory coolant cells. The original native graph checks the
/// common caloric MaterialId and each cell's balanced mass flow.
#[derive(Clone, Copy, Debug)]
pub struct PressureCoolantLoop<'a> {
    pub sections: &'a [PressureCoolantSection],
    pub links: &'a [PressureConvectiveLink],
    pub driving_pressure_pa: f64,
    pub pressure_budget: SeriesPressureBudget,
    pub maximum_links: usize,
}

/// Original closed-loop geometry and fixed-speed pump equipment. There is no
/// second prescribed pressure, volume flow or convective speed. The first
/// section is the pump inlet used for head conversion and viscosity admission.
#[derive(Clone, Copy, Debug)]
pub struct PumpCoolantLoop<'a> {
    pub sections: &'a [PressureCoolantSection],
    pub links: &'a [PressureConvectiveLink],
    pub pump: &'a PumpHeadCurve,
    pub pressure_budget: SeriesPressureBudget,
    pub operating_point: PumpOperatingControls,
    pub maximum_links: usize,
}

/// Pump equipment and operating point from the EXACT accepted thermal trial.
/// The reported hydraulic output is not shaft/electrical input or motor loss.
#[derive(Debug)]
pub struct CoolantPumpReceipt<'a> {
    pub curve: &'a PumpHeadCurve,
    pub supply: PumpSupplyState,
    pub network_pressure_pa: f64,
    pub pressure_residual_pa: f64,
    pub network_evaluations: usize,
    pub total_flow_evaluations: usize,
    /// Original pump output minus independently summed local hydraulic loss.
    /// Reported, never assigned to a compensating thermal source.
    pub hydraulic_power_residual_w: f64,
}

#[derive(Debug)]
pub struct PumpRelaxationConvection<'a, F> {
    pub coolant: PressureRelaxationConvection<F, ResolvedLiquidState>,
    pub pump: CoolantPumpReceipt<'a>,
}

#[derive(Debug)]
pub enum PressureCoolantError<E> {
    Invalid(&'static str),
    Source { section: usize, source: E },
    Material { section: usize, source: MaterialError },
    Hydraulic(SeriesPressureError),
    Pump(PumpOperatingError),
    Convection(ConvectionHeatError<E>),
    Thermal(DepositionError),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for PressureCoolantError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pressure-driven coolant: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for PressureCoolantError<E> {}

/// Actual accepted pressure/flow and heat-transfer materials. Bulk hydraulic
/// and declared film sources remain separate views of the same prescribed fluid
/// law, NOT anonymous property replacements. Nothing here advances after acceptance.
#[derive(Debug)]
pub struct PressureRelaxationConvection<F, L = F> {
    pub hydraulics: SeriesPressureFlow<L>,
    pub convection: FiniteRelaxationConvection<F, L>,
    /// All section dp*Q is deposited into its original finite coolant cell.
    /// It is already part of report.external_thermal_heat_j: the held pressure
    /// source supplies this work and has no finite mechanical state in this model.
    pub hydraulic_heat_w: BTreeMap<CellId, f64>,
}

fn checkpoint<E>(cx: &Cx<'_>, gate: &CancelGate) -> Result<(), PressureCoolantError<E>> {
    if cx.checkpoint().is_err() { gate.request(); }
    if gate.is_requested() { Err(PressureCoolantError::Cancelled) } else { Ok(()) }
}
fn finite<E>(x: f64) -> Result<f64, PressureCoolantError<E>> {
    if x.is_finite() { Ok(x) } else { Err(PressureCoolantError::Invalid("unrepresentable coolant pressure/heat")) }
}
fn same_density(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a-b).abs() <= 128.0*f64::EPSILON*a.max(b)
}

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Solve pressure-driven CLOSED circulation within the original nonlinear
    /// enthalpy/circuit/material residual, including all local hydraulic heat.
    /// At each trial: original bulk sources -> pressure/flow -> actual speeds
    /// -> original wall correlations -> the same upwind enthalpy transfers.
    ///
    /// Reuse the initial finite-link admission shown by the prescribed-flow
    /// API. `resolve_coolant` must retain each original caloric fluid's source
    /// pins, pressure/composition and applicability at the supplied temperature.
    /// It is called at bulk T for resistance and at the DECLARED bulk/film T for
    /// convection. A generic trait cannot prove caloric/source correspondence.
    ///
    /// This fixed-volume, fixed-mass loop requires one common density matching
    /// every cell's mass/volume. No density tolerance is a thermal-expansion
    /// model: only roundoff is admitted. Variable-density or compressible flow,
    /// evolving volume, phase change and unknown pressure-work closure refuse
    /// or require another owner. The pressure drop is a held drive, not a pump
    /// curve, shaft/electrical efficiency, fluid-inertia or pressure-EOS solve.
    ///
    /// All dp*Q is thermalized in the corresponding coolant cell. It is external
    /// pressure-source work, distinct from Joule and solid-relaxation loss. The
    /// generic entropy account leaves this explicitly authored power at-cell
    /// unclassified; no reservoir/photon/stream entropy claim is manufactured.
    /// Closed enthalpy circulation remains internal. Do not add hydraulic heat
    /// or advective watts again in any other source callback.
    ///
    /// The final original source payload comes from the accepted trial. No
    /// discarded flow/circuit/thermal evaluation mutates an accepted history.
    /// Refusal is atomic, including unsupported zero-flow correlation states.
    /// Existing first-order material/thermal splitting and upwind spatial error
    /// remain; fixed mechanical exterior loading is not retuned by coolant T.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_pressure_driven_relaxation_with_links<F: ThermalFluidTransport, L: ThermalFluidTransport, M, E>(
        &self, cx: &Cx<'_>,
        owner: &mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal: &mut DepositionSimulation<'_>, external: &[f64],
        electrical: Option<RelaxationElectricalInput<'_>>, controls: FrameThermalSolve,
        energy_guess_j: &[f64], gate: &CancelGate, path: PressureCoolantLoop<'_>,
        resolve_coolant: impl Fn(&PressureCoolantSection, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
    ) -> Result<RelaxationHeatReport<PressureRelaxationConvection<F, L>>,
        RelaxationTransportError<PressureCoolantError<E>>>
    where R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E> {
        checkpoint(cx, gate).map_err(RelaxationTransportError::Source)?;
        if !path.driving_pressure_pa.is_finite() {
            return Err(RelaxationTransportError::Source(PressureCoolantError::Invalid("finite pressure required")));
        }
        self.step_driven_relaxation(cx, owner, thermal, external, electrical, controls,
            energy_guess_j, gate, path.sections, path.links, path.pressure_budget.maximum_sections,
            path.maximum_links, resolve_coolant,
            |sections| Ok((SeriesPressureFlow::solve(sections, path.driving_pressure_pa,
                path.pressure_budget, gate).map_err(PressureCoolantError::Hydraulic)?, ())),
            |hydraulics, convection, hydraulic_heat_w, ()|
                PressureRelaxationConvection { hydraulics, convection, hydraulic_heat_w })
    }

    /// Select pump pressure AND flow inside each original thermal endpoint.
    /// Reuses the existing head-curve matcher, original source resistance and
    /// this same closed-loop heat/advection path. No pump law or second thermal
    /// solver is introduced. Each resolved liquid owner stays with the final
    /// candidate; intermediate pump pressure probes borrow checked snapshots.
    ///
    /// Source viscosity controls the network and the pump's declared validity.
    /// It does not imply an efficiency or viscosity correction to the head
    /// curve. Every probe must stay inside equipment and section support;
    /// unsupported branches, failed brackets or exhausted work limits refuse
    /// the ENTIRE mechanical/electrical/caloric transaction without fallback.
    ///
    /// Local loss dp*Q enters coolant heat exactly once as supplied hydraulic
    /// work. The pump-versus-loss power difference is bounded by the pressure
    /// residual and the interval energy tolerance, and visible in the report.
    /// It is not added as a compensating heater.
    /// No additional motor/shaft loss, speed evolution, compression, gravity
    /// circulation or dynamic fluid inertia is inferred. The original common-
    /// density/fixed-inventory and first-order thermal-split limits still apply.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_pump_driven_relaxation_with_links<'p, F: ThermalFluidTransport, M, E>(
        &self, cx: &Cx<'_>,
        owner: &mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal: &mut DepositionSimulation<'_>, external: &[f64],
        electrical: Option<RelaxationElectricalInput<'_>>, controls: FrameThermalSolve,
        energy_guess_j: &[f64], gate: &CancelGate, path: PumpCoolantLoop<'p>,
        resolve_coolant: impl Fn(&PressureCoolantSection, &DepositionSimulation<'_>, Temperature)
            -> Result<ResolvedLiquidState, E>,
    ) -> Result<RelaxationHeatReport<PumpRelaxationConvection<'p, F>>,
        RelaxationTransportError<PressureCoolantError<E>>>
    where R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E> {
        checkpoint(cx, gate).map_err(RelaxationTransportError::Source)?;
        let sample_period_s = owner.mechanical().sample_period_s();
        self.step_driven_relaxation(cx, owner, thermal, external, electrical, controls,
            energy_guess_j, gate, path.sections, path.links, path.pressure_budget.maximum_sections,
            path.maximum_links, resolve_coolant, |sections| {
                let PumpOperatingPoint { curve, state, supply, network_pressure_pa,
                    pressure_residual_pa, network_evaluations, total_flow_evaluations, .. } =
                    SeriesPressureFlow::solve_with_pump(sections, path.pump, path.pressure_budget,
                        path.operating_point, cx, gate).map_err(PressureCoolantError::Pump)?;
                let hydraulic_power_residual_w = finite(supply.hydraulic_power_w
                    - state.dissipated_power_w())?;
                if finite(sample_period_s * hydraulic_power_residual_w)?.abs()
                    > controls.total_energy_tolerance_j {
                    return Err(PressureCoolantError::Invalid(
                        "pump hydraulic work and original coolant loss exceed the interval energy tolerance"));
                }
                Ok((state, CoolantPumpReceipt { curve, supply, network_pressure_pa,
                    pressure_residual_pa, network_evaluations, total_flow_evaluations,
                    hydraulic_power_residual_w }))
            }, |hydraulics, convection, hydraulic_heat_w, pump| PumpRelaxationConvection {
                coolant: PressureRelaxationConvection { hydraulics, convection, hydraulic_heat_w }, pump,
            })
    }

    // Pressure selection changes; every physical heat/state owner below stays
    // shared. The payload constructor runs within the source evaluation and is
    // pure: it cannot re-query or advance an already accepted state.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_driven_relaxation<F: ThermalFluidTransport, L: ThermalFluidTransport, M, E, D, S>(
        &self, cx: &Cx<'_>, owner: &mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal: &mut DepositionSimulation<'_>, external: &[f64],
        electrical: Option<RelaxationElectricalInput<'_>>, controls: FrameThermalSolve,
        energy_guess_j: &[f64], gate: &CancelGate, sections_geometry: &[PressureCoolantSection],
        links_geometry: &[PressureConvectiveLink], maximum_sections: usize, maximum_links: usize,
        resolve_coolant: impl Fn(&PressureCoolantSection, &DepositionSimulation<'_>, Temperature) -> Result<L, E>,
        resolve_drive: impl Fn(Vec<PressureSection<L>>) -> Result<(SeriesPressureFlow<L>, D), PressureCoolantError<E>>,
        payload: impl Fn(SeriesPressureFlow<L>, FiniteRelaxationConvection<F, L>, BTreeMap<CellId, f64>, D) -> S,
    ) -> Result<RelaxationHeatReport<S>, RelaxationTransportError<PressureCoolantError<E>>>
    where R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E> {
        let admission = || -> Result<BTreeMap<CellId, usize>, PressureCoolantError<E>> {
            checkpoint(cx, gate)?;
            if sections_geometry.len() < 2 || sections_geometry.len() > maximum_sections
                || links_geometry.len() > maximum_links {
                return Err(PressureCoolantError::Invalid("bounded closed loop and links required"));
            }
            let mut cells = BTreeMap::new();
            for (index, section) in sections_geometry.iter().enumerate() {
                checkpoint(cx, gate)?;
                if cells.insert(section.cell, index).is_some() {
                    return Err(PressureCoolantError::Invalid("each loop section needs a distinct original coolant inventory"));
                }
            }
            let mut ids = BTreeSet::new();
            for link in links_geometry {
                checkpoint(cx, gate)?;
                if link.section >= sections_geometry.len() || !ids.insert(link.link) {
                    return Err(PressureCoolantError::Invalid("unique thermal links must name original hydraulic sections"));
                }
            }
            Ok(cells)
        };
        let sections_by_cell = admission().map_err(RelaxationTransportError::Source)?;
        owner.step_heat_network_with_resolved_advection(thermal, external,
            RelaxationHeatInput { network: self.heat, cx }, electrical, controls,
            energy_guess_j, gate, |trial| {
                checkpoint(cx, gate)?;
                let mut sections = Vec::with_capacity(sections_geometry.len());
                let mut thermal_states = Vec::with_capacity(sections_geometry.len());
                for (index, section) in sections_geometry.iter().enumerate() {
                    checkpoint(cx, gate)?;
                    let temperature = trial.temperature_k(section.cell).map_err(PressureCoolantError::Thermal)?;
                    let medium = resolve_coolant(section, trial, Temperature::new(temperature))
                        .map_err(|source| PressureCoolantError::Source { section: index, source })?;
                    let state = medium.thermal_transport_state()
                        .map_err(|source| PressureCoolantError::Material { section: index, source })?;
                    if state.temperature_k().to_bits() != temperature.to_bits() {
                        return Err(PressureCoolantError::Invalid("hydraulic source must use the actual coolant bulk temperature"));
                    }
                    thermal_states.push(state);
                    sections.push(PressureSection { resistance: section.resistance, medium });
                }
                let (pressure, drive) = resolve_drive(sections)?;
                let density = pressure.evaluations()[0].transport.density_kg_m3();
                let mut hydraulic_heat_w = BTreeMap::new();
                for ((section, state), loss) in sections_geometry.iter().zip(&thermal_states)
                    .zip(pressure.evaluations()) {
                    checkpoint(cx, gate)?;
                    let cell = trial.cells().get(&section.cell)
                        .ok_or(PressureCoolantError::Invalid("hydraulic coolant cell is inactive"))?;
                    if state.transport() != loss.transport
                        || !same_density(density, loss.transport.density_kg_m3())
                        || !same_density(density, cell.mass_kg/cell.geometry.volume_m3) {
                        return Err(PressureCoolantError::Invalid("fixed inventory pressure loop requires coherent common density and cell mass/volume"));
                    }
                    hydraulic_heat_w.insert(section.cell, loss.dissipated_power_w);
                }
                let mut convective = Vec::with_capacity(links_geometry.len());
                for link in links_geometry {
                    checkpoint(cx, gate)?;
                    let loss = &pressure.evaluations()[link.section];
                    let speed = finite(loss.volume_flow_m3_s.abs()/loss.law.flow_area().value())?;
                    convective.push(SourceConvectiveLink { link: link.link,
                        fluid_cell: sections_geometry[link.section].cell,
                        characteristic_length: link.characteristic_length,
                        property_temperature: link.property_temperature,
                        law: SourceCellConvectionLaw::Forced { correlation: link.correlation,
                            speed_m_s: speed, settings: link.settings },
                    });
                }
                let ContactThermalCoefficients { mut bulk, contacts } = self.resolve_with_links(
                    cx, trial, &convective, maximum_links, |link, thermal, temperature| {
                        let index = sections_by_cell[&link.fluid_cell];
                        resolve_coolant(&sections_geometry[index], thermal, temperature)
                    }).map_err(PressureCoolantError::Convection)?;
                let ConvectionHeatSample { heat, convection, finite_convection } = bulk.source;
                let finite_links = finite_convection.ok_or(PressureCoolantError::Invalid("missing pressure-resolved wall transfer"))?;
                for (&cell, &power) in &hydraulic_heat_w {
                    checkpoint(cx, gate)?;
                    let previous = bulk.sources_w.get(&cell).copied().unwrap_or(0.0);
                    let next = finite(previous + power)?;
                    if power != 0.0 && next == previous {
                        return Err(PressureCoolantError::Invalid("hydraulic heat disappears beside existing heat"));
                    }
                    if power != 0.0 { bulk.sources_w.insert(cell, next); }
                }
                let mass = pressure.mass_flow_kg_s();
                let mut flows = Vec::with_capacity(sections_geometry.len());
                for (index, section) in sections_geometry.iter().enumerate() {
                    checkpoint(cx, gate)?;
                    let id = u64::try_from(index).map_err(|_| PressureCoolantError::Invalid("loop edge identity overflow"))?;
                    flows.push(MassFlow { id, a: section.cell,
                        b: sections_geometry[(index+1)%sections_geometry.len()].cell, kg_s: mass });
                }
                let flow = Advected::new((), flows, maximum_sections)
                    .map_err(PressureCoolantError::Thermal)?;
                checkpoint(cx, gate)?;
                Ok(RelaxationFlowCoefficients { flow, heat: ContactThermalCoefficients { contacts,
                    bulk: CircuitThermalCoefficients { resistors: bulk.resistors,
                        sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
                        source: (heat, payload(pressure,
                            FiniteRelaxationConvection { reservoirs: convection, finite: finite_links }, hydraulic_heat_w, drive)),
                    } },
                })
            })
    }
}
