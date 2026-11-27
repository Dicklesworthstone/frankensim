//! Thermal solids whose material card names a registered caloric law,
//! lowered to the cells of one `fs_thermal` `ThermalTimeline`.
//!
//! A part is a cell when it is a `Solid` with only thermal physics and a
//! model on its first material card names a law [`CaloricLaws`] knows: the
//! card, not the kernel's lumped default, picks the path. A cell outside the
//! product law's `T` band is refused at lowering, or mid-run by the law
//! (`MaterialDomain`); nothing clamps.
//! A cell is a box, born at the start: its mass is density times volume, its
//! centreline the box's x extent at its declared pose, its specific energy
//! the law's at its declared temperature (or the ambient). A
//! `thermal-boundary` interface from the cell's face to a face of a
//! stationary frame (a `DrivenBody` part) gives its surface exchange at ambient (`convection_coefficient`, `reservoir_area`
//! and `reservoir_conductance` on its card); without one it exchanges
//! nothing. A `thermal-contact` between two cells is a link: its card's `R''`
//! in series with each cell's bulk conduction to the contact face.
//! A constant source on its [`HEAT_SOURCE_PORT`] is its applied
//! heat; a step source switches it at `t_switch` by one boundary event.
//! The timeline steps at most the time plan's step, checks energy to
//! the cells' shared absolute tolerance, and sends boundary heat through its
//! `environment` port to a sink of that name. A cell a `TemperatureSource`
//! names publishes its temperature at the start of each step as a signal.
//!
//! Under an implicit scheme it steps on fs-thermal's implicit conduction
//! graph, on each law's implicit form (the crate's CONTRACT.md has the rules).

use std::collections::BTreeMap;
use std::rc::Rc;

use fs_matdb::{ConstitutiveModelCard, MaterialCard};
use fs_mbd::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::MaterialRef;
use fs_scenario::assembly::draft::DraftInterfaceKind;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, DraftPartDecl, DraftPartKind, PhysicsDomain,
    PrimitiveShape, SourceTarget, StepPolicy, TimeHistory, TimeScheme,
};
use fs_scenario::entity::EntityId;
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::enthalpy::{ThermalLink, series_contact_conductance_w_k};

use super::{
    Scoped, THERMAL_BOUNDARY_FAMILY, family_entry, of_family, owning_part, query_env, slot,
};
use crate::lowering::ConstitutiveModelRef;
use crate::lowering::{
    CONDUCTANCE_DIMS, CardStore, GroupLoweringContext, GroupSelection, LoweredGroupOutput,
    LoweringError, LoweringGraph, LoweringRegistry, MASS_DIMS, MatDbQueryEnv, Observed, POWER_DIMS,
    ProvenanceRecorder, SourceLoweringReceipt, TEMPERATURE_DIMS, VOLUME_DIMS,
    record_time_history_provenance,
};
use crate::lowering::{CONTACT_RESISTANCE_DIMS, LENGTH_DIMS, query_thermal_conductance};
use crate::{NodeId, RouteTemperature, System, heat_sink};

/// Group name of the thermal timeline.
pub const GROUP: &str = "owners.thermal.timeline";
/// Law id of the product caloric law: the card's constant specific heat
/// from the model's `datum_temperature`, at which the specific energy is zero.
pub const SENSIBLE_HEAT_LAW: &str = "sensible-heat";
/// The port of a cell its applied heat source targets.
pub const HEAT_SOURCE_PORT: &str = "heat_source";
const TIMELINE: &str = "thermal-timeline";
/// Card families of a contact between two cells: a resistance, a conductance.
const CONTACT_FAMILIES: [&str; 2] = ["thermal-contact-resistance", "thermal-conductance"];
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const HEAT_TRANSFER: Dims = Dims([0, 1, -3, -1, 0, 0]);
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);

/// What a caloric-law factory is asked for: the law of one cell's material,
/// and the specific energy it holds at the cell's initial temperature.
pub struct CaloricRequest<'a> {
    /// Receipt scope of the cell (`part:<name>`).
    pub item: &'a str,
    /// Card store the material resolves in.
    pub cards: &'a CardStore,
    /// The cell's material.
    pub material: &'a MaterialRef,
    /// The constitutive model naming the law.
    pub model: &'a ConstitutiveModelCard,
    /// Ambient state and claim policy for card queries.
    pub env: MatDbQueryEnv,
    /// The cell's initial temperature [K].
    pub temperature_k: f64,
    /// Whether the timeline steps implicitly, so needs the law's implicit form.
    pub implicit: bool,
}

/// A cell's caloric law and its initial specific energy [J/kg].
pub struct CaloricCell {
    /// The law, held by the timeline itself.
    pub law: Rc<dyn ThermalLaw>,
    /// The same law object's implicit form, if it has one (fs-thermal's graph
    /// matches it to `law` by address); an implicit timeline needs it.
    pub implicit: Option<Rc<dyn ImplicitCaloricLaw>>,
    /// Specific energy at the initial temperature, on the law's own datum.
    pub specific_energy_j_kg: f64,
}

/// A caloric law with its implicit conduction form, as one object: the
/// implicit graph matches a law to the simulation's by address.
pub trait ImplicitCaloricLaw: ImplicitConductionLaw + ConductivityLaw + ThermalLaw {}

impl<T: ImplicitConductionLaw + ConductivityLaw + ThermalLaw> ImplicitCaloricLaw for T {}

type Factory =
    Box<dyn Fn(&CaloricRequest<'_>, &mut ProvenanceRecorder) -> Result<CaloricCell, LoweringError>>;

/// Caloric-law factories by constitutive law id. A factory records every
/// value it reads in the recorder it is given.
#[derive(Default)]
pub struct CaloricLaws {
    factories: BTreeMap<String, Factory>,
}

impl CaloricLaws {
    /// No laws: no part lowers to a cell.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The product laws: [`SENSIBLE_HEAT_LAW`].
    #[must_use]
    pub fn product() -> Self {
        let mut laws = Self::new();
        laws.register(SENSIBLE_HEAT_LAW, sensible_heat);
        laws
    }

    /// Register (or replace) the factory of constitutive law `law`.
    pub fn register<F>(&mut self, law: impl Into<String>, factory: F)
    where
        F: Fn(&CaloricRequest<'_>, &mut ProvenanceRecorder) -> Result<CaloricCell, LoweringError>
            + 'static,
    {
        self.factories.insert(law.into(), Box::new(factory));
    }

    /// The first model of `models` whose law has a factory, and the factory.
    fn model<'m>(
        &self,
        models: &'m [ConstitutiveModelCard],
    ) -> Option<(&'m ConstitutiveModelCard, &Factory)> {
        models
            .iter()
            .find_map(|model| Some((model, self.factories.get(&model.law.0)?)))
    }
}

/// A constant specific heat `c` [J/(kg K)] from a datum `t0` [K]:
/// `e = c (T - t0)`, within the model's `band` [K] when it has one.
struct SensibleHeat {
    c_j_kg_k: f64,
    datum_k: f64,
    band: Option<(f64, f64)>,
    /// Conductivity [W/(m K)], held for the run; only an implicit cell has it.
    k_w_m_k: Option<f64>,
}

impl ThermalLaw for SensibleHeat {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let temperature_k = self.datum_k + specific_h_j_kg / self.c_j_kg_k;
        let within = |(lo, hi): (f64, f64)| (lo..=hi).contains(&temperature_k);
        let banded = self.band.is_none_or(within);
        if !(temperature_k.is_finite() && temperature_k > 0.0 && banded) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(temperature_k)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.c_j_kg_k
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial: f64,
        end: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let ratio = self.temperature_k(end)? / self.temperature_k(initial)?;
        Ok(Some(self.c_j_kg_k * ratio.ln()))
    }
}

impl ConductivityLaw for SensibleHeat {
    fn conductivity_w_m_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.temperature_k(specific_h_j_kg)?;
        self.k_w_m_k.ok_or(DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for SensibleHeat {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    /// The band's ends (an implicit cell has a band; without one, NaN, which
    /// fs-thermal refuses).
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        let (lo, hi) = self.band.unwrap_or((f64::NAN, f64::NAN));
        [lo, hi].map(|t| self.c_j_kg_k * (t - self.datum_k))
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(specific_h_j_kg)?,
            temperature_slope_kg_k_j: 1.0 / self.c_j_kg_k,
            conductivity_w_m_k: self.conductivity_w_m_k(specific_h_j_kg)?,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The product law: the card's `specific_heat_capacity`, from the model's
/// `datum_temperature`, within its `T` validity band (receipted as
/// `validity.T.lower` and `validity.T.upper`) when it has one.
fn sensible_heat(
    request: &CaloricRequest<'_>,
    r: &mut ProvenanceRecorder,
) -> Result<CaloricCell, LoweringError> {
    let (item, name, model) = (request.item, "specific_heat_capacity", request.model);
    let (cards, material, env) = (request.cards, request.material, request.env);
    let (c, _) = r.query_material_property(item, cards, material, name, name, env)?;
    let datum = model.parameters.get("datum_temperature");
    let datum = datum.filter(|datum| datum.dims == TEMPERATURE_DIMS && datum.value.is_finite());
    let (Some(datum), true) = (datum, c.dims == SPECIFIC_HEAT && c.value > 0.0) else {
        let detail = format!("needs a positive {name} and a datum_temperature [K]");
        return Err(LoweringError::Invalid {
            item: item.to_owned(),
            detail,
        });
    };
    let reference = ConstitutiveModelRef {
        card_hash: material.card,
        model_hash: model.content_hash(),
        law: &model.law.0,
        law_version: model.law_version,
        parameter: "datum_temperature",
    };
    r.record_model_parameter("datum_temperature", datum.value, datum.dims, reference);
    let band = model.validity.bound("T");
    if let Some((lo, hi)) = band {
        let mut bound = reference;
        for (parameter, value) in [("validity.T.lower", lo), ("validity.T.upper", hi)] {
            bound.parameter = parameter;
            r.record_model_parameter(parameter, value, TEMPERATURE_DIMS, bound);
        }
        if !(lo..=hi).contains(&request.temperature_k) {
            let t = request.temperature_k;
            let detail = format!("initial {t} K is outside the law's band [{lo}, {hi}] K");
            return Err(LoweringError::Invalid {
                item: item.to_owned(),
                detail,
            });
        }
    }
    let refuse = |detail: &str| LoweringError::Invalid {
        item: item.to_owned(),
        detail: detail.to_owned(),
    };
    let k_w_m_k = if request.implicit {
        if band.is_none() {
            return Err(refuse("an implicit cell needs the law's `T` validity band"));
        }
        let (name, held) = ("thermal_conductivity", "thermal_conductivity.held");
        let (k, _) = r.query_material_property(item, cards, material, name, held, env)?;
        if !(k.dims == CONDUCTIVITY && k.value > 0.0) {
            return Err(refuse("needs a positive thermal_conductivity [W/(m K)]"));
        }
        Some(k.value)
    } else {
        None
    };
    let (c_j_kg_k, datum_k) = (c.value, datum.value);
    let law = Rc::new(SensibleHeat {
        c_j_kg_k,
        datum_k,
        band,
        k_w_m_k,
    });
    let implicit = Rc::clone(&law) as Rc<dyn ImplicitCaloricLaw>;
    Ok(CaloricCell {
        implicit: request.implicit.then_some(implicit),
        law,
        specific_energy_j_kg: c_j_kg_k * (request.temperature_k - datum_k),
    })
}

/// Register the thermal-timeline group, resolving caloric laws in `laws`.
pub(crate) fn register(registry: &mut LoweringRegistry, laws: CaloricLaws) {
    let laws = Rc::new(laws);
    let selecting = Rc::clone(&laws);
    registry.register_group(
        GROUP,
        move |draft: &DraftAssembly, cards: &CardStore| select(&selecting, draft, cards),
        move |ctx: &GroupLoweringContext<'_>, graph: &mut LoweringGraph| lower(&laws, ctx, graph),
    );
}

fn is_cell(laws: &CaloricLaws, cards: &CardStore, part: &DraftPartDecl) -> bool {
    let card = part
        .materials
        .first()
        .and_then(|material| cards.material(&material.card));
    part.kind == DraftPartKind::Solid
        && part.physics.iter().eq([&PhysicsDomain::Thermal])
        && card.is_some_and(|card| laws.model(card.models()).is_some())
}

/// The cells, their boundary interfaces, the contacts between them and their
/// heat sources, as one selection.
#[allow(clippy::unnecessary_wraps)] // Signature required by GroupSelector.
fn select(
    laws: &CaloricLaws,
    draft: &DraftAssembly,
    cards: &CardStore,
) -> Result<Vec<GroupSelection>, LoweringError> {
    let input = draft.input();
    let mut selection = of_family(draft, cards, THERMAL_BOUNDARY_FAMILY);
    selection.parts = (0..input.parts.len())
        .filter(|&i| is_cell(laws, cards, &input.parts[i]))
        .collect();
    let cell = |end| owning_part(input, end).is_some_and(|part| selection.parts.contains(&part));
    for family in CONTACT_FAMILIES {
        for index in of_family(draft, cards, family).interfaces {
            let iface = &input.interfaces[index];
            if cell(iface.feature_a) && cell(iface.feature_b) {
                selection.interfaces.push(index);
            }
        }
    }
    let heats_a_cell = |target: &SourceTarget| match target {
        SourceTarget::Port { entity, port } if port == HEAT_SOURCE_PORT => {
            owning_part(input, *entity).is_some_and(|part| selection.parts.contains(&part))
        }
        _ => false,
    };
    let sources = input
        .sources
        .iter()
        .enumerate()
        .filter(|(_, source)| heats_a_cell(&source.target));
    selection.sources = sources.map(|(index, _)| index).collect();
    let found = !selection.parts.is_empty();
    Ok(found.then_some(selection).into_iter().collect())
}

/// One timeline of the selected cells.
fn lower(
    laws: &CaloricLaws,
    ctx: &GroupLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredGroupOutput, LoweringError> {
    let input = ctx.draft.input();
    let ambient_k = input.environment.base.ambient_temperature.value;
    let mut out = LoweredGroupOutput::default();
    let boundaries = boundaries(ctx, &mut out)?;
    let (contacts, conductive) = contacts(ctx, &mut out)?;
    let sources = sources(ctx, &mut out)?;
    let first = &input.parts[ctx.selection.parts[0]];
    let budget = &first.fidelity.error_budget;
    let (StepPolicy::Fixed { dt: step } | StepPolicy::Adaptive { dt_max: step, .. }) =
        &input.time_plan.step_policy;
    let config = StepConfig {
        max_dt_s: step.value,
        energy_tolerance_j: budget.absolute_tol.value,
    };
    let work = budget
        .max_work_units
        .and_then(|work| usize::try_from(work).ok());
    let Some(work) = work else {
        return Err(Scoped::new("part", &first.name).invalid("a timeline needs max_work_units"));
    };
    let scheme = input.time_plan.scheme;
    let (mut schedule, mut cell_laws, mut implicit) = (Vec::new(), Vec::new(), Vec::new());
    for (cell, &index) in ctx.selection.parts.iter().enumerate() {
        let part = &input.parts[index];
        let mut it = Scoped::new("part", &part.name);
        if part.fidelity.error_budget != *budget {
            return Err(it.invalid("the cells of a timeline share one error budget"));
        }
        it.declare_on("time_plan", "step_policy.dt", *step);
        it.declare("fidelity.error_budget.absolute_tol", budget.absolute_tol);
        let work_units = QtyAny::dimensionless(work as f64);
        it.declare("fidelity.error_budget.max_work_units", work_units);
        let source_w = sources.get(&cell).map_or(0.0, |&(power, _)| power);
        let boundary = match boundaries.get(&cell) {
            Some(&[htc_w_m2_k, area_m2, conductance_w_k]) => {
                let reservoir = Reservoir {
                    area_m2,
                    conductance_w_k,
                    temperature_k: ambient_k,
                };
                Boundary {
                    ambient_k,
                    htc_w_m2_k,
                    // fs-thermal's implicit graph refuses an unbound `Some`
                    // reservoir; a zero one, all `boundaries` admits there,
                    // exchanges nothing.
                    reservoir: (scheme == TimeScheme::Explicit).then_some(reservoir),
                    source_w,
                }
            }
            None => Boundary {
                ambient_k,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w,
            },
        };
        let (born, caloric) = new_cell(laws, ctx, &mut it, part, cell as u64, boundary)?;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(born),
        });
        implicit.extend(caloric.implicit.map(|law| (MaterialId(cell as u64), law)));
        cell_laws.push((MaterialId(cell as u64), caloric.law));
        out.parts.push(it.part(TIMELINE));
    }
    let switches = switched(&schedule, &sources);
    let (cells, links, steps) = (cell_laws.len(), contacts.len(), switches.len());
    schedule.extend(contacts);
    schedule.extend(switches);
    let fail =
        |err: &dyn std::fmt::Debug| Scoped::new("node", TIMELINE).invalid(format!("{err:?}"));
    let limits = DepositionBudget {
        max_cells: cells,
        max_links: cells + links,
        max_events: cells + links + steps,
        max_substeps: work,
    };
    let simulation = DepositionSimulation::new_shared(schedule, &cell_laws, limits);
    let simulation = simulation.and_then(|start| start.advance_to(0.0, config, || Ok(())));
    let simulation = simulation.map_err(|err| fail(&err))?.0;
    // An implicit timeline always has its `environment` port.
    let environment =
        scheme != TimeScheme::Explicit || !boundaries.is_empty() || !sources.is_empty();
    let implicit = (implicit, conductive);
    let timeline = timeline(&simulation, &implicit, environment, config, work, scheme)?;
    let timeline = timeline.publishing(published(ctx));
    let timeline = timeline.map_err(|err| fail(&err))?;
    let start = timeline.start(&simulation).map_err(|err| fail(&err))?;
    let node = graph.add(TIMELINE, Some(first.entity), timeline, start)?;
    for (cell, &index) in ctx.selection.parts.iter().enumerate() {
        let entity = input.parts[index].entity;
        if cell > 0 {
            graph.entity_nodes.entry(entity).or_default().push(node);
        }
        expose_cell(graph, (node, entity, cell))?;
    }
    if environment {
        let port = graph.port(node, "environment")?;
        let sink = heat_sink("environment", ambient_k).map_err(|err| fail(&err))?;
        let temperature = RouteTemperature::Fixed(ambient_k);
        graph.add_heat_route(port, Vec::new(), Some(sink), temperature);
    }
    Ok(out)
}

/// The cells' laws' implicit forms and the contacts' links.
type ImplicitForms = (
    Vec<(MaterialId, Rc<dyn ImplicitCaloricLaw>)>,
    Vec<ConductiveLink>,
);

/// The cells' timeline under the plan's `scheme`. An implicit one steps on
/// fs-thermal's implicit conduction graph of the laws' implicit forms and
/// the contacts' links, with the plan's step, the cells' absolute tolerance
/// and work, and the graph's other controls at their defaults.
fn timeline(
    simulation: &DepositionSimulation<'static>,
    (laws, links): &ImplicitForms,
    environment: bool,
    config: StepConfig,
    work: usize,
    scheme: TimeScheme,
) -> Result<ThermalTimeline<'static>, LoweringError> {
    let fail =
        |err: &dyn std::fmt::Debug| Scoped::new("node", TIMELINE).invalid(format!("{err:?}"));
    let ids = simulation.cells().keys().copied().collect();
    let scheme = match scheme {
        TimeScheme::Explicit => {
            return ThermalTimeline::new(ids, environment, config).map_err(|err| fail(&err));
        }
        TimeScheme::BackwardEuler => ImplicitScheme::BackwardEuler,
        TimeScheme::Sdirk2 => ImplicitScheme::Sdirk2,
    };
    let (mut conducting, mut stepping) = (Vec::new(), Vec::new());
    for (id, law) in laws {
        conducting.push((*id, Rc::clone(law) as Rc<dyn ConductivityLaw>));
        stepping.push((*id, Rc::clone(law) as Rc<dyn ImplicitConductionLaw>));
    }
    let ok = || Ok(());
    let observer = ConductionObserver::new_shared(simulation, &conducting, links, &[], ok);
    let observer = observer.map_err(|err| fail(&err))?;
    let linear = ImplicitConductionLinearControls::default();
    let graph = ImplicitConductionGraph::new_shared(simulation, &observer, &stepping, linear, ok);
    let controls = ImplicitConductionControls {
        maximum_dt_s: config.max_dt_s,
        maximum_steps: work,
        cell_tolerance_j: config.energy_tolerance_j,
        total_tolerance_j: config.energy_tolerance_j,
        ..ImplicitConductionControls::default()
    };
    let graph = graph.map_err(|err| fail(&err))?;
    ThermalTimeline::implicit(ids, graph, controls, scheme).map_err(|err| fail(&err))
}

/// The cells the draft's temperature sources name, in cell order.
fn published(ctx: &GroupLoweringContext<'_>) -> Vec<CellId> {
    let input = ctx.draft.input();
    let sources = input.interfaces.iter();
    let sources = sources.filter(|iface| iface.kind == DraftInterfaceKind::TemperatureSource);
    let cells = sources.filter_map(|iface| slot(input, &ctx.selection.parts, iface.feature_b));
    let mut cells: Vec<CellId> = cells.map(|cell| CellId(cell as u64)).collect();
    cells.sort_unstable();
    cells.dedup();
    cells
}

/// Exposes cell `cell` on its `entity`: its heat port as `heat`, its published
/// temperature, if it is, as `temperature`, and its reader for probes.
fn expose_cell(
    graph: &mut LoweringGraph,
    (node, entity, cell): (NodeId, EntityId, usize),
) -> Result<(), LoweringError> {
    let port = graph.port(node, &format!("cell{cell}"))?;
    graph.expose_port(entity, "heat", port);
    if let Ok(port) = graph.port(node, &format!("temperature{cell}")) {
        graph.expose_port(entity, "temperature", port);
    }
    observe_cell(graph, node, entity, cell);
    Ok(())
}

/// Registers cell `cell`'s reader for probes: its temperature, from its slot in
/// the timeline's state.
fn observe_cell(graph: &mut LoweringGraph, node: NodeId, entity: EntityId, cell: usize) {
    let id = CellId(cell as u64);
    let read = move |system: &System| {
        let timeline = system.state::<TimelineState<'static>>(node)?.simulation();
        timeline.temperature_k(id).ok().map(Observed::Cell)
    };
    graph.observers.insert(entity, Box::new(read));
}

/// Each boundary's cell and exchange: the convection coefficient, and the
/// reservoir's area and conductance.
fn boundaries(
    ctx: &GroupLoweringContext<'_>,
    out: &mut LoweredGroupOutput,
) -> Result<BTreeMap<usize, [f64; 3]>, LoweringError> {
    let input = ctx.draft.input();
    let mut boundaries = BTreeMap::new();
    for &index in &ctx.selection.interfaces {
        let iface = &input.interfaces[index];
        let Some(card) = family_entry(ctx.cards, iface, THERMAL_BOUNDARY_FAMILY) else {
            continue; // A contact between cells: `contacts`.
        };
        let mut it = Scoped::new("interface", &iface.name);
        let ends =
            [iface.feature_a, iface.feature_b].map(|end| slot(input, &ctx.selection.parts, end));
        let frame = |end| {
            owning_part(input, end)
                .is_some_and(|p| input.parts[p].kind == DraftPartKind::DrivenBody)
        };
        let cell = match ends {
            [Some(cell), None] if frame(iface.feature_b) => cell,
            [None, Some(cell)] if frame(iface.feature_a) => cell,
            _ => return Err(it.invalid("a thermal boundary joins a cell to a `DrivenBody` frame")),
        };
        let exchange = [
            it.interface(ctx, card, "convection_coefficient", HEAT_TRANSFER)?,
            it.interface(ctx, card, "reservoir_area", AREA)?,
            it.interface(ctx, card, "reservoir_conductance", CONDUCTANCE_DIMS)?,
        ];
        let implicit = input.time_plan.scheme != TimeScheme::Explicit;
        if implicit && (exchange[1] != 0.0 || exchange[2] != 0.0) {
            let why =
                "an implicit thermal boundary's reservoir needs fs-thermal's reservoir binding";
            return Err(it.invalid(why));
        }
        let ambient = input.environment.base.ambient_temperature;
        it.declare_on("environment", "base.ambient_temperature", ambient);
        if boundaries.insert(cell, exchange).is_some() {
            return Err(it.invalid("a cell has one thermal boundary"));
        }
        out.interfaces.push(it.interface_output(Some(TIMELINE)));
    }
    Ok(boundaries)
}

/// Each contact between two cells, as a link from the start: the card's
/// `R''` at `:area` in series with each cell's bulk conduction from its centre
/// to the contact face, `G = A / (d_a/k_a + d_b/k_b + R'')`, each `k` held.
/// An implicit graph gets the same `d` and `R''` as a `ConductiveLink`.
fn contacts(
    ctx: &GroupLoweringContext<'_>,
    out: &mut LoweredGroupOutput,
) -> Result<(Vec<TimedEvent>, Vec<ConductiveLink>), LoweringError> {
    let input = ctx.draft.input();
    let (mut links, mut conductive) = (Vec::new(), Vec::new());
    for &index in &ctx.selection.interfaces {
        let iface = &input.interfaces[index];
        let card = CONTACT_FAMILIES.map(|family| family_entry(ctx.cards, iface, family));
        let (resistance, Some((entry, state))) = (card[0].is_some(), card[0].or(card[1])) else {
            continue;
        };
        let mut it = Scoped::new("interface", &iface.name);
        let DraftInterfaceKind::ThermalContact { third_body, .. } = &iface.kind else {
            continue;
        };
        if third_body.is_some() {
            return Err(it.invalid("a thermal contact with a third body has no owner form"));
        }
        let area = state.iter().find(|(key, _)| key.name == "area");
        let Some(&(_, area)) = area.filter(|(_, q)| q.dims == AREA && q.value > 0.0) else {
            return Err(it.invalid("a thermal contact needs a positive :area [m^2]"));
        };
        let (env, item, area) = (query_env(ctx), it.item.clone(), area.value);
        let g = query_thermal_conductance(&item, &mut it.r, entry, state, resistance, env)?;
        if !resistance {
            it.r.record_declared_field("area", area, AREA, &item, "state.area");
        }
        let (r, why) = (area / g, "A / thermal_conductance");
        it.r.record_declared_field("contact_resistance", r, CONTACT_RESISTANCE_DIMS, &item, why);
        let mut ends = [(CellId(0), 0.0, 0.0); 2];
        for (side, end) in [iface.feature_a, iface.feature_b].into_iter().enumerate() {
            let name = ["a", "b"][side];
            let cell = slot(input, &ctx.selection.parts, end);
            let cell = cell.ok_or_else(|| it.invalid("a contact joins two cells"))?;
            let part = &input.parts[ctx.selection.parts[cell]];
            let DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) = part.geometry
            else {
                return Err(it.invalid("a cell is a box"));
            };
            let Some(d) = half_depth([lx, ly, lz].map(|side| side.value), area) else {
                let why = format!("`{}`: a contact covers a whole box face", part.name);
                return Err(it.invalid(why));
            };
            let held = format!("{name}.thermal_conductivity.held");
            let (k, material) = ("thermal_conductivity", &part.materials[0]);
            let at = format!("{item}: the conductivity of `{}`", part.name);
            let (k, _) =
                it.r.query_material_property(&at, ctx.cards, material, k, &held, env)?;
            if !(k.dims == CONDUCTIVITY && k.value > 0.0) {
                return Err(it.invalid(format!("`{}` needs a positive {held}", part.name)));
            }
            let (half, rule) = (format!("{name}.half_distance"), "half the normal side");
            it.derive(part, &half, d, LENGTH_DIMS, rule);
            let bulk = format!("{name}.bulk_resistance");
            it.r.record_declared_field(bulk, d / k.value, CONTACT_RESISTANCE_DIMS, &item, "d / k");
            ends[side] = (CellId(cell as u64), d, k.value);
        }
        let [(a, d_a, k_a), (b, d_b, k_b)] = ends;
        let g = series_contact_conductance_w_k(area, [d_a, d_b], [k_a, k_b], r);
        let g = g.map_err(|err| it.invalid(format!("{err:?}")))?;
        let why = "A / (d_a/k_a + d_b/k_b + R'')";
        it.r.record_declared_field("contact_conductance", g, CONDUCTANCE_DIMS, &item, why);
        let link = ThermalLink {
            id: index as u64,
            a,
            b,
            area_m2: area,
            conductance_w_k: g,
        };
        links.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link),
        });
        conductive.push(ConductiveLink {
            link_id: index as u64,
            half_distances_m: [d_a, d_b],
            contact_resistance_m2_k_w: r,
        });
        out.interfaces.push(it.interface_output(Some(TIMELINE)));
    }
    Ok((links, conductive))
}

/// Half the side of a box of sides `sides` normal to its face of area `area`.
/// Faces of equal area have equal normal sides (a face is the volume over its
/// normal side), so it is unique. The match is to 1e-9 relative: a declared
/// area is one decimal, a face the product of two (0.05 x 0.05 != 0.0025).
fn half_depth([x, y, z]: [f64; 3], area: f64) -> Option<f64> {
    let faces = [(y * z, x), (z * x, y), (x * y, z)];
    let matches = |(face, _): &(f64, f64)| (face - area).abs() <= 1e-9 * area;
    faces.into_iter().find(matches).map(|(_, side)| 0.5 * side)
}

/// Each source's cell, its power from the start and, for a step, when it
/// switches and to what.
type Sources = BTreeMap<usize, (f64, Option<[f64; 2]>)>;

/// The draft's heat sources.
fn sources(
    ctx: &GroupLoweringContext<'_>,
    out: &mut LoweredGroupOutput,
) -> Result<Sources, LoweringError> {
    let input = ctx.draft.input();
    let mut sources = BTreeMap::new();
    for &index in &ctx.selection.sources {
        let source = &input.sources[index];
        let mut it = Scoped::new("source", &source.name);
        let SourceTarget::Port { entity, .. } = &source.target else {
            return Err(it.invalid("a cell's heat source targets its port"));
        };
        let cell = slot(input, &ctx.selection.parts, *entity);
        // The draft holds a step's `after` to its `before`'s units.
        let held = match &source.history {
            TimeHistory::Constant(power) => Some((power, None)),
            TimeHistory::Step {
                t_switch,
                before,
                after,
            } => Some((before, Some([t_switch.value, after.value]))),
            _ => None,
        };
        let (Some((power, switch)), Some(cell)) = (held, cell) else {
            return Err(it.invalid("a cell's heat source is constant or a step"));
        };
        if power.dims != POWER_DIMS || sources.insert(cell, (power.value, switch)).is_some() {
            return Err(it.invalid("a cell has at most one heat source, a power"));
        }
        record_time_history_provenance(&it.item, &source.history, &mut it.r);
        out.sources.push(SourceLoweringReceipt {
            source_name: source.name.clone(),
            target_description: format!("{}.{HEAT_SOURCE_PORT}", entity.token()),
            created_nodes: vec![TIMELINE.to_owned()],
            parameter_sources: it.r.finish().1,
        });
    }
    Ok(sources)
}

/// The stepping sources' events, in time order (fs-thermal's): each changes
/// its cell's boundary at birth, the first `cells` of `births`, in its
/// source alone.
fn switched(births: &[TimedEvent], sources: &Sources) -> Vec<TimedEvent> {
    let mut events = Vec::new();
    for (&cell, &(_, switch)) in sources {
        let (Some([time_s, source_w]), Event::Birth(born)) = (switch, &births[cell].event) else {
            continue;
        };
        let mut boundary = born.boundary;
        boundary.source_w = source_w;
        let event = Event::Boundary(born.id, boundary);
        events.push(TimedEvent { time_s, event });
    }
    events.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    events
}

/// A cell born at the start, and its law.
fn new_cell(
    laws: &CaloricLaws,
    ctx: &GroupLoweringContext<'_>,
    it: &mut Scoped,
    part: &DraftPartDecl,
    id: u64,
    boundary: Boundary,
) -> Result<(DepositedCell, CaloricCell), LoweringError> {
    let ambient = ctx.draft.input().environment.base.ambient_temperature;
    let initial = part
        .initial
        .as_ref()
        .ok_or_else(|| it.invalid("a cell declares its initial pose"))?;
    let DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) = part.geometry else {
        return Err(it.invalid("a cell is a box"));
    };
    let temperature = if let Some(temperature) = initial.temperature {
        it.declare("initial.temperature", temperature);
        temperature
    } else {
        it.declare_on("environment", "base.ambient_temperature", ambient);
        ambient
    };
    let pose = it.pose(initial)?;
    let material = part
        .materials
        .first()
        .ok_or_else(|| it.invalid("a cell needs a material"))?;
    let models = ctx
        .cards
        .material(&material.card)
        .map_or(&[][..], MaterialCard::models);
    let (model, factory) = laws
        .model(models)
        .ok_or_else(|| it.invalid("no caloric law"))?;
    let [x, y, z] = [lx, ly, lz].map(|side| side.value);
    let (volume, area) = (x * y * z, 2.0 * (x * y + y * z + z * x));
    let mass = it.material(ctx, material, "density", DENSITY)? * volume;
    it.derive(part, "volume", volume, VOLUME_DIMS, "volume = lx * ly * lz");
    it.derive(
        part,
        "surface_area",
        area,
        AREA,
        "area = 2 (lx ly + ly lz + lz lx)",
    );
    it.derive(part, "mass", mass, MASS_DIMS, "mass = density * volume");
    let (item, env, temperature_k) = (it.item.as_str(), query_env(ctx), temperature.value);
    let implicit = ctx.draft.input().time_plan.scheme != TimeScheme::Explicit;
    let request = CaloricRequest {
        item,
        cards: ctx.cards,
        material,
        model,
        env,
        temperature_k,
        implicit,
    };
    let caloric = factory(&request, &mut it.r)?;
    if implicit && caloric.implicit.is_none() {
        return Err(it.invalid("its caloric law has no implicit form"));
    }
    let half = pose
        .orientation()
        .rotate_body_to_world(Vec3::new(0.5 * x, 0.0, 0.0));
    let p = pose.position_world();
    let geometry = CellGeometry {
        start_m: [p.x - half.x, p.y - half.y, p.z - half.z],
        end_m: [p.x + half.x, p.y + half.y, p.z + half.z],
        volume_m3: volume,
        surface_area_m2: area,
    };
    let cell = DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id),
        mass_kg: mass,
        specific_h_j_kg: caloric.specific_energy_j_kg,
        geometry,
        boundary,
    };
    Ok((cell, caloric))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sensible_law_inverts_its_own_specific_energy() {
        let law = SensibleHeat {
            c_j_kg_k: 500.0,
            datum_k: 0.0,
            band: None,
            k_w_m_k: None,
        };
        assert_eq!(
            law.temperature_k(500.0 * 300.0).unwrap().to_bits(),
            300.0_f64.to_bits()
        );
        let ds = law
            .specific_entropy_change_j_kg_k(500.0 * 300.0, 500.0 * 600.0)
            .unwrap();
        assert!((ds.unwrap() - 500.0 * 2.0_f64.ln()).abs() <= 1e-12);
        assert!(law.temperature_k(0.0).is_err());
    }
}
