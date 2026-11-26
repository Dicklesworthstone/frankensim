//! Thermal solids whose material card names a registered caloric law,
//! lowered to the cells of one `fs_thermal` `ThermalTimeline`.
//!
//! A part is a cell when it is a `Solid` with only thermal physics and a
//! model on its first material card names a law [`CaloricLaws`] knows. A
//! cell is a box, born at the start: its mass is density times volume, its
//! centreline the box's x extent at its declared pose, its specific energy
//! the law's at its declared temperature (or the ambient). A
//! `thermal-boundary` interface from the cell's face to a face of a
//! stationary frame (a `DrivenBody` part) gives its surface exchange at ambient (`convection_coefficient`, `reservoir_area`
//! and `reservoir_conductance` on its card); without one it exchanges
//! nothing. A constant source on its [`HEAT_SOURCE_PORT`] is its applied
//! heat. The timeline steps at most the time plan's step, checks energy to
//! the cells' shared absolute tolerance, and sends boundary heat through its
//! `environment` port to a sink of that name.

use std::collections::BTreeMap;
use std::rc::Rc;

use fs_matdb::{ConstitutiveModelCard, MaterialCard};
use fs_mbd::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::MaterialRef;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, DraftPartDecl, DraftPartKind, PhysicsDomain,
    PrimitiveShape, SourceTarget, StepPolicy, TimeHistory,
};
use fs_thermal::component::ThermalTimeline;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};

use super::{
    Scoped, THERMAL_BOUNDARY_FAMILY, family_entry, of_family, owning_part, query_env, slot,
};
use crate::lowering::ConstitutiveModelRef;
use crate::lowering::{
    CONDUCTANCE_DIMS, CardStore, GroupLoweringContext, GroupSelection, LoweredGroupOutput,
    LoweringError, LoweringGraph, LoweringRegistry, MASS_DIMS, MatDbQueryEnv, POWER_DIMS,
    ProvenanceRecorder, SourceLoweringReceipt, TEMPERATURE_DIMS, VOLUME_DIMS,
    record_time_history_provenance,
};
use crate::{RouteTemperature, heat_sink};

/// Group name of the thermal timeline.
pub const GROUP: &str = "owners.thermal.timeline";
/// Law id of the product caloric law: the card's constant specific heat
/// from the model's `datum_temperature`, at which the specific energy is zero.
pub const SENSIBLE_HEAT_LAW: &str = "sensible-heat";
/// The port of a cell its applied heat source targets.
pub const HEAT_SOURCE_PORT: &str = "heat_source";
const TIMELINE: &str = "thermal-timeline";
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const HEAT_TRANSFER: Dims = Dims([0, 1, -3, -1, 0, 0]);
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
}

/// A cell's caloric law and its initial specific energy [J/kg].
pub struct CaloricCell {
    /// The law, held by the timeline itself.
    pub law: Rc<dyn ThermalLaw>,
    /// Specific energy at the initial temperature, on the law's own datum.
    pub specific_energy_j_kg: f64,
}

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
/// `e = c (T - t0)`.
struct SensibleHeat {
    c_j_kg_k: f64,
    datum_k: f64,
}

impl ThermalLaw for SensibleHeat {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let temperature_k = self.datum_k + specific_h_j_kg / self.c_j_kg_k;
        if !(temperature_k.is_finite() && temperature_k > 0.0) {
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

/// The product law: the card's `specific_heat_capacity`, from the model's
/// `datum_temperature`.
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
    let (c_j_kg_k, datum_k) = (c.value, datum.value);
    Ok(CaloricCell {
        law: Rc::new(SensibleHeat { c_j_kg_k, datum_k }),
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

/// The cells, their boundary interfaces and their heat sources, as one
/// selection.
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
    let (mut schedule, mut cell_laws) = (Vec::new(), Vec::new());
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
        let source_w = sources.get(&cell).copied().unwrap_or(0.0);
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
                    reservoir: Some(reservoir),
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
        let (born, law) = new_cell(laws, ctx, &mut it, part, cell as u64, boundary)?;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(born),
        });
        cell_laws.push((MaterialId(cell as u64), law));
        out.parts.push(it.part(TIMELINE));
    }
    let cells = cell_laws.len();
    let fail =
        |err: &dyn std::fmt::Debug| Scoped::new("node", TIMELINE).invalid(format!("{err:?}"));
    let limits = DepositionBudget {
        max_cells: cells,
        max_links: cells,
        max_events: cells,
        max_substeps: work,
    };
    let simulation = DepositionSimulation::new_shared(schedule, &cell_laws, limits);
    let simulation = simulation.and_then(|start| start.advance_to(0.0, config, || Ok(())));
    let simulation = simulation.map_err(|err| fail(&err))?.0;
    let environment = !boundaries.is_empty() || !sources.is_empty();
    let ids = (0..cells as u64).map(CellId).collect();
    let timeline = ThermalTimeline::new(ids, environment, config).map_err(|err| fail(&err))?;
    let start = timeline.start(&simulation).map_err(|err| fail(&err))?;
    let node = graph.add(TIMELINE, Some(first.entity), timeline, start)?;
    for (cell, &index) in ctx.selection.parts.iter().enumerate() {
        let entity = input.parts[index].entity;
        if cell > 0 {
            graph.entity_nodes.entry(entity).or_default().push(node);
        }
        let port = graph.port(node, &format!("cell{cell}"))?;
        graph.expose_port(entity, "heat", port);
    }
    if environment {
        let port = graph.port(node, "environment")?;
        let sink = heat_sink("environment", ambient_k).map_err(|err| fail(&err))?;
        graph.add_heat_route(
            port,
            Vec::new(),
            Some(sink),
            RouteTemperature::Fixed(ambient_k),
        );
    }
    Ok(out)
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
        let card = family_entry(ctx.cards, iface, THERMAL_BOUNDARY_FAMILY);
        let card = card.ok_or_else(|| it.invalid("not a thermal-boundary interface"))?;
        let exchange = [
            it.interface(ctx, card, "convection_coefficient", HEAT_TRANSFER)?,
            it.interface(ctx, card, "reservoir_area", AREA)?,
            it.interface(ctx, card, "reservoir_conductance", CONDUCTANCE_DIMS)?,
        ];
        let ambient = input.environment.base.ambient_temperature;
        it.declare_on("environment", "base.ambient_temperature", ambient);
        if boundaries.insert(cell, exchange).is_some() {
            return Err(it.invalid("a cell has one thermal boundary"));
        }
        out.interfaces.push(it.interface_output(Some(TIMELINE)));
    }
    Ok(boundaries)
}

/// Each source's cell and constant power.
fn sources(
    ctx: &GroupLoweringContext<'_>,
    out: &mut LoweredGroupOutput,
) -> Result<BTreeMap<usize, f64>, LoweringError> {
    let input = ctx.draft.input();
    let mut sources = BTreeMap::new();
    for &index in &ctx.selection.sources {
        let source = &input.sources[index];
        let mut it = Scoped::new("source", &source.name);
        let SourceTarget::Port { entity, .. } = &source.target else {
            return Err(it.invalid("a cell's heat source targets its port"));
        };
        let cell = slot(input, &ctx.selection.parts, *entity);
        let (TimeHistory::Constant(power), Some(cell)) = (&source.history, cell) else {
            return Err(it.invalid("a cell's heat source is constant"));
        };
        if power.dims != POWER_DIMS || sources.insert(cell, power.value).is_some() {
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

/// A cell born at the start, and its law.
fn new_cell(
    laws: &CaloricLaws,
    ctx: &GroupLoweringContext<'_>,
    it: &mut Scoped,
    part: &DraftPartDecl,
    id: u64,
    boundary: Boundary,
) -> Result<(DepositedCell, Rc<dyn ThermalLaw>), LoweringError> {
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
    let request = CaloricRequest {
        item,
        cards: ctx.cards,
        material,
        model,
        env,
        temperature_k,
    };
    let caloric = factory(&request, &mut it.r)?;
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
    Ok((cell, caloric.law))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sensible_law_inverts_its_own_specific_energy() {
        let law = SensibleHeat {
            c_j_kg_k: 500.0,
            datum_k: 0.0,
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
