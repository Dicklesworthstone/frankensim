//! Owner constructors: thin adapters that lower groups of draft items into
//! the physics owners' own Components, through their public constructors.
//!
//! - [`mbd`]: the rigid bodies joined by [`VISCOUS_BUSHING_FAMILY`]
//!   interfaces or card-less `ShearFilm`s lower to one `fs_mbd` `ViscousNetwork`,
//!   and those joined by [`RELAXATION_BUSHING_FAMILY`] ones to one `RelaxationNetwork`.
//! - [`thermal`]: thermal solids whose material card names a registered
//!   caloric law lower to the cells of one `fs_thermal` `ThermalTimeline`.
//! - Loss routes: a card-less `HeatRoute` interface from a part that one
//!   bushing joins to a cell carries its declared `fraction` of that
//!   bushing's heat; the rest goes to a sink named `unallocated`. A card-less
//!   `TemperatureSource` from such a part to a cell gives that bushing's
//!   temperature shift the cell's temperature at the start of each step.
//!
//! [`register_owner_defaults`] adds the three groups and overrides no kernel
//! default: a draft that names none of these families or laws lowers as
//! before.

use fs_matdb::PropertySample;
use fs_mbd::{Pose, UnitQuaternion, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftInitialState, DraftInterfaceDecl, DraftInterfaceKind,
    DraftPartDecl, ThirdBodyDecl,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey};
use fs_scenario::entity::{Entity, EntityId};

use super::{
    CardStore, GroupLoweringContext, GroupSelection, InterfaceCardEntry, LoweredGroupOutput,
    LoweredInterfaceOutput, LoweredPartOutput, LoweringError, LoweringGraph, LoweringRegistry,
    MatDbQueryEnv, ProvenanceRecorder, geometry_generator_label, resolve_geometry_fingerprint,
};
use crate::{HeatSink, PortRef, RouteTemperature, heat_sink};

pub mod gas;
pub mod mbd;
pub mod thermal;

/// Card family of a passive viscous bushing between two rigid bodies (a
/// `MechanicalContact` interface).
pub const VISCOUS_BUSHING_FAMILY: &str = "viscous-bushing";
/// Card family of a viscous bushing that also has Maxwell branches, each with
/// its own elastic memory (a `MechanicalContact` interface).
pub const RELAXATION_BUSHING_FAMILY: &str = "relaxation-bushing";
/// Card family of a cell's surface exchange with the ambient (a
/// `ThermalContact` interface between a cell's face and a face of a
/// stationary `DrivenBody` frame).
pub const THERMAL_BOUNDARY_FAMILY: &str = "thermal-boundary";
/// Group name of the loss routes.
pub const LOSS_ROUTES: &str = "owners.loss-routes";

/// A card's state, as declared on an interface.
type State<'a> = &'a [(ParameterKey, QtyAny)];
/// A bushing's interface entity, its heat port and its (cell, fraction) sinks.
type Route = (EntityId, PortRef, Vec<(PortRef, f64)>);

/// Register the owner groups with the product caloric laws.
pub fn register_owner_defaults(registry: &mut LoweringRegistry) {
    register_owners(registry, thermal::CaloricLaws::product());
}

/// Register the owner groups, resolving caloric laws through `caloric`.
/// They lower in this order: the viscous network, the thermal timeline, the
/// loss routes (which join the first two).
pub fn register_owners(registry: &mut LoweringRegistry, caloric: thermal::CaloricLaws) {
    registry.register_group(mbd::GROUP, mbd::select, mbd::lower);
    thermal::register(registry, caloric);
    gas::register(registry);
    registry.register_group(LOSS_ROUTES, select_loss_routes, lower_loss_routes);
}

/// The provenance of one draft item (`part:<name>`, `interface:<name>`).
pub(crate) struct Scoped {
    item: String,
    r: ProvenanceRecorder,
}

impl Scoped {
    pub(crate) fn new(kind: &str, name: &str) -> Self {
        let item = format!("{kind}:{name}");
        Self {
            item,
            r: ProvenanceRecorder::default(),
        }
    }

    pub(crate) fn invalid(&self, detail: impl Into<String>) -> LoweringError {
        LoweringError::Invalid {
            item: self.item.clone(),
            detail: detail.into(),
        }
    }

    /// Record a value declared on this item, named by its field.
    pub(crate) fn declare(&mut self, field: &str, value: QtyAny) {
        let item = self.item.as_str();
        self.r
            .record_declared_field(field, value.value, value.dims, item, field);
    }

    /// Record a value declared on another item (`environment`, `time_plan`).
    pub(crate) fn declare_on(&mut self, item: &str, field: &str, value: QtyAny) {
        self.r
            .record_declared_field(field, value.value, value.dims, item, field);
    }

    /// The declared initial pose, recorded.
    pub(crate) fn pose(&mut self, initial: &DraftInitialState) -> Result<Pose, LoweringError> {
        for (axis, value) in initial.position.iter().enumerate() {
            self.declare(&format!("initial.position[{axis}]"), *value);
        }
        for (axis, value) in initial.orientation.into_iter().enumerate() {
            let value = QtyAny::dimensionless(value);
            self.declare(&format!("initial.orientation[{axis}]"), value);
        }
        let [x, y, z] = initial.position.map(|value| value.value);
        let [w, i, j, k] = initial.orientation;
        let pose = UnitQuaternion::new(w, i, j, k).and_then(|q| Pose::new(Vec3::new(x, y, z), q));
        pose.map_err(|err| self.invalid(format!("initial pose: {err:?}")))
    }

    /// Record a value derived from `part`'s geometry by `formula`.
    pub(crate) fn derive(
        &mut self,
        part: &DraftPartDecl,
        name: &str,
        value: f64,
        dims: Dims,
        formula: &str,
    ) {
        let hash = resolve_geometry_fingerprint(&part.geometry);
        let label = geometry_generator_label(&part.geometry);
        self.r
            .record_geometry_derived(name, value, dims, hash, label, formula);
    }

    /// Property `name` of `material`, recorded, in dimensions `dims`.
    pub(crate) fn material(
        &mut self,
        ctx: &GroupLoweringContext<'_>,
        material: &MaterialRef,
        name: &str,
        dims: Dims,
    ) -> Result<f64, LoweringError> {
        let (r, env) = (&mut self.r, query_env(ctx));
        let (q, _) = r.query_material_property(&self.item, ctx.cards, material, name, name, env)?;
        self.checked(name, &q, dims)
    }

    /// Property `name` of an interface card at its state, recorded, in `dims`.
    pub(crate) fn interface(
        &mut self,
        ctx: &GroupLoweringContext<'_>,
        (entry, state): (&InterfaceCardEntry, State<'_>),
        name: &str,
        dims: Dims,
    ) -> Result<f64, LoweringError> {
        let (r, env) = (&mut self.r, query_env(ctx));
        let (q, _) = r.query_interface_property(&self.item, entry, state, name, name, env)?;
        self.checked(name, &q, dims)
    }

    /// Property `name` of a third body's card at its query point, recorded.
    pub(crate) fn third_body(
        &mut self,
        ctx: &GroupLoweringContext<'_>,
        tb: &ThirdBodyDecl,
        name: &str,
        dims: Dims,
    ) -> Result<f64, LoweringError> {
        let (r, env) = (&mut self.r, query_env(ctx));
        let (q, _) = r.query_third_body_property(&self.item, ctx.cards, tb, name, name, env)?;
        self.checked(name, &q, dims)
    }

    fn checked(&self, name: &str, q: &PropertySample, dims: Dims) -> Result<f64, LoweringError> {
        if q.dims == dims {
            return Ok(q.value);
        }
        Err(self.invalid(format!(
            "`{name}` has dimensions {:?}, not {dims:?}",
            q.dims
        )))
    }

    pub(crate) fn part(self, node: &str) -> LoweredPartOutput {
        let (property_receipts, parameter_sources) = self.r.finish();
        let created_nodes = vec![node.to_owned()];
        LoweredPartOutput {
            created_nodes,
            property_receipts,
            parameter_sources,
        }
    }

    pub(crate) fn interface_output(self, node: Option<&str>) -> LoweredInterfaceOutput {
        let (property_receipts, parameter_sources) = self.r.finish();
        let created_nodes = node.map(str::to_owned).into_iter().collect();
        LoweredInterfaceOutput {
            created_nodes,
            property_receipts,
            parameter_sources,
        }
    }
}

fn query_env(ctx: &GroupLoweringContext<'_>) -> MatDbQueryEnv {
    let base = &ctx.draft.input().environment.base;
    MatDbQueryEnv {
        ambient_temp_k: base.ambient_temperature.value,
        ambient_pressure_pa: base.ambient_pressure.value,
        policy: ctx.selection_policy,
    }
}

/// The index of the part that owns `entity`: the part itself, or the part a
/// feature of it belongs to.
pub(crate) fn owning_part(input: &DraftAssemblyInput, entity: EntityId) -> Option<usize> {
    let mut current = Some(entity);
    for _ in 0..=input.catalog.entities().len() {
        let id = current?;
        if let Some(index) = input.parts.iter().position(|part| part.entity == id) {
            return Some(index);
        }
        current = input.catalog.get(id).and_then(Entity::parent);
    }
    None
}

/// The position in `parts` of the part that owns `entity`.
pub(crate) fn slot(input: &DraftAssemblyInput, parts: &[usize], e: EntityId) -> Option<usize> {
    let part = owning_part(input, e)?;
    parts.iter().position(|&p| p == part)
}

/// The card entry and state of `iface` when its card is of `family`: a
/// `MechanicalContact` for a bushing, a `ThermalContact` otherwise.
pub(crate) fn family_entry<'a>(
    cards: &'a CardStore,
    iface: &'a DraftInterfaceDecl,
    family: &str,
) -> Option<(&'a InterfaceCardEntry, State<'a>)> {
    let bushing = [VISCOUS_BUSHING_FAMILY, RELAXATION_BUSHING_FAMILY].contains(&family);
    let (card, state) = match &iface.kind {
        DraftInterfaceKind::MechanicalContact {
            law_card, state, ..
        } if bushing => (law_card, state),
        DraftInterfaceKind::ThermalContact {
            law_card, state, ..
        } if !bushing => (law_card, state),
        _ => return None,
    };
    let entry = cards
        .interface(card)
        .filter(|entry| entry.family == family)?;
    Some((entry, state.as_slice()))
}

/// The state of `iface` when it is a heat route.
pub(crate) fn route_state(iface: &DraftInterfaceDecl) -> Option<State<'_>> {
    match &iface.kind {
        DraftInterfaceKind::HeatRoute { state } => Some(state.as_slice()),
        _ => None,
    }
}

/// Whether `iface` is a bushing: a `viscous-bushing` or `relaxation-bushing`
/// card, or a shear film.
pub(crate) fn is_bushing(cards: &CardStore, iface: &DraftInterfaceDecl) -> bool {
    matches!(iface.kind, DraftInterfaceKind::ShearFilm { .. })
        || family_entry(cards, iface, VISCOUS_BUSHING_FAMILY).is_some()
        || family_entry(cards, iface, RELAXATION_BUSHING_FAMILY).is_some()
}

/// The bushing a loss route draws from: the one bushing that joins the part
/// owning the route's `feature_a`.
pub(crate) fn route_bushing<'a>(
    input: &'a DraftAssemblyInput,
    cards: &CardStore,
    route: &DraftInterfaceDecl,
) -> Option<&'a DraftInterfaceDecl> {
    let part = owning_part(input, route.feature_a)?;
    let on_part = |end| owning_part(input, end) == Some(part);
    let mut joined = input
        .interfaces
        .iter()
        .filter(|b| (on_part(b.feature_a) || on_part(b.feature_b)) && is_bushing(cards, b));
    let bushing = joined.next()?;
    joined.next().is_none().then_some(bushing)
}

/// The draft's interfaces of `family`, as a selection.
pub(crate) fn of_family(draft: &DraftAssembly, cards: &CardStore, family: &str) -> GroupSelection {
    let interfaces = draft.input().interfaces.iter().enumerate();
    let interfaces = interfaces.filter(|(_, iface)| family_entry(cards, iface, family).is_some());
    let interfaces = interfaces.map(|(index, _)| index).collect();
    GroupSelection {
        interfaces,
        ..GroupSelection::default()
    }
}

/// The sink of a bushing's heat that no route takes.
pub(crate) fn unallocated(ambient_k: f64) -> Result<HeatSink, LoweringError> {
    heat_sink("unallocated", ambient_k).map_err(|err| LoweringError::Invalid {
        item: "environment".to_owned(),
        detail: format!("no unallocated sink at the ambient temperature: {err:?}"),
    })
}

#[allow(clippy::unnecessary_wraps)] // Signature required by GroupSelector.
fn select_loss_routes(
    draft: &DraftAssembly,
    _cards: &CardStore,
) -> Result<Vec<GroupSelection>, LoweringError> {
    let mut selection = GroupSelection::default();
    let routes = draft.input().interfaces.iter().enumerate();
    let source = |iface: &DraftInterfaceDecl| iface.kind == DraftInterfaceKind::TemperatureSource;
    let routes = routes.filter(|(_, iface)| route_state(iface).is_some() || source(iface));
    let routes = routes.map(|(index, _)| index);
    selection.interfaces = routes.collect();
    let found = !selection.interfaces.is_empty();
    Ok(found.then_some(selection).into_iter().collect())
}

/// One heat route per bushing, in the order of its first loss route: the
/// declared fractions to the cells, the rest to `unallocated` at ambient; and
/// each temperature source's signal.
fn lower_loss_routes(
    ctx: &GroupLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredGroupOutput, LoweringError> {
    let input = ctx.draft.input();
    let mut out = LoweredGroupOutput::default();
    let mut routes: Vec<Route> = Vec::new();
    for &index in &ctx.selection.interfaces {
        let iface = &input.interfaces[index];
        let mut it = Scoped::new("interface", &iface.name);
        let bushing = route_bushing(input, ctx.cards, iface)
            .ok_or_else(|| it.invalid("feature_a must be on a part exactly one bushing joins"))?;
        if iface.kind == DraftInterfaceKind::TemperatureSource {
            connect_temperature(ctx, graph, (&mut it, iface, bushing))?;
            out.interfaces.push(it.interface_output(None));
            continue;
        }
        let source = graph.entity_port(bushing.entity, "heat");
        let source =
            source.ok_or_else(|| it.invalid("its bushing must be a lowered viscous bushing"))?;
        let receiver = graph.resolve_entity_port(&input.catalog, iface.feature_b, "heat");
        let receiver = receiver
            .ok_or_else(|| it.invalid("feature_b must be a cell of the thermal timeline"))?;
        let fraction = route_state(iface)
            .unwrap_or_default()
            .iter()
            .find(|(key, _)| key.name == "fraction")
            .map(|(_, q)| *q);
        let Some(fraction) =
            fraction.filter(|q| q.dims == Dims::NONE && (0.0..=1.0).contains(&q.value))
        else {
            return Err(it.invalid("a loss route declares a dimensionless `fraction` in [0, 1]"));
        };
        it.declare("state.fraction", fraction);
        out.interfaces.push(it.interface_output(None));
        match routes
            .iter_mut()
            .find(|(entity, ..)| *entity == bushing.entity)
        {
            Some((_, _, sinks)) => sinks.push((receiver, fraction.value)),
            None => routes.push((bushing.entity, source, vec![(receiver, fraction.value)])),
        }
    }
    let ambient_k = input.environment.base.ambient_temperature.value;
    for (_, source, sinks) in routes {
        let sink = Some(unallocated(ambient_k)?);
        graph.add_heat_route(source, sinks, sink, RouteTemperature::SinkEquivalent);
    }
    Ok(out)
}

/// A temperature source's cell's published temperature as the signal its
/// bushing's shift reads, guessed at first as the cell's initial temperature
/// (the ambient when it declares none), which is recorded.
fn connect_temperature(
    ctx: &GroupLoweringContext<'_>,
    graph: &mut LoweringGraph,
    (it, iface, bushing): (&mut Scoped, &DraftInterfaceDecl, &DraftInterfaceDecl),
) -> Result<(), LoweringError> {
    let input = ctx.draft.input();
    let law = graph.entity_port(bushing.entity, "temperature");
    let law = law.ok_or_else(|| it.invalid("its bushing's card states no temperature shift"))?;
    let cell = graph.resolve_entity_port(&input.catalog, iface.feature_b, "temperature");
    let cell =
        cell.ok_or_else(|| it.invalid("feature_b must be a cell of the thermal timeline"))?;
    graph.builder_mut().signal(cell, law)?;
    let part = owning_part(input, iface.feature_b).map(|part| &input.parts[part]);
    let start = part.and_then(|part| Some((part, part.initial.as_ref()?.temperature?)));
    let start = if let Some((part, temperature)) = start {
        let item = format!("part:{}", part.name);
        it.declare_on(&item, "initial.temperature", temperature);
        temperature
    } else {
        let ambient = input.environment.base.ambient_temperature;
        it.declare_on("environment", "base.ambient_temperature", ambient);
        ambient
    };
    graph.set_initial_guess(law, vec![start.value]);
    Ok(())
}
