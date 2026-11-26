//! The rigid bodies joined by `viscous-bushing` interfaces or shear films,
//! lowered to one `fs_mbd` `ViscousNetwork`; or joined by `relaxation-bushing`
//! interfaces, lowered to one `RelaxationNetwork`. A network's bushings are all
//! relaxation bushings or none.
//!
//! A body's mass and principal inertia come from its primitive geometry and
//! the density on its card, about its centroid, which is its body origin;
//! its pose and motion are its declared initial state. A bushing's damping
//! is the diagonal SI form on its card (`damping_x`, `damping_y`,
//! `damping_z` in N s/m; `damping_rx`, `damping_ry`, `damping_rz` in
//! N m s/rad), acting at both bodies' origins. A card-less shear film's
//! normal is the bushing's z axis: `mu A / h` on x and y, `mu J / h` on rz
//! with `J = A^2 / (2 pi)` (a full disk), and 0 on z, rx and ry, since
//! squeeze is not linear. `mu` is its third body's at the declared query
//! temperature and stays there (one-way coupling). The network checks impulse,
//! work and speed to the bodies' shared absolute tolerance, relaxes at most
//! their `max_work_units` pairs, and takes one source substep per step. Its
//! body forces are bound to zero; a bushing's heat goes where its loss
//! routes send it, all of it to `unallocated` when it has none.
//!
//! A relaxation bushing has a viscous bushing's damping and, for each branch
//! `k = 0, 1, ...` while its card has `relaxation_time_<k>` (s), a diagonal
//! stiffness `stiffness_<k>_x` .. `stiffness_<k>_rz` (N/m, N m/rad), whose
//! factor `L` is its square root. A branch's memory has no default: the
//! interface's state declares its elastic displacement `displacement_<k>_x`
//! .. `displacement_<k>_rz` (m, rad), held as `a = L x`. A relaxation network
//! prepares at most `max_work_units` dense multiply-adds per step.

use core::f64::consts::PI;

use fs_mbd::component::network::ViscousNetwork;
use fs_mbd::component::relaxation::{RelaxationNetwork, RelaxationState};
use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::{
    RelaxationNetworkLimits, RelaxationNetworkLink, SpatialRelaxationBranch,
};
use fs_mbd::spatial_loads::internal_impulse::viscous::network::{
    ViscousNetworkLimits, ViscousNetworkLink,
};
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::{
    ViscousNetworkMotionControls, relaxation::RelaxationNetworkMotionControls,
};
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::ParameterKey;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl,
    DraftPartKind, PrimitiveShape, ThirdBodyDecl,
};
use fs_scenario::entity::EntityId;

use super::{
    RELAXATION_BUSHING_FAMILY, Scoped, State, VISCOUS_BUSHING_FAMILY, family_entry, is_bushing,
    owning_part, route_bushing, route_state, slot, unallocated,
};
use crate::lowering::{
    AREA_DIMS, CardStore, DAMPING_DIMS, GroupLoweringContext, GroupSelection, INERTIA_DIMS,
    InterfaceCardEntry, LENGTH_DIMS, LoweredGroupOutput, LoweringError, LoweringGraph, MASS_DIMS,
    Observed, STIFFNESS_DIMS, TEMPERATURE_DIMS, TIME_DIMS, TORSIONAL_DAMPING_DIMS,
    TORSIONAL_STIFFNESS_DIMS, VOLUME_DIMS, geometry_volume_m3,
};
use crate::{NodeId, RouteTemperature, System};

/// Group name of the bushing network.
pub const GROUP: &str = "owners.mbd.viscous-network";
const NETWORK: &str = "viscous-network";
const RELAXATION: &str = "relaxation-network";
/// Each connector's branches and their memory amplitudes [sqrt(J)].
type Memory = Vec<(Vec<SpatialRelaxationBranch>, Vec<[f64; 6]>)>;
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const VISCOSITY: Dims = Dims([-1, 1, -1, 0, 0, 0]);
const POLAR_MOMENT: Dims = Dims([4, 0, 0, 0, 0, 0]);
/// A bushing card's diagonal coefficients, in the matrix's order.
const DAMPING: [(&str, Dims); 6] = [
    ("damping_x", DAMPING_DIMS),
    ("damping_y", DAMPING_DIMS),
    ("damping_z", DAMPING_DIMS),
    ("damping_rx", TORSIONAL_DAMPING_DIMS),
    ("damping_ry", TORSIONAL_DAMPING_DIMS),
    ("damping_rz", TORSIONAL_DAMPING_DIMS),
];
/// A branch's coordinates: key suffix, stiffness and displacement dimensions.
const BRANCH: [(&str, Dims, Dims); 6] = [
    ("x", STIFFNESS_DIMS, LENGTH_DIMS),
    ("y", STIFFNESS_DIMS, LENGTH_DIMS),
    ("z", STIFFNESS_DIMS, LENGTH_DIMS),
    ("rx", TORSIONAL_STIFFNESS_DIMS, Dims::NONE),
    ("ry", TORSIONAL_STIFFNESS_DIMS, Dims::NONE),
    ("rz", TORSIONAL_STIFFNESS_DIMS, Dims::NONE),
];
/// Working memory one network step may hold: a resource cap, not physics.
const LIVE_BYTES: usize = 64 << 20;

/// The bushings, and every part one joins, as one selection.
pub(crate) fn select(
    draft: &DraftAssembly,
    cards: &CardStore,
) -> Result<Vec<GroupSelection>, LoweringError> {
    let input = draft.input();
    let mut selection = GroupSelection::default();
    let bushings = input.interfaces.iter().enumerate();
    let bushings = bushings.filter(|(_, iface)| is_bushing(cards, iface));
    selection.interfaces = bushings.map(|(index, _)| index).collect();
    for &index in &selection.interfaces {
        let iface = &input.interfaces[index];
        let ends = [iface.feature_a, iface.feature_b].map(|end| owning_part(input, end));
        let [Some(a), Some(b)] = ends else {
            let it = Scoped::new("interface", &iface.name);
            return Err(it.invalid("a viscous bushing joins two parts"));
        };
        selection.parts.extend([a, b]);
    }
    selection.parts.sort_unstable();
    selection.parts.dedup();
    let found = !selection.parts.is_empty();
    Ok(found.then_some(selection).into_iter().collect())
}

/// One network of the selected bodies and bushings: a `RelaxationNetwork` when
/// the bushings have memory, a `ViscousNetwork` otherwise.
pub(crate) fn lower(
    ctx: &GroupLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredGroupOutput, LoweringError> {
    let input = ctx.draft.input();
    let parts = &ctx.selection.parts;
    let first = &input.parts[parts[0]];
    let budget = &first.fidelity.error_budget;
    let tolerance = budget.absolute_tol.value;
    let work = budget
        .max_work_units
        .and_then(|work| usize::try_from(work).ok());
    let Some(work) = work.filter(|_| tolerance > 0.0 && tolerance.is_finite()) else {
        let it = Scoped::new("part", &first.name);
        return Err(it.invalid("a network needs a positive absolute_tol and max_work_units"));
    };
    let memory = with_memory(ctx)?;
    let label = if memory { RELAXATION } else { NETWORK };
    let mut out = LoweredGroupOutput::default();
    let (mut properties, mut states) = (Vec::new(), Vec::new());
    for &index in parts {
        let part = &input.parts[index];
        let mut it = Scoped::new("part", &part.name);
        if part.kind != DraftPartKind::RigidBody || part.fidelity.error_budget != *budget {
            return Err(it.invalid("a network's bodies are rigid bodies with one error budget"));
        }
        it.declare("fidelity.error_budget.absolute_tol", budget.absolute_tol);
        let work_units = QtyAny::dimensionless(work as f64);
        it.declare("fidelity.error_budget.max_work_units", work_units);
        for (axis, g) in input.environment.base.gravity.iter().enumerate() {
            it.declare_on("environment", &format!("base.gravity[{axis}]"), *g);
        }
        let mass = mass_properties(ctx, &mut it, part)?;
        states.push(initial_state(&mut it, part, mass)?);
        properties.push(mass);
        out.parts.push(it.part(label));
    }
    let (mut links, mut bank) = (Vec::new(), Vec::new());
    for &index in &ctx.selection.interfaces {
        let iface = &input.interfaces[index];
        let mut it = Scoped::new("interface", &iface.name);
        let mut coefficients = [[0.0; 6]; 6];
        if let DraftInterfaceKind::ShearFilm { third_body, state } = &iface.kind {
            film(ctx, &mut it, (third_body, state), &mut coefficients)?;
        } else {
            let family = [VISCOUS_BUSHING_FAMILY, RELAXATION_BUSHING_FAMILY][usize::from(memory)];
            let card = family_entry(ctx.cards, iface, family);
            let card = card.ok_or_else(|| it.invalid(format!("not a {family} interface")))?;
            for (axis, (name, dims)) in DAMPING.into_iter().enumerate() {
                coefficients[axis][axis] = it.interface(ctx, card, name, dims)?;
            }
            if memory {
                bank.push(branches(ctx, &mut it, card)?);
            }
        }
        // The card's coefficients are SI, so the coordinates are SI units.
        let scales = BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        };
        let damping = PassiveBushingMatrix::try_new(coefficients, scales);
        let damping = damping.map_err(|err| it.invalid(format!("not passive: {err:?}")))?;
        let ends = [iface.feature_a, iface.feature_b].map(|end| slot(input, parts, end));
        let [Some(a), Some(b)] = ends else {
            return Err(it.invalid("a viscous bushing joins two parts"));
        };
        let arms_body_m = [Vec3::ZERO; 2];
        links.push(ViscousNetworkLink {
            bodies: [a, b],
            arms_body_m,
            damping,
        });
        out.interfaces.push(it.interface_output(Some(label)));
    }
    let [x, y, z] = input.environment.base.gravity.map(|g| g.value);
    let gravity = Gravity::new(Vec3::new(x, y, z));
    let gravity =
        gravity.map_err(|err| Scoped::new("environment", "gravity").invalid(format!("{err:?}")))?;
    let controls = controls(tolerance, properties.len(), links.len(), work);
    let entities: Vec<EntityId> = parts.iter().map(|&p| input.parts[p].entity).collect();
    let (bodies, bank) = ((properties, states), memory.then_some(bank));
    let node = add_network(graph, &entities, bodies, (links, bank), (gravity, controls))?;
    for &entity in &entities[1..] {
        graph.entity_nodes.entry(entity).or_default().push(node);
    }
    for body in 0..entities.len() {
        let force = graph.port(node, &format!("force{body}"))?;
        graph.bind_zero_boundary(force);
    }
    let ambient_k = input.environment.base.ambient_temperature.value;
    for (link, &index) in ctx.selection.interfaces.iter().enumerate() {
        let bushing = input.interfaces[index].entity;
        let heat = graph.port(node, &format!("heat{link}"))?;
        graph.expose_port(bushing, "heat", heat);
        let routed = input.interfaces.iter().any(|route| {
            route_state(route).is_some()
                && route_bushing(input, ctx.cards, route).is_some_and(|b| b.entity == bushing)
        });
        if !routed {
            let sink = Some(unallocated(ambient_k)?);
            graph.add_heat_route(heat, Vec::new(), sink, RouteTemperature::SinkEquivalent);
        }
    }
    Ok(out)
}

/// The network as one node, named for its first body, with each body's
/// reader: a `RelaxationNetwork` when the links have memory, else a
/// `ViscousNetwork`. A relaxation network's live bytes are the viscous one's,
/// and its preparation terms what the viscous one has as pair relaxations.
fn add_network(
    graph: &mut LoweringGraph,
    entities: &[EntityId],
    (properties, states): (Vec<MassProperties>, Vec<RigidBodyState>),
    (links, memory): (Vec<ViscousNetworkLink>, Option<Memory>),
    (gravity, controls): (Gravity, ViscousNetworkMotionControls),
) -> Result<NodeId, LoweringError> {
    let masses = properties.clone();
    let Some(memory) = memory else {
        let network = ViscousNetwork::new(properties, links, gravity, controls);
        let network =
            network.map_err(|err| Scoped::new("network", NETWORK).invalid(format!("{err:?}")))?;
        let node = graph.add(NETWORK, Some(entities[0]), network, states)?;
        observe_bodies(
            graph,
            (node, entities),
            masses,
            Vec::<RigidBodyState>::as_slice,
        );
        return Ok(node);
    };
    let ViscousNetworkLimits {
        maximum_bodies,
        maximum_connectors,
        maximum_pair_relaxations,
        maximum_live_bytes,
    } = controls.limits;
    let source = DrivenViscousNetworkLimits {
        maximum_bodies,
        maximum_connectors,
        maximum_live_bytes,
        maximum_preparation_terms: maximum_pair_relaxations as u64,
    };
    let maximum_branches = memory.iter().map(|(branches, _)| branches.len()).sum();
    let limits = RelaxationNetworkLimits {
        source,
        maximum_branches,
    };
    let controls = RelaxationNetworkMotionControls {
        source: controls.source,
        limits,
    };
    let links = links.into_iter().zip(memory);
    let (links, memory) = links
        .map(|(viscous, (branches, memory))| (RelaxationNetworkLink { viscous, branches }, memory))
        .unzip();
    let network = RelaxationNetwork::new(properties, links, gravity, controls);
    let network =
        network.map_err(|err| Scoped::new("network", RELAXATION).invalid(format!("{err:?}")))?;
    let state = RelaxationState {
        bodies: states,
        memory,
    };
    let node = graph.add(RELAXATION, Some(entities[0]), network, state)?;
    observe_bodies(
        graph,
        (node, entities),
        masses,
        |state: &RelaxationState| &state.bodies,
    );
    Ok(node)
}

/// Whether the bushings have memory: all on `relaxation-bushing` cards, or
/// none. A network that mixes them is refused, naming both sets.
fn with_memory(ctx: &GroupLoweringContext<'_>) -> Result<bool, LoweringError> {
    let input = ctx.draft.input();
    let bushings = ctx.selection.interfaces.iter();
    let bushings = bushings.map(|&index| &input.interfaces[index]);
    let has = |iface: &&DraftInterfaceDecl| {
        family_entry(ctx.cards, iface, RELAXATION_BUSHING_FAMILY).is_some()
    };
    let (memory, none): (Vec<_>, Vec<_>) = bushings.partition(has);
    if memory.is_empty() || none.is_empty() {
        return Ok(!memory.is_empty());
    }
    let names = |set: Vec<&DraftInterfaceDecl>| {
        set.iter()
            .map(|i| format!("`{}`", i.name))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let (memory, none) = (names(memory), names(none));
    let [with, without] = [RELAXATION_BUSHING_FAMILY, VISCOUS_BUSHING_FAMILY];
    let mixed = format!("{memory} on `{with}`, {none} on `{without}` or a shear film");
    let it = Scoped::new("network", RELAXATION);
    Err(it.invalid(format!(
        "a network's bushings are all on `{with}` or none: {mixed}"
    )))
}

/// A relaxation bushing's branches and their memory, from its card and its
/// interface's state (see the module's notes).
fn branches(
    ctx: &GroupLoweringContext<'_>,
    it: &mut Scoped,
    (card, state): (&InterfaceCardEntry, State<'_>),
) -> Result<(Vec<SpatialRelaxationBranch>, Vec<[f64; 6]>), LoweringError> {
    let (mut branches, mut memory) = (Vec::new(), Vec::new());
    loop {
        let k = branches.len();
        let time = format!("relaxation_time_{k}");
        if card.card.claims_for(&time).is_empty() {
            return Ok((branches, memory));
        }
        let tau = it.interface(ctx, (card, state), &time, TIME_DIMS)?;
        let (mut rows, mut x) = ([[0.0; 6]; 6], [0.0; 6]);
        for (axis, (suffix, stiffness, length)) in BRANCH.into_iter().enumerate() {
            let name = format!("stiffness_{k}_{suffix}");
            rows[axis][axis] = it.interface(ctx, (card, state), &name, stiffness)?.sqrt();
            let key = format!("displacement_{k}_{suffix}");
            let declared = state
                .iter()
                .find(|(p, q)| p.name == key && q.dims == length);
            let no_default =
                || it.invalid(format!("a branch's memory has no default: declare `{key}`"));
            let &(_, q) = declared.ok_or_else(no_default)?;
            it.declare(&format!("state.{key}"), q);
            x[axis] = q.value;
        }
        let branch = SpatialRelaxationBranch::from_gram_rows(rows, tau);
        let branch = branch.map_err(|err| it.invalid(format!("branch {k}: {err:?}")))?;
        let a = branch.elastic_memory(x);
        memory.push(a.map_err(|err| it.invalid(format!("branch {k} memory: {err:?}")))?);
        branches.push(branch);
    }
}

/// Registers each body's reader for probes: its velocity and angular velocity
/// in the assembly frame, from its slot in the network's state `S`.
fn observe_bodies<S: 'static>(
    graph: &mut LoweringGraph,
    (node, entities): (NodeId, &[EntityId]),
    masses: Vec<MassProperties>,
    bodies: fn(&S) -> &[RigidBodyState],
) {
    for (body, (&entity, mass)) in entities.iter().zip(masses).enumerate() {
        let read = move |system: &System| {
            let state = *bodies(system.state::<S>(node)?).get(body)?;
            let velocity = state.center_of_mass_velocity_world(mass).ok()?;
            let omega = mass.angular_velocity_body_checked(state.angular_momentum_body());
            let omega = state.pose().orientation().rotate_body_to_world(omega.ok()?);
            Some(Observed::Body {
                velocity_m_s: [velocity.x, velocity.y, velocity.z],
                angular_velocity_rad_s: [omega.x, omega.y, omega.z],
            })
        };
        graph.observers.insert(entity, Box::new(read));
    }
}

/// The network's controls: every tolerance at the bodies' shared absolute
/// tolerance, one source substep, and limits sized to this network.
fn controls(
    tolerance: f64,
    bodies: usize,
    connectors: usize,
    work: usize,
) -> ViscousNetworkMotionControls {
    ViscousNetworkMotionControls {
        source: ViscousConnectorControls {
            impulse: InternalWrenchImpulseControls {
                angular_origin_world_m: Vec3::ZERO,
                linear_tolerance_ns: tolerance,
                angular_tolerance_nms: tolerance,
                work_tolerance_j: tolerance,
            },
            linear_speed_tolerance_m_s: tolerance,
            angular_speed_tolerance_rad_s: tolerance,
        },
        source_substeps: 1,
        limits: ViscousNetworkLimits {
            maximum_bodies: bodies,
            maximum_connectors: connectors,
            maximum_pair_relaxations: work,
            maximum_live_bytes: LIVE_BYTES,
        },
    }
}

/// A shear film's diagonal damping, each value receipted with its formula or reason.
fn film(
    ctx: &GroupLoweringContext<'_>,
    it: &mut Scoped,
    (fluid, state): (&ThirdBodyDecl, State<'_>),
    coefficients: &mut [[f64; 6]; 6],
) -> Result<(), LoweringError> {
    let mu = it.third_body(ctx, fluid, "dynamic_viscosity", VISCOSITY)?;
    // Only a temperature is held; the other keys are plain declared values.
    for (key, q) in &fluid.query_point {
        let held = (q.dims == TEMPERATURE_DIMS)
            .then_some("held: mu does not follow the film's temperature (one-way)");
        let why = held.unwrap_or(&key.name);
        it.r.record_declared_field(&key.name, q.value, q.dims, &it.item, why);
    }
    it.declare("third_body.film_thickness", fluid.film_thickness);
    let wetted = ParameterKey::new("area", AREA_DIMS);
    let area = state.iter().find(|(k, _)| *k == wetted).map(|&(_, q)| q);
    let area = area.ok_or_else(|| it.invalid("a film needs `area`"))?;
    if area.value <= 0.0 {
        return Err(it.invalid("a film's `area` must be positive"));
    }
    it.declare("state.area", area);
    let (a, h) = (area.value, fluid.film_thickness.value);
    let (j, disk) = (a * a / (2.0 * PI), "disk: J = A^2 / (2 pi)");
    it.r.record_declared_field("polar_moment", j, POLAR_MOMENT, &it.item, disk);
    let shear = (mu * a / h, "mu A / h");
    let twist = (mu * j / h, "mu J / h");
    let none = (0.0, "absent: squeeze is not linear");
    let axes = [shear, shear, none, none, none, twist];
    for (axis, ((name, dims), (c, why))) in DAMPING.into_iter().zip(axes).enumerate() {
        it.r.record_declared_field(name, c, dims, &it.item, why);
        coefficients[axis][axis] = c;
    }
    Ok(())
}

/// Mass and principal inertia about the centroid from the part's primitive
/// geometry and the density on its card.
fn mass_properties(
    ctx: &GroupLoweringContext<'_>,
    it: &mut Scoped,
    part: &DraftPartDecl,
) -> Result<MassProperties, LoweringError> {
    let shape = geometry_volume_m3(&part.geometry).zip(principal_inertia_per_mass(&part.geometry));
    let Some(((volume, volume_formula), (per_mass, formula))) = shape else {
        return Err(it.invalid("a network body is a box, cylinder, sphere or tube"));
    };
    let material = part
        .materials
        .first()
        .ok_or_else(|| it.invalid("a body needs a material card"))?;
    let mass = it.material(ctx, material, "density", DENSITY)? * volume;
    it.derive(part, "volume", volume, VOLUME_DIMS, volume_formula);
    it.derive(part, "mass", mass, MASS_DIMS, "mass = density * volume");
    let inertia = per_mass.map(|i| mass * i);
    for (name, value) in ["inertia_x", "inertia_y", "inertia_z"]
        .into_iter()
        .zip(inertia)
    {
        it.derive(part, name, value, INERTIA_DIMS, formula);
    }
    let principal = Vec3::new(inertia[0], inertia[1], inertia[2]);
    let mass = MassProperties::new(mass, Vec3::ZERO, principal);
    mass.map_err(|err| it.invalid(format!("mass properties: {err:?}")))
}

/// Principal moments of inertia per unit mass about the centroid, along the
/// primitive's own axes (a cylinder's or tube's axis is z), and their formula.
#[must_use]
pub fn principal_inertia_per_mass(
    geometry: &DraftGeometrySource,
) -> Option<([f64; 3], &'static str)> {
    let DraftGeometrySource::Primitive(shape) = geometry else {
        return None;
    };
    Some(match *shape {
        PrimitiveShape::Box { lx, ly, lz } => {
            let [x, y, z] = [lx.value, ly.value, lz.value].map(|side| side * side / 12.0);
            (
                [y + z, x + z, x + y],
                "box: I_x = m (ly^2 + lz^2) / 12, and so on",
            )
        }
        PrimitiveShape::Cylinder { radius, height } => {
            let (r2, h2) = (radius.value * radius.value, height.value * height.value);
            let side = (3.0 * r2 + h2) / 12.0;
            (
                [side, side, r2 / 2.0],
                "cylinder: I_x = m (3 r^2 + h^2) / 12, I_z = m r^2 / 2",
            )
        }
        PrimitiveShape::Sphere { radius } => (
            [0.4 * radius.value * radius.value; 3],
            "sphere: I = 2 m r^2 / 5",
        ),
        PrimitiveShape::Tube {
            inner_radius,
            outer_radius,
            length,
        } => {
            let sum = inner_radius.value.powi(2) + outer_radius.value.powi(2);
            let side = (3.0 * sum + length.value.powi(2)) / 12.0;
            (
                [side, side, sum / 2.0],
                "tube: I_x = m (3 (a^2 + b^2) + l^2) / 12, I_z = m (a^2 + b^2) / 2",
            )
        }
    })
}

/// The declared initial state as the network's state: the pose, linear
/// momentum `m v`, and body-frame angular momentum `I omega`, with `v` and
/// `omega` declared in the assembly frame.
fn initial_state(
    it: &mut Scoped,
    part: &DraftPartDecl,
    mass: MassProperties,
) -> Result<RigidBodyState, LoweringError> {
    let initial = part.initial.as_ref();
    let initial = initial.ok_or_else(|| it.invalid("a network body declares its initial state"))?;
    let mut vector = |field: &str, values: &[QtyAny; 3]| {
        for (axis, value) in values.iter().enumerate() {
            it.declare(&format!("initial.{field}[{axis}]"), *value);
        }
        let [x, y, z] = values.map(|value| value.value);
        Vec3::new(x, y, z)
    };
    let velocity = vector("linear_velocity", &initial.linear_velocity);
    let omega = vector("angular_velocity", &initial.angular_velocity);
    let pose = it.pose(initial)?;
    let omega_body = pose.orientation().rotate_world_to_body(omega);
    let state = mass
        .angular_momentum_body_checked(omega_body)
        .and_then(|angular| RigidBodyState::new(pose, velocity.scale(mass.mass()), angular));
    state.map_err(|err| it.invalid(format!("initial state: {err:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_inertia_is_the_closed_form() {
        // A solid box of sides a, b, c: I_x = m (b^2 + c^2) / 12, and so on.
        let side = |value| QtyAny::new(value, Dims([1, 0, 0, 0, 0, 0]));
        let shape = PrimitiveShape::Box {
            lx: side(0.3),
            ly: side(0.5),
            lz: side(0.7),
        };
        let geometry = DraftGeometrySource::Primitive(shape);
        let (per_mass, _) = principal_inertia_per_mass(&geometry).unwrap();
        let expected = [0.74 / 12.0, 0.58 / 12.0, 0.34 / 12.0];
        for (got, want) in per_mass.into_iter().zip(expected) {
            assert!((got - want).abs() <= 1e-15, "{got} vs {want}");
        }
    }
}
