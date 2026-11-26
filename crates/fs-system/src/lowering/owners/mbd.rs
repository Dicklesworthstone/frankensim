//! The rigid bodies joined by `viscous-bushing` interfaces, lowered to one
//! `fs_mbd` `ViscousNetwork`.
//!
//! A body's mass and principal inertia come from its primitive geometry and
//! the density on its card, about its centroid, which is its body origin;
//! its pose and motion are its declared initial state. A bushing's damping
//! is the diagonal SI form on its card (`damping_x`, `damping_y`,
//! `damping_z` in N s/m; `damping_rx`, `damping_ry`, `damping_rz` in
//! N m s/rad), acting at both bodies' origins. The network checks impulse,
//! work and speed to the bodies' shared absolute tolerance, relaxes at most
//! their `max_work_units` pairs, and takes one source substep per step. Its
//! body forces are bound to zero; a bushing's heat goes where its loss
//! routes send it, all of it to `unallocated` when it has none.

use fs_mbd::component::network::ViscousNetwork;
use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::{
    ViscousNetworkLimits, ViscousNetworkLink,
};
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::ViscousNetworkMotionControls;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, DraftPartDecl, DraftPartKind, PrimitiveShape,
};

use super::{
    LOSS_ROUTE_FAMILY, Scoped, VISCOUS_BUSHING_FAMILY, family_entry, of_family, owning_part,
    route_bushing, slot, unallocated,
};
use crate::RouteTemperature;
use crate::lowering::{
    CardStore, DAMPING_DIMS, GroupLoweringContext, GroupSelection, INERTIA_DIMS,
    LoweredGroupOutput, LoweringError, LoweringGraph, MASS_DIMS, TORSIONAL_DAMPING_DIMS,
    VOLUME_DIMS, geometry_volume_m3,
};

/// Group name of the viscous network.
pub const GROUP: &str = "owners.mbd.viscous-network";
const NETWORK: &str = "viscous-network";
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
/// A bushing card's diagonal coefficients, in the matrix's order.
const DAMPING: [(&str, Dims); 6] = [
    ("damping_x", DAMPING_DIMS),
    ("damping_y", DAMPING_DIMS),
    ("damping_z", DAMPING_DIMS),
    ("damping_rx", TORSIONAL_DAMPING_DIMS),
    ("damping_ry", TORSIONAL_DAMPING_DIMS),
    ("damping_rz", TORSIONAL_DAMPING_DIMS),
];
/// Working memory one network step may hold: a resource cap, not physics.
const LIVE_BYTES: usize = 64 << 20;

/// The bushings, and every part one joins, as one selection.
pub(crate) fn select(
    draft: &DraftAssembly,
    cards: &CardStore,
) -> Result<Vec<GroupSelection>, LoweringError> {
    let input = draft.input();
    let mut selection = of_family(draft, cards, VISCOUS_BUSHING_FAMILY);
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

/// One `ViscousNetwork` of the selected bodies and bushings.
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
        out.parts.push(it.part(NETWORK));
    }
    let mut links = Vec::new();
    for &index in &ctx.selection.interfaces {
        let iface = &input.interfaces[index];
        let mut it = Scoped::new("interface", &iface.name);
        let card = family_entry(ctx.cards, iface, VISCOUS_BUSHING_FAMILY);
        let card = card.ok_or_else(|| it.invalid("not a viscous-bushing interface"))?;
        let mut coefficients = [[0.0; 6]; 6];
        for (axis, (name, dims)) in DAMPING.into_iter().enumerate() {
            coefficients[axis][axis] = it.interface(ctx, card, name, dims)?;
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
        out.interfaces.push(it.interface_output(Some(NETWORK)));
    }
    let [x, y, z] = input.environment.base.gravity.map(|g| g.value);
    let gravity = Gravity::new(Vec3::new(x, y, z));
    let gravity =
        gravity.map_err(|err| Scoped::new("environment", "gravity").invalid(format!("{err:?}")))?;
    let controls = controls(tolerance, properties.len(), links.len(), work);
    let bodies = properties.len();
    let network = ViscousNetwork::new(properties, links, gravity, controls);
    let network =
        network.map_err(|err| Scoped::new("network", NETWORK).invalid(format!("{err:?}")))?;
    let node = graph.add(NETWORK, Some(first.entity), network, states)?;
    for &index in &parts[1..] {
        graph
            .entity_nodes
            .entry(input.parts[index].entity)
            .or_default()
            .push(node);
    }
    for body in 0..bodies {
        let force = graph.port(node, &format!("force{body}"))?;
        graph.bind_zero_boundary(force);
    }
    let ambient_k = input.environment.base.ambient_temperature.value;
    for (link, &index) in ctx.selection.interfaces.iter().enumerate() {
        let bushing = input.interfaces[index].entity;
        let heat = graph.port(node, &format!("heat{link}"))?;
        graph.expose_port(bushing, "heat", heat);
        let routed = input.interfaces.iter().any(|route| {
            family_entry(ctx.cards, route, LOSS_ROUTE_FAMILY).is_some()
                && route_bushing(input, ctx.cards, route).is_some_and(|b| b.entity == bushing)
        });
        if !routed {
            let sink = Some(unallocated(ambient_k)?);
            graph.add_heat_route(heat, Vec::new(), sink, RouteTemperature::SinkEquivalent);
        }
    }
    Ok(out)
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
