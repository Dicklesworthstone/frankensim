//! `conduction_body`: fs-conduction's backward-Euler finite-element body,
//! [`ConductionBody`], as a scene element, and `trace_contact`: the owner's
//! contact between two of its traces as the kernel's [`ConductanceMatrix`].
//!
//! The scene authors a box `[0, Lx] x [0, Ly] x [0, Lz]` meshed as an
//! `nx x ny x nz` Kuhn grid (fs-conduction's `fixtures::box_grid`, six
//! tetrahedra a cell), one declared conductivity and volumetric heat
//! capacity, the start temperature, and the box faces that are heat ports.
//! The faces not named are insulated. A face under `ports` receives the heat
//! rate into the body and returns its area-mean end temperature; a face under
//! `traces` receives one heat load for every vertex of the face and returns
//! the end temperature at each. The owner's own backward Euler steps the
//! field. Nothing here conducts, stores or integrates heat.
//!
//! A `trace_contact` names two traces of two bodies, where the second body's
//! origin lies in the first's frame, and one declared area-specific
//! resistance with its source. The owner computes the interface points
//! (`ConductionBody::matching_contact`, or `mapped_contact` when the scene
//! gives the refinement's budgets under `mapped`) and the kernel's element
//! carries the heat. The contact is computed from both bodies, so the bodies
//! a contact names are built before the kernel takes any element
//! ([`prepare`]); every other body is built where the scene lists it.
//!
//! The kernel keeps component states as `'static` values and the body
//! borrows its mesh, boundary and material, so they are leaked: they live
//! as long as the process that runs the scene.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::{ConductionBody, PatchPort};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::NonmatchingOptions;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RESISTANCE_PROPERTY, ConductionMesh,
    ConductionProblem, ConductivityModel, InterfaceResistance, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_ports::component::PortSpec;
use fs_system::elements::{ConductanceMatrix, InterfacePoint};
use fs_system::{NodeId, System, SystemBuilder};

use super::{J, Refusal, add, fields, number, number_or, positive};

/// The kernel holds the nodal temperature field [K] in vertex order.
pub(super) type State = Vec<f64>;

/// The box faces a port may name: the axis and whether the face is at the
/// box's far end.
const FACES: [(&str, usize, bool); 6] = [
    ("x_min", 0, false),
    ("x_max", 0, true),
    ("y_min", 1, false),
    ("y_max", 1, true),
    ("z_min", 2, false),
    ("z_max", 2, true),
];
/// Cells along one axis, and in the whole box: the body is stepped by a
/// sparse Krylov solve on every kernel step, so a scene stays small.
const MAX_CELLS_PER_AXIS: f64 = 64.0;
const MAX_CELLS: usize = 4096;
/// Bounds of the owner's Krylov budget a scene may give.
const MAX_ITERATIONS: f64 = 100_000.0;
const MAX_RESTART: f64 = 1000.0;
/// The largest budget of the common refinement a scene may give: its pair
/// tests and its retained triangles.
const MAX_REFINEMENT_BUDGET: f64 = 1.0e9;

/// A whole number from 1 to `max` at `key`, or `default`.
fn whole(value: &J, key: &str, default: usize, max: f64, what: &str) -> Result<usize, String> {
    let Some(given) = value.get(key) else {
        return Ok(default);
    };
    let number = given
        .as_f64()
        .filter(|n| (1.0..=max).contains(n) && n.fract() == 0.0)
        .ok_or_else(|| format!("{what}: {key} is a whole number from 1 to {max}"))?;
    // A whole number no larger than `max`, checked above.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    Ok(number as usize)
}

/// Three numbers at `key`.
fn three(definition: &J, key: &str, what: &str) -> Result<[f64; 3], String> {
    let values: Vec<f64> = definition
        .get(key)
        .and_then(J::as_array)
        .filter(|values| values.len() == 3)
        .ok_or_else(|| format!("{what}: {key} is an array of three numbers"))?
        .iter()
        .map(J::as_f64)
        .collect::<Option<_>>()
        .ok_or_else(|| format!("{what}: {key} is an array of three numbers"))?;
    Ok([values[0], values[1], values[2]])
}

/// Three positive finite numbers at `key`.
fn triple(definition: &J, key: &str, what: &str) -> Result<[f64; 3], String> {
    let values = three(definition, key, what)?;
    if values
        .iter()
        .any(|value| !(value.is_finite() && *value > 0.0))
    {
        return Err(format!("{what}: {key} must be positive and finite"));
    }
    Ok(values)
}

/// Three finite numbers at `key`, of either sign.
fn vector(definition: &J, key: &str, what: &str) -> Result<[f64; 3], String> {
    let values = three(definition, key, what)?;
    if values.iter().any(|value| !value.is_finite()) {
        return Err(format!("{what}: {key} must be finite"));
    }
    Ok(values)
}

/// The cell counts along x, y and z.
fn cells(definition: &J, what: &str) -> Result<[usize; 3], String> {
    let counts = triple(definition, "cells", what)?;
    let mut cells = [0; 3];
    for (cell, count) in cells.iter_mut().zip(counts) {
        if count > MAX_CELLS_PER_AXIS || count.fract() != 0.0 {
            return Err(format!(
                "{what}: cells are whole numbers from 1 to {MAX_CELLS_PER_AXIS}"
            ));
        }
        // A whole number from 1 to 64, checked above.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = count as usize;
        *cell = count;
    }
    if cells.iter().product::<usize>() > MAX_CELLS {
        return Err(format!("{what}: at most {MAX_CELLS} cells"));
    }
    Ok(cells)
}

/// The owner's step configuration: its energy gate (required: it is the
/// scene's scale) and its Krylov budget, by default the owner's own.
fn config(definition: &J, what: &str) -> Result<StepConfig, String> {
    let inner = format!("{what}.solve");
    let solve = definition
        .get("solve")
        .ok_or_else(|| format!("{what}: missing solve (energy_tolerance_j)"))?;
    fields(
        solve,
        &[
            "energy_tolerance_j",
            "tolerance",
            "max_iterations",
            "restart",
        ],
        &inner,
    )?;
    let default = LinearConfig::default();
    let tolerance = number_or(solve, "tolerance", default.tolerance, &inner)?;
    if !(tolerance > 0.0 && tolerance < 1.0) {
        return Err(format!("{inner}: tolerance must lie in (0, 1)"));
    }
    Ok(StepConfig {
        linear: LinearConfig {
            tolerance,
            max_iterations: whole(
                solve,
                "max_iterations",
                default.max_iterations,
                MAX_ITERATIONS,
                &inner,
            )?,
            restart: whole(solve, "restart", default.restart, MAX_RESTART, &inner)?,
        },
        energy_tolerance_j: positive(solve, "energy_tolerance_j", &inner)?,
    })
}

/// The faces named under `key`, each once: name, axis and far end. A body
/// with no `key` has none.
fn faces(
    definition: &J,
    key: &str,
    what: &str,
) -> Result<Vec<(&'static str, usize, bool)>, String> {
    let Some(given) = definition.get(key) else {
        return Ok(Vec::new());
    };
    let names = given
        .as_array()
        .filter(|names| !names.is_empty())
        .ok_or_else(|| format!("{what}: {key} is a nonempty array of box faces"))?;
    let mut faces = Vec::with_capacity(names.len());
    for name in names {
        let face = FACES
            .iter()
            .find(|(face, _, _)| Some(*face) == name.as_str())
            .ok_or_else(|| {
                format!("{what}: {key} name box faces: x_min, x_max, y_min, y_max, z_min, z_max")
            })?;
        if faces.contains(face) {
            return Err(format!("{what}: port {} is named twice", face.0));
        }
        faces.push(*face);
    }
    Ok(faces)
}

/// The faces that are ports, the mean ports first and then the traces, with
/// whether each is a trace. A face is one or the other.
fn ports(definition: &J, what: &str) -> Result<Vec<(&'static str, usize, bool, bool)>, String> {
    let means = faces(definition, "ports", what)?;
    let traces = faces(definition, "traces", what)?;
    if means.is_empty() && traces.is_empty() {
        return Err(format!(
            "{what}: ports is a nonempty array of box faces (or traces is)"
        ));
    }
    if let Some(both) = traces.iter().find(|face| means.contains(face)) {
        return Err(format!(
            "{what}: face {} is under ports and under traces; it is one port",
            both.0
        ));
    }
    let mean = means
        .into_iter()
        .map(|(face, axis, far)| (face, axis, far, false));
    let trace = traces
        .into_iter()
        .map(|(face, axis, far)| (face, axis, far, true));
    Ok(mean.chain(trace).collect())
}

/// A context of its own for an owner's preparation, which keeps no part of
/// it.
fn with_cx<R>(run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 0,
            kernel_id: 0x5953,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        run(&cx)
    })
}

/// A body before the kernel takes it.
struct Body {
    component: ConductionBody<'static>,
    start: State,
    temperature_k: f64,
    /// The body's ports in port order: the face and whether it is a trace.
    faces: Vec<(&'static str, bool)>,
}

/// The owner's body on its box mesh, at a uniform start temperature.
fn build(definition: &J, what: &str) -> Result<Body, Refusal> {
    fields(
        definition,
        &[
            "kind",
            "box_m",
            "cells",
            "conductivity_w_m_k",
            "volumetric_heat_capacity_j_m3_k",
            "temperature_k",
            "ports",
            "traces",
            "solve",
        ],
        what,
    )?;
    let extent = triple(definition, "box_m", what)?;
    let cells = cells(definition, what)?;
    let ports = ports(definition, what)?;
    let config = config(definition, what)?;
    let temperature_k = positive(definition, "temperature_k", what)?;
    let conductivity =
        ConductivityModel::isotropic_declared(number(definition, "conductivity_w_m_k", what)?)
            .map_err(|error| format!("{what}: {error}"))?;
    let capacity = VolumetricHeatCapacity::declared(number(
        definition,
        "volumetric_heat_capacity_j_m3_k",
        what,
    )?)
    .map_err(|error| format!("{what}: {error}"))?;
    let (complex, positions) = box_grid(cells, extent);
    let mesh: &'static ConductionMesh = Box::leak(Box::new(
        ConductionMesh::new(complex, positions).map_err(|error| format!("{what}: {error}"))?,
    ));
    let mut regions = ThermalBoundaryBuilder::new(mesh);
    for &(port, axis, far, _) in &ports {
        let at = if far { extent[axis] } else { 0.0 };
        regions = regions
            .region(
                port,
                |face| on_box_face(face.centroid[axis], at),
                ThermalBc::adiabatic(),
            )
            .map_err(|error| format!("{what}: {error}"))?;
    }
    let boundary = Box::leak(Box::new(
        regions
            .adiabatic_remainder()
            .finish()
            .map_err(|error| format!("{what}: {error}"))?,
    ));
    // The owner prepares its capacitance once.
    let solver = with_cx(|cx| BackwardEuler::uniform(cx, mesh, capacity))
        .map_err(|error| format!("{what}: {error}"))?;
    let problem = ConductionProblem {
        mesh,
        boundary,
        material: Box::leak(Box::new(conductivity)),
        element_materials: None,
        source: Box::leak(Box::new(ScalarField::Uniform(0.0))),
    };
    // A face under `ports` is a mean port: one heat rate in, its area-mean
    // temperature out. A face under `traces` carries a value a vertex.
    let patches: Vec<PatchPort<'_>> = ports
        .iter()
        .map(|&(port, _, _, trace)| {
            if trace {
                PatchPort::Trace(port)
            } else {
                PatchPort::Mean(port)
            }
        })
        .collect();
    let component = ConductionBody::new(solver, problem, &patches, config)
        .map_err(|error| format!("{what}: {error}"))?;
    Ok(Body {
        component,
        start: vec![temperature_k; mesh.vertex_count()],
        temperature_k,
        faces: ports
            .iter()
            .map(|&(face, _, _, trace)| (face, trace))
            .collect(),
    })
}

/// What a scene declared for a contact and what the owner made of it: kept
/// for the wiring check and the report.
#[derive(Debug, Clone)]
pub(super) struct Declared {
    /// The two traces, `element.face`, that the contact was computed for:
    /// its port `a` meets the first and its port `b` the second.
    pub(super) between: [String; 2],
    pub(super) resistance_m2_k_w: f64,
    /// Where the scene says the resistance comes from.
    pub(super) source: String,
    /// The sum of the interface points' conductances [W/K].
    conductance_w_k: f64,
    points: usize,
}

/// The bodies that a scene's contacts name and those contacts, built before
/// the kernel takes any element.
#[derive(Default)]
pub(super) struct Prepared {
    bodies: Vec<(String, Body)>,
    contacts: Vec<(String, ConductanceMatrix, Declared)>,
}

/// One end of a contact: the element and the face of `text`, the body among
/// `bodies`, built now if no earlier contact named it, and the face's port.
fn end(
    definitions: &[(String, J)],
    bodies: &mut Vec<(String, Body)>,
    text: &str,
    what: &str,
) -> Result<(usize, usize), Refusal> {
    let (name, face) = text
        .rsplit_once('.')
        .ok_or_else(|| format!("{what}: between names two traces as element.face, not {text:?}"))?;
    let definition = definitions
        .iter()
        .find(|(known, _)| known == name)
        .map(|(_, definition)| definition)
        .ok_or_else(|| format!("{what}: between names {text}, and no element is named {name:?}"))?;
    if definition.str_field("kind") != Some("conduction_body") {
        return Err(
            format!("{what}: between names {text}, and {name} is not a conduction_body").into(),
        );
    }
    let body = if let Some(body) = bodies.iter().position(|(known, _)| known == name) {
        body
    } else {
        let built = build(definition, &format!("elements.{name}"))?;
        bodies.push((name.to_owned(), built));
        bodies.len() - 1
    };
    let port = bodies[body]
        .1
        .faces
        .iter()
        .position(|&(known, trace)| known == face && trace)
        .ok_or_else(|| {
            format!("{what}: between names {text}, which is not one of {name}'s traces")
        })?;
    Ok((body, port))
}

/// The declared resistance as the owner's interface resistance, through a
/// card that says the value is the caller's declaration.
fn resistance(name: &str, declared: &Declared) -> Result<InterfaceResistance, String> {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RESISTANCE_PROPERTY, RESISTANCE_DIMS),
            value: PropertyValue::Scalar {
                value: declared.resistance_m2_k_w,
                dims: RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!(
                    "caller-declared system-scene contact {name}: {}",
                    declared.source
                ),
                license: "unspecified".into(),
                artifact: None,
            },
        })
        .map_err(|error| error.to_string())?;
    let side = |trace: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: trace.to_owned(),
            phase: "caller-unspecified".into(),
            process: "caller-declared contact side".into(),
            revision: 0,
        },
        texture_frame: "caller-unspecified".into(),
    };
    let card = InterfaceSystemCard::assemble(
        side(&declared.between[0]),
        side(&declared.between[1]),
        SystemContext {
            medium: "caller-unspecified".into(),
            third_body: None,
            environment: "caller-unspecified".into(),
            history: "caller-unspecified".into(),
        },
        claims,
        Vec::new(),
    )
    .map_err(|error| error.to_string())?;
    InterfaceResistance::from_card(
        name,
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .map_err(|error| error.to_string())
}

/// The budgets of the owner's common refinement, all four as the scene
/// gives them: none has a default here.
fn refinement(definition: &J, what: &str) -> Result<Option<NonmatchingOptions>, String> {
    let Some(mapped) = definition.get("mapped") else {
        return Ok(None);
    };
    let inner = format!("{what}.mapped");
    fields(
        mapped,
        &[
            "plane_tolerance_m",
            "coverage_relative_tolerance",
            "max_pair_tests",
            "max_overlap_triangles",
        ],
        &inner,
    )?;
    let budget = |key: &str| {
        number(mapped, key, &inner)?;
        whole(mapped, key, 0, MAX_REFINEMENT_BUDGET, &inner)
    };
    Ok(Some(NonmatchingOptions {
        plane_tolerance_m: number(mapped, "plane_tolerance_m", &inner)?,
        coverage_relative_tolerance: number(mapped, "coverage_relative_tolerance", &inner)?,
        max_pair_tests: budget("max_pair_tests")?,
        max_overlap_triangles: budget("max_overlap_triangles")?,
    }))
}

/// A `trace_contact` element: the owner's interface points between the two
/// traces as the kernel's element.
fn contact(
    definitions: &[(String, J)],
    bodies: &mut Vec<(String, Body)>,
    name: &str,
    definition: &J,
) -> Result<(ConductanceMatrix, Declared), Refusal> {
    let what = format!("elements.{name}");
    fields(
        definition,
        &[
            "kind",
            "between",
            "offset_m",
            "resistance_m2_k_w",
            "source",
            "mapped",
        ],
        &what,
    )?;
    let between = match definition.get("between").and_then(J::as_array) {
        Some([a, b]) => a.as_str().zip(b.as_str()),
        _ => None,
    }
    .map(|(a, b)| [a, b])
    .ok_or_else(|| {
        format!("{what}: between is [\"element.face\", \"element.face\"], two traces")
    })?;
    let offset_m = vector(definition, "offset_m", &what)?;
    let source = definition
        .str_field("source")
        .filter(|source| !source.trim().is_empty())
        .ok_or_else(|| format!("{what}: source says where resistance_m2_k_w comes from"))?;
    let options = refinement(definition, &what)?;
    let (body_a, port_a) = end(definitions, bodies, between[0], &what)?;
    let (body_b, port_b) = end(definitions, bodies, between[1], &what)?;
    if body_a == body_b {
        return Err(format!(
            "{what}: between names two traces of one body; a contact joins two bodies"
        )
        .into());
    }
    let mut declared = Declared {
        between: [between[0].to_owned(), between[1].to_owned()],
        resistance_m2_k_w: positive(definition, "resistance_m2_k_w", &what)?,
        source: source.to_owned(),
        conductance_w_k: 0.0,
        points: 0,
    };
    let resistance = resistance(name, &declared).map_err(|error| format!("{what}: {error}"))?;
    let (a, b) = (&bodies[body_a].1.component, &bodies[body_b].1.component);
    let joint = match options {
        // The scene names its own remedy: the owner knows nothing of scenes.
        None => a
            .matching_contact(port_a, b, port_b, offset_m, &resistance)
            .map_err(|error| {
                format!(
                    "{what}: {error} (without mapped the two faces must have the same points \
                     and triangles; mapped gives the refinement's budgets for faces that do not)"
                )
            }),
        Some(options) => {
            with_cx(|cx| a.mapped_contact(cx, port_a, b, port_b, offset_m, &resistance, options))
                .map_err(|error| format!("{what}: {error}"))
        }
    }?;
    declared.points = joint.points.len();
    declared.conductance_w_k = joint.points.iter().map(|point| point.conductance_w_k).sum();
    let points = joint
        .points
        .iter()
        .map(|point| InterfacePoint {
            a: point.a,
            b: point.b,
            conductance_w_k: point.conductance_w_k,
        })
        .collect();
    let wall = ConductanceMatrix::new(joint.width_a, joint.width_b, points)
        .map_err(|error| format!("{what}: {error}"))?;
    Ok((wall, declared))
}

/// Every `trace_contact` of the scene with the bodies it names. A scene
/// with no contact prepares nothing.
pub(super) fn prepare(definitions: &[(String, J)]) -> Result<Prepared, Refusal> {
    let mut prepared = Prepared::default();
    for (name, definition) in definitions {
        if definition.str_field("kind") == Some("trace_contact") {
            let (wall, declared) = contact(definitions, &mut prepared.bodies, name, definition)?;
            prepared.contacts.push((name.clone(), wall, declared));
        }
    }
    Ok(prepared)
}

/// A `conduction_body` element: the owner's body on its box mesh, at a
/// uniform start temperature, which it also returns.
pub(super) fn body(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    prepared: &mut Prepared,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>, f64), Refusal> {
    let body = match prepared.bodies.iter().position(|(known, _)| known == name) {
        Some(built) => prepared.bodies.swap_remove(built).1,
        None => build(definition, what)?,
    };
    let (node, specs) = add(builder, name, body.component, body.start)?;
    Ok((node, specs, body.temperature_k))
}

/// A `trace_contact` element, prepared with its bodies.
pub(super) fn wall(
    builder: &mut SystemBuilder,
    name: &str,
    prepared: &mut Prepared,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>, Declared), Refusal> {
    let built = prepared
        .contacts
        .iter()
        .position(|(known, _, _)| known == name)
        .ok_or_else(|| format!("{what}: the contact was not prepared"))?;
    let (_, wall, declared) = prepared.contacts.swap_remove(built);
    let (node, specs) = add(builder, name, wall, ())?;
    Ok((node, specs, declared))
}

/// Each mean port's area-mean temperature and each trace's least and
/// greatest temperature at the end of the last step, and the field's least
/// and greatest nodal temperatures.
pub(super) fn values(system: &System, node: NodeId, ports: &[PortSpec]) -> Vec<(String, f64)> {
    let mut values = Vec::with_capacity(ports.len() + 2);
    for port in ports {
        let id = port.id().as_str();
        let returned = system
            .port(node, id)
            .filter(|_| system.steps() > 0)
            .map(|port| system.last_output(port));
        if port.width() == 1 {
            let mean = returned.and_then(|values| values.first().copied());
            values.push((format!("{id}.temperature_k"), mean.unwrap_or(f64::NAN)));
        } else {
            let (low, high) = returned.map_or((f64::NAN, f64::NAN), span);
            values.push((format!("{id}.min_temperature_k"), low));
            values.push((format!("{id}.max_temperature_k"), high));
        }
    }
    let field = system.state::<State>(node);
    let (low, high) = field.map_or((f64::NAN, f64::NAN), |field| span(field));
    values.push(("min_temperature_k".to_owned(), low));
    values.push(("max_temperature_k".to_owned(), high));
    values
}

/// The least and the greatest of `values`.
fn span(values: &[f64]) -> (f64, f64) {
    values
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), &t| {
            (low.min(t), high.max(t))
        })
}

/// A contact's heat rate from its side `a` to its side `b` over the last
/// step, the sum of its interface points' conductances, the resistance the
/// scene declared and the number of interface points.
pub(super) fn wall_values(
    system: &System,
    node: NodeId,
    declared: &Declared,
) -> Vec<(String, f64)> {
    let rate = system
        .port(node, "a")
        .filter(|_| system.steps() > 0)
        .map_or(f64::NAN, |port| system.last_output(port).iter().sum());
    // A count of interface points, far below 2^53.
    #[allow(clippy::cast_precision_loss)]
    let points = declared.points as f64;
    vec![
        ("heat_rate_w".to_owned(), rate),
        ("conductance_w_k".to_owned(), declared.conductance_w_k),
        ("resistance_m2_k_w".to_owned(), declared.resistance_m2_k_w),
        ("interface_points".to_owned(), points),
    ]
}
