//! `conduction_body`: fs-conduction's backward-Euler finite-element body,
//! [`ConductionBody`], as a scene element, and `trace_contact`: the owner's
//! contact between two of its traces as the kernel's [`ConductanceMatrix`].
//!
//! The scene authors a box `[0, Lx] x [0, Ly] x [0, Lz]` meshed as an
//! `nx x ny x nz` Kuhn grid (fs-conduction's `fixtures::box_grid`, six
//! tetrahedra a cell), one conductivity, one declared volumetric heat
//! capacity, the start temperature, and the box faces that are heat ports.
//! The faces not named are insulated. A face under `ports` receives the heat
//! rate into the body and returns its area-mean end temperature; a face under
//! `traces` receives one heat load for every vertex of the face and returns
//! the end temperature at each. The owner's own backward Euler steps the
//! field. Nothing here conducts, stores or integrates heat.
//!
//! The conductivity is a number the scene declares, or the k(T) curve of a
//! material from the scene's `materials` section ([`conductivity`]): the
//! curve the material's card gives over its declared temperature span. A
//! body whose conductivity varies over that span steps on the owner's Newton
//! path (`ConductionBody::temperature_dependent`) under the scene's `newton`
//! policy, by default the owner's. Its heat capacity stays the scene's
//! constant, as the owner it wraps requires.
//!
//! A body that names a material and declares no heat capacity takes its
//! heat capacity and density from the material's card as well
//! ([`card_material`]): fs-conduction prepares its paired material from the
//! card (`MaterialEnthalpyTransport::from_isotropic_card_exact`), the chart
//! of the card's Cp(T) and density, its k(T), and the density read once at
//! the body's `reference_temperature_k`, which fixes the nodal masses. The
//! body is then the enthalpy owner behind the same ports
//! ([`EnthalpyBody`]): its state is the nodal specific enthalpy, and the
//! report reads temperatures through the card's chart as the owner does.
//!
//! A `trace_contact` names two traces of two bodies, where the second body's
//! origin lies in the first's frame, and one area-specific resistance: a
//! number the scene declares with its source, or an interface card of a
//! compiled pack read once at a declared state ([`card`]). The owner computes
//! the interface points (`ConductionBody::matching_contact`, or
//! `mapped_contact` when the scene gives the refinement's budgets under
//! `mapped`) and the kernel's element carries the heat. A card-backed
//! contact that names the card coordinates its traces' temperatures take
//! (`follows`) is the owner's [`TraceResistanceLaw`] on the element's hook:
//! the card is asked at every interface point of every step, and at every
//! pair of facing trace points, and a state it has no answer for refuses
//! the step. The contact is
//! computed from both bodies, so the bodies a contact names are built before
//! the kernel takes any element ([`prepare`]); every other body is built
//! where the scene lists it.
//!
//! The kernel keeps component states as `'static` values and the body
//! borrows its mesh, boundary and material, so they are leaked: they live
//! as long as the process that runs the scene.

use std::fmt::Write as _;
use std::path::Path;
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_bytes;
use fs_conduction::component::{
    ConductionBody, EnthalpyBody, PatchPort, TraceContact, TraceResistanceLaw,
};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::NonmatchingOptions;
use fs_conduction::interface::temperature::TemperatureContactCoordinates;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::transient::enthalpy::material::MaterialEnthalpyTransport;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RESISTANCE_PROPERTY, ConductionMesh,
    ConductionProblem, ConductivityModel, ConductivityTable, InterfaceResistance, LinearConfig,
    ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_ports::component::PortSpec;
use fs_solver::NewtonKrylovConfig;
use fs_system::elements::{
    ConductanceLaw, ConductanceMatrix, ElementError, GuardRow, InterfacePoint,
};
use fs_system::{NodeId, System, SystemBuilder};

use crate::cards::{MAX_CARD_PACK_BYTES, MAX_CARD_PACK_SOURCE_BYTES};
use crate::material_state_input::state_point;

use super::materials::{Materials, read_pack};
use super::{Element, J, Refusal, add, fields, num, number, number_or, positive, quote};

/// The kernel holds the nodal field in vertex order: the temperature [K],
/// or the specific enthalpy [J/kg] of a body whose heat capacity is its
/// material card's.
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
/// The largest Newton update budget of one step a scene may give.
const MAX_NEWTON_UPDATES: f64 = 1000.0;
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

/// The owner's Newton policy for a conductivity that depends on
/// temperature: its update budget and its two residual gates, each by
/// default the owner's own. The owner refuses a policy it cannot run.
fn newton(definition: &J, what: &str) -> Result<NonlinearStepConfig, String> {
    let default = NonlinearStepConfig::default();
    let Some(policy) = definition.get("newton") else {
        return Ok(default);
    };
    let inner = format!("{what}.newton");
    fields(
        policy,
        &["max_iterations", "residual_rtol", "residual_atol_j"],
        &inner,
    )?;
    Ok(NonlinearStepConfig {
        max_iterations: whole(
            policy,
            "max_iterations",
            default.max_iterations,
            MAX_NEWTON_UPDATES,
            &inner,
        )?,
        residual_rtol: number_or(policy, "residual_rtol", default.residual_rtol, &inner)?,
        residual_atol_j: number_or(policy, "residual_atol_j", default.residual_atol_j, &inner)?,
        ..default
    })
}

/// A body's conductivity, the scene material it came from, and the owner's
/// Newton policy when it varies with temperature.
struct Conductivity {
    model: ConductivityModel,
    material: Option<String>,
    nonlinear: Option<NonlinearStepConfig>,
}

/// The scene's `conductivity_w_m_k`, or the k(T) curve of the scene
/// material named by `material`, whose span must hold the start
/// temperature.
fn conductivity(
    definition: &J,
    materials: &Materials,
    temperature_k: f64,
    what: &str,
) -> Result<Conductivity, String> {
    let Some(name) = definition.get("material") else {
        if definition.get("newton").is_some() {
            return Err(format!(
                "{what}: newton is the policy of a material's temperature-dependent \
                 conductivity; a body with conductivity_w_m_k takes none"
            ));
        }
        let model =
            ConductivityModel::isotropic_declared(number(definition, "conductivity_w_m_k", what)?)
                .map_err(|error| format!("{what}: {error}"))?;
        return Ok(Conductivity {
            model,
            material: None,
            nonlinear: None,
        });
    };
    if definition.get("conductivity_w_m_k").is_some() {
        return Err(format!(
            "{what}: the conductivity is conductivity_w_m_k or a material's, not both"
        ));
    }
    let name = name
        .as_str()
        .ok_or_else(|| format!("{what}.material: expected a material name"))?;
    let curve = materials
        .get(name)
        .ok_or_else(|| format!("{what}.material: unknown material {name:?}"))?
        .conductivity
        .as_ref()
        .ok_or_else(|| {
            format!(
                "{what}.material: material {name:?} has no conductivity (its declaration names \
                 no conductivity_property)"
            )
        })?;
    let [lower_k, upper_k] = curve.temperature_range_k();
    if !(lower_k..=upper_k).contains(&temperature_k) {
        return Err(format!(
            "{what}: temperature_k {temperature_k} K lies outside material {name:?}'s \
             conductivity, {lower_k} K to {upper_k} K"
        ));
    }
    let knots = curve.knots().iter().map(|&[t_k, k_w_m_k]| (t_k, k_w_m_k));
    let table = ConductivityTable::declared_curve(knots.collect())
        .map_err(|error| format!("{what}.material: {error}"))?;
    let nonlinear = if table.is_temperature_dependent() {
        Some(newton(definition, what)?)
    } else if definition.get("newton").is_some() {
        return Err(format!(
            "{what}: material {name:?}'s conductivity does not vary over its span, so the body \
             steps on the owner's linear path and takes no newton"
        ));
    } else {
        None
    };
    Ok(Conductivity {
        model: ConductivityModel::isotropic(table),
        material: Some(name.to_owned()),
        nonlinear,
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

/// The owner behind a body's ports, boxed until the kernel takes it.
enum Owner {
    /// A declared volumetric heat capacity: the backward-Euler owner, whose
    /// state is the nodal temperature.
    Capacity(Box<ConductionBody<'static>>),
    /// The material card's heat capacity and density: the enthalpy owner,
    /// whose state is the nodal specific enthalpy, and the paired material
    /// whose chart reads it as temperature.
    Enthalpy(
        Box<EnthalpyBody<'static, 'static>>,
        &'static MaterialEnthalpyTransport,
    ),
}

impl Body {
    /// The body's backward-Euler owner, which fs-conduction's contacts join;
    /// `what` names the contact refused one of a card's heat capacity.
    fn capacity(&self, what: &str) -> Result<&ConductionBody<'static>, String> {
        match &self.owner {
            Owner::Capacity(body) => Ok(body),
            Owner::Enthalpy(..) => Err(format!(
                "{what}: between names a body whose heat capacity is its material card's, and \
                 fs-conduction's contacts join bodies with volumetric_heat_capacity_j_m3_k"
            )),
        }
    }
}

/// A body before the kernel takes it.
struct Body {
    owner: Owner,
    start: State,
    temperature_k: f64,
    /// The scene material whose conductivity the body takes, if any.
    material: Option<String>,
    /// The body's ports in port order: the face and whether it is a trace.
    faces: Vec<(&'static str, bool)>,
}

/// A body's material: a conductivity with a declared heat capacity, or the
/// whole of a material card's.
enum Source {
    Declared(Conductivity, VolumetricHeatCapacity),
    Card(CardMaterial),
}

/// A body whose heat capacity, density and conductivity are a material
/// card's: the paired material, the start state on its chart, the owner's
/// step policy and the scene material's name.
struct CardMaterial {
    transport: &'static MaterialEnthalpyTransport,
    start_j_kg: f64,
    step: EnthalpyStepConfig,
    material: String,
}

/// The paired material of the scene material `material` names, prepared by
/// fs-conduction from its card with the density read at the body's
/// `reference_temperature_k`; the start temperature's specific enthalpy on
/// its chart; and the owner's Newton policy: `newton`'s update budget and
/// residual gates (as for a k(T) body, by default fs-conduction's), the
/// rest of the Newton-Krylov policy fs-solver's default, and `config`'s
/// energy gate.
fn card_material(
    definition: &J,
    materials: &Materials,
    temperature_k: f64,
    config: StepConfig,
    what: &str,
) -> Result<CardMaterial, String> {
    if definition.get("conductivity_w_m_k").is_some() {
        return Err(format!(
            "{what}: the conductivity is conductivity_w_m_k or a material's, not both"
        ));
    }
    let name = definition
        .str_field("material")
        .ok_or_else(|| format!("{what}.material: expected a material name"))?;
    let material = materials
        .get(name)
        .ok_or_else(|| format!("{what}.material: unknown material {name:?}"))?;
    if definition.get("reference_temperature_k").is_none() {
        return Err(format!(
            "{what}: a body that declares no volumetric_heat_capacity_j_m3_k takes its heat \
             capacity and density from material {name:?}'s card, and needs \
             reference_temperature_k, the temperature at which box_m holds its mass"
        ));
    }
    let reference_k = positive(definition, "reference_temperature_k", what)?;
    let transport: &'static MaterialEnthalpyTransport =
        Box::leak(Box::new(material.enthalpy_transport(reference_k).map_err(
            |error| format!("{what}.material: material {name:?}: {error}"),
        )?));
    let start_j_kg = transport
        .caloric()
        .state_at_temperature(temperature_k)
        .map_err(|error| format!("{what}: temperature_k {temperature_k} K: {error}"))?
        .specific_enthalpy_j_kg();
    let policy = newton(definition, what)?;
    let step = EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: policy.residual_atol_j,
            relative_tolerance: policy.residual_rtol,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: policy.max_iterations,
        energy_tolerance_j: config.energy_tolerance_j,
    };
    Ok(CardMaterial {
        transport,
        start_j_kg,
        step,
        material: name.to_owned(),
    })
}

/// A body's material: a material card's whole when the body names a material
/// and declares no heat capacity, else a conductivity with the declared heat
/// capacity.
fn source(
    definition: &J,
    materials: &Materials,
    temperature_k: f64,
    config: StepConfig,
    what: &str,
) -> Result<Source, String> {
    if definition.get("material").is_some()
        && definition.get("volumetric_heat_capacity_j_m3_k").is_none()
    {
        return card_material(definition, materials, temperature_k, config, what).map(Source::Card);
    }
    if definition.get("reference_temperature_k").is_some() {
        return Err(format!(
            "{what}: reference_temperature_k fixes the mass of a body whose heat capacity is its \
             material card's; a body with volumetric_heat_capacity_j_m3_k takes none"
        ));
    }
    let conductivity = conductivity(definition, materials, temperature_k, what)?;
    let capacity = VolumetricHeatCapacity::declared(number(
        definition,
        "volumetric_heat_capacity_j_m3_k",
        what,
    )?)
    .map_err(|error| format!("{what}: {error}"))?;
    Ok(Source::Declared(conductivity, capacity))
}

/// The owner `source` makes on `mesh` and `boundary` behind the ports
/// `patches`, its start state at `temperature_k`, and the scene material it
/// names.
fn owner(
    source: Source,
    mesh: &'static ConductionMesh,
    boundary: &'static ThermalBoundary,
    patches: &[PatchPort<'_>],
    temperature_k: f64,
    config: StepConfig,
    what: &str,
) -> Result<(Owner, State, Option<String>), String> {
    let source_field = Box::leak(Box::new(ScalarField::Uniform(0.0)));
    match source {
        Source::Declared(
            Conductivity {
                model,
                material,
                nonlinear,
            },
            capacity,
        ) => {
            // The owner prepares its capacitance once.
            let solver = with_cx(|cx| BackwardEuler::uniform(cx, mesh, capacity))
                .map_err(|error| format!("{what}: {error}"))?;
            let problem = ConductionProblem {
                mesh,
                boundary,
                material: Box::leak(Box::new(model)),
                element_materials: None,
                source: source_field,
            };
            let component = match nonlinear {
                None => ConductionBody::new(solver, problem, patches, config),
                Some(policy) => {
                    ConductionBody::temperature_dependent(solver, problem, patches, config, policy)
                }
            }
            .map_err(|error| format!("{what}: {error}"))?;
            let start = vec![temperature_k; mesh.vertex_count()];
            Ok((Owner::Capacity(Box::new(component)), start, material))
        }
        Source::Card(card) => {
            // The owner prepares its reference masses and storage once.
            let budget = EnthalpyBudget {
                max_vertices: mesh.vertex_count(),
                max_elements: mesh.element_count(),
            };
            let solver = with_cx(|cx| card.transport.prepare_uniform(cx, mesh, budget))
                .map_err(|error| format!("{what}: {error}"))?;
            let problem = ConductionProblem {
                mesh,
                boundary,
                material: card.transport.conductivity(),
                element_materials: None,
                source: source_field,
            };
            // The tangent's transpose solves take `solve`'s Krylov policy.
            let component = EnthalpyBody::new(solver, problem, patches, card.step, config.linear)
                .map_err(|error| format!("{what}: {error}"))?;
            let start = vec![card.start_j_kg; mesh.vertex_count()];
            Ok((
                Owner::Enthalpy(Box::new(component), card.transport),
                start,
                Some(card.material),
            ))
        }
    }
}

/// The owner's body on its box mesh, at a uniform start temperature.
fn build(definition: &J, materials: &Materials, what: &str) -> Result<Body, Refusal> {
    fields(
        definition,
        &[
            "kind",
            "box_m",
            "cells",
            "conductivity_w_m_k",
            "material",
            "newton",
            "volumetric_heat_capacity_j_m3_k",
            "reference_temperature_k",
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
    let source = source(definition, materials, temperature_k, config, what)?;
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
    let (owner, start, material) = owner(
        source,
        mesh,
        boundary,
        &patches,
        temperature_k,
        config,
        what,
    )?;
    Ok(Body {
        owner,
        start,
        temperature_k,
        material,
        faces: ports
            .iter()
            .map(|&(face, _, _, trace)| (face, trace))
            .collect(),
    })
}

/// Where a contact's resistance comes from.
#[derive(Debug, Clone)]
enum Origin {
    /// The scene's own number, with the scene's sentence on its source.
    Declared(String),
    /// An interface card read at a declared state: the pack's name and the
    /// members of the report's `card` object that identify the pack, the
    /// card, the claim, the receipt and the state, and say which coordinates
    /// follow the traces if the contact follows its card.
    Card { pack_id: String, report: String },
}

/// What a scene declared for a contact and what the owner made of it: kept
/// for the wiring check and the report.
#[derive(Debug, Clone)]
pub(super) struct Declared {
    /// The two traces, `element.face`, that the contact was computed for:
    /// its port `a` meets the first and its port `b` the second.
    pub(super) between: [String; 2],
    /// The scene's number, or the value its card gave at the declared state.
    resistance_m2_k_w: f64,
    origin: Origin,
    /// The sum of the interface points' conductances [W/K].
    conductance_w_k: f64,
    points: usize,
    /// The law of a contact that follows its card, if it does.
    following: Option<Following>,
}

/// A card's resistance law as the kernel element's hook: the owner's refusal
/// is the element's reason.
#[derive(Debug)]
struct CardLaw(TraceResistanceLaw);

impl ConductanceLaw for CardLaw {
    fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), String> {
        self.0
            .factor(kelvin_a, kelvin_b)
            .map_err(|error| error.to_string())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        self.0.identity_words(words);
    }
}

/// What a contact that follows its card keeps for the report: the law and
/// the element's interface points.
#[derive(Debug, Clone)]
struct Following {
    law: TraceResistanceLaw,
    points: Vec<InterfacePoint>,
}

/// The names of the card coordinates that take the trace temperatures, the
/// first trace's first.
fn axes(coordinates: &TemperatureContactCoordinates) -> Vec<String> {
    match coordinates {
        TemperatureContactCoordinates::Mean(axis) => vec![axis.clone()],
        TemperatureContactCoordinates::TwoSided(axes) => axes.to_vec(),
    }
}

impl Following {
    /// The least and the greatest resistance [m^2 K/W] the card gives over
    /// the interface points at the trace temperatures `a` and `b`, or NaN
    /// where it gives none.
    fn span(&self, a: &[f64], b: &[f64]) -> (f64, f64) {
        let at = |side: &[(usize, f64); 3], trace: &[f64]| {
            side.iter()
                .map(|&(point, weight)| weight * trace.get(point).copied().unwrap_or(f64::NAN))
                .sum::<f64>()
        };
        let reference = self.law.reference_resistance().value_m2_k_per_w();
        let resistances: Vec<f64> = self
            .points
            .iter()
            .map(|point| {
                self.law
                    .factor(at(&point.a, a), at(&point.b, b))
                    .map_or(f64::NAN, |(factor, _)| reference / factor)
            })
            .collect();
        if resistances.iter().any(|value| value.is_nan()) {
            (f64::NAN, f64::NAN)
        } else {
            span(&resistances)
        }
    }
}

impl Declared {
    /// Whether the resistance is a card's value, not the scene's own number.
    pub(super) fn is_card_backed(&self) -> bool {
        matches!(self.origin, Origin::Card { .. })
    }

    /// Whether the contact follows its card at the traces' temperatures.
    pub(super) fn follows_its_card(&self) -> bool {
        self.following.is_some()
    }

    /// The contact's member of its element in the JSON report, with its
    /// leading comma: `declared` for the scene's own number and its source,
    /// `card` for a card's value with the identities of what it was read
    /// from.
    pub(super) fn json_member(&self) -> String {
        let head = format!(
            "\"between\":[{},{}],\"resistance_m2_k_w\":{}",
            quote(&self.between[0]),
            quote(&self.between[1]),
            num(self.resistance_m2_k_w)
        );
        match &self.origin {
            Origin::Declared(source) => {
                format!(",\"declared\":{{{head},\"source\":{}}}", quote(source))
            }
            Origin::Card { report, .. } => format!(",\"card\":{{{head},{report}}}"),
        }
    }
}

/// The text report's lines for the scene's contacts: the resistances the
/// scene declares, then those read from cards. Empty when a scene has none.
pub(super) fn text_report(elements: &[Element]) -> String {
    let mut out = String::new();
    let contacts = || {
        elements
            .iter()
            .filter_map(|element| Some((&element.name, element.contact.as_ref()?)))
    };
    for (group, heading) in [
        "\nContact resistances (declared in the scene, not measured):",
        "\nContact resistances from interface cards (each read once at its declared state; \
         --json gives the receipts):",
        "\nContact resistances that follow an interface card at the temperatures of the two \
         traces (each shown at its declared state; --json gives the receipts):",
    ]
    .into_iter()
    .enumerate()
    {
        let mut first = true;
        for (name, contact) in contacts() {
            let (from, own) = match (&contact.origin, &contact.following) {
                (Origin::Declared(source), _) => (source, 0),
                (Origin::Card { pack_id, .. }, None) => (pack_id, 1),
                (Origin::Card { pack_id, .. }, Some(_)) => (pack_id, 2),
            };
            if own != group {
                continue;
            }
            if std::mem::take(&mut first) {
                let _ = writeln!(out, "{heading}");
            }
            let following = match &contact.following {
                Some(following) => {
                    let axes = axes(following.law.coordinates());
                    format!(", following {}", axes.join(" and "))
                }
                None => String::new(),
            };
            let _ = writeln!(
                out,
                "{name}: {:e} m^2 K/W between {} and {} <- {from}{following}",
                contact.resistance_m2_k_w, contact.between[0], contact.between[1]
            );
        }
    }
    out
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
    materials: &Materials,
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
        let built = build(definition, materials, &format!("elements.{name}"))?;
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
fn resistance(
    name: &str,
    between: &[String; 2],
    resistance_m2_k_w: f64,
    source: &str,
) -> Result<InterfaceResistance, String> {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RESISTANCE_PROPERTY, RESISTANCE_DIMS),
            value: PropertyValue::Scalar {
                value: resistance_m2_k_w,
                dims: RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("caller-declared system-scene contact {name}: {source}"),
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
        side(&between[0]),
        side(&between[1]),
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

/// The card coordinates that take the temperatures of a contact's two
/// traces, when the scene names them under `follows`: one name for a
/// coordinate that takes their mean at each interface point, two for one
/// that takes side a's and one that takes side b's, in the order of
/// `between`.
fn follows(definition: &J, what: &str) -> Result<Option<TemperatureContactCoordinates>, String> {
    let Some(value) = definition.get("follows") else {
        return Ok(None);
    };
    let two = match value.as_array() {
        Some([a, b]) => a.as_str().zip(b.as_str()),
        _ => None,
    };
    match (value.as_str(), two) {
        (Some(axis), _) => Ok(Some(TemperatureContactCoordinates::Mean(axis.to_owned()))),
        (None, Some((a, b))) => Ok(Some(TemperatureContactCoordinates::TwoSided([
            a.to_owned(),
            b.to_owned(),
        ]))),
        (None, None) => Err(format!(
            "{what}.follows: the name of the card coordinate that takes the mean of the two \
             traces' temperatures, or two names, the coordinate that takes the first trace's \
             and the one that takes the second's"
        )),
    }
}

/// A contact's resistance from an interface card: the compiled pack at
/// `pack`, a path relative to the scene, read at the scene's `state` under
/// `selection` (one claim only, unless the scene names another policy). The
/// owner resolves the card, so every refusal about its claims and their
/// validity is matdb's own. With the resistance, where it came from for the
/// report, and the owner's law bound at that state when the scene says
/// which coordinates follow the traces: the state is then the law's
/// reference point and the resistance its value there.
fn card(
    scene: &Path,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(InterfaceResistance, Origin, Option<TraceResistanceLaw>), String> {
    let path = definition
        .str_field("pack")
        .filter(|path| !path.trim().is_empty() && path.len() <= MAX_CARD_PACK_SOURCE_BYTES)
        .ok_or_else(|| format!("{what}.pack: an interface pack path is required"))?;
    let path = scene.parent().unwrap_or(Path::new(".")).join(path);
    let bytes = read_pack(&path, MAX_CARD_PACK_BYTES, &format!("{what}.pack"))?;
    let pack = NormalizedInterfacePack::from_bytes(&bytes)
        .map_err(|error| format!("{what}.pack: {error}"))?;
    let point = state_point(definition).map_err(|error| format!("{what}.state: {error}"))?;
    let policy = match definition.get("selection") {
        None => SelectionPolicy::SingleClaimOnly,
        Some(tag) => tag
            .as_str()
            .map(SelectionPolicy::from_tag)
            .ok_or_else(|| format!("{what}.selection: expected a policy name"))?
            .map_err(|error| format!("{what}.selection: {error}"))?,
    };
    let coordinates = follows(definition, what)?;
    let mut followed = String::new();
    let (resistance, law) = match coordinates {
        None => (
            InterfaceResistance::from_card(name, pack.card(), &point, policy)
                .map_err(|error| format!("{what}: {error}"))?,
            None,
        ),
        Some(coordinates) => {
            let names: Vec<String> = axes(&coordinates).iter().map(|axis| quote(axis)).collect();
            followed = format!(",\"follows\":[{}]", names.join(","));
            let law = TraceResistanceLaw::bind(name, pack.card(), &point, coordinates, policy)
                .map_err(|error| format!("{what}: {error}"))?;
            (law.reference_resistance().clone(), Some(law))
        }
    };
    let receipt = resistance.receipt();
    let state = point
        .axes()
        .iter()
        .map(|(axis, value)| format!("{}:{}", quote(axis), num(*value)))
        .collect::<Vec<_>>()
        .join(",");
    let report = format!(
        "\"pack_id\":{},\"pack_hash\":{},\"pack_bytes_hash\":{},\"card_hash\":{},\"resistance_receipt_hash\":{},\"selected_claim\":{},\"state\":{{{state}}}{followed}",
        quote(pack.pack_id()),
        quote(&pack.content_hash().to_hex()),
        quote(&hash_bytes(&bytes).to_hex()),
        quote(&resistance.card_identity().to_hex()),
        quote(&receipt.content_hash().to_hex()),
        quote(&receipt.selected.0.to_hex()),
    );
    let origin = Origin::Card {
        pack_id: pack.pack_id().to_owned(),
        report,
    };
    Ok((resistance, origin, law))
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
    materials: &Materials,
    bodies: &mut Vec<(String, Body)>,
    name: &str,
    definition: &J,
    scene: &Path,
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
            "pack",
            "state",
            "selection",
            "follows",
            "mapped",
        ],
        &what,
    )?;
    // One source of the resistance: the scene's number, or a card.
    let given = |keys: &[&str]| keys.iter().any(|key| definition.get(key).is_some());
    let carded = given(&["pack", "state", "selection", "follows"]);
    if carded == given(&["resistance_m2_k_w", "source"]) {
        return Err(format!(
            "{what}: a contact takes its resistance from resistance_m2_k_w with its source, or \
             from a card (pack with state): one of the two, and not both"
        )
        .into());
    }
    let between = match definition.get("between").and_then(J::as_array) {
        Some([a, b]) => a.as_str().zip(b.as_str()),
        _ => None,
    }
    .map(|(a, b)| [a, b])
    .ok_or_else(|| {
        format!("{what}: between is [\"element.face\", \"element.face\"], two traces")
    })?;
    let between = [between[0].to_owned(), between[1].to_owned()];
    let offset_m = vector(definition, "offset_m", &what)?;
    let options = refinement(definition, &what)?;
    let (body_a, port_a) = end(definitions, materials, bodies, &between[0], &what)?;
    let (body_b, port_b) = end(definitions, materials, bodies, &between[1], &what)?;
    if body_a == body_b {
        return Err(format!(
            "{what}: between names two traces of one body; a contact joins two bodies"
        )
        .into());
    }
    let (resistance, resistance_m2_k_w, origin, law) = if carded {
        let (resistance, origin, law) = card(scene, name, definition, &what)?;
        let value = resistance.value_m2_k_per_w();
        (resistance, value, origin, law)
    } else {
        let source = definition
            .str_field("source")
            .filter(|source| !source.trim().is_empty())
            .ok_or_else(|| format!("{what}: source says where resistance_m2_k_w comes from"))?;
        let value = positive(definition, "resistance_m2_k_w", &what)?;
        let resistance = resistance(name, &between, value, source)
            .map_err(|error| format!("{what}: {error}"))?;
        (resistance, value, Origin::Declared(source.to_owned()), None)
    };
    let mut declared = Declared {
        between,
        resistance_m2_k_w,
        origin,
        conductance_w_k: 0.0,
        points: 0,
        following: None,
    };
    let (a, b) = (
        bodies[body_a].1.capacity(&what)?,
        bodies[body_b].1.capacity(&what)?,
    );
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
    let (wall, following) =
        kernel_element(&joint, law).map_err(|error| format!("{what}: {error}"))?;
    declared.following = following;
    Ok((wall, declared))
}

/// The kernel's element for the owner's contact `joint`: its interface
/// points, and, for a contact that follows its card, `law` on the element's
/// hook. The card is then asked at the interface points and, as the
/// owner's pointwise contact asks it, at every pair of trace points that
/// face each other: a mapped contact has none.
fn kernel_element(
    joint: &TraceContact,
    law: Option<TraceResistanceLaw>,
) -> Result<(ConductanceMatrix, Option<Following>), ElementError> {
    let points: Vec<InterfacePoint> = joint
        .points
        .iter()
        .map(|point| InterfacePoint {
            a: point.a,
            b: point.b,
            conductance_w_k: point.conductance_w_k,
        })
        .collect();
    let Some(law) = law else {
        let wall = ConductanceMatrix::new(joint.width_a, joint.width_b, points)?;
        return Ok((wall, None));
    };
    let alone = |point: usize| [(point, 1.0), (point, 0.0), (point, 0.0)];
    let rows = joint
        .facing
        .iter()
        .map(|&[a, b]| GuardRow {
            a: alone(a),
            b: alone(b),
        })
        .collect();
    let wall = ConductanceMatrix::new(joint.width_a, joint.width_b, points.clone())?
        .following(Arc::new(CardLaw(law.clone())))
        .guarded_at(rows)?;
    Ok((wall, Some(Following { law, points })))
}

/// Every `trace_contact` of the scene with the bodies it names. A scene
/// with no contact prepares nothing.
pub(super) fn prepare(
    definitions: &[(String, J)],
    materials: &Materials,
    scene: &Path,
) -> Result<Prepared, Refusal> {
    let mut prepared = Prepared::default();
    for (name, definition) in definitions {
        if definition.str_field("kind") == Some("trace_contact") {
            let (wall, declared) = contact(
                definitions,
                materials,
                &mut prepared.bodies,
                name,
                definition,
                scene,
            )?;
            prepared.contacts.push((name.clone(), wall, declared));
        }
    }
    Ok(prepared)
}

/// The parts of a `conduction_body` the scene keeps: its node and ports, its
/// start temperature, the scene material whose conductivity it takes, if
/// any, and the paired material whose chart reads its state, if that state
/// is specific enthalpy.
pub(super) struct Added {
    pub(super) node: NodeId,
    pub(super) ports: Vec<PortSpec>,
    pub(super) temperature_k: f64,
    pub(super) material: Option<String>,
    pub(super) enthalpy: Option<&'static MaterialEnthalpyTransport>,
}

/// A `conduction_body` element: the owner's body on its box mesh, at a
/// uniform start temperature.
pub(super) fn body(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    materials: &Materials,
    prepared: &mut Prepared,
    what: &str,
) -> Result<Added, Refusal> {
    let body = match prepared.bodies.iter().position(|(known, _)| known == name) {
        Some(built) => prepared.bodies.swap_remove(built).1,
        None => build(definition, materials, what)?,
    };
    let ((node, ports), enthalpy) = match body.owner {
        Owner::Capacity(component) => (add(builder, name, *component, body.start)?, None),
        Owner::Enthalpy(component, transport) => {
            (add(builder, name, *component, body.start)?, Some(transport))
        }
    };
    Ok(Added {
        node,
        ports,
        temperature_k: body.temperature_k,
        material: body.material,
        enthalpy,
    })
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
/// and greatest nodal temperatures. A body whose state is specific enthalpy
/// (`enthalpy`) reads each node's temperature through the card's chart, as
/// the owner reads it, and adds its frozen reference density.
pub(super) fn values(
    system: &System,
    node: NodeId,
    ports: &[PortSpec],
    enthalpy: Option<&MaterialEnthalpyTransport>,
) -> Vec<(String, f64)> {
    let mut values = Vec::with_capacity(ports.len() + 3);
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
    let (low, high) = match (field, enthalpy) {
        (None, _) => (f64::NAN, f64::NAN),
        (Some(field), None) => span(field),
        (Some(field), Some(transport)) => {
            let chart = transport.caloric().curve();
            field
                .iter()
                .map(|&h| chart.temperature_at_specific_enthalpy(h))
                .collect::<Result<Vec<f64>, _>>()
                .map_or((f64::NAN, f64::NAN), |temperatures| span(&temperatures))
        }
    };
    values.push(("min_temperature_k".to_owned(), low));
    values.push(("max_temperature_k".to_owned(), high));
    if let Some(transport) = enthalpy {
        values.push((
            "reference_density_kg_m3".to_owned(),
            transport.reference_density_kg_m3(),
        ));
    }
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
/// scene declared and the number of interface points. A contact that
/// follows its card adds the least and the greatest resistance the card
/// gave over the interface points at the temperatures the last step ended
/// with; the conductance and the resistance before them are the card's at
/// the declared state.
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
    let mut values = vec![
        ("heat_rate_w".to_owned(), rate),
        ("conductance_w_k".to_owned(), declared.conductance_w_k),
        ("resistance_m2_k_w".to_owned(), declared.resistance_m2_k_w),
        ("interface_points".to_owned(), points),
    ];
    if let Some(following) = &declared.following {
        let received = |id: &str| {
            system
                .port(node, id)
                .filter(|_| system.steps() > 0)
                .map(|port| system.last_input(port))
        };
        let (low, high) = received("a")
            .zip(received("b"))
            .map_or((f64::NAN, f64::NAN), |(a, b)| following.span(a, b));
        values.push(("min_resistance_m2_k_w".to_owned(), low));
        values.push(("max_resistance_m2_k_w".to_owned(), high));
    }
    values
}
