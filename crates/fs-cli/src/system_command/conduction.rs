//! `conduction_body`: fs-conduction's backward-Euler finite-element body,
//! [`ConductionBody`], as a scene element.
//!
//! The scene authors a box `[0, Lx] x [0, Ly] x [0, Lz]` meshed as an
//! `nx x ny x nz` Kuhn grid (fs-conduction's `fixtures::box_grid`, six
//! tetrahedra a cell), one declared conductivity and volumetric heat
//! capacity, the start temperature, and the box faces that are heat ports.
//! The faces not named are insulated. Each port receives the heat rate into
//! the body and returns its face's area-mean end temperature; the owner's
//! own backward Euler steps the field. Nothing here conducts, stores or
//! integrates heat.
//!
//! The kernel keeps component states as `'static` values and the body
//! borrows its mesh, boundary and material, so they are leaked: they live
//! as long as the process that runs the scene.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::component::ConductionBody;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField, ThermalBc,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::PortSpec;
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

/// Three positive finite numbers at `key`.
fn triple(definition: &J, key: &str, what: &str) -> Result<[f64; 3], String> {
    let values: Vec<f64> = definition
        .get(key)
        .and_then(J::as_array)
        .filter(|values| values.len() == 3)
        .ok_or_else(|| format!("{what}: {key} is an array of three numbers"))?
        .iter()
        .map(J::as_f64)
        .collect::<Option<_>>()
        .ok_or_else(|| format!("{what}: {key} is an array of three numbers"))?;
    if values
        .iter()
        .any(|value| !(value.is_finite() && *value > 0.0))
    {
        return Err(format!("{what}: {key} must be positive and finite"));
    }
    Ok([values[0], values[1], values[2]])
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

/// The faces that are ports, each once: name, axis and far end.
fn ports(definition: &J, what: &str) -> Result<Vec<(&'static str, usize, bool)>, String> {
    let names = definition
        .get("ports")
        .and_then(J::as_array)
        .filter(|names| !names.is_empty())
        .ok_or_else(|| format!("{what}: ports is a nonempty array of box faces"))?;
    let mut ports = Vec::with_capacity(names.len());
    for name in names {
        let face = FACES
            .iter()
            .find(|(face, _, _)| Some(*face) == name.as_str())
            .ok_or_else(|| {
                format!("{what}: ports name box faces: x_min, x_max, y_min, y_max, z_min, z_max")
            })?;
        if ports.contains(face) {
            return Err(format!("{what}: port {} is named twice", face.0));
        }
        ports.push(*face);
    }
    Ok(ports)
}

/// A `conduction_body` element: the owner's body on its box mesh, at a
/// uniform start temperature, which it also returns.
pub(super) fn body(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>, f64), Refusal> {
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
    let mut faces = ThermalBoundaryBuilder::new(mesh);
    for &(port, axis, far) in &ports {
        let at = if far { extent[axis] } else { 0.0 };
        faces = faces
            .region(
                port,
                |face| on_box_face(face.centroid[axis], at),
                ThermalBc::adiabatic(),
            )
            .map_err(|error| format!("{what}: {error}"))?;
    }
    let boundary = Box::leak(Box::new(
        faces
            .adiabatic_remainder()
            .finish()
            .map_err(|error| format!("{what}: {error}"))?,
    ));
    // The owner prepares its capacitance once, with a context of its own:
    // the prepared solver keeps no part of it.
    let gate = CancelGate::new_clock_free();
    let solver = ArenaPool::new(ArenaConfig::default())
        .scope(|arena| {
            let key = StreamKey {
                seed: 0,
                kernel_id: 0x5953,
                tile: 0,
                iteration: 0,
            };
            let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
            BackwardEuler::uniform(&cx, mesh, capacity)
        })
        .map_err(|error| format!("{what}: {error}"))?;
    let problem = ConductionProblem {
        mesh,
        boundary,
        material: Box::leak(Box::new(conductivity)),
        element_materials: None,
        source: Box::leak(Box::new(ScalarField::Uniform(0.0))),
    };
    let names: Vec<&str> = ports.iter().map(|&(port, _, _)| port).collect();
    let component = ConductionBody::new(solver, problem, &names, config)
        .map_err(|error| format!("{what}: {error}"))?;
    let (node, specs) = add(
        builder,
        name,
        component,
        vec![temperature_k; mesh.vertex_count()],
    )?;
    Ok((node, specs, temperature_k))
}

/// Each port's area-mean temperature at the end of the last step, and the
/// field's least and greatest nodal temperatures.
pub(super) fn values(system: &System, node: NodeId, ports: &[PortSpec]) -> Vec<(String, f64)> {
    let mut values = Vec::with_capacity(ports.len() + 2);
    for port in ports {
        let id = port.id().as_str();
        let mean = system
            .port(node, id)
            .filter(|_| system.steps() > 0)
            .and_then(|port| system.last_output(port).first().copied())
            .unwrap_or(f64::NAN);
        values.push((format!("{id}.temperature_k"), mean));
    }
    let field = system.state::<State>(node);
    let (low, high) = field.map_or((f64::NAN, f64::NAN), |field| {
        field
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), &t| {
                (low.min(t), high.max(t))
            })
    });
    values.push(("min_temperature_k".to_owned(), low));
    values.push(("max_temperature_k".to_owned(), high));
    values
}
