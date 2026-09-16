//! File-driven nominal coupled cooling. This is an explicit experimental JSON
//! request, not a .fsim migration or a ledger-backed `solve` result. The binary
//! adapter reuses the CLI's single JSON reader and existing numerical producers.

#[path = "json_read.rs"]
mod json;
mod design;
mod solid_data;
mod objective;
mod fan_drive;
mod convection;
mod fan_speed;
mod contacts;
mod transient;
mod acceleration;

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fmt::{self, Write as _};
use std::fs::File;
use std::io::Read;
use std::time::{Duration, Instant};

use acceleration::solve_coupled_transport;
use fs_airflow::conjugate::{AirSegment, ConjugateConfig, Relaxation, SolidRegionState};
use fs_airflow::graph::thermal::coupled_transport::CoupledTransportSolution;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::{CoupledGradient, CoupledLinearization, InterfaceSolveConfig};
use fs_airflow::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolution, GraphSolveConfig, LossGraph};
use fs_airflow::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_cli::{CommandOutput, exit};
use fs_conduction::adjoint::robin::RobinLinearization;
use fs_conduction::{ConductionMesh, ConductionProblem, InitialGuess, ScalarField, SolveConfig, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};
use fs_rep_mesh::TetComplex;
use json::JsonValue as J;

const MAX_INPUT_BYTES: u64 = 16 * 1024 * 1024;
const SCHEMA: &str = "frankensim.cooling-network.v1";
const RESULT_SCHEMA: &str = "frankensim.cooling-network.result.v1";
const NO_CLAIM: &str = "nominal fixed-geometry solid model; hydraulics and convection coefficients frozen within each thermal solve; caller-declared constant isotropic/anisotropic or bounded scalar k(T) materials and frozen fluid properties; k(T) transients require explicit Newton/Armijo settings and temperature-independent heat capacity; coefficients declared or derived from validity-gated duct correlations, without coupled boundary-layer evolution; explicit matching-P1 contacts have fixed caller-declared resistance; component sources use nodal P1 support; maxima concern the discrete field only; no CFD, recirculation, fan heating, nonmatching contact, radiation, uncertainty certification, mesh-convergence or experimental-validation claim; not a .fsim or ledger-backed solve";
const HELP: &str = "Usage: frankensim [--json] cooling-network <request.json>\n\nSolve a prescribed-pressure or fan-driven network and heterogeneous solid,\nincluding component heating, directional conductivity, bounded scalar k(T),\nfinite-resistance thermal contacts, downstream mixing and declared or\nflow-derived duct convection. Compute mean/peak temperatures, conditional\nthermal gradients (including contact resistance), and effective-h or full\nfan-speed target searches. All quantities use coherent SI. Transient k(T)\nrequires explicit transient.nonlinear settings; heat capacity remains constant.\nRequest schema: frankensim.cooling-network.v1.\n\nSee examples/cooling-network/README.md, MATERIAL_COOLING.md, FAN_COOLING.md,\nNONLINEAR_TRANSIENT_COOLING.md and CONTACT_COOLING.md. Results are nominal\nestimates, not validated hardware or ledger-backed .fsim runs.\n";

type Result<T> = std::result::Result<T, Failure>;
#[derive(Debug)]
struct Failure { code: &'static str, message: String }
impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}: {}", self.code, self.message) }
}
impl std::error::Error for Failure {}
fn bad(message: impl Into<String>) -> Failure { Failure { code: "cooling-network-input", message: message.into() } }
fn producer(error: impl fmt::Display) -> Failure { Failure { code: "cooling-network-solve", message: error.to_string() } }
fn poll(cx: &Cx<'_>) -> Result<()> { cx.checkpoint().map_err(|_| Failure { code: "cooling-network-cancelled", message: "work cancelled before publication".into() }) }

#[derive(Debug)]
struct Surface {
    name: String, faces: BTreeSet<[u32; 3]>, h: f64, area: f64,
    convection: Option<convection::Law>,
}
#[derive(Debug, Clone, Copy)]
struct Limits {
    graph: usize, coupling: usize, linear: usize, derivative: usize, wall_seconds: f64,
    flow: f64, heat: f64, temperature: f64, relative: f64, relaxation: f64,
}
#[derive(Debug)]
struct Request {
    seed: u64, graph: LossGraph, boundaries: Vec<FixedPressure>, inlets: Vec<TransportInlet>,
    region_paths: Vec<Vec<String>>, air: TransportAir, mesh: ConductionMesh,
    surfaces: Vec<Surface>, conductivity: f64, source: f64, adiabatic: bool,
    solid_data: solid_data::SolidData,
    contacts: Option<contacts::Contacts>,
    fan: Option<fan_drive::FanDrive>,
    fan_speed_design: Option<fan_speed::FanSpeedDesign>,
    transient: Option<transient::Schedule>,
    objective: objective::Objective, gradient: bool, limits: Limits, design: Option<design::DesignRequest>,
}
#[derive(Debug)]
struct Evaluation {
    coupled: CoupledTransportSolution, temperatures: Vec<f64>, gradient: Option<CoupledGradient>,
    objective: f64, objective_state: objective::ObjectiveState,
    robin_total_w: f64, source_total_w: f64, htc: Vec<f64>,
    convection: Vec<convection::Derived>,
    contact_fluxes: Vec<fs_conduction::InterfaceFlux>,
}

fn object<'a>(value: &'a J, allowed: &[&str], path: &str) -> Result<&'a J> {
    let members = value.as_object().ok_or_else(|| bad(format!("{path} must be an object")))?;
    for (name, _) in members {
        if !allowed.contains(&name.as_str()) { return Err(bad(format!("unknown field {path}.{name}"))); }
    }
    Ok(value)
}
fn get<'a>(value: &'a J, key: &str) -> Result<&'a J> { value.get(key).ok_or_else(|| bad(format!("missing required field {key}"))) }
fn string(value: &J, field: &str) -> Result<String> {
    let s = value.as_str().ok_or_else(|| bad(format!("{field} must be a string")))?;
    if s.is_empty() || s.len() > 256 || s.trim() != s || s.chars().any(char::is_control) {
        return Err(bad(format!("{field} must be a nonempty, trimmed, control-free string of at most 256 bytes")));
    }
    Ok(s.to_string())
}
fn number(value: &J, field: &str) -> Result<f64> {
    value.as_f64().filter(|n| n.is_finite()).ok_or_else(|| bad(format!("{field} must be a finite number")))
}
fn positive(value: &J, field: &str) -> Result<f64> {
    let n = number(value, field)?;
    if n <= 0.0 { Err(bad(format!("{field} must be positive"))) } else { Ok(n) }
}
fn integer(value: &J, field: &str, max: usize) -> Result<usize> {
    let n = integer_raw(value, field)?;
    if n > max { Err(bad(format!("{field} exceeds {max}"))) } else { Ok(n) }
}
fn integer_raw(value: &J, field: &str) -> Result<usize> {
    value.number_raw().and_then(|s| s.parse::<usize>().ok())
        .ok_or_else(|| bad(format!("{field} must use a nonnegative integer JSON spelling")))
}
fn count(value: &J, field: &str, max: usize) -> Result<usize> {
    let n = integer(value, field, max)?;
    if n == 0 { Err(bad(format!("{field} must be positive"))) } else { Ok(n) }
}
fn boolean(value: &J, field: &str) -> Result<bool> {
    match value { J::Bool(v) => Ok(*v), _ => Err(bad(format!("{field} must be boolean"))) }
}
fn array<'a>(value: &'a J, field: &str, max: usize) -> Result<&'a [J]> {
    let items = value.as_array().ok_or_else(|| bad(format!("{field} must be an array")))?;
    if items.len() > max { return Err(bad(format!("{field} has more than {max} entries"))); }
    Ok(items)
}
fn indices<const N: usize>(value: &J, field: &str, vertices: usize) -> Result<[u32; N]> {
    let items = array(value, field, N)?;
    if items.len() != N { return Err(bad(format!("{field} requires {N} vertex indices"))); }
    let mut out = [0_u32; N];
    for (slot, item) in out.iter_mut().zip(items) {
        let n = integer(item, field, vertices.saturating_sub(1))?;
        *slot = u32::try_from(n).map_err(|_| bad("vertex index exceeds u32"))?;
    }
    let mut unique = out;
    unique.sort_unstable();
    if unique.windows(2).any(|p| p[0] == p[1]) { return Err(bad(format!("{field} repeats a vertex"))); }
    Ok(out)
}

impl Request {
    fn parse(text: &str) -> Result<Self> {
        if text.len() as u64 > MAX_INPUT_BYTES { return Err(bad("request exceeds 16 MiB")); }
        let root = J::parse(text).map_err(bad_parse)?;
        object(&root, &["schema", "units", "seed", "budgets", "tolerances", "air", "hydraulics", "solid", "objective", "design", "fan_speed_design", "transient"], "request")?;
        if get(&root, "schema")?.as_str() != Some(SCHEMA) || get(&root, "units")?.as_str() != Some("SI") {
            return Err(bad("expected schema frankensim.cooling-network.v1 and units SI"));
        }
        let seed = string(get(&root, "seed")?, "seed")?.parse::<u64>().map_err(|_| bad("seed must be a decimal u64 string"))?;
        let b = object(get(&root, "budgets")?, &["graph_sweeps", "coupling_iterations", "linear_iterations", "derivative_iterations", "wall_seconds"], "budgets")?;
        let t = object(get(&root, "tolerances")?, &["flow_m3_s", "heat_w", "temperature_k", "linear_relative", "relaxation"], "tolerances")?;
        let limits = Limits {
            graph: count(get(b, "graph_sweeps")?, "graph_sweeps", 100_000)?,
            coupling: count(get(b, "coupling_iterations")?, "coupling_iterations", 10_000)?,
            linear: count(get(b, "linear_iterations")?, "linear_iterations", 100_000)?,
            derivative: count(get(b, "derivative_iterations")?, "derivative_iterations", 10_000)?,
            wall_seconds: positive(get(b, "wall_seconds")?, "wall_seconds")?,
            flow: positive(get(t, "flow_m3_s")?, "flow_m3_s")?,
            heat: positive(get(t, "heat_w")?, "heat_w")?,
            temperature: positive(get(t, "temperature_k")?, "temperature_k")?,
            relative: positive(get(t, "linear_relative")?, "linear_relative")?,
            relaxation: positive(get(t, "relaxation")?, "relaxation")?,
        };
        if limits.wall_seconds > 3600.0 || limits.relative >= 1.0 || limits.relaxation > 1.0 {
            return Err(bad("wall_seconds must be <=3600, linear_relative <1, and relaxation <=1"));
        }
        let a = object(get(&root, "air")?, &["density_kg_m3", "specific_heat_j_kg_k"], "air")?;
        let air = TransportAir { density: Density::new(positive(get(a, "density_kg_m3")?, "density_kg_m3")?),
            specific_heat_j_kg_k: positive(get(a, "specific_heat_j_kg_k")?, "specific_heat_j_kg_k")? };
        let s = object(get(&root, "solid")?, &["vertices_m", "tetrahedra", "conductivity_w_m_k", "materials", "element_materials", "source_w_m3", "component_power", "adiabatic_remainder", "surfaces", "contacts"], "solid")?;
        let mut positions = Vec::new();
        for point in array(get(s, "vertices_m")?, "vertices_m", 20_000)? {
            let xyz = array(point, "vertex", 3)?;
            if xyz.len() != 3 { return Err(bad("each vertex requires three metre coordinates")); }
            positions.push([number(&xyz[0], "x")?, number(&xyz[1], "y")?, number(&xyz[2], "z")?]);
        }
        if positions.len() < 4 { return Err(bad("at least four vertices required")); }
        let mut tets = Vec::new();
        let mut unique_tets = BTreeSet::new();
        let mut face_counts = BTreeMap::new();
        for entry in array(get(s, "tetrahedra")?, "tetrahedra", 100_000)? {
            let tet = indices::<4>(entry, "tetrahedron", positions.len())?;
            let mut key = tet;
            key.sort_unstable();
            if !unique_tets.insert(key) { return Err(bad("duplicate tetrahedron")); }
            for omitted in 0..4 {
                let face: Vec<_> = key.iter().enumerate().filter_map(|(i, &v)| (i != omitted).then_some(v)).collect();
                let count = face_counts.entry(face).or_insert(0_usize);
                *count += 1;
                if *count > 2 { return Err(bad("a tetrahedral face has more than two incident cells")); }
            }
            tets.push(tet);
        }
        if tets.is_empty() { return Err(bad("at least one tetrahedron required")); }
        let used: BTreeSet<_> = tets.iter().flatten().copied().collect();
        if used.len() != positions.len() { return Err(bad("unused solid vertices would create unanchored degrees of freedom")); }
        let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).map_err(producer)?;
        let exterior: BTreeMap<_, _> = mesh.boundary().iter().map(|f| (f.vertices, f.area)).collect();
        let mut names = BTreeSet::new();
        let mut owned_faces = BTreeSet::new();
        let mut surfaces = Vec::new();
        for entry in array(get(s, "surfaces")?, "surfaces", 4096)? {
            object(entry, &["name", "faces", "htc_w_m2_k", "convection"], "surface")?;
            let name = string(get(entry, "name")?, "surface.name")?;
            if !names.insert(name.clone()) { return Err(bad(format!("duplicate surface {name}"))); }
            let (h, convection) = convection::parse(entry)?;
            let mut faces = BTreeSet::new();
            for entry in array(get(entry, "faces")?, "surface.faces", 200_000)? {
                let mut face = indices::<3>(entry, "surface face", mesh.vertex_count())?;
                face.sort_unstable();
                if !exterior.contains_key(&face) || !owned_faces.insert(face) {
                    return Err(bad(format!("surface {name} contains a non-exterior or multiply owned face {face:?}")));
                }
                faces.insert(face);
            }
            if faces.is_empty() { return Err(bad(format!("surface {name} has no faces"))); }
            // Same order as the production boundary integrator and Robin ports.
            let area: f64 = mesh.boundary().iter().filter(|f| faces.contains(&f.vertices)).map(|f| f.area).sum();
            surfaces.push(Surface { name, faces, h, area, convection });
        }
        if surfaces.is_empty() { return Err(bad("at least one heat-exchanging surface required")); }
        let adiabatic = boolean(get(s, "adiabatic_remainder")?, "adiabatic_remainder")?;
        let contacts = contacts::Contacts::parse(s.get("contacts"), &mesh, &surfaces, adiabatic)?;
        let (solid_data, conductivity, source) = solid_data::SolidData::parse(s, &mesh)?;
        let hyd = object(get(&root, "hydraulics")?, &["node_count", "boundaries", "fan", "branches"], "hydraulics")?;
        let node_count = count(get(hyd, "node_count")?, "node_count", 4096)?;
        let (boundaries, inlets, fan) = fan_drive::parse(hyd, node_count)?;
        let mut branches = Vec::new();
        let mut region_paths = Vec::new();
        let mut owned_regions = BTreeSet::new();
        for entry in array(get(hyd, "branches")?, "branches", 16_384)? {
            object(entry, &["name", "from", "to", "resistance_pa_s2_m6", "source", "regions"], "branch")?;
            let name = string(get(entry, "name")?, "branch.name")?;
            let from = integer(get(entry, "from")?, "branch.from", node_count - 1)?;
            let to = integer(get(entry, "to")?, "branch.to", node_count - 1)?;
            let resistance = positive(get(entry, "resistance_pa_s2_m6")?, "resistance_pa_s2_m6")?;
            let source = string(get(entry, "source")?, "branch.source")?;
            // Zero declared coefficient uncertainty is not a physical certificate.
            let loss = LossElement::new(name.clone(), LossResistance::new(resistance), 0.0,
                SourceProvenance::new(source, format!("{SCHEMA}:{name}")), ToleranceBasis::EngineeringAllowance).map_err(producer)?;
            branches.push(GraphBranch { from, to, loss });
            let mut path = Vec::new();
            for row in array(get(entry, "regions")?, "branch.regions", surfaces.len())? {
                let region = string(row, "branch region")?;
                if !names.contains(&region) || !owned_regions.insert(region.clone()) {
                    return Err(bad(format!("unknown or multiply owned thermal region {region}")));
                }
                path.push(region);
            }
            region_paths.push(path);
        }
        if owned_regions != names { return Err(bad("every solid surface must belong to exactly one branch")); }
        let graph = LossGraph::new(node_count, branches).map_err(producer)?;
        let o = get(&root, "objective")?;
        let objective = objective::Objective::parse(o, &surfaces, &mesh)?;
        let gradient = boolean(get(o, "gradient")?, "gradient")?;
        if root.get("design").is_some() && root.get("fan_speed_design").is_some() {
            return Err(bad("choose effective-h design or fan-speed design, not both"));
        }
        let fan_speed_design = root.get("fan_speed_design").map(|value| {
            let drive = fan.as_ref().ok_or_else(|| bad("fan_speed_design requires hydraulics.fan"))?;
            fan_speed::FanSpeedDesign::parse(value, drive)
        }).transpose()?;
        let design = root.get("design").map(|value| design::DesignRequest::parse(value, &names, objective.is_mean())).transpose()?;
        if design.is_some() && !gradient { return Err(bad("design requires objective.gradient=true")); }
        if let Some(target) = root.get("design").and_then(|d| d.str_field("surface")) {
            if surfaces.iter().any(|s| s.name == target && s.convection.is_some()) {
                return Err(bad("effective-h design cannot override a flow-derived convection law"));
            }
        }
        let transient = root.get("transient").map(|value| transient::Schedule::parse(
            value, mesh.vertex_count(), mesh.element_count(), fan.as_ref(),
        )).transpose()?;
        if transient.is_some() && (gradient || design.is_some() || fan_speed_design.is_some()) {
            return Err(bad("transient requires gradient=false and no steady design search"));
        }
        Ok(Self { seed, graph, boundaries, inlets, region_paths, air, mesh, surfaces,
            conductivity, source, adiabatic, solid_data, contacts, fan, fan_speed_design, transient, objective, gradient, limits, design })
    }

    fn flow(&self, cx: &Cx<'_>) -> Result<GraphSolution> {
        if let Some(fan) = &self.fan {
            return fan.solve(cx, &self.graph, self.limits, fan.speed_ratio);
        }
        self.graph.solve(&self.boundaries, GraphSolveConfig { max_sweeps: self.limits.graph,
            max_node_iterations: 80, absolute_flow_tolerance: VolumetricFlowRate::new(self.limits.flow),
            relative_flow_tolerance: 0.0 }, cx).map_err(producer)
    }

    fn boundary(&self, names: &[&str], references: &[f64], htc: &BTreeMap<String, f64>) -> Result<ThermalBoundary> {
        if names.len() != references.len() { return Err(bad("internal reference arity mismatch")); }
        let refs: BTreeMap<_, _> = names.iter().copied().zip(references.iter().copied()).collect();
        let mut builder = ThermalBoundaryBuilder::new(&self.mesh);
        for surface in &self.surfaces {
            let reference = *refs.get(surface.name.as_str()).ok_or_else(|| bad("missing surface reference"))?;
            let h = *htc.get(&surface.name).ok_or_else(|| bad("missing surface coefficient"))?;
            builder = builder.region(&surface.name, |face| surface.faces.contains(&face.vertices),
                ThermalBc::robin(h, reference).map_err(producer)?).map_err(producer)?;
        }
        // Contacts own the untagged paired faces. Parse already checked every
        // OTHER face when an adiabatic remainder was not requested.
        if self.adiabatic || self.contacts.is_some() { builder = builder.adiabatic_remainder(); }
        builder.finish().map_err(producer)
    }

    fn transport<'f>(&self, cx: &Cx<'_>, flow: &'f GraphSolution, htc: &BTreeMap<String, f64>) -> Result<TransportNetwork<'f>> {
        poll(cx)?;
        let by_name: BTreeMap<_, _> = self.surfaces.iter().map(|s| (s.name.as_str(), s)).collect();
        let models = self.region_paths.iter().map(|path| {
            if path.is_empty() { return Ok(BranchThermalModel::Adiabatic); }
            path.iter().map(|name| {
                let surface = by_name[name.as_str()];
                AirSegment::new(name, surface.area, htc[name]).map_err(producer)
            }).collect::<Result<Vec<_>>>().map(BranchThermalModel::Exchange)
        }).collect::<Result<Vec<_>>>()?;
        TransportNetwork::new(cx, flow, self.air, models, &self.inlets, TransportConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(self.limits.flow), relative_flow_tolerance: 0.0,
            absolute_heat_tolerance_w: self.limits.heat, relative_heat_tolerance: 0.0,
        }).map_err(producer)
    }

    fn evaluate(&self, cx: &Cx<'_>, flow: &GraphSolution, htc: &BTreeMap<String, f64>, want_gradient: bool) -> Result<Evaluation> {
        poll(cx)?;
        let (coefficients, convection) = convection::resolve(self, cx, flow, htc)?;
        let htc = &coefficients;
        let network = self.transport(cx, flow, htc)?;
        let names = network.regions();
        let material = fs_conduction::ConductivityModel::isotropic_declared(self.conductivity).map_err(producer)?;
        let uniform_source = ScalarField::Uniform(self.source);
        let source = self.solid_data.nodal_source.as_ref().unwrap_or(&uniform_source);
        let gate = ConjugateConfig { max_iterations: self.limits.coupling,
            temperature_tolerance_k: self.limits.temperature, balance_tolerance_w: self.limits.heat,
            balance_relative_tolerance: 0.0, relaxation: Relaxation::Fixed { omega: self.limits.relaxation } };
        let mut final_linear = None;
        let mut failure = None;
        let coupled = solve_coupled_transport(cx, &network, &gate, |cx, references| {
            let evaluated = (|| -> Result<Vec<SolidRegionState>> {
                let boundary = self.boundary(&names, references, htc)?;
                let mut config = SolveConfig::default();
                config.initial = InitialGuess::Uniform(references.iter().sum::<f64>() / references.len() as f64);
                config.linear.tolerance = self.limits.relative;
                config.linear.max_iterations = self.limits.linear;
                config.stop.residual_rtol = self.limits.relative;
                config.stop.step_atol = 0.0;
                let problem = ConductionProblem { mesh: &self.mesh, boundary: &boundary, material: &material,
                    element_materials: self.solid_data.element_materials.as_ref(), source };
                let linear = match &self.contacts {
                    Some(contacts) => RobinLinearization::new_with_interfaces(cx, problem, &contacts.interfaces, config, &names),
                    None => RobinLinearization::new(cx, problem, config, &names),
                }.map_err(producer)?;
                let states = names.iter().map(|name| linear.primal().report.robin_fluxes.iter()
                    .find(|flux| flux.region == *name).map(SolidRegionState::from_robin_flux)
                    .ok_or_else(|| bad(format!("solid report lacks surface {name}"))))
                    .collect::<Result<Vec<_>>>()?;
                final_linear = Some(linear);
                Ok(states)
            })();
            evaluated.map_err(|error| {
                failure = Some(error);
                AirflowError::Cancelled { iteration: 0, references_k: references.to_vec() }
            })
        });
        if let Some(error) = failure { return Err(error); }
        let coupled = coupled.map_err(producer)?;
        for derived in &convection { derived.check_direction(&coupled.solid, self.limits.heat)?; }
        let linear = final_linear.ok_or_else(|| bad("coupling produced no solid field"))?;
        let objective_state = self.objective.evaluate(cx, &linear.primal().temperature, &coupled.solid)?;
        let objective = objective_state.value;
        let gradient = if want_gradient {
            let binding = CoupledLinearization::new(cx, &network, &linear, &gate).map_err(producer)?;
            let mut weights = binding.zero_objective();
            objective_state.seed(&mut weights);
            Some(binding.pullback_iqn(cx, &weights, InterfaceSolveConfig { max_iterations: self.limits.derivative,
                absolute_tolerance: self.limits.relative, relative_tolerance: self.limits.relative,
                relaxation: self.limits.relaxation }, acceleration::POLICY).map_err(producer)?)
        } else { None };
        let total_solid: f64 = coupled.solid.iter().map(|s| s.heat_rate_w).sum();
        let robin = linear.primal().report.energy.robin_out_w;
        let source_w = linear.primal().report.energy.source_w;
        if !total_solid.is_finite() || !robin.is_finite() || !source_w.is_finite()
            || (total_solid - robin).abs() > self.limits.heat || (robin - source_w).abs() > self.limits.heat {
            return Err(producer("whole-domain solid energy/decomposition balance missed the declared watt tolerance"));
        }
        if let Some(audit) = &self.solid_data.power {
            if (source_w - audit.delivered_total_w()).abs() > self.limits.heat {
                return Err(producer("assembled solid source disagrees with the component power map"));
            }
        }
        let contact_fluxes = self.contacts.as_ref().map(|contacts|
            contacts.interfaces.fluxes(&linear.primal().temperature).map_err(producer))
            .transpose()?.unwrap_or_default();
        poll(cx)?;
        Ok(Evaluation { coupled, temperatures: linear.primal().temperature.clone(), gradient,
            objective, objective_state, robin_total_w: robin, source_total_w: source_w,
            htc: linear.ports().iter().map(|port| port.htc_w_m2_k).collect(), convection, contact_fluxes })
    }
}
fn bad_parse(error: json::JsonReadError) -> Failure { bad(error.to_string()) }

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""), '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"), '\r' => out.push_str("\\r"), '\t' => out.push_str("\\t"),
            c if c < '\u{20}' => { let _ = write!(out, "\\u{:04x}", c as u32); }
            c => out.push(c),
        }
    }
    out.push('"'); out
}
fn num(n: f64) -> Result<String> {
    if n.is_finite() { Ok(n.to_string()) } else { Err(producer("nonfinite result cannot be published as JSON")) }
}
fn numbers(values: &[f64]) -> Result<String> {
    Ok(format!("[{}]", values.iter().map(|&v| num(v)).collect::<Result<Vec<_>>>()?.join(",")))
}
fn optional(n: Option<f64>) -> Result<String> { n.map_or_else(|| Ok("null".into()), num) }

fn render(request: &Request, flow: &GraphSolution, evaluated: &Evaluation) -> Result<String> {
    let mean = request.objective.is_mean();
    let mut walls = Vec::new();
    for (i, state) in evaluated.coupled.solid.iter().enumerate() {
        walls.push(format!("{{\"region\":{},\"area_m2\":{},\"htc_w_m2_k\":{},\"mean_temperature_k\":{},\"reference_k\":{},\"outward_heat_w\":{},\"dmean_dlog_htc\":{},\"dobjective_dlog_htc\":{}}}",
            quote(&state.region), num(state.area_m2)?, num(evaluated.htc[i])?, num(state.mean_wall_temperature_k)?,
            num(evaluated.coupled.reference_temperatures_k[i])?, num(state.heat_rate_w)?,
            optional(evaluated.gradient.as_ref().filter(|_| mean).map(|g| g.log_htc[i]))?,
            optional(evaluated.gradient.as_ref().map(|g| g.log_htc[i]))?));
    }
    let mut branches = Vec::new();
    for (hydraulic, thermal) in flow.branches.iter().zip(&evaluated.coupled.transport.branches) {
        branches.push(format!("{{\"name\":{},\"flow_m3_s\":{},\"inlet_k\":{},\"outlet_k\":{}}}",
            quote(&hydraulic.loss.name), num(hydraulic.flow.value())?, optional(thermal.inlet_temperature_k)?, optional(thermal.outlet_temperature_k)?));
    }
    let nodes = evaluated.coupled.transport.node_temperatures_k.iter().copied().map(optional).collect::<Result<Vec<_>>>()?.join(",");
    Ok(format!("{{\"schema\":{},\"authority\":\"nominal-estimate\",\"no_claim\":{},\"seed\":{},\"objective_region\":{},\"objective_mean_k\":{},\"coupling_iterations\":{},\"graph_sweeps\":{},\"source_w\":{},\"robin_out_w\":{},\"air_heat_imbalance_w\":{},\"walls\":[{}],\"branches\":[{}],\"node_temperatures_k\":[{}],\"solid_temperatures_k\":{},\"dmean_dinlet_k\":{},\"adjoint_residual\":{}}}\n",
        quote(RESULT_SCHEMA), quote(NO_CLAIM), quote(&request.seed.to_string()),
        request.objective.region().map_or_else(|| "null".into(), quote), optional(mean.then_some(evaluated.objective))?,
        evaluated.coupled.iterations, flow.sweeps, num(evaluated.source_total_w)?, num(evaluated.robin_total_w)?,
        num(evaluated.coupled.transport.heat_imbalance_w)?, walls.join(","), branches.join(","), nodes,
        numbers(&evaluated.temperatures)?, evaluated.gradient.as_ref().filter(|_| mean).map(|g| numbers(&g.inlets)).transpose()?.unwrap_or_else(|| "null".into()),
        optional(evaluated.gradient.as_ref().map(|g| g.interface_residual))?))
        .and_then(|result| {
            let prefix = result.strip_suffix("}\n").ok_or_else(|| bad("internal result framing mismatch"))?;
            Ok(format!("{prefix},\"solid_inputs\":{},\"objective\":{},\"dobjective_dinlet_k\":{},\"convection\":[{}],\"contacts\":{},\"contact_sensitivities\":{},\"coupling_solver\":{},\"gradient_scope\":\"steady thermal inlet, effective-coefficient and named contact-resistance sensitivities at fixed hydraulics, geometry, material laws and fluid properties; not fan-speed or channel-geometry derivatives; transient output has no adjoint\"}}\n",
                request.solid_data.render(request.conductivity, request.source)?,
                request.objective.render(&evaluated.objective_state, &request.mesh)?,
                evaluated.gradient.as_ref().map(|g| numbers(&g.inlets)).transpose()?.unwrap_or_else(|| "null".into()),
                evaluated.convection.iter().map(convection::Derived::render).collect::<Result<Vec<_>>>()?.join(","),
                request.contacts.as_ref().map(|contacts| contacts.render(&evaluated.contact_fluxes)).transpose()?.unwrap_or_else(|| "[]".into()),
                request.contacts.as_ref().map(|contacts| contacts.sensitivity_json(&evaluated.temperatures, evaluated.gradient.as_ref())).transpose()?.unwrap_or_else(|| "null".into()),
                acceleration::render(request.limits.relaxation, evaluated.gradient.is_some())?))
        })
}

fn execute(request: &Request, gate: &CancelGate) -> Result<String> {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(gate, arena, StreamKey { seed: request.seed, kernel_id: 717, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        poll(&cx)?;
        if let Some(schedule) = &request.transient {
            return transient::solve(request, &cx, schedule);
        }
        if let Some(design) = &request.fan_speed_design {
            return fan_speed::solve(request, &cx, design);
        }
        let flow = request.flow(&cx)?;
        let output = if let Some(design) = &request.design {
            let designed = design::solve(request, &cx, &flow, design)?;
            design::attach(render(request, &flow, &designed.passing)?, &designed)?
        } else {
            let coefficients = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let evaluated = request.evaluate(&cx, &flow, &coefficients, request.gradient)?;
            render(request, &flow, &evaluated)?
        };
        poll(&cx)?;
        match &request.fan {
            Some(fan) => fan.attach(output, &flow, fan.speed_ratio),
            None => Ok(output),
        }
    })
}

fn diagnostic(code: u8, failure: Failure, json_mode: bool) -> CommandOutput {
    let stderr = if json_mode {
        format!("{{\"schema\":\"frankensim.cooling-network.diagnostic.v1\",\"code\":{},\"message\":{}}}\n", quote(failure.code), quote(&failure.message))
    } else { format!("{failure}\n") };
    CommandOutput { exit_code: code, stdout: String::new(), stderr }
}

pub(super) fn run(args: &[OsString], json_mode: bool) -> CommandOutput {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return CommandOutput { exit_code: exit::SUCCESS, stdout: if json_mode {
            format!("{{\"schema\":{},\"help\":{}}}\n", quote(RESULT_SCHEMA), quote(HELP))
        } else { HELP.into() }, stderr: String::new() };
    }
    if args.len() != 1 { return diagnostic(exit::USAGE, bad(HELP), json_mode); }
    let read = (|| -> Result<String> {
        let mut text = String::new();
        File::open(std::path::Path::new(&args[0])).map_err(|e| bad(e.to_string()))?.take(MAX_INPUT_BYTES + 1)
            .read_to_string(&mut text).map_err(|e| bad(e.to_string()))?;
        if text.len() as u64 > MAX_INPUT_BYTES { return Err(bad("request exceeds 16 MiB")); }
        Ok(text)
    })();
    let text = match read { Ok(t) => t, Err(e) => return diagnostic(exit::INPUT, e, json_mode) };
    let request = match Request::parse(&text) {
        Ok(r) => r,
        Err(e) => {
            let class = if e.code == "cooling-network-transient-budget" { exit::BUDGET } else { exit::REFUSED };
            return diagnostic(class, e, json_mode);
        }
    };
    let gate = CancelGate::new();
    let duration = Duration::from_secs_f64(request.limits.wall_seconds);
    let started = Instant::now();
    let result = std::thread::scope(|scope| {
        let (stop, stopped) = std::sync::mpsc::channel::<()>();
        let gate_ref = &gate;
        scope.spawn(move || {
            if matches!(stopped.recv_timeout(duration), Err(std::sync::mpsc::RecvTimeoutError::Timeout)) { gate_ref.request(); }
        });
        let result = execute(&request, &gate);
        let _ = stop.send(());
        result
    });
    if started.elapsed() >= duration || gate.is_requested() {
        return diagnostic(exit::BUDGET, Failure { code: "cooling-network-time-budget", message: "wall-time budget exhausted; no partial result published".into() }, json_mode);
    }
    match result {
        Ok(stdout) => CommandOutput { exit_code: exit::SUCCESS, stdout, stderr: String::new() },
        Err(e) => {
            let class = if matches!(e.code, "cooling-network-design-budget" | "cooling-network-transient-budget") { exit::BUDGET } else { exit::REFUSED };
            diagnostic(class, e, json_mode)
        }
    }
}

#[cfg(test)]
mod tests;
