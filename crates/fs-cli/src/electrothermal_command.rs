//! `frankensim electrothermal`: canonical material packs -> coupled spatial fields.
//! Fixed-geometry isotropic conduction with optional prescribed fluid transport.

mod interfaces;
mod source;
mod radiation;
#[cfg(test)]
mod radiation_tests;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod transient_tests;
#[cfg(test)]
mod interface_tests;
#[cfg(test)]
mod thermal_interface_tests;

use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::Path;

use crate::{CommandOutput, exit};
use fs_exec::CancelGate;
use fs_lbm::conjugate::material_coupling::MaterialCouplingConfig;
use fs_lbm::conjugate::{
    ChtError, ContactResistance, ElectricSetup, Electrode, ElectrodeDrive,
    EnergyConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use crate::material_state_input::{JsonValue as J, quantity, state_point};
use source::{Problem, ResultField, Source};

const MAX_CELLS: usize = 65_536;
const MAX_MATERIALS: usize = 64;
const MAX_SCENE_BYTES: u64 = 4 * 1024 * 1024;
const MAX_PACK_BYTES: u64 = 4 * 1024 * 1024;
const MAX_TOTAL_PACK_BYTES: u64 = 64 * 1024 * 1024;
const NO_CLAIM: &str = "steady fixed-geometry all-solid isotropic ohmic/thermal conduction; source nominal values, not propagated uncertainty or experimental qualification; optional source-backed exterior gray radiation to black isothermal unit-view reservoirs; zero-thickness electrical interfaces use explicit heat partition, thermal interfaces use independently sourced resistance; no phase change, thermal deformation, fluid flow, microscopic contact-area model or AC effects";
const HELP: &str = "Usage: frankensim [--json] electrothermal <scene.json> [--transient <time.json>]\n\nLoad canonical FSMCDPK material-card packs and solve spatial current, Joule\nheating and temperature together. Both resistivity and conductivity are\nre-resolved at every changed local temperature; unsupported states refuse.\nMaterial substitutions change pack data, not solver code. Optional\nelectrical_interfaces and thermal_interfaces bind canonical FSINTPK ordered\ninterface cards to voxel faces. Electrical heating and thermal resistance\nremain independent responses, evaluated at the same local temperatures.\n\nAn optional fluid object assigns one material to prescribed conservative flow,\nwith sourced frozen reference transport. See docs/MATERIAL_ELECTROTHERMAL_FLUID.md.\n\n--transient additionally requires sourced density and heat capacity over an\nexplicit temperature interval, and advances stored energy with scheduled\nelectrode drives and authored heat sources. Cell mass remains fixed.\n\nScene schema: frankensim.electrothermal.v1. Time request schema:\nfrankensim.electrothermal.transient.v1. See docs/MATERIAL_ELECTROTHERMAL.md,\ndocs/MATERIAL_ELECTROTHERMAL_TRANSIENT.md, docs/ELECTRICAL_CONTACTS.md and\ndocs/THERMAL_FACE_CONTACTS.md. For source-backed exterior radiation see docs/ELECTROTHERMAL_RADIATION.md.\nAll lengths are metres, temperatures kelvin,\ncurrent amperes and power watts. --json includes fields, source identities\nand coupling/energy defects.\n";

pub fn run(args: &[OsString], json: bool) -> CommandOutput {
    if args.len() == 1 && matches!(args[0].to_str(), Some("--help" | "-h")) {
        return CommandOutput { exit_code: exit::SUCCESS, stdout: HELP.into(), stderr: String::new() };
    }
    let transient = args.len() == 3 && args[1] == "--transient";
    if args.len() != 1 && !transient {
        return refusal(exit::USAGE, "electrothermal-usage", HELP, json);
    }
    let mut problem = match load(Path::new(&args[0])) {
        Ok(problem) => problem,
        Err(error) => return refusal(exit::INPUT, "electrothermal-input", &error, json),
    };
    if transient {
        let request = match source::transient::load(Path::new(&args[2]), problem.domain.cell_count()) {
            Ok(request) => request,
            Err(error) => return refusal(exit::INPUT, "electrothermal-transient-input", &error, json),
        };
        return match source::transient::march(&problem, &request, &CancelGate::new()) {
            Ok(result) => CommandOutput {
                exit_code: exit::SUCCESS, stdout: source::transient::render(&problem, &result, json), stderr: String::new(),
            },
            Err(error) => {
                let code = if matches!(error, ChtError::Cancelled) { exit::CANCELLED } else { exit::REFUSED };
                refusal(code, "electrothermal-transient-solve", &error.to_string(), json)
            }
        };
    }
    let result = match problem.solve(&CancelGate::new()) {
        Ok(result) => result,
        Err(error) => {
            let code = if matches!(error, ChtError::Cancelled) { exit::CANCELLED } else { exit::REFUSED };
            return refusal(code, "electrothermal-solve", &error.to_string(), json);
        }
    };
    let balance = &result.coupled.thermal.report.balance;
    if [result.coupled.electric.total_joule_w, result.coupled.electric.delivered_w,
        balance.residual_w, balance.relative_residual, balance.boundary_outflow_w]
        .iter().any(|v| !v.is_finite())
    {
        return refusal(exit::REFUSED, "electrothermal-result", "nonfinite aggregate result", json);
    }
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout: render(&problem, &result, json),
        stderr: String::new(),
    }
}

fn read(path: &Path, cap: u64) -> Result<Vec<u8>, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let metadata = file.metadata().map_err(|e| e.to_string())?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err(format!("{}: expected a regular file of at most {cap} bytes", path.display()));
    }
    let mut bytes = Vec::new();
    file.take(cap + 1).read_to_end(&mut bytes).map_err(|e| e.to_string())?;
    if bytes.len() as u64 > cap { return Err("input grew beyond its byte limit".into()); }
    Ok(bytes)
}

fn fields(value: &J, allowed: &[&str]) -> Result<(), String> {
    let object = value.as_object().ok_or("expected an object")?;
    if let Some((key, _)) = object.iter().find(|(key, _)| !allowed.contains(&key.as_str())) {
        return Err(format!("unsupported field {key:?}"));
    }
    Ok(())
}
fn required<'a>(value: &'a J, key: &str) -> Result<&'a J, String> {
    value.get(key).ok_or_else(|| format!("{key} is required"))
}
fn string<'a>(value: &'a J, key: &str) -> Result<&'a str, String> {
    required(value, key)?.as_str().filter(|s| !s.trim().is_empty())
        .ok_or_else(|| format!("{key} must be a nonempty string"))
}
fn number(value: &J) -> Result<f64, String> {
    value.as_f64().filter(|x| x.is_finite()).ok_or_else(|| "expected a finite number".into())
}
fn scalar(value: &J, key: &str) -> Result<f64, String> { number(required(value, key)?) }
fn positive(value: &J, key: &str) -> Result<f64, String> {
    let v = scalar(value, key)?;
    if v <= 0.0 { Err(format!("{key} must be positive")) } else { Ok(v) }
}
fn integer(value: &J) -> Result<usize, String> {
    value.number_raw().and_then(|s| s.parse().ok()).ok_or_else(|| "expected an unsigned integer".into())
}
fn array<'a>(value: &'a J, key: &str, cap: usize) -> Result<&'a [J], String> {
    let items = required(value, key)?.as_array().ok_or_else(|| format!("{key} must be an array"))?;
    if items.len() > cap { Err(format!("{key} exceeds {cap} entries")) } else { Ok(items) }
}
fn triple(value: &J) -> Result<[usize; 3], String> {
    let items = value.as_array().filter(|a| a.len() == 3).ok_or("expected three integer coordinates")?;
    Ok([integer(&items[0])?, integer(&items[1])?, integer(&items[2])?])
}
fn budget(value: &J, key: &str, cap: usize) -> Result<usize, String> {
    let n = integer(required(value, key)?)?;
    if n == 0 || n > cap { Err(format!("{key} must be in 1..={cap}")) } else { Ok(n) }
}

fn load(path: &Path) -> Result<Problem, String> {
    let bytes = read(path, MAX_SCENE_BYTES)?;
    let root = J::parse(std::str::from_utf8(&bytes).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    fields(&root, &["schema", "cells", "cell_size_m", "materials", "material_indices",
        "initial_temperature_k", "thermal_faces", "power_w", "contacts", "electrodes", "numerics",
        "fluid", "electrical_interfaces", "thermal_interfaces", "radiation"])?;
    if string(&root, "schema")? != "frankensim.electrothermal.v1" { return Err("unsupported electrothermal schema".into()); }
    let [nx, ny, nz] = triple(required(&root, "cells")?)?;
    let cells = nx.checked_mul(ny).and_then(|n| n.checked_mul(nz))
        .filter(|n| *n > 0 && *n <= MAX_CELLS).ok_or("domain exceeds its positive cell budget")?;
    let dx = positive(&root, "cell_size_m")?;
    let mut domain = VoxelDomain::new(nx, ny, nz, dx).map_err(|e| e.to_string())?;
    let definitions = array(&root, "materials", MAX_MATERIALS)?;
    if definitions.is_empty() { return Err("materials must not be empty".into()); }
    let mut materials = Vec::with_capacity(definitions.len());
    let mut total_bytes = 0u64;
    for definition in definitions {
        fields(definition, &["pack", "conducting", "temperature_axis", "temperature_kind", "state"])?;
        let pack_path = path.parent().unwrap_or(Path::new(".")).join(string(definition, "pack")?);
        let pack_bytes = read(&pack_path, MAX_PACK_BYTES.min(MAX_TOTAL_PACK_BYTES - total_bytes))?;
        total_bytes += pack_bytes.len() as u64;
        let pack = NormalizedMaterialCardPack::from_bytes(&pack_bytes).map_err(|e| e.to_string())?;
        let conducting = match required(definition, "conducting")? {
            J::Bool(value) => *value,
            _ => return Err("conducting must explicitly be true or false".into()),
        };
        let point = state_point(definition)?;
        materials.push(Source::new(pack, point, string(definition, "temperature_axis")?.into(),
            quantity(string(definition, "temperature_kind")?, Temperature::DIMS)?, conducting)
            .map_err(|e| e.to_string())?);
    }
    let indices = array(&root, "material_indices", MAX_CELLS)?;
    if indices.len() != cells { return Err("material_indices needs exactly one assignment per cell".into()); }
    let mut material_indices = Vec::with_capacity(cells);
    for (cell, value) in indices.iter().enumerate() {
        let m = integer(value)?;
        if m >= materials.len() { return Err(format!("cell {cell} names an unknown material")); }
        let [x, y, z] = domain.coords(cell);
        domain.set(x, y, z, Voxel::Solid(m as u16));
        material_indices.push(m);
    }
    let fluid = source::load_fluid(&root, &mut domain, &material_indices, &materials)?;
    let electrical_interfaces = interfaces::load(
        &root, path, &domain, &materials, &material_indices, &mut total_bytes,
    )?;
    let thermal_interfaces = interfaces::load_thermal(
        &root, path, &domain, &materials, &material_indices, &mut total_bytes,
    )?;
    let initial_value = required(&root, "initial_temperature_k")?;
    let initial_temperature = if let Some(values) = initial_value.as_array() {
        if values.len() != cells { return Err("initial temperature field has wrong length".into()); }
        values.iter().map(number).collect::<Result<Vec<_>, _>>()?
    } else { vec![number(initial_value)?; cells] };
    if initial_temperature.iter().any(|t| *t <= 0.0) { return Err("initial temperatures must be positive kelvin".into()); }
    let face_inputs = array(&root, "thermal_faces", 6)?;
    if face_inputs.len() != 6 { return Err("declare all six thermal faces in x-,x+,y-,y+,z-,z+ order".into()); }
    let mut faces = [ThermalFace::Adiabatic; 6];
    for (face, value) in faces.iter_mut().zip(face_inputs) {
        *face = if value.as_str() == Some("adiabatic") { ThermalFace::Adiabatic }
        else if value.get("temperature_k").is_some() {
            fields(value, &["temperature_k"])?;
            ThermalFace::Temperature(positive(value, "temperature_k")?)
        } else if value.get("inflow_temperature_k").is_some() {
            fields(value, &["inflow_temperature_k"])?;
            ThermalFace::Inflow { temperature: positive(value, "inflow_temperature_k")? }
        } else if value.get("backflow_temperature_k").is_some() {
            fields(value, &["backflow_temperature_k"])?;
            ThermalFace::Outflow { backflow_temperature: positive(value, "backflow_temperature_k")? }
        } else if value.get("heat_flux_w_m2").is_some() {
            fields(value, &["heat_flux_w_m2"])?;
            ThermalFace::HeatFlux(scalar(value, "heat_flux_w_m2")?)
        } else {
            fields(value, &["h_w_m2_k", "ambient_k"])?;
            ThermalFace::Convective { h: positive(value, "h_w_m2_k")?, ambient: positive(value, "ambient_k")? }
        };
    }
    let mut thermal = ThermalSetup::new(faces);
    if root.get("power_w").is_some() {
        thermal.power_w = array(&root, "power_w", MAX_CELLS)?.iter().map(number).collect::<Result<Vec<_>, _>>()?;
        if thermal.power_w.len() != cells { return Err("power_w needs one value per cell when supplied".into()); }
    }
    if root.get("contacts").is_some() {
        let mut pairs = std::collections::BTreeSet::new();
        for contact in array(&root, "contacts", 4096)? {
            fields(contact, &["materials", "resistance_m2_k_w"])?;
            let pair = array(contact, "materials", 2)?;
            if pair.len() != 2 { return Err("contact needs two material indices".into()); }
            let (a, b) = (integer(&pair[0])?, integer(&pair[1])?);
            if a == b || a.max(b) >= materials.len() || !pairs.insert((a.min(b), a.max(b))) {
                return Err("contact pair must be distinct, declared and unique".into());
            }
            if fluid.as_ref().is_some_and(|f| f.material_index == a || f.material_index == b) {
                return Err("solid contact resistance cannot be assigned to the fluid material".into());
            }
            let resistance = scalar(contact, "resistance_m2_k_w")?;
            if resistance < 0.0 { return Err("contact resistance cannot be negative".into()); }
            thermal.contacts.push(ContactResistance { materials: (a as u16, b as u16), resistance_m2_k_w: resistance });
        }
    }
    let mut electrodes = Vec::new();
    for electrode in array(&root, "electrodes", 256)? {
        fields(electrode, &["lo", "hi", "voltage_v", "current_a"])?;
        let drive = match (electrode.get("voltage_v"), electrode.get("current_a")) {
            (Some(v), None) => ElectrodeDrive::Voltage(number(v)?),
            (None, Some(i)) => ElectrodeDrive::Current(number(i)?),
            _ => return Err("each electrode needs exactly one voltage_v or current_a".into()),
        };
        electrodes.push(Electrode { lo: triple(required(electrode, "lo")?)?, hi: triple(required(electrode, "hi")?)?, drive });
    }
    let n = required(&root, "numerics")?;
    fields(n, &["max_passes", "max_thermal_passes", "temperature_tolerance_k", "power_tolerance_w",
        "relative_tolerance", "conductivity_relative_tolerance", "energy_tolerance", "energy_max_iterations",
        "thermal_contact_relative_tolerance"])?;
    let coupling = MaterialCouplingConfig {
        max_passes: budget(n, "max_passes", 1000)?,
        temperature_tolerance_k: positive(n, "temperature_tolerance_k")?,
        power_tolerance_w: scalar(n, "power_tolerance_w")?,
        relative_tolerance: scalar(n, "relative_tolerance")?,
    };
    let conductivity_tolerance = positive(n, "conductivity_relative_tolerance")?;
    let thermal_contact_tolerance = if !thermal_interfaces.is_empty()
        || n.get("thermal_contact_relative_tolerance").is_some()
    {
        positive(n, "thermal_contact_relative_tolerance")?
    } else {
        // No thermal interface is present, so this tolerance is not consumed.
        conductivity_tolerance
    };
    let energy = EnergyConfig {
        tolerance: positive(n, "energy_tolerance")?,
        max_iterations: budget(n, "energy_max_iterations", 1_000_000)?,
        ..EnergyConfig::default()
    };
    if coupling.power_tolerance_w < 0.0 || !(0.0..1.0).contains(&coupling.relative_tolerance)
        || conductivity_tolerance >= 1.0 || energy.tolerance >= 1.0 || thermal_contact_tolerance >= 1.0
    { return Err("invalid coupling or solver tolerance".into()); }
    let radiation = radiation::load(&root, &domain, &materials, &material_indices, &thermal)?;
    Ok(Problem { domain, material_indices, materials, fluid, electrical_interfaces, radiation,
        thermal_interfaces, thermal_contact_tolerance,
        electric: ElectricSetup { conductors: Vec::new(), electrodes }, thermal, initial_temperature,
        coupling, energy, max_thermal_passes: budget(n, "max_thermal_passes", 1000)?, conductivity_tolerance })
}

fn quote(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""), '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"), '\r' => out.push_str("\\r"), '\t' => out.push_str("\\t"),
            c if c < '\u{20}' => { let _ = write!(out, "\\u{:04x}", c as u32); }
            c => out.push(c),
        }
    }
    out.push('"'); out
}
fn values(values: impl IntoIterator<Item = f64>) -> String {
    format!("[{}]", values.into_iter().map(|v| {
        if v.is_finite() { v.to_string() } else { "null".into() }
    }).collect::<Vec<_>>().join(","))
}
fn refusal(exit_code: u8, code: &str, message: &str, json: bool) -> CommandOutput {
    CommandOutput { exit_code, stdout: String::new(), stderr: if json {
        format!("{{\"schema\":\"frankensim.electrothermal.diagnostic.v1\",\"code\":{},\"message\":{}}}\n", quote(code), quote(message))
    } else { format!("{code}: {message}\n") } }
}
fn render(problem: &Problem, result: &ResultField, json: bool) -> String {
    let r = &result.coupled;
    let no_claim = if problem.fluid.is_some() { source::FLUID_NO_CLAIM } else { NO_CLAIM };
    if !json {
        return format!("Electrothermal solve: {} cells, {} coupling passes, {} thermal solves\nPeak temperature: {} K\nJoule power: {} W; electrode power: {} W\nTemperature defect: {} K; Joule-source defect: {} W/cell; conductivity relative defect: {}; thermal-contact relative defect: {}\n{}\n",
            problem.domain.cell_count(), r.report.passes, result.thermal_solves,
            r.thermal.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            r.electric.total_joule_w, r.electric.delivered_w, r.report.temperature_change_k,
            r.report.joule_source_change_w, result.conductivity_relative_defect,
            result.thermal_contact_relative_defect, no_claim);
    }
    let materials = problem.materials.iter().map(|m| format!(
        "{{\"pack\":{},\"pack_identity\":{},\"card_identity\":{}}}",
        quote(m.pack.pack_id()), quote(&m.pack.content_hash().to_hex()), quote(&m.pack.card().content_hash().to_hex()),
    )).collect::<Vec<_>>().join(",");
    let sources = result.sources.iter().enumerate().map(|(cell, state)| format!(
        "{{\"material_index\":{},\"thermal_state\":{},\"electrical_state\":{}}}",
        problem.material_indices[cell], quote(&state.thermal.identity().to_hex()),
        state.conductor.as_ref().map_or_else(|| "null".into(), |c| quote(&c.material().identity().to_hex())),
    )).collect::<Vec<_>>().join(",");
    let electrodes = r.electric.electrodes.iter().map(|e| format!(
        "{{\"voltage_v\":{},\"current_a\":{}}}",
        if e.voltage_v.is_finite() { e.voltage_v.to_string() } else { "null".into() },
        if e.current_a.is_finite() { e.current_a.to_string() } else { "null".into() },
    )).collect::<Vec<_>>().join(",");
    format!(concat!(
        "{{\"schema\":\"frankensim.electrothermal.result.v1\",\"authority\":\"estimated-numerical\",",
        "\"cells\":{:?},\"cell_size_m\":{},\"temperature_k\":{},\"potential_v\":{},\"joule_w\":{},",
        "\"solved_conductivity_w_m_k\":{},\"resolved_conductivity_w_m_k\":{},",
        "\"total_joule_w\":{},\"electrode_power_w\":{},\"thermal_balance_residual_w\":{},",
        "\"thermal_balance_relative_residual\":{},\"passes\":{},\"thermal_solves\":{},",
        "\"temperature_defect_k\":{},\"joule_source_defect_w\":{},\"conductivity_relative_defect\":{},",
        "\"materials\":[{}],\"source_states\":[{}],\"electrodes\":[{}],\"fluid_reference\":{},",
        "\"electrical_interfaces\":{},\"thermal_interfaces\":{},",
        "\"thermal_contact_relative_defect\":{},\"radiation\":{},\"no_claim\":{}}}\n"),
        problem.domain.dims(), problem.domain.dx(), values(r.thermal.temperature.iter().copied()),
        values(r.electric.potential_v.iter().copied()), values(r.electric.joule_w.iter().copied()),
        values(r.thermal.conductivity.iter().map(|k| k[0])), values(result.sources.iter().map(|s| s.k)),
        r.electric.total_joule_w, r.electric.delivered_w, r.thermal.report.balance.residual_w,
        r.thermal.report.balance.relative_residual, r.report.passes, result.thermal_solves,
        r.report.temperature_change_k, r.report.joule_source_change_w, result.conductivity_relative_defect,
        materials, sources, electrodes, problem.fluid.as_ref().map_or_else(|| "null".into(), |f| f.render()),
        interfaces::render(problem, &result.interfaces), interfaces::render_thermal(problem, result),
        result.thermal_contact_relative_defect, radiation::render(&result.radiation), quote(no_claim))
}
