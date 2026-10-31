//! File-driven wind sources through the shared pressure and WAV owners.
#[path = "wind_input/attachments.rs"]
mod attachments;
#[path = "wind_input/gas_sources.rs"]
mod gas_sources;
#[path = "wind_input/finite_gas_heat.rs"]
mod finite_gas_heat;
use std::io::{Read, Write};
use std::path::Path;
use fs_blake3::ContentHash;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, ObservationPolicy, PressureRenderer};
use fs_couple::render::schedule::reed::{ReedPerformance, REED_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::performance::{ApertureObservation, AperturePerformance, CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::{PlateValvePerformance, MAX_PLATE_VALVE_PERFORMANCE_BYTES, PLATE_VALVE_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::bernoulli_aperture::radiation::environment::multiport::input::MAX_MIMO_RADIATION_BYTES;
use fs_couple::flue_loop::network::performance::file::{
    FlueFilePerformance, FLUE_PERFORMANCE_SCHEMA, MAX_FLUE_WALL_PACK_BYTES,
};
use fs_matdb::NormalizedMaterialCardPack;
use fs_exec::CancelGate;
use fs_material::gas_species::input::ThermochemicalSources;
use super::{RATE, create_outputs, json_string, stream_output};

fn opening_loss_provenance<'a>(
    bindings: impl Iterator<Item = (&'static str,
        &'a fs_couple::bernoulli_aperture::opening_loss::BoundOpeningLoss)>,
) -> String {
    let records = bindings.map(|(side, binding)| {
        let input = binding.input();
        let gas = binding.gas();
        let law = binding.loss();
        format!("{{\"side\":\"{side}\",\"node\":{},\"radius_m\":{:e},\"area_ratio\":{:e},\"maximum_contracted_mach\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"coefficient_pa_s2_m6\":{:e},\"maximum_flow_m3_s\":{:e},\"source\":{}}}",
            input.node, binding.radius_m(), input.contraction.area_ratio,
            input.contraction.maximum_contracted_mach, gas.density, gas.sound_speed,
            law.coefficient_pa_s2_m6(), law.maximum_flow_m3_s(),
            json_string(&input.contraction.provenance))
    }).collect::<Vec<_>>().join(",");
    format!(",\"opening_losses\":[{records}],\"opening_loss_scope\":\"same-sample symmetric quasi-steady separation in series with original local radiation; no vortex memory, directional pressure recovery or measured instrument coefficient\"")
}

// The aperture-flow path, and whether `--labium` adds the labium displacement to it.
type FlowsOption<'a> = Option<(&'a str, bool)>;

fn options(args: &[String]) -> Result<(&str, &str, usize, bool, ObservationPolicy, FlowsOption<'_>, bool, Vec<&str>), String> {
    let mut paths = Vec::new();
    let mut flows = None;
    let mut labium = false;
    let mut block = None;
    let mut decimate = false;
    let mut self_contained = false;
    let mut attached_sources = Vec::new();
    let mut policy = ObservationPolicy::Diagnostic;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--self-contained" => {
                if self_contained { return Err("--self-contained may only be supplied once".into()); }
                self_contained = true;
            }
            "--attached-source" => {
                let path = iter.next().map(String::as_str).filter(|value| !value.starts_with('-'))
                    .ok_or_else(|| "--attached-source needs a relative source path".to_string())?;
                if attached_sources.len() == attachments::MAX_FILES {
                    return Err("--attached-source allows at most 32 files".into());
                }
                attached_sources.push(path);
            }
            "--acceptance" => {
                if policy == ObservationPolicy::Acceptance { return Err("--acceptance may only be supplied once".into()); }
                policy = ObservationPolicy::Acceptance;
            }
            "--aperture-flows" => {
                if flows.is_some() { return Err("--aperture-flows may only be supplied once".into()); }
                flows = Some(iter.next().map(String::as_str).filter(|v| !v.starts_with('-'))
                    .ok_or_else(|| "--aperture-flows needs an output path".to_string())?);
            }
            "--labium" => {
                if labium { return Err("--labium may only be supplied once".into()); }
                labium = true;
            }
            "--decimate" => {
                if decimate { return Err("--decimate may only be supplied once".into()); }
                decimate = true;
            }
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                block = Some(iter.next().and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?);
            }
            value if value.starts_with('-') => return Err(format!(
                "unsupported wind option {value:?}; geometry, pressure, duration and full-scale belong to the source file"
            )),
            value => paths.push(value),
        }
    }
    if !attached_sources.is_empty() && !self_contained {
        return Err("--attached-source requires --self-contained".into());
    }
    if labium && flows.is_none() {
        return Err("--labium adds the labium jet displacement to --aperture-flows".into());
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render wind INPUT.performance OUT.wav [--block N] [--decimate] [--acceptance] [--aperture-flows FLOWS [--labium]] [--self-contained [--attached-source RELATIVE_PATH]...]".into());
    };
    Ok((*input, *output, block.unwrap_or(512), decimate, policy, flows.map(|f| (f, labium)), self_contained, attached_sources))
}

// A physical propagation delay is source geometry, not an output/filter latency
// to be removed by ensemble alignment. Keep it on the mechanical sample clock.
pub(super) fn outlet_provenance(p: &AperturePerformance, radiation: Option<BaffledRadiationLoad>) -> String {
    let load_json=radiation_provenance(radiation);
    let (c,radius,_density,_speed)=match (p.observation(),p.system()) {
        (ApertureObservation::TubeBaffled(receiver), CoupledAperture::Dg(runtime)) => {
            return super::duct_metadata::dg_baffled_observation_provenance(
                runtime.load(), receiver, p.baffled_receiver().expect("admitted DG microphone"),
            );
        }
        (ApertureObservation::TubeBaffled(c),CoupledAperture::Tube(t))=>(c,t.spec().radius_m,t.aperture().spec().density_kg_m3,t.spec().sound_speed_m_s),
        (ApertureObservation::NetworkBaffled {node,receiver},CoupledAperture::Network(n))=> {
            let (index,section)=n.spec().sections.iter().enumerate().find(|(_,s)|s.nodes.contains(&node)).expect("admitted physical outlet section");
            let (density,speed)=n.section_medium(index).expect("admitted outlet medium");
            (receiver,section.radius_m,density,speed)
        }
        _=>return if radiation.is_some() {
            format!(",\"observation_scope\":\"internal coupled tube pressure; not an exterior microphone\"{load_json}")
        } else {String::new()},
    };
    let medium=p.baffled_medium().expect("admitted receiver medium");
    let (density,speed)=(medium.density,medium.sound_speed);
    let mic=p.baffled_receiver().expect("admitted receiver");
    let scope=if radiation.is_some() {
        "exterior Rayleigh pressure from radiation-loaded terminal flow; compact passive load and uniform circular outlet in infinite rigid baffle"
    } else {
        "one-way exterior Rayleigh pressure; uniform circular outlet in infinite rigid baffle; terminal load is independently supplied"
    };
    format!(",\"observation_scope\":\"{scope}\",\"outlet_receiver\":{{\"position_m\":[{:e},{:e},{:e}],\"radius_m\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"propagation_delay_mechanical_samples\":[{},{}],\"radiation_feedback_added\":{}}}{load_json}",
        c.position_m[0],c.position_m[1],c.position_m[2],radius,
        density,speed,
        c.radial_rings,c.angular_points,c.maximum_frequency_hz,mic.delay_samples.0,mic.delay_samples.1,radiation.is_some())
}

fn radiation_provenance(radiation: Option<BaffledRadiationLoad>) -> String {
    if let Some(load)=radiation {
        let spec=load.load();let term=spec.terms()[0];
        format!(",\"radiation_load\":{{\"model\":\"compact-baffled-piston-positive-real-v1\",\"maximum_frequency_hz\":{:e},\"resistance_pa_s_m3\":{:e},\"pole_rate_per_s\":{:e},\"checked_max_complex_relative_error\":{:e},\"checked_max_resistance_relative_error\":{:e},\"replaces_memoryless_reflection\":true,\"extra_end_correction_m\":0,\"scope\":\"finite-sample low-ka analogue and bilinear checks, not broadband or measured calibration\"}}",
            load.maximum_frequency_hz(),term.resistance_pa_s_m3,term.rate_per_s,
            load.max_complex_relative_error(),load.max_resistance_relative_error())
    } else {String::new()}
}

// The flue jet's declared speed edge (ruling 8132): `jet_limits` MIN_SPEED and its state below.
fn jet_speed_edge_provenance(voice: &fs_couple::flue_loop::FlueVoice) -> String {
    format!(",\"jet_speed_edge\":{{\"minimum_jet_speed_m_s\":{:e},\"scope\":\"declared validity edge, not a stability limit; below it the onset reads the rest state before the first motion, and a lookup beyond the retained history carries no jet deflection, so the windway flow splits at zero deflection\"}}",
        voice.minimum_jet_speed_m_s())
}

fn flue_wall_provenance(p: &FlueFilePerformance) -> String {
    use fs_couple::bernoulli_aperture::wall::material::AxialConstraint;
    let walls = p.renderer().voice().network().expect("file admits graph").material_walls();
    let records = p.wall_records().iter().map(|record| {
        let wall = walls.iter().find(|wall| wall.source().node == record.node())
            .expect("retained material source has its original physical wall");
        let assignment = wall.source();
        let material = &assignment.material;
        let (model, directions) = graph_wall_constitutive_metadata(assignment);
        let state = material.material();
        let coordinate = |name: &str| material.query_point().iter()
            .find(|(axis, _)| axis == name).expect("file supplies T and omega").1;
        let axial = match assignment.axial {
            AxialConstraint::StressFree => "stress-free",
            AxialConstraint::StrainFixed => "strain-fixed",
        };
        let (loss, relaxation) = graph_wall_loss_metadata(assignment);
        let patch = wall.patch();
        format!("{{\"node\":{},\"model\":\"{model}\",\"pack_path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"resolved_state_hash\":\"{}\",\"source_artifact\":\"{}\",\"material_chemistry\":{},\"material_phase\":{},\"material_process\":{},\"material_revision\":{},\"temperature_k\":{:e},\"angular_frequency_rad_s\":{:e},\"band_hz\":[{:e},{:e}],\"radius_m\":{:e},\"covered_length_m\":{:e},\"area_m2\":{:e},\"thickness_m\":{:e},\"axial_constraint\":\"{axial}\",\"loss\":\"{loss}\",\"surface_density_kg_m2\":{:e},\"stiffness_pa_m\":{:e},\"resistance_pa_s_m\":{:e}{relaxation}{directions}}}",
            record.node(),json_string(record.path()),record.pack().content_hash().to_hex(),
            material.card_identity().to_hex(),material.identity().to_hex(),record.pack().source_artifact().to_hex(),
            json_string(&state.chemistry),json_string(&state.phase),json_string(&state.process),state.revision,
            coordinate("T"),coordinate("omega"),wall.band_hz()[0],wall.band_hz()[1],wall.radius_m(),
            wall.length_m(),patch.area_m2,assignment.thickness_m,
            patch.wall.surface_density,patch.wall.stiffness_per_area,patch.wall.resistance)
    }).collect::<Vec<_>>().join(",");
    let scope = if walls.is_empty() { "rigid" } else {
        "original interior dual strips; unselected and endpoint half-strips rigid; fixed-reference small radial motion; no axial or bending shell coupling or exterior wall radiation"
    };
    format!(",\"material_walls\":[{records}],\"wall_strips\":{},\"wall_model_scope\":\"{scope}\"",walls.len())
}

fn wall_constitutive_metadata(
    material: &fs_material::state_point::ResolvedMaterialStatePoint,
    axes: Option<[usize; 2]>, time_property: Option<&str>,
    normal_properties: Option<[&str; 3]>,
) -> (&'static str, String) {
    match axes {
        None => ("isotropic-local-ring-membrane", String::new()),
        Some([hoop, axial]) => {
            let mut viscosity = time_property.map_or_else(String::new, |name| {
                let time = material.property(name).expect("admitted original viscosity time").value_si();
                format!(",\"proportional_viscosity\":{{\"time_property\":{},\"time_s\":{time:e},\"law\":\"C_viscous=tau*C_elastic\"}}", json_string(name))
            });
            if let Some(names) = normal_properties {
                let entries = ["hh", "ha", "aa"].into_iter().zip(names).map(|(axis, name)| {
                    let value = material.property(name).expect("admitted original normal viscosity").value_si();
                    format!("\"{axis}\":{{\"property\":{},\"value_pa_s\":{value:e}}}", json_string(name))
                }).collect::<Vec<_>>().join(",");
                viscosity.push_str(&format!(",\"normal_viscosity\":{{{entries},\"law\":\"symmetric normal-plane Kelvin-Voigt; exact constrained C+sD reduction\"}}"));
            }
            ("orthotropic-principal-local-ring-membrane", format!(
                ",\"material_axes\":{{\"hoop\":{},\"axial\":{},\"basis\":\"one-based principal axes following local cylindrical directions\"}}{viscosity}", hoop + 1, axial + 1))
        }
    }
}

fn graph_wall_constitutive_metadata(
    source: &fs_couple::bernoulli_aperture::wall::material::MaterialWallAssignment,
) -> (&'static str, String) {
    use fs_couple::bernoulli_aperture::wall::material::WallMaterialModel;
    match &source.model {
        WallMaterialModel::Isotropic => wall_constitutive_metadata(&source.material, None, None, None),
        WallMaterialModel::Orthotropic { axes, proportional_viscosity_time_property } =>
            wall_constitutive_metadata(&source.material, Some(*axes), proportional_viscosity_time_property.as_deref(), None),
        WallMaterialModel::OrthotropicNormal { axes, properties, .. } =>
            wall_constitutive_metadata(&source.material, Some(*axes), None, Some(properties.each_ref().map(String::as_str))),
    }
}

fn graph_wall_loss_metadata(
    source: &fs_couple::bernoulli_aperture::wall::material::MaterialWallAssignment,
) -> (&'static str, String) {
    use fs_couple::bernoulli_aperture::wall::material::{WallMaterialLoss, WallMaterialModel};
    match &source.loss {
        WallMaterialLoss::Elastic => ("elastic", String::new()),
        WallMaterialLoss::KelvinVoigt => (match &source.model {
            WallMaterialModel::Isotropic => "kelvin-voigt",
            WallMaterialModel::Orthotropic { .. } => "proportional-kelvin-voigt",
            WallMaterialModel::OrthotropicNormal { .. } => "normal-kelvin-voigt",
        }, String::new()),
        WallMaterialLoss::Maxwell(spectrum) => {
            let material = spectrum.material();
            let terms = material.terms.iter().map(|(modulus, time)| format!(
                "{{\"modulus_pa\":{modulus:e},\"relaxation_time_s\":{time:e}}}"))
                .collect::<Vec<_>>().join(",");
            ("maxwell", format!(",\"equilibrium_modulus_pa\":{:e},\"maxwell_terms\":[{terms}]", material.e_inf))
        }
    }
}

fn dg_wall_source_metadata(
    source: &fs_couple::dg_port::input::wall::DgMaterialWallRecord,
) -> String {
    use fs_couple::bernoulli_aperture::wall::material::{AxialConstraint, WallMaterialLoss};
    let wall = source.wall();
    let (model, directions) = wall_constitutive_metadata(wall.material(),
        wall.orthotropic_axes(), wall.proportional_viscosity_time_property(), source.normal_viscosity_properties());
    let axial = match wall.axial_constraint() {
        AxialConstraint::StressFree => "stress-free",
        AxialConstraint::StrainFixed => "strain-fixed",
    };
    let (loss, relaxation) = match wall.loss() {
        WallMaterialLoss::Elastic => ("elastic", String::new()),
        WallMaterialLoss::KelvinVoigt => (
            if wall.proportional_viscosity_time_property().is_some() { "proportional-kelvin-voigt" }
                else if source.normal_viscosity_properties().is_some() { "normal-kelvin-voigt" }
                else { "kelvin-voigt" }, String::new()),
        WallMaterialLoss::Maxwell(spectrum) => {
            let terms = spectrum.material().terms.iter().map(|(modulus, time)| format!(
                "{{\"modulus_pa\":{modulus:e},\"relaxation_time_s\":{time:e}}}"))
                .collect::<Vec<_>>().join(",");
            ("maxwell", format!(",\"maxwell_terms\":[{terms}]"))
        }
    };
    format!("{{\"model\":\"{model}\",\"pack_path\":{},\"pack_hash\":\"{}\",\"resolved_state_hash\":\"{}\",\"thickness_m\":{:e},\"axial_constraint\":\"{axial}\",\"band_hz\":[{:e},{:e}],\"loss\":\"{loss}\",\"propagator\":\"original-exact-dg\"{relaxation}{directions}}}",
        json_string(source.path()), source.pack().content_hash().to_hex(),
        wall.material().identity().to_hex(), wall.thickness_m(), wall.band_hz()[0], wall.band_hz()[1])
}

fn dg_wall_provenance(
    layout: &fs_couple::dg_port::input::wall::regions::DgMaterialWallLayout,
) -> String {
    fn source(layout: &fs_couple::dg_port::input::wall::regions::DgMaterialWallLayout) -> String {
        if let Some(source) = layout.uniform_record() {
        dg_wall_source_metadata(source)
    } else if layout.regions().is_empty() {
        "null".into()
    } else {
        let regions = layout.regions().iter().map(|region| {
            let interval = region.interval_m();
            format!("{{\"interval_m\":[{:e},{:e}],\"source\":{}}}",
                interval[0], interval[1], dg_wall_source_metadata(region.source()))
        }).collect::<Vec<_>>().join(",");
        let band = layout.band_hz().expect("resolved nonempty regional wall");
        let model = if layout.regions().iter().any(|region| region.source().wall().orthotropic_axes().is_some()) {
            "regional-local-ring-membranes"
        } else { "regional-isotropic-local-ring-membranes" };
        format!("{{\"model\":\"{model}\",\"regions\":[{regions}],\"band_hz\":[{:e},{:e}],\"propagator\":\"original-exact-dg\"}}",
            band[0], band[1])
        }
    }
    let main = source(layout);
    if main == "null" && layout.chimneys().is_empty() {
        return ",\"solid_wall\":null,\"solid_wall_scope\":\"rigid\"".into();
    }
    let chimneys = if layout.chimneys().is_empty() { String::new() } else {
        let records = layout.chimneys().iter().map(|chimney| format!(
            "{{\"port\":{},\"source\":{}}}", chimney.port().0, source(chimney.layout())
        )).collect::<Vec<_>>().join(",");
        format!(",\"chimney_walls\":[{records}]")
    };
    format!(",\"solid_wall\":{main}{chimneys},\"solid_wall_scope\":\"original cylindrical main-bore and explicitly assigned finite-chimney sources; physical regions aligned to each duct's DG faces; reciprocal frozen-state radial inertia and material response in each original gas propagator; omitted chimneys rigid; no axial shell waves or exterior wall radiation\"")
}

fn flue_boundary_provenance(p: &FlueFilePerformance) -> String {
    let network = p.renderer().voice().network().expect("file admits graph");
    let records = p.boundary_records().iter().map(|record| {
        let aperture = network.boundary_apertures().expect("source-bound apertures")
            .aperture(record.node()).expect("original boundary address");
        let drive = network.boundary_drive(record.node()).expect("original prescribed drive");
        material_boundary_record(record, aperture, drive)
    }).collect::<Vec<_>>().join(",");
    moving_boundary_provenance(&records)
}

fn material_boundary_record(
    record: &fs_couple::bernoulli_aperture::material_boundary::MaterialBoundaryRecord,
    aperture: &fs_couple::bernoulli_aperture::dynamic::DynamicAperture,
    drive: &fs_couple::bernoulli_aperture::network::boundaries::BoundaryDrive,
) -> String {
        let plate = aperture.plate_reduction().expect("source-bound plate");
        let chart = plate.material_chart().expect("retained material chart");
        let material = &chart.regions()[0].input.material;
        let spec = aperture.spec();
        let loss = if chart.bending_viscosity().is_some() { "kelvin-voigt" }
            else if aperture.relaxation().is_some() { "maxwell" } else { "elastic" };
        let bending_viscosity = chart.bending_viscosity().map_or_else(|| "null".into(), |laws|
            match (laws[0].isotropic(), laws[0].omega_band_rad_s()) {
                (Some(law), (lo, hi)) => format!("{{\"viscosity_pa_s\":{:e},\"omega_band_rad_s\":[{lo:e},{hi:e}]}}",
                    law.viscosity_pa_s),
                (None, (lo, hi)) => format!("{{\"plane_stress_tensor\":true,\"omega_band_rad_s\":[{lo:e},{hi:e}]}}"),
            });
        let terms = aperture.relaxation().map_or_else(String::new, |memory|
            memory.spec().regions[0].material.terms.iter().map(|(modulus, time)|
                format!("{{\"modulus_pa\":{modulus:e},\"relaxation_time_s\":{time:e}}}"))
                .collect::<Vec<_>>().join(","));
        let slit = aperture.viscous_slit().map_or_else(|| "null".into(), |law|
            format!("{{\"length_m\":{:e},\"dynamic_viscosity_pa_s\":{:e},\"source\":\"local-chimney-gas\"}}",
                law.length_m(), law.dynamic_viscosity_pa_s()));
        format!("{{\"node\":{},\"model\":\"material-cantilever-slit-v1\",\"propagator\":\"exact-zoh\",\"pack_path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"resolved_state_hash\":\"{}\",\"force_entity\":{},\"force_parameter\":\"closing_force\",\"force_dimensions\":[1,1,-2,0,0,0],\"initial_closing_force_n\":{:e},\"exterior_pressure_pa\":{:e},\"mass_kg\":{:e},\"stiffness_n_m\":{:e},\"pressure_area_m2\":{:e},\"slit_width_m\":{:e},\"rest_gap_m\":{:e},\"local_density_kg_m3\":{:e},\"local_impedance_pa_s_m3\":{:e},\"material_loss\":\"{loss}\",\"damping_ratio\":{:e},\"bending_viscosity\":{bending_viscosity},\"maxwell_terms\":[{terms}],\"viscous_slit\":{slit}}}",
            record.node(), json_string(record.path()), record.pack().content_hash().to_hex(),
            material.card_identity().to_hex(), material.identity().to_hex(),
            json_string(&record.closing_force().entity.token()), -drive.opening_force_n,
            drive.exterior_pressure_pa, spec.mass_kg, spec.stiffness_n_m, plate.pressure_area_m2(),
            plate.width_m(), spec.aperture.rest_opening_m, spec.density_kg_m3, spec.impedance_pa_s_m3,
            spec.damping_ratio)
}

fn moving_boundary_provenance(records: &str) -> String {
    format!(",\"moving_boundaries\":[{records}],\"boundary_model_scope\":\"one material plate opening mode with contact, swept volume and retained chimney state; prescribed exterior pressure; not a circular pad or resolved key linkage\"")
}

fn modal_wetted_provenance(
    valve: &fs_couple::bernoulli_aperture::modal::ModalAperture,
) -> String {
    use fs_plate::shell::pressure::ShellFace;
    let surfaces = valve.material_plates().iter().enumerate().filter_map(|(blade, plate)| {
        let port = plate.wetted_surface()?;
        let surface = port.surface();
        let facets = surface.upstream_facets().iter().map(|facet| {
            let side = match facet.side {
                ShellFace::Positive => "positive",
                ShellFace::Negative => "negative",
            };
            format!("{{\"triangle\":{},\"upstream_face\":\"{side}\"}}", facet.triangle)
        }).collect::<Vec<_>>().join(",");
        let areas = port.modal_area_m2().iter().map(|area| format!("{area:e}"))
            .collect::<Vec<_>>().join(",");
        Some(format!("{{\"blade\":{blade},\"specimen_identity\":\"{}\",\"source\":{},\"area_m2\":{:e},\"upstream_facets\":[{facets}],\"modal_pressure_areas_m2\":[{areas}]}}",
            surface.specimen_identity().to_hex(), json_string(surface.provenance()), port.area_m2()))
    }).collect::<Vec<_>>().join(",");
    format!(",\"wetted_surfaces\":[{surfaces}],\"wetted_surface_scope\":\"explicit original triangle faces share one signed projection for pressure traction and swept volume; omitted blades retain the complete pressure surface\"")
}

fn modal_boundary_provenance(
    p: &fs_couple::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance,
) -> String {
    use fs_couple::bernoulli_aperture::two_sided::ValveSide;
    let renderer = p.renderer();
    let records = p.boundary_records().iter().map(|(side, record)| {
        let aperture = renderer.source().system().boundary_apertures(*side)
            .expect("source-bound material keys").aperture(record.node())
            .expect("original boundary address");
        let drive = renderer.held_boundary_drives(*side).iter()
            .find(|drive| drive.node == record.node()).expect("complete explicit initial key drives");
        let name = match side { ValveSide::Upstream => "upstream", ValveSide::Downstream => "downstream" };
        let source = material_boundary_record(record, aperture, drive);
        let controls = renderer.pending_boundary_controls(*side).iter()
            .filter(|control| control.drive.node == record.node()).map(|control| format!(
                "{{\"mechanical_sample\":{},\"closing_force_n\":{:e},\"exterior_pressure_pa\":{:e}}}",
                control.sample, -control.drive.opening_force_n, control.drive.exterior_pressure_pa))
            .collect::<Vec<_>>().join(",");
        format!("{{\"side\":\"{name}\",\"source\":{source},\"controls\":[{controls}]}}")
    }).collect::<Vec<_>>().join(",");
    format!("{},\"boundary_control_clock_hz\":{}",
        moving_boundary_provenance(&records), p.info().mechanics_rate_hz)
}

struct Loaded {
    source: Box<dyn PressureRenderer>, rate:u32, samples:u64, full_scale_pa:f64,
    hash:ContentHash, fixture:&'static str, source_json:String,
    // A source that already owns multirate physical observation reports THAT
    // clock/filter history; the identity PCM wrapper must not call it 48k mechanics.
    observation_json:Option<String>,
}

pub(super) fn resolved_provenance(
    p: &fs_couple::nsf_dg::performance::file::ResolvedDgFilePerformance,
) -> String {
    let info = p.info();
    let renderer = p.renderer();
    let system = renderer.system();
    let law = system.zone().law();
    let limits = system.limits();
    let load = system.load();
    let duct = load.duct();
    let (point, mut observer) = match renderer.observation() {
        ApertureObservation::Inlet => ("mean-inlet",
            ",\"observation_scope\":\"internal interval-mean DG pressure; not an exterior microphone\"".into()),
        ApertureObservation::TubeBaffled(receiver) => ("baffled-outlet",
            super::duct_metadata::dg_baffled_projected_observation_provenance(load, receiver,
                renderer.baffled_receiver().expect("admitted resolved microphone"), info.sample_rate_hz)),
        _ => unreachable!("resolved input admits only its actual DG pressure observers"),
    };
    let gas = super::duct_metadata::gas_state(system.reference_gas());
    let gas_sources = super::thermochemical::provenance(p.gas_records());
    let solid_wall = dg_wall_provenance(p.wall_layout());
    let radiation = super::duct_metadata::dg_radiation_provenance(load);
    let propagation = super::duct_metadata::observer_gas_provenance(renderer.observer_gas());
    let sgs = info.sgs.map_or_else(|| ",\"sgs\":null".into(), |config| {
        let model = match config.closure.model {
            fs_lbm::sgs3::SgsModel::Smagorinsky { .. } => "smagorinsky",
            fs_lbm::sgs3::SgsModel::Wale { .. } => "wale",
            fs_lbm::sgs3::SgsModel::Vreman { .. } => "vreman",
        };
        format!(",\"sgs\":{{\"usage\":\"exploration\",\"model\":\"{model}\",\"coefficient\":{:e},\"filter_width_lu\":{:e},\"filter_width_m\":{:e},\"maximum_eddy_to_molecular\":{:e},\"turbulent_prandtl\":{:e},\"active\":{},\"scope\":\"original NSF deviatoric stress and conservative heat transport; open-face gradients use the original characteristic or reservoir state; no synthetic noise or turbulence qualification\"}}",
            config.closure.model.coefficient(), config.closure.filter_width_lu,
            config.closure.filter_width_lu * law.dx_m, config.closure.maximum_eddy_to_molecular,
            config.turbulent_prandtl, system.zone().sgs_config().is_some())
    });
    let aperture = [
        (system.coupled_face_is_high(), "coupled_patch",
            "selected aperture pressure, volume flow and acoustic work feed the bore; complete face mass and energy include escaping flow; no exterior radiation observer on the remainder"),
        (!system.coupled_face_is_high(), "source_patch",
            "selected aperture receives the declared source pressure or finite upstream duct; source work and upstream duct flow use only this aperture; complete face mass and energy include the independently driven remainder; no extra microphone source"),
    ].into_iter().filter_map(|(high, name, scope)| {
        let patch = system.zone().characteristic_patch_on(high)?;
        let (kind, value, units, temperature) = match system.patch_background_drive_on(high).expect("explicit patch background") {
            fs_lbm::nsf3::NsfBoundaryDrive::PressureOffset(value) => ("prescribed-pressure-offset", value, "Pa", String::new()),
            fs_lbm::nsf3::NsfBoundaryDrive::IncomingCharacteristic(value) => ("incoming-characteristic", value, "m/s", String::new()),
            fs_lbm::nsf3::NsfBoundaryDrive::PressureTemperature { pressure_offset_pa, inflow_temperature_k } =>
                ("prescribed-pressure-and-inflow-temperature", pressure_offset_pa, "Pa",
                    format!(",\"inflow_temperature_k\":{inflow_temperature_k:e},\"temperature_scope\":\"prescribed only at inward flow; outgoing entropy remains fluid-owned\"")),
        };
        let open = system.zone().characteristic_axis().expect("original characteristic axis");
        Some(format!(",\"{name}\":{{\"normal_axis\":{},\"high_face\":{},\"start\":[{},{}],\"end\":[{},{}],\"bounds\":\"half-open original cells on increasing tangential axes\",\"remainder_boundary\":{{\"kind\":\"{kind}\",\"value\":{value:e},\"units\":\"{units}\"{temperature}}},\"scope\":\"{scope}\"}}",
            open.axis, high, patch.start[0], patch.start[1], patch.end[0], patch.end[1]))
    }).collect::<String>();
    let reservoir = if system.zone().reservoir_aperture_count() == 0 { String::new() } else {
        let openings = (0..system.zone().reservoir_aperture_count()).map(|aperture| {
            let opening = system.zone().reservoir_aperture_at(aperture).expect("declared original aperture");
            let state = opening.state;
            let kind = if system.finite_reservoir_aperture() == Some(aperture) {
                "finite-thermochemical-chamber"
            } else { "prescribed-exterior-ghost-state" };
            format!("{{\"aperture\":{aperture},\"kind\":\"{kind}\",\"normal_axis\":{},\"high_face\":{},\"cell_count\":{},\"area_m2\":{:e},\"density_kg_m3\":{:e},\"velocity_m_s\":[{:e},{:e},{:e}],\"temperature_k\":{:e},\"mass_account\":{{\"field\":\"reservoir_exchanges[aperture].mass_out_kg\",\"unit\":\"kg\",\"positive\":\"out of resolved fluid\"}},\"energy_account\":{{\"field\":\"reservoir_exchanges[aperture].total_energy_out_j\",\"unit\":\"J\",\"positive\":\"out of resolved fluid\"}},\"scope\":\"original material-owned reservoir on independent transverse wall tiles; per-owner streamed mass and conservative FV energy; observation and acoustic reaction select an explicit owner\"}}",
                opening.axis, opening.high, opening.cells.len(), opening.cells.len() as f64 * law.dx_m.powi(2),
                state.density_kg_m3, state.velocity_m_s[0], state.velocity_m_s[1], state.velocity_m_s[2], state.temperature_k)
        }).collect::<Vec<_>>().join(",");
        format!(",\"reservoir_apertures\":[{openings}]")
    };
    let reservoir_observation = if renderer.reservoir_observation_count() == 0 {
        String::new()
    } else {
        // The bore helper describes one component of the complete microphone.
        observer = observer.replacen("\"observation_scope\":", "\"bore_observation_scope\":", 1);
        let records = (0..system.zone().reservoir_aperture_count()).filter_map(|aperture| {
            let input = renderer.reservoir_observation_at(aperture)?;
            let microphone = renderer.reservoir_receiver_at(aperture).expect("bound reservoir receiver");
            let material = super::duct_metadata::observer_gas_provenance(renderer.reservoir_observer_gas_at(aperture));
            let tiles = system.zone().last_reservoir_flows().iter()
                .filter(|flow| flow.aperture == aperture).map(|flow|
                format!("[{},{},{}]", flow.cell[0], flow.cell[1], flow.cell[2])).collect::<Vec<_>>().join(",");
            Some(format!("{{\"aperture\":{aperture},\"kind\":\"original-tile-continuum-flow-Rayleigh\",\"origin_m\":[{:e},{:e}],\"angle_rad\":{:e},\"grid_cells\":[{tiles}],\"tile_area_m2\":{:e},\"flow_account\":\"per-tile outward continuum volume flow; fluid/ghost face midpoint and original time trapezoid, distinct from streamed mass\",\"placement\":\"cyclic tangential axes, second reversed on low faces, rotated into common baffle\",\"kernel_support_observation_samples\":[{},{}],\"filter_delay_observation_samples\":{:e},\"observer_gas\":{material},\"radiation_feedback_added\":{}}}",
                input.origin_m[0], input.origin_m[1], input.angle_rad, law.dx_m.powi(2),
                microphone.delay_samples.0, microphone.delay_samples.1, microphone.filter_delay_samples(),
                system.reservoir_radiation_at(aperture).is_some()))

        }).collect::<Vec<_>>();
        let surfaces = if records.len() == 1 {
            format!(",\"reservoir_observation\":{}", records[0])
        } else {
            format!(",\"reservoir_observations\":[{}]", records.join(","))
        };
        format!(",\"observation_scope\":\"signed pressure sum of original DG outlets and individual resolved reservoir tiles in one explicit rigid baffle and material exterior; observation adds no reaction or flow/noise source; no free-field body or room scattering\"{surfaces}")
    };
    let mut reservoir_feedback = system.reservoir_radiation().map_or_else(String::new, |input| {
        let aperture = system.reservoir_radiation_aperture().expect("bound reaction aperture");
        format!(",\"reservoir_feedback\":{{\"aperture\":{aperture},\"kind\":\"same-interval-original-joint-DG-mouth-radiation\",\"origin_m\":[{:e},{:e}],\"mouth_quadrature\":[{},{}],\"ghost_pressure_bracket_pa\":[{:e},{:e}],\"flow_coordinates\":[\"into bore\",\"outward through resolved reservoir\"],\"scope\":\"independent original flows; full computed self and mutual exterior response; pressure and acoustic work use the original cached conservative FV predictors over the complete physical interval; endpoint tile traces remain observation-only; pressure and work checked before fluid, bore and radiation publish; uniform pressure projection on the original reservoir ghost boundary; no transparent-boundary or free-field qualification\"}}",
            input.origin_m[0], input.origin_m[1], input.cells[0], input.cells[1],
            input.pressure_bracket_pa[0], input.pressure_bracket_pa[1])
    });
    if let Some(solve) = system.radiation_solve() {
        let openings = system.reservoir_radiations().iter().enumerate()
            .map(|(index, (aperture, input))| format!(
                "{{\"aperture\":{aperture},\"joint_flow_coordinate\":{},\"origin_m\":[{:e},{:e}],\"mouth_quadrature\":[{},{}],\"ghost_pressure_bracket_pa\":[{:e},{:e}]}}",
                index + 1, input.origin_m[0], input.origin_m[1],
                input.cells[0], input.cells[1],
                input.pressure_bracket_pa[0], input.pressure_bracket_pa[1]))
            .collect::<Vec<_>>().join(",");
        let config = solve.newton;
        let globalization = match config.globalization {
            fs_solver::Globalization::LineSearch(policy) => format!(
                "{{\"kind\":\"line-search\",\"armijo\":{:e},\"contraction\":{:e},\"minimum_step\":{:e}}}",
                policy.armijo, policy.contraction, policy.minimum_step),
            fs_solver::Globalization::TrustRegion(policy) => format!(
                "{{\"kind\":\"trust-region\",\"initial_radius\":{:e},\"minimum_radius\":{:e},\"maximum_radius\":{:e},\"acceptance_ratio\":{:e},\"expansion_ratio\":{:e},\"shrink\":{:e},\"expansion\":{:e}}}",
                policy.initial_radius, policy.minimum_radius, policy.maximum_radius,
                policy.acceptance_ratio, policy.expansion_ratio, policy.shrink, policy.expansion),
        };
        reservoir_feedback = format!(
            ",\"reservoir_feedbacks\":[{openings}],\"reservoir_feedback_solve\":{{\"kind\":\"original-Newton-Krylov-shared-exterior\",\"maximum_iterations\":{},\"normalized_difference_step\":{:e},\"absolute_tolerance_pa\":{:e},\"relative_tolerance\":{:e},\"linear_restart\":{},\"maximum_linear_cycles\":{},\"forcing_minimum\":{:e},\"forcing_maximum\":{:e},\"forcing_gamma\":{:e},\"forcing_exponent\":{:e},\"globalization\":{globalization},\"flow_coordinates\":\"bore first, then outward through each listed original reservoir aperture\",\"scope\":\"one original joint self/mutual radiation owner; every resolved preview charged to the shared evaluation and cell-step bounds; out-of-bracket proposals refuse without clipping; each original pressure and work balance checked before fluid, both ducts, finite source and radiation publish; uniform pressure projection, not transparent or free-field qualification\"}}",
            solve.maximum_iterations, solve.difference_step,
            config.absolute_tolerance, config.relative_tolerance,
            config.linear_restart, config.max_linear_cycles,
            config.forcing_minimum, config.forcing_maximum,
            config.forcing_gamma, config.forcing_exponent);
    }
    let boundary = match renderer.drive_history().first().expect("finite native drive history") {
        fs_couple::nsf_dg::NsfDgDrive::Boundary(fs_lbm::nsf3::NsfBoundaryDrive::PressureTemperature { .. }) =>
            ",\"external_boundary\":{\"kind\":\"prescribed-pressure-and-inflow-temperature\",\"pressure_units\":\"Pa offset from source reference\",\"inflow_temperature_units\":\"K\",\"flow\":\"determined by the live outgoing characteristic at every original face reconstruction\",\"temperature\":\"held physical source temperature for inward flow; outgoing fluid entropy is retained\",\"material\":\"original source equation of state, caloric and transport laws; fixed composition\"}",
        fs_couple::nsf_dg::NsfDgDrive::Boundary(fs_lbm::nsf3::NsfBoundaryDrive::PressureOffset(_)) =>
            ",\"external_boundary\":{\"kind\":\"prescribed-pressure-offset\",\"units\":\"Pa\",\"reference\":\"original source-face material pressure\",\"flow\":\"determined by the live outgoing characteristic at every original face reconstruction\"}",
        fs_couple::nsf_dg::NsfDgDrive::Boundary(fs_lbm::nsf3::NsfBoundaryDrive::IncomingCharacteristic(_)) =>
            ",\"external_boundary\":{\"kind\":\"incoming-characteristic\",\"units\":\"m/s\",\"definition\":\"delta outward normal velocity minus delta pressure over reference material impedance\"}",
        fs_couple::nsf_dg::NsfDgDrive::UpstreamPressureRise(_) =>
            ",\"external_boundary\":{\"kind\":\"upstream-series-pressure-rise\",\"units\":\"Pa\",\"definition\":\"resolved face pressure minus finite upstream duct pressure\",\"flow\":\"both original acoustic loads and resolved faces solved on the same physical interval\",\"source_work\":\"minus pressure rise times outward upstream volume flow times physical interval\"}",
        fs_couple::nsf_dg::NsfDgDrive::FiniteReservoirActuation { .. } =>
            ",\"external_boundary\":{\"kind\":\"finite-material-reservoir-actuation\",\"volume_rate_units\":\"m3/s; positive expands chamber\",\"heat_power_units\":\"W; positive adds heat to chamber\",\"definition\":\"accepted finite inventory supplies aperture pressure and temperature; authored chamber-volume and heat rates change that same original named-mixture or thermochemical inventory\",\"partition\":\"explicit first order: original NSF flow exchange, then heat/volume solve over the same physical interval; refine interval for source feedback accuracy\"}",
        fs_couple::nsf_dg::NsfDgDrive::FiniteReservoir { .. } =>
            ",\"external_boundary\":{\"kind\":\"finite-thermochemical-reservoir\",\"pressure_units\":\"Pa absolute\",\"temperature_units\":\"K\",\"definition\":\"accepted finite inventory supplies held pressure and temperature to the original aperture; its streamed mass and complete FV energy update the same inventory atomically with fluid and ducts\",\"partition\":\"explicit first order over each original physical interval; refine interval for source feedback accuracy\"}",
        fs_couple::nsf_dg::NsfDgDrive::Reservoir { .. }
        | fs_couple::nsf_dg::NsfDgDrive::ReservoirAperture { .. } =>
            ",\"external_boundary\":{\"kind\":\"transverse-reservoir-pressure-temperature\",\"pressure_units\":\"Pa absolute\",\"temperature_units\":\"K\",\"definition\":\"original gas derives density and energy at held pressure and temperature; declared ghost velocity retained\",\"flow\":\"resolved fluid and both attached acoustic loads respond together with zero axial series effort; an unloaded other face receives its unperturbed incoming characteristic\"}",
    };
    let boundary = match renderer.drive_history().first().expect("finite native drive history") {
        fs_couple::nsf_dg::NsfDgDrive::ReservoirAperture { aperture, .. }
        | fs_couple::nsf_dg::NsfDgDrive::FiniteReservoir { aperture }
        | fs_couple::nsf_dg::NsfDgDrive::FiniteReservoirActuation { aperture, .. } =>
            format!("{boundary},\"driven_reservoir_aperture\":{aperture}"),
        fs_couple::nsf_dg::NsfDgDrive::Reservoir { .. } =>
            format!("{boundary},\"driven_reservoir_aperture\":0"),
        _ => boundary.to_string(),
    };
    let finite_reservoir_actuation = renderer.drive_history().iter().find_map(|drive| match drive {
        fs_couple::nsf_dg::NsfDgDrive::FiniteReservoirActuation {
            maximum_relative_volume_change, ..
        } => Some(*maximum_relative_volume_change),
        _ => None,
    });
    let finite_reservoir = system.finite_reservoir_inventory().map_or_else(String::new, |parcel| {
        let aperture = parcel.aperture;
        let limits = system.finite_reservoir_limits().expect("bound finite source limits");
        let gas = super::duct_metadata::gas_state(&parcel.gas);
        let source_kind = match parcel.source_kind {
            fs_couple::nsf_dg::NsfFiniteReservoirSource::Thermochemical => "thermochemical",
            fs_couple::nsf_dg::NsfFiniteReservoirSource::NamedMixture => "named-mixture",
        };
        let (scope, actuation) = if let Some(limit) = finite_reservoir_actuation {
            ("original named-mixture or thermochemical law; well-mixed chamber thermalizes incoming kinetic energy; original substep exchange variations bound hidden withdrawal and energy excursions; prescribed heat and chamber-volume rates update the same retained inventory; no stored source momentum, anatomical player model, or mean transport into acoustic ducts",
                format!(",\"actuation\":{{\"maximum_relative_volume_change\":{limit:e},\"volume_rate_units\":\"m3/s; positive expands chamber\",\"heat_power_units\":\"W; positive adds heat\",\"increments\":\"physical rates times original NSF/DG interval; independent of output sample clock\",\"mechanical_work\":\"midpoint pressure times volume change; positive work by gas on expansion; separate from complete FV energy already exchanged at aperture\",\"partition\":\"first-order original flow exchange then heat/volume solve; all owners publish atomically\",\"scope\":\"prescribed well-mixed chamber actuation; fixed composition; original NASA9 heat-volume owner or analytic conservative midpoint with original named constant-Cv law\"}}"))
        } else {
            ("same original frozen named-mixture or thermochemical law; rigid well-mixed chamber thermalizes incoming kinetic energy; original substep exchange variations bound hidden withdrawal and energy excursions; complete FV energy includes work already; no stored source momentum, mechanical breathing effort, pump, or mean transport into acoustic ducts", String::new())
        };
        format!(",\"finite_reservoir\":{{\"aperture\":{aperture},\"source_kind\":\"{source_kind}\",\"volume_m3\":{:e},\"mass_kg\":{:e},\"gas\":{gas},\"sensible_energy_j\":{:e},\"energy_reference_j_kg\":{:e},\"maximum_relative_mass_exchange\":{:e},\"maximum_relative_energy_exchange\":{:e},\"mass_absolute_tolerance_kg\":{:e},\"energy_absolute_tolerance_j\":{:e},\"relative_tolerance\":{:e},\"scope\":\"{scope}\"{actuation}}}",
            parcel.volume_m3, parcel.mass_kg, parcel.sensible_energy_j,
            parcel.energy_reference_j_kg, limits.maximum_relative_mass_exchange,
            limits.maximum_relative_energy_exchange, limits.mass_absolute_tolerance_kg,
            limits.energy_absolute_tolerance_j, limits.relative_tolerance)
    });
    let upstream = system.upstream_load().map_or_else(String::new, |upstream| {
        let duct = upstream.duct();
        let gas = super::duct_metadata::gas_state(duct.inlet_gas());
        let wall = dg_wall_provenance(p.upstream_wall_layout().expect("native upstream source retains its wall layout"));
        let radiation = super::duct_metadata::dg_radiation_provenance(upstream);
        format!(",\"upstream_acoustic_load\":{{\"operator_identity\":\"{}\",\"total_dofs\":{},\"length_m\":{:e},\"inlet_area_m2\":{:e},\"inlet_gas\":{gas},\"scope\":\"original finite material DG duct and terminal histories; acoustic coupling without mean mass or enthalpy advection; pressure observer remains on the downstream duct\"{wall}{radiation}}}",
            duct.operator_identity().to_hex(), duct.total_dofs(), duct.total_length(), duct.inlet_area_m2())
    });
    let timing = format!(",\"time_projection\":{{\"physical_interval_s\":{:e},\"declared_physical_intervals\":{},\"declared_physical_duration_s\":{:e},\"observation_rate_hz\":{},\"declared_observation_samples\":{},\"declared_observation_duration_s\":{:e},\"method\":\"conservative overlap of original interval means; direct same-clock path\",\"boundary_control\":\"declared hold or linear ramp sampled at each physical interval start; event error at most one physical interval\",\"bandlimited_antialias_filter\":false,\"scope\":\"box-average temporal projection requires clock refinement; physical prefix may lead observation by less than one physical interval; shared PCM decimation may truncate an incomplete output interval\"}}",
        load.propagator().dt(), info.physical_intervals,
        load.propagator().dt() * info.physical_intervals as f64,
        info.sample_rate_hz, info.samples, info.samples as f64 / f64::from(info.sample_rate_hz));
    format!("\"resolved_input\":{{\"schema\":\"frankensim-resolved-dg-performance-v1\",\"blake3\":\"{}\",\"source_owner\":\"original-NsfGrid3-D3Q27\",\"dimensions\":[{},{},{}],\"solid_cells\":{},\"initial_state_cells\":{},\"boundary_segments\":{},\"boundary_ramps\":{},\"dx_m\":{:e},\"nsf_dt_s\":{:e},\"substeps\":{},\"source_gas\":{gas},\"bulk_to_shear_ratio\":{:e},\"bulk_ratio_source\":{},\"coupled_area_m2\":{:e},\"maximum_root_evaluations\":{},\"maximum_cell_steps\":{},\"pressure_absolute_tolerance_pa\":{:e},\"pressure_relative_tolerance\":{:e},\"projection_absolute_tolerance_j\":{:e},\"projection_relative_tolerance\":{:e},\"duct_operator_identity\":\"{}\",\"duct_total_dofs\":{},\"pressure_point\":\"{point}\",\"observer_gas\":{propagation},\"model_scope\":\"original resolved fixed-composition fluid and live material DG bore; uniform acoustic interface projection with separate mass/energy flux; stationary voxel walls and two-cell fluid collars; source, geometry, boundary and coupling convergence remain required; no measured flute/reed qualification\"{sgs}{aperture}{reservoir}{boundary}{upstream}{timing}{observer}{radiation}{solid_wall}{gas_sources}{reservoir_observation}{reservoir_feedback}{finite_reservoir}}}",
        info.input_hash.to_hex(), info.dimensions[0], info.dimensions[1], info.dimensions[2],
        info.solid_cells, info.initial_state_cells, info.boundary_segments, info.boundary_ramps, law.dx_m, law.dt_s,
        limits.substeps, law.bulk_to_shear_ratio, json_string(system.zone().bulk_ratio_source()),
        system.area_m2(), limits.maximum_evaluations, limits.maximum_cell_steps,
        limits.pressure_absolute_tolerance_pa, limits.pressure_relative_tolerance,
        limits.projection_absolute_tolerance_j, limits.projection_relative_tolerance,
        duct.operator_identity().to_hex(), duct.total_dofs())
}

fn flue_duct_provenance(p: &fs_couple::flue_loop::network::performance::FlueOutletPerformance) -> String {
    use fs_couple::bernoulli_aperture::network::NetworkNode;
    let source = p.voice().network().expect("physical flue network").source();
    let inlet = source.graph.inlet_section_index().expect("unique inlet");
    let nodes = source.graph.nodes.iter().map(|node| match node {
        NetworkNode::Inlet => "{\"kind\":\"inlet\"}".to_string(),
        NetworkNode::Driven => "{\"kind\":\"driven\"}".to_string(),
        NetworkNode::Junction => "{\"kind\":\"junction\"}".to_string(),
        NetworkNode::Termination { reflection } => format!("{{\"kind\":\"reflection\",\"reflection\":{reflection:e}}}"),
        NetworkNode::Relaxation { load } => {
            let base = load.base();
            let compliance = base.compliance_m3_pa.map_or_else(|| "null".to_string(), |c| format!("{c:e}"));
            let terms = load.terms().iter().map(|term| format!(
                "{{\"resistance_pa_s_m3\":{:e},\"rate_per_s\":{:e}}}",
                term.resistance_pa_s_m3,term.rate_per_s)).collect::<Vec<_>>().join(",");
            format!("{{\"kind\":\"relaxation\",\"resistance_pa_s_m3\":{:e},\"inertance_pa_s2_m3\":{:e},\"compliance_m3_pa\":{compliance},\"terms\":[{terms}]}}",
                base.resistance_pa_s_m3,base.inertance_pa_s2_m3)
        }
        node => {
            let (kind, load) = match node {
                NetworkNode::Impedance { load } => ("impedance", *load),
                NetworkNode::Series { load } => ("series", load.base()),
                NetworkNode::Shunt { load } => ("shunt", load.base()),
                _ => unreachable!("file declares only explicit R-L-C loads"),
            };
            let compliance = load.compliance_m3_pa.map_or_else(|| "null".to_string(), |c| format!("{c:e}"));
            format!("{{\"kind\":\"{kind}\",\"resistance_pa_s_m3\":{:e},\"inertance_pa_s2_m3\":{:e},\"compliance_m3_pa\":{compliance}}}",
                load.resistance_pa_s_m3, load.inertance_pa_s2_m3)
        }
    }).collect::<Vec<_>>().join(",");
    let sections = source.graph.sections.iter().zip(&source.gases).map(|(s,g)| format!(
        "{{\"nodes\":[{},{}],\"length_m\":{:e},\"radius_m\":{:e},\"temperature_k\":{:e},\"pressure_pa\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"dynamic_viscosity_pa_s\":{:e}}}",
        s.nodes[0],s.nodes[1],s.length_m,s.radius_m,g.temperature,g.pressure,g.density,g.sound_speed,g.dynamic_viscosity
    )).collect::<Vec<_>>().join(",");
    let radiation = p.radiation_input().expect("explicit flue exterior");
    let mut outlets = format!("{{\"node\":{},\"center_m\":[0,0]}}",p.config().outlet_node);
    for outlet in &radiation.additional_outlets {
        outlets.push_str(&format!(",{{\"node\":{},\"center_m\":[{:e},{:e}]}}",
            outlet.node,outlet.center_m[0],outlet.center_m[1]));
    }
    let mut text = format!(",\"duct_network\":{{\"inlet_section\":{inlet},\"nodes\":[{nodes}],\"sections\":[{sections}],\"outlets\":[{outlets}]}}");
    if source.graph.sections.len() == 1 {
        let section = source.graph.sections[0];
        text.push_str(&format!(",\"tube_length_m\":{:e},\"tube_radius_m\":{:e}",section.length_m,section.radius_m));
        if let NetworkNode::Termination { reflection } = source.graph.nodes[p.config().outlet_node] {
            text.push_str(&format!(",\"terminal_reflection\":{reflection:e}"));
        }
    }
    text
}

/// Original material/DG source metadata as one named JSON member (no leading comma).
pub(super) fn material_modal_dg_provenance(
    p: &fs_couple::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance,
) -> String {
    let i = p.info();
    let source = p.renderer();
    let load = source.system().load();
    let (pressure_point, observation) = match source.observation() {
        ApertureObservation::Inlet => ("mean-inlet", String::new()),
        ApertureObservation::TubeBaffled(receiver) => (
            "baffled-outlet",
            super::duct_metadata::dg_baffled_observation_provenance(
                load, receiver, source.baffled_receiver().expect("admitted DG receiver"),
            ),
        ),
        _ => unreachable!("authored material DG observes its original inlet or outlet"),
    };
    let materials = if let Some(solid) = p.solid_record() {
        solid.materials().iter().map(|record| format!(
            "{{\"section\":{},\"path\":{},\"pack_hash\":\"{}\"}}",
            record.section(), json_string(record.path()), record.pack().content_hash().to_hex(),
        )).collect::<Vec<_>>().join(",")
    } else {
        p.material_records().iter().map(|record| format!(
            "{{\"blade\":{},\"section\":{},\"region\":{},\"path\":{},\"pack_hash\":\"{}\"}}",
            record.blade(), record.section(), record.region(), json_string(record.path()),
            record.pack().content_hash().to_hex(),
        )).collect::<Vec<_>>().join(",")
    };
    let inlet = super::duct_metadata::gas_state(load.duct().inlet_gas());
    let outlet = super::duct_metadata::gas_state(load.duct().outlet_gas());
    let gas_sources = super::thermochemical::provenance(p.gas_records());
    let radiation = super::duct_metadata::dg_radiation_provenance(load);
    let solid_wall = dg_wall_provenance(p.wall_layout());
    let observation = format!("{observation},\"observer_gas\":{}",
        super::duct_metadata::observer_gas_provenance(source.observer_gas()));
    let upstream = source.system().upstream_load().map_or_else(String::new, |tract| {
        let duct = tract.duct();
        let mesh = duct.mesh();
        let field = duct.gas_field().map_or_else(String::new, super::duct_metadata::dg_field_provenance);
        let radiation = super::duct_metadata::dg_radiation_provenance(tract);
        let wall = dg_wall_provenance(p.upstream_wall_layout()
            .expect("source-bound upstream DG wall layout"));
        format!(",\"upstream_tract\":{{\"dg_operator_identity\":\"{}\",\"length_m\":{:e},\"inlet_gas\":{},\"outlet_gas\":{},\"elements\":{},\"order\":{},\"total_dofs\":{},\"scope\":\"same-step original upstream and downstream mean-port laws in the sole material-valve solve; blowing pressure is a series source; matching inlet gases; original upstream terminal exchange stays separate from the downstream microphone\"{field}{radiation}{wall}}}",
            duct.operator_identity().to_hex(), duct.total_length(),
            super::duct_metadata::gas_state(duct.inlet_gas()),
            super::duct_metadata::gas_state(duct.outlet_gas()),
            mesh.elements, mesh.order, duct.total_dofs())
    });
    let relaxation_json = match source.system().aperture() {
        fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(valve) =>
            modal_relaxation_provenance(valve),
        _ => String::new(),
    };
    let (structure, model_scope) = if let Some(record) = p.solid_record() {
        let fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(valve) =
            source.system().aperture() else { unreachable!("solid source owns a modal aperture") };
        if let Some(solid) = valve.hereditary_solid() {
            let image = solid.source();
            let band = image.frequency_band_rad_s();
            (format!(",\"solid\":{{\"element_order\":{},\"tetrahedra\":{},\"source_identity\":\"{}\",\"coordinate_mass_kg\":{:e},\"material_response\":\"hereditary\",\"material_quadrature_points\":{},\"memory_channels\":{},\"basis_error\":{:e},\"frequency_band_rad_s\":[{:e},{:e}],\"instantaneous_frequency_bound_rad_s\":{:e},\"selected_zero_opening_modes_retained\":true}}",
                json_string(record.element_order()), record.tetrahedra(), image.identity().to_hex(),
                solid.coordinate_mass_kg(), image.material_points().len(), image.memory_channel_count(),
                image.basis_error(), band.0, band.1, image.instantaneous_frequency_bound_rad_s()),
                "original three-dimensional regional tensor material, complete selected equilibrium basis and causal Maxwell history in the shared exact pressure/contact/DG owner; explicit source clock, strain, basis and capacity gates; mesh, mode, time and DG refinement remain separate; no measured-instrument qualification")
        } else {
        let solid = valve.material_solid().expect("native solid retains its original material image");
        let image = solid.source();
        let approximation = image.approximation();
        (format!(",\"solid\":{{\"element_order\":{},\"tetrahedra\":{},\"source_identity\":\"{}\",\"coordinate_mass_kg\":{:e},\"reference_omega_rad_s\":{:e},\"maximum_storage_sample_error\":{:e},\"maximum_loss_sample_error\":{:e},\"full_modal_viscosity\":true,\"selected_zero_opening_modes_retained\":true}}",
            json_string(record.element_order()), record.tetrahedra(), image.identity().to_hex(),
            solid.coordinate_mass_kg(), approximation.reference_omega_rad_s,
            approximation.maximum_storage_sample_error, approximation.maximum_loss_sample_error),
            "original three-dimensional regional material solid with all selected modes, signed pressure/force maps and full coupled viscosity in the shared DG aperture solve; frozen-frequency material image with sampled dispersion checks, not hereditary solid dynamics; mesh, mode, time and DG refinement remain separate; no measured-instrument qualification")
        }
    } else {
        (String::new(), "original regional material plates, signed modal pressure/force projections and retained material memory coupled to the shared exact DG bore; flat DKT and finite-mode interface; nonlinear time and DG spatial convergence require separate refinement; no measured-instrument qualification")
    };
    let surface_closure = match source.system().aperture() {
        fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(valve) =>
            valve.solid_surface_closure().map_or_else(String::new, |spec| {
                let contact = &spec.contact;
                let faces = contact.lay_triangles.iter().map(usize::to_string)
                    .collect::<Vec<_>>().join(",");
                let edges = spec.slit_edges.iter().map(|edge| format!("[{},{}]", edge[0], edge[1]))
                    .collect::<Vec<_>>().join(",");
                let gaps = contact.nodal_rest_gap_m.iter().map(|gap| format!("{gap:e}"))
                    .collect::<Vec<_>>().join(",");
                format!(",\"solid_surface_closure\":{{\"original_exterior_face_indices\":[{faces}],\"slit_edge_original_nodes\":[{edges}],\"opening_sign\":{:e},\"normal_convention\":\"gap = nodal_rest_gap + opening_sign * original_outward_normal dot original_nodal_displacement\",\"nodal_rest_gap_m\":[{gaps}],\"contact_law\":{{\"stiffness_pa_per_m_alpha\":{:e},\"alpha\":{:e},\"internal_loss_s_per_m\":{:e},\"max_penetration_m\":{:e},\"provenance\":{}}},\"scope\":\"fixed-reference Tet4 parallel-normal patches; original area-lumped contact and positive P1 slit integral on all retained modes; original pressure and swept-volume projections remain the shared conjugate work port; no resolved moving-channel flow\"}}",
                    spec.opening_sign, contact.stiffness_pa_per_m_alpha, contact.alpha,
                    contact.internal_loss_s_per_m, contact.max_penetration_m,
                    json_string(&contact.provenance))
            }),
        _ => String::new(),
    };
    let source_json = format!("\"material_modal_dg_valve_input\":{{\"schema\":\"frankensim-material-modal-dg-valve-v1\",\"blake3\":\"{}\",\"blades\":{},\"retained_modes\":{},\"nodes\":{},\"triangles\":{},\"mechanics_rate_hz\":{},\"pressure_point\":\"{pressure_point}\",\"material_sources\":[{materials}],\"structural_propagator\":\"exact-zoh\",\"dg_operator_identity\":\"{}\",\"length_m\":{:e},\"inlet_gas\":{inlet},\"outlet_gas\":{outlet},\"model_scope\":\"{model_scope}\"{structure}{surface_closure}{observation}{radiation}{gas_sources}{solid_wall}{upstream}{relaxation_json}}}",
        i.input_hash.to_hex(), i.blades, i.modes, i.nodes, i.triangles,
        i.mechanics_rate_hz, load.duct().operator_identity().to_hex(),
        load.duct().total_length());
    source_json
}

fn modal_relaxation_provenance(
    valve: &fs_couple::bernoulli_aperture::modal::ModalAperture,
) -> String {
    valve.material_relaxation().map_or_else(String::new, |memory| {
        let initial = match memory.initial_memory() {
            fs_couple::bernoulli_aperture::modal::relaxation::InitialModalMemory::Relaxed => "relaxed",
            fs_couple::bernoulli_aperture::modal::relaxation::InitialModalMemory::Unrelaxed => "unrelaxed",
        };
        let arms: usize = memory.memory_sqrt_j().iter().map(Vec::len).sum();
        if let Some(solid) = valve.hereditary_solid() {
            let branches = solid.source().branches().iter().map(|branch| format!(
                "{{\"relaxation_time_s\":{:e},\"channels\":{}}}",
                branch.relaxation_time_s(), branch.rows().len(),
            )).collect::<Vec<_>>().join(",");
            return format!(",\"solid_material_relaxation\":{{\"initial_memory\":\"{initial}\",\"live_arms\":{arms},\"branches\":[{branches}],\"scope\":\"original volume material laws and strain maps; equal-time factors retain every loss direction; Kelvin and Maxwell losses share the original exact structural owner\"}}");
        }
        let points = memory.tensor_material_points();
        if !points.is_empty() {
            let mut regions = std::collections::BTreeSet::new();
            let matrix_json = |matrix: &[[f64; 6]; 6]| matrix.iter().map(|row|
                format!("[{}]", row.iter().map(|v| format!("{v:e}"))
                    .collect::<Vec<_>>().join(","))).collect::<Vec<_>>().join(",");
            let laws = points.iter().filter_map(|point| {
                if !regions.insert((point.blade(), point.region())) { return None; }
                let law = point.image().point().material();
                let equilibrium = matrix_json(law.equilibrium_mandel());
                let branches = law.branches().iter().map(|branch| format!(
                    "{{\"relaxation_time_s\":{:e},\"stiffness_mandel_pa\":[{}]}}",
                    branch.relaxation_time_s(), matrix_json(branch.stiffness_mandel())))
                    .collect::<Vec<_>>().join(",");
                Some(format!("{{\"blade\":{},\"region\":{},\"band_rad_s\":[{:e},{:e}],\"equilibrium_mandel_pa\":[{equilibrium}],\"branches\":[{branches}]}}",
                    point.blade(), point.region(), law.band_rad_s().0, law.band_rad_s().1))
            }).collect::<Vec<_>>().join(",");
            return format!(",\"tensor_material_relaxation\":{{\"initial_memory\":\"{initial}\",\"live_arms\":{arms},\"material_quadrature_points\":{},\"basis\":\"spatial Mandel xx yy zz sqrt2_xy sqrt2_yz sqrt2_zx\",\"regional_laws\":[{laws}],\"scope\":\"complete 3D plane-stress history constraint; all source coordinates retained by an invertible passive transform; original DKT strains and outer-face material limits; fixed-state homogeneous linear bending\"}}", points.len());
        }
        let spectra = memory.material_spectra().iter().enumerate().flat_map(|(blade, regions)| {
            regions.iter().enumerate().map(move |(section, law)| {
                let terms = law.terms.iter().map(|(modulus, time)| format!(
                    "{{\"modulus_pa\":{modulus:e},\"relaxation_time_s\":{time:e}}}"))
                    .collect::<Vec<_>>().join(",");
                format!("{{\"blade\":{blade},\"section\":{section},\"equilibrium_young_modulus_pa\":{:e},\"terms\":[{terms}]}}", law.e_inf)
            })
        }).collect::<Vec<_>>().join(",");
        let kind = if memory.is_regional() { "regional_material_relaxation" }
            else { "proportional_material_relaxation" };
        format!(",\"{kind}\":{{\"initial_memory\":\"{initial}\",\"live_arms\":{arms},\"regional_spectra\":[{spectra}]}}")
    })
}

/// Original two-sided material source metadata shared by wind and ensemble.
pub(super) fn material_modal_provenance(
    p: &fs_couple::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance,
) -> String {
    let i = p.info();
    let valve = match p.renderer().source().system().valve() {
        fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(a) => a,
        _ => unreachable!("material-modal file owns the requested modal valve"),
    };
    let materials = p.material_records().iter().map(|record| {
        let source = valve.material_plates()[record.blade()].source();
        let region = &source.regions()[record.region()];
        let material = &region.input.material;
        format!("{{\"blade\":{},\"section\":{},\"region\":{},\"path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"state_hash\":\"{}\",\"specimen_hash\":\"{}\",\"mass_kg\":{:e},\"thickness_m\":{:e},\"material_angle_rad\":{:e}}}",
            record.blade(), record.section(), record.region(), json_string(record.path()),
            record.pack().content_hash().to_hex(), material.card_identity().to_hex(),
            material.identity().to_hex(), source.specimen_identity().to_hex(),
            region.mass_kg, region.section.thickness, region.input.material_angle_rad)
    }).collect::<Vec<_>>().join(",");
    let gas_sources = super::thermochemical::provenance(p.gas_records());
    let system = p.renderer().source().system();
    let walls = p.wall_records().iter().map(|record| {
        let side = match record.side() {
            fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream => "upstream",
            fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream => "downstream",
        };
        let wall = system.material_walls(record.side()).iter()
            .find(|wall| wall.source().node == record.node()).expect("source-bound duct wall");
        let material = &wall.source().material;
        let (model, directions) = graph_wall_constitutive_metadata(wall.source());
        let (loss, relaxation) = graph_wall_loss_metadata(wall.source());
        let patch = wall.patch();
        format!("{{\"side\":\"{side}\",\"node\":{},\"model\":\"{model}\",\"path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"state_hash\":\"{}\",\"thickness_m\":{:e},\"radius_m\":{:e},\"length_m\":{:e},\"band_hz\":[{:e},{:e}],\"loss\":\"{loss}\",\"surface_density_kg_m2\":{:e},\"stiffness_pa_m\":{:e},\"resistance_pa_s_m\":{:e}{relaxation}{directions}}}",
            record.node(), json_string(record.path()), record.pack().content_hash().to_hex(),
            material.card_identity().to_hex(), material.identity().to_hex(),
            wall.source().thickness_m, wall.radius_m(), wall.length_m(), wall.band_hz()[0], wall.band_hz()[1],
            patch.wall.surface_density, patch.wall.stiffness_per_area, patch.wall.resistance)
    }).collect::<Vec<_>>().join(",");
    let directional_walls = [fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream,
        fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream].into_iter().any(|side|
            system.material_walls(side).iter().any(|wall|
                !matches!(&wall.source().model, fs_couple::bernoulli_aperture::wall::material::WallMaterialModel::Isotropic)));
    let wall_model = if directional_walls { "isotropic or principal-aligned orthotropic" } else { "isotropic" };
    let wall_json = format!(",\"material_duct_walls\":[{walls}],\"duct_wall_scope\":\"fixed-state local {wall_model} membrane; original regional duct gas; no full shell or exterior wall radiation\"");
    let supply = system.supply_channel().map_or_else(|| "null".into(), |channel| {
        match channel {
            fs_couple::bernoulli_aperture::slit_channel::FixedChannel::Slit(channel) => {
                let spec = channel.spec();
                format!("{{\"length_m\":{:e},\"width_m\":{:e},\"height_m\":{:e},\"transverse_cells\":{},\"density_kg_m3\":{:e},\"viscosity_pa_s\":{:e}}}",
                    spec.length_m, spec.width_m, spec.height_m, spec.transverse_cells,
                    channel.gas().density, channel.gas().dynamic_viscosity)
            }
            fs_couple::bernoulli_aperture::slit_channel::FixedChannel::Circular(channel) => {
                let spec = channel.spec();
                format!("{{\"kind\":\"circular\",\"length_m\":{:e},\"radius_m\":{:e},\"radial_cells\":{},\"density_kg_m3\":{:e},\"viscosity_pa_s\":{:e}}}",
                    spec.length_m, spec.radius_m, spec.radial_cells,
                    channel.gas().density, channel.gas().dynamic_viscosity)
            }
        }
    });
    let (propagator, linear_scope) = match i.propagator {
        fs_couple::bernoulli_aperture::paired::BladePropagator::Midpoint =>
            ("midpoint", "original material structure with loaded one-cent midpoint admission"),
        fs_couple::bernoulli_aperture::paired::BladePropagator::ExactZoh =>
            ("exact-zoh", if valve.material_relaxation().is_some_and(|memory| !memory.tensor_material_points().is_empty()) {
                "original M/K and complete 3D plane-stress tensor histories on retained DKT strains; pressure and contact held per mechanical sample; integrated material dissipation; nonlinear time error requires refinement"
            } else if valve.material_relaxation().is_some_and(|memory| memory.is_regional()) {
                "original M/K and full cross-mode regional Maxwell internal strains; pressure and contact held per mechanical sample; integrated material dissipation; nonlinear time error requires refinement"
            } else {
                "original M/K and selected full C or proportional Maxwell internal strains; pressure and contact held per mechanical sample; integrated material dissipation; nonlinear time error requires refinement"
            }),
    };
    let relaxation_json = modal_relaxation_provenance(valve);
    let (radiation_scope, mimo_json) = if p.mimo_radiation_records().is_empty()
        && p.computed_radiation_records().is_empty() {
        ("compact independent baffled reaction and shared exterior pressure", String::new())
    } else {
        let groups = p.mimo_radiation_records().iter().map(|(side, record)| {
            let side = match side {
                fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream => "upstream",
                fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream => "downstream",
            };
            format!("{{\"side\":\"{side}\",\"source\":{}}}", super::duct_metadata::mimo_radiation_provenance(record))
        }).chain(p.computed_radiation_records().iter().map(|(side, record)| {
            let side = match side {
                fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream => "upstream",
                fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream => "downstream",
            };
            format!("{{\"side\":\"{side}\",\"source\":{}}}", super::duct_metadata::computed_baffled_provenance(record))
        })).collect::<Vec<_>>().join(",");
        ("source-bound mutual baffled reaction within each selected side and matched exterior pressure; unselected outlets retain authored loads; no cross-side exterior coupling",
            format!(",\"mimo_radiation_groups\":[{groups}]"))
    };
    let terminals = p.radiation_bindings().len()
        + p.mimo_radiation_records().iter().map(|(_, r)| r.apertures().len()).sum::<usize>()
        + p.computed_radiation_records().iter().map(|(_, r)| r.apertures().len()).sum::<usize>();
    let boundary_json = modal_boundary_provenance(p);
    let wetted_json = modal_wetted_provenance(valve);
    let opening_loss_json = opening_loss_provenance(
        [fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream,
         fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream]
            .into_iter().flat_map(|side| {
                let label = match side {
                    fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream => "upstream",
                    fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream => "downstream",
                };
                system.opening_losses(side).iter().map(move |binding| (label, binding))
            }));
    let force_programs = valve.material_plates().iter().enumerate().filter_map(|(blade, plate)| {
        let events = p.renderer().pending_plate_force_controls(blade);
        if events.is_empty() { return None; }
        let weights = plate.options().reduction.unit_force_weights.iter()
            .map(|weight| format!("{weight:e}")).collect::<Vec<_>>().join(",");
        let projection = valve.spec().modes[plate.mode_range()].iter()
            .map(|mode| format!("{:e}", mode.force_projection[blade]))
            .collect::<Vec<_>>().join(",");
        let events = events.iter().map(|event| format!(
            "{{\"mechanical_sample\":{},\"force_n\":{:e}}}", event.sample, event.force_n))
            .collect::<Vec<_>>().join(",");
        Some(format!("{{\"blade\":{blade},\"node_force_weights\":[{weights}],\"modal_force_projection\":[{projection}],\"controls\":[{events}]}}"))
    }).collect::<Vec<_>>().join(",");
    let force_json = format!(",\"source_force_programs\":[{force_programs}],\"source_force_clock_hz\":{},\"source_force_scope\":\"prescribed physical newtons through original plate footprints and all retained modal coordinates; no tongue or lip contact law\"", i.mechanics_rate_hz);
    let compiled_controls = p.renderer().pending_controls().len()
        + [fs_couple::bernoulli_aperture::two_sided::ValveSide::Upstream,
           fs_couple::bernoulli_aperture::two_sided::ValveSide::Downstream].into_iter()
            .map(|side| p.renderer().pending_boundary_controls(side).len()).sum::<usize>()
        + (0..i.blades).map(|blade| p.renderer().pending_plate_force_controls(blade).len()).sum::<usize>();
    let closure_scope = if valve.plate_closure().is_some() {
        "source-mesh distributed contact and positive-P1 slit flow on every retained mode"
    } else { "effective uniform slit/contact" };
    let source_json = format!("\"material_modal_valve_input\":{{\"schema\":\"frankensim-material-modal-valve-v1\",\"blake3\":\"{}\",\"blades\":{},\"retained_modes\":{},\"nodes\":{},\"triangles\":{},\"compiled_controls\":{},\"full_tensor_viscosity\":{},\"radiation_terminals\":{},\"structural_propagator\":\"{propagator}\",\"linear_scope\":\"{linear_scope}\",\"material_sources\":[{materials}],\"supply_channel\":{supply},\"model_scope\":\"regional material flat DKT walls; original signed pressure and force projections; {closure_scope}; both gas-bound ducts; {radiation_scope}; no measured reed, unbaffled or resolved moving-channel qualification\"{gas_sources}{mimo_json}{boundary_json}{wall_json}{force_json}{relaxation_json}{opening_loss_json}{wetted_json}}}",
        i.input_hash.to_hex(), i.blades, i.modes, i.nodes, i.triangles,
        compiled_controls, valve.viscosity_matrix_n_s_m().is_some(),
        terminals);
    source_json
}

// Share the original source receipts with ensemble composition without rebuilding
// any jet, material, gas, control or microphone state.
pub(super) fn flue_dg_provenance(
    p: &fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance,
) -> String {
    let i = p.info();
    let voice = p.renderer().voice();
    let load = voice.dg().expect("file retains original DG source").load();
    let duct = load.duct();
    let geometry = p.geometry();
    let band = duct.viscothermal_band_rad_s().map_or_else(|| "null".into(), |band|
        format!("[{:e},{:e}]", band[0], band[1]));
    let field = duct.gas_field().map_or_else(String::new,
        super::duct_metadata::dg_field_provenance);
    let inlet_gas = super::duct_metadata::gas_state(duct.inlet_gas());
    let outlet_gas = super::duct_metadata::gas_state(duct.outlet_gas());
    let gas_sources = super::thermochemical::provenance(p.gas_records());
    let radiation = super::duct_metadata::dg_radiation_provenance(load);
    let (pressure_point, pressure_observation) =
        if let Some(receiver) = p.renderer().baffled_receiver_config() {
            let terminal = p.renderer().baffled_receiver().expect("admitted DG microphone");
            if let Some(mouth) = p.renderer().mouth_radiation() {
                let source = load.radiation_record().expect("bound DG exterior");
                let exterior = super::duct_metadata::gas_state(source.gas_record().state());
                let microphone = p.renderer().mouth_receiver().expect("admitted DG mouth");
                let support = p.renderer().propagation_delay_samples().expect("exterior kernels");
                let mouth_reaction = voice.mouth_radiation_load();
                let mutual_mouth = load.has_driven_mouth();
                let terminal_json = if source.apertures().len() == 1 {
                    let aperture = &source.apertures()[0];
                    format!(",\"outlet_receiver\":{{\"radius_m\":{:e},\"radial_rings\":{},\"angular_points\":{},\"kernel_support_mechanical_samples\":[{},{}],\"filter_delay_mechanical_samples\":{:e}}}",
                        aperture.radius_m, receiver.radial_rings, receiver.angular_points,
                        terminal.delay_samples.0, terminal.delay_samples.1, terminal.filter_delay_samples())
                } else {
                    let apertures = load.terminal_ports().iter().zip(source.apertures())
                        .map(|(port, aperture)| format!("{{\"port\":{},\"center_m\":[{:e},{:e}],\"radius_m\":{:e}}}",
                            port.0, aperture.center_m[0], aperture.center_m[1], aperture.radius_m))
                        .collect::<Vec<_>>().join(",");
                    format!(",\"terminal_receiver\":{{\"apertures\":[{apertures}],\"flow_order\":\"original DG terminal channels\",\"radial_rings\":{},\"angular_points\":{},\"kernel_support_mechanical_samples\":[{},{}],\"filter_delay_mechanical_samples\":{:e}}}",
                        receiver.radial_rings, receiver.angular_points, terminal.delay_samples.0,
                        terminal.delay_samples.1, terminal.filter_delay_samples())
                };
                let mouth_json = format!(",\"mouth_receiver\":{{\"shape\":\"rectangle\",\"center_m\":[{:e},{:e}],\"size_m\":[{:e},{:e}],\"cells\":[{},{}],\"source\":\"negative accepted DG inlet volume flow\",\"kernel_support_mechanical_samples\":[{},{}],\"filter_delay_mechanical_samples\":{:e},\"reaction_added\":{}}}",
                    mouth.center_m[0], mouth.center_m[1], mouth.size_m(&geometry)[0], mouth.size_m(&geometry)[1],
                    mouth.cells[0], mouth.cells[1], microphone.delay_samples.0,
                    microphone.delay_samples.1, microphone.filter_delay_samples(), mouth_reaction.is_some() || mutual_mouth);
                let reaction_json = if mutual_mouth {
                    ",\"mouth_reaction\":{\"model\":\"joint-mixed-baffled-v1\",\"legacy_load_replaced\":true,\"memory\":\"one original MIMO mean-output ZOH\",\"mutual_mouth_terminal_reaction\":true,\"source_effort\":\"actual inlet pressure minus exterior mouth pressure\",\"work\":\"signed original joint-mouth exchange; not heat or certified rational storage\"}".to_owned()
                } else {
                    mouth_reaction.map_or_else(|| ",\"mouth_reaction\":null".to_owned(), |reaction| {
                        let spec = reaction.load();
                        let term = spec.terms()[0];
                        format!(",\"mouth_reaction\":{{\"model\":\"compact-rectangular-baffled-v1\",\"legacy_load_replaced\":true,\"maximum_frequency_hz\":{:e},\"resistance_pa_s_m3\":{:e},\"rate_per_s\":{:e},\"static_quadrature_relative_error\":{:e},\"max_complex_relative_error\":{:e},\"max_resistance_relative_error\":{:e},\"memory\":\"shared passive midpoint relaxation\",\"mutual_mouth_terminal_reaction\":false}}",
                            reaction.maximum_frequency_hz(), term.resistance_pa_s_m3, term.rate_per_s,
                            reaction.static_quadrature_relative_error(), reaction.max_complex_relative_error(),
                            reaction.max_resistance_relative_error())
                    })
                };
                let point = if source.apertures().len() == 1 { "baffled-mouth-and-outlet" }
                    else { "baffled-mouth-and-terminals" };
                let scope = if mutual_mouth {
                    "signed Rayleigh pressure sum of the actual rectangular mouth and circular terminals in one bound material exterior and infinite rigid baffle; complete same-step mutual mouth/terminal reaction; no resolved jet, free-field body or room radiation"
                } else if source.apertures().len() == 1 {
                    "signed Rayleigh pressure sum of the actual rectangular mouth and computed circular outlet in one material exterior and infinite rigid baffle; independent self loads, no mutual mouth/outlet reaction, resolved jet, free-field body or room radiation"
                } else {
                    "signed Rayleigh pressure sum of the actual rectangular mouth and computed circular terminals in one material exterior and infinite rigid baffle; no mutual reaction between mouth and terminals, resolved jet, free-field body or room radiation"
                };
                (point, format!(",\"observation_scope\":\"{scope}\",\"receiver_position_m\":[{:e},{:e},{:e}],\"receiver_band_hz\":{:e},\"exterior_gas\":{exterior},\"reconstruction\":\"same-clock-bandlimited-differentiating-Rayleigh\",\"kernel_support_mechanical_samples\":[{},{}],\"delay_compensated\":false{terminal_json}{mouth_json}{reaction_json}",
                    receiver.position_m[0], receiver.position_m[1], receiver.position_m[2],
                    receiver.maximum_frequency_hz, support.0, support.1))
            } else {
                let point = if load.terminal_ports().len() == 1 { "baffled-outlet" }
                    else { "baffled-terminals" };
                (point, super::duct_metadata::dg_baffled_observation_provenance(
                    load, receiver, terminal,
                ))
            }
        } else {
            ("mean-inlet", ",\"observation_scope\":\"internal interval-mean inlet pressure; not an exterior microphone\"".into())
        };
    let wall = duct.spatial_viscothermal_admission().map_or_else(
        || "\"uniform-consistent-mass\"".to_owned(),
        |wall| format!("{{\"model\":\"local-cylinder-positive-quadrature\",\"quadrature_points\":{},\"poles_per_field\":{},\"checked_max_complex_relative_error\":{:e},\"checked_max_real_loss_relative_error\":{:e}}}",
            wall.resolution.quadrature_points, wall.resolution.poles_per_field,
            wall.max_complex_error, wall.max_real_loss_error));
    let solid_wall = dg_wall_provenance(p.wall_layout());
    let pressure_observation = format!("{pressure_observation},\"observer_gas\":{},\"mouth_observer_gas\":{}",
        super::duct_metadata::observer_gas_provenance(p.renderer().observer_gas()),
        super::duct_metadata::observer_gas_provenance(p.renderer().mouth_observer_gas()));
    let edge_json = jet_speed_edge_provenance(voice);
    let source_json = format!("\"flue_dg_input\":{{\"schema\":\"frankensim-flue-dg-performance-v1\",\"blake3\":\"{}\",\"jet_card_source\":\"{}\",\"jet_card_hash\":\"{:016x}\",\"jet_card_authority\":\"{}\",\"compiled_controls\":{},\"propagator\":\"original-exact-dg\",\"operator_identity\":\"{}\",\"total_dofs\":{},\"length_m\":{:e},\"inlet_area_m2\":{:e},\"outlet_area_m2\":{:e},\"inlet_gas\":{inlet_gas},\"outlet_gas\":{outlet_gas},\"wall_band_rad_s\":{band},\"wall_owner\":{wall},\"flue_width_m\":{:e},\"flue_height_m\":{:e},\"cut_up_m\":{:e},\"pressure_point\":\"{pressure_point}\",\"model_scope\":\"original estimated jet and causal material-gas DG bore; solid wall image explicitly recorded; taper uses local-cylinder Webster approximation; source, spatial and nonlinear time convergence remain separate; no resolved turbulence or measured flute qualification\"{pressure_observation}{radiation}{field}{gas_sources}{solid_wall}{edge_json}}}",
        i.input_hash.to_hex(), i.jet_card_source, i.jet_card_hash, p.card().authority.label(),
        i.compiled_controls, duct.operator_identity().to_hex(), duct.total_dofs(),
        duct.total_length(), duct.inlet_area_m2(), duct.outlet_area_m2(),
        geometry.flue_width_m, geometry.flue_height_m, geometry.cut_up_m);
    source_json
}

pub(super) fn flue_provenance(p: &FlueFilePerformance) -> String {
    let i = p.info();
    let renderer = p.renderer();
    let network = renderer.voice().network().expect("file admits the physical network");
    let inlet = network.source().graph.inlet_section_index().expect("admitted unique inlet");
    let gas = network.source().gases[inlet];
    let exterior = renderer.radiation_input().expect("file binds exterior explicitly").exterior;
    let geometry = p.geometry();
    let island = renderer.voice().island();
    // Named only when the opt-in forced jet replaces the authored-gain chain.
    let jet_json = renderer.voice().forced_jet().map_or_else(String::new, |jet| format!(
        ",\"jet_response\":{{\"kind\":\"forced-bickley\",\"spacing_b\":{:e},\"window_b_over_u0\":{:e},\"half_width_m\":{:e},\"centre_over_bulk\":{:e},\"reynolds\":{:e}}}",
        jet.spec().spacing, jet.spec().window, jet.half_width_m(), jet.centre_over_bulk(), jet.reynolds()));
    // Named only when the opt-in feeds the card the exit's centre speed instead of the bulk speed.
    let card_speed_json = renderer.voice().card_centreline_speed().map_or_else(String::new, |speed| format!(
        ",\"jet_card_speed\":{{\"definition\":\"centreline\",\"law\":{},\"held_pressure_pa\":{:e},\"held_bulk_speed_m_s\":{:e},\"held_card_speed_m_s\":{:e}}}",
        json_string(speed.law), speed.held_pressure_pa, speed.held_bulk_speed_m_s, speed.held_card_speed_m_s));
    let edge_json = jet_speed_edge_provenance(renderer.voice());
    let windway_json = if let Some(flow) = renderer.voice().unsteady_windway() {
        let channel = flow.channel();
        let spec = channel.spec();
        format!(",\"windway_model\":\"unsteady-no-slip-slit-v1\",\"windway_length_m\":{:e},\"windway_dynamic_viscosity_pa_s\":{:e},\"windway_transverse_cells\":{},\"jet_convection_cells\":{}",
            spec.length_m, channel.gas().dynamic_viscosity, spec.transverse_cells, flow.convection_cells())
    } else if let Some((law, held)) = renderer.voice().developing_windway() {
        // The held note's validity measures; the attack and release are not covered.
        format!(",\"windway_model\":\"developing-slit-shah1978-v1\",\"windway_length_m\":{:e},\"windway_dynamic_viscosity_pa_s\":{:e},\"windway_held\":{{\"pressure_pa\":{:e},\"bulk_speed_m_s\":{:e},\"core_speed_m_s\":{:e},\"x_plus\":{:e},\"reynolds_hydraulic\":{:e},\"reynolds_length\":{:e},\"displacement_thickness_m\":{:e},\"reynolds_displacement\":{:e},\"chen_inlet_shift\":{:e},\"developed\":{},\"attack_covered\":false,\"geometry\":\"straight parallel plates; a converging windway declared as its exit gap is most likely overstated\"}}",
            law.length_m(), law.gas().dynamic_viscosity, held.pressure_pa, held.bulk_speed_m_s,
            held.core_speed_m_s, held.x_plus, held.reynolds_hydraulic, held.reynolds_length,
            held.displacement_thickness_m, held.reynolds_displacement, held.chen_inlet_shift, held.developed)
    } else {
        match renderer.voice().viscous_windway() {
            Some(law) => format!(",\"windway_model\":\"quasi-steady-slit-v1\",\"windway_length_m\":{:e},\"windway_dynamic_viscosity_pa_s\":{:e}", law.length_m(), law.dynamic_viscosity_pa_s()),
            None => ",\"windway_model\":\"ideal\",\"windway_length_m\":null".to_string(),
        }
    };
    let mouth_loss_json = match renderer.voice().mouth_contraction() {
        Some(loss) => format!(",\"mouth_contraction\":{{\"area_ratio\":{:e},\"maximum_contracted_mach\":{:e},\"source\":{},\"effective_area_m2\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e}}}",
            loss.spec().area_ratio, loss.spec().maximum_contracted_mach,
            json_string(&loss.spec().provenance), loss.effective_area_m2(),
            loss.gas().density, loss.gas().sound_speed),
        None => ",\"mouth_contraction\":null".to_string(),
    };
    let mouth_reaction = renderer.voice().mouth_radiation_load();
    let mouth_radiation_json = match renderer.mouth_radiation() {
        Some(mouth) => format!(",\"mouth_radiation\":{{\"shape\":\"rectangle\",\"center_m\":[{:e},{:e}],\"size_m\":[{:e},{:e}],\"cells\":[{},{}],\"source\":\"signed outward inlet volume flow\",\"reaction_added\":{}}}",
            mouth.center_m[0], mouth.center_m[1], mouth.size_m(&geometry)[0], mouth.size_m(&geometry)[1],
            mouth.cells[0], mouth.cells[1], mouth_reaction.is_some()),
        None => ",\"mouth_radiation\":null".to_string(),
    };
    let mouth_reaction_json = match mouth_reaction {
        Some(load) => {
            let spec = load.load();
            let term = spec.terms()[0];
            format!(",\"mouth_reaction\":{{\"model\":\"compact-rectangular-baffled-v1\",\"legacy_load_replaced\":true,\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e},\"maximum_frequency_hz\":{:e},\"resistance_pa_s_m3\":{:e},\"rate_per_s\":{:e},\"static_quadrature_relative_error\":{:e},\"max_complex_relative_error\":{:e},\"max_resistance_relative_error\":{:e},\"memory\":\"shared passive midpoint relaxation\",\"mutual_mouth_terminal_reaction\":false}}",
                load.exterior().density, load.exterior().sound_speed, load.maximum_frequency_hz(),
                term.resistance_pa_s_m3, term.rate_per_s, load.static_quadrature_relative_error(),
                load.max_complex_relative_error(), load.max_resistance_relative_error())
        }
        None => ",\"mouth_reaction\":null".to_string(),
    };
    let window_radiation = renderer.voice().window_radiation_load();
    let window_radiation_json = match window_radiation {
        Some(load) => {
            let cells = load.load().cells().iter().map(|c| format!("[{:e},{:e},{:e}]",
                c.resistance_pa_s_m3, c.inertance_pa_s2_m3, c.compliance_m3_pa)).collect::<Vec<_>>().join(",");
            let [band_lo, band_hi] = fs_couple::flue_loop::window_radiation::WINDOW_RADIATION_BAND_HZ;
            format!(",\"window_radiation\":{{\"model\":\"window-on-rigid-outer-cylinder-v1\",\"flux\":\"{}\",\"legacy_load_replaced\":\"mouth_loss_rel\",\"outer_radius_m\":{:e},\"opening_cells\":{},\"open_area_m2\":{:e},\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e},\"band_hz\":[{band_lo:e},{band_hi:e}],\"ka_band\":[{:e},{:e}],\"static_mass_kramers_kronig_pa_s2_m3\":{:e},\"static_mass_realized_pa_s2_m3\":{:e},\"static_mass_tail_share\":{:e},\"high_frequency_resistance_pa_s_m3\":{:e},\"band_error\":{:e},\"discrete_band_error\":{:e},\"tail_error\":{:e},\"tail_top_hz\":{:e},\"cells_r_l_c\":[{cells}],\"memory\":\"fs-vfit parallel-cell midpoint chain\"}}",
                load.flux().label(), load.outer_radius_m(), load.cells(), load.open_area_m2(),
                load.exterior().density, load.exterior().sound_speed, load.ka_band()[0], load.ka_band()[1],
                load.static_mass_pa_s2_m3(), load.load().static_inertance_pa_s2_m3(),
                load.static_mass_tail_share(), load.high_frequency_resistance_pa_s_m3(),
                load.band_error(), load.discrete_band_error(), load.tail_error(), load.tail_top_hz())
        }
        None => ",\"window_radiation\":null".to_string(),
    };
    let window_loss = if window_radiation.is_some() { "replaced by window_radiation" } else { "unchanged" };
    let window_json = match p.window() {
        Some(record) => {
            let (w, l) = (&record.window, &record.ladder);
            let pocket = match w.upstream {
                fs_bem::wall_window::WindowUpstream::ClosedEnd { distance }
                | fs_bem::wall_window::WindowUpstream::Port { distance } => distance,
            };
            let bore_area = core::f64::consts::PI * w.bore_radius * w.bore_radius;
            let fine = &l.rungs[l.rungs.len() - 1];
            let cheeks = match w.cheeks { fs_bem::wall_window::WindowCheeks::Radial => "radial", fs_bem::wall_window::WindowCheeks::Parallel => "parallel" };
            let ramp = w.ramp.map_or_else(|| "null".to_string(), |r| format!("{{\"edge_height_m\":{:e},\"angle_rad\":{:e}}}", r.edge_height, r.angle));
            let cover = w.cover.map_or_else(|| "null".to_string(), |c| match c.side {
                fs_bem::wall_window::CoverSide::Centred => format!("{{\"open_width_m\":{:e},\"side\":\"centred\",\"face_setback_m\":null}}", c.open_width),
                fs_bem::wall_window::CoverSide::OneSided { face } => format!("{{\"open_width_m\":{:e},\"side\":\"one-sided\",\"face_setback_m\":{}}}", c.open_width, face.map_or_else(|| "null".to_string(), |s| format!("{s:e}"))),
            });
            let rungs = l.rungs.iter().map(|r| format!("{:e}", r.bore_length)).collect::<Vec<_>>().join(",");
            let iterations = l.rungs.iter().map(|r| r.iterations.to_string()).collect::<Vec<_>>().join(",");
            let order_source = if w.is_prism() { "re-entrant edges" } else { "observed over three rungs" };
            format!(",\"window\":{{\"model\":\"wall-cut-potential-flow-v2\",\"legacy_load_replaced\":\"0.61 a_eq mouth inertance\",\"bore_radius_m\":{:e},\"wall_thickness_m\":{:e},\"arc_width_m\":{:e},\"resolved_arc_width_m\":{:e},\"axial_length_m\":{:e},\"cheeks\":\"{cheeks}\",\"ramp\":{ramp},\"cover\":{cover},\"resolved_outer_axial_length_m\":{:e},\"pocket_end_m\":{:e},\"cells\":[{},{},{}],\"rung_bore_lengths_m\":[{rungs}],\"rung_iterations\":[{iterations}],\"extrapolation_order\":{:e},\"order_source\":\"{order_source}\",\"bore_length_m\":{:e},\"network_channel_line_m\":{:e},\"riser_mass_per_density_per_m\":{:e},\"mouth_inertance_pa_s2_m3\":{:e},\"jet_source_coefficient\":\"unchanged\",\"mouth_loss_rel\":\"{window_loss}\"}}",
                w.bore_radius, w.wall_thickness, w.arc_width, fine.resolved_arc_width, w.axial_length,
                fine.outer_axial_length, pocket, record.cells[0], record.cells[1], record.cells[2], l.order,
                l.mass_per_density * bore_area, record.channel_line_m, record.riser_mass_per_density,
                renderer.voice().mouth_inertance_pa_s2_m3())
        }
        None => ",\"window\":null".to_string(),
    };
    let (phase, details) = network.fractional_sections().expect("file selects fractional travel");
    let actual_length: f64 = details.iter().map(|s|s.realized_dc_length_m).sum();
    let phase_error: f64 = details.iter().map(|s|s.phase_error_bound_rad).sum();
    let input_hash = i.input_hash.to_hex();
    let mut gas_sources = super::thermochemical::provenance(p.gas_records());
    if let Some(report) = renderer.observer_gas() {
        gas_sources.push_str(&format!(
            ",\"observer_gas\":{},\"mouth_observer_gas\":{}",
            super::duct_metadata::observer_gas_provenance(Some(report)),
            super::duct_metadata::observer_gas_provenance(renderer.mouth_observer_gas())
        ));
    }
    let receiver = renderer.config().receiver;
    let duct_json = flue_duct_provenance(renderer);
    let wall_json = flue_wall_provenance(p);
    let boundary_json = flue_boundary_provenance(p);
    let opening_loss_json = opening_loss_provenance(
        network.opening_losses().iter().map(|binding| ("network", binding)));
    let mimo_json = if let Some(record) = p.computed_radiation_record() {
        format!(",\"computed_radiation\":{}", super::duct_metadata::computed_baffled_provenance(record))
    } else {
        p.mimo_radiation_record().map_or_else(String::new, |record|
            format!(",\"mimo_radiation\":{}", super::duct_metadata::mimo_radiation_provenance(record)))
    };
    let mutual = p.mimo_radiation_record().is_some() || p.computed_radiation_record().is_some();
    let terminal_reaction = !p.radiation_bindings().is_empty() || mutual;
    let reaction = terminal_reaction || mouth_reaction.is_some() || window_radiation.is_some();
    let radiation_json = p.radiation_bindings().iter().map(|binding| format!(
        "{{\"node\":{},\"original_reflection\":{:e},\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e}{}}}",
        binding.opening().node,binding.previous_reflection(),binding.load().density_kg_m3(),
        binding.load().sound_speed_m_s(),radiation_provenance(Some(*binding.load()))
    )).collect::<Vec<_>>().join(",");
    let observation_scope = if p.computed_radiation_record().is_some() {
        "in-loop geometry/gas-derived mutual baffled terminal reaction and matched-medium ordered outlet pressure"
    } else if p.mimo_radiation_record().is_some() {
        "in-loop source-bound mutual baffled terminal reaction and matched-medium ordered outlet pressure"
    } else if terminal_reaction {
        "in-loop independent compact baffled terminal reaction and matched-medium outlet pressure"
    } else { "one-way infinite-baffle outlet pressure; supplied passive loads" };
    let radiation_exclusion = if mutual { "broadband unbaffled radiation" }
        else { "broadband unbaffled or mutual radiation" };
    let mouth_scope = if mouth_reaction.is_some() {
        "same-sample compact rectangular mouth reaction and matched exterior observation"
    } else if window_radiation.is_some() {
        "same-sample window radiation from the head's rigid outer cylinder; optional rectangular mouth is observation only"
    } else { "optional rectangular mouth is observation only" };
    let source_json = format!("\"flue_input\":{{\"schema\":\"{FLUE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\"jet_card_source\":\"{}\",\"jet_card_hash\":\"{:016x}\",\"jet_card_claim\":\"{}\",\"jet_card_authority\":\"{}\",\"composition_authority\":\"EstimateOnly\",\"compiled_controls\":{},\"represented_length_m\":{actual_length:e},\"flue_width_m\":{:e},\"flue_height_m\":{:e},\"cut_up_m\":{:e},\"labium_offset_m\":{:e},\"jet_angle_rad\":{:e},\"mouth_loss_rel\":{:e},\"authored_displacement_gain\":{:e}{jet_json}{card_speed_json}{edge_json},\"authored_seed_amplitude_rel\":{:e},\"gas_temperature_k\":{:e},\"gas_pressure_pa\":{:e},\"gas_density_kg_m3\":{:e},\"gas_sound_speed_m_s\":{:e},\"gas_dynamic_viscosity_pa_s\":{:e},\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e},\"viscothermal_sections\":{},\"phase_band_hz\":{:e},\"summed_phase_error_bound_rad\":{phase_error:e},\"receiver_position_m\":[{:e},{:e},{:e}],\"receiver_band_hz\":{:e},\"radiation_feedback_added\":{reaction},\"radiation_bindings\":[{radiation_json}],\"model_scope\":\"estimated tonal jet and supplied cylindrical graph; {observation_scope}; {mouth_scope}; wall image separately declared; no full shell dynamics, mouth-to-terminal mutual reaction, {radiation_exclusion}, resolved turbulence or measured flute/SPL claim\"{gas_sources}{windway_json}{mouth_loss_json}{mouth_radiation_json}{mouth_reaction_json}{window_json}{window_radiation_json}{duct_json}{wall_json}{boundary_json}{mimo_json}{opening_loss_json}}}",
        i.jet_card_source,i.jet_card_hash,p.card().claim.kind(),p.card().authority.label(),i.compiled_controls,
        geometry.flue_width_m,geometry.flue_height_m,geometry.cut_up_m,geometry.labium_offset_m,
        geometry.jet_angle_rad,geometry.mouth_loss_rel,island.displacement_gain,island.seed_amplitude_rel,
        gas.temperature,gas.pressure,gas.density,gas.sound_speed,gas.dynamic_viscosity,
        exterior.density,exterior.sound_speed,network.viscothermal_sections().len(),phase.maximum_frequency_hz,
        receiver.position_m[0],receiver.position_m[1],receiver.position_m[2],receiver.maximum_frequency_hz);
    source_json
}

fn load(bytes: &[u8], block: usize, policy: ObservationPolicy) -> Result<Loaded, String> {
    load_with_gas_sources(bytes, block, policy, None)
}
fn load_with_gas_sources(
    bytes: &[u8],
    block: usize,
    policy: ObservationPolicy,
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<Loaded, String> {
    load_with_sources(bytes, block, policy, sources,
        |_| Err("material-card packs require an explicit source loader".into()))
}
fn load_with_sources(
    bytes: &[u8],
    block: usize,
    policy: ObservationPolicy,
    sources: Option<&ThermochemicalSources<'_>>,
    load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
) -> Result<Loaded, String> {
    load_with_all_sources(bytes, block, policy, sources, load_material,
        |_| Err("MIMO radiation artifacts require an explicit source loader".into()), None)
}
fn load_with_all_sources(
    bytes: &[u8],
    block: usize,
    policy: ObservationPolicy,
    sources: Option<&ThermochemicalSources<'_>>,
    load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    load_radiation: impl FnMut(&str) -> Result<Vec<u8>, String>,
    flow_tap: Option<(Box<dyn Write + Send>, bool)>,
) -> Result<Loaded, String> {
    if flow_tap.is_some() && !bytes.starts_with(b"frankensim-flue-performance-v1\n") {
        return Err("--aperture-flows observes flue files only".into());
    }
    if bytes.starts_with(b"frankensim-resolved-dg-performance-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("resolved flow still requires source geometry, boundary and coupling convergence plus measured specimen qualification before --acceptance".into());
        }
        let p = fs_couple::nsf_dg::performance::file::ResolvedDgFilePerformance::from_bytes_with_sources(
            bytes, block, std::sync::Arc::new(CancelGate::new()), sources, load_material,
        ).map_err(|e| e.to_string())?;
        let info = p.info();
        let source_json = resolved_provenance(&p);
        return Ok(Loaded {
            source: Box::new(p.into_renderer()), rate: info.sample_rate_hz, samples: info.samples,
            full_scale_pa: info.full_scale_pa, hash: info.input_hash,
            fixture: "resolved-dg-input", source_json, observation_json: None,
        });
    }
    if bytes.starts_with(b"frankensim-material-modal-dg-valve-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("material-modal DG files retain flat DKT plates and a finite modal interface; exact propagation does not qualify a measured instrument".into());
        }
        let p = fs_couple::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance::from_bytes_with_sources(
            bytes, block, &CancelGate::new(), sources, load_material,
        ).map_err(|e| e.to_string())?;
        let i = p.info();
        let source_json = material_modal_dg_provenance(&p);
        return Ok(Loaded {
            source: Box::new(p.into_renderer()), rate: i.output_rate_hz,
            samples: i.samples, full_scale_pa: i.full_scale_pa, hash: i.input_hash,
            fixture: "material-modal-dg-input", source_json, observation_json: None,
        });
    }
    if bytes.starts_with(b"frankensim-material-modal-valve-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("material-modal files retain flat DKT walls and a finite modal interface; matched baffled observation is not a measured instrument qualification".into());
        }
        let p = fs_couple::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance::from_bytes_with_all_sources(
            bytes, block, &CancelGate::new(), sources, load_material, load_radiation).map_err(|e| e.to_string())?;
        let i = p.info();
        if i.output_rate_hz != RATE {
            return Err("material-modal wind output must declare 48000 Hz; mechanical refinement belongs to its separate source clock".into());
        }
        let source_json = material_modal_provenance(&p);
        let ratio = i.mechanics_rate_hz / i.output_rate_hz;
        let observation_json = format!("{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"output_sample_rate_hz\":{},\"ratio\":{ratio},\"filter\":\"shared multirate differentiating Rayleigh receiver with bandlimited propagation\",\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
            i.mechanics_rate_hz, i.samples * u64::from(ratio), i.output_rate_hz);
        return Ok(Loaded { source: Box::new(p.into_renderer()), rate: i.output_rate_hz,
            samples: i.samples, full_scale_pa: i.full_scale_pa, hash: i.input_hash,
            fixture: "material-modal-valve-input", source_json, observation_json: Some(observation_json) });
    }
    // Exact versioned source header selects the existing owner; filenames and
    // instrument labels never choose geometry or material constants.
    // The exact DG source reuses the same jet and material input owners. Its
    // selected same-clock pressure follows the explicit decimation path below.
    if bytes.starts_with(b"frankensim-flue-dg-performance-v1\n") {
        use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
        if policy == ObservationPolicy::Acceptance {
            return Err("DG flue files retain an exploration-only jet card; --acceptance requires resolved qualified source flow".into());
        }
        if flow_tap.is_some() {
            return Err("DG flue aperture-flow export is not supported by this pressure input".into());
        }
        let p = FlueDgFilePerformance::from_bytes_with_sources(
            bytes, block, &CancelGate::new(), sources, load_material).map_err(|e| e.to_string())?;
        let i = p.info();
        let source_json = flue_dg_provenance(&p);
        return Ok(Loaded { source: Box::new(p.into_renderer()), rate: i.sample_rate_hz,
            samples: i.samples, full_scale_pa: i.full_scale_pa, hash: i.input_hash,
            fixture: "flue-dg-input", source_json, observation_json: None });
    }
    if bytes.starts_with(b"frankensim-flue-performance-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("flue files retain an explicitly estimated tonal-jet/cylindrical-graph image; --acceptance requires an independently qualified physical model".into());
        }
        let p = FlueFilePerformance::from_bytes_with_all_sources(bytes, block,
            &CancelGate::new(), sources, load_material, load_radiation).map_err(|e| e.to_string())?;
        let i = p.info();
        if i.sample_rate_hz != RATE {
            return Err("flue wind output must declare 48000 Hz; its file separately declares the faster mechanical clock".into());
        }
        let renderer = p.renderer();
        let source_json = flue_provenance(&p);
        let ratio = renderer.mechanics_per_output();
        let mechanical_samples = i.samples*ratio as u64; // already checked by the performance owner
        let propagation = renderer.propagation_delay_samples();
        let observation_json = format!("{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{mechanical_samples},\"output_sample_rate_hz\":{RATE},\"ratio\":{ratio},\"filter\":\"shared multirate differentiating Rayleigh receiver with bandlimited propagation\",\"delay_output_samples\":{},\"propagation_fir_support_output_samples\":[{},{}],\"first_output_source_index\":{},\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
            i.mechanics_rate_hz,renderer.receiver().filter_delay_samples(),propagation.0,
            propagation.1,ratio-1);
        let renderer = match flow_tap {
            Some((sink, false)) => p.into_renderer().with_aperture_flow_tap(sink).map_err(|e| e.to_string())?,
            Some((sink, true)) => p.into_renderer().with_aperture_flow_and_labium_tap(sink)
                .map_err(|e| e.to_string())?,
            None => p.into_renderer(),
        };
        Ok(Loaded { source:Box::new(renderer),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"flue-input",source_json,
            observation_json:Some(observation_json) })
    } else if bytes.starts_with(b"frankensim-plate-valve-performance-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("--acceptance selects the explicit ReedPerformance flanged-development path; plate-valve files require explicit matched-load builder admission".into());
        }
        let p=PlateValvePerformance::from_bytes_with_all_sources(bytes,block,&CancelGate::new(),sources,load_material).map_err(|e|e.to_string())?;
        let i=p.info();let a=p.renderer().system().aperture();
        let plate=a.plate_reduction().expect("plate file retains its specimen");
        let observation=match p.renderer().observation() {
            ApertureObservation::Inlet=>"inlet",
            ApertureObservation::TubeTerminal=>"terminal",
            ApertureObservation::NetworkNode(node)=>match p.renderer().system() {
                CoupledAperture::Network(n) if i.duct_sections==1 && matches!(n.spec().nodes[node],
                    fs_couple::bernoulli_aperture::network::NetworkNode::Termination {..}
                    | fs_couple::bernoulli_aperture::network::NetworkNode::Impedance {..}
                    | fs_couple::bernoulli_aperture::network::NetworkNode::Relaxation {..})=>"terminal",
                _=>"network-node",
            },
            ApertureObservation::TubeBaffled(_)|ApertureObservation::NetworkBaffled {..}=>"baffled-outlet",
        };
        let receiver_json=format!("{}{}",outlet_provenance(p.renderer(),i.radiation_load),super::duct_metadata::provenance(&p));
        let wall_json=p.dg_wall_layout().map_or_else(String::new,dg_wall_provenance);
        let computed_dg = matches!(p.renderer().system(), CoupledAperture::Dg(runtime)
            if runtime.load().radiation_record().is_some());
        let scope=if observation!="baffled-outlet" { "internal pressure, not an exterior microphone" }
            else if computed_dg { "exterior baffled-outlet pressure from the accepted causal DG terminal flow; computed full-clock-band reaction with separate continuous and sampled GridOnly checks; no measured-instrument qualification" }
            else if i.radiation_load.is_some() { "exterior baffled-outlet pressure with compact passive radiation feedback; not broadband matched radiation or measured calibration" }
            else { "one-way exterior baffled-outlet pressure; no matched radiation load or measured-instrument claim" };
        let requested=match p.renderer().system() {
            CoupledAperture::Dg(d)=>d.load().duct().total_length(),
            CoupledAperture::Tube(t)=>t.spec().length_m,
            CoupledAperture::Network(n)=>n.spec().sections.iter().map(|s|s.length_m).sum(),
        };
        let propagation=if matches!(p.renderer().system(),CoupledAperture::Dg(_)) {"gas-derived exact-Bessel viscous/thermal memory in the original DG exponential"}
            else if p.viscothermal_losses().is_empty() {"lossless propagating sections with explicit local loads"}
            else {"gas-derived viscothermal series/thermal memory distributed along selected sections"};
        let gas_json=format!("{}{}{}", super::thermochemical::provenance(p.gas_records()),
            finite_gas_heat::provenance(p.renderer(), i.sample_rate_hz),
            finite_gas_heat::wall_provenance(p.renderer(), i.sample_rate_hz));
        // Each section's pack, hash, state and query point, or its fixture label (pexfp.32).
        let section_sources=p.section_provenance().iter().map(|s|format!("\"{s}\"")).collect::<Vec<_>>().join(",");
        let source_json=format!("\"plate_valve_input\":{{\"schema\":\"{PLATE_VALVE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{}\",\"nodes\":{},\"triangles\":{},\"sections\":{},\"section_sources\":[{section_sources}],\"memory_branches\":{},\"compiled_controls\":{},\"pressure_point\":\"{observation}\",\"requested_tube_length_m\":{:e},\"represented_tube_length_m\":{:e},\"effective_mass_kg\":{:e},\"stiffness_n_m\":{:e},\"pressure_area_m2\":{:e},\"model_scope\":\"one linear plate mode; spatial lay and supplied material history; {propagation}; {scope}\"{receiver_json}{gas_json}{wall_json}}}",
            i.input_hash.to_hex(),i.nodes,i.triangles,i.sections,i.memory_branches,i.compiled_controls,
            requested,i.represented_tube_length_m,plate.mass_kg(),plate.stiffness_n_m(),plate.pressure_area_m2());
        Ok(Loaded{source:Box::new(p.into_renderer()),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"plate-valve-input",source_json,observation_json:None})
    } else {
        let p=ReedPerformance::from_bytes_with_gas_sources(bytes,block,policy,sources).map_err(|e|e.to_string())?;
        let i=p.info();
        let observer_receipt=super::duct_metadata::reed_observer_json(&p).map_err(|e|e.to_string())?;
        let observer_json=format!("{}{}",super::duct_metadata::reed_provenance(&p),super::thermochemical::provenance(p.gas_records()));
        let slit_json = p.viscous_slit().map_or_else(String::new, |slit| {
            let gas = slit.gas().expect("reed channel retains its actual inlet gas");
            format!(",\"slit_channel\":{{\"model\":\"quasisteady parallel-plate resistance plus separated jet\",\"length_m\":{:e},\"inlet_temperature_k\":{:e},\"inlet_pressure_pa\":{:e},\"inlet_density_kg_m3\":{:e},\"dynamic_viscosity_pa_s\":{:e}}}",
                slit.length_m(), gas.temperature, gas.pressure, gas.density, gas.dynamic_viscosity)
        });
        let policy_name=if policy==ObservationPolicy::Acceptance {"acceptance"} else {"diagnostic"};
        let source_json=format!("\"reed_input\":{{\"schema\":\"{REED_PERFORMANCE_SCHEMA}\",\"observation_policy\":\"{policy_name}\",\"blake3\":\"{}\",\"segments\":{},\"tone_holes\":{},\"massive_reed\":{},\"gesture_events\":{},\"compiled_controls\":{}{observer_json}{slit_json}{observer_receipt}}}",
            i.input_hash.to_hex(),i.segments,i.tone_holes,i.massive_reed,i.gesture_events,i.compiled_controls);
        Ok(Loaded{source:Box::new(p),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"reed-input",source_json,observation_json:None})
    }
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    run_with_gas_sources(args, None)
}

/// Native file options use the same snapshots for gas, solids and radiation.
pub(super) fn run_with_gas_options(
    args: &[String],
    gas_options: &super::thermochemical::Options,
) -> Result<(), String> {
    gas_sources::run_with_gas_options(args, gas_options, None)
}

pub(super) fn run_resolved_with_gas_options(
    args: &[String],
    gas_options: &super::thermochemical::Options,
) -> Result<(), String> {
    if options(args)?.5.is_some() {
        return Err("resolved flow does not provide the flue aperture-flow export".into());
    }
    gas_sources::run_with_gas_options(args, gas_options,
        Some(fs_couple::nsf_dg::performance::file::RESOLVED_DG_PERFORMANCE_SCHEMA))
}

pub(super) fn run_with_gas_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    run_with_prepared_sources(args, sources, None, None)
}

fn run_with_prepared_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
    prepared: Option<attachments::AttachedSources>,
    required_schema: Option<&'static str>,
) -> Result<(), String> {
    let (input, output, block, decimate, policy, flows, self_contained, attached_names) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let flow_tap: Option<(Box<dyn Write + Send>, bool)> = match flows.map(|(p, l)| (Path::new(p), l)) {
        Some((path, _)) if path == output || path == sidecar.as_path() => {
            return Err("--aperture-flows must name a file distinct from the WAV and its sidecar".into());
        }
        Some((path, labium)) => Some((Box::new(std::io::BufWriter::new(std::fs::File::create_new(path)
            .map_err(|e| format!("aperture flow file {} create failed: {e}", path.display()))?)), labium)),
        None => None,
    };
    let file = std::fs::File::open(input).map_err(|e| format!("wind input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_PLATE_VALVE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("wind input read failed: {e}"))?;
    if let Some(schema) = required_schema {
        if !bytes.strip_prefix(schema.as_bytes()).is_some_and(|rest| rest.starts_with(b"\n")) {
            return Err("resolved command requires frankensim-resolved-dg-performance-v1 source bytes".into());
        }
    }
    let parent = Path::new(input).parent().unwrap_or(Path::new("."));
    let attached_sources = match prepared {
        Some(snapshot) => snapshot,
        None => attachments::AttachedSources::load(parent, &attached_names)?,
    };
    let loaded = load_with_all_sources(&bytes, block, policy, sources, |name| {
        // Restricted callbacks can only read eagerly admitted input snapshots.
        // A missing name never falls through to ambient filesystem access.
        if self_contained {
            let pack_bytes = attached_sources.get(name).ok_or_else(||
                format!("--self-contained refuses external material-card pack {name:?}"))?;
            return NormalizedMaterialCardPack::from_bytes(pack_bytes)
                .map_err(|e| format!("material-card pack {name:?}: {e}"));
        }
        let file = std::fs::File::open(parent.join(name))
            .map_err(|e| format!("material-card pack {name:?} open failed: {e}"))?;
        let mut pack_bytes = Vec::new();
        file.take((MAX_FLUE_WALL_PACK_BYTES + 1) as u64).read_to_end(&mut pack_bytes)
            .map_err(|e| format!("material-card pack {name:?} read failed: {e}"))?;
        if pack_bytes.len() > MAX_FLUE_WALL_PACK_BYTES {
            return Err("material-card pack exceeds 8 MiB".into());
        }
        NormalizedMaterialCardPack::from_bytes(&pack_bytes)
            .map_err(|e| format!("material-card pack {name:?}: {e}"))
    }, |name| {
        if self_contained {
            return attached_sources.get(name).map(<[u8]>::to_vec).ok_or_else(||
                format!("--self-contained refuses external radiation artifact {name:?}"));
        }
        let file = std::fs::File::open(parent.join(name))
            .map_err(|e| format!("MIMO radiation artifact {name:?} open failed: {e}"))?;
        let mut bytes = Vec::new();
        file.take((MAX_MIMO_RADIATION_BYTES + 1) as u64).read_to_end(&mut bytes)
            .map_err(|e| format!("MIMO radiation artifact {name:?} read failed: {e}"))?;
        if bytes.len() > MAX_MIMO_RADIATION_BYTES {
            return Err("MIMO radiation artifact exceeds 8 MiB".into());
        }
        Ok(bytes)
    }, flow_tap)?;
    if decimate && loaded.observation_json.is_some() {
        return Err("this source already declares mechanical and output clocks; omit --decimate to avoid a second observer".into());
    }
    if !decimate && loaded.rate != RATE {
        return Err("non-48000-Hz mechanics require explicit --decimate; no implicit resampling".into());
    }
    if decimate && loaded.rate <= RATE {
        return Err("--decimate requires mechanics above 48000 Hz; ordinary 48000-Hz input needs no conversion".into());
    }
    // The entire finite physical runtime stays inside this observer, including
    // spatial contact and material memory. Only its output pressure is filtered.
    let mut renderer = DecimatedRenderer::new(loaded.source, loaded.rate, RATE, block)
        .map_err(|e| e.to_string())?;
    let samples = renderer.output_samples_for(loaded.samples).map_err(|e| e.to_string())?;
    if samples == 0 { return Err("wind source needs at least one complete output interval".into()); }
    renderer.validate_sample_count(samples).map_err(|e| e.to_string())?;
    let samples = usize::try_from(samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let observation = renderer.info();
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, samples, block, loaded.full_scale_pa,
    )?;
    let input_hash = loaded.hash.to_hex();
    let wav_hash = rendered.hash.to_hex();
    let fixture=loaded.fixture;
    let source_clock = if fixture == "resolved-dg-input" { "observation" } else { "mechanics" };
    let observation_json = loaded.observation_json.unwrap_or_else(|| format!(
        "{{\"{source_clock}_sample_rate_hz\":{},\"{source_clock}_samples\":{},\"output_sample_rate_hz\":{RATE},\"ratio\":{},\"filter\":\"{}\",\"delay_output_samples\":{},\"first_output_source_index\":{},\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
        loaded.rate,loaded.samples,observation.ratio,observation.filter_profile,
        observation.delay_output_samples,observation.first_output_source_index));
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"{fixture}\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",{},\
         \"observation\":{observation_json}}}",
        loaded.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,loaded.source_json,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!("{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"{fixture}\",\
        \"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
        \"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa, rendered.rms_pa);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wind_observation_policy_is_explicit_and_duplicate_flag_refuses() {
        let args = vec!["source.performance".into(), "out.wav".into()];
        assert_eq!(options(&args).unwrap().4, ObservationPolicy::Diagnostic);
        let mut accepted = args;
        accepted.push("--acceptance".into());
        assert_eq!(options(&accepted).unwrap().4, ObservationPolicy::Acceptance);
        accepted.push("--acceptance".into());
        assert!(options(&accepted).is_err());
    }

    #[test]
    fn acceptance_rejects_legacy_listener_before_creating_audio() {
        let physical = include_str!("../../../examples/reed-flanged-baffled.performance");
        let legacy = physical.replace("baffled_outlet 0 0 0.5 4 16 1000", "listener_m 0.5");
        let error = load(legacy.as_bytes(), 37, ObservationPolicy::Acceptance)
            .err()
            .unwrap();
        assert!(error.contains("diagnostic-only"), "{error}");
        assert!(load(physical.as_bytes(), 37, ObservationPolicy::Acceptance).is_ok());
    }
}

