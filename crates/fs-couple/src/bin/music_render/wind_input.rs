//! File-driven wind sources through the shared pressure and WAV owners.
use std::io::{Read, Write};
use std::path::Path;
use fs_blake3::ContentHash;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, ObservationPolicy, PressureRenderer};
use fs_couple::render::schedule::reed::{ReedPerformance, REED_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::performance::{ApertureObservation, AperturePerformance, CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::{PlateValvePerformance, MAX_PLATE_VALVE_PERFORMANCE_BYTES, PLATE_VALVE_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::flue_loop::network::performance::file::{
    FlueFilePerformance, FLUE_PERFORMANCE_SCHEMA, MAX_FLUE_WALL_PACK_BYTES,
};
use fs_matdb::NormalizedMaterialCardPack;
use fs_exec::CancelGate;
use fs_material::gas_species::input::ThermochemicalSources;
use super::{RATE, create_outputs, json_string, stream_output};

fn options(args: &[String]) -> Result<(&str, &str, usize, bool, ObservationPolicy), String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut decimate = false;
    let mut policy = ObservationPolicy::Diagnostic;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--acceptance" => {
                if policy == ObservationPolicy::Acceptance { return Err("--acceptance may only be supplied once".into()); }
                policy = ObservationPolicy::Acceptance;
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
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render wind INPUT.performance OUT.wav [--block N] [--decimate] [--acceptance]".into());
    };
    Ok((*input, *output, block.unwrap_or(512), decimate, policy))
}

// A physical propagation delay is source geometry, not an output/filter latency
// to be removed by ensemble alignment. Keep it on the mechanical sample clock.
pub(super) fn outlet_provenance(p: &AperturePerformance, radiation: Option<BaffledRadiationLoad>) -> String {
    let load_json=radiation_provenance(radiation);
    let (c,radius,_density,_speed)=match (p.observation(),p.system()) {
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

fn flue_wall_provenance(p: &FlueFilePerformance) -> String {
    use fs_couple::bernoulli_aperture::wall::material::{AxialConstraint, WallMaterialLoss};
    let walls = p.renderer().voice().network().expect("file admits graph").material_walls();
    let records = p.wall_records().iter().map(|record| {
        let wall = walls.iter().find(|wall| wall.source().node == record.node())
            .expect("retained material source has its original physical wall");
        let assignment = wall.source();
        let material = &assignment.material;
        let state = material.material();
        let coordinate = |name: &str| material.query_point().iter()
            .find(|(axis, _)| axis == name).expect("file supplies T and omega").1;
        let axial = match assignment.axial {
            AxialConstraint::StressFree => "stress-free",
            AxialConstraint::StrainFixed => "strain-fixed",
        };
        let loss = match &assignment.loss {
            WallMaterialLoss::Elastic => "elastic",
            WallMaterialLoss::KelvinVoigt => "kelvin-voigt",
            WallMaterialLoss::Maxwell(_) => unreachable!("file profile selects elastic or Kelvin-Voigt"),
        };
        let patch = wall.patch();
        format!("{{\"node\":{},\"model\":\"isotropic-local-ring-membrane\",\"pack_path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"resolved_state_hash\":\"{}\",\"source_artifact\":\"{}\",\"material_chemistry\":{},\"material_phase\":{},\"material_process\":{},\"material_revision\":{},\"temperature_k\":{:e},\"angular_frequency_rad_s\":{:e},\"band_hz\":[{:e},{:e}],\"radius_m\":{:e},\"covered_length_m\":{:e},\"area_m2\":{:e},\"thickness_m\":{:e},\"axial_constraint\":\"{axial}\",\"loss\":\"{loss}\",\"surface_density_kg_m2\":{:e},\"stiffness_pa_m\":{:e},\"resistance_pa_s_m\":{:e}}}",
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

struct Loaded {
    source: Box<dyn PressureRenderer>, rate:u32, samples:u64, full_scale_pa:f64,
    hash:ContentHash, fixture:&'static str, source_json:String,
    // A source that already owns multirate physical observation reports THAT
    // clock/filter history; the identity PCM wrapper must not call it 48k mechanics.
    observation_json:Option<String>,
}
fn flue_duct_provenance(p: &fs_couple::flue_loop::network::performance::FlueOutletPerformance) -> String {
    use fs_couple::bernoulli_aperture::network::NetworkNode;
    let source = p.voice().network().expect("physical flue network").source();
    let inlet = source.graph.inlet_section_index().expect("unique inlet");
    let nodes = source.graph.nodes.iter().map(|node| match node {
        NetworkNode::Inlet => "{\"kind\":\"inlet\"}".to_string(),
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
    if bytes.starts_with(b"frankensim-material-modal-valve-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("material-modal files retain flat DKT walls and an effective uniform slit; matched baffled observation is not a measured instrument qualification".into());
        }
        let p = fs_couple::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance::from_bytes_with_sources(
            bytes, block, &CancelGate::new(), sources, load_material).map_err(|e| e.to_string())?;
        let i = p.info();
        if i.output_rate_hz != RATE {
            return Err("material-modal wind output must declare 48000 Hz; mechanical refinement belongs to its separate source clock".into());
        }
        let valve = match p.renderer().source().system().valve() {
            fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(a) => a,
            _ => unreachable!("material-modal file owns the requested modal valve"),
        };
        let materials = p.material_records().iter().map(|record| {
            let source = valve.material_plates()[record.blade()].source();
            let region = &source.regions()[record.section()];
            let material = &region.input.material;
            format!("{{\"blade\":{},\"section\":{},\"path\":{},\"pack_hash\":\"{}\",\"card_hash\":\"{}\",\"state_hash\":\"{}\",\"specimen_hash\":\"{}\",\"mass_kg\":{:e},\"thickness_m\":{:e},\"material_angle_rad\":{:e}}}",
                record.blade(), record.section(), json_string(record.path()),
                record.pack().content_hash().to_hex(), material.card_identity().to_hex(),
                material.identity().to_hex(), source.specimen_identity().to_hex(),
                region.mass_kg, region.section.thickness, region.input.material_angle_rad)
        }).collect::<Vec<_>>().join(",");
        let gas_sources = super::thermochemical::provenance(p.gas_records());
        let system = p.renderer().source().system();
        let supply = system.supply_channel().map_or_else(|| "null".into(), |channel| {
            let spec = channel.spec();
            format!("{{\"length_m\":{:e},\"width_m\":{:e},\"height_m\":{:e},\"transverse_cells\":{},\"density_kg_m3\":{:e},\"viscosity_pa_s\":{:e}}}",
                spec.length_m, spec.width_m, spec.height_m, spec.transverse_cells,
                channel.gas().density, channel.gas().dynamic_viscosity)
        });
        let source_json = format!("\"material_modal_valve_input\":{{\"schema\":\"frankensim-material-modal-valve-v1\",\"blake3\":\"{}\",\"blades\":{},\"retained_modes\":{},\"nodes\":{},\"triangles\":{},\"compiled_controls\":{},\"full_tensor_viscosity\":{},\"radiation_terminals\":{},\"material_sources\":[{materials}],\"supply_channel\":{supply},\"model_scope\":\"regional material flat DKT walls; original signed pressure and force projections; effective uniform slit/contact; both gas-bound ducts; compact independent baffled reaction and shared exterior pressure; no measured reed, unbaffled or resolved moving-channel qualification\"{gas_sources}}}",
            i.input_hash.to_hex(), i.blades, i.modes, i.nodes, i.triangles,
            p.renderer().pending_controls().len(), valve.viscosity_matrix_n_s_m().is_some(),
            p.radiation_bindings().len());
        let ratio = i.mechanics_rate_hz / i.output_rate_hz;
        let observation_json = format!("{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"output_sample_rate_hz\":{},\"ratio\":{ratio},\"filter\":\"shared multirate differentiating Rayleigh receiver with bandlimited propagation\",\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
            i.mechanics_rate_hz, i.samples * u64::from(ratio), i.output_rate_hz);
        return Ok(Loaded { source: Box::new(p.into_renderer()), rate: i.output_rate_hz,
            samples: i.samples, full_scale_pa: i.full_scale_pa, hash: i.input_hash,
            fixture: "material-modal-valve-input", source_json, observation_json: Some(observation_json) });
    }
    // Exact versioned source header selects the existing owner; filenames and
    // instrument labels never choose geometry or material constants.
    if bytes.starts_with(b"frankensim-flue-performance-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("flue files retain an explicitly estimated tonal-jet/cylindrical-graph image; --acceptance requires an independently qualified physical model".into());
        }
        let p = FlueFilePerformance::from_bytes_with_sources(bytes, block,
            &CancelGate::new(), sources, load_material).map_err(|e| e.to_string())?;
        let i = p.info();
        if i.sample_rate_hz != RATE {
            return Err("flue wind output must declare 48000 Hz; its file separately declares the faster mechanical clock".into());
        }
        let renderer = p.renderer();
        let network = renderer.voice().network().expect("file admits the physical network");
        let inlet = network.source().graph.inlet_section_index().expect("admitted unique inlet");
        let gas = network.source().gases[inlet];
        let exterior = renderer.radiation_input().expect("file binds exterior explicitly").exterior;
        let geometry = p.geometry();
        let island = renderer.voice().island();
        let windway_json = if let Some(flow) = renderer.voice().unsteady_windway() {
            let channel = flow.channel();
            let spec = channel.spec();
            format!(",\"windway_model\":\"unsteady-no-slip-slit-v1\",\"windway_length_m\":{:e},\"windway_dynamic_viscosity_pa_s\":{:e},\"windway_transverse_cells\":{},\"jet_convection_cells\":{}",
                spec.length_m, channel.gas().dynamic_viscosity, spec.transverse_cells, flow.convection_cells())
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
        let (phase, details) = network.fractional_sections().expect("file selects fractional travel");
        let actual_length: f64 = details.iter().map(|s|s.realized_dc_length_m).sum();
        let phase_error: f64 = details.iter().map(|s|s.phase_error_bound_rad).sum();
        let input_hash = i.input_hash.to_hex();
        let gas_sources = super::thermochemical::provenance(p.gas_records());
        let receiver = renderer.config().receiver;
        let duct_json = flue_duct_provenance(renderer);
        let wall_json = flue_wall_provenance(&p);
        let reaction = !p.radiation_bindings().is_empty();
        let radiation_json = p.radiation_bindings().iter().map(|binding| format!(
            "{{\"node\":{},\"original_reflection\":{:e},\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e}{}}}",
            binding.opening().node,binding.previous_reflection(),binding.load().density_kg_m3(),
            binding.load().sound_speed_m_s(),radiation_provenance(Some(*binding.load()))
        )).collect::<Vec<_>>().join(",");
        let observation_scope = if reaction {
            "in-loop independent compact baffled terminal reaction and matched-medium outlet pressure"
        } else { "one-way infinite-baffle outlet pressure; supplied passive loads" };
        let source_json = format!("\"flue_input\":{{\"schema\":\"{FLUE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\"jet_card_source\":\"{}\",\"jet_card_hash\":\"{:016x}\",\"jet_card_claim\":\"{}\",\"jet_card_authority\":\"{}\",\"composition_authority\":\"EstimateOnly\",\"compiled_controls\":{},\"represented_length_m\":{actual_length:e},\"flue_width_m\":{:e},\"flue_height_m\":{:e},\"cut_up_m\":{:e},\"labium_offset_m\":{:e},\"jet_angle_rad\":{:e},\"mouth_loss_rel\":{:e},\"authored_displacement_gain\":{:e},\"authored_seed_amplitude_rel\":{:e},\"gas_temperature_k\":{:e},\"gas_pressure_pa\":{:e},\"gas_density_kg_m3\":{:e},\"gas_sound_speed_m_s\":{:e},\"gas_dynamic_viscosity_pa_s\":{:e},\"exterior_density_kg_m3\":{:e},\"exterior_sound_speed_m_s\":{:e},\"viscothermal_sections\":{},\"phase_band_hz\":{:e},\"summed_phase_error_bound_rad\":{phase_error:e},\"receiver_position_m\":[{:e},{:e},{:e}],\"receiver_band_hz\":{:e},\"radiation_feedback_added\":{reaction},\"radiation_bindings\":[{radiation_json}],\"model_scope\":\"estimated tonal jet and supplied cylindrical graph; {observation_scope}; wall image separately declared; no full shell dynamics, broadband unbaffled or mutual radiation, resolved turbulence or measured flute/SPL claim\"{gas_sources}{windway_json}{mouth_loss_json}{duct_json}{wall_json}}}",
            i.jet_card_source,i.jet_card_hash,p.card().claim.kind(),p.card().authority.label(),i.compiled_controls,
            geometry.flue_width_m,geometry.flue_height_m,geometry.cut_up_m,geometry.labium_offset_m,
            geometry.jet_angle_rad,geometry.mouth_loss_rel,island.displacement_gain,island.seed_amplitude_rel,
            gas.temperature,gas.pressure,gas.density,gas.sound_speed,gas.dynamic_viscosity,
            exterior.density,exterior.sound_speed,network.viscothermal_sections().len(),phase.maximum_frequency_hz,
            receiver.position_m[0],receiver.position_m[1],receiver.position_m[2],receiver.maximum_frequency_hz);
        let ratio = renderer.mechanics_per_output();
        let mechanical_samples = i.samples*ratio as u64; // already checked by the performance owner
        let observation_json = format!("{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{mechanical_samples},\"output_sample_rate_hz\":{RATE},\"ratio\":{ratio},\"filter\":\"shared multirate differentiating Rayleigh receiver with bandlimited propagation\",\"delay_output_samples\":{},\"propagation_fir_support_output_samples\":[{},{}],\"first_output_source_index\":{},\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
            i.mechanics_rate_hz,renderer.receiver().filter_delay_samples(),renderer.receiver().delay_samples.0,
            renderer.receiver().delay_samples.1,ratio-1);
        Ok(Loaded { source:Box::new(p.into_renderer()),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"flue-input",source_json,
            observation_json:Some(observation_json) })
    } else if bytes.starts_with(b"frankensim-plate-valve-performance-v1\n") {
        if policy == ObservationPolicy::Acceptance {
            return Err("--acceptance selects the explicit ReedPerformance flanged-development path; plate-valve files require explicit matched-load builder admission".into());
        }
        let p=PlateValvePerformance::from_bytes_with_gas_sources(bytes,block,&CancelGate::new(),sources).map_err(|e|e.to_string())?;
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
        let scope=if observation!="baffled-outlet" { "internal pressure, not an exterior microphone" }
            else if i.radiation_load.is_some() { "exterior baffled-outlet pressure with compact passive radiation feedback; not broadband matched radiation or measured calibration" }
            else { "one-way exterior baffled-outlet pressure; no matched radiation load or measured-instrument claim" };
        let requested=match p.renderer().system() {
            CoupledAperture::Tube(t)=>t.spec().length_m,
            CoupledAperture::Network(n)=>n.spec().sections.iter().map(|s|s.length_m).sum(),
        };
        let propagation=if p.viscothermal_losses().is_empty() {"lossless propagating sections with explicit local loads"}
            else {"gas-derived viscothermal series/thermal memory distributed along selected sections"};
        let gas_json=super::thermochemical::provenance(p.gas_records());
        let source_json=format!("\"plate_valve_input\":{{\"schema\":\"{PLATE_VALVE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{}\",\"nodes\":{},\"triangles\":{},\"sections\":{},\"memory_branches\":{},\"compiled_controls\":{},\"pressure_point\":\"{observation}\",\"requested_tube_length_m\":{:e},\"represented_tube_length_m\":{:e},\"effective_mass_kg\":{:e},\"stiffness_n_m\":{:e},\"pressure_area_m2\":{:e},\"model_scope\":\"one linear plate mode; spatial lay and supplied material history; {propagation}; {scope}\"{receiver_json}{gas_json}}}",
            i.input_hash.to_hex(),i.nodes,i.triangles,i.sections,i.memory_branches,i.compiled_controls,
            requested,i.represented_tube_length_m,plate.mass_kg(),plate.stiffness_n_m(),plate.pressure_area_m2());
        Ok(Loaded{source:Box::new(p.into_renderer()),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"plate-valve-input",source_json,observation_json:None})
    } else {
        let p=ReedPerformance::from_bytes_with_gas_sources(bytes,block,policy,sources).map_err(|e|e.to_string())?;
        let i=p.info();
        let observer_json=format!("{}{}",super::duct_metadata::reed_provenance(&p),super::thermochemical::provenance(p.gas_records()));
        let policy_name=if policy==ObservationPolicy::Acceptance {"acceptance"} else {"diagnostic"};
        let source_json=format!("\"reed_input\":{{\"schema\":\"{REED_PERFORMANCE_SCHEMA}\",\"observation_policy\":\"{policy_name}\",\"blake3\":\"{}\",\"segments\":{},\"tone_holes\":{},\"massive_reed\":{},\"gesture_events\":{},\"compiled_controls\":{}{observer_json}}}",
            i.input_hash.to_hex(),i.segments,i.tone_holes,i.massive_reed,i.gesture_events,i.compiled_controls);
        Ok(Loaded{source:Box::new(p),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"reed-input",source_json,observation_json:None})
    }
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    run_with_gas_sources(args, None)
}
pub(super) fn run_with_gas_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    let (input, output, block, decimate, policy) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("wind input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_PLATE_VALVE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("wind input read failed: {e}"))?;
    let parent = Path::new(input).parent().unwrap_or(Path::new("."));
    let loaded = load_with_sources(&bytes, block, policy, sources, |name| {
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
    })?;
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
    let observation_json = loaded.observation_json.unwrap_or_else(|| format!(
        "{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"output_sample_rate_hz\":{RATE},\"ratio\":{},\"filter\":\"{}\",\"delay_output_samples\":{},\"first_output_source_index\":{},\"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}",
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
