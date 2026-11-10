//! Glottal-source plus tract-area-function command; reuses the streaming output owner.

use super::{RATE, create_outputs, json_string, stream_output};
use fs_blake3::DomainHasher;
use fs_couple::voice_performance::{
    MAX_VOICE_PERFORMANCE_BYTES, NasalSinus, VOICE_PERFORMANCE_SCHEMA, VoiceDescription,
};
use fs_material::gas_species::input::ThermochemicalSources;
use std::io::{Read, Write};
use std::path::Path;

const INPUT_HASH_DOMAIN: &str = "org.frankensim.fs-couple.voice-performance.v1";

/// A nasal card's `model_scope` without the `head_sphere` row, in place of the series sentence; with it, see
/// [`nasal_head_scope`].
const NASAL_SCOPE: &str = "lumped fold pHS on a branched plane-wave tract (the oral tract; from the velopharyngeal port \
    the nasal passage and the two passages to the nostrils; each declared maxillary sinus a lumped cavity through its \
    ostium) with TMM viscothermal loss and the declared radiation load at the lips and at each nostril; the mouth's and \
    both nostrils' flows radiated through infinite-baffle circular openings on one plane (no nostril in front of the \
    mouth's plane, no loading of one opening by another), no head diffraction; the supraglottal pressure is never the \
    output";

/// A nasal card's `model_scope` with the `head_sphere` row (N3b): each opening is a cap on its own ray from the head's
/// centre, the three pressures summed after their FIRs; `head_load` and `torso_sphere` as in the series sentences.
fn nasal_head_scope(head_load: bool, torso: bool) -> String {
    let body = if torso { "a rigid sphere with a rigid torso sphere under it" } else { "a rigid sphere" };
    let fir = if torso {
        "the two spheres' causal FIR, the torso's change below 1.25 times the resolved band"
    } else {
        "the sphere's causal FIR"
    };
    let load = match (head_load, torso) {
        (true, true) => "the radiation load at the lips and at each nostril is the declared flanged load plus that \
            opening's cap change on the same head and torso",
        (true, false) => "the radiation load at the lips and at each nostril is the declared flanged load plus that \
            opening's cap change on the same sphere",
        (false, _) => "the radiation load at the lips and at each nostril is still the declared one",
    };
    let absent = if torso { "no neck" } else { "no torso" };
    format!(
        "lumped fold pHS on a branched plane-wave tract (the oral tract; from the velopharyngeal port the nasal passage \
         and the two passages to the nostrils; each declared maxillary sinus a lumped cavity through its ostium) with \
         TMM viscothermal loss; the mouth's and each nostril's flow observed as a cap of that opening's area on {body}, \
         each cap centred on the ray from the sphere's centre through the opening's baffle centre (a modelling choice), \
         each opening's baffled pressure through {fir}, the baffle's image kept above the resolved band; the three \
         pressures summed after the FIRs (no loading of one opening by another, no scattering off another opening's \
         cap); {load}; {absent}; the supraglottal pressure is never the output"
    )
}

fn options(args: &[String]) -> Result<(&str, &str, usize), String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                block = Some(iter.next().and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?);
            }
            value if value.starts_with('-') => return Err(format!(
                "unsupported voice option {value:?}; the fold, tract, gas, duration, pressure and microphone belong to the file"
            )),
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render voice INPUT.performance OUT.wav [--block N]".into());
    };
    Ok((*input, *output, block.unwrap_or(512)))
}

pub(super) fn run_with_gas_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    let (input, output, block) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("voice input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_VOICE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("voice input read failed: {e}"))?;
    let description = VoiceDescription::from_bytes_with_gas_sources(&bytes, sources)?;
    if description.output_rate_hz != RATE {
        return Err("music_render requires an explicitly declared 48000 Hz output clock; no implicit resampling".into());
    }
    let samples = usize::try_from(description.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let mut performance = description.build()?;
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut performance, &mut audio, samples, block, description.full_scale_pa,
    )?;
    let mut hasher = DomainHasher::new(INPUT_HASH_DOMAIN);
    hasher.update(&bytes);
    let input_hash = hasher.finalize().to_hex();
    let wav_hash = rendered.hash.to_hex();
    let info = performance.info();
    let stats = performance.stats();
    let mic = description.microphone.position_m;
    let gas_sources =
        super::thermochemical::provenance(core::slice::from_ref(description.gas_record()));
    let tract_gas = tract_gas_json(&description);
    let tract_motion = tract_motion_json(&description, info);
    let rest_gap = rest_gap_json(&description.rest_gap);
    let observer = performance.observer_gas_report().map_or_else(String::new, |report| {
        format!(",\"observer_gas\":{}",
            super::duct_metadata::observer_gas_provenance(Some(report)))
    });
    // Rigid-end limit of this tract, mouth load and wall; see VoiceInfo::passive_formants_hz.
    let passive = match &info.passive_formants_hz {
        Ok(formants) => format!(
            "[{}]",
            formants.iter().map(|f| format!("{f:.1}")).collect::<Vec<_>>().join(",")
        ),
        Err(why) => format!("{{\"refused\":{}}}", json_string(why)),
    };
    // Only a performance with the `glottal_orifice` row (a flow rig) writes this; without it the sidecar is byte for byte
    // what it was. The mean runs over the render's last 2 s, the rig's read window.
    let rig = stats.rig_mean_glottal_flow_m3_s.map_or_else(String::new, |u| {
        format!("\"flow_rig\":{{\"mean_glottal_flow_last_2_s_m3_s\":{u:e}}},")
    });
    // Only a performance with the `jet_noise` row writes this; without it the sidecar is byte for byte what it was.
    let jet_band = performance.jet_noise_coverage().map_or_else(String::new, |c| {
        format!(
            "\"jet_noise_band\":{{\"open_gate_steps\":{},\"whole\":{},\"truncated\":{},\"silent\":{}}},",
            c.open_gate(), c.whole, c.truncated, c.silent
        )
    });
    // Only a performance with the `subglottal_airway` row writes this block; without it the sidecar is byte for byte what
    // it was. The wall is the card's `tract_wall` on every tube's perimeter, or rigid without one.
    let airway = match (&info.subglottal_airway, description.subglottal_airway) {
        (Some(receipt), Some(generations)) => format!(
            "\"subglottal_airway\":{{\"tree\":\"{}\",\"geometry\":\"{}\",\"generations\":{generations},\
             \"wall\":\"{}\",\"port_impedance_pa_s_m3\":{:e},\"dc_resistance_pa_s_m3\":{:e},\
             \"same_sample_reflection\":{:e},\"series_impedance_pa_s_m3\":{:e},\"reflectance_taps\":{},\
             \"source_taps\":{},\"schedule\":\"lung pressure at the tree's far end\"}},",
            fs_couple::voice_performance::SUBGLOTTAL_AIRWAY_TREE,
            fs_couple::voice_performance::SUBGLOTTAL_AIRWAY_GEOMETRY,
            if description.wall.is_some() { "the card's tract_wall pin on every tube's perimeter" } else { "rigid" },
            receipt.port_impedance, receipt.dc_resistance, receipt.same_sample_reflection,
            receipt.series_impedance, receipt.reflectance.n_fft, receipt.source_taps,
        ),
        _ => String::new(),
    };
    // Only a performance with the `head_sphere` row writes this block and the sphere's scope; without it the sidecar is byte
    // for byte what it was. `filter_delay_samples` is the microphone's whole filter latency, the head's FIR included.
    // Only one with `head_load` as well names the lips' load inside the block and in the scope. Only one with
    // `torso_sphere` as well names the torso inside the block, in the lips' load and in the scope.
    let snow = description.torso_sphere.as_ref();
    let torso = snow.map_or_else(String::new, |t| {
        format!(
            ",\"torso_sphere\":{{\"radius_m\":{:e},\"neck_gap_m\":{:e},\"source\":{}}}",
            t.torso.radius_m, t.torso.neck_gap_m, json_string(&t.source)
        )
    });
    // A nasal card (N3b) names its three caps inside the block as well; the lips' load form is each opening's.
    let openings = match (&description.nasal, description.head_load) {
        (None, _) => "",
        (Some(_), true) => ",\"openings\":\"the lips and both nostrils, each a cap of its own area on the ray from the \
            sphere's centre through its baffle centre (a modelling choice), each loaded in the lips_load form on its own \
            cap; the three pressures summed after the FIRs\"",
        (Some(_), false) => ",\"openings\":\"the lips and both nostrils, each a cap of its own area on the ray from the \
            sphere's centre through its baffle centre (a modelling choice); the three pressures summed after the FIRs\"",
    };
    let (head, scope) = match &description.head_sphere {
        Some(card) if description.head_load => (
            format!(
                "\"head_sphere\":{{\"radius_m\":{:e},\"source\":{},\"filter_delay_samples\":{:e},\
                 \"lips_load\":\"{}\"{openings}{torso}}},",
                card.radius_m, json_string(&card.source), performance.filter_delay_samples(),
                if snow.is_some() {
                    "flanged + (cap on this sphere and torso - baffled piston), the torso's part below 1.25 times the resolved band"
                } else {
                    "flanged + (cap on this sphere - baffled piston)"
                }
            ),
            if snow.is_some() {
                "lumped fold pHS on a plane-wave series tract with TMM viscothermal loss; mouth flow observed as a cap of the opening's area on a rigid sphere with a rigid torso sphere under it (the baffled opening's pressure through the two spheres' causal FIR, the torso's change below 1.25 times the resolved band, the baffle's image kept above the resolved band), and the lips' radiation load is the declared flanged load plus that cap's change on the same head and torso; no neck, nasal branch or side cavity; the supraglottal pressure is never the output"
            } else {
                "lumped fold pHS on a plane-wave series tract with TMM viscothermal loss; mouth flow observed as a cap of the opening's area on a rigid sphere (the baffled opening's pressure through the sphere's causal FIR, the baffle's image kept above the resolved band), and the lips' radiation load is the declared flanged load plus that cap's change on the same sphere; no torso, nasal branch or side cavity; the supraglottal pressure is never the output"
            },
        ),
        Some(card) => (
            format!(
                "\"head_sphere\":{{\"radius_m\":{:e},\"source\":{},\"filter_delay_samples\":{:e}{openings}{torso}}},",
                card.radius_m, json_string(&card.source), performance.filter_delay_samples()
            ),
            if snow.is_some() {
                "lumped fold pHS on a plane-wave series tract with TMM viscothermal loss and a declared radiation load; mouth flow observed as a cap of the opening's area on a rigid sphere with a rigid torso sphere under it (the baffled opening's pressure through the two spheres' causal FIR, the torso's change below 1.25 times the resolved band, the baffle's image kept above the resolved band; the lips' radiation load is still the declared one), no neck, nasal branch or side cavity; the supraglottal pressure is never the output"
            } else {
                "lumped fold pHS on a plane-wave series tract with TMM viscothermal loss and a declared radiation load; mouth flow observed as a cap of the opening's area on a rigid sphere (the baffled opening's pressure through the sphere's causal FIR, the baffle's image kept above the resolved band; the lips' radiation load is still the declared one), no torso, nasal branch or side cavity; the supraglottal pressure is never the output"
            },
        ),
        None => (
            String::new(),
            "lumped fold pHS on a plane-wave series tract with TMM viscothermal loss and a declared radiation load; mouth flow radiated through an infinite-baffle circular opening, no head diffraction, nasal branch or side cavity; the supraglottal pressure is never the output",
        ),
    };
    // Only a performance with the nasal rows writes this block and the nasal scope; without them the sidecar is byte for
    // byte what it was, and a nasal card without `head_sphere` keeps N3's scope byte for byte.
    let nasal = nasal_json(&description, info);
    let scope: std::borrow::Cow<'_, str> = match (&description.nasal, &description.head_sphere) {
        (None, _) => scope.into(),
        (Some(_), None) => NASAL_SCOPE.into(),
        (Some(_), Some(_)) => nasal_head_scope(description.head_load, snow.is_some()).into(),
    };
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"voice-input\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
         \"voice_input\":{{\"schema\":\"{VOICE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\
         \"mechanics_rate_hz\":{},\"two_mass\":{},\"sections\":{},\"tract_length_m\":{:e},\
         \"glottal_radius_m\":{:e},\"mouth_radius_m\":{:e},\"sound_speed_m_s\":{:e},\
         \"yielding_wall\":{},\"tract_source\":{},\"tract_licence\":{},\
         \"microphone_m\":[{:e},{:e},{:e}],\"mouth_flow_taps\":{},\"mouth_flow_rungs\":{},\
         \"passive_tract_formants_hz\":{passive},\
         \"peak_glottal_flow_m3_s\":{:e},\"min_glottal_flow_m3_s\":{:e},\"closed_fraction\":{:e},\
         \"max_force_passes\":{},\"mean_force_passes\":{:e},{rig}{jet_band}{head}{airway}{nasal}\
         \"model_scope\":\"{scope}\"{observer}}}{gas_sources}{tract_gas}{tract_motion}{rest_gap}}}",
        description.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        info.mechanics_rate_hz, info.two_mass, info.sections, info.tract_length_m,
        info.glottal_radius_m, info.mouth_radius_m, info.sound_speed_m_s,
        info.has_wall, json_string(&description.tract_source), json_string(&description.tract_licence),
        mic[0], mic[1], mic[2], info.mouth_flow_taps, info.mouth_flow_rungs,
        stats.peak_glottal_flow_m3_s, stats.min_glottal_flow_m3_s, stats.closed_fraction,
        stats.max_force_passes, stats.mean_force_passes,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"voice-input\",\
         \"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"peak_glottal_flow_m3_s\":{:e},\"min_glottal_flow_m3_s\":{:e},\"closed_fraction\":{:e},\
         \"max_force_passes\":{},\"mean_force_passes\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa,
        rendered.rms_pa, stats.peak_glottal_flow_m3_s, stats.min_glottal_flow_m3_s, stats.closed_fraction,
        stats.max_force_passes, stats.mean_force_passes,
    );
    Ok(())
}

/// The `tract_gas` row's block, `,"tract_gas":{...}`, or nothing without the row, so a sidecar without it is byte for
/// byte what it was. The voice's `sound_speed_m_s` is the tract's; the record outside this block stays the exterior's.
pub(super) fn tract_gas_json(description: &VoiceDescription) -> String {
    description.tract_gas_record().map_or_else(String::new, |record| {
        let (tract, exterior) = (record.state(), description.ambient_gas());
        format!(
            ",\"tract_gas\":{{\"sound_speed_m_s\":{:e},\"density_kg_m3\":{:e},\"dynamic_viscosity_pa_s\":{:e},\
             \"exterior_sound_speed_m_s\":{:e},\"exterior_density_kg_m3\":{:e}{}}}",
            tract.sound_speed, tract.density, tract.dynamic_viscosity, exterior.sound_speed, exterior.density,
            super::thermochemical::provenance(core::slice::from_ref(record))
        )
    })
}

/// The nasal rows' block, `"nasal":{...},`, or nothing without them, so a sidecar without them is byte for byte what it
/// was. `after_section_0_based` counts from 0: the port opens at the junction after that many plus one oral sections.
/// The first passage section's own area is recorded; the line takes `port_area_m2` in its place.
pub(super) fn nasal_json(
    description: &VoiceDescription,
    info: &fs_couple::voice_performance::VoiceInfo,
) -> String {
    let Some(card) = description.nasal.as_ref() else {
        return String::new();
    };
    let sections = card
        .sections
        .iter()
        .map(|s| format!("[{:e},{:e}]", s.area_m2, s.length_m))
        .collect::<Vec<_>>()
        .join(",");
    let pairs = card
        .pairs
        .iter()
        .map(|p| format!("[{:e},{:e},{:e}]", p.left_area_m2, p.right_area_m2, p.length_m))
        .collect::<Vec<_>>()
        .join(",");
    let sinus = |s: Option<&NasalSinus>| {
        s.map_or_else(
            || "null".to_string(),
            |s| {
                format!(
                    "{{\"from_nostril_m\":{:e},\"ostium_area_m2\":{:e},\"ostium_length_m\":{:e},\"volume_m3\":{:e}}}",
                    s.from_nostril_m, s.ostium_area_m2, s.ostium_length_m, s.volume_m3
                )
            },
        )
    };
    let taps = info.nostril_flow_taps.map_or_else(|| "null".to_string(), |[left, right]| format!("[{left},{right}]"));
    format!(
        "\"nasal\":{{\"source\":{},\"licence\":{},\"velum\":{{\"after_section_0_based\":{},\"port_area_m2\":{:e},\
         \"attach_source\":{}}},\"sections_area_length_m\":[{sections}],\"pairs_left_right_length_m\":[{pairs}],\
         \"sinus_left\":{},\"sinus_right\":{},\"nostrils\":{{\"x_m\":{:e},\"y_m\":{:e},\"source\":{}}},\
         \"nostril_flow_taps\":{taps}}},",
        json_string(&card.source), json_string(&card.licence), card.after_section, card.port_area_m2,
        json_string(&card.attach_source), sinus(card.sinus_left.as_ref()), sinus(card.sinus_right.as_ref()),
        card.nostril_x_m, card.nostril_y_m, json_string(&card.nostril_source),
    )
}

/// The `rest_gap` rows' block, `,"rest_gap_time_s_gap_m":[[time,gap],...]`, as the card declares them, or nothing
/// without the rows, so a sidecar without them is byte for byte what it was.
fn rest_gap_json(rows: &[(f64, f64)]) -> String {
    if rows.is_empty() {
        return String::new();
    }
    let rows = rows
        .iter()
        .map(|(time, gap)| format!("[{time:e},{gap:e}]"))
        .collect::<Vec<_>>()
        .join(",");
    format!(",\"rest_gap_time_s_gap_m\":[{rows}]")
}

/// The moving tract's block, `,"tract_motion":{...}`, or nothing without `tract_move` rows, so a sidecar without them
/// is byte for byte what it was.
pub(super) fn tract_motion_json(
    description: &VoiceDescription,
    info: &fs_couple::voice_performance::VoiceInfo,
) -> String {
    if description.tract_moves.is_empty() {
        return String::new();
    }
    let name = |target: Option<usize>| target.map_or("card", |at| description.keyframes[at].name.as_str());
    let moves = description
        .tract_moves
        .iter()
        .map(|m| {
            format!(
                "{{\"start_s\":{:e},\"end_s\":{:e},\"target\":{},\"steps\":{}}}",
                m.start_s, m.end_s, json_string(name(m.target)), m.steps
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let keyframes = description
        .keyframes
        .iter()
        .map(|k| {
            format!(
                "{{\"name\":{},\"source\":{},\"licence\":{},\"sections\":{},\"first_area_m2\":{:e}}}",
                json_string(&k.name), json_string(&k.source), json_string(&k.licence), k.sections.len(),
                k.sections[0].area_m2
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let formants = info
        .keyframe_formants_hz
        .iter()
        .map(|(target, f)| {
            let value = match f {
                Ok(f) => format!("[{}]", f.iter().map(|x| format!("{x:.1}")).collect::<Vec<_>>().join(",")),
                Err(why) => format!("{{\"refused\":{}}}", json_string(why)),
            };
            format!("{{\"target\":{},\"passive_tract_formants_hz\":{value}}}", json_string(target))
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        ",\"tract_motion\":{{\"shapes\":{},\"moves\":[{moves}],\"keyframes\":[{keyframes}],\"targets\":[{formants}],\
         \"limits\":\"frozen linear shapes crossfaded with the outgoing history carried (no volume source from moving \
         walls); the glottal section held at the card's first area; the card's mouth radius at the microphone\"}}",
        info.tract_shapes
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_string()).collect() }
    #[test]
    fn voice_file_owns_physics_and_clock_while_block_size_remains_host_owned() {
        for flag in ["--seconds", "--schedule", "--full-scale-pa", "--frequency", "--vowel"] {
            assert!(options(&args(&["input", "output", flag, "1"])).is_err());
        }
        for value in ["0", "65537", "-1", "NaN"] {
            assert!(options(&args(&["input", "output", "--block", value])).is_err());
        }
        assert!(options(&args(&["input", "output", "--block", "37", "--block", "64"])).is_err());
        assert!(options(&args(&["input"])).is_err());
        let a = args(&["input", "output", "--block", "37"]);
        assert_eq!(options(&a).unwrap(), ("input", "output", 37));
    }
    #[test]
    fn rest_gap_rows_reach_the_sidecar_and_a_card_without_them_adds_nothing() {
        assert_eq!(rest_gap_json(&[]), "");
        assert_eq!(
            rest_gap_json(&[(0.0, 2.5e-3), (0.2375, 3.57e-4)]),
            ",\"rest_gap_time_s_gap_m\":[[0e0,2.5e-3],[2.375e-1,3.57e-4]]"
        );
    }
    #[test]
    fn a_nasal_card_on_a_head_names_each_cap_its_placement_and_the_sum() {
        for (load, torso) in [(false, false), (false, true), (true, false), (true, true)] {
            let s = nasal_head_scope(load, torso);
            assert!(s.contains("(a modelling choice)") && s.contains("summed after the FIRs"), "{s}");
            assert!(s.contains("no loading of one opening by another, no scattering off another opening's cap"), "{s}");
            assert!(!s.contains("infinite-baffle") && !s.contains(['"', '\\']), "{s}");
            assert_eq!(s.contains("torso sphere under it"), torso, "{s}");
            assert_eq!(s.contains("plus that opening's cap change"), load, "{s}");
        }
    }
}
