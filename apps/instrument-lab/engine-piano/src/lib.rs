//! `fs-live-piano`: the REDUCED REAL-TIME piano model for the browser.
//!
//! It is explicitly separate from the offline reference
//! (`crates/fs-couple/examples/grand_piano`). It reuses the reference's own
//! source modules for the instrument description and every cold computation
//! (see [`refmods`]) and steps the same equations with the same discretization
//! (see [`kernel`]), but only for the keys that are sounding.
//!
//! Inputs are SI: strikes are post-escapement hammer velocities in m/s, the
//! output is pressure in Pa at a Rayleigh half-space receiver, mapped to float
//! with a FIXED declared full scale. There is no sample, no assigned pitch,
//! no normalization, no AGC, no EQ and no reverb anywhere in this crate.
pub mod audio;
pub mod dense;
pub mod fastfelt;
pub mod kernel;
pub mod preparation;
pub mod refmods;
pub mod stretch;
// The reference engine's radiation module names the bank as `crate::linear`.
#[cfg(not(target_arch = "wasm32"))]
#[allow(unused_imports)]
use refmods::linear;

pub use audio::PianoAudio;
pub use kernel::{Config, Counters, Description, Ledger, LivePiano};

#[cfg(target_arch = "wasm32")]
mod wasm;

/// MIDI velocity 127 maps to this POST-ESCAPEMENT hammer speed [m/s]; the map
/// is linear and explicitly UNCALIBRATED, exactly as `grand_piano/midi.rs`.
pub const MIDI_MAXIMUM_VELOCITY_M_S: f64 = 4.5;

/// Uncalibrated linear MIDI-velocity -> hammer-speed map of the reference.
pub fn midi_velocity_to_m_s(velocity: u8) -> f64 {
    MIDI_MAXIMUM_VELOCITY_M_S * (f64::from(velocity.min(127)) / 127.0)
}

/// Machine-readable statement of what the live model is (for the site).
pub fn describe(piano: &LivePiano) -> String {
    let (travel, filter) = piano.microphone_delay_samples();
    describe_clock(piano, piano.config().sample_rate, 1,
        (travel as f64, filter), piano.microphone_multiplies())
}

/// Host-clock description of the same physical engine, including observation
/// latency and microphone work per emitted frame, not per mechanical substep.
pub fn describe_audio(audio: &PianoAudio) -> String {
    describe_clock(audio.piano(), audio.output_rate(), audio.substeps(),
        audio.microphone_delay_frames(), audio.multiplies_per_output_frame())
}

fn describe_clock(piano: &LivePiano, output_rate: u32, substeps: usize,
    (travel, filter): (f64, f64), microphone_work: usize) -> String {
    let c = piano.config();
    let gas = piano.acoustic_gas();
    let rate = f64::from(output_rate);
    let freqs = piano.board_frequencies_hz();
    let provenance: String = piano.board_provenance().chars()
        .map(|ch| if ch == '"' || ch == '\\' { '\'' } else if ch.is_control() { ' ' } else { ch }).collect();
    format!(
        concat!(
            "{{\"model\":\"reduced real-time piano (not the offline reference)\",",
            "\"sampleRate\":{},\"mechanicalRateHz\":{},\"substeps\":{},",
            "\"keys\":{},\"speakingStrings\":{},\"tablePartials\":{},",
            "\"boardCoordinates\":{},\"boardExactLowModes\":{},\"boardCoordinateHz\":[{:.3},{:.3}],",
            "\"boardProvenance\":\"{}\",",
            "\"acousticGas\":{{\"temperatureK\":{},\"pressurePa\":{},\"densityKgM3\":{},\"soundSpeedMS\":{},\"gamma\":{}}},",
            "\"partialBandHz\":{},\"maxPartialsPerString\":{},\"courseModeBudget\":{},",
            "\"maxSoundingKeys\":{},\"maxDynamicModes\":{},\"damperPoints\":{},\"duplex\":{},\"staticAttachment\":{},",
            "\"microphoneM\":[{},{},{}],\"microphoneTravelSamples\":{},\"microphoneFilterDelaySamples\":{},",
            "\"microphoneLatencyMs\":{:.3},\"microphoneMultipliesPerSample\":{},",
            "\"fullScalePa\":{},\"midiMaxVelocityMS\":{},\"stations\":{},",
            "\"silentBankHz\":{},\"silentBankModes\":{},",
            "\"retireRingingJ\":{:e},\"retireDampedJ\":{:e},\"stringDamping\":\"{}\",\"hammerLaw\":\"{}\",",
            "\"stringGasLoading\":\"{}\",\"stringMaterials\":{}}}"
        ),
        output_rate, c.sample_rate, substeps, piano.courses().len(), piano.total_strings(), piano.table_modes(),
        piano.board_modes(), piano.board_exact_low_modes(), freqs.first().copied().unwrap_or(0.0), freqs.last().copied().unwrap_or(0.0),
        provenance,
        gas.temperature, gas.pressure, gas.density, gas.sound_speed, gas.gamma,
        c.partial_band_hz, c.max_partials, c.course_mode_budget,
        c.max_courses, c.max_dynamic_modes, c.damper_points, c.duplex, c.static_attachment,
        c.microphone_m[0], c.microphone_m[1], c.microphone_m[2], travel, filter,
        (travel + filter) / rate * 1000.0, microphone_work,
        c.full_scale_pa, MIDI_MAXIMUM_VELOCITY_M_S, LivePiano::stations(),
        c.silent_bank_hz, piano.bank_modes(),
        c.retire_ringing_j, c.retire_damped_j,
        if let Some(materials) = piano.string_materials() {
            match materials.bending_loss_law() {
                crate::preparation::materials::BendingLossLaw::KelvinVoigt => "source-kelvin-voigt",
                crate::preparation::materials::BendingLossLaw::LossAngle => "source-loss-angle",
            }
        }
            else if c.rt0425_string_damping { "rt0425" } else { "authored" },
        if c.rt0425_hammer_dissipation { "rt0425 K_H and R_H" }
            else if c.rt0425_hammer_stiffness { "rt0425 K_H, estimated Prony" } else { "legacy" },
        if piano.string_materials().is_some() { "sampled-cylinder-retained-modes" }
            else { "unresolved-effective-scale" },
        piano.string_materials().map_or_else(|| "null".into(), |materials| materials.json()),
    )
}
