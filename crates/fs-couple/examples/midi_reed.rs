//! Authored four-note SMF -> Performance -> existing reed schedule seam.
//! Test fixture only; the fixed bore repeats one configuration, never converts
//! score pitch to frequency. The onset observable matches convergence_ladder.
use fs_blake3::{ContentHash, hash_bytes};
use fs_couple::description::per_note::{
    ThresholdBudget, ThresholdModel, record_resolution, resolve_threshold,
};
use fs_couple::{render::ReedBoreVoice, thin_plate::PlateBank};
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::*;
use fs_scenario::gesture::{GestureTarget, GestureValue};
use fs_scenario::performance::*;
use fs_scenario::{BeatingReed, EntityCatalog, EntityDeclaration, Environment, Scenario};
use std::io::Write;
use std::path::Path;

pub const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
pub fn limits() -> SmfLimits {
    SmfLimits {
        max_input_bytes: 4096,
        max_tracks: 2,
        max_events: 64,
        max_control_ticks: 2000,
    }
}
pub fn budget() -> ThresholdBudget {
    ThresholdBudget {
        pressure_tolerance_pa: 18.75,
        max_evaluations: 10,
    }
}
pub fn fixture() -> (AssemblyDescription, ImportPolicy, DynamicsTable) {
    let scenario = Scenario::new(
        "test-fixture:authored-midi-reed",
        1,
        Environment::earth_lab(),
    );
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("test-fixture:reed-player"))
        .unwrap();
    let target = ParameterPath {
        entity: root,
        parameter: ParameterKey::new("player_pressure", PA),
    };
    let desc = AssemblyDescription::new(
        &scenario,
        AssemblyInput {
            root,
            catalog,
            parts: vec![],
            interfaces: vec![],
            configurations: vec![ConfigurationState {
                name: "fixed-open-bore".into(),
                assignments: vec![],
                provenance: StateProvenance::Declared {
                    rationale: "owner-authored M-CL0 test fixture; not a measured specimen".into(),
                },
            }],
            observers: vec![],
            reference_state: Some(ReferenceState::Assembled),
        },
        DEFAULT_ASSEMBLY_BUDGET,
    )
    .unwrap();
    let source = hash_bytes(b"owner-authored midi-reed test fixture v1");
    let policy = ImportPolicy {
        chart: NoteChart {
            source,
            entries: vec![NoteChartEntry {
                channel: 0,
                note: 60,
                configuration: "fixed-open-bore".into(),
            }],
        },
        velocities: VelocityBins {
            source,
            bins: vec![VelocityBin {
                first: 1,
                last: 127,
                marking: "mf".into(),
            }],
        },
        bindings: vec![CommandBinding {
            path: target.clone(),
            track_id: "blow".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.0),
        }],
        controllers: vec![],
        control_rate_hz: 1000,
        initial_tempo_us: None,
        observation_prior: hash_bytes(
            b"prior:fixed on-axis free-field listener 1 m; no recording input",
        ),
    };
    let table = DynamicsTable {
        source,
        rows: vec![DynamicsRow {
            marking: "mf".into(),
            controls: vec![DynamicControl {
                target,
                domain: (QtyAny::new(0.0, PA), QtyAny::new(4800.0, PA)),
                choice: ControlChoice::PerNote {
                    objective: DeclaredObjective::MinimumThresholdPressure,
                },
            }],
        }],
    };
    (desc, policy, table)
}
pub fn authored_midi() -> Vec<u8> {
    // 480 PPQ, 500000 us/qn; four 96-tick notes, each followed by 24 ticks
    // of release. Duration 0.5 s. One fixed-bore state, authored velocities.
    let mut track = vec![0, 0xff, 0x51, 3, 7, 0xa1, 0x20];
    for (index, velocity) in [60, 72, 84, 96].into_iter().enumerate() {
        track.extend_from_slice(&[
            if index == 0 { 0 } else { 24 },
            0x90,
            60,
            velocity,
            96,
            60,
            0,
        ]);
    }
    track.extend_from_slice(&[24, 0xff, 0x2f, 0]);
    let mut b = b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk".to_vec();
    b.extend_from_slice(&(track.len() as u32).to_be_bytes());
    b.extend(track);
    b
}

/// Existing reed fixture, strict junction solver, no recordings. Finite-window
/// sustained RMS is a declared model observable, not an asymptotic threshold.
pub struct ReedOnset;
impl ThresholdModel for ReedOnset {
    fn identity(&self) -> ContentHash {
        hash_bytes(concat!("M-CL0:duct=cylinder(0.0022m,0.5m);gas=dry_air_ussa1976(293.15K,101325Pa);",
            "reed=(0.0004m,0.013m,6000Pa,attack=0.008s,mass=0,stiffness=0,damping=0.35);",
            "termination=UnflangedOpen;plates=[];listener=1m;rate=48000;line_samples=24000;",
            "fft=4096;solver=Strict(adjacent-float);configured_fast_tolerance=1e-9;horizon=0.25s;block=256;",
            "onset=RMS(last-quarter)>0.001Pa && RMS(last-quarter)>=0.999*RMS(previous-quarter);",
            "no asymptotic-bifurcation or realism claim").as_bytes())
    }
    fn oscillates(&mut self, pressure_pa: f64) -> Result<bool, String> {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0)
            .map_err(|e| format!("{e:?}"))?;
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0022,
                length: 0.5,
            }],
        };
        let reed = BeatingReed {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: pressure_pa,
            attack_s: 0.008,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        };
        let mut voice = ReedBoreVoice::new_with_fft_length(
            &duct,
            &gas,
            reed,
            Termination::UnflangedOpen,
            PlateBank::default(),
            1.0,
            48_000,
            24_000,
            None,
            4096,
        )
        .map_err(|e| format!("{e:?}"))?;
        voice
            .set_junction_step_tolerance(1e-9)
            .map_err(|e| format!("{e:?}"))?;
        let mut signal = vec![0.0; 12_000];
        for block in signal.chunks_mut(256) {
            voice.step_block(block).map_err(|e| format!("{e:?}"))?;
        }
        if signal.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite onset trajectory".into());
        }
        let rms = |x: &[f64]| det::sqrt(x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64);
        let early = rms(&signal[6000..9000]);
        let late = rms(&signal[9000..]);
        Ok(late > 1e-3 && late >= 0.999 * early)
    }
}
fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?
        .write_all(bytes)?;
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let dir = args
        .next()
        .ok_or("usage: midi_reed OUT_DIR [--refuse-unmapped]")?;
    let refuse = args.next();
    if refuse.as_deref().is_some_and(|s| s != "--refuse-unmapped") || args.next().is_some() {
        return Err("unknown MIDI fixture argument".into());
    }
    let dir = Path::new(&dir);
    let (d, mut p, t) = fixture();
    if refuse.is_some() {
        p.chart.entries.clear();
    }
    let bytes = authored_midi();
    let mut perf = import_smf(&bytes, &d, &p, &t, limits())?;
    let mut receipts = Vec::new();
    for i in 0..perf.per_note().len() {
        let receipt = resolve_threshold(&mut ReedOnset, &perf.per_note()[i], budget())?;
        record_resolution(&mut perf, i, &receipt)?;
        receipts.push(receipt);
    }
    let out = lower(&perf, &d, &p, &t)?;
    // The serialized artifact is independently decoded and lowered before write.
    let decoded = Performance::from_bytes(&perf.to_bytes(), &d, &p, &t, limits())?;
    if lower(&decoded, &d, &p, &t)?.schedule.to_canonical_bytes()
        != out.schedule.to_canonical_bytes()
    {
        return Err("performance replay differs".into());
    }
    std::fs::create_dir(dir)?;
    write_new(&dir.join("score.mid"), &bytes)?;
    write_new(&dir.join("score.performance"), &perf.to_bytes())?;
    write_new(
        &dir.join("score.gesture"),
        &out.schedule.to_canonical_bytes(),
    )?;
    write_new(
        &dir.join("score.assembly"),
        &d.to_bytes(DEFAULT_ASSEMBLY_BUDGET)?,
    )?;
    let mut per_note = Vec::new();
    for (c, r) in out.provenance.per_note.iter().zip(&receipts) {
        let resolved = c.resolved.ok_or("lowered control has no resolved value")?;
        let receipt = c.receipt.ok_or("lowered control has no receipt")?;
        let probes: Vec<_> = r
            .probes
            .iter()
            .map(|(pressure, on)| format!("[{pressure},{on}]"))
            .collect();
        per_note.push(format!(
            "{{\"note_index\":{},\"target\":\"{}:{}\",\"objective\":\"MinimumThresholdPressure\",\"domain\":[{},{}],\"units\":\"Pa\",\"resolved\":{},\"receipt\":\"{}\",\"model\":\"{}\",\"tolerance_pa\":{},\"max_evaluations\":{},\"bracket_pa\":[{},{}],\"probes_pa\":[{}]}}",
            c.note_index,c.target.entity.token(),c.target.parameter.name,c.domain.0.value,c.domain.1.value,resolved.value,receipt.to_hex(),r.model.to_hex(),r.tolerance_pa,r.max_evaluations,r.bracket_pa.0,r.bracket_pa.1,probes.join(",")
        ));
    }
    let metadata = format!(
        "{{\"midi_blake3\":\"{}\",\"chart_blake3\":\"{}\",\"dynamics_table_blake3\":\"{}\",\"schedule_blake3\":\"{}\",\"performance_blake3\":\"{}\",\"description_blake3\":\"{}\",\"observation_prior_blake3\":\"{}\",\"per_note\":[{}],\"source\":\"owner-authored test fixture\",\"no_claim\":\"score-to-render plumbing on legacy reed; no physical validation or realism\"}}\n",
        out.provenance.midi.to_hex(),
        out.provenance.chart.to_hex(),
        out.provenance.dynamics.to_hex(),
        out.schedule.content_hash().to_hex(),
        out.provenance.performance.to_hex(),
        out.provenance.description.to_hex(),
        out.provenance.observation_prior.to_hex(),
        per_note.join(",")
    );
    write_new(
        &dir.join("performance-provenance.json"),
        metadata.as_bytes(),
    )?;
    println!("{metadata}");
    Ok(())
}
