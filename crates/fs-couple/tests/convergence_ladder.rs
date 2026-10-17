//! G1/G3 measured ladders on the existing M-CL0 ReedBoreVoice.
//! This fixture verifies the harness. It does not validate this legacy reed
//! model or claim that finite-horizon onset equals an asymptotic bifurcation.
use fs_couple::{reed_bore::ReedSolverMode, render::ReedBoreVoice, thin_plate::PlateBank};
use fs_duct::{Duct, Segment, Termination};
use fs_ladder::axis::*;
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_scenario::BeatingReed;
use std::{
    collections::BTreeMap,
    io::Write,
    path::Path,
    sync::{Mutex, OnceLock},
    time::Instant,
};

const LENGTH: f64 = 0.5;
const SECONDS: f64 = 0.25;
const ONSET_RMS_PA: f64 = 1e-3;
const ONSET_STEPS: usize = 8;

#[derive(Clone)]
struct Measurement {
    record: RungRecord,
    signal: Vec<f64>,
    rate: u32,
    fft_request: usize,
    mode: ReedSolverMode,
    tolerance: f64,
}
type Key = (u32, usize, u64, bool);
static CACHE: OnceLock<Mutex<BTreeMap<Key, Measurement>>> = OnceLock::new();

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn duct() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: LENGTH,
        }],
    }
}
fn voice(rate: u32, fft: usize, tol: f64, mode: ReedSolverMode, pressure: f64) -> ReedBoreVoice {
    let mut v = ReedBoreVoice::new_with_fft_length(
        &duct(),
        &gas(),
        BeatingReed {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: pressure,
            attack_s: 0.008,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        },
        // The 96 and 192 kHz rungs put Nyquist at ka 1.9 and 3.9 for this bore.
        // The unflanged low-ka fit stops at ka 1 and the exact thin-wall load
        // at the 3.83 cutoff; the multimodal flanged load holds to ka 20.
        Termination::FlangedOpen,
        PlateBank::default(),
        1.0,
        rate,
        rate as usize / 2,
        None,
        fft,
    )
    .unwrap();
    v.set_junction_step_tolerance(tol).unwrap();
    v.set_solver_mode(mode);
    v
}
fn rms(x: &[f64]) -> f64 {
    det::sqrt(x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64)
}
fn render(v: &mut ReedBoreVoice, rate: u32) -> Vec<f64> {
    let mut signal = vec![0.0; (SECONDS * f64::from(rate)) as usize];
    for block in signal.chunks_mut(256) {
        v.step_block(block).unwrap();
    }
    assert!(signal.iter().all(|x| x.is_finite()));
    signal
}
fn sustained(rate: u32, fft: usize, tol: f64, mode: ReedSolverMode, pressure: f64) -> bool {
    let signal = render(&mut voice(rate, fft, tol, mode, pressure), rate);
    let n = signal.len();
    let early = rms(&signal[n / 2..3 * n / 4]);
    let late = rms(&signal[3 * n / 4..]);
    late > ONSET_RMS_PA && late >= 0.999 * early
}
fn onset(rate: u32, fft: usize, tol: f64, mode: ReedSolverMode) -> f64 {
    let mut low = 0.0;
    let mut high = 4800.0;
    assert!(
        !sustained(rate, fft, tol, mode, low),
        "finite-horizon onset lower endpoint already sustains"
    );
    assert!(
        sustained(rate, fft, tol, mode, high),
        "finite-horizon onset is not bracketed; never substitute a threshold"
    );
    for _ in 0..ONSET_STEPS {
        let middle = 0.5 * (low + high);
        if sustained(rate, fft, tol, mode, middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}
fn pitch_hz(signal: &[f64], rate: u32) -> f64 {
    let x = &signal[signal.len() / 2..];
    let expected = gas().sound_speed / (4.0 * LENGTH);
    let min = (f64::from(rate) / (1.6 * expected)).floor() as usize;
    let max = (f64::from(rate) / (0.6 * expected)).ceil() as usize;
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    let corr = |lag: usize| {
        let mut cross = 0.0;
        let mut a = 0.0;
        let mut b = 0.0;
        for i in lag..x.len() {
            let u = x[i] - mean;
            let v = x[i - lag] - mean;
            cross += u * v;
            a += u * u;
            b += v * v;
        }
        cross / det::sqrt((a * b).max(f64::MIN_POSITIVE))
    };
    let peak = (min + 1..max)
        .max_by(|a, b| corr(*a).total_cmp(&corr(*b)))
        .unwrap();
    let (a, b, c) = (corr(peak - 1), corr(peak), corr(peak + 1));
    let curvature = a - 2.0 * b + c;
    let offset = if curvature.abs() > f64::EPSILON {
        0.5 * (a - c) / curvature
    } else {
        0.0
    };
    f64::from(rate) / (peak as f64 + offset.clamp(-0.5, 0.5))
}
fn obs(name: impl Into<String>, value: f64, unit: &str) -> Observable {
    Observable {
        name: name.into(),
        value,
        unit: unit.into(),
    }
}
fn measure(rate: u32, fft: usize, tol: f64, mode: ReedSolverMode) -> Measurement {
    let key = (rate, fft, tol.to_bits(), mode == ReedSolverMode::Strict);
    let cache = CACHE.get_or_init(|| Mutex::new(BTreeMap::new()));
    if let Some(m) = cache.lock().unwrap().get(&key) {
        return m.clone();
    }
    let start = Instant::now();
    let mut v = voice(rate, fft, tol, mode, 2800.0);
    let resolution = v.realized_resolution();
    let (tail, _total) = v.initial_ir_tail_energy().unwrap();
    let signal = render(&mut v, rate);
    let hz = pitch_hz(&signal, rate);
    let reference = gas().sound_speed / (4.0 * LENGTH);
    let mut qoi = vec![
        obs(
            "pitch_cents",
            1200.0 * det::ln(hz / reference) / det::ln(2.0),
            "cent",
        ),
        obs("threshold_pa", onset(rate, fft, tol, mode), "Pa"),
    ];
    qoi.push(obs("threshold_percent", 100.0 * qoi[1].value / 6000.0, "%"));
    let x = &signal[signal.len() / 2..];
    for partial in 1..=8 {
        let mut re = 0.0;
        let mut im = 0.0;
        for (k, p) in x.iter().enumerate() {
            let phase =
                std::f64::consts::TAU * f64::from(partial) * hz * k as f64 / f64::from(rate);
            re += p * det::cos(phase);
            im -= p * det::sin(phase);
        }
        let amplitude =
            (2.0 * det::sqrt(re * re + im * im) / x.len() as f64).max(f64::MIN_POSITIVE);
        qoi.push(obs(
            format!("partial_{partial}_db"),
            20.0 * det::ln(amplitude / 2e-5) / det::ln(10.0),
            "dB re 20 uPa",
        ));
    }
    let residual = v
        .junction_residual_max()
        .expect("massless samples measured");
    qoi.push(obs("retained_ir_tail_l2", det::sqrt(tail), "1"));
    qoi.push(obs("junction_flow_defect", residual, "m3/s"));
    let mut dimensions = vec![
        Dimension {
            name: "dt".into(),
            value: 1.0 / resolution.sample_rate_hz,
            unit: "s".into(),
        },
        Dimension {
            name: "fft".into(),
            value: resolution.n_fft as f64,
            unit: "samples".into(),
        },
    ];
    if let Some(value) = resolution.step_tolerance {
        dimensions.push(Dimension {
            name: "tolerance".into(),
            value,
            unit: "1".into(),
        });
    }
    let record = RungRecord {
        requested: 1.0 / resolution.sample_rate_hz,
        dimensions,
        qoi,
        bounds: vec![
            obs("retained_ir_tail_l2", det::sqrt(tail), "1"),
            obs("junction_flow_defect", residual, "m3/s"),
        ],
        regime: format!("register-{}", (hz / reference).round() as u32),
        solver_residual: residual,
        wall_s: start.elapsed().as_secs_f64(),
    };
    let m = Measurement {
        record,
        signal,
        rate,
        fft_request: fft,
        mode,
        tolerance: tol,
    };
    cache.lock().unwrap().insert(key, m.clone());
    m
}
fn policy() -> AnalysisPolicy {
    AnalysisPolicy {
        maximum_order: 16.0,
        order_tolerance: 0.5,
        relative_fit_tolerance: 0.02,
        asymptotic_tolerance: 0.1,
        minimum_refinement_ratio: 1.1,
    }
}
fn axis(name: &str) -> RefinableAxis {
    let (unit, model, coordinate, targets) = match name {
        // The legacy observer differentiates jet force with a backward difference
        // (smooth-signal first order). This is a declared hypothesis to test,
        // not a claim of a global order for a nonlinear oscillation/bifurcation.
        "dt" => (
            "s",
            ErrorModel::Algebraic { formal_order: 1.0 },
            Coordinate::Direct,
            vec![],
        ),
        "fft" => (
            "samples",
            ErrorModel::SpectralTail,
            Coordinate::Reciprocal,
            vec!["retained_ir_tail_l2".into()],
        ),
        "tolerance" => (
            "1",
            ErrorModel::ResidualBound,
            Coordinate::Direct,
            vec!["junction_flow_defect".into()],
        ),
        _ => panic!("undeclared axis"),
    };
    RefinableAxis {
        name: name.into(),
        unit: unit.into(),
        model,
        group: None,
        coordinate,
        bound_targets: targets,
    }
}
fn axis_report(name: &str, measurements: &[Measurement]) -> AxisReport {
    let a = axis(name);
    let mut rs: Vec<_> = measurements.iter().map(|m| m.record.clone()).collect();
    let finest_pitch = rs.last().unwrap().qoi[0].value;
    for (r, m) in rs.iter_mut().zip(measurements) {
        r.requested = match name {
            "dt" => 1.0 / f64::from(m.rate),
            "fft" => m.fft_request as f64,
            _ => m.tolerance,
        };
        r.qoi[0].value -= finest_pitch;
    }
    // A measured waveform difference demonstrates sensitivity independently of
    // whether a main-study scalar QoI happens to be exactly invariant.
    let waveform_rms = |m: &Measurement| obs("observer_rms", rms(&m.signal), "Pa");
    let mut s=SensitivityFixture {fixture:"M-CL0: 0.5m/2.2mm cylinder, dry air 293.15K/101325Pa, published music_render reed card; raw observer RMS".into(),
        realized:[rs[0].dimensions.iter().find(|d|d.name==name).unwrap().value,
            rs[1].dimensions.iter().find(|d|d.name==name).unwrap().value],
        qoi:[waveform_rms(&measurements[0]),waveform_rms(&measurements[1])],minimum_change:0.0};
    if name == "tolerance" {
        // Main-study tolerances can agree to roundoff. The independent liveness
        // fixture deliberately spans a coarse and tight solve on this same path.
        s.realized = [0.1, 1e-12];
        s.fixture
            .push_str("; independent Newton sensitivity tolerances 0.1/1e-12");
        for (i, t) in s.realized.into_iter().enumerate() {
            let signal = render(
                &mut voice(48000, 2048, t, ReedSolverMode::FastNewton, 2800.0),
                48000,
            );
            s.qoi[i] = obs("observer_rms", rms(&signal), "Pa");
        }
    }
    let report = evaluate_axis(&a, &rs, &s, &policy());
    emit_axis(&report);
    assert!(report.live_rungs >= 3, "{report:?}");
    assert!(report.complete(), "{report:?}");
    report
}
fn emit_axis(report: &AxisReport) {
    for (i, r) in report.rungs.iter().enumerate() {
        let realized = r
            .dimensions
            .iter()
            .find(|d| d.name == report.axis.name)
            .unwrap()
            .value;
        let partials = r.qoi[3..11]
            .iter()
            .map(|q| q.value.to_string())
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "IRPLADDER {{\"stage\":\"rung\",\"lane\":\"reed-ladder\",\"axis\":{:?},\"model\":{:?},\"requested\":{},\"realized\":{},\"live\":{},\"regime\":{:?},\"qoi\":{{\"pitch_cents\":{},\"threshold_pa\":{},\"threshold_percent\":{},\"partial_db\":[{}]}},\"solver_residual\":{},\"wall_s\":{}}}",
            report.axis.name,
            format!("{:?}", report.axis.model),
            r.requested,
            realized,
            !report.dead_rungs.contains(&i),
            r.regime,
            r.qoi[0].value,
            r.qoi[1].value,
            r.qoi[2].value,
            partials,
            r.solver_residual,
            r.wall_s
        );
    }
    for q in &report.qois {
        let (order, gci, bound) = match q.verdict {
            AxisVerdict::Asymptotic {
                observed_order,
                gci,
            } => (Some(observed_order), gci, None),
            AxisVerdict::PreAsymptotic { observed_order } => (Some(observed_order), None, None),
            AxisVerdict::Bounded { bound } => (None, None, Some(bound)),
            _ => (None, None, None),
        };
        let number = |v: Option<f64>| v.map_or_else(|| "null".into(), |x| x.to_string());
        let formal = if let ErrorModel::Algebraic { formal_order } = report.axis.model {
            Some(formal_order)
        } else {
            None
        };
        println!(
            "IRPLADDER {{\"stage\":\"axis\",\"axis\":{:?},\"qoi\":{:?},\"verdict\":{:?},\"observed_order\":{},\"formal_order\":{},\"gci\":{},\"bound\":{},\"spread\":{}}}",
            report.axis.name,
            q.name,
            format!("{:?}", q.verdict),
            number(order),
            number(formal),
            number(gci),
            number(bound),
            q.spread
        );
    }
}

#[test]
fn reed_dt_ladder() {
    let mut ms: Vec<_> = [48000, 96000, 192000]
        .into_iter()
        .map(|rate| measure(rate, 8192, 1e-12, ReedSolverMode::Strict))
        .collect();
    let report = axis_report("dt", &ms);
    if report
        .qois
        .iter()
        .any(|q| matches!(q.verdict, AxisVerdict::PreAsymptotic { .. }))
    {
        ms.push(measure(384000, 8192, 1e-12, ReedSolverMode::Strict));
        axis_report("dt", &ms);
    }
}
#[test]
fn reed_fft_ladder_with_liveness() {
    let ms: Vec<_> = [2048, 4096, 8192, 16384]
        .into_iter()
        .map(|n| measure(48000, n, 1e-12, ReedSolverMode::Strict))
        .collect();
    let report = axis_report("fft", &ms);
    assert_eq!(report.dead_rungs, vec![3]);
    assert_eq!(report.rungs[2].dimensions[1].value, 8192.0);
    assert_eq!(report.rungs[3].dimensions[1].value, 8192.0);
}
#[test]
fn reed_solver_tolerance_ladder() {
    let ms: Vec<_> = [1e-6, 1e-9, 1e-12]
        .into_iter()
        .map(|t| measure(48000, 8192, t, ReedSolverMode::FastNewton))
        .collect();
    let report = axis_report("tolerance", &ms);
    let strict = measure(48000, 8192, 1e-12, ReedSolverMode::Strict);
    for (m, r) in ms.iter().zip(&report.rungs) {
        assert_eq!(r.dimensions[2].value, m.tolerance);
        let discrepancy = m
            .signal
            .iter()
            .zip(&strict.signal)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(discrepancy.is_finite());
        println!(
            "IRPLADDER {{\"stage\":\"strict-reference\",\"tolerance\":{},\"max_pressure_difference_pa\":{}}}",
            m.tolerance, discrepancy
        );
    }
}
#[test]
fn reed_joint_ladder() {
    let start = Instant::now();
    let dt: Vec<_> = [48000, 96000, 192000]
        .into_iter()
        .map(|r| measure(r, 8192, 1e-12, ReedSolverMode::Strict))
        .collect();
    let fft: Vec<_> = [2048, 4096, 8192]
        .into_iter()
        .map(|n| measure(192000, n, 1e-12, ReedSolverMode::Strict))
        .collect();
    let tol: Vec<_> = [1e-6, 1e-9, 1e-12]
        .into_iter()
        .map(|t| measure(192000, 8192, t, ReedSolverMode::FastNewton))
        .collect();
    let mut axes = vec![
        axis_report("dt", &dt),
        axis_report("fft", &fft),
        axis_report("tolerance", &tol),
    ];
    axes[0].axis.group = Some("physical-ir-duration".into());
    axes[1].axis.group = Some("physical-ir-duration".into());
    let joint: Vec<_> = [
        (48000, 2048, 1e-6),
        (96000, 4096, 1e-9),
        (192000, 8192, 1e-12),
    ]
    .into_iter()
    .map(|(r, n, t)| measure(r, n, t, ReedSolverMode::FastNewton))
    .collect();
    let group = CoupledGroup {
        name: "physical-ir-duration".into(),
        members: vec![("dt".into(), 1.0), ("fft".into(), 1.0)],
        invariant: 2048.0 / 48000.0,
        relative_tolerance: 1e-12,
    };
    let report = LadderReport::new(
        policy(),
        axes,
        vec![group],
        joint.iter().map(|m| m.record.clone()).collect(),
        EnsembleDecl::Deterministic {
            reason: "not applicable: deterministic fixture, no stochastic source".into(),
        },
    )
    .unwrap();
    assert_eq!(report.joint.interactions.len(), 3 * 13);
    for i in &report.joint.interactions {
        assert!(i.interaction.is_finite());
        println!(
            "IRPLADDER {{\"stage\":\"interaction\",\"rung\":{},\"qoi\":{:?},\"joint_error\":{},\"single_axis_sum\":{},\"interaction\":{}}}",
            i.rung, i.joint_error.name, i.joint_error.value, i.single_axis_sum, i.interaction
        );
    }
    retain(&report, joint.last().unwrap());
    println!(
        "IRPLADDER {{\"stage\":\"terminal\",\"verdict\":\"measured-no-realism-claim\",\"report_blake3\":{:?},\"git_sha\":{:?},\"host\":{:?},\"command\":{:?},\"ensemble\":\"not applicable: deterministic fixture, no stochastic source\",\"elapsed_s\":{}}}",
        report.hash.to_hex(),
        std::env::var("FRANKENSIM_LADDER_GIT_SHA")
            .unwrap_or_else(|_| "unrecorded-direct-test".into()),
        std::env::var("FRANKENSIM_LADDER_HOST").unwrap_or_else(|_| "unrecorded-direct-test".into()),
        std::env::var("FRANKENSIM_LADDER_COMMAND")
            .unwrap_or_else(|_| "cargo test -p fs-couple --test convergence_ladder".into()),
        start.elapsed().as_secs_f64()
    );
}
fn write_new(path: &Path, bytes: &[u8]) {
    let mut f = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .unwrap();
    f.write_all(bytes).unwrap();
    f.sync_all().unwrap();
}
fn retain(report: &LadderReport, m: &Measurement) {
    let Ok(out) = std::env::var("FRANKENSIM_LADDER_OUT_DIR") else {
        return;
    };
    let dir = Path::new(&out);
    let bytes = report.canonical_bytes();
    write_new(&dir.join("reed-ladder.report"), &bytes);
    let (wav, clipped) = fs_couple::pcm_wav::encode_pcm16_wav(&m.signal, m.rate, 100000.0).unwrap();
    assert_eq!(clipped, 0);
    write_new(&dir.join("reed.wav"), &wav);
    let sidecar = format!(
        "{{\"report_blake3\":{:?},\"wav_blake3\":{:?},\"sample_rate_hz\":{},\"n_fft_request\":{},\"solver_mode\":{:?},\"step_tolerance\":{},\"full_scale_pa\":100000,\"onset_horizon_s\":{},\"onset_rms_pa\":{},\"onset_bisections\":{},\"source_fixture\":\"music_render::reed_context / M-CL0\",\"inputs\":{{\"duct_length_m\":0.5,\"duct_radius_m\":0.0022,\"gas_spec\":\"dry_air_ussa1976\",\"temperature_k\":293.15,\"ambient_pressure_pa\":101325,\"rest_opening_m\":0.0004,\"reed_width_m\":0.013,\"closing_pressure_pa\":6000,\"blowing_pressure_pa\":2800,\"attack_s\":0.008,\"reed_mass_kg\":0,\"stiffness_n_m\":0,\"damping_ratio\":0.35,\"termination\":\"FlangedOpen\",\"plates\":[],\"listener_m\":1}},\"no_claim\":\"retained-periodic-IR tail and observed flow defect do not bound acoustic QoIs; finite-horizon onset is not a bifurcation threshold; legacy M-CL0 is not realism validation\"}}\n",
        report.hash.to_hex(),
        fs_blake3::hash_bytes(&wav).to_hex(),
        m.rate,
        m.fft_request,
        format!("{:?}", m.mode),
        m.tolerance,
        SECONDS,
        ONSET_RMS_PA,
        ONSET_STEPS
    );
    write_new(&dir.join("reed.provenance.json"), sidecar.as_bytes());
}

#[test]
fn resolution_and_tolerance_refusals_are_explicit() {
    let mut v = voice(48000, 8192, 1e-9, ReedSolverMode::Strict, 2800.0);
    assert_eq!(v.realized_resolution().n_fft, 8192);
    assert_eq!(
        voice(48000, 8193, 1e-9, ReedSolverMode::Strict, 2800.0)
            .realized_resolution()
            .n_fft,
        8192
    );
    assert_eq!(v.realized_resolution().step_tolerance, None);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0 + f64::EPSILON] {
        assert!(v.set_junction_step_tolerance(invalid).is_err());
    }
    assert!(v.set_junction_step_tolerance(1.0).is_ok());
    v.step_block(&mut [0.0; 1]).unwrap();
    assert!(v.initial_ir_tail_energy().is_err());
}
