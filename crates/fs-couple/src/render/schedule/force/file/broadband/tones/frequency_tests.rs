use super::tests::{driver, observer};
use super::*;
use fs_math::det;

const RATE: u32 = 8000;
const NATURAL: f64 = core::f64::consts::TAU * 700.0;
const LIMITS: &str = "tone_limits 8 1048576 10000000 100000\n";
const FREQUENCIES: [(u64, f64); 5] = [
    (0, core::f64::consts::TAU * 600.0),
    (128, core::f64::consts::TAU * 800.0),
    (192, core::f64::consts::TAU * 800.0),
    (256, core::f64::consts::TAU * 500.0),
    (359, NATURAL),
];
const GAINS: [(u64, f64); 8] = [
    (0, 0.0),
    (41, 0.0),
    (97, 1.0),
    (128, 1.0),
    (191, 0.0),
    (223, 0.0),
    (307, -0.5),
    (359, -0.5),
];

fn source_text() -> String {
    let gain = GAINS
        .iter()
        .map(|(sample, gain)| format!("gain {sample} {gain}\n"))
        .collect::<String>();
    let frequency = FREQUENCIES
        .iter()
        .map(|(sample, omega)| format!("omega {sample} {omega:.17e}\n"))
        .collect::<String>();
    format!(
        "{LIMITS}tones 1\ntone 0 0 {NATURAL:.17e} 0.03 -0.01\nenvelopes 1\nenvelope 0 {}\n{gain}frequency_sweeps 1\nsweep 0 {} 0.0002\n{frequency}",
        GAINS.len(),
        FREQUENCIES.len()
    )
}

// Independent whole-segment integral of the authored instantaneous frequency.
fn phase(sample: f64) -> f64 {
    let mut result = 0.0;
    for pair in FREQUENCIES.windows(2) {
        let [(start, a), (end, b)] = [pair[0], pair[1]];
        if sample <= start as f64 {
            break;
        }
        let elapsed = sample.min(end as f64) - start as f64;
        let slope = (b - a) / (end - start) as f64;
        result += (a * elapsed + 0.5 * slope * elapsed * elapsed) / f64::from(RATE);
    }
    result
}

fn gain(sample: f64) -> f64 {
    for pair in GAINS.windows(2) {
        let [(start, a), (end, b)] = [pair[0], pair[1]];
        if sample <= end as f64 {
            return a + (b - a) * (sample - start as f64) / (end - start) as f64;
        }
    }
    GAINS.last().unwrap().1
}

fn derivative(sample: f64, state: [f64; 4]) -> [f64; 4] {
    let angle = phase(sample);
    let force = -2.0 * gain(sample) * (0.03 * det::cos(angle) - 0.01 * det::sin(angle));
    let [q, v, _, _] = state;
    let loss = 0.04 * NATURAL;
    [
        v,
        force - loss * v - NATURAL * NATURAL * q,
        force * v,
        loss * v * v,
    ]
}

#[test]
fn swept_source_follows_integrated_phase_through_silence_and_physical_ode() {
    let text = source_text();
    let mut observed = observer(0.02);
    let ports = [vec![-2.0]];
    let program = ToneProgram::read(&text, 0, &[&ports], RATE, 1024, 0.9, &mut 10000).unwrap();
    let bound = program.maximum_sweep_force_error_bound_n;
    assert!(bound > 0.0 && bound <= 0.0002);
    let mut source = program
        .sources
        .into_iter()
        .next()
        .unwrap()
        .unwrap()
        .prepare(&observed, program.budget)
        .unwrap();
    let mut reference = [0.0; 4];
    let (mut work, mut heat) = (0.0, 0.0);
    for sample in 0..359 {
        // RK4 integrates the original chirp independently of the affine source
        // reconstruction and its cached mechanics/observer operators.
        for sub in 0..32 {
            let at = sample as f64 + sub as f64 / 32.0;
            let h = 1.0 / (f64::from(RATE) * 32.0);
            let k1 = derivative(at, reference);
            let k2 = derivative(
                at + 1.0 / 64.0,
                core::array::from_fn(|i| reference[i] + h * k1[i] / 2.0),
            );
            let k3 = derivative(
                at + 1.0 / 64.0,
                core::array::from_fn(|i| reference[i] + h * k2[i] / 2.0),
            );
            let k4 = derivative(
                at + 1.0 / 32.0,
                core::array::from_fn(|i| reference[i] + h * k3[i]),
            );
            reference = core::array::from_fn(|i| {
                reference[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
        let pressure = source.step(&mut observed, &[0.0]).unwrap();
        let angle = phase((sample + 1) as f64);
        let expected_c = -0.06 * det::cos(angle) + 0.02 * det::sin(angle);
        let expected_s = 0.02 * det::cos(angle) + 0.06 * det::sin(angle);
        assert!((source.cosine[1] - expected_c).abs() < 2e-14);
        assert!((source.sine[1] - expected_s).abs() < 2e-14);
        // The velocity Green function has norm <= 1/sqrt(1-zeta^2).
        // Include the signed port magnitude two and pressure gain two.
        let pressure_bound = 4.01 * bound * (sample + 1) as f64 / f64::from(RATE);
        assert!((pressure - 2.0 * reference[1]).abs() < pressure_bound + 1e-13);
        if sample < 41 {
            assert_eq!(pressure, 0.0);
        }
        work += observed.last_modal_frame().input_work_j;
        heat += observed.last_modal_frame().viscous_dissipation_j;
        assert!(observed.last_modal_frame().viscous_dissipation_j >= 0.0);
    }
    assert!(heat > 0.0 && reference[3] > 0.0);
    assert!((work - heat - observed.last_modal_frame().total_modal_energy_j).abs() < 1e-15);
}

#[test]
fn rejected_sweep_step_preserves_phase_and_both_curve_clocks_for_retry() {
    let text = source_text();
    let (mut a, mut b) = (observer(0.02), observer(0.02));
    let (mut sa, mut sb) = (
        driver(&text, &[vec![1.0]], &a),
        driver(&text, &[vec![1.0]], &b),
    );
    for _ in 0..223 {
        sa.step(&mut a, &[0.0]).unwrap();
        sb.step(&mut b, &[0.0]).unwrap();
    }
    let phase = (sa.cosine.clone(), sa.sine.clone());
    let states = a.model().states().to_vec();
    let frame = a.last_modal_frame().clone();
    assert!(sa.step(&mut a, &[1e20]).is_err());
    assert_eq!(sa.envelope.as_ref().unwrap().sample, 223);
    assert_eq!((sa.cosine.clone(), sa.sine.clone()), phase);
    assert_eq!(a.model().states(), states);
    assert_eq!(a.last_modal_frame(), &frame);
    for _ in 223..359 {
        assert_eq!(
            sa.step(&mut a, &[0.0]).unwrap().to_bits(),
            sb.step(&mut b, &[0.0]).unwrap().to_bits()
        );
        assert_eq!(a.last_modal_frame(), b.last_modal_frame());
    }
}

#[test]
fn original_sweep_ids_guards_and_resources_survive_voice_grouping() {
    let text = format!(
        "{LIMITS}tones 2\ntone 1 0 {NATURAL} 0.03 -0.01\ntone 0 0 {NATURAL} 0.03 -0.01\nfrequency_sweeps 1\nsweep 0 2 0.0002\nomega 0 4000\nomega 32 4500\n"
    );
    let ports = [vec![1.0]];
    let program =
        ToneProgram::read(&text, 0, &[&ports, &ports], RATE, 32, 0.9, &mut 10000).unwrap();
    assert_eq!(program.frequency_sweep_count, 1);
    assert_eq!(program.envelope_count, 0);
    assert!(program.sources[0].as_ref().unwrap().envelopes.is_none());
    assert!(program.sources[0].as_ref().unwrap().sweeps.is_none());
    assert!(program.sources[1].as_ref().unwrap().envelopes.is_some());
    assert!(program.sources[1].as_ref().unwrap().sweeps.is_some());
    for bad in [
        text.replace("1048576", "1000"),
        text.replace("100000\n", "16\n"),
        text.replace("0.0002", "0.000000000001"),
        text.replace("omega 32 4500", "omega 32 23000"),
    ] {
        assert!(ToneProgram::read(&bad, 0, &[&ports, &ports], RATE, 32, 0.9, &mut 10000).is_err());
    }
    assert!(ToneProgram::read(&text, 0, &[&ports, &ports], RATE, 32, 0.1, &mut 10000).is_err());
}
