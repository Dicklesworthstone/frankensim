use super::tests::{driver, observer};
use super::*;
use fs_math::{c64::C64, det};

const RATE: u32 = 8000;
const NATURAL: f64 = core::f64::consts::TAU * 700.0;
const DRIVE: f64 = core::f64::consts::TAU * 2600.0;
const LIMITS: &str = "tone_limits 8 1048576 10000000 100000\n";
const KNOTS: [(u64, f64); 8] = [
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
    let knots = KNOTS
        .iter()
        .map(|(sample, gain)| format!("gain {sample} {gain}\n"))
        .collect::<String>();
    format!(
        "{LIMITS}tones 1\ntone 0 0 {DRIVE:.17e} 0.03 -0.01\nenvelopes 1\nenvelope 0 {}\n{knots}",
        KNOTS.len()
    )
}

// Independent polynomial-times-carrier particular solution of q''+w^2*q=F.
// Whole knot segments are evaluated directly, independently of numerical steps.
fn expected(sample: u64) -> (f64, f64) {
    let mut q = 0.0;
    let mut v = 0.0;
    for pair in KNOTS.windows(2) {
        let [(start, a), (stop, b)] = [pair[0], pair[1]];
        if start >= sample {
            break;
        }
        let end = sample.min(stop);
        let duration = (end - start) as f64 / f64::from(RATE);
        let end_gain = a + (b - a) * (end - start) as f64 / (stop - start) as f64;
        let phase = DRIVE * start as f64 / f64::from(RATE);
        let c = -0.06 * det::cos(phase) + 0.02 * det::sin(phase);
        let s = 0.02 * det::cos(phase) + 0.06 * det::sin(phase);
        let den = C64::from_re(NATURAL * NATURAL - DRIVE * DRIVE);
        let f0 = C64::new(a * c, -a * s);
        let slope = C64::new(
            (end_gain - a) * c / duration,
            -(end_gain - a) * s / duration,
        );
        let linear = slope / den;
        let constant = f0 / den - C64::new(0.0, 2.0 * DRIVE) * linear / den;
        let phase_end = C64::new(det::cos(DRIVE * duration), det::sin(DRIVE * duration));
        let end_coefficient = constant + linear * C64::from_re(duration);
        let initial_velocity = (linear + C64::new(0.0, DRIVE) * constant).re;
        let cw = det::cos(NATURAL * duration);
        let sw = det::sin(NATURAL * duration);
        let forced_q =
            (end_coefficient * phase_end).re - constant.re * cw - initial_velocity / NATURAL * sw;
        let forced_v = ((linear + C64::new(0.0, DRIVE) * end_coefficient) * phase_end).re
            + NATURAL * constant.re * sw
            - initial_velocity * cw;
        (q, v) = (
            q * cw + v / NATURAL * sw + forced_q,
            -NATURAL * q * sw + v * cw + forced_v,
        );
    }
    (q, v)
}

#[test]
fn attacks_releases_and_signed_gain_follow_continuous_force_work_without_phase_reset() {
    let mut observed = observer(0.0);
    let mut source = driver(&source_text(), &[vec![-2.0]], &observed);
    let mut work = 0.0;
    for sample in 1..=359 {
        let pressure = source.step(&mut observed, &[0.0]).unwrap();
        let (q, v) = expected(sample);
        assert!(
            (pressure - 2.0 * v).abs() < 3e-14,
            "sample {sample}: {pressure} vs {}",
            2.0 * v
        );
        assert!((observed.model().states()[0].displacement_m_sqrt_kg - q).abs() < 3e-18);
        if sample <= 41 {
            assert_eq!(pressure, 0.0);
        }
        // A held start gain would apply zero force over this entire first
        // attack interval. The affine owner integrates its actual rising load.
        if sample == 42 {
            assert!(pressure.abs() > 1e-12);
        }
        work += observed.last_modal_frame().input_work_j;
        assert_eq!(observed.last_modal_frame().viscous_dissipation_j, 0.0);
    }
    assert!((work - observed.last_modal_frame().total_modal_energy_j).abs() < 1e-17);
    assert_eq!(source.envelope.as_ref().unwrap().sample, 359);
}

#[test]
fn refused_affine_sample_retains_carrier_curve_cursor_and_physical_frame() {
    let text = source_text();
    let (mut a, mut b) = (observer(0.02), observer(0.02));
    let (mut sa, mut sb) = (
        driver(&text, &[vec![-2.0]], &a),
        driver(&text, &[vec![-2.0]], &b),
    );
    for _ in 0..96 {
        sa.step(&mut a, &[0.0]).unwrap();
        sb.step(&mut b, &[0.0]).unwrap();
    }
    let phase = (sa.cosine[1..].to_vec(), sa.sine[1..].to_vec());
    let states = a.model().states().to_vec();
    let frame = a.last_modal_frame().clone();
    assert!(sa.step(&mut a, &[1e20]).is_err());
    assert_eq!(sa.envelope.as_ref().unwrap().sample, 96);
    assert_eq!((sa.cosine[1..].to_vec(), sa.sine[1..].to_vec()), phase);
    assert_eq!(a.model().states(), states);
    assert_eq!(a.last_modal_frame(), &frame);
    for _ in 96..359 {
        assert_eq!(
            sa.step(&mut a, &[0.0]).unwrap().to_bits(),
            sb.step(&mut b, &[0.0]).unwrap().to_bits()
        );
        assert_eq!(a.last_modal_frame(), b.last_modal_frame());
    }
}

#[test]
fn original_tone_indices_survive_voice_grouping_and_endpoint_buffers_are_budgeted() {
    let text = format!(
        "{LIMITS}tones 2\ntone 1 0 {DRIVE} 0.03 -0.01\ntone 0 0 {DRIVE} 0.03 -0.01\nenvelopes 1\nenvelope 0 2\ngain 0 0\ngain 8 1\n"
    );
    let ports = [vec![-2.0]];
    let program = ToneProgram::read(&text, 0, &[&ports, &ports], RATE, 32, 0.9, &mut 1000).unwrap();
    assert_eq!(program.envelope_count, 1);
    assert!(program.sources[0].as_ref().unwrap().envelopes.is_none());
    assert!(program.sources[1].as_ref().unwrap().envelopes.is_some());
    for bad in [
        text.replace("1048576", "1000"),
        text.replace("100000\n", "8\n"),
    ] {
        assert!(ToneProgram::read(&bad, 0, &[&ports, &ports], RATE, 32, 0.9, &mut 1000).is_err());
    }
}
