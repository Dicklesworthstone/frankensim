use super::*;
use crate::graph::tests::with_cx;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}

fn pipe(cx: &Cx<'_>) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap(),
        LaminarPipeSpec {
            length: Length::new(0.4),
            radius: Length::new(0.001),
            maximum_reynolds: 1000.0,
            band_rad_s: [0.01, 900.0],
            maximum_relative_wall_error: 1e-4,
        },
        &ZkDiffusiveBasis::new(16).unwrap(),
    )
    .unwrap()
}

fn tone(cosine: f64, sine: f64) -> LaminarPressureTone {
    LaminarPressureTone {
        cosine: Pressure::new(cosine),
        sine: Pressure::new(sine),
    }
}

// Closure of one prepared exact step's energy balance, in units of the step's own scale. The fs-phs
// admission (linear_exact.rs, balance_allowance; bead rc-1109 .1.5, ruling R1) is max(1e-12, 2 * 2^s * eps)
// of the scale, where s is the number of doublings of the plain-f64 positive-heat chain that builds the
// dissipated side: that chain carries about 2^s eps of rounding (measured at most 1.18 * 2^s eps over 57
// chains with s = 14..19), and 2 is the leading factor of 2 (2^s - 1) gamma_n with gamma_n -> u.
// The pulsatile-start-up test closes the same quantity, so they carry the same allowance. s is READ from each fixture's run
// (the chain's doublings printed by a diagnostic build), not derived: a change to a fixture that adds
// doublings must move its s with it.
fn closure_tolerance(doublings: i32) -> f64 {
    2.0 * 2f64.powi(doublings) * f64::EPSILON
}

fn closes(actual: f64, expected: f64, scale: f64, tolerance: f64) -> bool {
    (actual - expected).abs() <= tolerance * scale
}

fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        closes(actual, expected, scale, tolerance),
        "{actual:e} != {expected:e} (scale {scale:e})"
    );
}

// Independent radial J0 eigenfunctions of the Stokes equation, with
// tau=rho*r^2/mu=1 s. Each mode solves q'_n + j_0n^2*q_n =
// (32/R/j_0n^4)*j_0n^2*dp(t). The exact Rayleigh sum gives the
// instantaneous part; the remaining truncated tail decays as j_0n^-6.
fn radial_flow(time: f64, omega: f64, cosine: f64, sine: f64, resistance: f64) -> f64 {
    const ZEROS: [f64; 16] = [
        2.404_825_557_695_772_4,
        5.520_078_110_286_311,
        8.653_727_912_911_013,
        11.791_534_439_014_281,
        14.930_917_708_487_787,
        18.071_063_967_910_924,
        21.211_636_629_879_258,
        24.352_471_530_749_302,
        27.493_479_132_040_253,
        30.634_606_468_431_976,
        33.775_820_213_573_57,
        36.917_098_353_664_045,
        40.058_425_764_628_24,
        43.199_791_713_176_73,
        46.341_188_371_661_815,
        49.482_609_897_397_815,
    ];
    let (sn, cs) = (omega * time).sin_cos();
    let pressure = cosine * cs + sine * sn;
    let correction: f64 = ZEROS
        .into_iter()
        .map(|zero| {
            let rate = zero * zero;
            let decay = (-rate * time).exp();
            let denominator = rate * rate + omega * omega;
            let c = rate * (rate * cs + omega * sn - rate * decay) / denominator;
            let s = rate * (rate * sn - omega * cs + omega * decay) / denominator;
            (pressure - cosine * c - sine * s) / rate.powi(2)
        })
        .sum();
    (pressure - 32.0 * correction) / resistance
}

#[test]
fn pulsatile_startup_matches_independent_radial_stokes_modes() {
    with_cx(|cx| {
        let pipe = pipe(cx);
        for (omega, c, s) in [(3.0, 20.0, 0.0), (17.0, 0.0, 20.0), (17.0, -12.0, 16.0)] {
            let mut state = pipe
                .steady_state(VolumetricFlowRate::new(0.0), 0.0)
                .unwrap();
            let mut source = [tone(c, s)];
            let report = pipe
                .prepare_tones(cx, &[omega], 0.25, budget())
                .unwrap()
                .advance(cx, &mut state, &mut source)
                .unwrap();
            let resistance = pipe.steady_resistance_pa_s_m3();
            let expected = radial_flow(0.25, omega, c, s, resistance);
            near(
                report.final_volume_flow_m3_s,
                expected,
                20.0 / resistance,
                2e-6,
            );
            assert!(report.native.dissipated > 0.0);
            // The balance residual of the 0.25 s tone step, in the gate's own scale (supplied or dissipated, the
            // larger). s = 19: the larger of the work chain (n = 23, 18 doublings) and the heat chain (n = 38,
            // 19 doublings) of this step, the value the tone gate takes. Measured 5.730560851424842e-12 of the
            // scale against the 2.3283064365386963e-10 allowed here.
            let scale = report.native.supplied.abs().max(report.native.dissipated);
            let tolerance = closure_tolerance(19);
            near(report.native.balance_residual(), 0.0, scale, tolerance);
            // Planted negative: a residual moved by three tolerances is refused.
            assert!(!closes(
                report.native.balance_residual() + 3.0 * tolerance * scale,
                0.0,
                scale,
                tolerance
            ));
            let (sn, cs) = (omega * 0.25).sin_cos();
            near(source[0].cosine.value(), c * cs + s * sn, 20.0, 1e-14);
            near(source[0].sine.value(), s * cs - c * sn, 20.0, 1e-14);
        }
    });
}

#[test]
fn continuous_pressure_subdivision_preserves_flow_displacement_work_and_heat() {
    with_cx(|cx| {
        let pipe = pipe(cx);
        let frequencies = [0.0, 3.0, 17.0];
        let initial = pipe
            .steady_state(VolumetricFlowRate::new(-1e-8), 2.0)
            .unwrap();
        let mut whole_state = initial.clone();
        let mut divided_state = initial;
        let mut whole_source = [tone(5.0, 0.0), tone(20.0, -4.0), tone(-12.0, 16.0)];
        let mut divided_source = whole_source;
        let whole = pipe
            .prepare_tones(cx, &frequencies, 0.25, budget())
            .unwrap()
            .advance(cx, &mut whole_state, &mut whole_source)
            .unwrap();
        let mut prepared = pipe
            .prepare_tones(cx, &frequencies, 0.05, budget())
            .unwrap();
        let mut volume = 0.0;
        let mut heat = 0.0;
        let mut work = 0.0;
        for _ in 0..5 {
            let step = prepared
                .advance(cx, &mut divided_state, &mut divided_source)
                .unwrap();
            volume += step.displaced_volume_m3;
            heat += step.native.dissipated;
            work += step.native.supplied;
        }
        for (&a, &b) in whole_state
            .coordinates()
            .iter()
            .zip(divided_state.coordinates())
        {
            near(a, b, whole_state.kinetic_energy_j().sqrt(), 2e-10);
        }
        near(volume, whole.displaced_volume_m3, 1e-8, 2e-10);
        near(heat, whole.native.dissipated, whole.native.dissipated, 2e-9);
        near(
            work,
            whole.native.supplied,
            whole.native.supplied.abs(),
            2e-9,
        );
        for (a, b) in whole_source.iter().zip(&divided_source) {
            near(a.cosine.value(), b.cosine.value(), 20.0, 1e-14);
            near(a.sine.value(), b.sine.value(), 20.0, 1e-14);
        }
    });
}

#[test]
fn dc_recovers_original_held_pressure_with_unchanged_source_phase() {
    with_cx(|cx| {
        let pipe = pipe(cx);
        let mut held = pipe
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut continuous = held.clone();
        let reference = pipe
            .prepare(cx, 0.1)
            .unwrap()
            .advance(cx, &mut held, [Pressure::new(20.0), Pressure::new(0.0)])
            .unwrap();
        let mut source = [tone(20.0, 0.0)];
        let result = pipe
            .prepare_tones(cx, &[0.0], 0.1, budget())
            .unwrap()
            .advance(cx, &mut continuous, &mut source)
            .unwrap();
        near(
            result.final_volume_flow_m3_s,
            reference.final_volume_flow_m3_s,
            1e-8,
            2e-10,
        );
        near(
            result.displaced_volume_m3,
            reference.displaced_volume_m3,
            1e-9,
            2e-10,
        );
        near(
            result.native.dissipated,
            reference.native.dissipated,
            reference.native.dissipated,
            2e-9,
        );
        assert_eq!(source, [tone(20.0, 0.0)]);
    });
}

#[test]
fn pressure_band_regime_and_history_refusals_preserve_flow_and_phase() {
    with_cx(|cx| {
        let pipe = pipe(cx);
        for frequencies in [vec![], vec![-1.0], vec![3.0, 3.0], vec![900.1], vec![0.001]] {
            assert!(pipe.prepare_tones(cx, &frequencies, 0.1, budget()).is_err());
        }
        let mut prepared = pipe
            .prepare_tones(cx, &[0.0, 17.0], 0.25, budget())
            .unwrap();
        let mut state = pipe
            .steady_state(VolumetricFlowRate::new(0.0), 2.0)
            .unwrap();
        let before = state.clone();
        // The sinusoidal pressure is zero at the interval's start, but its
        // interior amplitude cannot bypass the original Reynolds ceiling.
        for mut source in [
            [tone(0.0, 0.0), tone(0.0, 1e8)],
            [tone(1.0, 1.0), tone(2.0, 0.0)],
        ] {
            let original_source = source;
            assert!(prepared.advance(cx, &mut state, &mut source).is_err());
            assert_eq!(state.coordinates(), before.coordinates());
            assert_eq!(state.accepted_time_s(), before.accepted_time_s());
            assert_eq!(source, original_source);
        }
        let foreign = self::pipe(cx);
        let mut foreign_state = foreign
            .steady_state(VolumetricFlowRate::new(0.0), 2.0)
            .unwrap();
        let mut source = [tone(1.0, 0.0), tone(2.0, 0.0)];
        assert!(
            prepared
                .advance(cx, &mut foreign_state, &mut source)
                .is_err()
        );
        assert_eq!(source, [tone(1.0, 0.0), tone(2.0, 0.0)]);
        assert!(prepared.advance(cx, &mut state, &mut source[..1]).is_err());
        // Rejected scratch is reusable without contaminating a later proposal.
        assert!(prepared.advance(cx, &mut state, &mut source).is_ok());
    });
}
