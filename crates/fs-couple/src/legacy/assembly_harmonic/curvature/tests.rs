//! G1 original admitted Maxwell source, complete response, and physical power peak.
use super::*;
use crate::assembly_harmonic::source_frequency_derivative;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_modal::harmonic::resonance::{
    HarmonicPowerCurvatureSample, ResonanceControls, refine_power_extremum_with_cancel,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicPortControls, HarmonicResponseControls,
    PreparedHarmonicResponse,
};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 39,
                kernel_id: 39,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicFrequencyCurvatureControls {
    HarmonicFrequencyCurvatureControls {
        maximum_additional_payload_bytes: 1 << 22,
        maximum_work: 10_000_000,
    }
}
fn source() -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        &MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -1.0,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(11.0 / 8.0)]],
            d: vec![0.0],
            e: vec![0.0],
            reciprocal: true,
        },
        1.0 / 48_000.0,
        MimoLoadLimits {
            band_rad_s: [0.1, 160_000.0],
            sampled_passivity_points: 33,
            max_states: 1,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn prepared(source: &PreparedMimoLoad, omega: f64, cx: &Cx<'_>) -> PreparedHarmonicResponse {
    let z = source.source_model().eval(C64::new(0.0, -omega));
    PreparedHarmonicResponse::try_with_port_impedance(
        &[1.0],
        &[0.0],
        &[1.0],
        1,
        &[1.0],
        1,
        &z,
        omega,
        HarmonicResponseControls {
            maximum_coordinates: 1,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 1e-10,
        },
        HarmonicPortControls {
            maximum_ports: 1,
            maximum_projection_terms: 100_000,
            maximum_payload_bytes: 1 << 20,
        },
        cx,
    )
    .unwrap()
}
fn sample(
    source: &PreparedMimoLoad,
    omega: f64,
    delay: f64,
    cx: &Cx<'_>,
) -> HarmonicPowerCurvatureSample {
    let (zw, zww, remaining) = source_frequency_curvature(source, omega, controls(), cx).unwrap();
    let phase = C64::new((omega * delay).cos(), (omega * delay).sin());
    let f = C64::new(1.2, -0.7) * phase;
    let fw = C64::new(0.0, delay) * f;
    let fww = f.scale(-delay * delay);
    prepared(source, omega, cx)
        .solve_with_frequency_curvature(
            &[f],
            Some(&[fw]),
            Some(&[fww]),
            Some(&zw),
            Some(&zww),
            remaining,
            cx,
        )
        .unwrap()
        .port_power_curvature_sample()
        .unwrap()
}
fn analytic(omega: f64) -> [f64; 3] {
    // Cleared original Maxwell memory equation:
    // ((1-w²)(1-iw)-iw*11/8)q=(1-iw)F.
    // Thus P=N*x/D(x), x=w², D=1+233x/64-15x²/4+x³.
    let x = omega * omega;
    let d = 1.0 + (233.0 / 64.0) * x - (15.0 / 4.0) * x * x + x * x * x;
    let first = 233.0 / 64.0 - (15.0 / 2.0) * x + 3.0 * x * x;
    let second = -15.0 / 2.0 + 6.0 * x;
    let numerator = 0.5 * (11.0 / 8.0) * (1.2 * 1.2 + 0.7 * 0.7);
    let px = numerator * (d - x * first) / (d * d);
    let pxx = numerator * (-x * second / (d * d) - 2.0 * (d - x * first) * first / (d * d * d));
    [
        numerator * x / d,
        2.0 * omega * px,
        2.0 * px + 4.0 * x * pxx,
    ]
}
#[test]
fn g1_admitted_dispersive_source_drives_actual_power_peak_and_curvature() {
    let mut source = source();
    source.step_mean_into(&[0.7], &mut [0.0]).unwrap();
    let history = source.pole_states().to_vec();
    let continuous_class = source.reports().0.class;
    let sampled_class = source.reports().1.class;
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for w in [0.4, 1.1, 2.0_f64.sqrt(), 1.7, 2.4] {
            for delay in [0.0, 0.35] {
                let actual = sample(&source, w, delay, cx);
                for (a, b) in [
                    actual.sample.average_power_w,
                    actual.sample.angular_derivative,
                    actual.angular_second_derivative,
                ]
                .into_iter()
                .zip(analytic(w))
                {
                    assert!((a - b).abs() < 3e-11 * b.abs().max(1.0), "{w}: {a} != {b}");
                }
            }
        }
        let peak = refine_power_extremum_with_cancel(
            (1.2, 1.6),
            ResonanceControls {
                absolute_angular_tolerance_rad_s: 1e-12,
                relative_angular_tolerance: 1e-12,
                maximum_evaluations: 128,
            },
            || false,
            |w| -> Result<_, String> {
                let (zw, remaining) = source_frequency_derivative(
                    &source,
                    w,
                    HarmonicFrequencyControls {
                        maximum_additional_payload_bytes: 1 << 22,
                        maximum_work: 10_000_000,
                    },
                    cx,
                )
                .map_err(|e| e.to_string())?;
                prepared(&source, w, cx)
                    .solve_with_frequency_derivative(
                        &[C64::new(1.2, -0.7)],
                        None,
                        Some(&zw),
                        remaining,
                        cx,
                    )
                    .map_err(|e| e.to_string())?
                    .port_power_sample()
                    .map_err(|e| e.to_string())
            },
            |w| -> Result<_, String> { Ok(sample(&source, w, 0.0, cx)) },
        )
        .unwrap();
        assert!((peak.angular_frequency_rad_s - 2.0_f64.sqrt()).abs() < 5e-11);
        assert!(peak.angular_second_derivative < 0.0);
    });
    assert_eq!(source.pole_states(), history);
    assert_eq!(source.reports().0.class, continuous_class);
    assert_eq!(source.reports().1.class, sampled_class);
}
#[test]
fn source_curvature_retains_band_resources_and_reusable_history() {
    let source = source();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (_, _, remaining) = source_frequency_curvature(&source, 1.0, controls(), cx).unwrap();
        assert_eq!(remaining.maximum_additional_payload_bytes, (1 << 22) - 48);
        assert_eq!(remaining.maximum_work, 10_000_000 - 2048);
        for w in [0.0, 0.05, 160_001.0, f64::NAN] {
            assert!(matches!(
                source_frequency_curvature(&source, w, controls(), cx),
                Err(AssemblyMimoHarmonicError::Source(_))
            ));
        }
        for limited in [
            HarmonicFrequencyCurvatureControls {
                maximum_additional_payload_bytes: 47,
                ..controls()
            },
            HarmonicFrequencyCurvatureControls {
                maximum_work: 2047,
                ..controls()
            },
        ] {
            assert!(source_frequency_curvature(&source, 1.0, limited, cx).is_err());
        }
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            source_frequency_curvature(&source, 1.0, controls(), cx),
            Err(AssemblyMimoHarmonicError::Cancelled)
        ));
    });
}
