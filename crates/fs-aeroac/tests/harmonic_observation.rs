//! G0/G1/G3 finite-window interference against analytic and time-domain oracles.

use fs_aeroac::curle2d::PreparedDipoleObserver;
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_aeroac::observation::{HarmonicFieldLine, HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn prepare<const D: usize>(lines: &[HarmonicFieldLine<D>]) -> PreparedHarmonicField<D> {
    PreparedHarmonicField::new(
        lines,
        HarmonicObservationBudget {
            max_lines: lines.len(),
            max_line_pairs: lines.len() * lines.len(),
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn phase(theta: f64) -> C64 {
    let (sin, cos) = theta.sin_cos();
    C64::new(cos, sin)
}

// Independent composite Simpson integration of the actual real time signals.
// This does not use frequency-product identities or the production window API.
fn time_oracle<const D: usize>(
    lines: &[HarmonicFieldLine<D>],
    start: f64,
    duration: f64,
) -> (f64, f64, [f64; D], [f64; D]) {
    let n = 20_000_usize;
    let mut mean_p = 0.0;
    let mut mean_p2 = 0.0;
    let mut mean_u = [0.0; D];
    let mut mean_pu = [0.0; D];
    for step in 0..=n {
        let time = start + duration * step as f64 / n as f64;
        let mut p = 0.0;
        let mut u = [0.0; D];
        for line in lines {
            let rotation = phase(-TAU * line.frequency_hz * time);
            p += (line.pressure_pa * rotation).re;
            for (component, phasor) in u.iter_mut().zip(line.particle_velocity_m_s) {
                *component += (phasor * rotation).re;
            }
        }
        let weight = if step == 0 || step == n {
            1.0
        } else if step % 2 == 0 {
            2.0
        } else {
            4.0
        } / (3.0 * n as f64);
        mean_p += weight * p;
        mean_p2 += weight * p * p;
        for axis in 0..D {
            mean_u[axis] += weight * u[axis];
            mean_pu[axis] += weight * p * u[axis];
        }
    }
    (mean_p, mean_p2, mean_u, mean_pu)
}

#[test]
fn subperiod_window_retains_sum_frequency_and_pressure_velocity_phase() {
    let line = HarmonicFieldLine {
        frequency_hz: 1000.0,
        pressure_pa: C64::from_re(2.0),
        particle_velocity_m_s: [C64::new(0.0, 3.0)],
    };
    let window = HarmonicWindow::new(0.0, 1.0 / 8000.0).unwrap();
    let result = prepare(&[line]).observe_window(window).unwrap();
    // p=2 cos(wt), u=3 sin(wt), wt in [0, pi/4]. The full-period
    // p*u average is zero, but the requested short-window flux is not.
    close(result.mean_pressure_pa, 4.0 * 2.0_f64.sqrt() / PI, 2.0);
    close(result.mean_square_pressure_pa2, 2.0 + 4.0 / PI, 4.0);
    close(result.mean_intensity_w_m2[0], 6.0 / PI, 3.0);
    close(result.pressure_rms_pa(), (2.0 + 4.0 / PI).sqrt(), 2.0);
    assert_eq!(result.window.start_s(), 0.0);
    assert_eq!(result.window.duration_s(), 1.0 / 8000.0);
}

#[test]
fn same_frequency_coalescing_preserves_nulls_and_small_residuals() {
    let mut lines = [HarmonicFieldLine {
        frequency_hz: 1000.0,
        pressure_pa: C64::new(1e16, -2e16),
        particle_velocity_m_s: [C64::new(-1e16, 1e16)],
    }; 3];
    lines[1].pressure_pa = C64::new(1.0, 2.0);
    lines[1].particle_velocity_m_s = [C64::new(3.0, -4.0)];
    lines[2].pressure_pa = -lines[0].pressure_pa;
    lines[2].particle_velocity_m_s = [-lines[0].particle_velocity_m_s[0]];
    let prepared = prepare(&lines);
    assert_eq!(prepared.frequency_count(), 1);
    let result = prepared
        .observe_window(HarmonicWindow::new(0.0, 0.001).unwrap())
        .unwrap();
    close(result.mean_square_pressure_pa2, 2.5, 2.5);
    close(result.mean_intensity_w_m2[0], -2.5, 2.5);

    let null = prepare(&[lines[0], lines[2]])
        .observe_window(HarmonicWindow::new(0.00031, 0.00002).unwrap())
        .unwrap();
    assert_eq!(null.mean_square_pressure_pa2, 0.0);
    assert_eq!(null.mean_intensity_w_m2, [0.0]);
}

#[test]
fn unresolved_tones_beat_in_observer_time_instead_of_adding_powers() {
    let lines = [
        HarmonicFieldLine {
            frequency_hz: 1000.0,
            pressure_pa: C64::ONE,
            particle_velocity_m_s: [C64::new(0.1, 0.2), C64::new(-0.2, 0.3)],
        },
        HarmonicFieldLine {
            frequency_hz: 1005.0,
            pressure_pa: C64::ONE,
            particle_velocity_m_s: [C64::new(-0.3, 0.2), C64::new(0.4, -0.1)],
        },
    ];
    let prepared = prepare(&lines);
    assert_eq!(prepared.frequency_count(), 2);
    let mut powers = Vec::new();
    for start in [0.0, 0.098, 0.2] {
        let duration = 0.004;
        let result = prepared
            .observe_window(HarmonicWindow::new(start, duration).unwrap())
            .unwrap();
        let (p, p2, u, pu) = time_oracle(&lines, start, duration);
        close(result.mean_pressure_pa, p, 2.0);
        close(result.mean_square_pressure_pa2, p2, 2.0);
        for axis in 0..2 {
            close(result.mean_particle_velocity_m_s[axis], u[axis], 1.0);
            close(result.mean_intensity_w_m2[axis], pu[axis], 1.0);
        }
        powers.push(result.mean_square_pressure_pa2);
    }
    // Incoherent line-power addition would give 1 Pa² in all three windows.
    assert!(powers[0] > 1.9 && powers[2] > 1.9);
    assert!(powers[1] < 0.01);
}

#[test]
fn nearly_cancelling_carriers_retain_positive_weak_beat_energy() {
    for (separation_hz, duration) in [(2.0_f64.powi(-20), 0.0001), (2.0_f64.powi(-30), 0.0011)] {
        let lines = [
            HarmonicFieldLine {
                frequency_hz: 1000.0,
                pressure_pa: C64::ONE,
                particle_velocity_m_s: [C64::from_re(1.0 / 400.0)],
            },
            HarmonicFieldLine {
                frequency_hz: 1000.0 + separation_hz,
                pressure_pa: -C64::ONE,
                particle_velocity_m_s: [C64::from_re(-1.0 / 400.0)],
            },
        ];
        // Integrate the stable real beat waveform rather than subtracting two
        // almost equal cosines or using the production frequency-pair kernel.
        let n = 20_000;
        let mut expected = 0.0;
        for step in 0..=n {
            let t = duration * f64::from(step) / f64::from(n);
            let pressure = 2.0
                * (TAU * (1000.0 + 0.5 * separation_hz) * t).sin()
                * (PI * separation_hz * t).sin();
            let weight = if step == 0 || step == n {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            };
            expected += weight * pressure * pressure / (3.0 * f64::from(n));
        }
        let actual = prepare(&lines)
            .observe_window(HarmonicWindow::new(0.0, duration).unwrap())
            .unwrap();
        assert!(actual.mean_square_pressure_pa2 > 0.0);
        assert!(
            (actual.mean_square_pressure_pa2 / expected - 1.0).abs() < 5e-7,
            "beat energy {:?} versus {expected:e}",
            actual.mean_square_pressure_pa2
        );
        assert!((actual.mean_intensity_w_m2[0] / (expected / 400.0) - 1.0).abs() < 5e-7);
    }
}

#[test]
fn common_clock_shift_preserves_the_finite_observation() {
    let lines = [
        HarmonicFieldLine {
            frequency_hz: 257.0,
            pressure_pa: C64::new(1.2, -0.7),
            particle_velocity_m_s: [C64::new(0.4, 0.9)],
        },
        HarmonicFieldLine {
            frequency_hz: 263.0,
            pressure_pa: C64::new(-0.2, 0.8),
            particle_velocity_m_s: [C64::new(-0.3, 0.2)],
        },
    ];
    let shift = 0.041;
    let shifted = lines.map(|line| {
        let rotation = phase(-TAU * line.frequency_hz * shift);
        HarmonicFieldLine {
            pressure_pa: line.pressure_pa * rotation,
            particle_velocity_m_s: line.particle_velocity_m_s.map(|v| v * rotation),
            ..line
        }
    });
    let a = prepare(&lines)
        .observe_window(HarmonicWindow::new(0.13, 0.017).unwrap())
        .unwrap();
    let b = prepare(&shifted)
        .observe_window(HarmonicWindow::new(0.13 - shift, 0.017).unwrap())
        .unwrap();
    close(a.mean_pressure_pa, b.mean_pressure_pa, 2.0);
    close(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2, 2.0);
    close(a.mean_intensity_w_m2[0], b.mean_intensity_w_m2[0], 2.0);
}

#[test]
fn real_3d_source_propagation_feeds_multitone_observation() {
    let medium = FluidAcousticState::try_new(1.3, 327.0).unwrap();
    let positions = [[-0.2, 0.1, 0.05], [0.4, -0.3, 0.2]];
    let sources = [
        CompactSourcePhasor::monopole(C64::new(0.0002, -0.0001)),
        CompactSourcePhasor::force_dipole([
            C64::new(0.3, 0.2),
            C64::new(-0.1, 0.4),
            C64::new(0.1, -0.2),
        ]),
    ];
    let lines = [211.0, 223.0, 419.0].map(|frequency| {
        let observer =
            PreparedCompactSourceObserver::new(medium, frequency, [0.2, 0.3, 1.0], &positions, 2)
                .unwrap();
        let line = observer.evaluate_harmonic_line(&sources).unwrap();
        let field = observer.evaluate(&sources).unwrap();
        assert_eq!(line.pressure_pa, field.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
        line
    });
    let window = HarmonicWindow::new(0.007, 0.012).unwrap();
    let result = prepare(&lines).observe_window(window).unwrap();
    let (p, p2, u, pu) = time_oracle(&lines, window.start_s(), window.duration_s());
    close(result.mean_pressure_pa, p, 1.0);
    close(result.mean_square_pressure_pa2, p2, 1.0);
    for axis in 0..3 {
        close(result.mean_particle_velocity_m_s[axis], u[axis], 0.01);
        close(result.mean_intensity_w_m2[axis], pu[axis], 0.01);
    }
}

#[test]
fn real_2d_observer_recovers_existing_cycle_average() {
    let frequency = 347.0;
    let observer = PreparedDipoleObserver::new(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        frequency,
        [0.3, 0.7],
        &[[0.0, 0.0], [-0.2, 0.1]],
        2,
    )
    .unwrap();
    let forces = [
        [C64::new(0.2, -0.1), C64::new(0.5, 0.1)],
        [C64::new(-0.1, 0.3), C64::new(0.2, 0.4)],
    ];
    let field = observer.evaluate(&forces).unwrap();
    let line = observer.evaluate_harmonic_line(&forces).unwrap();
    let result = prepare(&[line])
        .observe_window(HarmonicWindow::new(0.00017, 3.0 / frequency).unwrap())
        .unwrap();
    close(
        result.mean_square_pressure_pa2,
        field.pressure_pa.norm_sq() / 2.0,
        1.0,
    );
    for axis in 0..2 {
        close(
            result.mean_intensity_w_m2[axis],
            field.active_intensity_w_m2[axis],
            0.01,
        );
    }
}

#[test]
fn admission_budgets_invalid_lines_and_silence() {
    let line = HarmonicFieldLine {
        frequency_hz: 100.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::ONE],
    };
    for budget in [
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 4,
        },
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 3,
        },
    ] {
        assert!(PreparedHarmonicField::new(&[line, line], budget).is_err());
    }
    for frequency_hz in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::MAX] {
        assert!(
            PreparedHarmonicField::new(
                &[HarmonicFieldLine {
                    frequency_hz,
                    ..line
                }],
                HarmonicObservationBudget {
                    max_lines: 1,
                    max_line_pairs: 1
                }
            )
            .is_err()
        );
    }
    assert!(
        PreparedHarmonicField::new(
            &[HarmonicFieldLine {
                pressure_pa: C64::new(f64::NAN, 0.0),
                ..line
            }],
            HarmonicObservationBudget {
                max_lines: 1,
                max_line_pairs: 1
            }
        )
        .is_err()
    );
    let silence = prepare::<3>(&[])
        .observe_window(HarmonicWindow::new(-1.0, 0.5).unwrap())
        .unwrap();
    assert_eq!(silence.mean_pressure_pa, 0.0);
    assert_eq!(silence.mean_square_pressure_pa2, 0.0);
    assert_eq!(silence.mean_particle_velocity_m_s, [0.0; 3]);
    assert_eq!(silence.mean_intensity_w_m2, [0.0; 3]);
}
