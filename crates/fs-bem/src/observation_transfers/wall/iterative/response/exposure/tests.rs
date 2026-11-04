use super::*;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_fft::harmonic_window::HarmonicWindowProfile as WindowProfile;

fn parts(value: C64) -> [Dd; 2] {
    [Dd::from_f64(value.re), Dd::from_f64(value.im)]
}

fn line(omega: f64, pressure: &[C64], drive: &[C64], admittance: &[C64]) -> Line {
    Line {
        omega,
        fields: pressure
            .iter()
            .zip(drive)
            .zip(admittance)
            .map(|((&p, &v), &y)| [parts(p), parts(v - y * p), parts(v), parts(y * p)])
            .collect(),
    }
}

fn bank(statistics: HarmonicExposureFactorStatistics, shift: f64) -> Bank {
    let mut lines = vec![
        line(
            3.1,
            &[C64::new(2.0, -0.4), C64::new(-0.7, 0.9)],
            &[C64::new(0.5, 0.2), C64::new(-0.3, -0.1)],
            &[C64::new(0.08, 0.3), C64::new(0.2, -0.15)],
        ),
        line(
            4.7,
            &[C64::new(-0.9, 0.8), C64::new(1.2, 0.4)],
            &[C64::new(-0.1, 0.4), C64::new(0.6, -0.2)],
            &[C64::new(0.06, -0.2), C64::new(0.13, 0.25)],
        ),
    ];
    if shift != 0.0 {
        for line in &mut lines {
            let rotation = C64::new((line.omega * shift).cos(), (line.omega * shift).sin());
            for field in &mut line.fields {
                for component in field {
                    *component =
                        parts(C64::new(component[0].to_f64(), component[1].to_f64()) * rotation);
                }
            }
        }
    }
    Bank { lines, statistics }
}

fn observer(mean: Bank, factors: Vec<Bank>) -> PreparedBoundaryWorkExposure<'static> {
    PreparedBoundaryWorkExposure {
        source_clock: "retained-field-reference",
        sources: &[],
        areas: &[0.3, 0.7],
        operator_hash: ContentHash([4; 32]),
        mean,
        factors,
        workspace_bytes: 0,
        panel_frequency_pairs: 0,
    }
}
fn zero() -> Bank {
    Bank {
        lines: vec![],
        statistics: HarmonicExposureFactorStatistics::Real,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn powers(value: WindowedBoundaryWork) -> [f64; 3] {
    [
        value.input_power_w,
        value.outward_power_w,
        value.wall_power_w,
    ]
}

// Independent time/phase integration of real boundary traces. Four equally
// weighted phases integrate exactly the proper coordinate's quadratic law.
fn quadrature(
    owner: &PreparedBoundaryWorkExposure<'_>,
    window: HarmonicWindow,
    triangular: bool,
) -> [f64; 3] {
    let steps = 16384;
    let mut result = [0.0; 3];
    for j in 0..=steps {
        let x = j as f64 / steps as f64;
        let t = window.start_s() + window.duration_s() * x;
        let simpson = if j == 0 || j == steps {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let weight = simpson / (3.0 * steps as f64)
            * if triangular {
                4.0 * x.min(1.0 - x)
            } else {
                1.0
            };
        for bank in core::iter::once(&owner.mean).chain(owner.factors.iter()) {
            let phases = if bank.statistics == HarmonicExposureFactorStatistics::Real {
                1
            } else {
                4
            };
            for phase in 0..phases {
                for (panel, area) in owner.areas.iter().enumerate() {
                    let mut field = [0.0; 4];
                    for line in &bank.lines {
                        let angle = line.omega * t - TAU * phase as f64 / phases as f64;
                        for (c, value) in field.iter_mut().enumerate() {
                            *value += line.fields[panel][c][0].to_f64() * angle.cos()
                                + line.fields[panel][c][1].to_f64() * angle.sin();
                        }
                    }
                    for (sum, c) in result.iter_mut().zip([2, 1, 3]) {
                        *sum += weight * area * field[0] * field[c] / phases as f64;
                    }
                }
            }
        }
    }
    result
}

#[test]
fn g1_mixed_boundary_work_matches_independent_time_and_phase_quadrature() {
    let real = HarmonicExposureFactorStatistics::Real;
    let proper = HarmonicExposureFactorStatistics::ProperComplex;
    let owner = observer(bank(real, 0.0), vec![bank(real, 0.31), bank(proper, -0.17)]);
    for triangular in [false, true] {
        let window = HarmonicWindow::with_profile(
            0.23,
            0.71,
            if triangular {
                WindowProfile::Triangular
            } else {
                WindowProfile::Uniform
            },
        )
        .unwrap();
        let actual = owner
            .observe_window(window, &CancelGate::new_clock_free())
            .unwrap();
        for (a, b) in powers(actual)
            .into_iter()
            .zip(quadrature(&owner, window, triangular))
        {
            close(a, b, 5e-12);
        }
        close(actual.balance_w, 0.0, 5e-16);
    }
}

#[test]
fn g1_reactive_wall_releases_stored_energy_without_becoming_negative_heat() {
    let mut owner = observer(zero(), vec![]);
    owner.mean.lines.push(line(
        1.0,
        &[C64::ONE; 2],
        &[C64::ZERO; 2],
        &[C64::new(0.0, 1.0); 2],
    ));
    let gate = CancelGate::new_clock_free();
    let quarter = owner
        .observe_window(HarmonicWindow::new(FRAC_PI_2, FRAC_PI_2).unwrap(), &gate)
        .unwrap();
    close(quarter.wall_power_w, -1.0 / PI, 3e-16);
    close(quarter.outward_power_w, 1.0 / PI, 3e-16);
    close(quarter.input_power_w, 0.0, 0.0);
    close(quarter.balance_w, 0.0, 0.0);
    let cycle = owner
        .observe_window(HarmonicWindow::new(0.0, TAU).unwrap(), &gate)
        .unwrap();
    close(cycle.wall_power_w, 0.0, 2e-16);
    // A proper coordinate removes carrier-sum storage exchange at one tone.
    let mut proper_bank = core::mem::replace(&mut owner.mean, zero());
    proper_bank.statistics = HarmonicExposureFactorStatistics::ProperComplex;
    owner.factors.push(proper_bank);
    let stationary = owner
        .observe_window(HarmonicWindow::new(FRAC_PI_2, FRAC_PI_2).unwrap(), &gate)
        .unwrap();
    close(stationary.wall_power_w, 0.0, 0.0);
}

#[test]
fn g3_window_partition_and_common_clock_shift_preserve_boundary_work() {
    let real = HarmonicExposureFactorStatistics::Real;
    let proper = HarmonicExposureFactorStatistics::ProperComplex;
    let owner = observer(bank(real, 0.0), vec![bank(proper, 0.0)]);
    let gate = CancelGate::new_clock_free();
    let start = 0.4;
    let split = 0.23;
    let duration = 0.79;
    let whole = powers(
        owner
            .observe_window(HarmonicWindow::new(start, duration).unwrap(), &gate)
            .unwrap(),
    );
    let left = powers(
        owner
            .observe_window(HarmonicWindow::new(start, split).unwrap(), &gate)
            .unwrap(),
    );
    let right = powers(
        owner
            .observe_window(
                HarmonicWindow::new(start + split, duration - split).unwrap(),
                &gate,
            )
            .unwrap(),
    );
    for c in 0..3 {
        close(
            whole[c],
            (left[c] * split + right[c] * (duration - split)) / duration,
            5e-16,
        );
    }
    let shift = 0.7;
    let moved = observer(bank(real, shift), vec![bank(proper, shift)]);
    let shifted = powers(
        moved
            .observe_window(HarmonicWindow::new(start + shift, duration).unwrap(), &gate)
            .unwrap(),
    );
    for c in 0..3 {
        close(shifted[c], whole[c], 2e-15);
    }
}

#[test]
fn g3_wall_work_is_accumulated_directly_beside_large_canceling_powers() {
    let p = parts(C64::from_re(1e20));
    let v = parts(C64::from_re(1e20));
    let wall = parts(C64::from_re(1e-20));
    let owner = observer(
        zero(),
        vec![Bank {
            lines: vec![Line {
                omega: 1.0,
                fields: vec![[p, v, v, wall]; 2],
            }],
            statistics: HarmonicExposureFactorStatistics::ProperComplex,
        }],
    );
    let work = owner
        .observe_window(
            HarmonicWindow::new(0.0, 0.01).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(work.input_power_w, work.outward_power_w);
    close(work.wall_power_w, 0.5, 1e-16);
    // The retained traces intentionally omit the weak wall velocity: expose
    // that mismatch, rather than declaring closure by subtracting totals.
    close(work.balance_w, -0.5, 1e-16);
}

#[test]
fn g4_cancellation_leaves_prepared_boundary_work_reusable() {
    let owner = observer(bank(HarmonicExposureFactorStatistics::Real, 0.0), vec![]);
    let window = HarmonicWindow::new(0.2, 0.3).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = owner.observe_window(window, &gate).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(owner.observe_window(window, &cancelled).is_err());
    assert_eq!(owner.observe_window(window, &gate).unwrap(), expected);
    let mut sum = [Dd::ZERO; 2];
    for amplitude in [1e16, 1.0, -1e16] {
        add_product(&mut sum, C64::ONE, C64::from_re(amplitude));
    }
    assert_eq!(sum[0].to_f64(), 1.0);
}
