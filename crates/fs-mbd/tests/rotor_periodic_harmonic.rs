//! G1/G3 continuous unbalance, exact resonant growth and switched-coefficient forcing.
use fs_mbd::rotordynamics::periodic::{
    MAX_ROTOR_FORCE_HARMONIC, RotorFloquetBudget, RotorHarmonicForce, RotorPeriodicSegment,
    periodic_floquet, prepare_harmonic_cycle,
};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::f64::consts::{PI, TAU};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn oscillator(mass: f64, damping: f64, stiffness: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[mass], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

#[test]
fn resonant_continuous_drive_matches_analytic_orbit_across_scales_and_subdivision() {
    let model = oscillator(1.0, 0.02, 1.0);
    for amplitude in [1e-80, 1.0, 1e80] {
        let cosine = [amplitude];
        let force = RotorHarmonicForce {
            cosine: &cosine,
            sine: &[0.0],
        };
        for fractions in [&[1.0][..], &[0.17, 0.23, 0.60][..]] {
            let segments: Vec<_> = fractions
                .iter()
                .map(|fraction| RotorPeriodicSegment {
                    model: &model,
                    duration_seconds: fraction * TAU,
                })
                .collect();
            let loads = vec![force; segments.len()];
            let plan = prepare_harmonic_cycle(
                1.0,
                &segments,
                &loads,
                1,
                1e-12,
                RotorFloquetBudget::default(),
                || false,
            )
            .unwrap();
            let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
            // At omega=1, q=F/C*sin(t), v=F/C*cos(t), without a held-load approximation.
            for (&t, state) in orbit
                .trajectory
                .boundary_times_seconds
                .iter()
                .zip(&orbit.trajectory.boundary_states)
            {
                near(state.displacement[0] / amplitude, t.sin() / 0.02, 2e-9);
                near(state.velocity[0] / amplitude, t.cos() / 0.02, 2e-9);
            }
            assert_eq!(plan.floquet_report().multipliers.len(), 2);
            near(
                plan.floquet_report().spectral_radius,
                (-0.01 * TAU).exp(),
                1e-11,
            );
            assert!(orbit.relative_fixed_point_residual < 1e-9);
            assert!(orbit.relative_closure_residual < 1e-9);
        }
    }
}

#[test]
fn rotating_unbalance_retains_the_signed_gyroscopic_shift_and_quadrature() {
    let (mass, damping, stiffness, gyro, frequency, force) = (2.0, 0.2, 8.0, 0.4, 2.0, 0.5);
    let model = RotorBearingModel::new(
        2,
        &[mass, 0.0, 0.0, mass],
        &[damping, 0.0, 0.0, damping],
        &[stiffness, 0.0, 0.0, stiffness],
        &[0.0, gyro, -gyro, 0.0],
        || false,
    )
    .unwrap();
    for spin in [-3.0, 3.0] {
        let plan = prepare_harmonic_cycle(
            spin,
            &[RotorPeriodicSegment {
                model: &model,
                duration_seconds: TAU / frequency,
            }],
            &[RotorHarmonicForce {
                cosine: &[force, 0.0],
                sine: &[0.0, force],
            }],
            1,
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
        // [qx,qy]=Re([Q,-iQ]exp(i wt)); D=k-m*w^2+Omega*g*w+i*c*w.
        let a = stiffness - mass * frequency * frequency + spin * gyro * frequency;
        let b = damping * frequency;
        let in_phase = force * a / (a * a + b * b);
        let quadrature = force * b / (a * a + b * b);
        for state in &orbit.trajectory.boundary_states {
            near(state.displacement[0], in_phase, 2e-10);
            near(state.displacement[1], -quadrature, 2e-10);
            near(state.velocity[0], frequency * quadrature, 2e-10);
            near(state.velocity[1], frequency * in_phase, 2e-10);
        }
        assert_eq!(plan.floquet_report().multipliers.len(), 4);
    }
}

#[test]
fn full_mutual_inertia_and_free_coordinates_follow_the_integrated_physical_waveform() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.0],
        &[0.0; 4],
        &[0.0; 4],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    let frequency = 1.7;
    let period = TAU / frequency;
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: period / 4.0,
    }; 4];
    let forces = [RotorHarmonicForce {
        cosine: &[1.3, 0.4],
        sine: &[0.2, -0.7],
    }; 4];
    let plan = prepare_harmonic_cycle(
        0.0,
        &segments,
        &forces,
        1,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let motion = plan
        .advance(
            &RotorTransientState {
                displacement: vec![0.0; 2],
                velocity: vec![0.0; 2],
            },
            1e-10,
            || false,
        )
        .unwrap();
    // Original M^-1 f from the exact two-by-two inverse, then integrate twice.
    let cosine_acceleration = [(1.3 - 0.5 * 0.4) / 1.75, (-0.5 * 1.3 + 2.0 * 0.4) / 1.75];
    let sine_acceleration = [(0.2 + 0.5 * 0.7) / 1.75, (-0.5 * 0.2 - 2.0 * 0.7) / 1.75];
    for (&t, state) in motion
        .boundary_times_seconds
        .iter()
        .zip(&motion.boundary_states)
    {
        let (s, c) = (frequency * t).sin_cos();
        for i in 0..2 {
            let q = cosine_acceleration[i] * (1.0 - c) / frequency.powi(2)
                + sine_acceleration[i] * (t / frequency - s / frequency.powi(2));
            let v = cosine_acceleration[i] * s / frequency
                + sine_acceleration[i] * (1.0 - c) / frequency;
            near(state.displacement[i], q, 2e-11);
            near(state.velocity[i], v, 2e-11);
        }
    }
}

#[test]
fn lossless_resonance_keeps_secular_motion_and_refuses_a_false_periodic_solution() {
    let model = oscillator(1.0, 0.0, 4.0);
    let plan = prepare_harmonic_cycle(
        2.0,
        &[RotorPeriodicSegment {
            model: &model,
            duration_seconds: PI,
        }],
        &[RotorHarmonicForce {
            cosine: &[1.3],
            sine: &[-0.4],
        }],
        1,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let motion = plan
        .advance(
            &RotorTransientState {
                displacement: vec![0.0],
                velocity: vec![0.0],
            },
            1e-10,
            || false,
        )
        .unwrap();
    // Resonance: q_cos=Fc*t*sin(wt)/(2w),
    // q_sin=Fs*(sin(wt)/w^2-t*cos(wt)/w)/2, w=2.
    let end = motion.boundary_states.last().unwrap();
    near(end.displacement[0], 0.4 * PI / 4.0, 2e-11);
    near(end.velocity[0], 1.3 * PI / 2.0, 2e-11);
    assert!(matches!(
        plan.periodic_orbit(1e-9, 1e8, || false),
        Err(RotorError::PeriodicOrbitUnresolved { .. } | RotorError::PeriodicLinear(_))
    ));
}

// Independent scalar oscillator solution with harmonic particular motion and
// the homogeneous underdamped closed form. No production exponential is used.
fn analytic_segment(
    mass: f64,
    damping: f64,
    stiffness: f64,
    frequency: f64,
    cosine: f64,
    sine: f64,
    start: f64,
    duration: f64,
    state: [f64; 2],
) -> [f64; 2] {
    let a = stiffness - mass * frequency.powi(2);
    let b = damping * frequency;
    let real = (cosine * a - sine * b) / (a * a + b * b);
    let imaginary = -(cosine * b + sine * a) / (a * a + b * b);
    let particular = |t: f64| {
        let (s, c) = (frequency * t).sin_cos();
        [
            real * c - imaginary * s,
            -frequency * (real * s + imaginary * c),
        ]
    };
    let initial = particular(start);
    let end = particular(start + duration);
    let alpha = damping / (2.0 * mass);
    let wd = (stiffness / mass - alpha.powi(2)).sqrt();
    let (s, c) = (wd * duration).sin_cos();
    let e = (-alpha * duration).exp();
    let q = state[0] - initial[0];
    let v = state[1] - initial[1];
    [
        end[0] + e * ((c + alpha * s / wd) * q + s * v / wd),
        end[1] + e * (-stiffness * s * q / (mass * wd) + (c - alpha * s / wd) * v),
    ]
}

#[test]
fn switched_bearings_and_loads_keep_one_global_phase_at_every_boundary() {
    let soft = oscillator(1.3, 0.07, 4.0);
    let stiff = oscillator(1.3, 0.07, 9.0);
    let segments = [
        RotorPeriodicSegment {
            model: &soft,
            duration_seconds: 0.41,
        },
        RotorPeriodicSegment {
            model: &stiff,
            duration_seconds: 0.69,
        },
    ];
    let loads = [
        RotorHarmonicForce {
            cosine: &[0.7],
            sine: &[0.3],
        },
        RotorHarmonicForce {
            cosine: &[-0.2],
            sine: &[0.8],
        },
    ];
    let plan = prepare_harmonic_cycle(
        -5.0,
        &segments,
        &loads,
        3,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let motion = plan
        .advance(
            &RotorTransientState {
                displacement: vec![0.2],
                velocity: vec![-0.3],
            },
            1e-10,
            || false,
        )
        .unwrap();
    let frequency = 3.0 * TAU / 1.1;
    let middle = analytic_segment(1.3, 0.07, 4.0, frequency, 0.7, 0.3, 0.0, 0.41, [0.2, -0.3]);
    let end = analytic_segment(1.3, 0.07, 9.0, frequency, -0.2, 0.8, 0.41, 0.69, middle);
    for (actual, expected) in motion.boundary_states[1..].iter().zip([middle, end]) {
        near(actual.displacement[0], expected[0], 3e-10);
        near(actual.velocity[0], expected[1], 3e-10);
    }
}

#[test]
fn zero_harmonic_load_shares_homogeneous_floquet_and_admits_no_invalid_schedule() {
    let model = oscillator(1.7, 0.2, 3.1);
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.4,
    }; 2];
    let zero = [RotorHarmonicForce {
        cosine: &[0.0],
        sine: &[0.0],
    }; 2];
    let budget = RotorFloquetBudget::default();
    let homogeneous = periodic_floquet(0.0, &segments, 1e-12, budget, || false).unwrap();
    let plan = prepare_harmonic_cycle(0.0, &segments, &zero, 7, 1e-12, budget, || false).unwrap();
    assert_eq!(plan.floquet_report(), &homogeneous);
    for order in [0, MAX_ROTOR_FORCE_HARMONIC + 1] {
        assert!(
            prepare_harmonic_cycle(0.0, &segments, &zero, order, 1e-12, budget, || false).is_err()
        );
    }
    assert!(
        prepare_harmonic_cycle(0.0, &segments, &zero[..1], 1, 1e-12, budget, || false).is_err()
    );
    let invalid = [RotorHarmonicForce {
        cosine: &[f64::NAN],
        sine: &[0.0],
    }; 2];
    assert!(prepare_harmonic_cycle(0.0, &segments, &invalid, 1, 1e-12, budget, || false).is_err());
    let mismatch = [RotorHarmonicForce {
        cosine: &[0.0],
        sine: &[],
    }; 2];
    assert!(prepare_harmonic_cycle(0.0, &segments, &mismatch, 1, 1e-12, budget, || false).is_err());
    assert!(matches!(
        prepare_harmonic_cycle(
            0.0,
            &segments,
            &zero,
            1,
            1e-12,
            RotorFloquetBudget {
                max_dense_bytes: 1,
                ..budget
            },
            || false
        ),
        Err(RotorError::PeriodicBudget(_))
    ));
    assert!(matches!(
        prepare_harmonic_cycle(0.0, &segments, &zero, 1, 1e-12, budget, || true),
        Err(RotorError::Cancelled)
    ));
}
