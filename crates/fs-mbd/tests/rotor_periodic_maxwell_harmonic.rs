//! G1/G3 continuous physical loads on every original Maxwell memory state.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce, RotorMaxwellState,
    RotorPeriodicSegment, prepare_maxwell_cycle, prepare_maxwell_forced_cycle,
    prepare_maxwell_fourier_cycle, prepare_maxwell_harmonic_cycle,
};
use fs_mbd::rotordynamics::relaxation::MaxwellDampingBranch;
use std::f64::consts::{PI, TAU};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn arm(strain: &[f64], stiffness: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: strain.to_vec(),
        stiffness,
        relaxation_time_s: tau,
    }
}

fn oscillator(m: f64, c: f64, k: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap()
}

// Independent original scalar constitutive algebra, not a period map or
// reduced viscous approximation: sigma_hat=k*i*w*tau/(1+i*w*tau)*b*q_hat.
fn scalar_phasors(
    mechanical: (f64, f64, f64),
    branches: &[MaxwellDampingBranch],
    omega: f64,
    force: C64,
) -> (C64, C64, Vec<C64>) {
    let (m, c, k) = mechanical;
    let transfer: Vec<_> = branches
        .iter()
        .map(|branch| {
            C64::new(0.0, omega * branch.relaxation_time_s)
                / C64::new(1.0, omega * branch.relaxation_time_s)
                * C64::from_re(branch.stiffness)
        })
        .collect();
    let mut dynamic = C64::new(k - m * omega * omega, c * omega);
    for (branch, h) in branches.iter().zip(&transfer) {
        dynamic = dynamic + h.scale(branch.strain[0] * branch.strain[0]);
    }
    let q = force / dynamic;
    let v = C64::new(0.0, omega) * q;
    let stress = branches
        .iter()
        .zip(transfer)
        .map(|(branch, h)| (h * q).scale(branch.strain[0]))
        .collect();
    (q, v, stress)
}

fn physical(value: C64, omega: f64, time: f64) -> f64 {
    let (sin, cos) = (omega * time).sin_cos();
    (value * C64::new(cos, sin)).re
}

#[test]
fn harmonic_orbit_preserves_original_stress_phase_and_passive_heat_across_force_units() {
    let (m, c, omega) = (1.3_f64, 0.04_f64, 1.7_f64);
    let branches = [
        arm(&[1.0], 0.45, 0.7),
        arm(&[-0.6], 0.3, 1.1),
        arm(&[0.0], 0.8, 1.3),
    ];
    // Cancel the real dynamic stiffness at this declared drive frequency.
    let storage: f64 = branches
        .iter()
        .map(|branch| {
            let a = omega * branch.relaxation_time_s;
            branch.stiffness * a * a / (1.0 + a * a) * branch.strain[0].powi(2)
        })
        .sum();
    let k = m * omega * omega - storage;
    let model = oscillator(m, c, k);
    let period = TAU / omega;
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: period / 8.0,
    }; 8];
    let (q, v, stress) = scalar_phasors((m, c, k), &branches, omega, C64::new(0.6, -0.8));
    let expected_heat = 0.5 * c * v.norm_sq()
        + branches
            .iter()
            .zip(&stress)
            .map(|(branch, sigma)| sigma.norm_sq() / (2.0 * branch.stiffness * branch.relaxation_time_s))
            .sum::<f64>();
    for scale in [1e-80, 1.0, 1e80] {
        let cosine = [0.6 * scale];
        let sine = [0.8 * scale];
        let forces = [RotorHarmonicForce { cosine: &cosine, sine: &sine }; 8];
        let plan = prepare_maxwell_harmonic_cycle(
            -2.7, &segments, &branches, &forces, 1, 1e-12,
            RotorFloquetBudget::default(), || false,
        )
        .unwrap();
        let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
        let mut heat = 0.0;
        let mut supply = 0.0;
        for (index, (&time, state)) in orbit.trajectory.boundary_times_seconds
            .iter().zip(&orbit.trajectory.boundary_states).enumerate()
        {
            near(state.displacement[0] / scale, physical(q, omega, time), 3e-9);
            near(state.velocity[0] / scale, physical(v, omega, time), 3e-9);
            for (actual, sigma) in state.branch_stress.iter().zip(&stress) {
                near(actual / scale, physical(*sigma, omega, time), 3e-9);
            }
            if index < 8 {
                // Uniform periodic quadrature exactly resolves the constant
                // and second-harmonic terms of this single-tone power oracle.
                let velocity = state.velocity[0] / scale;
                heat += c * velocity * velocity;
                for (branch, sigma) in branches.iter().zip(&state.branch_stress) {
                    heat += (sigma / scale).powi(2)
                        / (branch.stiffness * branch.relaxation_time_s);
                }
                supply += physical(C64::new(0.6, -0.8), omega, time) * velocity;
            }
        }
        near(heat / 8.0, expected_heat, 5e-9);
        near(supply / 8.0, expected_heat, 5e-9);
        assert_eq!(plan.floquet_report().multipliers.len(), 5);
        assert_eq!(plan.floquet_report().memory_states, 3);
        assert!(plan.floquet_report().maximum_memory_residual < 1e-12);
        assert!(orbit.relative_fixed_point_residual < 1e-9);
        assert!(orbit.relative_closure_residual < 1e-9);
    }
}

#[test]
fn dc_and_simultaneous_tones_recover_every_material_stress_without_extra_modes() {
    let mechanical = (1.5, 0.13, 4.0);
    let model = oscillator(mechanical.0, mechanical.1, mechanical.2);
    let branches = [arm(&[1.2], 0.7, 0.8), arm(&[-0.4], 1.1, 1.3)];
    let fundamental = 1.9;
    let segments: Vec<_> = [0.13, 0.27, 0.6].iter().map(|fraction| RotorPeriodicSegment {
        model: &model,
        duration_seconds: fraction * TAU / fundamental,
    }).collect();
    let orders = [1, 3, 7];
    let cosine = [[0.4], [-0.5], [0.2]];
    let sine = [[0.2], [-0.3], [-0.7]];
    let tones: Vec<_> = (0..3).map(|i| RotorHarmonicForce {
        cosine: &cosine[i], sine: &sine[i],
    }).collect();
    let loads = [RotorFourierForce { steady: &[0.9], harmonics: &tones }; 3];
    let plan = prepare_maxwell_fourier_cycle(
        4.2, &segments, &branches, &loads, &orders, 1e-12,
        RotorFloquetBudget::default(), || false,
    ).unwrap();
    let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
    for (&time, state) in orbit.trajectory.boundary_times_seconds
        .iter().zip(&orbit.trajectory.boundary_states)
    {
        let (mut q, mut v) = (0.9 / mechanical.2, 0.0);
        let mut expected_stress = [0.0; 2];
        for (i, &order) in orders.iter().enumerate() {
            let omega = order as f64 * fundamental;
            let (qh, vh, sh) = scalar_phasors(
                mechanical, &branches, omega, C64::new(cosine[i][0], -sine[i][0]),
            );
            q += physical(qh, omega, time);
            v += physical(vh, omega, time);
            for j in 0..2 {
                expected_stress[j] += physical(sh[j], omega, time);
            }
        }
        near(state.displacement[0], q, 3e-9);
        near(state.velocity[0], v, 3e-9);
        for (&actual, expected) in state.branch_stress.iter().zip(expected_stress) {
            near(actual, expected, 3e-9);
        }
    }
    let homogeneous = prepare_maxwell_cycle(
        4.2, &segments, &branches, 1e-12, RotorFloquetBudget::default(), || false,
    ).unwrap();
    assert_eq!(plan.floquet_report().multipliers.len(), 4);
    for (&actual, &expected) in plan.floquet_report().monodromy.iter()
        .zip(&homogeneous.floquet_report().monodromy)
    {
        near(actual, expected, 3e-12);
    }
}

// Independent RK4 of the ORIGINAL four scalar ODEs, retaining each sigma.
// Time is global across a coefficient/amplitude switch.
#[allow(clippy::too_many_arguments)]
fn integrate_segment(
    mut state: [f64; 4],
    mechanical: (f64, f64, f64),
    branches: &[MaxwellDampingBranch; 2],
    force: impl Fn(f64) -> f64,
    start: f64,
    duration: f64,
    steps: usize,
) -> [f64; 4] {
    let (m, c, k) = mechanical;
    let derivative = |time: f64, x: [f64; 4]| [
        x[1],
        (force(time) - c * x[1] - k * x[0]
            - branches[0].strain[0] * x[2] - branches[1].strain[0] * x[3]) / m,
        branches[0].stiffness * branches[0].strain[0] * x[1] - x[2] / branches[0].relaxation_time_s,
        branches[1].stiffness * branches[1].strain[0] * x[1] - x[3] / branches[1].relaxation_time_s,
    ];
    let h = duration / steps as f64;
    for step in 0..steps {
        let t = start + step as f64 * h;
        let a = derivative(t, state);
        let b = derivative(t + h / 2.0, std::array::from_fn(|i| state[i] + h * a[i] / 2.0));
        let c = derivative(t + h / 2.0, std::array::from_fn(|i| state[i] + h * b[i] / 2.0));
        let d = derivative(t + h, std::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn switched_fourier_forces_keep_global_phase_and_arbitrary_original_memory() {
    let first = oscillator(1.2, 0.07, 2.4);
    let second = oscillator(1.2, 0.13, 5.7);
    let segments = [
        RotorPeriodicSegment { model: &first, duration_seconds: 0.37 },
        RotorPeriodicSegment { model: &second, duration_seconds: 0.63 },
    ];
    let branches = [arm(&[1.0], 0.5, 0.5), arm(&[-0.4], 0.9, 1.1)];
    let dc = [[0.2], [-0.3]];
    let cosine = [[[0.4], [-0.7]], [[-0.2], [0.6]]];
    let sine = [[[0.3], [0.1]], [[0.5], [-0.4]]];
    let tones: Vec<Vec<_>> = (0..2).map(|i| (0..2).map(|j| RotorHarmonicForce {
        cosine: &cosine[i][j], sine: &sine[i][j],
    }).collect()).collect();
    let loads: Vec<_> = (0..2).map(|i| RotorFourierForce {
        steady: &dc[i], harmonics: &tones[i],
    }).collect();
    let orders = [1, 3];
    let plan = prepare_maxwell_fourier_cycle(
        -3.0, &segments, &branches, &loads, &orders, 1e-12,
        RotorFloquetBudget::default(), || false,
    ).unwrap();
    let initial = RotorMaxwellState {
        displacement: vec![0.2], velocity: vec![-0.3], branch_stress: vec![0.7, -0.2],
    };
    let trajectory = plan.advance(&initial, 1e-10, || false).unwrap();
    let mut reference = [0.2, -0.3, 0.7, -0.2];
    let mut start = 0.0;
    for i in 0..2 {
        let force = |time: f64| dc[i][0] + (0..2).map(|j| {
            let (sin, cos) = (TAU * orders[j] as f64 * time).sin_cos();
            cosine[i][j][0] * cos + sine[i][j][0] * sin
        }).sum::<f64>();
        let mechanics = if i == 0 { (1.2, 0.07, 2.4) } else { (1.2, 0.13, 5.7) };
        let h = segments[i].duration_seconds;
        let coarse = integrate_segment(reference, mechanics, &branches, force, start, h, 1024);
        reference = integrate_segment(reference, mechanics, &branches, force, start, h, 2048);
        let state = &trajectory.boundary_states[i + 1];
        let actual = [state.displacement[0], state.velocity[0], state.branch_stress[0], state.branch_stress[1]];
        for j in 0..4 {
            near(coarse[j], reference[j], 1e-9);
            near(actual[j], reference[j], 1e-9);
        }
        start += h;
    }
    assert_eq!(trajectory.boundary_states[0], initial);
    assert_eq!(plan.floquet_report().multipliers.len(), 4);
}

#[test]
fn harmonic_full_inertia_signed_gyro_and_coupled_stresses_match_original_complex_equations() {
    let mass = [2.0, 0.3, 0.3, 1.4];
    let damping = [0.17, 0.04, -0.02, 0.22];
    let stiffness = [3.2, 0.12, 0.08, 4.1];
    let gyro = [0.0, 0.5, -0.5, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    let branches = [arm(&[1.0, -0.3], 0.4, 0.8), arm(&[-0.2, 0.8], 0.7, 1.3)];
    let omega = 1.9;
    let force = [C64::new(0.6, -0.2), C64::new(-0.3, 0.5)];
    let segments: Vec<_> = [0.21, 0.34, 0.45].iter().map(|fraction| RotorPeriodicSegment {
        model: &model, duration_seconds: fraction * TAU / omega,
    }).collect();
    let sources = [RotorHarmonicForce { cosine: &[0.6, -0.3], sine: &[0.2, -0.5] }; 3];
    for spin in [-3.0, 3.0] {
        let mut dynamic = std::array::from_fn::<_, 4, _>(|i| C64::new(
            stiffness[i] - omega * omega * mass[i],
            omega * (damping[i] + spin * gyro[i]),
        ));
        let transfer: Vec<_> = branches.iter().map(|branch| {
            C64::new(0.0, omega * branch.relaxation_time_s)
                / C64::new(1.0, omega * branch.relaxation_time_s)
                * C64::from_re(branch.stiffness)
        }).collect();
        for (branch, h) in branches.iter().zip(&transfer) {
            for i in 0..2 {
                for j in 0..2 {
                    dynamic[i * 2 + j] = dynamic[i * 2 + j] + h.scale(branch.strain[i] * branch.strain[j]);
                }
            }
        }
        let determinant = dynamic[0] * dynamic[3] - dynamic[1] * dynamic[2];
        let q = [
            (force[0] * dynamic[3] - dynamic[1] * force[1]) / determinant,
            (dynamic[0] * force[1] - force[0] * dynamic[2]) / determinant,
        ];
        let v = q.map(|value| C64::new(0.0, omega) * value);
        let stress: Vec<_> = branches.iter().zip(transfer).map(|(branch, h)| {
            h * (q[0].scale(branch.strain[0]) + q[1].scale(branch.strain[1]))
        }).collect();
        let plan = prepare_maxwell_harmonic_cycle(
            spin, &segments, &branches, &sources, 1, 1e-12,
            RotorFloquetBudget::default(), || false,
        ).unwrap();
        let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
        for (&time, state) in orbit.trajectory.boundary_times_seconds
            .iter().zip(&orbit.trajectory.boundary_states)
        {
            for i in 0..2 {
                near(state.displacement[i], physical(q[i], omega, time), 3e-9);
                near(state.velocity[i], physical(v[i], omega, time), 3e-9);
                near(state.branch_stress[i], physical(stress[i], omega, time), 3e-9);
            }
        }
        assert_eq!(plan.floquet_report().multipliers.len(), 6);
        assert!(plan.floquet_report().maximum_mass_residual < 1e-12);
        assert!(plan.floquet_report().maximum_memory_residual < 1e-12);
    }
}

#[test]
fn lossless_secular_motion_and_uncoupled_material_decay_survive_together() {
    let model = oscillator(1.0, 0.0, 4.0);
    let branches = [arm(&[0.0], 2.0, 0.4)];
    let segments = [RotorPeriodicSegment { model: &model, duration_seconds: PI / 4.0 }; 4];
    let forces = [RotorHarmonicForce { cosine: &[1.3], sine: &[-0.4] }; 4];
    let plan = prepare_maxwell_harmonic_cycle(
        2.0, &segments, &branches, &forces, 1, 1e-12,
        RotorFloquetBudget::default(), || false,
    ).unwrap();
    let initial = RotorMaxwellState {
        displacement: vec![0.2], velocity: vec![-0.1], branch_stress: vec![0.6],
    };
    let motion = plan.advance(&initial, 1e-9, || false).unwrap();
    for (&t, state) in motion.boundary_times_seconds.iter().zip(&motion.boundary_states) {
        let (sin, cos) = (2.0 * t).sin_cos();
        let q = 0.2 * cos - 0.05 * sin
            + 1.3 * t * sin / 4.0 - 0.2 * (sin / 4.0 - t * cos / 2.0);
        let v = -0.4 * sin - 0.1 * cos
            + 1.3 * (sin / 4.0 + t * cos / 2.0) - 0.2 * t * sin;
        near(state.displacement[0], q, 3e-10);
        near(state.velocity[0], v, 3e-10);
        near(state.branch_stress[0], 0.6 * (-t / 0.4).exp(), 3e-10);
    }
    assert_eq!(plan.floquet_report().multipliers.len(), 3);
    // Both mechanical increments are resolved because Fc AND Fs are nonzero.
    // A tiny rounded I-P must not manufacture a finite periodic resonance.
    assert!(matches!(
        plan.periodic_orbit(1e-9, 1e8, || false),
        Err(RotorError::PeriodicOrbitUnresolved { .. } | RotorError::PeriodicLinear(_))
    ));
}

#[test]
fn zero_tones_reuse_held_forcing_and_admission_cancellation_publish_no_partial_cycle() {
    let model = oscillator(1.0, 0.2, 4.0);
    let branches = [arm(&[1.0], 0.7, 0.9)];
    let segments = [RotorPeriodicSegment { model: &model, duration_seconds: 0.4 }];
    let rows = [RotorFourierForce { steady: &[0.3], harmonics: &[] }];
    let budget = RotorFloquetBudget::default();
    let fourier = prepare_maxwell_fourier_cycle(
        0.0, &segments, &branches, &rows, &[], 1e-12, budget, || false,
    ).unwrap();
    let held = prepare_maxwell_forced_cycle(
        0.0, &segments, &branches, &[&[0.3]], 1e-12, budget, || false,
    ).unwrap();
    assert_eq!(fourier.forcing_offset(), held.forcing_offset());
    assert_eq!(fourier.floquet_report(), held.floquet_report());
    let tones = [RotorHarmonicForce { cosine: &[0.4], sine: &[-0.2] }];
    assert!(matches!(
        prepare_maxwell_harmonic_cycle(0.0, &segments, &branches, &tones, 0, 1e-12, budget, || false),
        Err(RotorError::Invalid(_))
    ));
    let wrong = [RotorFourierForce { steady: &[0.0], harmonics: &tones }];
    assert!(matches!(
        prepare_maxwell_fourier_cycle(0.0, &segments, &branches, &wrong, &[2, 1], 1e-12, budget, || false),
        Err(RotorError::Invalid(_))
    ));
    let small = RotorFloquetBudget { max_multiply_terms: 1, ..budget };
    assert!(matches!(
        prepare_maxwell_harmonic_cycle(0.0, &segments, &branches, &tones, 1, 1e-12, small, || false),
        Err(RotorError::PeriodicBudget(_))
    ));
    assert!(matches!(
        prepare_maxwell_harmonic_cycle(0.0, &segments, &branches, &tones, 1, 1e-12, budget, || true),
        Err(RotorError::Cancelled)
    ));
    let mut checks = 0;
    prepare_maxwell_harmonic_cycle(0.0, &segments, &branches, &tones, 1, 1e-12, budget, || {
        checks += 1;
        false
    }).unwrap();
    let mut late = 0;
    assert!(matches!(
        prepare_maxwell_harmonic_cycle(0.0, &segments, &branches, &tones, 1, 1e-12, budget, || {
            late += 1;
            late == checks
        }),
        Err(RotorError::Cancelled)
    ));
    let initial = RotorMaxwellState {
        displacement: vec![0.2], velocity: vec![-0.1], branch_stress: vec![0.6],
    };
    assert!(matches!(fourier.advance(&initial, 1e-10, || true), Err(RotorError::Cancelled)));
    assert_eq!(
        fourier.advance(&initial, 1e-10, || false).unwrap(),
        held.advance(&initial, 1e-10, || false).unwrap()
    );
}
