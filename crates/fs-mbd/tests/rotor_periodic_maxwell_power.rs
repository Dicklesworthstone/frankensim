//! G1/G3 original-memory work, heat, switched storage and source interference.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::periodic::power::{RotorCyclePowerError, prepare_forced_power_cycle};
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce, RotorMaxwellState,
    RotorPeriodicSegment, prepare_maxwell_forced_power_cycle, prepare_maxwell_fourier_power_cycle,
    prepare_maxwell_harmonic_power_cycle,
};
use fs_mbd::rotordynamics::relaxation::MaxwellDampingBranch;
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::f64::consts::TAU;

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

fn state(q: f64, v: f64, stress: &[f64]) -> RotorMaxwellState {
    RotorMaxwellState {
        displacement: vec![q],
        velocity: vec![v],
        branch_stress: stress.to_vec(),
    }
}

fn budget() -> RotorFloquetBudget {
    RotorFloquetBudget {
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 500_000_000,
        ..RotorFloquetBudget::default()
    }
}

#[test]
fn dark_memory_releases_its_own_original_storage_as_individual_heat() {
    let model = oscillator(1.0, 0.3, 2.0);
    let branches = [arm(&[0.0], 2.0, 0.25), arm(&[0.0], 3.0, 0.5)];
    let segments = [
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.15,
        },
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.25,
        },
    ];
    let prepared = prepare_maxwell_forced_power_cycle(
        4.0,
        &segments,
        &branches,
        &[&[0.0], &[0.0]],
        1e-12,
        budget(),
        || false,
    )
    .unwrap();
    let initial = state(0.0, 0.0, &[4.0, -3.0]);
    let (trajectory, power) = prepared.advance(&initial, 1e-10, 1e-11, || false).unwrap();
    for (j, branch) in branches.iter().enumerate() {
        let energy = initial.branch_stress[j].powi(2) / (2.0 * branch.stiffness);
        let expected = energy * (1.0 - (-0.8 / branch.relaxation_time_s).exp());
        near(power.branch_heat_j[j], expected, 3e-12);
        near(power.branch_elastic_energy_change_j[j], -expected, 3e-12);
        near(power.average_branch_heat_power_w[j], expected / 0.4, 3e-12);
        for (&time, current) in trajectory
            .boundary_times_seconds
            .iter()
            .zip(&trajectory.boundary_states)
        {
            near(
                current.branch_stress[j],
                initial.branch_stress[j] * (-time / branch.relaxation_time_s).exp(),
                3e-12,
            );
        }
    }
    assert_eq!(power.applied_work_j, 0.0);
    assert_eq!(power.velocity_force_work_j, 0.0);
    assert_eq!(power.kinetic_energy_change_j, 0.0);
    assert_eq!(power.relaxed_elastic_energy_change_j, 0.0);
    near(power.balance_residual_j, 0.0, 1e-11);
    assert_eq!(initial, state(0.0, 0.0, &[4.0, -3.0]));
}

#[test]
fn independent_decay_mode_resolves_viscous_and_memory_dissipation_separately() {
    // Exact physical mode of M=1,C=2,K=1.5,k=1.5,tau=.25:
    // [q,v,sigma]=[1,-1,-.5]*exp(-t). Initial energy is4/3 J;
    // integrated C work=(1-e^-2T), branch heat=(1-e^-2T)/3.
    let model = oscillator(1.0, 2.0, 1.5);
    let branches = [arm(&[1.0], 1.5, 0.25)];
    let segments = [
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.3,
        },
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.4,
        },
    ];
    let prepared = prepare_maxwell_forced_power_cycle(
        0.0,
        &segments,
        &branches,
        &[&[0.0], &[0.0]],
        1e-12,
        budget(),
        || false,
    )
    .unwrap();
    let (trajectory, power) = prepared
        .advance(&state(1.0, -1.0, &[-0.5]), 1e-10, 1e-11, || false)
        .unwrap();
    let released = -(-1.4_f64).exp_m1();
    near(power.velocity_force_work_j, released, 3e-12);
    near(power.branch_heat_j[0], released / 3.0, 3e-12);
    near(power.kinetic_energy_change_j, -0.5 * released, 3e-12);
    near(
        power.relaxed_elastic_energy_change_j,
        -0.75 * released,
        3e-12,
    );
    near(
        power.branch_elastic_energy_change_j[0],
        -released / 12.0,
        3e-12,
    );
    near(power.balance_residual_j, 0.0, 1e-11);
    for (&time, current) in trajectory
        .boundary_times_seconds
        .iter()
        .zip(&trajectory.boundary_states)
    {
        near(current.displacement[0], (-time).exp(), 3e-12);
        near(current.velocity[0], -(-time).exp(), 3e-12);
        near(current.branch_stress[0], -0.5 * (-time).exp(), 3e-12);
    }
}

#[test]
fn resonant_periodic_consumption_matches_original_phasor_heat_and_source_units() {
    let (m, c, omega) = (1.3_f64, 0.04_f64, 1.7_f64);
    let branches = [
        arm(&[1.0], 0.45, 0.7),
        arm(&[-0.6], 0.3, 1.1),
        arm(&[0.0], 0.8, 1.3),
    ];
    let transfer: Vec<_> = branches
        .iter()
        .map(|branch| {
            C64::new(0.0, omega * branch.relaxation_time_s)
                / C64::new(1.0, omega * branch.relaxation_time_s)
                * C64::from_re(branch.stiffness)
        })
        .collect();
    let k = m * omega * omega
        - branches
            .iter()
            .zip(&transfer)
            .map(|(branch, h)| h.re * branch.strain[0].powi(2))
            .sum::<f64>();
    let model = oscillator(m, c, k);
    let period = TAU / omega;
    let segments = [
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.23 * period,
        },
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.77 * period,
        },
    ];
    let force = C64::new(0.6, -0.8);
    let dynamic = transfer.iter().zip(&branches).fold(
        C64::new(k - m * omega * omega, c * omega),
        |d, (h, branch)| d + h.scale(branch.strain[0].powi(2)),
    );
    let q = force / dynamic;
    let velocity = C64::new(0.0, omega) * q;
    let heat: Vec<_> = branches
        .iter()
        .zip(&transfer)
        .map(|(branch, h)| {
            let stress = (*h * q).scale(branch.strain[0]);
            stress.norm_sq() / (2.0 * branch.stiffness * branch.relaxation_time_s)
        })
        .collect();
    let viscous = 0.5 * c * velocity.norm_sq();
    let supply = 0.5 * (force * velocity.conj()).re;
    near(supply, viscous + heat.iter().sum::<f64>(), 1e-13);
    for scale in [1e-60, 1.0, 1e60] {
        let cosine = [0.6 * scale];
        let sine = [0.8 * scale];
        let forces = [RotorHarmonicForce {
            cosine: &cosine,
            sine: &sine,
        }; 2];
        let prepared = prepare_maxwell_harmonic_power_cycle(
            -2.7,
            &segments,
            &branches,
            &forces,
            1,
            1e-12,
            budget(),
            || false,
        )
        .unwrap();
        let (orbit, power) = prepared
            .periodic_orbit(1e-9, 1e8, 2e-9 * scale * scale, || false)
            .unwrap();
        near(power.average_applied_power_w / scale / scale, supply, 3e-9);
        near(
            power.average_velocity_force_power_w / scale / scale,
            viscous,
            3e-9,
        );
        for (j, expected) in heat.iter().enumerate() {
            near(
                power.average_branch_heat_power_w[j] / scale / scale,
                *expected,
                3e-9,
            );
            near(
                power.branch_elastic_energy_change_j[j] / scale / scale,
                0.0,
                3e-9,
            );
        }
        near(power.balance_residual_j / scale / scale, 0.0, 2e-9);
        assert!(orbit.relative_fixed_point_residual < 1e-9);
        assert!(orbit.relative_closure_residual < 1e-9);
        assert_eq!(
            prepared.motion_cycle().floquet_report().multipliers.len(),
            5
        );
    }
}

// Independent original scalar ODE and separately integrated physical work.
// This reference never calls a production exponential, Gramian or matrix solve.
fn rhs(
    t: f64,
    x: &[f64; 8],
    c: f64,
    k: f64,
    force: [f64; 5],
    branches: &[MaxwellDampingBranch],
) -> [f64; 8] {
    let f = force[0]
        + force[1] * (TAU * t).cos()
        + force[2] * (TAU * t).sin()
        + force[3] * (3.0 * TAU * t).cos()
        + force[4] * (3.0 * TAU * t).sin();
    let memory_force = branches[0].strain[0] * x[2] + branches[1].strain[0] * x[3];
    [
        x[1],
        (f - c * x[1] - k * x[0] - memory_force) / 1.2,
        branches[0].stiffness * branches[0].strain[0] * x[1] - x[2] / branches[0].relaxation_time_s,
        branches[1].stiffness * branches[1].strain[0] * x[1] - x[3] / branches[1].relaxation_time_s,
        f * x[1],
        c * x[1] * x[1],
        x[2] * x[2] / (branches[0].stiffness * branches[0].relaxation_time_s),
        x[3] * x[3] / (branches[1].stiffness * branches[1].relaxation_time_s),
    ]
}

fn reference_segment(
    mut x: [f64; 8],
    start: f64,
    duration: f64,
    c: f64,
    k: f64,
    force: [f64; 5],
    branches: &[MaxwellDampingBranch],
) -> [f64; 8] {
    let step = duration / 4096.0;
    for index in 0..4096 {
        let t = start + index as f64 * step;
        let a = rhs(t, &x, c, k, force, branches);
        let b = rhs(
            t + 0.5 * step,
            &std::array::from_fn(|i| x[i] + 0.5 * step * a[i]),
            c,
            k,
            force,
            branches,
        );
        let c1 = rhs(
            t + 0.5 * step,
            &std::array::from_fn(|i| x[i] + 0.5 * step * b[i]),
            c,
            k,
            force,
            branches,
        );
        let d = rhs(
            t + step,
            &std::array::from_fn(|i| x[i] + step * c1[i]),
            c,
            k,
            force,
            branches,
        );
        for i in 0..8 {
            x[i] += step * (a[i] + 2.0 * b[i] + 2.0 * c1[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn switched_fourier_memory_heat_keeps_all_cross_terms_and_modulator_work() {
    let models = [oscillator(1.2, 0.1, 2.0), oscillator(1.2, -0.03, 3.4)];
    let branches = [arm(&[1.0], 0.7, 0.3), arm(&[-0.6], 0.35, 0.8)];
    let segments = [
        RotorPeriodicSegment {
            model: &models[0],
            duration_seconds: 0.37,
        },
        RotorPeriodicSegment {
            model: &models[1],
            duration_seconds: 0.63,
        },
    ];
    let law = [[0.3, 0.4, 0.2, 0.1, -0.1], [-0.2, -0.3, 0.5, 0.2, 0.07]];
    let tones = [
        [
            RotorHarmonicForce {
                cosine: &[0.4],
                sine: &[0.2],
            },
            RotorHarmonicForce {
                cosine: &[0.1],
                sine: &[-0.1],
            },
        ],
        [
            RotorHarmonicForce {
                cosine: &[-0.3],
                sine: &[0.5],
            },
            RotorHarmonicForce {
                cosine: &[0.2],
                sine: &[0.07],
            },
        ],
    ];
    let forces = [
        RotorFourierForce {
            steady: &[0.3],
            harmonics: &tones[0],
        },
        RotorFourierForce {
            steady: &[-0.2],
            harmonics: &tones[1],
        },
    ];
    let prepared = prepare_maxwell_fourier_power_cycle(
        0.0,
        &segments,
        &branches,
        &forces,
        &[1, 3],
        1e-12,
        budget(),
        || false,
    )
    .unwrap();
    let initial = state(0.2, -0.3, &[0.4, -0.2]);
    let (trajectory, power) = prepared.advance(&initial, 1e-9, 1e-9, || false).unwrap();
    let first = reference_segment(
        [0.2, -0.3, 0.4, -0.2, 0.0, 0.0, 0.0, 0.0],
        0.0,
        0.37,
        0.1,
        2.0,
        law[0],
        &branches,
    );
    let last = reference_segment(first, 0.37, 0.63, -0.03, 3.4, law[1], &branches);
    for (current, expected) in trajectory.boundary_states[1..].iter().zip([first, last]) {
        near(current.displacement[0], expected[0], 2e-9);
        near(current.velocity[0], expected[1], 2e-9);
        near(current.branch_stress[0], expected[2], 2e-9);
        near(current.branch_stress[1], expected[3], 2e-9);
    }
    near(power.applied_work_j, last[4], 2e-9);
    near(power.velocity_force_work_j, last[5], 2e-9);
    near(power.branch_heat_j[0], last[6], 2e-9);
    near(power.branch_heat_j[1], last[7], 2e-9);
    near(power.segments[0].applied_work_j, first[4], 2e-9);
    near(power.segments[1].applied_work_j, last[4] - first[4], 2e-9);
    assert!(power.segments[1].velocity_force_work_j < 0.0);
    let switch1 = 0.5 * (3.4 - 2.0) * first[0].powi(2);
    let seam = 0.5 * (2.0 - 3.4) * last[0].powi(2);
    near(power.segments[0].stiffness_switch_work_j, switch1, 2e-9);
    near(power.segments[1].stiffness_switch_work_j, seam, 2e-9);
    near(power.stiffness_switch_work_j, switch1 + seam, 2e-9);
    near(power.balance_residual_j, 0.0, 1e-9);
    assert!(power.absolute_segment_residual_j < 1e-9);
}

#[test]
fn full_mass_gyro_and_circulatory_work_are_conjugate_coordinate_invariant() {
    let mass = [2.0, 0.4, 0.4, 1.5];
    let damping = [0.12, -0.02, 0.07, 0.18];
    let gyro = [0.0, 0.8, -0.8, 0.0];
    let stiffness = [[4.0, 0.7, -0.2, 3.0], [5.0, -0.3, 0.6, 2.5]];
    let models: [RotorBearingModel; 2] = std::array::from_fn(|i| {
        RotorBearingModel::new(2, &mass, &damping, &stiffness[i], &gyro, || false).unwrap()
    });
    let segments = [
        RotorPeriodicSegment {
            model: &models[0],
            duration_seconds: 0.17,
        },
        RotorPeriodicSegment {
            model: &models[1],
            duration_seconds: 0.29,
        },
    ];
    let branches = [arm(&[1.2, -0.4], 1.7, 0.7), arm(&[-0.3, 0.8], 2.3, 0.2)];
    let forces = [[0.7, -0.4], [-0.1, 0.9]];
    let initial = RotorMaxwellState {
        displacement: vec![0.2, -0.1],
        velocity: vec![-0.3, 0.4],
        branch_stress: vec![0.3, -0.2],
    };
    let prepared = prepare_maxwell_forced_power_cycle(
        -2.3,
        &segments,
        &branches,
        &[&forces[0], &forces[1]],
        1e-12,
        budget(),
        || false,
    )
    .unwrap();
    let (_, base) = prepared.advance(&initial, 1e-9, 1e-9, || false).unwrap();
    let units = [1e3, 1e-2];
    let congruence =
        |a: &[f64; 4]| std::array::from_fn::<_, 4, _>(|i| a[i] * units[i / 2] * units[i % 2]);
    let changed: [RotorBearingModel; 2] = std::array::from_fn(|i| {
        RotorBearingModel::new(
            2,
            &congruence(&mass),
            &congruence(&damping),
            &congruence(&stiffness[i]),
            &congruence(&gyro),
            || false,
        )
        .unwrap()
    });
    let new_segments = [
        RotorPeriodicSegment {
            model: &changed[0],
            duration_seconds: 0.17,
        },
        RotorPeriodicSegment {
            model: &changed[1],
            duration_seconds: 0.29,
        },
    ];
    let new_branches: Vec<_> = branches
        .iter()
        .map(|branch| {
            arm(
                &[branch.strain[0] * units[0], branch.strain[1] * units[1]],
                branch.stiffness,
                branch.relaxation_time_s,
            )
        })
        .collect();
    let new_forces: [[f64; 2]; 2] =
        std::array::from_fn(|i| std::array::from_fn(|j| forces[i][j] * units[j]));
    let new_initial = RotorMaxwellState {
        displacement: initial
            .displacement
            .iter()
            .zip(units)
            .map(|(x, scale)| x / scale)
            .collect(),
        velocity: initial
            .velocity
            .iter()
            .zip(units)
            .map(|(x, scale)| x / scale)
            .collect(),
        branch_stress: initial.branch_stress.clone(),
    };
    let mapped = prepare_maxwell_forced_power_cycle(
        -2.3,
        &new_segments,
        &new_branches,
        &[&new_forces[0], &new_forces[1]],
        1e-12,
        budget(),
        || false,
    )
    .unwrap();
    let (_, equivalent) = mapped.advance(&new_initial, 1e-9, 1e-9, || false).unwrap();
    for (a, b) in [
        base.applied_work_j,
        base.velocity_force_work_j,
        base.nonconservative_stiffness_work_j,
        base.stiffness_switch_work_j,
        base.kinetic_energy_change_j,
        base.relaxed_elastic_energy_change_j,
    ]
    .into_iter()
    .zip([
        equivalent.applied_work_j,
        equivalent.velocity_force_work_j,
        equivalent.nonconservative_stiffness_work_j,
        equivalent.stiffness_switch_work_j,
        equivalent.kinetic_energy_change_j,
        equivalent.relaxed_elastic_energy_change_j,
    ]) {
        near(a, b, 3e-10);
    }
    for (a, b) in base.branch_heat_j.iter().zip(&equivalent.branch_heat_j) {
        near(*a, *b, 3e-10);
    }
    assert!(base.nonconservative_stiffness_work_j.abs() > 1e-5);
}

#[test]
fn zero_branch_owner_equivalence_and_transactional_refusals() {
    let model = oscillator(1.3, 0.2, 2.4);
    let segments = [
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.2,
        },
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.3,
        },
    ];
    let forces: [&[f64]; 2] = [&[0.7], &[-0.2]];
    let initial = state(0.1, -0.4, &[]);
    let prepared =
        prepare_maxwell_forced_power_cycle(0.0, &segments, &[], &forces, 1e-12, budget(), || false)
            .unwrap();
    let expected = prepare_forced_power_cycle(0.0, &segments, &forces, 1e-12, budget(), || false)
        .unwrap()
        .advance(
            &RotorTransientState {
                displacement: initial.displacement.clone(),
                velocity: initial.velocity.clone(),
            },
            1e-9,
            1e-10,
            || false,
        )
        .unwrap()
        .1;
    let complete = prepared.advance(&initial, 1e-9, 1e-10, || false).unwrap();
    near(complete.1.applied_work_j, expected.applied_work_j, 1e-12);
    near(
        complete.1.velocity_force_work_j,
        expected.velocity_force_work_j,
        1e-12,
    );
    assert!(complete.1.branch_heat_j.is_empty());
    assert_eq!(
        complete,
        prepared.advance(&initial, 1e-9, 1e-10, || false).unwrap()
    );
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            prepared
                .advance(&initial, 1e-9, tolerance, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    prepared
        .advance(&initial, 1e-9, 1e-10, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            prepared.advance(&initial, 1e-9, 1e-10, || {
                count += 1;
                count == stop
            }),
            Err(RotorCyclePowerError::Dynamics(RotorError::Cancelled))
        ));
    }
    let tiny = RotorFloquetBudget {
        max_dense_bytes: 1,
        ..budget()
    };
    assert!(
        prepare_maxwell_forced_power_cycle(0.0, &segments, &[], &forces, 1e-12, tiny, || false)
            .is_err()
    );
    assert!(
        prepare_maxwell_forced_power_cycle(0.0, &[], &[], &[], 1e-12, budget(), || false).is_err()
    );
    assert_eq!(initial, state(0.1, -0.4, &[]));
}
