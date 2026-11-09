//! Original material clocks drive rotor memory, resonance and physical heat.
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_mbd::rotordynamics::periodic::{
    MAX_ROTOR_PERIODIC_SEGMENTS, RotorFloquetBudget, RotorHarmonicForce, RotorMaxwellState,
    RotorPeriodicForcing, RotorPeriodicSegment,
};
use fs_mbd::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use fs_mbd_material::rotor_kinetics::{RotorTemperatureError, RotorTemperatureKinetics};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use std::f64::consts::{LN_2, TAU};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn branch(strain: f64, k: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: vec![strain],
        stiffness: k,
        relaxation_time_s: tau,
    }
}

fn wlf() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 1.0,
            c2_k: 100.0,
        },
        300.0,
        [250.0, 650.0],
    )
    .unwrap()
}

fn arrhenius() -> TemperatureShift {
    // At 600 K, ln(a_T) = 600 ln(2) * (1/600 - 1/300) = -ln(2).
    TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 600.0 * LN_2 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        },
        300.0,
        [250.0, 650.0],
    )
    .unwrap()
}

#[test]
fn heterogeneous_calibrations_always_use_the_original_reference() {
    let branches = [branch(1.0, 2.0, 0.4), branch(-0.5, 3.0, 0.8)];
    let shifts = [wlf(), arrhenius()];
    let mut temperatures = vec![vec![300.0, 300.0], vec![400.0, 600.0], vec![300.0, 450.0]];
    let kinetics =
        RotorTemperatureKinetics::try_new(&branches, &shifts, &temperatures, || false).unwrap();
    let scales = kinetics.as_kinetics().relaxation_time_scales;
    assert_eq!(scales[0], [1.0, 1.0]);
    near(scales[1][0], 10.0_f64.powf(-0.5), 2e-15);
    near(scales[1][1], 0.5, 2e-15);
    assert_eq!(scales[2][0].to_bits(), 1.0_f64.to_bits());
    near(scales[2][1], 2.0_f64.powf(-2.0 / 3.0), 2e-15);
    assert!(std::ptr::eq(kinetics.branches(), branches.as_slice()));
    assert_eq!(kinetics.temperature_shifts(), shifts);
    assert_eq!(kinetics.temperatures_k(), temperatures);
    temperatures[0][0] = 500.0;
    assert_eq!(kinetics.temperatures_k()[0][0], 300.0);
}

#[test]
fn switched_material_clocks_preserve_memory_and_release_original_branch_storage() {
    // The mechanical mode is independently q=exp(-t), v=-exp(-t) for
    // M=2,C=3,K=1. Dark material channels cannot alter it or have C scaled.
    let model = RotorBearingModel::new(1, &[2.0], &[3.0], &[1.0], &[0.0], || false).unwrap();
    let branches = [branch(0.0, 2.0, 0.4), branch(0.0, 3.0, 0.8)];
    let shifts = [wlf(), arrhenius()];
    let temperatures = [vec![300.0, 300.0], vec![400.0, 600.0], vec![300.0, 450.0]];
    let durations = [0.15, 0.2, 0.25];
    let segments = durations.map(|duration_seconds| RotorPeriodicSegment {
        model: &model,
        duration_seconds,
    });
    let kinetics =
        RotorTemperatureKinetics::try_new(&branches, &shifts, &temperatures, || false).unwrap();
    let floquet = kinetics
        .floquet(0.0, &segments, 1e-11, RotorFloquetBudget::default(), || {
            false
        })
        .unwrap();
    let prepared = kinetics
        .prepare_power_cycle(
            0.0,
            &segments,
            RotorPeriodicForcing::Homogeneous,
            1e-11,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
    let initial = RotorMaxwellState {
        displacement: vec![1.0],
        velocity: vec![-1.0],
        branch_stress: vec![4.0, -3.0],
    };
    let (trajectory, power) = prepared.advance(&initial, 1e-10, 1e-10, || false).unwrap();
    // Independent closed-form factors, not values read back from the adapter.
    let ratios = [
        [1.0, 1.0],
        [10.0_f64.powf(-0.5), 0.5],
        [1.0, 2.0_f64.powf(-2.0 / 3.0)],
    ];
    let mut decay = [0.0_f64; 2];
    let mut time = 0.0;
    for (i, &dt) in durations.iter().enumerate() {
        time += dt;
        let state = &trajectory.boundary_states[i + 1];
        near(state.displacement[0], (-time).exp(), 2e-11);
        near(state.velocity[0], -(-time).exp(), 2e-11);
        for j in 0..2 {
            let old_energy = initial.branch_stress[j].powi(2) * (-2.0 * decay[j]).exp()
                / (2.0 * branches[j].stiffness);
            decay[j] += dt / (branches[j].relaxation_time_s * ratios[i][j]);
            let expected_stress = initial.branch_stress[j] * (-decay[j]).exp();
            let energy = expected_stress.powi(2) / (2.0 * branches[j].stiffness);
            near(state.branch_stress[j], expected_stress, 2e-11);
            near(
                power.segments[i].branch_heat_j[j],
                old_energy - energy,
                2e-11,
            );
            near(
                power.segments[i].branch_elastic_energy_change_j[j],
                energy - old_energy,
                2e-11,
            );
        }
        assert_eq!(power.segments[i].stiffness_switch_work_j, 0.0);
    }
    for j in 0..2 {
        let energy = initial.branch_stress[j].powi(2) / (2.0 * branches[j].stiffness);
        near(
            power.branch_heat_j[j],
            energy * (1.0 - (-2.0 * decay[j]).exp()),
            2e-11,
        );
        near(
            floquet.monodromy[(2 + j) * 4 + 2 + j],
            (-decay[j]).exp(),
            2e-11,
        );
    }
    near(
        power.velocity_force_work_j,
        1.5 * (1.0 - (-2.0 * time).exp()),
        2e-11,
    );
    near(power.balance_residual_j, 0.0, 1e-10);
    assert_eq!(power.stiffness_switch_work_j, 0.0);
    assert_eq!(power.applied_work_j, 0.0);
    assert_eq!(initial.branch_stress, [4.0, -3.0]);
}

#[test]
fn temperature_shifted_resonance_matches_analytic_maxwell_dispersion_and_heat() {
    let (m, c, k, branch_k, tau_reference, omega) = (1.2, 0.35, 2.0, 1.7, 0.4, 1.25);
    let model = RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap();
    let branches = [branch(1.0, branch_k, tau_reference)];
    let shifts = [arrhenius()];
    let period = TAU / omega;
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: period,
    }];
    let (cosine, sine) = (1.3, 0.2);
    let rows = [RotorHarmonicForce {
        cosine: &[cosine],
        sine: &[sine],
    }];
    let forcing = RotorPeriodicForcing::Harmonic {
        rows: &rows,
        order: 1,
    };
    for (temperature, ratio) in [(300.0, 1.0), (600.0, 0.5)] {
        let kinetics =
            RotorTemperatureKinetics::try_new(&branches, &shifts, &[vec![temperature]], || false)
                .unwrap();
        let motion = kinetics
            .prepare_cycle(
                0.0,
                &segments,
                forcing,
                1e-11,
                RotorFloquetBudget::default(),
                || false,
            )
            .unwrap();
        let orbit = motion.periodic_orbit(1e-10, 1e8, || false).unwrap();
        let power_cycle = kinetics
            .prepare_power_cycle(
                0.0,
                &segments,
                forcing,
                1e-11,
                RotorFloquetBudget::default(),
                || false,
            )
            .unwrap();
        let (_, power) = power_cycle
            .periodic_orbit(1e-10, 1e8, 1e-9, || false)
            .unwrap();

        // D = K - M*w² + i*C*w + k_b*(i*w*tau)/(1+i*w*tau).
        let tau = tau_reference * ratio;
        let x = omega * tau;
        let memory_re = branch_k * x * x / (1.0 + x * x);
        let memory_im = branch_k * x / (1.0 + x * x);
        let d_re = k - m * omega * omega + memory_re;
        let d_im = c * omega + memory_im;
        let denominator = d_re * d_re + d_im * d_im;
        // F = cosine - i*sine for exp(+i*w*t).
        let q_re = (cosine * d_re - sine * d_im) / denominator;
        let q_im = (-sine * d_re - cosine * d_im) / denominator;
        let stress_re = memory_re * q_re - memory_im * q_im;
        let stress_im = memory_re * q_im + memory_im * q_re;
        let initial = &orbit.trajectory.boundary_states[0];
        near(initial.displacement[0], q_re, 2e-10);
        near(initial.velocity[0], -omega * q_im, 2e-10);
        near(initial.branch_stress[0], stress_re, 2e-10);
        let expected_heat =
            period * (stress_re * stress_re + stress_im * stress_im) / (2.0 * branch_k * tau);
        let expected_viscous_work = period * c * omega * omega * (q_re * q_re + q_im * q_im) / 2.0;
        near(power.branch_heat_j[0], expected_heat, 2e-10);
        near(power.velocity_force_work_j, expected_viscous_work, 2e-10);
        near(
            power.applied_work_j,
            expected_heat + expected_viscous_work,
            2e-10,
        );
    }
    assert_eq!(branches[0].relaxation_time_s, tau_reference);
}

#[test]
fn unsupported_temperatures_invalid_clocks_and_cancellation_publish_no_schedule() {
    let branches = [branch(1.0, 2.0, 0.4)];
    let shifts = [wlf()];
    for bad in [249.0, 651.0, f64::NAN] {
        assert!(matches!(
            RotorTemperatureKinetics::try_new(
                &branches,
                &shifts,
                &[vec![300.0], vec![bad]],
                || false,
            ),
            Err(RotorTemperatureError::Shift {
                segment: 1,
                branch: 0,
                ..
            })
        ));
    }
    for rows in [
        vec![],
        vec![vec![]],
        vec![vec![300.0]; MAX_ROTOR_PERIODIC_SEGMENTS + 1],
    ] {
        assert!(matches!(
            RotorTemperatureKinetics::try_new(&branches, &shifts, &rows, || false,),
            Err(RotorTemperatureError::Invalid(_))
        ));
    }
    assert!(matches!(
        RotorTemperatureKinetics::try_new(&branches, &[], &[vec![300.0]], || false,),
        Err(RotorTemperatureError::Invalid(_))
    ));
    assert!(matches!(
        RotorTemperatureKinetics::try_new(
            &vec![branches[0].clone(); MAX_ROTOR_MAXWELL_BRANCHES + 1],
            &vec![shifts[0]; MAX_ROTOR_MAXWELL_BRANCHES + 1],
            &[vec![300.0; MAX_ROTOR_MAXWELL_BRANCHES + 1]],
            || false,
        ),
        Err(RotorTemperatureError::Invalid(_))
    ));
    for (time, temperature) in [(f64::from_bits(1), 300.0), (f64::MAX, 250.0)] {
        assert!(matches!(
            RotorTemperatureKinetics::try_new(
                &[branch(1.0, 2.0, time)],
                &shifts,
                &[vec![temperature]],
                || false,
            ),
            Err(RotorTemperatureError::Rotor(RotorError::Invalid(_)))
        ));
    }
    let rows = [vec![300.0], vec![400.0]];
    let mut polls = 0;
    RotorTemperatureKinetics::try_new(&branches, &shifts, &rows, || {
        polls += 1;
        false
    })
    .unwrap();
    for cancel_at in [1, polls / 2, polls] {
        let mut observed = 0;
        assert!(matches!(
            RotorTemperatureKinetics::try_new(&branches, &shifts, &rows, || {
                observed += 1;
                observed == cancel_at
            },),
            Err(RotorTemperatureError::Rotor(RotorError::Cancelled))
        ));
    }
    assert_eq!(branches[0].strain, [1.0]);
    assert_eq!(branches[0].stiffness, 2.0);
    assert_eq!(branches[0].relaxation_time_s, 0.4);
}
