use super::*;
use fs_material::fiber::felt_rate::{CreepBranch, MAX_FELT_RATE_BRANCHES, RateFelt};
use fs_material::fiber::{Uniaxial, WoolFelt};

fn material(crush: bool, count: usize) -> RateFelt {
    RateFelt::new(
        WoolFelt::new(
            30000.0,
            0.2,
            2.0,
            if crush { 3.0 } else { 2.0 },
            if crush { 0.15 } else { 0.0 },
            0.8,
        )
        .unwrap(),
        (0..count)
            .map(|i| CreepBranch {
                modulus_pa: 60000.0 * (i + 1) as f64,
                viscosity_pa_s: 1200.0 * (i + 1) as f64,
            })
            .collect(),
    )
    .unwrap()
}
fn pad(law: &RateFelt, compression: f64) -> FeltPad {
    FeltPad::from_material(0.001, 0.006, compression, vec![1.0], 0.0, law).unwrap()
}

#[test]
fn g1_all_material_arms_recover_with_the_original_exact_decay_and_positive_heat() {
    let gate = CancelGate::new_clock_free();
    let law = material(false, MAX_FELT_RATE_BRANCHES);
    let pad = pad(&law, -0.01); // original contact is separated throughout
    let clock = config(0.02);
    // This rate lies outside the old midpoint phase/decay gate.
    assert!(matches!(
        ImpactSystem::new(
            vec![linear(1.0, 0.0)],
            vec![],
            vec![pad.clone()],
            vec![],
            clock
        ),
        Err(ImpactError::TimeResolution(_))
    ));
    let mut run = ExponentialImpactSystem::new_with_felt(
        vec![linear(1.0, 0.0)],
        vec![],
        vec![pad.clone()],
        vec![],
        vec![],
        clock,
        &gate,
    )
    .unwrap();
    assert_eq!(run.state().len(), 2 + MAX_FELT_RATE_BRANCHES);
    let initial: Vec<_> = pad
        .creep
        .iter()
        .enumerate()
        .map(|(i, branch)| branch.stiffness_n_m.sqrt() * 1e-5 * (i + 1) as f64)
        .collect();
    run.inner.x[2..].copy_from_slice(&initial);
    let energy = run.stored_energy_j();
    let mut heat = 0.0;
    for tick in 1..=12 {
        let frame = run.step(&[0.0], &gate).unwrap();
        heat += frame.dissipated_energy_j;
        assert_eq!(frame.felt_crush_loss_j, 0.0);
        assert_eq!(&run.state()[..2], &[0.0, 0.0]);
        for (i, branch) in pad.creep.iter().enumerate() {
            let expected = initial[i]
                * (-(branch.stiffness_n_m / branch.viscosity_n_s_m) * tick as f64 * clock.dt_s)
                    .exp();
            assert!((run.state()[2 + i] - expected).abs() < 3e-10 * initial[i]);
        }
        assert!(frame.dissipated_energy_j >= 0.0);
    }
    let expected_heat = initial
        .iter()
        .map(|z| 0.5 * z * z * (-(-2.0 * 50.0 * 12.0 * clock.dt_s).exp_m1()))
        .sum::<f64>();
    assert!((heat - expected_heat).abs() < 3e-10 * energy);
    assert!((run.stored_energy_j() + heat - energy).abs() < 3e-10 * energy);
}

fn rhs(x: [f64; 4], omega: f64, p: &FeltPad) -> [f64; 4] {
    let compression = p.precompression_m + x[0] - x[2] - x[3];
    let strain = (compression / p.thickness_m).max(0.0);
    // Independent ORIGINAL conservative power law and Kelvin equations.
    let force = p.area_m2 * p.law.f_ref * (strain / p.law.eps_ref).powi(2);
    [
        x[1],
        -omega * omega * x[0] - force,
        (force - p.creep[0].stiffness_n_m * x[2]) / p.creep[0].viscosity_n_s_m,
        (force - p.creep[1].stiffness_n_m * x[3]) / p.creep[1].viscosity_n_s_m,
    ]
}
fn reference(horizon: f64, omega: f64, p: &FeltPad) -> [f64; 4] {
    let mut x = [0.0; 4];
    let dt = horizon / 16384.0;
    for _ in 0..16384 {
        let k1 = rhs(x, omega, p);
        let k2 = rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]), omega, p);
        let k3 = rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]), omega, p);
        let k4 = rhs(core::array::from_fn(|i| x[i] + dt * k3[i]), omega, p);
        for i in 0..4 {
            x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    x
}

#[test]
fn g1_nonlinear_material_contact_keeps_work_heat_and_temporal_refinement() {
    let gate = CancelGate::new_clock_free();
    let pad = pad(&material(false, 2), 0.0012);
    let (omega, horizon) = (200.0, 0.003);
    let expected = reference(horizon, omega, &pad);
    let mut errors = Vec::new();
    for ticks in [32, 64, 128] {
        let mut run = ExponentialImpactSystem::new_with_felt(
            vec![linear(omega, 0.0)],
            vec![],
            vec![pad.clone()],
            vec![],
            vec![],
            config(horizon / ticks as f64),
            &gate,
        )
        .unwrap();
        let energy = run.stored_energy_j();
        let mut heat = 0.0;
        for _ in 0..ticks {
            let frame = run.step(&[0.0], &gate).unwrap();
            assert!(frame.dissipated_energy_j >= -1e-15);
            heat += frame.dissipated_energy_j;
        }
        assert!(heat > 0.0);
        assert!((run.stored_energy_j() + heat - energy).abs() < 2e-9 * energy);
        let actual = [
            run.state()[0],
            run.state()[1],
            run.state()[2] / pad.creep[0].stiffness_n_m.sqrt(),
            run.state()[3] / pad.creep[1].stiffness_n_m.sqrt(),
        ];
        let error = ((actual[0] - expected[0]) * omega).powi(2)
            + (actual[1] - expected[1]).powi(2)
            + ((actual[2] - expected[2]) * omega).powi(2)
            + ((actual[3] - expected[3]) * omega).powi(2);
        errors.push(error.sqrt());
    }
    assert!(
        errors[0] / errors[1] > 3.2 && errors[1] / errors[2] > 3.2,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-5, "{errors:?}");
}

#[test]
fn g4_conditioning_and_creep_commit_together_and_failed_trials_retry_identically() {
    let gate = CancelGate::new_clock_free();
    let make = || {
        let (body, weight) = ImpactBody::free_mass(1.0, 0.0, 0.3).unwrap();
        let mut p = pad(&material(true, 2), 0.0001);
        p.weights[0] = weight;
        ExponentialImpactSystem::new_with_felt(
            vec![body],
            vec![],
            vec![p],
            vec![],
            vec![],
            config(1e-5),
            &gate,
        )
        .unwrap()
    };
    let (mut run, mut fresh) = (make(), make());
    let initial = run.stored_energy_j();
    let mut heat = 0.0;
    let mut crush = 0.0;
    for _ in 0..100 {
        let frame = run.step(&[0.0], &gate).unwrap();
        fresh.step(&[0.0], &gate).unwrap();
        heat += frame.dissipated_energy_j;
        crush += frame.felt_crush_loss_j;
    }
    assert!(crush > 0.0 && heat >= crush);
    assert!((run.stored_energy_j() + heat - initial).abs() < 2e-9 * initial);
    let before = run.state().to_vec();
    let history = run.felt_history(0).unwrap();
    let samples = run.samples();
    let energy = run.stored_energy_j();
    run.set_iteration_limit(0).unwrap();
    assert!(run.step(&[0.0], &gate).is_err());
    assert_eq!(run.state(), before);
    assert_eq!(run.felt_history(0).unwrap(), history);
    assert_eq!(run.samples(), samples);
    assert_eq!(run.stored_energy_j(), energy);
    let stop = CancelGate::new_clock_free();
    stop.request();
    assert!(matches!(
        run.step(&[0.0], &stop),
        Err(ImpactError::Cancelled)
    ));
    run.set_iteration_limit(50).unwrap();
    let actual = run.step(&[0.0], &gate).unwrap();
    let expected = fresh.step(&[0.0], &gate).unwrap();
    assert_eq!(run.state(), fresh.state());
    assert_eq!(run.felt_history(0), fresh.felt_history(0));
    assert_eq!(actual.stored_energy_j, expected.stored_energy_j);
    assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
}

#[test]
fn g0_complete_kelvin_state_counts_toward_native_capacity_and_densification_refuses() {
    let gate = CancelGate::new_clock_free();
    let law = material(false, MAX_FELT_RATE_BRANCHES);
    let modes = 48;
    assert!(PreparedNonlinearZoh::fits_reference_capacity(
        2 * modes,
        modes
    ));
    assert!(!PreparedNonlinearZoh::fits_reference_capacity(
        2 * modes + 8,
        modes
    ));
    let mut pad = pad(&law, 0.0);
    pad.weights.resize(modes, 0.0);
    assert!(
        ExponentialImpactSystem::new_with_felt(
            vec![linear(1.0, 0.0); modes],
            vec![],
            vec![pad],
            vec![],
            vec![],
            config(1e-5),
            &gate,
        )
        .is_err()
    );
    let mut p = FeltPad::from_material(0.001, 0.006, 0.01, vec![1.0], 0.0, &law).unwrap();
    assert!(
        ExponentialImpactSystem::new_with_felt(
            vec![linear(1.0, 0.0)],
            vec![],
            vec![p.clone()],
            vec![],
            vec![],
            config(1e-5),
            &gate,
        )
        .is_err()
    );
    p.precompression_m = 0.0;
    p.creep.push(p.creep[0]);
    assert!(
        ExponentialImpactSystem::new_with_felt(
            vec![linear(1.0, 0.0)],
            vec![],
            vec![p],
            vec![],
            vec![],
            config(1e-5),
            &gate,
        )
        .is_err()
    );
}
