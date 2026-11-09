//! Original material and acoustic memory inside the complete exponential flow.
use super::super::super::{felt::FeltPad, radiation, relaxation::InitialMemory};
use super::*;
use fs_material::fiber::WoolFelt;
use fs_material::fiber::felt_rate::{CreepBranch, RateFelt};
use fs_phs::RelaxationBranch;

fn arm(dimension: usize, stiffness: f64, tau: f64) -> RelaxationBranch {
    let mut projection = vec![0.0; dimension];
    projection[0] = 1.0;
    RelaxationBranch {
        projection,
        stiffness,
        relaxation_time_s: tau,
    }
}

#[test]
fn g1_original_maxwell_eigenmode_has_exact_decay_and_independent_heat() {
    let gate = CancelGate::new_clock_free();
    let (amplitude, rate, tau, stiffness) = (0.01_f64, 2.0_f64, 0.1_f64, 116.0_f64);
    let mut run = ExponentialImpactSystem::new(
        vec![linear(5.0, amplitude)],
        vec![],
        vec![],
        config(0.2),
        &gate,
    )
    .unwrap()
    .with_relaxation_branches(
        vec![arm(2, stiffness, tau)],
        InitialMemory::Unrelaxed,
        1,
        &gate,
    )
    .unwrap();
    // An independently manufactured eigenmode: lambda=-2 solves
    // lambda²+25+116*lambda*tau/(1+lambda*tau)=0.
    // This is an explicit physical initial viscous strain, not a solver seed.
    run.inner.x[1] = -rate * amplitude;
    run.inner.x[2] = stiffness.sqrt() * amplitude / (1.0 - rate * tau);
    let initial = run.stored_energy_j();
    let mut heat = 0.0;
    for tick in 1..=12 {
        let frame = run.step(&[0.0], &gate).unwrap();
        let time = tick as f64 * 0.2;
        let q = amplitude * (-rate * time).exp();
        let expected = [q, -rate * q, stiffness.sqrt() * q / (1.0 - rate * tau)];
        for (&a, &b) in run.state().iter().zip(&expected) {
            assert!((a - b).abs() < 2e-10 * amplitude, "{tick}: {a:e} vs {b:e}");
        }
        heat += frame.dissipated_energy_j;
        assert!(frame.dissipated_energy_j >= 0.0);
        assert_eq!(frame.felt_crush_loss_j, 0.0);
        assert_eq!(run.relaxation_branch_count(), 1);
        let observation = run.relaxation_observation().unwrap();
        let effort = stiffness.sqrt() * q - expected[2];
        assert!((observation.stored_energy_j - 0.5 * effort * effort).abs() < 2e-10 * initial);
        assert!(
            (observation.dissipated_power_w - effort * effort / tau).abs() < 2e-10 * initial / tau
        );
        let expected_heat = stiffness * amplitude * amplitude / 16.0 / tau
            * (-(-2.0 * rate * time).exp_m1())
            / (2.0 * rate);
        assert!((heat - expected_heat).abs() < 2e-10 * initial);
        assert!((run.stored_energy_j() + heat - initial).abs() < 2e-10 * initial);
    }
}

fn cross_arms(dimension: usize) -> Vec<RelaxationBranch> {
    [(120.0, 0.005, [1.0, -0.4]), (80.0, 0.06, [-0.2, 0.7])]
        .into_iter()
        .map(|(stiffness, relaxation_time_s, row)| {
            let mut projection = vec![0.0; dimension];
            projection[0] = row[0];
            projection[2] = row[1];
            RelaxationBranch {
                projection,
                stiffness,
                relaxation_time_s,
            }
        })
        .collect()
}
fn acoustic() -> radiation::Model {
    radiation::Model {
        ports: 2,
        poles: vec![radiation::Pole {
            omega: 60.0,
            zeta: 0.15,
            coupling: vec![1.7, -2.4],
        }],
    }
}
fn cross_model(acoustic_first: bool) -> ExponentialImpactSystem {
    let gate = CancelGate::new_clock_free();
    let mut bodies = vec![linear(20.0, 0.002), linear(37.0, -0.001)];
    bodies[0].initial[0].velocity_m_sqrt_kg_per_s = 0.015;
    bodies[1].initial[0].velocity_m_sqrt_kg_per_s = -0.007;
    bodies[0].damping_per_s[0] = 2.0;
    bodies[1].damping_per_s[0] = 3.0;
    let run = ExponentialImpactSystem::new(bodies, vec![], vec![], config(0.01), &gate).unwrap();
    if acoustic_first {
        run.with_radiation_load(&acoustic(), &[1, 0], 1, &gate)
            .unwrap()
            .with_relaxation_branches(cross_arms(6), InitialMemory::Relaxed, 2, &gate)
            .unwrap()
    } else {
        run.with_relaxation_branches(cross_arms(4), InitialMemory::Relaxed, 2, &gate)
            .unwrap()
            .with_radiation_load(&acoustic(), &[1, 0], 1, &gate)
            .unwrap()
    }
}

// Independent original equations in canonical physical order:
// [q0,p0,q1,p1,z0,z1,acoustic_q,acoustic_p,heat,external_work].
fn cross_rhs(x: [f64; 10]) -> [f64; 10] {
    let rows = [
        [120.0_f64.sqrt(), -0.4 * 120.0_f64.sqrt()],
        [-0.2 * 80.0_f64.sqrt(), 0.7 * 80.0_f64.sqrt()],
    ];
    let tau = [0.005, 0.06];
    let effort =
        core::array::from_fn::<_, 2, _>(|i| rows[i][0] * x[0] + rows[i][1] * x[2] - x[4 + i]);
    let c = [-2.4, 1.7];
    [
        x[1],
        -400.0 * x[0] - 2.0 * x[1] - rows[0][0] * effort[0] - rows[1][0] * effort[1] - c[0] * x[7]
            + 0.1,
        x[3],
        -1369.0 * x[2]
            - 3.0 * x[3]
            - rows[0][1] * effort[0]
            - rows[1][1] * effort[1]
            - c[1] * x[7]
            - 0.05,
        effort[0] / tau[0],
        effort[1] / tau[1],
        x[7],
        -3600.0 * x[6] - 18.0 * x[7] + c[0] * x[1] + c[1] * x[3],
        2.0 * x[1] * x[1]
            + 3.0 * x[3] * x[3]
            + effort[0] * effort[0] / tau[0]
            + effort[1] * effort[1] / tau[1]
            + 18.0 * x[7] * x[7],
        0.1 * x[1] - 0.05 * x[3],
    ]
}
fn cross_reference() -> [f64; 10] {
    let mut x = [0.002, 0.015, -0.001, -0.007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    x[4] = 120.0_f64.sqrt() * (x[0] - 0.4 * x[2]);
    x[5] = 80.0_f64.sqrt() * (-0.2 * x[0] + 0.7 * x[2]);
    let dt = 0.12 / 32768.0;
    for _ in 0..32768 {
        let a = cross_rhs(x);
        let b = cross_rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = cross_rhs(core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = cross_rhs(core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..10 {
            x[i] += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    x
}

#[test]
fn g1_cross_mode_material_and_acoustic_reaction_share_one_continuous_flow_in_either_order() {
    let gate = CancelGate::new_clock_free();
    let (mut first, mut second) = (cross_model(false), cross_model(true));
    assert_eq!(first.state().len(), 8);
    assert_eq!(second.state().len(), 8);
    let initial = first.stored_energy_j();
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..12 {
        let a = first.step(&[0.1, -0.05], &gate).unwrap();
        let b = second.step(&[0.1, -0.05], &gate).unwrap();
        heat += a.dissipated_energy_j;
        work += a.supplied_work_j;
        assert!(a.dissipated_energy_j >= 0.0 && b.dissipated_energy_j >= 0.0);
        let mapping = [0, 1, 2, 3, 6, 7, 4, 5];
        for (i, j) in mapping.into_iter().enumerate() {
            assert!((first.state()[i] - second.state()[j]).abs() < 1e-10);
        }
        assert!((a.stored_energy_j - b.stored_energy_j).abs() < 1e-10 * initial);
        assert!((first.radiation_force(0).unwrap() - 2.4 * first.state()[7]).abs() < 1e-13);
        assert!((first.radiation_force(1).unwrap() + 1.7 * first.state()[7]).abs() < 1e-13);
    }
    let reference = cross_reference();
    for (&a, &b) in first.state().iter().zip(&reference[..8]) {
        assert!((a - b).abs() < 2e-10, "{a:e} vs {b:e}");
    }
    assert!((heat - reference[8]).abs() < 2e-9 * initial);
    assert!((work - reference[9]).abs() < 2e-9 * initial);
    assert!((first.stored_energy_j() + heat - work - initial).abs() < 2e-9 * initial);
    assert_eq!(first.relaxation_observation().unwrap().states, 2);
    assert_eq!(first.radiation_observation().unwrap().poles, 1);
}

fn nonlinear_model() -> ExponentialImpactSystem {
    let gate = CancelGate::new_clock_free();
    let source = RateFelt::new(
        WoolFelt::new(30000.0, 0.2, 2.0, 3.0, 0.15, 0.8).unwrap(),
        vec![CreepBranch {
            modulus_pa: 60000.0,
            viscosity_pa_s: 1200.0,
        }],
    )
    .unwrap();
    let pad = FeltPad::from_material(0.001, 0.006, 0.0001, vec![1.0], 0.0, &source).unwrap();
    let mut body = linear(10.0, 0.0);
    body.initial[0].velocity_m_sqrt_kg_per_s = 0.3;
    ExponentialImpactSystem::new_with_felt(
        vec![body],
        vec![],
        vec![pad],
        vec![],
        vec![],
        config(1e-5),
        &gate,
    )
    .unwrap()
    .with_relaxation_branches(
        vec![arm(3, 12000.0, 0.01)],
        InitialMemory::Relaxed,
        1,
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &radiation::Model {
            ports: 1,
            poles: vec![radiation::Pole {
                omega: 200.0,
                zeta: 0.1,
                coupling: vec![10.0],
            }],
        },
        &[0],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn g4_felt_and_both_memories_keep_the_same_storage_and_joint_commit_after_retry() {
    let gate = CancelGate::new_clock_free();
    let (mut run, mut fresh) = (nonlinear_model(), nonlinear_model());
    assert_eq!(run.state().len(), 6); // mechanical(2), felt(1), Maxwell(1), air(2)
    let initial = run.stored_energy_j();
    let (mut heat, mut crush) = (0.0, 0.0);
    for _ in 0..100 {
        let frame = run.step(&[0.0], &gate).unwrap();
        fresh.step(&[0.0], &gate).unwrap();
        heat += frame.dissipated_energy_j;
        crush += frame.felt_crush_loss_j;
        assert!((frame.stored_energy_j - run.stored_energy_j()).abs() < 1e-12);
    }
    assert!(heat >= crush && crush > 0.0);
    assert!((run.stored_energy_j() + heat - initial).abs() < 2e-9 * initial);
    assert!(run.relaxation_observation().unwrap().stored_energy_j > 0.0);
    assert!(run.radiation_observation().unwrap().stored_energy_j > 0.0);
    let before = run.state().to_vec();
    let history = run.felt_history(0).unwrap();
    let energy = run.stored_energy_j();
    let samples = run.samples();
    run.set_iteration_limit(0).unwrap();
    assert!(run.step(&[0.0], &gate).is_err());
    assert_eq!(run.state(), before);
    assert_eq!(run.felt_history(0).unwrap(), history);
    assert_eq!(run.stored_energy_j(), energy);
    assert_eq!(run.samples(), samples);
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
    assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
}

#[test]
fn g0_full_memory_capacity_source_mapping_clock_and_cold_lifecycle_are_admitted() {
    let gate = CancelGate::new_clock_free();
    let make = || {
        ExponentialImpactSystem::new(
            vec![linear(5.0, 0.001)],
            vec![],
            vec![],
            config(0.01),
            &gate,
        )
        .unwrap()
    };
    let native = make();
    assert!(
        native
            .with_relaxation_branches(
                vec![arm(2, 10.0, 0.1); 126],
                InitialMemory::Relaxed,
                126,
                &gate,
            )
            .is_err()
    ); // the old two-state image fits, the COMPLETE 128-state flow does not
    assert!(
        make()
            .with_relaxation_branches(vec![arm(2, 1e8, 0.1)], InitialMemory::Relaxed, 1, &gate,)
            .is_err()
    ); // original combined instantaneous spectrum, not the bare-body band
    let mut momentum = arm(2, 10.0, 0.1);
    momentum.projection = [0.0, 1.0].to_vec();
    assert!(
        make()
            .with_relaxation_branches(vec![momentum], InitialMemory::Relaxed, 1, &gate,)
            .is_err()
    );
    assert!(
        make()
            .with_relaxation_branches(vec![arm(2, 10.0, -0.1)], InitialMemory::Relaxed, 1, &gate,)
            .is_err()
    );
    let attached = make()
        .with_relaxation_branches(vec![arm(2, 10.0, 0.1)], InitialMemory::Unrelaxed, 1, &gate)
        .unwrap();
    assert!(
        attached
            .with_relaxation_branches(vec![arm(3, 10.0, 0.1)], InitialMemory::Relaxed, 1, &gate,)
            .is_err()
    );
    let mut late = make();
    late.step(&[0.0], &gate).unwrap();
    assert!(
        late.with_relaxation_branches(vec![arm(2, 10.0, 0.1)], InitialMemory::Relaxed, 1, &gate,)
            .is_err()
    );
    let zero = radiation::Model {
        ports: 1,
        poles: vec![radiation::Pole {
            omega: 10.0,
            zeta: 0.3,
            coupling: vec![0.0],
        }],
    };
    let mut untouched = make();
    let mut no_op = make()
        .with_radiation_load(&zero, &[0], 1, &gate)
        .unwrap()
        .with_relaxation_branches(vec![], InitialMemory::Relaxed, 0, &gate)
        .unwrap();
    assert_eq!(no_op.state().len(), 2);
    untouched.step(&[0.1], &gate).unwrap();
    no_op.step(&[0.1], &gate).unwrap();
    assert_eq!(untouched.state(), no_op.state());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        make().with_relaxation_branches(
            vec![arm(2, 10.0, 0.1)],
            InitialMemory::Relaxed,
            1,
            &stopped,
        ),
        Err(ImpactError::Cancelled)
    ));
}
