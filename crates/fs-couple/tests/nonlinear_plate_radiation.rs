//! G1/G3: full mutual radiation, nonlinear mechanics and Maxwell history share time.
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::thin_plate::VkBody;
use fs_exec::CancelGate;
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_phs::RelaxationBranch;

fn body(nonlinear: bool) -> VkBody {
    VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![2.0, 3.0],
            channels: if nonlinear {
                vec![StressChannel {
                    coefficient: 30.0,
                    coupling: vec![1.0, 0.2, 0.2, 0.6],
                }]
            } else {
                vec![]
            },
        },
        &[0.4, -0.1, -0.1, 0.2],
        &[1.0, 0.3],
        vec![0.02, -0.01],
    )
    .unwrap()
    .with_relaxation_branches(
        vec![RelaxationBranch {
            projection: vec![1.0, 0.0, -0.5, 0.0],
            stiffness: 3.5,
            relaxation_time_s: 0.15,
        }],
        1,
    )
    .unwrap()
}
fn model() -> Model {
    Model {
        ports: 2,
        poles: vec![
            Pole {
                omega: 4.0,
                zeta: 0.2,
                coupling: vec![1.7, -0.9],
            },
            Pole {
                omega: 6.0,
                zeta: 0.6,
                coupling: vec![0.4, 1.2],
            },
        ],
    }
}

// Independently assembled physical q/v/Maxwell/acoustic equations. Last two
// entries independently integrate force work and every physical damper's loss.
fn derivative(x: &[f64; 11], force: f64, nonlinear: bool) -> [f64; 11] {
    let q = [x[0], x[2]];
    let v = [x[1], x[3]];
    let eq = [q[0] + 0.2 * q[1], 0.2 * q[0] + 0.6 * q[1]];
    let cubic = if nonlinear {
        30.0 * (q[0] * eq[0] + q[1] * eq[1])
    } else {
        0.0
    };
    let memory = 3.5_f64.sqrt() * (q[0] - 0.5 * q[1]) - x[4];
    [
        v[0],
        -4.0 * q[0] - cubic * eq[0] - 0.4 * v[0] + 0.1 * v[1]
            - 3.5_f64.sqrt() * memory
            - 1.7 * x[6]
            - 0.4 * x[8]
            + force,
        v[1],
        -9.0 * q[1] - cubic * eq[1] + 0.1 * v[0] - 0.2 * v[1]
            + 0.5 * 3.5_f64.sqrt() * memory
            + 0.9 * x[6]
            - 1.2 * x[8]
            + 0.3 * force,
        memory / 0.15,
        x[6],
        -16.0 * x[5] - 1.6 * x[6] + 1.7 * v[0] - 0.9 * v[1],
        x[8],
        -36.0 * x[7] - 7.2 * x[8] + 0.4 * v[0] + 1.2 * v[1],
        force * (v[0] + 0.3 * v[1]),
        0.4 * v[0] * v[0] - 0.2 * v[0] * v[1]
            + 0.2 * v[1] * v[1]
            + memory * memory / 0.15
            + 1.6 * x[6] * x[6]
            + 7.2 * x[8] * x[8],
    ]
}
fn rk4(x: &mut [f64; 11], force: f64, h: f64, nonlinear: bool) {
    let a = derivative(x, force, nonlinear);
    let b = derivative(
        &core::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
        force,
        nonlinear,
    );
    let c = derivative(
        &core::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
        force,
        nonlinear,
    );
    let d = derivative(&core::array::from_fn(|i| x[i] + h * c[i]), force, nonlinear);
    for i in 0..11 {
        x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
    }
}
fn energy(x: &[f64], nonlinear: bool) -> f64 {
    let strain = x[0] * x[0] + 0.4 * x[0] * x[2] + 0.6 * x[2] * x[2];
    let memory = 3.5_f64.sqrt() * (x[0] - 0.5 * x[2]) - x[4];
    0.5 * (x[1] * x[1]
        + x[3] * x[3]
        + 4.0 * x[0] * x[0]
        + 9.0 * x[2] * x[2]
        + memory * memory
        + x[6] * x[6]
        + 16.0 * x[5] * x[5]
        + x[8] * x[8]
        + 36.0 * x[7] * x[7])
        + if nonlinear {
            7.5 * strain * strain
        } else {
            0.0
        }
}

#[test]
fn complete_loaded_plate_matches_independent_physics_and_work() {
    let gate = CancelGate::new();
    for nonlinear in [false, true] {
        // The nonlinear Gonzalez trajectory has second-order time error; use
        // a resolved clock for this absolute heat/work comparison. The test
        // below separately measures that order over three coarser clocks.
        let steps = if nonlinear { 800 } else { 200 };
        let dt = 0.2 / steps as f64;
        let mut p = body(nonlinear)
            .with_radiation_load(&model(), 2, &gate)
            .unwrap()
            .into_prepared(dt)
            .unwrap();
        assert_eq!(p.state().len(), 9);
        let mut reference = [0.0; 11];
        let (mut work, mut heat) = (0.0, 0.0);
        for k in 0..steps {
            let force = if k < steps / 2 { 4.0 } else { -0.7 };
            let q = [p.state()[0], p.state()[2]];
            let record = p.drive_and_radiate(force, 1.2, 1.5).unwrap();
            for _ in 0..16 {
                rk4(&mut reference, force, dt / 16.0, nonlinear);
            }
            let independent_work = force * ((p.state()[0] - q[0]) + 0.3 * (p.state()[2] - q[1]));
            assert!((record.work.supplied - independent_work).abs() < 2e-12);
            work += record.work.supplied;
            heat += record.work.dissipated;
            assert!((p.stored_energy_j() - energy(p.state(), nonlinear)).abs() < 1e-13);
            assert!((energy(p.state(), nonlinear) + heat - work).abs() < 1e-10);
        }
        let error = p
            .state()
            .iter()
            .zip(reference)
            .take(9)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            error < if nonlinear { 3e-8 } else { 2e-11 },
            "nonlinear={nonlinear}, error={error}"
        );
        assert!(
            (heat - reference[10]).abs() < if nonlinear { 3e-8 } else { 2e-11 },
            "nonlinear={nonlinear}: heat={heat}, reference={}, error={}",
            reference[10],
            heat - reference[10]
        );
        assert!(
            (work - reference[9]).abs() < if nonlinear { 3e-8 } else { 2e-11 },
            "nonlinear={nonlinear}: work={work}, reference={}, error={}",
            reference[9],
            work - reference[9]
        );
        assert!(p.state()[6].abs() > 1e-3 && p.state()[8].abs() > 1e-3);
    }
}

#[test]
fn nonlinear_radiation_trajectory_converges_with_the_original_second_order_owner() {
    let gate = CancelGate::new();
    let mut reference = [0.; 11];
    for _ in 0..12800 {
        rk4(&mut reference, 4.0, 0.2 / 12800.0, true);
    }
    let mut errors = Vec::new();
    for steps in [100, 200, 400] {
        let dt = 0.2 / steps as f64;
        let mut p = body(true)
            .with_radiation_load(&model(), 2, &gate)
            .unwrap()
            .into_prepared(dt)
            .unwrap();
        for _ in 0..steps {
            p.drive_and_radiate(4.0, 1.2, 1.5).unwrap();
        }
        errors.push(
            p.state()
                .iter()
                .zip(reference)
                .take(9)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}

#[test]
fn load_keeps_existing_motion_and_memory_and_refuses_incomplete_physics() {
    let gate = CancelGate::new();
    let mut original = body(true);
    for _ in 0..10 {
        original.drive_and_radiate(2.0, 0.001, 1.2, 1.5).unwrap();
    }
    let before = original.state().to_vec();
    let h = original.stored_energy_j();
    let loaded = original.with_radiation_load(&model(), 2, &gate).unwrap();
    assert_eq!(&loaded.state()[..before.len()], before);
    assert_eq!(&loaded.state()[before.len()..], &[0.; 4]);
    assert_eq!(loaded.stored_energy_j(), h);
    assert!(body(true).with_radiation_load(&model(), 1, &gate).is_err());
    let mut wrong = model();
    wrong.ports = 1;
    for p in &mut wrong.poles {
        p.coupling.truncate(1);
    }
    assert!(body(true).with_radiation_load(&wrong, 2, &gate).is_err());
    let mut active = model();
    active.poles[0].zeta = -0.1;
    assert!(body(true).with_radiation_load(&active, 2, &gate).is_err());
    let mut large = model();
    large.poles = vec![large.poles[0].clone(); 64];
    assert!(body(true).with_radiation_load(&large, 64, &gate).is_err());
    let mut zero = model();
    for p in &mut zero.poles {
        p.coupling.fill(0.0);
    }
    let no_op = body(true).with_radiation_load(&zero, 2, &gate).unwrap();
    assert_eq!(no_op.state().len(), 5);
    gate.request();
    assert!(body(true).with_radiation_load(&model(), 2, &gate).is_err());
}
