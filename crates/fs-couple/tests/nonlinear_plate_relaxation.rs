//! G1/G3: nonlinear plate motion retains causal cross-mode material history.
use fs_couple::thin_plate::VkBody;
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_phs::RelaxationBranch;

fn body(nonlinear: bool, force: &[f64]) -> VkBody {
    let storage = SosModalStorage {
        omegas: vec![2.0, 3.0],
        channels: if nonlinear {
            vec![StressChannel {
                coefficient: 30.0,
                coupling: vec![1.0, 0.2, 0.2, 0.6],
            }]
        } else {
            vec![]
        },
    };
    VkBody::from_modal_storage_with_damping(
        storage,
        &[0.4, -0.1, -0.1, 0.2],
        force,
        vec![0.02, -0.01],
    )
    .unwrap()
}

fn branch() -> RelaxationBranch {
    RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 3.5,
        relaxation_time_s: 0.15,
    }
}

// Independent continuous mechanical/internal-strain equations. No pHS helper,
// finite-difference force, modal eigensolver or production time step is used.
fn derivative(x: [f64; 5], force: f64, nonlinear: bool) -> [f64; 5] {
    let [q0, v0, q1, v1, z] = x;
    let root = 3.5_f64.sqrt();
    let memory_effort = root * (q0 - 0.5 * q1) - z;
    let form = q0 * q0 + 0.4 * q0 * q1 + 0.6 * q1 * q1;
    let weight = if nonlinear { 30.0 * form } else { 0.0 };
    [
        v0,
        force - 4.0 * q0 - weight * (q0 + 0.2 * q1) - 0.4 * v0 + 0.1 * v1 - root * memory_effort,
        v1,
        0.3 * force - 9.0 * q1 - weight * (0.2 * q0 + 0.6 * q1) + 0.1 * v0 - 0.2 * v1
            + 0.5 * root * memory_effort,
        memory_effort / 0.15,
    ]
}

fn reference(nonlinear: bool, force: f64, duration: f64) -> [f64; 5] {
    let steps = 12_000;
    let h = duration / steps as f64;
    let mut x = [0.0; 5];
    for _ in 0..steps {
        let a = derivative(x, force, nonlinear);
        let b = derivative(
            std::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
            force,
            nonlinear,
        );
        let c = derivative(
            std::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
            force,
            nonlinear,
        );
        let d = derivative(std::array::from_fn(|i| x[i] + h * c[i]), force, nonlinear);
        x = std::array::from_fn(|i| x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

fn energy(x: &[f64]) -> f64 {
    let form = x[0] * x[0] + 0.4 * x[0] * x[2] + 0.6 * x[2] * x[2];
    let memory = 3.5_f64.sqrt() * (x[0] - 0.5 * x[2]) - x[4];
    0.5 * (x[1] * x[1] + x[3] * x[3] + 4.0 * x[0] * x[0] + 9.0 * x[2] * x[2])
        + 7.5 * form * form
        + 0.5 * memory * memory
}

#[test]
fn linear_material_memory_is_exact_under_subdivision_and_matches_continuous_equations() {
    let make = || {
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![branch()], 1)
            .unwrap()
    };
    let mut whole = make();
    let mut split = make();
    whole.drive_and_radiate(2.0, 0.6, 1.2, 1.5).unwrap();
    for _ in 0..12 {
        split.drive_and_radiate(2.0, 0.05, 1.2, 1.5).unwrap();
    }
    let expected = reference(false, 2.0, 0.6);
    assert_eq!(whole.state().len(), 5);
    assert!(
        whole.state()[4].abs() > 0.01,
        "material history must actually evolve"
    );
    for ((a, b), truth) in whole.state().iter().zip(split.state()).zip(expected) {
        assert!((a - b).abs() < 2e-12);
        assert!((a - truth).abs() < 2e-11, "{a} versus {truth}");
    }
    assert!((whole.volume_velocity() - (0.02 * expected[1] - 0.01 * expected[3])).abs() < 2e-12);
}

#[test]
fn localized_material_memory_keeps_the_unstrained_combination_lossless() {
    let make = |memory| {
        let result = VkBody::from_modal_storage_with_damping(
            SosModalStorage {
                omegas: vec![6.0, 6.0],
                channels: vec![],
            },
            &[0.0; 4],
            &[1.0, 1.0],
            vec![0.02, 0.02],
        )
        .unwrap();
        if memory {
            result
                .with_relaxation_branches(
                    vec![RelaxationBranch {
                        projection: vec![1.0, 0.0, -1.0, 0.0],
                        stiffness: 9.0,
                        relaxation_time_s: 0.07,
                    }],
                    1,
                )
                .unwrap()
        } else {
            result
        }
    };
    let mut dry = make(false);
    let mut material = make(true);
    dry.drive_and_radiate(0.8, 0.2, 1.2, 1.0).unwrap();
    material.drive_and_radiate(0.8, 0.2, 1.2, 1.0).unwrap();
    for (actual, expected) in material.state()[..4].iter().zip(dry.state()) {
        assert!((actual - expected).abs() < 2e-12);
    }
    assert!(material.state()[4].abs() < 2e-13);
    assert!((material.stored_energy_j() - dry.stored_energy_j()).abs() < 2e-12);
}

#[test]
fn nonlinear_memory_converges_and_dissipates_the_complete_storage() {
    let expected = reference(true, 2.0, 0.6);
    let mut errors = Vec::new();
    for steps in [200, 400, 800] {
        let mut plate = body(true, &[1.0, 0.3])
            .with_relaxation_branches(vec![branch()], 1)
            .unwrap();
        let dt = 0.6 / steps as f64;
        for _ in 0..steps {
            let previous_flow = plate.volume_velocity();
            let (pressure, flow) = plate.drive_and_radiate(2.0, dt, 1.2, 1.5).unwrap();
            let mean_pressure =
                1.2 * (flow - previous_flow) / (2.0 * std::f64::consts::PI * 1.5 * dt);
            assert!((pressure - mean_pressure).abs() < 2e-12);
        }
        let error = plate
            .state()
            .iter()
            .zip(expected)
            .map(|(actual, truth)| (actual - truth).powi(2))
            .sum::<f64>()
            .sqrt();
        errors.push(error);
        assert!((plate.stored_energy_j() - energy(plate.state())).abs() < 2e-12);
        let initial_energy = plate.stored_energy_j();
        let mut previous_energy = initial_energy;
        for _ in 0..40 {
            plate.drive_and_radiate(0.0, 0.002, 1.2, 1.5).unwrap();
            let next_energy = plate.stored_energy_j();
            assert!(next_energy <= previous_energy + 2e-12);
            previous_energy = next_energy;
        }
        assert!(previous_energy < 0.99 * initial_energy);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 2e-5, "{errors:?}");
}

#[test]
fn memory_admission_preserves_original_state_and_never_loads_momentum() {
    let mut started = body(false, &[1.0, 0.3]);
    started.drive_and_radiate(1.0, 0.01, 1.2, 1.0).unwrap();
    assert!(started.with_relaxation_branches(vec![branch()], 1).is_err());
    let attached = body(false, &[1.0, 0.3])
        .with_relaxation_branches(vec![branch()], 1)
        .unwrap();
    assert!(
        attached
            .with_relaxation_branches(vec![branch()], 1)
            .is_err()
    );
    assert!(
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![branch()], 0)
            .is_err()
    );
    let mut invalid = branch();
    invalid.projection[1] = 1.0;
    assert!(
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![invalid], 1)
            .is_err()
    );
    let mut invalid = branch();
    invalid.projection.pop();
    assert!(
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![invalid], 1)
            .is_err()
    );
    let mut invalid = branch();
    invalid.stiffness = -1.0;
    assert!(
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![invalid], 1)
            .is_err()
    );
    let mut invalid = branch();
    invalid.relaxation_time_s = f64::NAN;
    assert!(
        body(false, &[1.0, 0.3])
            .with_relaxation_branches(vec![invalid], 1)
            .is_err()
    );
    let mut raw = body(false, &[1.0, 0.3]);
    let mut empty = body(false, &[1.0, 0.3])
        .with_relaxation_branches(vec![], 0)
        .unwrap();
    raw.drive_and_radiate(1.0, 0.02, 1.2, 1.0).unwrap();
    empty.drive_and_radiate(1.0, 0.02, 1.2, 1.0).unwrap();
    assert_eq!(raw.state(), empty.state());
    let mut nonlinear = body(true, &[1.0, 0.3])
        .with_relaxation_branches(vec![branch()], 1)
        .unwrap();
    let before = nonlinear.state().to_vec();
    assert!(nonlinear.drive_and_radiate(1.0, 100.0, 1.2, 1.0).is_err());
    assert_eq!(nonlinear.state(), before);
}
