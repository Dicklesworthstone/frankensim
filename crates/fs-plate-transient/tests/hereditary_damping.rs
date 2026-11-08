use fs_plate::viscosity::{PlateBendingViscosity, PlateViscosityBudget, ViscousPlateSection};
use fs_plate::{PlateMesh, PlateModel};
use fs_plate_transient::PlateDynamicsBudget;
use fs_plate_transient::relaxation::{
    PlateRelaxationBranch, PlateRelaxationBudget, PlateRelaxationError, RelaxingPlate,
};
use fs_sparse::{Coo, Csr, LdltOptions};

fn matrix(n: usize, values: &[f64]) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[n * i + j] != 0.0 {
                coo.push(i, j, values[n * i + j]);
            }
        }
    }
    coo.assemble()
}
fn mesh() -> PlateMesh {
    PlateMesh {
        nodes: vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        tris: vec![[0, 1, 2]],
    }
}
fn plate() -> PlateModel {
    let mut dof_map = vec![None; 9];
    dof_map[0] = Some(0);
    PlateModel {
        m: matrix(1, &[1.3]),
        k: matrix(1, &[17.0]),
        dof_map,
        free: 1,
    }
}
fn viscosity_budget(plate: &PlateModel) -> PlateViscosityBudget {
    PlateViscosityBudget {
        max_elements: 1,
        max_dofs: plate.free,
        max_full_dofs: 9,
        max_assembly_entries: 81,
    }
}
fn viscosity<'a>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    tensor: [f64; 9],
) -> PlateBendingViscosity<'a> {
    PlateBendingViscosity::assemble(
        mesh,
        plate,
        vec![ViscousPlateSection::from_resultant(tensor).unwrap()],
        viscosity_budget(plate),
        &mut || false,
    )
    .unwrap()
}
fn branch<'a>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    tau: f64,
    scale: f64,
) -> PlateRelaxationBranch<'a> {
    let tensor = [1.0, 0.2, 0.1, 0.2, 2.0, -0.05, 0.1, -0.05, 0.8].map(|x| scale * x);
    PlateRelaxationBranch::assemble(
        mesh,
        plate,
        vec![tensor],
        tau,
        viscosity_budget(plate),
        &mut || false,
    )
    .unwrap()
}
fn budget(plate: &PlateModel, branches: usize) -> PlateRelaxationBudget {
    PlateRelaxationBudget {
        dynamics: PlateDynamicsBudget {
            max_dofs: plate.free,
            max_nonzeros: plate.m.nnz() + plate.k.nnz(),
            max_full_dofs: 9,
            max_load_parameters: 0,
        },
        max_branches: branches,
        max_effective_entries: (2 + branches) * plate.free * plate.free,
        max_state_scalars: (7 + 2 * branches) * plate.free,
    }
}
fn quadratic(matrix: &Csr, x: &[f64]) -> f64 {
    (0..x.len())
        .map(|i| {
            (0..x.len())
                .map(|j| x[i] * matrix.get(i, j) * x[j])
                .sum::<f64>()
        })
        .sum()
}

// G1: independent continuous ODE, with separate relaxation/viscous heat states.
// This does not use the production ramp, effective matrix, factor or integrator.
fn exact_reference(stiffness: [f64; 2], taus: [f64; 2], viscosity: f64, duration: f64) -> [f64; 6] {
    let count = 20000;
    let h = duration / count as f64;
    let derivative = |x: [f64; 6]| {
        [
            x[1],
            -(17.0 * x[0] + viscosity * x[1] + stiffness[0] * x[2] + stiffness[1] * x[3]) / 1.3,
            x[1] - x[2] / taus[0],
            x[1] - x[3] / taus[1],
            stiffness[0] * x[2] * x[2] / taus[0] + stiffness[1] * x[3] * x[3] / taus[1],
            viscosity * x[1] * x[1],
        ]
    };
    let mut x = [0.02, 0.1, 0.02, 0.02, 0.0, 0.0];
    for _ in 0..count {
        let a = derivative(x);
        let b = derivative(core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = derivative(core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = derivative(core::array::from_fn(|i| x[i] + h * c[i]));
        x = core::array::from_fn(|i| x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn original_plate_memory_matches_continuous_maxwell_and_kelvin_motion() {
    let mesh = mesh();
    let plate = plate();
    let taus = [0.035, 0.6];
    let source = vec![
        branch(&mesh, &plate, taus[0], 0.4),
        branch(&mesh, &plate, taus[1], 0.7),
    ];
    let stiffness = core::array::from_fn(|i| source[i].excess_stiffness().get(0, 0));
    let viscosity = viscosity(
        &mesh,
        &plate,
        [0.08, 0.01, 0.0, 0.01, 0.16, 0.0, 0.0, 0.0, 0.04],
    );
    for with_viscosity in [false, true] {
        let damping = if with_viscosity {
            viscosity.matrix().get(0, 0)
        } else {
            0.0
        };
        let reference = exact_reference(stiffness, taus, damping, 0.8);
        let mut errors = Vec::new();
        for steps in [400, 800, 1600] {
            let branches = vec![
                branch(&mesh, &plate, taus[0], 0.4),
                branch(&mesh, &plate, taus[1], 0.7),
            ];
            let owner = if with_viscosity {
                let mut caps = budget(&plate, 2);
                caps.max_effective_entries += viscosity.matrix().nnz();
                RelaxingPlate::prepare_with_bending_viscosity(
                    &mesh,
                    &viscosity,
                    branches,
                    0.8 / steps as f64,
                    1e-12,
                    caps,
                    viscosity_budget(&plate),
                    &LdltOptions::default(),
                    &mut || false,
                )
            } else {
                RelaxingPlate::prepare(
                    &plate,
                    branches,
                    0.8 / steps as f64,
                    1e-12,
                    budget(&plate, 2),
                    &LdltOptions::default(),
                    &mut || false,
                )
            }
            .unwrap();
            assert_eq!(owner.bending_viscosity().is_some(), with_viscosity);
            let mut state = owner.state(&[0.02], &[0.1], &[0.02, 0.02]).unwrap();
            let initial = state.mechanical_energy_j();
            for _ in 0..steps {
                let energy = state.mechanical_energy_j();
                let step = state.advance(&[0.0], &mut || false).unwrap();
                assert!(step.dissipated_j >= 0.0);
                assert!(step.relaxation_dissipated_j >= 0.0 && step.viscous_dissipated_j >= 0.0);
                assert_eq!(
                    step.dissipated_j,
                    step.relaxation_dissipated_j + step.viscous_dissipated_j
                );
                assert!(step.balance_residual_j.abs() < 2e-16);
                assert!(state.mechanical_energy_j() <= energy + 2e-16);
            }
            let actual = [
                state.displacement()[0],
                state.velocity()[0],
                state.branch_strains()[0],
                state.branch_strains()[1],
                state.relaxation_dissipated_j(),
                state.viscous_dissipated_j(),
            ];
            errors.push(
                actual
                    .iter()
                    .zip(reference)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0_f64, f64::max),
            );
            assert!((state.mechanical_energy_j() + state.dissipated_j() - initial).abs() < 2e-15);
            assert!(
                state.dissipated_j() > 0.2 * initial,
                "material must measurably damp the trajectory"
            );
        }
        assert!(
            errors[0] > 3.8 * errors[1] && errors[1] > 3.8 * errors[2],
            "second-order convergence: {errors:?}"
        );
        assert!(errors[2] < 2e-7, "continuous ODE mismatch: {errors:?}");
    }
}

#[test]
fn held_shape_relaxes_exactly_even_when_tau_is_far_below_the_time_step() {
    let mesh = mesh();
    let plate = plate();
    let dt = 0.125;
    for rate in [1e-5, 1.0, 1e4] {
        let branch = branch(&mesh, &plate, dt / rate, 0.4);
        let k = branch.excess_stiffness().get(0, 0);
        let owner = RelaxingPlate::prepare(
            &plate,
            vec![branch],
            dt,
            1e-12,
            budget(&plate, 1),
            &LdltOptions::default(),
            &mut || false,
        )
        .unwrap();
        let mut state = owner.state(&[0.03], &[0.0], &[0.01]).unwrap();
        let mean_force = 17.0 * 0.03 + k * (-(-rate).exp_m1() / rate) * 0.01;
        let expected_heat = 0.5 * k * 0.01_f64.powi(2) * (-(-2.0 * rate).exp_m1());
        let step = state.advance(&[mean_force], &mut || false).unwrap();
        assert!((state.displacement()[0] - 0.03).abs() < 2e-17);
        assert!(state.velocity()[0].abs() < 2e-16);
        assert!((state.branch_strains()[0] - 0.01 * (-rate).exp()).abs() < 2e-17);
        assert!((step.dissipated_j - expected_heat).abs() < 2e-11 * expected_heat);
        assert!((step.stored_change_j + expected_heat).abs() < 2e-16);
    }
}

#[test]
fn full_consistent_inertia_and_anisotropic_memory_close_independent_forced_work() {
    let mesh = mesh();
    let mut dof_map = vec![None; 9];
    dof_map[0] = Some(0);
    dof_map[4] = Some(1);
    let plate = PlateModel {
        m: matrix(2, &[2.0, 0.4, 0.4, 1.0]),
        k: matrix(2, &[7.0, -2.0, -2.0, 5.0]),
        dof_map,
        free: 2,
    };
    let branch = branch(&mesh, &plate, 0.17, 0.3);
    let stiffness = branch.excess_stiffness().clone();
    assert!(
        stiffness.get(0, 1).abs() > 1e-3,
        "the physical memory must couple coordinates"
    );
    let viscosity = viscosity(
        &mesh,
        &plate,
        [0.08, 0.01, 0.0, 0.01, 0.16, 0.0, 0.0, 0.0, 0.04],
    );
    assert!(viscosity.matrix().get(0, 1).abs() > 1e-3);
    let mut caps = budget(&plate, 1);
    caps.max_effective_entries += viscosity.matrix().nnz();
    let owner = RelaxingPlate::prepare_with_bending_viscosity(
        &mesh,
        &viscosity,
        vec![branch],
        0.005,
        1e-12,
        caps,
        viscosity_budget(&plate),
        &LdltOptions::default(),
        &mut || false,
    )
    .unwrap();
    let mut state = owner
        .state(&[0.02, -0.03], &[0.1, 0.05], &[0.01, 0.02])
        .unwrap();
    let energy = |u: &[f64], v: &[f64], q: &[f64]| {
        0.5 * (quadratic(&plate.k, u) + quadratic(&plate.m, v) + quadratic(&stiffness, q))
    };
    let initial = energy(
        state.displacement(),
        state.velocity(),
        state.branch_strains(),
    );
    let mut work = 0.0;
    for step in 0..200 {
        let old = state.displacement().to_vec();
        let old_velocity = state.velocity().to_vec();
        let t = (step as f64 + 0.5) * owner.step_s();
        let force = [0.03 * (3.0 * t).sin(), -0.04 * (2.0 * t).cos()];
        let record = state.advance(&force, &mut || false).unwrap();
        let mean_velocity: Vec<_> = old_velocity
            .iter()
            .zip(state.velocity())
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let viscous_heat = owner.step_s() * quadratic(viscosity.matrix(), &mean_velocity);
        assert!((record.viscous_dissipated_j - viscous_heat).abs() < 1e-17);
        work += force
            .iter()
            .zip(old.iter().zip(state.displacement()))
            .map(|(f, (a, b))| f * (b - a))
            .sum::<f64>();
        assert!(record.balance_residual_j.abs() < 3e-16);
    }
    let final_energy = energy(
        state.displacement(),
        state.velocity(),
        state.branch_strains(),
    );
    assert!((final_energy + state.dissipated_j() - initial - work).abs() < 5e-15);
    assert!((state.supplied_work_j() - work).abs() < 1e-16);
}

#[test]
fn rigid_translation_stays_undamped_and_refused_advances_preserve_history() {
    let mesh = mesh();
    let n = 9;
    let mut identity = vec![0.0; n * n];
    for i in 0..n {
        identity[n * i + i] = 1.0;
    }
    let plate = PlateModel {
        m: matrix(n, &identity),
        k: matrix(n, &vec![0.0; n * n]),
        dof_map: (0..n).map(Some).collect(),
        free: n,
    };
    let branch = PlateRelaxationBranch::assemble(
        &mesh,
        &plate,
        vec![[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
        0.03,
        PlateViscosityBudget {
            max_elements: 1,
            max_dofs: 9,
            max_full_dofs: 9,
            max_assembly_entries: 81,
        },
        &mut || false,
    )
    .unwrap();
    let mut caps = budget(&plate, 1);
    caps.max_state_scalars -= 1;
    assert!(matches!(
        RelaxingPlate::prepare(
            &plate,
            vec![],
            0.01,
            1e-12,
            PlateRelaxationBudget {
                max_state_scalars: 0,
                ..caps
            },
            &LdltOptions::default(),
            &mut || false
        ),
        Err(PlateRelaxationError::Budget(_))
    ));
    let viscosity = viscosity(&mesh, &plate, [0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let mut caps = budget(&plate, 1);
    caps.max_effective_entries += viscosity.matrix().nnz();
    let owner = RelaxingPlate::prepare_with_bending_viscosity(
        &mesh,
        &viscosity,
        vec![branch],
        0.01,
        1e-12,
        caps,
        viscosity_budget(&plate),
        &LdltOptions::default(),
        &mut || false,
    )
    .unwrap();
    let mut velocity = vec![0.0; n];
    for i in [0, 3, 6] {
        velocity[i] = 1.0;
    }
    let mut state = owner
        .state(&vec![0.0; n], &velocity, &vec![0.0; n])
        .unwrap();
    for _ in 0..30 {
        state.advance(&vec![0.0; n], &mut || false).unwrap();
    }
    assert!(state.dissipated_j() < 1e-28);
    for (&a, &b) in state.velocity().iter().zip(&velocity) {
        assert!((a - b).abs() < 2e-14);
    }
    let before = (
        state.displacement().to_vec(),
        state.velocity().to_vec(),
        state.branch_strains().to_vec(),
        state.time_s(),
        state.dissipated_j(),
        state.supplied_work_j(),
    );
    let mut polls = 0;
    let error = state
        .advance(&vec![0.0; n], &mut || {
            polls += 1;
            polls > n + 5
        })
        .unwrap_err();
    assert!(matches!(error, PlateRelaxationError::Cancelled));
    assert_eq!(state.displacement(), before.0);
    assert_eq!(state.velocity(), before.1);
    assert_eq!(state.branch_strains(), before.2);
    assert_eq!(
        (
            state.time_s(),
            state.dissipated_j(),
            state.supplied_work_j()
        ),
        (before.3, before.4, before.5)
    );
    assert!(state.advance(&[f64::NAN; 9], &mut || false).is_err());
    assert_eq!(state.branch_strains(), before.2);
}

#[test]
fn energy_ledger_exposes_rounding_of_the_actual_accepted_endpoints() {
    // An ordinary small-deflection equilibrium is a large background relative
    // to this perturbation. Source work follows the intended affine path, but
    // the actual displacement cannot acquire its sub-ulp elastic increment.
    let plate = plate();
    let owner = RelaxingPlate::prepare(
        &plate,
        vec![],
        0.125,
        1e-12,
        budget(&plate, 0),
        &LdltOptions::default(),
        &mut || false,
    )
    .unwrap();
    let mut state = owner.state(&[0.01], &[1e-20], &[]).unwrap();
    let initial_velocity = state.velocity()[0];
    let record = state.advance(&[17.0 * 0.01], &mut || false).unwrap();
    let velocity = state.velocity()[0];
    let kinetic = 0.5 * 1.3 * (velocity - initial_velocity) * (velocity + initial_velocity);
    assert_eq!(state.displacement(), &[0.01]);
    assert_eq!(record.stored_change_j, 0.0);
    assert!((record.kinetic_change_j - kinetic).abs() < 1e-54);
    assert!(record.supplied_work_j > 1e-23);
    assert!((record.balance_residual_j + record.supplied_work_j - kinetic).abs() < 1e-36);

    // A real Maxwell branch can release resolvable heat while its rounded
    // spring-strain endpoint remains unchanged. Do not book the ideal storage
    // decrement into an accepted state which did not actually lose it.
    let mesh = mesh();
    let mut supported = plate.clone();
    supported.m = matrix(1, &[1.0]);
    supported.k = matrix(1, &[1.0]);
    let branch = PlateRelaxationBranch::assemble(
        &mesh,
        &supported,
        vec![[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
        0.125 / 1e-20,
        PlateViscosityBudget {
            max_elements: 1,
            max_dofs: 1,
            max_full_dofs: 9,
            max_assembly_entries: 81,
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(branch.excess_stiffness().get(0, 0), 3.0);
    let owner = RelaxingPlate::prepare(
        &supported,
        vec![branch],
        0.125,
        1e-12,
        budget(&supported, 1),
        &LdltOptions::default(),
        &mut || false,
    )
    .unwrap();
    let mut state = owner.state(&[0.25], &[0.0], &[0.125]).unwrap();
    let record = state.advance(&[0.625], &mut || false).unwrap();
    assert_eq!(state.displacement(), &[0.25]);
    assert_eq!(state.velocity(), &[0.0]);
    assert_eq!(state.branch_strains(), &[0.125]);
    assert_eq!(record.stored_change_j, 0.0);
    assert_eq!(record.kinetic_change_j, 0.0);
    assert!(record.dissipated_j > 0.0);
    assert_eq!(record.balance_residual_j, record.dissipated_j);

    // The intended velocity increment can straddle a binary exponent boundary
    // and round on publication, even when the midpoint solve is exact.
    let mut free = plate.clone();
    free.m = matrix(1, &[1.0]);
    free.k = matrix(1, &[0.0]);
    let owner = RelaxingPlate::prepare(
        &free,
        vec![],
        1.0,
        1e-12,
        budget(&free, 0),
        &LdltOptions::default(),
        &mut || false,
    )
    .unwrap();
    let initial_velocity = 1.0 - 0.5 * f64::EPSILON;
    let target_midpoint = 1.0 + f64::EPSILON;
    let force = 2.0 * (target_midpoint - initial_velocity);
    let mut state = owner.state(&[0.0], &[initial_velocity], &[]).unwrap();
    let record = state.advance(&[force], &mut || false).unwrap();
    let velocity = state.velocity()[0];
    let kinetic = 0.5 * (velocity - initial_velocity) * (velocity + initial_velocity);
    assert!((record.kinetic_change_j - kinetic).abs() < 1e-30);
    assert!((kinetic - record.supplied_work_j).abs() > 1e-17);
    assert!((record.balance_residual_j - (kinetic - record.supplied_work_j)).abs() < 1e-30);
}
