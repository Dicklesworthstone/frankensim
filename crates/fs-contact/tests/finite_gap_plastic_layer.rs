//! Physical integration oracles for normal indentation plus independent layer shear.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticDrive, FiniteGapPlasticLayerContact,
    FiniteGapPlasticLayerError, FiniteGapPlasticLayerShearAdvance, FiniteGapPlasticLayerSite,
    FiniteGapPlasticLayerState, FiniteGapPlasticLoadControls, FiniteGapPlasticModel,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_tribo::layer_contact::tensor::{ShearStiffness, ShearViscosity, TensorLayerRelaxationArm};
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerState};
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, MultipleTangentialError, PreparedTangentialContacts,
};
use fs_tribo::{TangentialLaw, TangentialSolveConfig};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 83,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn normal(cells: &[usize], cx: &Cx<'_>) -> FiniteGapPlasticModel {
    let mut gap = vec![0.01; 25];
    for &cell in cells {
        gap[cell] = 0.0;
    }
    FiniteGapPlasticModel::new(
        FiniteGapPlasticConfig {
            grid: FiniteGapGrid {
                cells_x: 5,
                cells_y: 5,
                cell_width_m: 1e-3,
                cell_depth_m: 1e-3,
            },
            undeformed_gap_m: gap,
            reduced_modulus_pa: 1e9,
            hardness_pa: vec![1e6; 25],
            maximum_active_set_iterations: 512,
            complementarity_tolerance_m: 1e-12,
            boundary_clearance_cells: 1,
        },
        cx,
    )
    .unwrap()
}

fn response(count: usize, matrix: Vec<f64>) -> PreparedTangentialContacts {
    PreparedTangentialContacts::new(
        count,
        matrix,
        MultipleTangentialConfig {
            maximum_contacts: 8,
            maximum_matrix_terms: 256,
            maximum_sweeps: 256,
            maximum_law_evaluations: 65536,
            local: MatrixTangentialConfig {
                maximum_evaluations: 512,
                displacement_tolerance_m: 2e-13,
                energy_tolerance_j: 2e-12,
                matrix_relative_tolerance: 1e-12,
            },
            total_energy_tolerance_j: 4e-12,
        },
        || false,
    )
    .unwrap()
}

fn controls() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 128,
        maximum_force_n: 1000.0,
        force_absolute_tolerance_n: 1e-13,
        force_relative_tolerance: 1e-12,
    }
}

fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1e-5,
        maximum_approach_m: 1e-4,
        absolute_force_tolerance_n: 1e-11,
        relative_force_tolerance: 1e-10,
        approach_tolerance_m: 1e-14,
        maximum_equilibrium_evaluations: 128,
    }
}

fn elastic(k: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(k, 100.0, 1.0, 300.0, 0.02).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

// Complete Gaussian elimination is independent of the native block/constitutive solve.
fn solve(mut matrix: Vec<f64>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| matrix[i * n + k].abs().total_cmp(&matrix[j * n + k].abs()))
            .unwrap();
        for j in 0..n {
            matrix.swap(k * n + j, pivot * n + j);
        }
        rhs.swap(k, pivot);
        for i in k + 1..n {
            let factor = matrix[i * n + k] / matrix[k * n + k];
            for j in k..n {
                matrix[i * n + j] -= factor * matrix[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] -= matrix[i * n + j] * rhs[j];
        }
        rhs[i] /= matrix[i * n + i];
    }
    rhs
}

#[test]
fn full_cell_recoil_matches_independent_elastic_solution_and_patch_wrench_work() {
    with_cx(|cx, _| {
        let matrix = vec![
            0.0004, 0.00005, 0.00012, -0.00003, 0.00005, 0.0003, 0.00002, 0.00008, 0.00012,
            0.00002, 0.0005, -0.00004, -0.00003, 0.00008, -0.00004, 0.00035,
        ];
        let stiffness = [1000.0, 1500.0];
        let model = FiniteGapPlasticLayerContact::new(
            normal(&[11, 13], cx),
            vec![
                FiniteGapPlasticLayerSite {
                    cell_index: 11,
                    law: elastic(stiffness[0]),
                },
                FiniteGapPlasticLayerSite {
                    cell_index: 13,
                    law: elastic(stiffness[1]),
                },
            ],
            response(2, matrix.clone()),
            cx,
        )
        .unwrap();
        let loaded = model
            .advance_normal(&model.initial_state(), 1e-5, cx)
            .unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        // A rigid patch translation plus spin about the centre, x = +/- 1 mm.
        let translation = [0.002, 0.001];
        let angle = 0.4;
        let free = [
            [translation[0], translation[1] - 1e-3 * angle],
            [translation[0], translation[1] + 1e-3 * angle],
        ];
        let step = model
            .advance_shear(&before, &free, 0.002, controls(), cx)
            .unwrap();
        let released: Vec<_> = free.into_iter().flatten().collect();
        let mut system = matrix.clone();
        for i in 0..4 {
            system[4 * i + i] += 2.0 / stiffness[i / 2];
        }
        let expected = solve(system.clone(), released.clone());
        let mut diagonalized = system;
        for i in 0..4 {
            for j in 0..4 {
                if i / 2 != j / 2 {
                    diagonalized[4 * i + j] = 0.0;
                }
            }
        }
        let wrong = solve(diagonalized, released.clone());
        assert!(
            expected
                .iter()
                .zip(wrong)
                .any(|(a, b)| (a - b).abs() > 1e-3)
        );
        let actual: Vec<_> = step
            .shear
            .contacts
            .iter()
            .flat_map(|s| s.traction_n)
            .collect();
        let mut work = 0.0;
        let mut storage = 0.0;
        for i in 0..4 {
            close(actual[i], expected[i], 2e-9);
            let recoil: f64 = (0..4).map(|j| matrix[4 * i + j] * actual[j]).sum();
            let slip = released[i] - recoil;
            close(
                step.shear.contacts[i / 2].slip_increment_m[i % 2],
                slip,
                3e-13,
            );
            work += actual[i] * slip;
            storage += 2.0 * actual[i] * actual[i] / stiffness[i / 2];
        }
        close(step.shear.stored_change_j, storage, 4e-12);
        close(step.shear.stored_change_j, work, 4e-12);
        assert_eq!(step.shear.dissipated_j, 0.0);
        assert_eq!(step.pressure_pa, vec![1e6, 1e6]);
        close(step.shear_force_n[0], actual[0] + actual[2], 1e-15);
        close(step.shear_force_n[1], actual[1] + actual[3], 1e-15);
        close(step.spin_moment_n_m, 1e-3 * (actual[3] - actual[1]), 1e-17);
        let released_work: f64 = actual.iter().zip(&released).map(|(f, u)| f * u).sum();
        close(
            released_work,
            step.shear_force_n[0] * translation[0]
                + step.shear_force_n[1] * translation[1]
                + step.spin_moment_n_m * angle,
            1e-17,
        );
        assert_eq!(step.candidate.normal_state(), before.normal_state());
        assert_eq!(before, saved);
    });
}

fn memory_layer() -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite([[1000.0, 70.0], [70.0, 700.0]]).unwrap(),
        vec![TensorLayerRelaxationArm {
            stiffness: ShearStiffness::from_gram_rows(vec![[10.0, 2.0], [0.0, 5.0]]).unwrap(),
            relaxation_time_s: 0.02,
        }],
        0.1,
        0.4,
        300.0,
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity(ShearViscosity::from_gram_rows(vec![[0.2, -0.1]]).unwrap())
    .unwrap()
}

fn independent_shear_energy(state: &ThirdBodyLayerState) -> f64 {
    let e = state.elastic_slip_m();
    let r = state.relaxation_slip_m()[0];
    0.5 * (1000.0 * e[0] * e[0] + 140.0 * e[0] * e[1] + 700.0 * e[1] * e[1])
        + 0.5 * ((10.0 * r[0] + 2.0 * r[1]).powi(2) + (5.0 * r[1]).powi(2))
}

fn check_memory_step(
    before: &FiniteGapPlasticLayerState,
    step: &FiniteGapPlasticLayerShearAdvance,
    released: [f64; 2],
    matrix: &[f64],
    dt: f64,
) -> (f64, f64) {
    let old = &before.layer_states()[0];
    let new = &step.candidate.layer_states()[0];
    let dp = [
        new.plastic_slip_m()[0] - old.plastic_slip_m()[0],
        new.plastic_slip_m()[1] - old.plastic_slip_m()[1],
    ];
    let distance = dp[0].hypot(dp[1]);
    let plastic = (0.1 + 0.4 * distance / dt) * distance;
    let r = [
        0.5 * (old.relaxation_slip_m()[0][0] + new.relaxation_slip_m()[0][0]),
        0.5 * (old.relaxation_slip_m()[0][1] + new.relaxation_slip_m()[0][1]),
    ];
    let relaxation = dt / 0.02 * ((10.0 * r[0] + 2.0 * r[1]).powi(2) + (5.0 * r[1]).powi(2));
    let delta = [
        new.elastic_slip_m()[0] - old.elastic_slip_m()[0],
        new.elastic_slip_m()[1] - old.elastic_slip_m()[1],
    ];
    let kelvin = (0.2 * delta[0] - 0.1 * delta[1]).powi(2) / dt;
    close(step.dissipation[0].plastic_j, plastic, 2e-16);
    close(step.dissipation[0].relaxation_j, relaxation, 2e-16);
    close(step.dissipation[0].kelvin_j, kelvin, 2e-16);
    let heat = plastic + relaxation + kelvin;
    close(step.shear.dissipated_j, heat, 2e-16);
    let force = step.shear.contacts[0].traction_n;
    let actual = [
        released[0] - matrix[0] * force[0] - matrix[1] * force[1],
        released[1] - matrix[2] * force[0] - matrix[3] * force[1],
    ];
    let work = force[0] * actual[0] + force[1] * actual[1];
    close(
        independent_shear_energy(new) - independent_shear_energy(old) + heat,
        work,
        4e-12,
    );
    assert_eq!(step.candidate.normal_state(), before.normal_state());
    (work, heat)
}

#[test]
fn indentation_unload_and_shear_reversal_keep_both_plastic_histories_and_total_work() {
    with_cx(|cx, _| {
        // Exact self-compliance integral of the one loaded square cell.
        let compliance = 4.0 * 1e-3 * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * 1e9);
        let yield_approach = compliance * 1e6;
        let first_dent = 3e-6;
        let normal = normal(&[12], cx);
        let original = normal
            .advance(&normal.initial_state(), yield_approach + first_dent, cx)
            .unwrap();
        let matrix = vec![0.0002, 0.00003, 0.00003, 0.0003];
        let model = FiniteGapPlasticLayerContact::new(
            normal,
            vec![FiniteGapPlasticLayerSite {
                cell_index: 12,
                law: memory_layer(),
            }],
            response(1, matrix.clone()),
            cx,
        )
        .unwrap();
        let mut state = model.state_from_normal(original.candidate, cx).unwrap();
        close(
            state.normal_state().plastic_depth_m()[12],
            first_dent,
            1e-17,
        );
        close(
            state.normal_state().cumulative_plastic_work_j(),
            first_dent,
            1e-17,
        );
        let mut normal_work = 0.5 * yield_approach + first_dent; // H*A = 1 N.
        let (mut shear_work, mut shear_heat) = (0.0, 0.0);
        let first_free = [0.003, -0.001];
        let first = model
            .advance_shear(&state, &[first_free], 0.002, controls(), cx)
            .unwrap();
        let (work, heat) = check_memory_step(&state, &first, first_free, &matrix, 0.002);
        shear_work += work;
        shear_heat += heat;
        let forward_plastic = first.candidate.layer_states()[0].plastic_slip_m();
        assert!(first.dissipation[0].plastic_j > 0.0);
        assert!(first.dissipation[0].relaxation_j > 0.0);
        assert!(first.dissipation[0].kelvin_j > 0.0);
        state = first.candidate;

        let unloaded = model
            .advance_normal(&state, first_dent + 0.5 * yield_approach, cx)
            .unwrap();
        assert_eq!(unloaded.candidate.layer_states(), state.layer_states());
        assert_eq!(unloaded.normal.plastic_work_j, 0.0);
        close(unloaded.normal.normal.normal_force_n, 0.5, 1e-10);
        normal_work -= 0.5 * (1.0 + 0.5) * 0.5 * yield_approach;
        state = unloaded.candidate;
        let reverse_free = [-0.004, 0.002];
        let reverse = model
            .advance_shear(&state, &[reverse_free], 0.002, controls(), cx)
            .unwrap();
        let (work, heat) = check_memory_step(&state, &reverse, reverse_free, &matrix, 0.002);
        shear_work += work;
        shear_heat += heat;
        let reversed = reverse.candidate.layer_states()[0].plastic_slip_m();
        let dp = [
            reversed[0] - forward_plastic[0],
            reversed[1] - forward_plastic[1],
        ];
        assert!(dp[0] * forward_plastic[0] + dp[1] * forward_plastic[1] < 0.0);
        close(reverse.pressure_pa[0], 0.5e6, 1e-4);
        state = reverse.candidate;

        let second_dent = 2e-6;
        let reloaded = model
            .advance_normal(&state, first_dent + yield_approach + second_dent, cx)
            .unwrap();
        assert_eq!(reloaded.candidate.layer_states(), state.layer_states());
        close(reloaded.normal.plastic_work_j, second_dent, 1e-17);
        normal_work += 0.5 * (0.5 + 1.0) * 0.5 * yield_approach + second_dent;
        state = reloaded.candidate;
        let held = model
            .advance_shear(&state, &[[0.0; 2]], 0.002, controls(), cx)
            .unwrap();
        let (work, heat) = check_memory_step(&state, &held, [0.0; 2], &matrix, 0.002);
        shear_work += work;
        shear_heat += heat;
        state = held.candidate;
        close(
            state.normal_state().plastic_depth_m()[12],
            first_dent + second_dent,
            1e-17,
        );
        let shear_storage = independent_shear_energy(&state.layer_states()[0]);
        close(
            model.sites()[0]
                .law
                .stored_energy_j(&state.layer_states()[0]),
            shear_storage,
            2e-16,
        );
        let stored = 0.5 * yield_approach + shear_storage;
        let normal_plastic_work = state.normal_state().cumulative_plastic_work_j();
        close(
            stored + normal_plastic_work + shear_heat,
            normal_work + shear_work,
            2e-11,
        );
    });
}

#[test]
fn force_and_spring_normal_boundaries_keep_the_same_layer_memory() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticLayerContact::new(
            normal(&[12], cx),
            vec![FiniteGapPlasticLayerSite {
                cell_index: 12,
                law: elastic(1000.0),
            }],
            response(1, vec![0.0; 4]),
            cx,
        )
        .unwrap();
        let load = model
            .advance_normal_at_load(&model.initial_state(), 0.25, load_controls(), cx)
            .unwrap();
        close(load.boundary.achieved_force_n, 0.25, 1e-10);
        let shear = model
            .advance_shear(&load.candidate, &[[0.0001, 0.0]], 0.002, controls(), cx)
            .unwrap();
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1e5,
            position_m: 2e-5,
        };
        let driven = model
            .advance_normal_with_drive(&shear.candidate, drive, load_controls(), cx)
            .unwrap();
        close(driven.boundary.achieved_force_n, 1.0, 1e-10);
        close(
            driven.boundary.approach_m,
            drive.position_m - 1.0 / drive.stiffness_n_per_m,
            1e-15,
        );
        close(
            driven.boundary.spring_energy_j,
            0.5 / drive.stiffness_n_per_m,
            1e-15,
        );
        assert!(driven.normal.plastic_work_j > 0.0);
        assert_eq!(
            driven.candidate.layer_states(),
            shear.candidate.layer_states()
        );
        let dent = driven.candidate.normal_state().plastic_depth_m()[12];
        let unloaded = model
            .advance_normal_at_load(&driven.candidate, 0.5, load_controls(), cx)
            .unwrap();
        close(unloaded.boundary.achieved_force_n, 0.5, 1e-10);
        close(
            unloaded.candidate.normal_state().plastic_depth_m()[12],
            dent,
            1e-17,
        );
        assert_eq!(
            unloaded.candidate.layer_states(),
            shear.candidate.layer_states()
        );
        assert_eq!(unloaded.normal.plastic_work_j, 0.0);
    });
}

#[test]
fn opening_foreign_registration_late_cell_failure_and_cancellation_are_atomic() {
    with_cx(|cx, gate| {
        let sites = vec![
            FiniteGapPlasticLayerSite {
                cell_index: 11,
                law: elastic(1000.0),
            },
            FiniteGapPlasticLayerSite {
                cell_index: 13,
                law: elastic(1000.0),
            },
        ];
        let model = FiniteGapPlasticLayerContact::new(
            normal(&[11, 13], cx),
            sites.clone(),
            response(2, vec![0.0; 16]),
            cx,
        )
        .unwrap();
        assert!(matches!(
            model.advance_shear(
                &model.initial_state(),
                &[[0.0; 2]; 2],
                0.002,
                controls(),
                cx
            ),
            Err(FiniteGapPlasticLayerError::DetachedLayer { cell_index: 11 })
        ));
        let loaded = model
            .advance_normal(&model.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let saved = loaded.clone();
        let mut limited = controls();
        limited.maximum_force_n = 1.0;
        assert!(matches!(
            model.advance_shear(&loaded, &[[0.001, 0.0], [0.01, 0.0]], 0.002, limited, cx),
            Err(FiniteGapPlasticLayerError::Shear(
                MultipleTangentialError::Contact { index: 1, .. }
            ))
        ));
        assert_eq!(loaded, saved);
        let sheared = model
            .advance_shear(
                &loaded,
                &[[0.0001, 0.0], [0.001, 0.0]],
                0.002,
                controls(),
                cx,
            )
            .unwrap()
            .candidate;
        let saved = sheared.clone();
        assert!(matches!(
            model.advance_normal(&sheared, -1e-6, cx),
            Err(FiniteGapPlasticLayerError::OpeningWithStoredShear { cell_index: 11, .. })
        ));
        assert_eq!(sheared, saved);
        let mut foreign_sites = sites.clone();
        foreign_sites[0].cell_index = 12;
        let foreign = FiniteGapPlasticLayerContact::new(
            model.normal_model().clone(),
            foreign_sites,
            response(2, vec![0.0; 16]),
            cx,
        )
        .unwrap();
        assert!(matches!(
            foreign.advance_shear(&sheared, &[[0.0; 2]; 2], 0.002, controls(), cx),
            Err(FiniteGapPlasticLayerError::Invalid(_))
        ));
        assert!(
            model
                .advance_shear(&sheared, &[], 0.002, controls(), cx)
                .is_err()
        );
        assert!(
            model
                .advance_shear(&sheared, &[[0.0; 2]; 2], 0.0, controls(), cx)
                .is_err()
        );
        let mut reversed = sites;
        reversed.reverse();
        assert!(
            FiniteGapPlasticLayerContact::new(
                model.normal_model().clone(),
                reversed,
                response(2, vec![0.0; 16]),
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(
            model
                .advance_shear(&sheared, &[[0.0; 2]; 2], 0.002, controls(), cx)
                .is_err()
        );
        assert!(model.advance_normal(&sheared, 1e-5, cx).is_err());
        assert_eq!(sheared, saved);
    });
}

#[test]
fn opening_and_recontact_preserve_permanent_shear_when_storage_is_zero() {
    with_cx(|cx, _| {
        // Binary coefficients make the two native midpoint endpoints exact:
        // e_1=epsilon, dp_1=3*epsilon/8; e_2=0, dp_2=3*epsilon/8.
        let epsilon = 1.0 / 1048576.0;
        let layer = ThirdBodyLayerContact::new(1.0 / epsilon, 0.125, 1.0, 300.0, 1e-4).unwrap();
        let model = FiniteGapPlasticLayerContact::new(
            normal(&[12], cx),
            vec![FiniteGapPlasticLayerSite {
                cell_index: 12,
                law: layer,
            }],
            response(1, vec![0.0; 4]),
            cx,
        )
        .unwrap();
        let loaded = model
            .advance_normal(&model.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let first = model
            .advance_shear(
                &loaded,
                &[[11.0 * epsilon / 8.0, 0.0]],
                epsilon,
                controls(),
                cx,
            )
            .unwrap();
        close(
            first.candidate.layer_states()[0].elastic_slip_m()[0],
            epsilon,
            1e-20,
        );
        let second = model
            .advance_shear(
                &first.candidate,
                &[[-5.0 * epsilon / 8.0, 0.0]],
                epsilon,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            second.candidate.layer_states()[0].elastic_slip_m(),
            [0.0; 2]
        );
        assert_eq!(
            second.candidate.layer_states()[0].plastic_slip_m(),
            [0.75 * epsilon, 0.0]
        );
        let opened = model.advance_normal(&second.candidate, -1e-6, cx).unwrap();
        assert_eq!(opened.normal.normal.normal_force_n, 0.0);
        assert_eq!(
            opened.candidate.layer_states(),
            second.candidate.layer_states()
        );
        assert_eq!(
            opened.candidate.normal_state().plastic_depth_m(),
            second.candidate.normal_state().plastic_depth_m()
        );
        let reloaded = model.advance_normal(&opened.candidate, 1e-5, cx).unwrap();
        assert_eq!(
            reloaded.candidate.layer_states(),
            second.candidate.layer_states()
        );
        assert_eq!(reloaded.normal.plastic_work_j, 0.0);
    });
}

fn coefficient_successor(
    before: &FiniteGapPlasticLayerContact,
    modulus: f64,
    hardness: f64,
    laws: Vec<ThirdBodyLayerContact>,
    cx: &Cx<'_>,
) -> FiniteGapPlasticLayerContact {
    assert_eq!(laws.len(), before.sites().len());
    let normal = before
        .normal_model()
        .with_material_properties(
            modulus,
            vec![hardness; before.normal_model().config().hardness_pa.len()],
            cx,
        )
        .unwrap();
    let sites = before
        .sites()
        .iter()
        .zip(laws)
        .map(|(site, law)| FiniteGapPlasticLayerSite {
            cell_index: site.cell_index,
            law,
        })
        .collect();
    FiniteGapPlasticLayerContact::new(normal, sites, before.mechanical_response().clone(), cx)
        .unwrap()
}

fn assert_retained_layer_strains(
    before: &FiniteGapPlasticLayerState,
    after: &FiniteGapPlasticLayerState,
) {
    assert_eq!(before.layer_states().len(), after.layer_states().len());
    for (old, new) in before.layer_states().iter().zip(after.layer_states()) {
        assert_eq!(old.elastic_slip_m(), new.elastic_slip_m());
        assert_eq!(old.plastic_slip_m(), new.plastic_slip_m());
        assert_eq!(old.relaxation_slip_m(), new.relaxation_slip_m());
    }
}

// Independent quadratic potentials. Kelvin and kinetic coefficients store no energy.
fn scaled_tensor_energy(state: &ThirdBodyLayerState, equilibrium: f64, arm: f64) -> f64 {
    let e = state.elastic_slip_m();
    let r = state.relaxation_slip_m()[0];
    equilibrium * 0.5 * (1000.0 * e[0].powi(2) + 140.0 * e[0] * e[1] + 700.0 * e[1].powi(2))
        + arm * 0.5 * ((10.0 * r[0] + 2.0 * r[1]).powi(2) + (5.0 * r[1]).powi(2))
}

#[test]
fn coefficient_softening_and_cooling_preserve_memory_and_account_signed_work() {
    with_cx(|cx, _| {
        let compliance = 4.0 * 1e-3 * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * 1e9);
        let yield_approach = compliance * 1e6;
        let matrix = vec![0.0002, 0.00003, 0.00003, 0.0003];
        let cold = FiniteGapPlasticLayerContact::new(
            normal(&[12], cx),
            vec![FiniteGapPlasticLayerSite {
                cell_index: 12,
                law: memory_layer(),
            }],
            response(1, matrix.clone()),
            cx,
        )
        .unwrap();
        let loaded = cold
            .advance_normal(&cold.initial_state(), yield_approach + 3e-6, cx)
            .unwrap();
        let before = cold
            .advance_shear(&loaded.candidate, &[[0.003, -0.001]], 0.002, controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        assert!(before.layer_states()[0].plastic_slip_m()[0] > 0.0);
        assert!(before.layer_states()[0].relaxation_slip_m()[0][0] > 0.0);
        let scale = 0.75_f64.sqrt();
        let hot_law = ThirdBodyLayerContact::new_tensor(
            ShearStiffness::positive_definite([[500.0, 35.0], [35.0, 350.0]]).unwrap(),
            vec![TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![
                    [10.0 * scale, 2.0 * scale],
                    [0.0, 5.0 * scale],
                ])
                .unwrap(),
                relaxation_time_s: 0.05,
            }],
            0.04,
            0.2,
            420.0,
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity(ShearViscosity::from_gram_rows(vec![[0.3, -0.05]]).unwrap())
        .unwrap();
        let hot = coefficient_successor(&cold, 1e9, 0.5e6, vec![hot_law], cx);
        let refreshed = hot.refresh_from(&cold, &before, cx).unwrap();
        assert_retained_layer_strains(&before, &refreshed.candidate);
        assert_eq!(before, saved);
        close(
            refreshed.candidate.normal_state().pressure_pa()[12],
            0.5e6,
            1e-7,
        );
        close(
            refreshed.candidate.normal_state().plastic_depth_m()[12],
            3e-6 + 0.5 * yield_approach,
            1e-17,
        );
        assert_eq!(
            refreshed.candidate.normal_state().recession_m(),
            before.normal_state().recession_m()
        );
        let old_normal_energy = 0.5 * yield_approach; // H*A = 1 N.
        close(
            refreshed.normal.elastic_energy_change_j,
            -0.75 * old_normal_energy,
            1e-17,
        );
        close(
            refreshed.normal.response.plastic_work_j,
            0.5 * old_normal_energy,
            1e-17,
        );
        close(
            refreshed.normal.parameter_work_j,
            -0.25 * old_normal_energy,
            1e-17,
        );
        assert_eq!(refreshed.normal.normal_work_j, 0.0);
        let old_layer_energy = scaled_tensor_energy(&before.layer_states()[0], 1.0, 1.0);
        let hot_layer_energy = scaled_tensor_energy(&before.layer_states()[0], 0.5, 0.75);
        close(
            refreshed.layer_coefficient_free_energy_change_j[0],
            hot_layer_energy - old_layer_energy,
            2e-18,
        );
        close(
            refreshed.total_layer_coefficient_free_energy_change_j,
            hot_layer_energy - old_layer_energy,
            2e-18,
        );

        // Returning coefficients restores shear potential, never permanent indentation.
        let cooled = cold.refresh_from(&hot, &refreshed.candidate, cx).unwrap();
        assert_retained_layer_strains(&before, &cooled.candidate);
        assert_eq!(
            cooled.candidate.normal_state().plastic_depth_m(),
            refreshed.candidate.normal_state().plastic_depth_m()
        );
        close(
            cooled.candidate.normal_state().pressure_pa()[12],
            0.5e6,
            1e-7,
        );
        assert_eq!(cooled.normal.response.plastic_work_j, 0.0);
        close(
            cooled.total_layer_coefficient_free_energy_change_j
                + refreshed.total_layer_coefficient_free_energy_change_j,
            0.0,
            2e-18,
        );
        let unchanged = cold.refresh_from(&cold, &before, cx).unwrap();
        assert_eq!(unchanged.candidate, before);
        assert_eq!(unchanged.total_layer_coefficient_free_energy_change_j, 0.0);
        assert_eq!(unchanged.normal.parameter_work_j, 0.0);

        // Only the next positive-time interval performs recoil and dissipative evolution.
        let free = [-0.004, 0.001];
        let dt = 0.002;
        let next = hot
            .advance_shear(&refreshed.candidate, &[free], dt, controls(), cx)
            .unwrap();
        let old = &refreshed.candidate.layer_states()[0];
        let new = &next.candidate.layer_states()[0];
        let dp = [
            new.plastic_slip_m()[0] - old.plastic_slip_m()[0],
            new.plastic_slip_m()[1] - old.plastic_slip_m()[1],
        ];
        let distance = dp[0].hypot(dp[1]);
        let plastic = (0.04 + 0.2 * distance / dt) * distance;
        let r = [
            0.5 * (old.relaxation_slip_m()[0][0] + new.relaxation_slip_m()[0][0]),
            0.5 * (old.relaxation_slip_m()[0][1] + new.relaxation_slip_m()[0][1]),
        ];
        let relaxation =
            dt / 0.05 * 0.75 * ((10.0 * r[0] + 2.0 * r[1]).powi(2) + (5.0 * r[1]).powi(2));
        let de = [
            new.elastic_slip_m()[0] - old.elastic_slip_m()[0],
            new.elastic_slip_m()[1] - old.elastic_slip_m()[1],
        ];
        let kelvin = (0.3 * de[0] - 0.05 * de[1]).powi(2) / dt;
        close(next.dissipation[0].plastic_j, plastic, 2e-16);
        close(next.dissipation[0].relaxation_j, relaxation, 2e-16);
        close(next.dissipation[0].kelvin_j, kelvin, 2e-16);
        let f = next.shear.contacts[0].traction_n;
        let work = f[0] * (free[0] - matrix[0] * f[0] - matrix[1] * f[1])
            + f[1] * (free[1] - matrix[2] * f[0] - matrix[3] * f[1]);
        close(
            scaled_tensor_energy(new, 0.5, 0.75) - hot_layer_energy + plastic + relaxation + kelvin,
            work,
            4e-12,
        );
        assert!(
            cold.advance_shear(&refreshed.candidate, &[free], dt, controls(), cx)
                .is_err()
        );
    });
}

#[test]
fn coefficient_refresh_at_load_and_held_spring_use_original_plastic_history() {
    with_cx(|cx, _| {
        let compliance = 4.0 * 1e-3 * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * 1e9);
        let cold = FiniteGapPlasticLayerContact::new(
            normal(&[12], cx),
            vec![FiniteGapPlasticLayerSite {
                cell_index: 12,
                law: memory_layer(),
            }],
            response(1, vec![0.0002, 0.00003, 0.00003, 0.0003]),
            cx,
        )
        .unwrap();
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1e5,
            position_m: 2e-5,
        };
        let driven = cold
            .advance_normal_with_drive(&cold.initial_state(), drive, load_controls(), cx)
            .unwrap();
        let sheared = cold
            .advance_shear(&driven.candidate, &[[0.003, -0.001]], 0.002, controls(), cx)
            .unwrap()
            .candidate;
        let unloaded = cold
            .advance_normal_at_load(&sheared, 0.5, load_controls(), cx)
            .unwrap()
            .candidate;
        let old_energy = 0.5 * 1e-6 * compliance * 0.5e6_f64.powi(2);
        for (modulus_scale, hardness) in [(0.5, 0.8e6), (2.0, 0.75e6)] {
            let successor = coefficient_successor(
                &cold,
                modulus_scale * 1e9,
                hardness,
                vec![memory_layer()],
                cx,
            );
            let step = successor
                .refresh_from_at_load(&cold, &unloaded, 0.5, load_controls(), cx)
                .unwrap();
            assert_retained_layer_strains(&unloaded, &step.candidate);
            assert_eq!(
                step.candidate.normal_state().plastic_depth_m(),
                unloaded.normal_state().plastic_depth_m()
            );
            assert_eq!(step.normal.response.plastic_work_j, 0.0);
            close(
                step.normal.elastic_energy_change_j,
                old_energy * (1.0 / modulus_scale - 1.0),
                2e-15,
            );
            close(
                step.normal.normal_work_j,
                2.0 * old_energy * (1.0 / modulus_scale - 1.0),
                3e-15,
            );
            close(
                step.normal.parameter_work_j,
                -old_energy * (1.0 / modulus_scale - 1.0),
                3e-15,
            );
            assert_eq!(step.total_layer_coefficient_free_energy_change_j, 0.0);
            // At the old approach the stiffer model would exceed its new cap.
            // A preliminary fixed-approach refresh would manufacture a dent.
            close(step.boundary.achieved_force_n, 0.5, 1e-10);
        }
        let softened = coefficient_successor(&cold, 1e9, 0.5e6, vec![memory_layer()], cx);
        let spring = softened
            .refresh_from_with_drive(&cold, &sheared, drive, load_controls(), cx)
            .unwrap();
        assert_retained_layer_strains(&sheared, &spring.candidate);
        close(spring.boundary.achieved_force_n, 0.5, 1e-10);
        close(spring.boundary.approach_m, 1.5e-5, 1e-14);
        close(spring.boundary.spring_energy_j, 1.25e-6, 1e-15);
        close(spring.normal.normal_work_j, 3.75e-6, 1e-15);
        let plastic_work = 0.5 * (5e-6 + 0.5 * compliance * 1e6);
        close(spring.normal.response.plastic_work_j, plastic_work, 1e-14);
        let storage_change = -0.75 * 0.5 * compliance * 1e6;
        close(spring.normal.elastic_energy_change_j, storage_change, 1e-15);
        close(
            spring.normal.parameter_work_j,
            storage_change + plastic_work - 3.75e-6,
            1e-14,
        );
        assert_eq!(spring.total_layer_coefficient_free_energy_change_j, 0.0);
    });
}

#[test]
fn coefficient_refresh_refuses_incompatible_memory_and_opening_without_partial_adoption() {
    with_cx(|cx, gate| {
        let mut config = normal(&[12, 13], cx).config().clone();
        config.undeformed_gap_m[13] = 2e-7;
        let cold = FiniteGapPlasticLayerContact::new(
            FiniteGapPlasticModel::new(config, cx).unwrap(),
            vec![
                FiniteGapPlasticLayerSite {
                    cell_index: 12,
                    law: memory_layer(),
                },
                FiniteGapPlasticLayerSite {
                    cell_index: 13,
                    law: memory_layer(),
                },
            ],
            response(2, vec![0.0; 16]),
            cx,
        )
        .unwrap();
        let loaded = cold
            .advance_normal(&cold.initial_state(), 5e-7, cx)
            .unwrap();
        let force = loaded.normal.normal.normal_force_n;
        assert!(loaded.candidate.normal_state().pressure_pa()[13] > 0.0);
        let before = cold
            .advance_shear(
                &loaded.candidate,
                &[[0.0, 0.0], [0.0001, 0.0]],
                0.002,
                controls(),
                cx,
            )
            .unwrap()
            .candidate;
        let saved = before.clone();
        let incompatible =
            coefficient_successor(&cold, 1e9, 1e6, vec![memory_layer(), elastic(500.0)], cx);
        assert!(matches!(
            incompatible.refresh_from(&cold, &before, cx),
            Err(FiniteGapPlasticLayerError::Cell { cell_index: 13, .. })
        ));
        assert_eq!(before, saved);

        // Stiffening at held resultant narrows the footprint: the neighbouring
        // shear-bearing cell opens. Its stored shear cannot be discarded as heat.
        let stiff =
            coefficient_successor(&cold, 1e10, 1e6, vec![memory_layer(), memory_layer()], cx);
        let normal_only = stiff
            .normal_model()
            .refresh_from_at_load(
                cold.normal_model(),
                before.normal_state(),
                force,
                load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(normal_only.continuation.candidate.pressure_pa()[13], 0.0);
        assert!(
            matches!(stiff.refresh_from_at_load(&cold, &before, force, load_controls(), cx), Err(FiniteGapPlasticLayerError::OpeningWithStoredShear { cell_index: 13, stored_energy_j }) if stored_energy_j > 0.0)
        );
        assert_eq!(before, saved);

        let changed_response = FiniteGapPlasticLayerContact::new(
            cold.normal_model().clone(),
            cold.sites().to_vec(),
            response(
                2,
                vec![
                    1e-4, 0.0, 0.0, 0.0, 0.0, 1e-4, 0.0, 0.0, 0.0, 0.0, 1e-4, 0.0, 0.0, 0.0, 0.0,
                    1e-4,
                ],
            ),
            cx,
        )
        .unwrap();
        assert!(matches!(
            changed_response.refresh_from(&cold, &before, cx),
            Err(FiniteGapPlasticLayerError::Invalid(_))
        ));
        gate.request();
        assert!(cold.refresh_from(&cold, &before, cx).is_err());
        assert_eq!(before, saved);
    });
}

#[path = "finite_gap_plastic_layer/drive.rs"]
mod drive;
