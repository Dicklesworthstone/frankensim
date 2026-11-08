//! Progressive native shell bending and motion. Coefficients are synthetic.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::{
    DuctileJ2, DuctileJ2Conditions, DuctileJ2Law, DuctileJ2Settings,
};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};
use fs_plate::shell::{DRILLING_ALPHA, ShellMesh};
use fs_solid::IncrementSettings;
use fs_solid::plastic_shell::drilling::PlasticShellReferenceDrilling;
use fs_solid::plastic_shell::ductile::{DuctileShellEnergy, DuctileShellEnergyTolerances};
use fs_solid::plastic_shell::dynamics::{PlasticShellDensity, PlasticShellReferenceInertia};
use fs_solid::plastic_shell::{PlasticShellBudget, PlasticShellSections, SmallStrainShellProblem};
use fs_solid::plastic3::SmallStrainTetError;

const WIDTH: f64 = 0.1;
const GP: f64 = 0.00015;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 137,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}

fn model(damage_limit: f64) -> DuctileJ2<TabulatedDuctileLocus> {
    DuctileJ2::new(
        J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(), 1.0, 0.0).unwrap(),
        TabulatedDuctileLocus::new(
            vec![(-2.0, 0.0002), (2.0, 0.0002)],
            DuctileDomain {
                triaxiality: [-2.0, 2.0],
                strain_rate_s_inv: [1e-10, 10.0],
                temperature_k: [300.0, 300.0],
            },
        )
        .unwrap(),
        DuctileJ2Settings {
            characteristic_length_m: WIDTH,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: damage_limit,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-10,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn section(
    model: &DuctileJ2<TabulatedDuctileLocus>,
    durations: [f64; 2],
) -> LayeredPlaneStressSection<DuctileJ2Law<'_, TabulatedDuctileLocus>> {
    LayeredPlaneStressSection::new(
        [(-0.05, 0.0), (0.0, 0.05)]
            .into_iter()
            .enumerate()
            .map(|(i, (bottom_z_m, top_z_m))| PlaneStressLayer {
                bottom_z_m,
                top_z_m,
                subdivisions: 1,
                law: model
                    .small_strain_law(DuctileJ2Conditions {
                        duration_s: durations[i],
                        temperature_k: 300.0,
                    })
                    .unwrap(),
            })
            .collect(),
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 64,
            stress_tolerance_pa: 1e-11,
            strain_limit: 0.05,
        },
        SectionBudget {
            maximum_layers: 2,
            maximum_quadrature_points: 6,
        },
    )
    .unwrap()
}

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

/// Quadratic w=t*x*y, theta_x=t*x, theta_y=-t*y: a pure twisting
/// curvature kxy=-2*t, with zero membrane strain at all three DKT points.
fn twist(mesh: &ShellMesh, value: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            [
                0.0,
                0.0,
                value * p[0] * p[1],
                value * p[0],
                -value * p[1],
                0.0,
            ]
            .into_iter()
            .enumerate()
            .map(move |(c, value)| (6 * node + c, value))
        })
        .collect()
}

fn problem<'a, 'm>(
    mesh: &'a ShellMesh,
    section: &'a LayeredPlaneStressSection<DuctileJ2Law<'m, TabulatedDuctileLocus>>,
    prescribed: &'a [(usize, f64)],
    loads: &'a [(usize, f64)],
) -> SmallStrainShellProblem<'a, DuctileJ2Law<'m, TabulatedDuctileLocus>> {
    SmallStrainShellProblem {
        mesh,
        sections: PlasticShellSections::Uniform(section),
        coordinate_length_m: 1.0,
        rotation_limit_rad: 0.05,
        drilling_alpha: DRILLING_ALPHA,
        prescribed,
        nodal_loads: loads,
        load_path: &[1.0],
        budget: PlasticShellBudget::default(),
        settings: IncrementSettings {
            max_newton_iterations: 64,
            ..IncrementSettings::default()
        },
    }
}

fn tolerances() -> DuctileShellEnergyTolerances {
    DuctileShellEnergyTolerances {
        section_work_j_m2: 1e-11,
        shell_work_j: 1e-10,
    }
}

fn check_energy(energy: &DuctileShellEnergy) {
    near(energy.equilibrium_work_residual_j, 0.0, 1e-10);
    near(energy.balance_residual_j, 0.0, 1e-10);
    near(
        energy.material_work_j,
        energy.stored_change_j
            + energy.plastic_dissipation_j
            + energy.damage_release_j
            + energy.material_numerical_dissipation_j,
        1e-10,
    );
    near(
        energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            + energy.transverse_residual_work_j,
        energy.stored_change_j
            + energy.plastic_dissipation_j
            + energy.damage_release_j
            + energy.material_numerical_dissipation_j
            + energy.drilling.final_j
            - energy.drilling.initial_j
            + energy.drilling.numerical_loss_j
            + energy.inertial_work_j,
        1e-10,
    );
}

#[test]
fn native_bending_softens_each_layer_with_zero_membrane_resultant_and_separate_work() {
    with_cx(|_, cx| {
        let model = model(0.95);
        let section = section(&model, [1.0; 2]);
        let mesh = mesh();
        let low = twist(&mesh, 0.01);
        let p = problem(&mesh, &section, &low, &[]);
        let reference = p.prepare_reference_drilling(cx).unwrap();
        let elastic = p
            .solve_ductile_with_energy(None, &reference, tolerances(), cx)
            .unwrap();
        assert_eq!(elastic.solution, p.solve(cx).unwrap());
        assert_eq!(elastic.energy.damage_release_j, 0.0);
        let prescribed = twist(&mesh, 0.06);
        let loaded = problem(&mesh, &section, &prescribed, &[]);
        let accepted = loaded
            .solve_ductile_with_energy(Some(&elastic.solution), &reference, tolerances(), cx)
            .unwrap();
        let direct = loaded
            .solve_from_with_reference_drilling(
                &elastic.solution.coordinates,
                &elastic.solution.section_states,
                &reference,
                cx,
            )
            .unwrap();
        assert_eq!(accepted.solution, direct.solution);
        assert_eq!(accepted.energy.drilling, direct.drilling);
        assert_eq!(accepted.duration_s, 1.0);
        check_energy(&accepted.energy);
        assert!(accepted.energy.plastic_dissipation_j > 0.0);
        assert!(accepted.energy.damage_release_j > 0.0);
        assert!(accepted.energy.material_numerical_dissipation_j > 0.0);
        let mut damaged = 0;
        let mut diagnostic = 0.0;
        for (surface, state) in accepted.solution.section_states.iter().enumerate() {
            near(state.generalized_strain[5], -0.12, 1e-14);
            for force in &accepted.solution.section_resultants[surface][..3] {
                near(*force, 0.0, 1e-13);
            }
            for (point, native) in section.points().iter().zip(&state.points) {
                let shear = 0.5 * point.z_m * state.generalized_strain[5];
                let stress = native.history.integrity()
                    * 1000.0
                    * (native.strain[3] - native.history.native_state().plastic_strain[3]);
                assert!(stress * point.z_m < 0.0);
                near(native.strain[3], shear, 1e-15);
                damaged += usize::from(native.history.damage() > 0.0);
                diagnostic += (0.06 / 3.0) * point.weight_m * GP * native.history.damage() / WIDTH;
            }
            let bottom = &state.points[0].history;
            let top = &state.points[5].history;
            near(bottom.damage(), top.damage(), 1e-12);
            assert!(bottom.damage() > 0.5);
            assert_eq!(state.points[2].history.damage(), 0.0);
            assert_eq!(state.points[3].history.damage(), 0.0);
        }
        assert_eq!(damaged, 12);
        near(
            accepted.energy.degraded_equivalent_plastic_work_j,
            diagnostic,
            1e-15,
        );
    });
}

#[test]
fn unloading_reversal_and_restart_keep_complete_irreversible_shell_history() {
    with_cx(|_, cx| {
        let model = model(0.95);
        let section = section(&model, [1.0; 2]);
        let mesh = mesh();
        let prescribed = twist(&mesh, 0.06);
        let p = problem(&mesh, &section, &prescribed, &[]);
        let drilling = p.prepare_reference_drilling(cx).unwrap();
        let first = p
            .solve_ductile_with_energy(None, &drilling, tolerances(), cx)
            .unwrap();
        let target = twist(&mesh, 0.0595);
        let unloaded = problem(&mesh, &section, &target, &[])
            .solve_ductile_with_energy(Some(&first.solution), &drilling, tolerances(), cx)
            .unwrap();
        assert_eq!(unloaded.energy.damage_release_j, 0.0);
        assert_eq!(unloaded.energy.plastic_dissipation_j, 0.0);
        for (before, after) in first
            .solution
            .section_states
            .iter()
            .zip(&unloaded.solution.section_states)
        {
            for (before, after) in before.points.iter().zip(&after.points) {
                assert_eq!(before.history.damage(), after.history.damage());
                assert_eq!(before.history.native_state(), after.history.native_state());
            }
        }
        let reverse = twist(&mesh, -0.06);
        let p = problem(&mesh, &section, &reverse, &[]);
        let saved = unloaded.solution.clone();
        let reversed = p
            .solve_ductile_with_energy(Some(&unloaded.solution), &drilling, tolerances(), cx)
            .unwrap();
        let restarted = p
            .solve_ductile_with_energy(Some(&saved), &drilling, tolerances(), cx)
            .unwrap();
        assert_eq!(reversed.solution, restarted.solution);
        assert_eq!(reversed.energy, restarted.energy);
        assert_eq!(unloaded.solution, saved);
        assert!(reversed.energy.damage_release_j > 0.0);
        for (before, after) in saved
            .section_states
            .iter()
            .zip(&reversed.solution.section_states)
        {
            for (before, after) in before.points.iter().zip(&after.points) {
                assert!(after.history.damage() >= before.history.damage());
                assert!(after.history.native_state().alpha >= before.history.native_state().alpha);
            }
        }
        check_energy(&unloaded.energy);
        check_energy(&reversed.energy);
    });
}

#[test]
fn backward_euler_retains_native_motion_and_distinct_material_drilling_and_inertial_work() {
    with_cx(|_, cx| {
        let model = model(0.95);
        let dt = 0.2;
        let section = section(&model, [dt; 2]);
        let mesh = mesh();
        let mut prescribed = twist(&mesh, 0.06);
        prescribed[5].1 = 0.004;
        let p = problem(&mesh, &section, &prescribed, &[]);
        let drilling = p.prepare_reference_drilling(cx).unwrap();
        let inertia =
            PlasticShellReferenceInertia::prepare(&mesh, &p, PlasticShellDensity::Uniform(2.0), cx)
                .unwrap();
        // Damage releases previously stored energy. A first jump from virgin
        // zero strain has no old storage to release in the native endpoint
        // account, so start this motion from a physically elastic bend.
        let low = twist(&mesh, 0.01);
        let initial = problem(&mesh, &section, &low, &[])
            .solve_ductile_with_energy(None, &drilling, tolerances(), cx)
            .unwrap();
        let velocity = vec![[0.0; 6]; mesh.nodes.len()];
        let saved = velocity.clone();
        let direct = p
            .advance_backward_euler_with_reference_drilling(
                Some(&initial.solution),
                &velocity,
                &inertia,
                &drilling,
                dt,
                cx,
            )
            .unwrap();
        let accepted = p
            .advance_ductile_backward_euler_with_energy(
                Some(&initial.solution),
                &velocity,
                &inertia,
                &drilling,
                dt,
                tolerances(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.solution, direct.solution);
        assert_eq!(accepted.velocity, direct.velocity);
        assert_eq!(accepted.inertia, direct.inertia);
        assert_eq!(accepted.energy.drilling, direct.drilling);
        assert_eq!(accepted.duration_s, dt);
        assert_eq!(velocity, saved);
        assert!(accepted.energy.damage_release_j > 0.0);
        assert!(accepted.energy.plastic_dissipation_j > 0.0);
        assert!(accepted.energy.drilling.final_j > 0.0);
        assert!(accepted.energy.drilling.numerical_loss_j > 0.0);
        assert!(accepted.inertia.physical.kinetic_after_j > 0.0);
        assert!(accepted.inertia.drilling_gauge.kinetic_after_j > 0.0);
        assert!(accepted.inertia.physical.numerical_loss_j > 0.0);
        check_energy(&accepted.energy);
        near(accepted.balance_residual_j(), 0.0, 1e-10);
        assert!(
            p.advance_ductile_backward_euler_with_energy(
                None,
                &velocity,
                &inertia,
                &drilling,
                0.5 * dt,
                tolerances(),
                cx,
            )
            .is_err()
        );
        assert_eq!(velocity, saved);
    });
}

#[test]
fn caller_cutback_refreshes_actual_intervals_without_resetting_damage_or_reference_gauge() {
    with_cx(|_, cx| {
        let model = model(0.2);
        let original = section(&model, [1.0; 2]);
        let mesh = mesh();
        let zero = twist(&mesh, 0.0);
        let base = problem(&mesh, &original, &zero, &[]);
        let drilling = PlasticShellReferenceDrilling::prepare(&mesh, &base, cx).unwrap();
        let target = twist(&mesh, 0.06);
        assert!(
            problem(&mesh, &original, &target, &[])
                .solve_ductile_with_energy(None, &drilling, tolerances(), cx)
                .is_err()
        );
        let mut current = None;
        let mut duration = 0.0;
        for step in 1..=24 {
            let interval = section(&model, [1.0 / 24.0; 2]);
            let target = twist(&mesh, 0.06 * step as f64 / 24.0);
            let accepted = problem(&mesh, &interval, &target, &[])
                .solve_ductile_with_energy(current.as_ref(), &drilling, tolerances(), cx)
                .unwrap();
            check_energy(&accepted.energy);
            if let Some(old) = &current {
                for (a, b) in old
                    .section_states
                    .iter()
                    .zip(&accepted.solution.section_states)
                {
                    for (a, b) in a.points.iter().zip(&b.points) {
                        let increment = b.history.damage() - a.history.damage();
                        assert!(increment >= 0.0 && increment <= 0.2);
                    }
                }
            }
            duration += accepted.duration_s;
            current = Some(accepted.solution);
        }
        near(duration, 1.0, 1e-15);
        assert!(
            current.unwrap().section_states[0].points[0]
                .history
                .damage()
                > 0.5
        );
    });
}

#[test]
fn damage_owner_time_gauge_and_cancellation_refuse_without_mutating_checkpoints() {
    with_cx(|gate, cx| {
        let model = model(0.01);
        let material = section(&model, [1.0; 2]);
        let mesh = mesh();
        let low = twist(&mesh, 0.01);
        let p = problem(&mesh, &material, &low, &[]);
        let drilling = p.prepare_reference_drilling(cx).unwrap();
        let old = p
            .solve_ductile_with_energy(None, &drilling, tolerances(), cx)
            .unwrap();
        let saved = old.solution.clone();
        let high = twist(&mesh, 0.06);
        assert!(
            problem(&mesh, &material, &high, &[])
                .solve_ductile_with_energy(Some(&old.solution), &drilling, tolerances(), cx)
                .is_err()
        );
        let equivalent = model.clone();
        let foreign = section(&equivalent, [1.0; 2]);
        assert!(
            problem(&mesh, &foreign, &low, &[])
                .solve_ductile_with_energy(Some(&old.solution), &drilling, tolerances(), cx)
                .is_err()
        );
        let mixed = section(&model, [1.0, 0.5]);
        assert!(
            problem(&mesh, &mixed, &low, &[])
                .solve_ductile_with_energy(Some(&old.solution), &drilling, tolerances(), cx)
                .is_err()
        );
        let mut altered = problem(&mesh, &material, &low, &[]);
        altered.drilling_alpha *= 2.0;
        assert!(
            altered
                .solve_ductile_with_energy(Some(&old.solution), &drilling, tolerances(), cx)
                .is_err()
        );
        assert!(
            p.solve_ductile_with_energy(
                Some(&old.solution),
                &drilling,
                DuctileShellEnergyTolerances {
                    shell_work_j: 0.0,
                    ..tolerances()
                },
                cx
            )
            .is_err()
        );
        assert_eq!(old.solution, saved);
        gate.request();
        assert!(matches!(
            p.solve_ductile_with_energy(Some(&old.solution), &drilling, tolerances(), cx),
            Err(SmallStrainTetError::Cancelled)
        ));
        assert_eq!(old.solution, saved);
    });
}
