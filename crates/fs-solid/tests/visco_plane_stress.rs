//! G0/G1/G3/G4: actual P1/Q1 plane-stress trajectories retain material history.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{
    BranchEnergy, FiniteStrainBranch, ViscoHyperelastic, ViscoRemainderState, ViscoState,
};
use fs_solid::hyper2d::visco::{
    FiniteViscoPlaneEnergy, FiniteViscoPlaneError, FiniteViscoPlaneMode, FiniteViscoPlaneSolid,
    FiniteViscoStepControls, ViscoelasticLawField,
};
use fs_solid::hyper2d::{PlaneMaterialPoint, material_point_kinematics};
use fs_solid::{Mesh2, NewtonSettings, Patch, SolidError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 53,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn law() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 0.8,
                c01: 0.2,
                kappa: 25.0,
            },
            4.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 4.0 },
            tau_s: 0.2,
        }],
    )
    .unwrap()
}

fn mode(exact: bool) -> FiniteViscoPlaneMode {
    if exact {
        FiniteViscoPlaneMode::Incompressible
    } else {
        FiniteViscoPlaneMode::FiniteBulk
    }
}

fn local() -> PlaneStressSettings {
    PlaneStressSettings::default()
}

fn newton() -> NewtonSettings {
    NewtonSettings {
        tol: 1e-9,
        max_iters: 40,
        load_steps: 1,
        max_backtracks: 30,
    }
}

fn controls() -> FiniteViscoStepControls {
    FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 10.0,
    }
}

type Field<'a> = &'a dyn Fn(f64, f64) -> [f64; 2];

fn boundary<'a>(field: Field<'a>) -> [(Patch, Field<'a>); 4] {
    [
        (Patch::Left, field),
        (Patch::Right, field),
        (Patch::Bottom, field),
        (Patch::Top, field),
    ]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn accounts(e: &FiniteViscoPlaneEnergy) -> [f64; 12] {
    [
        e.stored_before_j,
        e.stored_after_j,
        e.dissipated_j,
        e.native_loading_work_j,
        e.incremental_potential_j,
        e.material_balance_error_j,
        e.internal_endpoint_work_j,
        e.dead_load_endpoint_work_j,
        e.support_endpoint_work_j,
        e.inertial_endpoint_work_j,
        e.equilibrium_work_residual_j,
        e.endpoint_work_defect_j,
    ]
}

fn check_energy(e: &FiniteViscoPlaneEnergy) {
    assert!(accounts(e).iter().all(|value| value.is_finite()));
    assert!(e.dissipated_j >= 0.0);
    close(e.material_balance_error_j, 0.0, 1e-9);
    close(e.equilibrium_work_residual_j, 0.0, 1e-9);
    close(
        e.native_loading_work_j,
        e.stored_after_j - e.stored_before_j + e.dissipated_j,
        1e-9,
    );
    close(
        e.endpoint_work_defect_j,
        e.internal_endpoint_work_j - e.native_loading_work_j,
        1e-9,
    );
    assert_eq!(e.inertial_endpoint_work_j, 0.0);
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    displacement: Vec<[f64; 2]>,
    points: Vec<PlaneMaterialPoint>,
    full: Option<Vec<ViscoState>>,
    remainder: Option<Vec<ViscoRemainderState>>,
    piola: Vec<[f64; 9]>,
    pressure: Option<Vec<f64>>,
    stored: f64,
    time: f64,
}

fn snapshot(body: &FiniteViscoPlaneSolid<'_>) -> Snapshot {
    Snapshot {
        displacement: body.displacement_m().to_vec(),
        points: body.material_points().to_vec(),
        full: body.point_states().map(<[_]>::to_vec),
        remainder: body.point_remainder_states().map(<[_]>::to_vec),
        piola: body.point_first_piola_pa().to_vec(),
        pressure: body.point_pressure_pa().map(<[_]>::to_vec),
        stored: body.stored_energy_j(),
        time: body.time_s(),
    }
}

#[test]
fn g0_nonaffine_q1_keeps_four_independent_histories_in_original_quadrature_order() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let law = law();
        let shape = |x: f64, y: f64| [0.12 * x + 0.09 * x * y, -0.025 * y + 0.05 * x * y];
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| shape(*x, *y)).collect();
        let thickness = 0.03;
        for exact in [false, true] {
            let body = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                thickness,
                &initial,
                local(),
                cx,
            )
            .unwrap();
            assert_eq!(body.material_points().len(), 4);
            assert_eq!(body.point_states().is_some(), !exact);
            assert_eq!(body.point_remainder_states().is_some(), exact);
            let before = snapshot(&body);
            let dt = 0.07;
            let step = body
                .advance(dt, &boundary(&shape), &[], newton(), controls(), cx)
                .unwrap();
            assert_eq!(snapshot(&body), before);
            assert_eq!(step.solid.displacement_m(), initial.as_slice());
            close(step.solid.time_s(), dt, 1e-15);
            let points = material_point_kinematics(&mesh, step.solid.displacement_m(), cx).unwrap();
            assert_eq!(step.solid.material_points(), points.as_slice());
            close(
                points.iter().map(|point| point.reference_area_m2).sum(),
                1.0,
                1e-14,
            );
            let mut stored = 0.0;
            let mut loss = 0.0;
            let mut work = 0.0;
            let mut potential = 0.0;
            for (index, point) in points.iter().enumerate() {
                assert_eq!((point.element, point.local_point), (0, index));
                let volume = point.reference_area_m2 * thickness;
                let expected = if exact {
                    let old = &body.point_remainder_states().unwrap()[index];
                    let f = old.deformation_gradient;
                    let native = law
                        .step_incompressible_plane_stress_response(
                            old,
                            &point.in_plane,
                            dt,
                            Some([f[2], f[5]]),
                            local(),
                        )
                        .unwrap();
                    assert_eq!(
                        step.solid.point_remainder_states().unwrap()[index],
                        native.material.state
                    );
                    close(
                        step.solid.point_pressure_pa().unwrap()[index],
                        native.plane_stress.pressure_pa,
                        1e-12,
                    );
                    (
                        native.material.state.stored,
                        native.material.dissipated_pa,
                        native.material.loading_work_pa,
                        native.material.incremental_potential_pa,
                        native.plane_stress.full_piola,
                    )
                } else {
                    let old = &body.point_states().unwrap()[index];
                    let f = old.deformation_gradient;
                    let native = law
                        .step_plane_stress_response(
                            old,
                            &point.in_plane,
                            dt,
                            Some([f[2], f[5], f[8]]),
                            local(),
                        )
                        .unwrap();
                    assert_eq!(
                        step.solid.point_states().unwrap()[index],
                        native.material.state
                    );
                    (
                        native.material.state.stored,
                        native.material.dissipated_pa,
                        native.material.loading_work_pa,
                        native.material.incremental_potential_pa,
                        native.material.first_piola_pa,
                    )
                };
                stored += volume * expected.0;
                loss += volume * expected.1;
                work += volume * expected.2;
                potential += volume * expected.3;
                for j in 0..9 {
                    close(
                        step.solid.point_first_piola_pa()[index][j],
                        expected.4[j],
                        1e-12,
                    );
                }
                for j in [2, 5, 8] {
                    close(step.solid.point_first_piola_pa()[index][j], 0.0, 1e-9);
                }
            }
            close(step.energy.stored_before_j, body.stored_energy_j(), 1e-13);
            close(step.energy.stored_after_j, stored, 1e-13);
            close(step.solid.stored_energy_j(), stored, 1e-13);
            close(step.energy.dissipated_j, loss, 1e-13);
            close(step.energy.native_loading_work_j, work, 1e-13);
            close(step.energy.incremental_potential_j, potential, 1e-13);
            check_energy(&step.energy);
            let metrics: Vec<_> = if exact {
                step.solid
                    .point_remainder_states()
                    .unwrap()
                    .iter()
                    .map(|state| state.inverse_viscous_cauchy_green[0])
                    .collect()
            } else {
                step.solid
                    .point_states()
                    .unwrap()
                    .iter()
                    .map(|state| state.inverse_viscous_cauchy_green[0])
                    .collect()
            };
            for i in 0..4 {
                for j in i + 1..4 {
                    assert_ne!(metrics[i], metrics[j]);
                }
            }
        }
    });
}

#[test]
fn g1_loaded_hold_and_unload_retain_memory_time_and_signed_endpoint_work_defects() {
    with_cx(|cx, _| {
        let law = law();
        for mesh in [
            Mesh2::triangles(1.0, 0.5, 1, 1),
            Mesh2::quads(1.0, 0.5, 1, 1),
        ] {
            for exact in [false, true] {
                let mut solid = FiniteViscoPlaneSolid::try_new_from(
                    &mesh,
                    ViscoelasticLawField::Uniform(&law),
                    mode(exact),
                    0.03,
                    &vec![[0.0; 2]; mesh.nodes.len()],
                    local(),
                    cx,
                )
                .unwrap();
                let mut time = 0.0;
                let mut held_stress = Vec::new();
                let mut finite_hold_work = 0.0;
                for (index, (stretch, dt)) in [
                    (0.15, 0.03),
                    (0.15, 0.07),
                    (0.15, 0.07),
                    (0.0, 0.03),
                    (0.0, 0.07),
                ]
                .into_iter()
                .enumerate()
                {
                    let shape = move |x: f64, _: f64| [stretch * x, 0.0];
                    let before = snapshot(&solid);
                    let step = solid
                        .advance(dt, &boundary(&shape), &[], newton(), controls(), cx)
                        .unwrap();
                    assert_eq!(snapshot(&solid), before);
                    time += dt;
                    close(step.solid.time_s(), time, 1e-14);
                    check_energy(&step.energy);
                    if index < 3 {
                        held_stress.push(step.solid.point_first_piola_pa()[0][0]);
                    }
                    if index == 1 || index == 2 {
                        close(step.energy.internal_endpoint_work_j, 0.0, 1e-12);
                        close(
                            step.energy.endpoint_work_defect_j,
                            -step.energy.native_loading_work_j,
                            1e-10,
                        );
                        assert!(step.energy.dissipated_j > 0.0);
                        finite_hold_work += step.energy.native_loading_work_j.abs();
                    }
                    solid = step.solid;
                }
                assert!(held_stress[2] < held_stress[1] && held_stress[1] < held_stress[0]);
                if !exact {
                    assert!(finite_hold_work > 1e-9);
                }
                let total_loss: f64 = if exact {
                    solid
                        .point_remainder_states()
                        .unwrap()
                        .iter()
                        .map(|state| state.dissipated)
                        .sum()
                } else {
                    solid
                        .point_states()
                        .unwrap()
                        .iter()
                        .map(|state| state.dissipated)
                        .sum()
                };
                assert!(total_loss > 0.0);
            }
        }
    });
}

#[test]
fn g1_constant_dead_traction_causes_actual_finite_and_exact_creep_displacement() {
    with_cx(|cx, _| {
        let law = law();
        let mesh = Mesh2::quads(1.0, 0.4, 2, 1);
        let thickness = 0.03;
        let clamp = |_: f64, _: f64| [0.0, 0.0];
        let pull = |_: f64, _: f64| [0.4, 0.0];
        let right = mesh.patch_nodes(Patch::Right);
        for exact in [false, true] {
            let mut solid = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                thickness,
                &vec![[0.0; 2]; mesh.nodes.len()],
                local(),
                cx,
            )
            .unwrap();
            let mut extension = Vec::new();
            let mut time = 0.0;
            for dt in [0.04, 0.08, 0.08] {
                let before = snapshot(&solid);
                let step = solid
                    .advance(
                        dt,
                        &[(Patch::Left, &clamp)],
                        &[(Patch::Right, &pull)],
                        newton(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(snapshot(&solid), before);
                extension.push(
                    right
                        .iter()
                        .map(|&node| step.solid.displacement_m()[node][0])
                        .sum::<f64>()
                        / right.len() as f64,
                );
                let reaction: f64 = step
                    .reactions_n
                    .iter()
                    .filter(|(dof, _)| dof % 2 == 0)
                    .map(|(_, force)| force)
                    .sum();
                close(reaction, -0.4 * 0.4 * thickness, 2e-10);
                check_energy(&step.energy);
                assert!(step.energy.dissipated_j > 0.0);
                assert_eq!(step.equilibrium.histories.len(), 1);
                assert!(step.equilibrium.histories[0].last().unwrap() <= &newton().tol);
                time += dt;
                close(step.solid.time_s(), time, 1e-14);
                solid = step.solid;
            }
            assert!(extension[0] > 0.0);
            assert!(extension[1] > extension[0] + 1e-6);
            assert!(extension[2] > extension[1] + 1e-6);
        }
    });
}

#[test]
fn g3_original_thickness_scales_all_physical_accounts_and_reactions_without_changing_history() {
    with_cx(|cx, _| {
        let law = law();
        let mesh = Mesh2::quads(1.0, 0.7, 1, 1);
        let zero = vec![[0.0; 2]; mesh.nodes.len()];
        let shape = |x: f64, y: f64| [0.12 * x + 0.03 * x * y, -0.02 * y + 0.01 * x * y];
        for exact in [false, true] {
            let thin = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                0.02,
                &zero,
                local(),
                cx,
            )
            .unwrap();
            let thick = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                0.1,
                &zero,
                local(),
                cx,
            )
            .unwrap();
            let a = thin
                .advance(0.06, &boundary(&shape), &[], newton(), controls(), cx)
                .unwrap();
            let b = thick
                .advance(0.06, &boundary(&shape), &[], newton(), controls(), cx)
                .unwrap();
            assert_eq!(a.solid.displacement_m(), b.solid.displacement_m());
            assert_eq!(a.solid.point_states(), b.solid.point_states());
            assert_eq!(
                a.solid.point_remainder_states(),
                b.solid.point_remainder_states()
            );
            assert_eq!(
                a.solid.point_first_piola_pa(),
                b.solid.point_first_piola_pa()
            );
            assert_eq!(a.solid.point_pressure_pa(), b.solid.point_pressure_pa());
            assert_eq!(a.solid.material_points(), b.solid.material_points());
            for (first, second) in accounts(&a.energy).into_iter().zip(accounts(&b.energy)) {
                close(second, 5.0 * first, 2e-12);
            }
            assert_eq!(a.reactions_n.len(), b.reactions_n.len());
            for ((i, p), (j, q)) in a.reactions_n.iter().zip(&b.reactions_n) {
                assert_eq!(i, j);
                close(*q, 5.0 * p, 2e-12);
            }
            close(
                b.solid.stored_energy_j(),
                5.0 * a.solid.stored_energy_j(),
                2e-12,
            );
        }
    });
}

#[test]
fn g0_late_endpoint_work_gate_refuses_the_complete_proposal_without_changing_accepted_state() {
    with_cx(|cx, _| {
        let law = law();
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let shape = |x: f64, _: f64| [0.18 * x, 0.0];
        for exact in [false, true] {
            let solid = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                0.03,
                &vec![[0.0; 2]; mesh.nodes.len()],
                local(),
                cx,
            )
            .unwrap();
            let before = snapshot(&solid);
            let proposal = solid
                .advance(0.07, &boundary(&shape), &[], newton(), controls(), cx)
                .unwrap();
            let defect = proposal.energy.endpoint_work_defect_j.abs();
            assert!(defect > 1e-8);
            let strict = FiniteViscoStepControls {
                maximum_endpoint_work_defect_j: 0.5 * defect,
                ..controls()
            };
            let error = solid
                .advance(0.07, &boundary(&shape), &[], newton(), strict, cx)
                .err()
                .expect("endpoint defect must refuse");
            match error {
                FiniteViscoPlaneError::EnergyGate {
                    quantity,
                    value_j,
                    maximum_abs_j,
                } => {
                    assert_eq!(quantity, "signed endpoint integration defect");
                    close(value_j, proposal.energy.endpoint_work_defect_j, 1e-12);
                    close(maximum_abs_j, 0.5 * defect, 1e-15);
                }
                other => panic!("expected the late endpoint-work gate, got {other:?}"),
            }
            assert_eq!(snapshot(&solid), before);
        }
    });
}

#[test]
fn g4_cancellation_during_load_query_leaves_accepted_strain_history_energy_and_time_unchanged() {
    for exact in [false, true] {
        with_cx(|cx, gate| {
            let law = law();
            let mesh = Mesh2::quads(1.0, 0.4, 2, 1);
            let initial: Vec<_> = mesh.nodes.iter().map(|[x, _]| [0.06 * x, 0.0]).collect();
            let solid = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                0.03,
                &initial,
                local(),
                cx,
            )
            .unwrap();
            let before = snapshot(&solid);
            let clamp = |_: f64, _: f64| [0.0, 0.0];
            let cancelling_load = |_: f64, _: f64| {
                gate.request();
                [0.4, 0.0]
            };
            let result = solid.advance(
                0.05,
                &[(Patch::Left, &clamp)],
                &[(Patch::Right, &cancelling_load)],
                newton(),
                controls(),
                cx,
            );
            assert!(matches!(
                result,
                Err(FiniteViscoPlaneError::Native(SolidError::Cancelled))
            ));
            assert_eq!(snapshot(&solid), before);
        });
    }
}

#[path = "visco_plane_stress/reference_change.rs"]
mod reference_change;
