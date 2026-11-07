//! G1/G3/G4: actual rigid tooling through original finite plastic and thermal history.
//! Diagonal constitutive roots and consistent-mass work are independent oracles.

use super::*;
use fs_contact::sustained::solve::SolveControls;
use fs_solid::hyper3::plastic::thermal::contact::{
    FiniteThermoplasticPlaneContactError, PlaneContactControls, RigidPlaneTool, TetPlaneContact,
};
use fs_solver::NewtonKrylovConfig;

fn contact_controls() -> PlaneContactControls {
    PlaneContactControls {
        solve: SolveControls {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-7,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            difference_step_n: 0.01,
            force_residual_tolerance_n: 1e-6,
            ..SolveControls::default()
        },
        maximum_normal_work_defect_j: 2e-6,
        work_balance_tolerance_j: 2e-6,
        ..PlaneContactControls::default()
    }
}

fn tool(axis: usize, before: f64, after: f64) -> RigidPlaneTool {
    let mut normal = [0.0; 3];
    let mut point_before_m = [0.0; 3];
    let mut translation_m = [0.0; 3];
    normal[axis] = -1.0;
    point_before_m[axis] = before;
    translation_m[axis] = after - before;
    RigidPlaneTool {
        normal,
        point_before_m,
        translation_m,
    }
}

fn row(initial_force_n: f64) -> TetPlaneContact {
    TetPlaneContact {
        node: 1,
        plane: 0,
        initial_force_n,
        scale_n_per_m: 1e4,
    }
}

fn supports(axis: usize) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3 + axis)
        .map(|dof| (dof, 0.0))
        .collect()
}

fn released(old: Point, temperature: f64) -> (Point, f64) {
    let (mut lower, mut upper) = (-0.5, 0.3);
    assert!(Point::at(lower, temperature, Some(old)).0.piola[0] < 0.0);
    assert!(Point::at(upper, temperature, Some(old)).0.piola[0] > 0.0);
    for _ in 0..80 {
        let a = 0.5 * (lower + upper);
        if Point::at(a, temperature, Some(old)).0.piola[0] < 0.0 {
            lower = a;
        } else {
            upper = a;
        }
    }
    Point::at(0.5 * (lower + upper), temperature, Some(old))
}

#[test]
fn moving_tool_heat_dwell_and_release_retain_permanent_shape_and_rotate_together() {
    with_cx(|cx| {
        for axis in [0, 1] {
            // A proper cyclic rotation moves original x to world y, including
            // the mesh, supports and tool. No synthetic compliance is supplied.
            let nodes = if axis == 0 {
                NODES
            } else {
                NODES.map(|v| [v[2], v[0], v[1]])
            };
            let mut solid = FiniteThermoplasticTetSolid::try_new(
                &nodes,
                &TETS,
                vec![family()],
                &[300.0],
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let (mut old, _) = Point::at(0.0, 300.0, None);
            let mut plane = 1.0;
            let mut guess = 0.0;
            for (index, (temperature, duration)) in
                [(300.0, 0.07), (450.0, 0.11), (450.0, 0.19), (300.0, 0.13)]
                    .into_iter()
                    .enumerate()
            {
                let (next, dissipation) = if index == 3 {
                    released(old, temperature)
                } else {
                    Point::at(-0.18, temperature, Some(old))
                };
                let end = if index == 3 { 1.1 } else { next.a.exp() };
                let candidate = solid
                    .advance_with_plane_contact(
                        &[temperature],
                        &supports(axis),
                        &[],
                        [0.0; 3],
                        duration,
                        &[tool(axis, plane, end)],
                        &[row(guess)],
                        settings(),
                        controls(),
                        contact_controls(),
                        cx,
                    )
                    .unwrap();
                let step = &candidate.step;
                let contact = &candidate.contact;
                let force = if index == 3 { 0.0 } else { -V * next.piola[0] };
                close(
                    step.solid.displacement_m()[1][axis],
                    next.a.exp() - 1.0,
                    3e-9,
                );
                close(contact.solve.forces_n[0][0], force, 2e-5);
                assert_eq!(contact.solve.forces_n[0][1..], [0.0, 0.0]);
                assert!(contact.solve.paths_m[0][0] >= 0.0);
                assert!(contact.minimum_gap_m >= 0.0);
                assert!(contact.path_windows > 0);
                close(contact.force_on_tool_n[0][axis], force, 2e-5);
                close(
                    contact.body_work_j,
                    -force * (next.a.exp() - old.a.exp()),
                    2e-6,
                );
                close(contact.tool_work_j, -force * (end - plane), 2e-6);
                close(contact.dead_load_work_j, 0.0, 1e-14);
                close(contact.load_decomposition_error_j, 0.0, 2e-6);
                close(contact.work_balance_error_j, 0.0, 2e-6);
                close(step.solid.element_states()[0].alpha, next.alpha, 3e-9);
                for diagonal in 0..3 {
                    let metric = (if diagonal == axis {
                        -2.0 * next.plastic
                    } else {
                        next.plastic
                    })
                    .exp();
                    close(
                        step.solid.element_states()[0].inverse_plastic_cauchy_green[diagonal],
                        metric,
                        3e-9,
                    );
                }
                let [_, _, _, _, reversible] = next.account(old, dissipation);
                close(step.energy.plastic_dissipation_j, V * dissipation, 2e-6);
                close(
                    step.energy.caloric_heat_j,
                    V * (dissipation + reversible),
                    2e-5,
                );
                assert_eq!(step.solid.temperatures_k(), &[temperature]);
                assert_eq!(step.solid.families(), &[family()]);
                if index == 1 {
                    assert!(next.alpha > old.alpha);
                }
                if index == 2 {
                    close(next.alpha, old.alpha, 1e-13);
                    close(step.energy.plastic_dissipation_j, 0.0, 1e-7);
                }
                if index == 3 {
                    assert!(next.a.exp() < 0.99);
                    assert!(contact.solve.paths_m[0][0] > 0.1);
                    close(contact.tool_work_j, 0.0, 1e-12);
                }
                guess = contact.solve.forces_n[0][0];
                solid = candidate.step.solid;
                old = next;
                plane = end;
            }
        }
    });
}

#[test]
fn reference_inertia_and_contact_arrival_defect_are_separate_from_heat_and_retry() {
    with_cx(|cx| {
        let solid = initial(cx)
            .with_inertia(&[6.0], &[[0.0; 3]; 4], cx)
            .unwrap();
        let (old, _) = Point::at(0.0, 300.0, None);
        let (next, dissipation) = Point::at(0.9_f64.ln(), 350.0, Some(old));
        let duration: f64 = 0.05;
        let velocity = (next.a.exp() - 1.0) / duration;
        // The original consistent P1 mass entry at the only free coordinate is .1 kg.
        let expected_force = -(V * next.piola[0] + 0.1 * velocity / duration);
        let tools = [tool(0, 1.02, 0.9)];
        let before_u = solid.displacement_m().to_vec();
        let before_f = solid.element_deformation().to_vec();
        let before_p = solid.element_states().to_vec();
        let before_v = solid.velocity_mps().unwrap().to_vec();
        let mut strict = contact_controls();
        strict.maximum_normal_work_defect_j = 0.0;
        let refused = solid
            .advance_with_plane_contact(
                &[350.0],
                &supports(0),
                &[],
                [0.0; 3],
                duration,
                &tools,
                &[row(0.0)],
                settings(),
                controls(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FiniteThermoplasticPlaneContactError::Energy {
                quantity: "absolute normal endpoint work defect",
                ..
            }
        ));
        assert_eq!(solid.displacement_m(), before_u);
        assert_eq!(solid.element_deformation(), before_f);
        assert_eq!(solid.element_states(), before_p);
        assert_eq!(solid.velocity_mps().unwrap(), before_v);
        assert_eq!(solid.temperatures_k(), &[300.0]);

        let mut admitted = contact_controls();
        admitted.maximum_normal_work_defect_j = 100.0;
        let candidate = solid
            .advance_with_plane_contact(
                &[350.0],
                &supports(0),
                &[],
                [0.0; 3],
                duration,
                &tools,
                &[row(0.0)],
                settings(),
                controls(),
                admitted,
                cx,
            )
            .unwrap();
        let step = &candidate.step;
        let contact = &candidate.contact;
        close(contact.solve.forces_n[0][0], expected_force, 2e-5);
        close(step.solid.displacement_m()[1][0], -0.1, 3e-9);
        close(step.solid.velocity_mps().unwrap()[1][0], velocity, 1e-7);
        close(step.solid.reference_mass_kg().unwrap(), 1.0, 1e-12);
        let inertia = step.inertia.as_ref().unwrap();
        close(inertia.kinetic_change_j, 0.05 * velocity.powi(2), 1e-7);
        close(inertia.numerical_loss_j, 0.05 * velocity.powi(2), 1e-7);
        close(contact.body_work_j, expected_force * 0.1, 2e-6);
        close(contact.tool_work_j, expected_force * 0.12, 2e-6);
        close(contact.normal_work_defect_j, expected_force * 0.02, 2e-6);
        close(
            contact.absolute_normal_work_defect_j,
            contact.normal_work_defect_j,
            1e-12,
        );
        let [_, _, _, _, reversible] = next.account(old, dissipation);
        close(
            step.energy.caloric_heat_j,
            V * (dissipation + reversible),
            2e-5,
        );
        assert!(contact.normal_work_defect_j > 1.0);
        assert!(inertia.numerical_loss_j > 0.0);
        close(contact.work_balance_error_j, 0.0, 2e-6);
        close(
            contact.tool_work_j + step.energy.support_endpoint_work_j,
            step.energy.internal_energy_change_j
                + step.energy.caloric_heat_j
                + step.energy.endpoint_work_defect_j
                + inertia.kinetic_change_j
                + inertia.numerical_loss_j
                + contact.normal_work_defect_j,
            2e-6,
        );
    });
}

#[test]
fn missing_contact_rows_and_cancelled_work_cannot_publish_penetrating_or_partial_history() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let solid = initial(&cx);
        // The only candidate is node 1, while prescribed node 2 crosses the
        // infinite plane. A zero candidate force is not proof of a safe mesh.
        let prescribed: Vec<_> = (0..12)
            .filter(|&dof| dof != 4)
            .map(|dof| (dof, if dof == 7 { 0.1 } else { 0.0 }))
            .collect();
        let refused = solid
            .advance_with_plane_contact(
                &[300.0],
                &prescribed,
                &[],
                [0.0; 3],
                0.1,
                &[tool(1, 1.0, 1.0)],
                &[row(0.0)],
                settings(),
                controls(),
                contact_controls(),
                &cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FiniteThermoplasticPlaneContactError::Penetration {
                plane: 0,
                node: 2,
                at_end: true,
                ..
            }
        ));
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.element_states(), &[family().initial_state()]);

        // A duplicate normal row would create a redundant, ambiguous force.
        assert!(
            solid
                .advance_with_plane_contact(
                    &[300.0],
                    &supports(0),
                    &[],
                    [0.0; 3],
                    0.1,
                    &[tool(0, 1.0, 0.9)],
                    &[row(0.0), row(0.0)],
                    settings(),
                    controls(),
                    contact_controls(),
                    &cx,
                )
                .is_err()
        );
        gate.request();
        assert!(matches!(
            solid.advance_with_plane_contact(
                &[300.0],
                &supports(0),
                &[],
                [0.0; 3],
                0.1,
                &[tool(0, 1.0, 0.9)],
                &[row(0.0)],
                settings(),
                controls(),
                contact_controls(),
                &cx,
            ),
            Err(FiniteThermoplasticPlaneContactError::Cancelled)
        ));
        assert_eq!(solid.temperatures_k(), &[300.0]);
        assert_eq!(solid.element_states(), &[family().initial_state()]);
    });
}
