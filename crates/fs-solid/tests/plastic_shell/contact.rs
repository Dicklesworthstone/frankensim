//! Native sustained contact: plastic history, physical moments and BE inertia.

use super::*;
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::SolveError;
use fs_solid::plastic_shell::contact::{
    PlasticShellContacts, ShellContactControls, ShellContactError, ShellContactStep,
    ShellContactWitness,
};
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn witness(gap: f64) -> ShellContactWitness {
    ShellContactWitness {
        id: 91,
        triangle: 0,
        barycentric: [0.0, 1.0, 0.0],
        offset_m: [0.0, 0.07, 0.5 * H],
        reference_obstacle_point_m: [1.0 + gap, 0.07, 0.5 * H],
        axes_world: [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        area_m2: 0.02,
    }
}

fn interval(end: [f64; 3]) -> ShellContactStep {
    ShellContactStep {
        obstacle_start_displacement_m: [0.0; 3],
        obstacle_end_displacement_m: end,
        friction_coefficient: 0.25,
        initial_force_n: [0.0; 3],
        scales: NaturalMapScales {
            normal_n_per_m: 100.0,
            tangent_n_per_m: 100.0,
        },
    }
}

fn controls() -> ShellContactControls {
    let mut controls = ShellContactControls::default();
    controls.solve.newton.absolute_tolerance = 1e-11;
    controls.solve.newton.relative_tolerance = 0.0;
    controls.solve.difference_step_n = 1e-5;
    controls.solve.force_residual_tolerance_n = 1e-9;
    controls.work_tolerance_j = 1e-10;
    controls
}

#[test]
fn g1_contact_press_yields_native_shell_and_sliding_work_keeps_drilling_separate() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = section(metal());
        let zero: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let old = problem(&mesh, &material, &zero, &[], &[1.0], 0.37)
            .solve(cx)
            .unwrap();
        let saved = old.clone();
        let depth = 0.008;
        let physical_rotation = 0.002;
        let drilling_rotation = 0.004;
        let target_x = -depth - 0.5 * H * physical_rotation;
        let prescribed: Vec<_> = (0..18)
            .filter(|d| *d != 6)
            .map(|d| {
                (
                    d,
                    match d {
                        10 => physical_rotation,
                        11 => drilling_rotation,
                        _ => 0.0,
                    },
                )
            })
            .collect();
        let p = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37);
        let gauge = p.prepare_reference_drilling(cx).unwrap();
        let contacts = PlasticShellContacts::prepare(&mesh, &p, &[witness(0.0)], 1, cx).unwrap();

        // Independent original-shell prescribed endpoint: same native local
        // returns, with no contact iterations or independently coded plastic law.
        let target: Vec<_> = (0..18)
            .map(|d| {
                (
                    d,
                    match d {
                        6 => target_x,
                        10 => physical_rotation,
                        11 => drilling_rotation,
                        _ => 0.0,
                    },
                )
            })
            .collect();
        let expected = problem(&mesh, &material, &target, &[], &[1.0], 0.37)
            .solve_from_with_reference_drilling(&old.coordinates, &old.section_states, &gauge, cx)
            .unwrap();
        let expected_normal = -expected.solution.increments[0].reactions[6];
        assert!(expected_normal > 0.0);
        let mut step = interval([-depth, 0.002, 0.0]);
        step.initial_force_n[0] = 0.8 * expected_normal;
        let solved = p
            .solve_from_with_contacts(
                &old.coordinates,
                &old.section_states,
                &gauge,
                &contacts,
                &[step],
                controls(),
                cx,
            )
            .unwrap();
        let report = &solved.contacts[0];
        near(solved.solution.coordinates[1][0], target_x, 2e-9);
        near(report.force_n[0], expected_normal, 2e-8);
        near(report.force_n[1], 0.25 * report.force_n[0], 2e-9);
        near(report.slip_increment_m[0], -0.002, 2e-12);
        near(report.shell_work_j, report.force_n[0] * depth, 2e-10);
        near(
            report.friction_loss_j,
            0.25 * report.force_n[0] * 0.002,
            2e-12,
        );
        near(report.normal_pressure_pa, report.force_n[0] / 0.02, 1e-12);
        near(report.work_residual_j, 0.0, controls().work_tolerance_j);
        assert_eq!(report.normal_numerical_loss_j, 0.0);
        assert_eq!(report.id, 91);
        assert!(report.final_gap_m >= 0.0);
        assert!(solved.drilling.final_j > 0.0);
        assert!(solved.motion.is_none());
        assert!(solved.contact_solve.response_evaluations > 1);
        for (actual, reference) in solved
            .solution
            .section_states
            .iter()
            .zip(&expected.solution.section_states)
        {
            for (a, b) in actual.points.iter().zip(&reference.points) {
                near(a.history.alpha, b.history.alpha, 2e-8);
            }
        }
        assert!(
            solved
                .solution
                .section_states
                .iter()
                .flat_map(|s| &s.points)
                .any(|p| p.history.alpha > 0.0)
        );
        assert_eq!(old, saved);

        // The same prepared witness survives release. The contact load vanishes
        // while the native material retains permanent compression and history.
        let free: Vec<_> = zero.iter().copied().filter(|(d, _)| *d != 6).collect();
        let release_problem = problem(&mesh, &material, &free, &[], &[1.0], 0.37);
        let release = ShellContactStep {
            obstacle_start_displacement_m: step.obstacle_end_displacement_m,
            obstacle_end_displacement_m: [0.01, 0.002, 0.0],
            initial_force_n: [0.0; 3],
            ..step
        };
        let released = release_problem
            .solve_from_with_contacts(
                &solved.solution.coordinates,
                &solved.solution.section_states,
                &gauge,
                &contacts,
                &[release],
                controls(),
                cx,
            )
            .unwrap();
        near(released.contacts[0].force_n[0], 0.0, 1e-10);
        near(released.contacts[0].friction_loss_j, 0.0, 1e-12);
        assert!(released.contacts[0].final_gap_m > 0.0);
        assert!(released.solution.coordinates[1][0] < -0.001);
    });
}

#[test]
fn g1_loaded_contact_restarts_with_identical_gaps_through_elastic_unload_and_reload() {
    with_cx(|_, cx| {
        let mesh = triangle();
        // Use the native plastic law below yield: a closed elastic cycle must
        // not invent plastic history during repeated contact response trials.
        let material = section(metal());
        let fixed: Vec<_> = (0..18).filter(|d| *d != 6).map(|d| (d, 0.0)).collect();
        let p = problem(&mesh, &material, &fixed, &[], &[1.0], 0.37);
        let gauge = p.prepare_reference_drilling(cx).unwrap();
        let weighted = ShellContactWitness {
            barycentric: [0.2, 0.3, 0.5],
            reference_obstacle_point_m: [0.3, 0.5 + 0.07, 0.5 * H],
            ..witness(0.0)
        };
        let contacts = PlasticShellContacts::prepare(&mesh, &p, &[weighted], 1, cx).unwrap();
        let mut old = p.solve(cx).unwrap();
        let virgin = old.section_states[0].points[0].history.clone();
        let mut policy = controls();
        // Representable positive closure residuals have an explicit small
        // budget. The signed floating-point gap is never clamped.
        policy.maximum_normal_numerical_loss_j = 1e-12;
        let mut start = [0.0; 3];
        let mut force = [0.0; 3];
        let mut previous_gap = 0.0_f64;
        let stiffness = E * H * 0.5 / (1.0 - NU * NU);
        for depth in [
            0.00010352029791283241,
            0.000032590353908501725,
            0.00006321098765,
            0.00014789012345,
            0.0,
            0.00010352029791283241,
            0.00006321098765,
            0.0,
        ] {
            let step = ShellContactStep {
                obstacle_start_displacement_m: start,
                friction_coefficient: 0.0,
                initial_force_n: force,
                ..interval([-depth, 0.0, 0.0])
            };
            let saved = old.clone();
            let next = p
                .solve_from_with_contacts(
                    &old.coordinates,
                    &old.section_states,
                    &gauge,
                    &contacts,
                    &[step],
                    policy,
                    cx,
                )
                .unwrap();
            let report = &next.contacts[0];
            assert_eq!(report.initial_gap_m.to_bits(), previous_gap.to_bits());
            assert!(report.final_gap_m >= 0.0);
            near(next.solution.coordinates[1][0], -depth / 0.3, 2e-9);
            near(report.force_n[0], stiffness * depth / (0.3 * 0.3), 2e-8);
            near(report.work_residual_j, 0.0, policy.work_tolerance_j);
            assert_eq!(report.friction_loss_j, 0.0);
            assert!(report.normal_numerical_loss_j <= policy.maximum_normal_numerical_loss_j);
            assert_eq!(next.solution.increments[0].evolving_elements, 0);
            for point in next
                .solution
                .section_states
                .iter()
                .flat_map(|state| &state.points)
            {
                assert_eq!(point.history, virgin);
            }
            assert_eq!(old, saved);
            start = step.obstacle_end_displacement_m;
            force = report.force_n;
            previous_gap = report.final_gap_m;
            old = next.solution;
        }
        near(old.coordinates[1][0], 0.0, 2e-9);
        for resultant in old.section_resultants.iter().flatten() {
            near(*resultant, 0.0, 2e-7);
        }
    });
}

#[test]
fn g1_dynamic_contact_uses_original_mass_and_separates_open_gap_approach_from_heat() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = section(IsotropicElastic::new(E, NU, 0.05).unwrap());
        let fixed: Vec<_> = (0..18).filter(|d| *d != 6).map(|d| (d, 0.0)).collect();
        let p = problem(&mesh, &material, &fixed, &[], &[1.0], 0.37);
        let gauge = p.prepare_reference_drilling(cx).unwrap();
        let density = 10.0;
        let inertia = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(density), cx)
            .unwrap();
        // An initial positive gap is admitted only with an explicit numerical
        // approach budget. Nothing in the closing work is deposited as heat.
        let contacts = PlasticShellContacts::prepare(&mesh, &p, &[witness(0.001)], 1, cx).unwrap();
        let step = ShellContactStep {
            friction_coefficient: 0.0,
            ..interval([-0.0018, 0.0, 0.0])
        };
        let velocity = vec![[0.0; 6]; 3];
        let dt = 0.1;
        let mut allowed = controls();
        allowed.maximum_normal_numerical_loss_j = 1e-3;
        let solved = p
            .advance_backward_euler_with_contacts(
                None,
                &velocity,
                &inertia,
                &gauge,
                &contacts,
                &[step],
                dt,
                allowed,
                cx,
            )
            .unwrap();
        let report = &solved.contacts[0];
        let depth = 0.0018 - report.initial_gap_m;
        let k = E * H * 0.5 / (1.0 - NU * NU);
        let mass_entry = density * H * 0.5 / 6.0;
        let expected_force = (k + mass_entry / (dt * dt)) * depth;
        near(report.force_n[0], expected_force, 2e-8);
        near(solved.solution.coordinates[1][0], -depth, 2e-10);
        near(
            report.normal_numerical_loss_j,
            report.force_n[0] * report.initial_gap_m,
            2e-15,
        );
        assert!(report.normal_numerical_loss_j > 0.0);
        assert_eq!(report.friction_loss_j, 0.0);
        near(
            report.shell_work_j + report.normal_numerical_loss_j,
            report.obstacle_work_j,
            1e-10,
        );
        let motion = solved.motion.as_ref().unwrap();
        near(motion.velocity[1][0], -depth / dt, 2e-9);
        near(
            motion.inertia.physical.kinetic_after_j,
            0.5 * mass_entry * (depth / dt).powi(2),
            2e-10,
        );
        assert!(motion.inertia.physical.numerical_loss_j > 0.0);
        assert!(matches!(
            p.advance_backward_euler_with_contacts(
                None,
                &velocity,
                &inertia,
                &gauge,
                &contacts,
                &[step],
                dt,
                controls(),
                cx,
            ),
            Err(ShellContactError::NormalApproachBudget { .. })
        ));
    });
}

#[test]
fn g3_contact_geometry_source_budget_and_cancellation_refuse_without_changing_history() {
    with_cx(|gate, cx| {
        let mesh = triangle();
        let material = section(metal());
        let fixed: Vec<_> = (0..18).filter(|d| *d != 6).map(|d| (d, 0.0)).collect();
        let p = problem(&mesh, &material, &fixed, &[], &[1.0], 0.37);
        let gauge = p.prepare_reference_drilling(cx).unwrap();
        let contacts = PlasticShellContacts::prepare(&mesh, &p, &[witness(0.0)], 1, cx).unwrap();
        assert_eq!(contacts.witnesses(), &[witness(0.0)]);
        assert!(
            PlasticShellContacts::prepare(&mesh, &p, &[witness(0.0), witness(0.0)], 2, cx).is_err()
        );
        let old = p.solve(cx).unwrap();
        let saved = old.clone();
        let step = interval([-0.008, 0.002, 0.0]);
        let mut limited = controls();
        limited.solve.maximum_response_evaluations = 1;
        assert!(matches!(
            p.solve_from_with_contacts(
                &old.coordinates,
                &old.section_states,
                &gauge,
                &contacts,
                &[step],
                limited,
                cx,
            ),
            Err(ShellContactError::Contact(
                SolveError::ResponseBudget { .. }
            ))
        ));
        assert!(matches!(
            p.solve_from_with_contacts(
                &old.coordinates,
                &old.section_states[..2],
                &gauge,
                &contacts,
                &[step],
                controls(),
                cx,
            ),
            Err(ShellContactError::Contact(SolveError::Response(
                SmallStrainTetError::InvalidInput { .. }
            )))
        ));
        let mut changed = mesh.clone();
        changed.nodes[1][0] *= 1.1;
        assert!(
            problem(&changed, &material, &fixed, &[], &[1.0], 0.37)
                .solve_from_with_contacts(
                    &old.coordinates,
                    &old.section_states,
                    &gauge,
                    &contacts,
                    &[step],
                    controls(),
                    cx,
                )
                .is_err()
        );
        let penetrating = ShellContactStep {
            obstacle_start_displacement_m: [-0.001, 0.0, 0.0],
            ..step
        };
        assert!(matches!(
            p.solve_from_with_contacts(
                &old.coordinates,
                &old.section_states,
                &gauge,
                &contacts,
                &[penetrating],
                controls(),
                cx,
            ),
            Err(ShellContactError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            p.solve_from_with_contacts(
                &old.coordinates,
                &old.section_states,
                &gauge,
                &contacts,
                &[step],
                controls(),
                cx,
            ),
            Err(ShellContactError::Solid(SmallStrainTetError::Cancelled))
        ));
        assert_eq!(old, saved);
    });
}
