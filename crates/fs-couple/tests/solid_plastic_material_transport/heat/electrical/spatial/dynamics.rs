//! Actual inertia, native deformed current and finite material feedback together.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

#[test]
fn rigid_acceleration_keeps_spatial_current_and_separate_kinetic_work() {
    let dt = 2e-6;
    let rho = 1e-5;
    let acceleration = 1e5;
    let initial_speed = 0.25;
    with_spatial(rho, false, |cx, _, thermal, solid, problem, prepared| {
        // Exact P1 consistent mass row sums: each adjacent tet contributes m/4.
        // Uniform force per mass gives rigid acceleration without any supports.
        let masses = [MASS / 2.0, MASS / 2.0, MASS / 2.0, MASS / 4.0, MASS / 4.0];
        let forces: Vec<_> = masses
            .iter()
            .enumerate()
            .map(|(node, mass)| (3 * node, mass * acceleration))
            .collect();
        let moving = ThermoplasticDynamics::try_new(
            solid.clone(),
            vec![[initial_speed, 0.0, 0.0]; NODES.len()],
            cx,
        )
        .unwrap();
        let boundary = electrical_boundary(problem);
        let first = moving
            .advance_spatial_electrical_candidate(
                thermal,
                problem,
                prepared,
                &boundary,
                dt,
                &[],
                &forces,
                IncrementSettings::default(),
                spatial_controls(),
                &[0.0; 2],
                cx,
            )
            .unwrap();
        let speed = initial_speed + dt * acceleration;
        let displacement = dt * speed;
        for node in 0..NODES.len() {
            close(first.dynamics().velocity_mps()[node][0], speed, 2e-10);
            close(
                first.problem().positions_m()[node][0],
                NODES[node][0] + displacement,
                2e-15,
            );
            for axis in 1..3 {
                close(first.dynamics().velocity_mps()[node][axis], 0.0, 2e-10);
            }
        }
        let c = dt * (DRIVE / L).powi(2) / (1000.0 * 5.0 * rho);
        let temperature = 300.0 + 2.0 * c / (1.0 + (1.0 + 4.0 * RHO_SLOPE * c).sqrt());
        for id in CELLS {
            close(
                first.thermal().temperature_k(id).unwrap(),
                temperature,
                3e-7,
            );
        }
        // Current is O(1e7 A/m²); 5e-6 A/m² is <1e-12 relative here.
        check_field_with_current_tolerance(
            first.problem(),
            &first.report().material.circuit_thermal.electrical.field,
            CELLS.map(|id| first.thermal().temperature_k(id).unwrap()),
            rho,
            [[-DRIVE / L, 0.0, 0.0]; 2],
            [1.0; 2],
            5e-6,
        );
        let report = first.report();
        let mass = 2.0 * MASS;
        let kinetic_change = 0.5 * mass * (speed * speed - initial_speed * initial_speed);
        let numerical_loss = 0.5 * mass * (dt * acceleration).powi(2);
        close(report.inertia.kinetic_change_j, kinetic_change, 1e-12);
        close(report.inertia.numerical_loss_j, numerical_loss, 1e-12);
        close(
            report.material.mechanical_work_j,
            mass * acceleration * displacement,
            1e-12,
        );
        close(
            report.material.mechanical_work_j,
            kinetic_change + numerical_loss,
            1e-12,
        );
        assert!(numerical_loss > 1e-6);
        close(report.material.deposited_mechanical_heat_j, 0.0, 1e-12);
        close(
            report
                .material
                .circuit_thermal
                .thermal
                .stored_energy_change_j,
            report.material.circuit_thermal.deposited_joule_j,
            1e-9,
        );
        close(report.material.combined_residual_j, 0.0, 1e-8);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(problem.positions_m(), &NODES);
        assert_eq!(
            moving.velocity_mps(),
            &vec![[initial_speed, 0.0, 0.0]; NODES.len()]
        );
    });
}

fn loaded_inertial_reference(rho: f64, dt: f64) -> [f64; 4] {
    // Eliminate the native free-apex DOF independently of the FE/Newton owner.
    // M_aa=m/10; M_a,y=m/20 couples the moving prescribed shear vertex.
    let apex = |temperature: f64| {
        let fy = 3e6 - 1e5 * (temperature - 300.0);
        let (g, h) = (5e8, 1e8);
        let q = |gamma: f64| 2.0 * g * 3.0_f64.sqrt() * (SHEAR * SHEAR + gamma * gamma).sqrt();
        let force = |gamma: f64| {
            let stress = 2.0 * g * gamma * (h + 3.0 * g * fy / q(gamma)) / (3.0 * g + h);
            VOLUME / L * stress
                + MASS / 10.0 * (2.0 * L * gamma) / (dt * dt)
                + MASS / 20.0 * (2.0 * SHEAR * L) / (dt * dt)
        };
        let (mut lo, mut hi) = (0.0, 0.03);
        let applied = SHEAR_STRESS_ZX * VOLUME / L;
        assert!(force(lo) < applied && force(hi) > applied);
        for _ in 0..80 {
            let gamma = 0.5 * (lo + hi);
            if force(gamma) < applied {
                lo = gamma;
            } else {
                hi = gamma;
            }
        }
        let gamma = 0.5 * (lo + hi);
        (2.0 * L * gamma, (q(gamma) - fy) / (3.0 * g + h), fy)
    };
    let capacity = MASS * 5.0;
    let bulk = LINK_AREA / (L / (2.0 * 0.01));
    let conductance = 1.0 / (1.0 / bulk + THERMAL_CONTACT / LINK_AREA);
    let radiation = 5.670374419e-8 * L * L / (2.0 / EMISSIVITY - 1.0);
    let exchange = |a: f64, b: f64| conductance * (a - b) + radiation * (a.powi(4) - b.powi(4));
    let joule = |t: f64, u: f64| {
        dt * VOLUME * (DRIVE / L).powi(2) * (1.0 + 4.0 * SHEAR * SHEAR + (u / L).powi(2))
            / (rho * (1.0 + RHO_SLOPE * (t - 300.0)))
    };
    let lower = |ta: f64| {
        let (mut lo, mut hi) = (300.0, 320.0);
        for _ in 0..80 {
            let tb = 0.5 * (lo + hi);
            if capacity * (tb - 300.0) - joule(tb, 0.0) - dt * exchange(ta, tb) < 0.0 {
                lo = tb;
            } else {
                hi = tb;
            }
        }
        0.5 * (lo + hi)
    };
    let residual = |ta: f64| {
        let (u, alpha, fy) = apex(ta);
        capacity * (ta - 300.0) - VOLUME * fy * alpha - joule(ta, u) + dt * exchange(ta, lower(ta))
    };
    let (mut lo, mut hi) = (300.0, 319.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let ta = 0.5 * (lo + hi);
        if residual(ta) < 0.0 {
            lo = ta;
        } else {
            hi = ta;
        }
    }
    let ta = 0.5 * (lo + hi);
    let (u, alpha, _) = apex(ta);
    [ta, lower(ta), u, alpha]
}

#[test]
fn resistivity_alone_changes_inertial_plastic_motion_and_warm_spatial_continuation() {
    let dt = 2e-5;
    let mut endpoints = Vec::new();
    for rho in [2e-5, 8e-5] {
        with_spatial(rho, false, |cx, _, thermal, solid, problem, prepared| {
            let moving =
                ThermoplasticDynamics::try_new(solid.clone(), vec![[0.0; 3]; NODES.len()], cx)
                    .unwrap();
            let boundary = electrical_boundary(problem);
            let first = moving
                .advance_spatial_electrical_candidate(
                    thermal,
                    problem,
                    prepared,
                    &boundary,
                    dt,
                    &supports(SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    spatial_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let actual = [
                first.thermal().temperature_k(CELLS[0]).unwrap(),
                first.thermal().temperature_k(CELLS[1]).unwrap(),
                first.dynamics().solid().solution().displacement_m[3][0],
                first.dynamics().solid().solution().element_states[0].alpha,
            ];
            let expected = loaded_inertial_reference(rho, dt);
            for i in 0..2 {
                close(actual[i], expected[i], 3e-6);
            }
            close(actual[2], expected[2], 3e-11);
            close(actual[3], expected[3], 2e-9);
            close(first.dynamics().velocity_mps()[3][0], actual[2] / dt, 2e-9);
            let scale = DRIVE / L;
            check_field_with_current_tolerance(
                first.problem(),
                &first.report().material.circuit_thermal.electrical.field,
                [actual[0], actual[1]],
                rho,
                [
                    [-scale, 2.0 * SHEAR * scale, actual[2] / L * scale],
                    [-scale, 2.0 * SHEAR * scale, 0.0],
                ],
                [1.0; 2],
                5e-6,
            );
            let report = &first.report().material;
            assert!(report.physical_plastic_heat_j > 0.0);
            assert!(first.report().inertia.numerical_loss_j > 1e-8);
            close(
                report.circuit_thermal.thermal.stored_energy_change_j,
                report.deposited_mechanical_heat_j + report.circuit_thermal.deposited_joule_j,
                1e-8,
            );
            close(report.combined_residual_j, 0.0, 1e-8);
            let kinetic = first.report().inertia.kinetic_after_j;
            let previous = first.dynamics().solid().solution().displacement_m.clone();
            let (heat, motion, current, _) = first.into_parts();
            let second = motion
                .advance_spatial_electrical_candidate(
                    &heat,
                    &current,
                    prepared,
                    &boundary,
                    2.0 * dt,
                    &supports(SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    spatial_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            close(second.report().inertia.kinetic_before_j, kinetic, 1e-12);
            for node in 0..NODES.len() {
                for axis in 0..3 {
                    let u = second.dynamics().solid().solution().displacement_m[node][axis];
                    close(
                        second.dynamics().velocity_mps()[node][axis],
                        (u - previous[node][axis]) / dt,
                        1e-12,
                    );
                    close(
                        second.problem().positions_m()[node][axis],
                        NODES[node][axis] + u,
                        1e-15,
                    );
                }
            }
            close(second.report().material.combined_residual_j, 0.0, 1e-8);
            assert_eq!(second.dynamics().time_s(), second.thermal().time_s());
            assert_eq!(moving.time_s(), 0.0);
            assert_eq!(problem.positions_m(), &NODES);
            endpoints.push(actual);
        });
    }
    assert!(endpoints[0][0] > endpoints[1][0] + 0.5);
    assert!(endpoints[0][2] > endpoints[1][2] + 1e-8);
    assert!(endpoints[0][3] > endpoints[1][3] + 1e-7);
}

#[test]
fn source_exit_and_late_work_budget_preserve_current_heat_and_velocity() {
    for short in [false, true] {
        with_spatial(2e-5, short, |cx, _, thermal, solid, problem, prepared| {
            let moving =
                ThermoplasticDynamics::try_new(solid.clone(), vec![[0.0; 3]; NODES.len()], cx)
                    .unwrap();
            let original_cells = thermal.cells().clone();
            let original_solid = moving.solid().solution().clone();
            let original_sources = problem.conductors().to_vec();
            let boundary = electrical_boundary(problem);
            let mut controls = spatial_controls();
            if !short {
                let good = moving
                    .advance_spatial_electrical_candidate(
                        thermal,
                        problem,
                        prepared,
                        &boundary,
                        2e-5,
                        &supports(SHEAR),
                        &applied(),
                        IncrementSettings::default(),
                        controls,
                        &[0.0; 2],
                        cx,
                    )
                    .unwrap();
                controls.coupled.maximum_mechanical_trials =
                    good.report().material.mechanical_trials - 1;
            }
            let result = moving.advance_spatial_electrical_candidate(
                thermal,
                problem,
                prepared,
                &boundary,
                2e-5,
                &supports(SHEAR),
                &applied(),
                IncrementSettings::default(),
                controls,
                &[0.0; 2],
                cx,
            );
            if short {
                assert!(result.is_err());
            } else {
                assert!(matches!(
                    result,
                    Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                        ElectricalSolidSourceError::MechanicalBudget
                    )))
                ));
            }
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(moving.solid().solution(), &original_solid);
            assert_eq!(moving.velocity_mps(), &vec![[0.0; 3]; NODES.len()]);
            assert_eq!(problem.conductors(), &original_sources);
            assert_eq!(problem.positions_m(), &NODES);
            assert_eq!(moving.time_s(), 0.0);
        });
    }
}
