//! Original current/heat histories and a separately checked inertial trajectory.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

#[test]
fn finite_circuit_and_material_motion_continue_together_without_heating_numerical_loss() {
    let mut temperatures = Vec::new();
    for resistivity in [1e-5, 2e-5] {
        with_powered(resistivity, false, false, CHARGE, |cx, _, t, s, c, p| {
            let dt = 2e-6;
            let mut velocity = vec![[0.0; 3]; NODES.len()];
            velocity[3][0] = 2.0;
            let moving = ThermoplasticDynamics::try_new(s.clone(), velocity, cx).unwrap();
            let first = moving
                .advance_electrical_candidate(
                    t,
                    c,
                    p,
                    &[],
                    dt,
                    &supports(0.0),
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let mass = MASS / 10.0;
            let stiffness = (1.2e9 / 2.4) * VOLUME / (L * L);
            let u = mass * dt * 2.0 / (mass + stiffness * dt * dt);
            close(
                first.dynamics().solid().solution().displacement_m[3][0],
                u,
                1e-13,
            );
            let report = first.report();
            let r = &report.material;
            assert!(r.circuit_thermal.deposited_joule_j > 0.0);
            assert!(report.inertia.numerical_loss_j > 1e-8);
            close(r.physical_plastic_heat_j, 0.0, 1e-12);
            close(r.reversible_heat_j, 0.0, 1e-12);
            close(
                r.circuit_thermal.thermal.stored_energy_change_j,
                r.circuit_thermal.deposited_joule_j,
                1e-9,
            );
            close(
                r.combined_residual_j,
                0.0,
                controls().coupled.total_energy_tolerance_j,
            );
            let endpoint = p.resolve(cx, first.thermal()).unwrap();
            let native = c
                .step_with_contacts_candidate(
                    endpoint.bulk.resistors,
                    endpoint.contacts,
                    &[],
                    dt,
                    controls().coupled.circuit,
                    || false,
                )
                .unwrap();
            assert_eq!(first.circuit(), native.circuit());
            temperatures.push(first.thermal().temperature_k(CELLS[0]).unwrap());
            let previous_kinetic = report.inertia.kinetic_after_j;
            let previous_charge = first.circuit().charges_c()[0];
            let (warm, moving_next, circuit_next, _) = first.into_parts();
            let second = moving_next
                .advance_electrical_candidate(
                    &warm,
                    &circuit_next,
                    p,
                    &[],
                    2.0 * dt,
                    &supports(0.0),
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            close(
                second.report().inertia.kinetic_before_j,
                previous_kinetic,
                1e-14,
            );
            let u2 = 2.0 * mass * u / (mass + stiffness * dt * dt);
            close(
                second.dynamics().solid().solution().displacement_m[3][0],
                u2,
                1e-13,
            );
            assert!(second.circuit().charges_c()[0] < previous_charge);
            assert_eq!(c.charges_c(), &[CHARGE]);
            assert_eq!(moving.time_s(), 0.0);
            assert_eq!(t.time_s(), 0.0);
        });
    }
    assert!(temperatures[0] > temperatures[1]); // Material data alone changes the heat.
}

#[test]
fn electrical_source_exit_rolls_back_motion_circuit_and_caloric_state() {
    with_powered(1e-5, false, true, CHARGE, |cx, _, t, s, c, p| {
        let mut velocity = vec![[0.0; 3]; NODES.len()];
        velocity[3][0] = 2.0;
        let moving = ThermoplasticDynamics::try_new(s.clone(), velocity.clone(), cx).unwrap();
        let original_cells = t.cells().clone();
        let result = moving.advance_electrical_candidate(
            t,
            c,
            p,
            &[],
            2e-6,
            &supports(0.0),
            &[],
            IncrementSettings::default(),
            controls(),
            &[0.0; 2],
            cx,
        );
        assert!(result.is_err()); // Heating exceeds this original resistivity source's 300.05 K limit.
        assert_eq!(moving.velocity_mps(), velocity);
        assert_eq!(moving.solid().solution(), s.solution());
        assert_eq!(t.cells(), &original_cells);
        assert_eq!(c.charges_c(), &[CHARGE]);
        assert_eq!(moving.time_s(), 0.0);
    });
}
