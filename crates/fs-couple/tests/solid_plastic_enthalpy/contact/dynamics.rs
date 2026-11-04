//! Independent one-DOF inertial pressure and finite-support heat balances.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

#[test]
fn inertia_controls_source_contact_pressure_and_numerical_loss_does_not_heat() {
    context(|cx| {
        let dt = 2e-6;
        let heat = 0.2;
        let a = Fixture::new(500.0);
        let b = Fixture::new(1000.0);
        let thermal = initial(&[&a, &b], &[(0, 1, 0.0)], &[heat / dt, 0.0]);
        let solid = bind(&thermal, &[&a], &NODES, &TETS, cx);
        let mut velocity = vec![[0.0; 3]; 4];
        velocity[3][2] = -1.0;
        let moving = ThermoplasticDynamics::try_new(solid, velocity, cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 11).map(|d| (d, 0.0)).collect();
        let (card, pin) = interface(&a, &b, 1.0, 8e6);
        let point = point();
        let contacts = [binding(&a, &b, &card, pin, &point)];
        let first = moving
            .advance_with_thermal_contacts(
                &thermal,
                &contacts,
                dt,
                &fixed,
                &[],
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                other,
                cx,
            )
            .unwrap();
        // P1 consistent M_33 = rho*V/10. The only free displacement sees C_zzzz.
        let mass = RHO * V / 10.0;
        let modulus = E * 0.8 / (1.2 * 0.6);
        let stiffness = modulus * V / (L * L);
        let u = -mass * dt / (mass + stiffness * dt * dt);
        let pressure = -modulus * u / L;
        let g = 0.5 * L * L / (0.002 - 0.0016 * pressure / 8e6);
        let ca = RHO * V * 500.0;
        let cb = RHO * V * 1000.0;
        let ta = 300.0 + heat / (ca + dt * g * cb / (cb + dt * g));
        let tb = 300.0 + dt * g * (ta - 300.0) / (cb + dt * g);
        let report = first.report();
        close(
            first.dynamics().solid().solution().displacement_m[3][2],
            u,
            1e-13,
        );
        close(first.dynamics().velocity_mps()[3][2], u / dt, 1e-8);
        close(report.material.contacts[0].compression_pa, pressure, 1e-4);
        close(
            first.thermal().temperature_k(CellId(100)).unwrap(),
            ta,
            1e-8,
        );
        close(
            first.thermal().temperature_k(CellId(101)).unwrap(),
            tb,
            1e-8,
        );
        close(
            report.material.balance.thermal_audit.stored_energy_change_j,
            heat,
            1e-9,
        );
        close(report.material.balance.physical_heat_j, 0.0, 1e-12);
        close(report.material.balance.reversible_heat_j, 0.0, 1e-12);
        assert!(report.inertia.numerical_loss_j > 1e-8);
        close(
            report.material.balance.combined_residual_j,
            0.0,
            controls().energy_tolerance_j,
        );
        let old_kinetic = report.inertia.kinetic_after_j;
        let (warm, moving_next, _) = first.into_parts();
        let second = moving_next
            .advance_with_thermal_contacts(
                &warm,
                &contacts,
                2.0 * dt,
                &fixed,
                &[],
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                other,
                cx,
            )
            .unwrap();
        let u2 = 2.0 * mass * u / (mass + stiffness * dt * dt);
        close(
            second.dynamics().solid().solution().displacement_m[3][2],
            u2,
            1e-13,
        );
        close(second.report().inertia.kinetic_before_j, old_kinetic, 1e-14);
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(moving.velocity_mps()[3][2], -1.0);
    });
}

#[test]
fn inertial_series_contact_refusal_preserves_velocity_and_all_original_heat() {
    context(|cx| {
        let a = Fixture::new(500.0);
        let b = Fixture::new(1000.0);
        let thermal = initial(&[&a, &b], &[(0, 1, 0.0)], &[0.0, 0.0]);
        let mut velocity = vec![[0.0; 3]; 4];
        velocity[3][2] = -1.0;
        let moving = ThermoplasticDynamics::try_new(
            bind(&thermal, &[&a], &NODES, &TETS, cx),
            velocity.clone(),
            cx,
        )
        .unwrap();
        let (card, pin) = interface(&a, &b, 1.0, 100.0);
        let point = point();
        let contacts = [binding(&a, &b, &card, pin, &point)];
        let fixed: Vec<_> = (0..12).filter(|d| *d != 11).map(|d| (d, 0.0)).collect();
        let original_cells = thermal.cells().clone();
        let result = moving.advance_with_series_thermal_contacts(
            &thermal,
            &contacts,
            2e-6,
            &fixed,
            &[],
            IncrementSettings::default(),
            controls(),
            &[0.0; 2],
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::from([(0, 0.1)]),
                    payload: (),
                })
            },
            cx,
        );
        assert!(result.is_err()); // Actual inertial compression exceeds the pinned source.
        assert_eq!(moving.velocity_mps(), velocity);
        assert_eq!(moving.solid().solution().displacement_m, vec![[0.0; 3]; 4]);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(moving.time_s(), 0.0);
    });
}
