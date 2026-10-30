//! G1/G3: original heat data affect continuing inertial plastic tetrahedra.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

#[test]
fn complete_heat_preserves_inertia_and_plastic_history_across_physical_intervals() {
    with_heat(0.7, 0.1, |cx, _, f, thermal, solid, prepared| {
        let mut velocity = vec![[0.0; 3]; NODES.len()];
        velocity[3][0] = 1.0;
        let moving = ThermoplasticDynamics::try_new(solid.clone(), velocity.clone(), cx).unwrap();
        let recorder = Recorder::default();
        let (first, audit) = moving
            .advance_heat_network_candidate(
                thermal,
                prepared,
                &recorder,
                1.0,
                &supports(SHEAR),
                &applied(),
                IncrementSettings::default(),
                thermal_config(),
                control(),
                cx,
            )
            .unwrap();
        let report = first.report();
        close(report.inertia.reference_mass_kg, 2.0 * MASS, 1e-18);
        // Only the upper apex initially moves; its P1 consistent diagonal mass
        // is one tenth of the actual upper tetrahedron inventory.
        close(report.inertia.kinetic_before_j, 0.5 * MASS / 10.0, 1e-16);
        assert!(report.inertia.kinetic_after_j > 0.0);
        assert!(report.inertia.numerical_loss_j > 0.0);
        assert!(report.material.physical_heat_j > 0.0);
        assert!(first.dynamics().solid().solution().element_states[0].alpha > 0.0);
        close(audit.external_heat_j, 0.0, 1e-10);
        // Closed internal radiation redistributes heat. Neither mechanical
        // backward-Euler loss is allowed to appear in caloric storage.
        close(
            report.material.thermal_audit.stored_energy_change_j,
            report.material.physical_heat_j
                + report.material.reversible_heat_j
                + audit.external_heat_j,
            control().energy_tolerance_j,
        );
        let closed = report.material.stored_mechanical_change_j
            + report.material.thermoelastic_storage_change_j
            + report.material.thermal_audit.stored_energy_change_j
            + report.material.numerical_loss_j
            + report.inertia.kinetic_change_j
            + report.inertia.numerical_loss_j
            - report.material.mechanical_work_j
            - audit.external_heat_j
            - report.material.thermal_audit.incoming_energy_j;
        close(closed, audit.combined_residual_j, 1e-12);
        assert!(audit.combined_residual_j.abs() <= control().energy_tolerance_j);
        assert!(audit.thermal_residual_j.abs() <= control().energy_tolerance_j);
        let sample = first.observer().first.sample().unwrap();
        for surface in sample.radiation.as_ref().unwrap().surfaces() {
            let side = CELLS.iter().position(|id| *id == surface.cell).unwrap();
            assert_eq!(
                surface.material.card_identity(),
                f.cards[side].content_hash()
            );
            assert_eq!(
                surface.phase.temperature_k(),
                first.thermal().temperature_k(surface.cell).unwrap()
            );
        }

        let kinetic = report.inertia.kinetic_after_j;
        let history = first.dynamics().solid().solution().element_states.clone();
        let displacement = first.dynamics().solid().solution().displacement_m.clone();
        let (second, second_audit) = first
            .dynamics()
            .advance_heat_network_candidate(
                first.thermal(),
                prepared,
                &first.observer().second,
                2.0,
                &supports(SHEAR - 0.0005),
                &[],
                IncrementSettings::default(),
                thermal_config(),
                control(),
                cx,
            )
            .unwrap();
        close(second.report().inertia.kinetic_before_j, kinetic, 1e-16);
        assert_ne!(
            second.dynamics().solid().solution().displacement_m,
            displacement
        );
        assert_eq!(second.dynamics().solid().solution().element_states, history);
        close(second.observer().second.durations.iter().sum(), 2.0, 1e-12);
        assert!(second_audit.combined_residual_j.abs() <= control().energy_tolerance_j);
        assert!(second_audit.thermal_residual_j.abs() <= control().energy_tolerance_j);
        assert_eq!(first.dynamics().time_s(), 1.0);
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(recorder.durations.is_empty());
    });
}

#[test]
fn original_emissivity_and_interface_data_change_actual_moving_solid_response() {
    let run = |emissivity, resistance| {
        with_heat(
            emissivity,
            resistance,
            |cx, _, _, thermal, solid, prepared| {
                let moving =
                    ThermoplasticDynamics::try_new(solid.clone(), vec![[0.0; 3]; NODES.len()], cx)
                        .unwrap();
                let (next, audit) = moving
                    .advance_heat_network_candidate(
                        thermal,
                        prepared,
                        &(),
                        1.0,
                        &supports(SHEAR),
                        &applied(),
                        IncrementSettings::default(),
                        thermal_config(),
                        control(),
                        cx,
                    )
                    .unwrap();
                assert!(audit.combined_residual_j.abs() <= control().energy_tolerance_j);
                [
                    next.thermal().temperature_k(CELLS[0]).unwrap(),
                    next.thermal().temperature_k(CELLS[1]).unwrap(),
                    next.dynamics().solid().solution().displacement_m[3][0],
                    next.dynamics().velocity_mps()[3][0],
                    next.dynamics().solid().solution().element_states[0].alpha,
                ]
            },
        )
    };
    let baseline = run(0.2, 0.01);
    for changed in [run(0.9, 0.01), run(0.2, 2.0)] {
        assert!((baseline[0] - changed[0]).abs() > 1e-4);
        assert!((baseline[1] - changed[1]).abs() > 1e-4);
        assert!((baseline[2] - changed[2]).abs() > 1e-10);
        assert!((baseline[3] - changed[3]).abs() > 1e-10);
        assert!((baseline[4] - changed[4]).abs() > 1e-8);
    }
}
