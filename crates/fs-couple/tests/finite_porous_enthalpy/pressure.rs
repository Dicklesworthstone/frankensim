//! G1/G3: pressure-driven thermal compaction, exact work and atomic feedback.

use super::*;
use fs_couple::iqn_ils::{IqnIlsConfig, driver::CouplingMethod};
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;

fn surface(cx: &Cx<'_>) -> ClosedTetPressureSurface {
    ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 16, cx).unwrap()
}

#[test]
fn pressure_drives_free_compaction_and_original_heat_feedback_with_both_methods() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let thermal = initial(&[&fixture], None);
        let solid = bind(&fixture, &thermal, cx);
        let surface = surface(cx);
        let stretch: f64 = 0.8;
        let temperature = adiabatic_temperature(stretch, CP);
        let expected = hydrostatic(stretch, temperature);
        // sigma = P / s^2. This independent Cauchy pressure must recover the
        // free displacement and heat together using CURRENT triangle areas.
        let pressure = -expected.first_piola_pa / stretch.powi(2);
        let supports: Vec<_> = prescribed(stretch)
            .into_iter()
            .filter(|(dof, _)| *dof != 3)
            .collect();
        let exact_work = pressure * V * (1.0 - stretch.powi(3));
        let endpoint_work = V * expected.work_j_m3;
        for method in [
            CouplingMethod::RelaxedPicard,
            CouplingMethod::IqnIls(IqnIlsConfig::default()),
        ] {
            let next = solid
                .advance_with_pressure_and_coupling(
                    &thermal,
                    &(),
                    1.0,
                    &supports,
                    &[],
                    [0.0; 3],
                    &surface,
                    pressure,
                    settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    method,
                    cx,
                )
                .unwrap();
            let mechanics = next.solid().mechanics();
            let report = next.report();
            let energy = &report.material_energy;
            close(
                next.thermal().temperature_k(CELL).unwrap(),
                temperature,
                3e-6,
            );
            close(mechanics.temperatures_k()[0], temperature, 3e-6);
            close(mechanics.displacement_m()[1][0], (stretch - 1.0) * L, 2e-10);
            close(mechanics.element_states()[0].alpha, expected.alpha, 3e-8);
            close(
                mechanics.element_plastic_jacobians()[0],
                expected.plastic_jacobian,
                3e-8,
            );
            close(
                next.solid().current_density_kg_m3(cx).unwrap()[0],
                RHO / stretch.powi(3),
                2e-4,
            );
            close(energy.pressure.work_j, exact_work, 3e-7);
            close(energy.pressure.endpoint_work_j, endpoint_work, 3e-7);
            close(
                energy.pressure.endpoint_work_correction_j,
                exact_work - endpoint_work,
                3e-7,
            );
            close(
                energy.pressure.external_work_defect_j,
                V * expected.endpoint_work_defect_j_m3 + exact_work - endpoint_work,
                4e-7,
            );
            close(energy.material.dead_load_endpoint_work_j, 0.0, 0.0);
            close(energy.material.support_endpoint_work_j, 0.0, 3e-7);
            let heat = V * expected.caloric_heat_j_m3();
            close(energy.material.caloric_heat_j, heat, 3e-7);
            close(report.thermal_audit.stored_energy_change_j, heat, 3e-7);
            close(report.applied_heat_j, heat, 3e-7);
            close(report.external_heat_j, 0.0, 2e-9);
            assert!(energy.pressure.endpoint_work_correction_j.abs() > 0.1);
            // The sizeable pressure quadrature correction must not become heat.
            close(
                exact_work,
                mechanics.internal_energy_j() + heat + energy.pressure.external_work_defect_j,
                5e-7,
            );
            close(
                energy.material.balance_error_j,
                0.0,
                controls().material.equilibrium_work_tolerance_j,
            );
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            assert!(std::ptr::eq(
                mechanics.families()[0].caloric_source(),
                &fixture.caloric
            ));
            close(next.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.mechanics().element_states()[0].alpha, 0.0);
    });
}

#[test]
fn thermal_pressure_retains_inertia_and_rolls_back_late_observer_and_heat_refusals() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let thermal = initial(&[&fixture], None);
        let solid = bind(&fixture, &thermal, cx)
            .with_inertia(&[[0.01, -0.005, 0.003]; 4], cx)
            .unwrap();
        let surface = surface(cx);
        // With prescribed motion, zero pressure must recover the original
        // physical heat and motion path, including the off-diagonal mass.
        let original = solid
            .advance(
                &thermal,
                &(),
                1.0,
                &prescribed(0.9),
                &[],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let zero = solid
            .advance_with_pressure(
                &thermal,
                &(),
                1.0,
                &prescribed(0.9),
                &[],
                [0.0; 3],
                &surface,
                0.0,
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            zero.solid().mechanics().element_states(),
            original.solid().mechanics().element_states()
        );
        assert_eq!(
            zero.solid().mechanics().velocity_mps(),
            original.solid().mechanics().velocity_mps()
        );
        assert_eq!(zero.thermal().cells(), original.thermal().cells());
        close(
            zero.report().material_energy.material.caloric_heat_j,
            original.report().material_energy.caloric_heat_j,
            1e-12,
        );
        close(
            zero.report()
                .material_energy
                .material
                .endpoint_work_defect_j,
            original.report().material_energy.endpoint_work_defect_j,
            1e-12,
        );
        close(
            zero.report().combined_residual_j,
            original.report().combined_residual_j,
            1e-12,
        );
        assert_eq!(zero.report().material_energy.pressure.work_j, 0.0);
        let first = solid
            .advance_with_pressure(
                &thermal,
                &(),
                1.0,
                &prescribed(0.9),
                &[],
                [0.0; 3],
                &surface,
                1e6,
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, solid, (), _) = first.into_parts();
        let old_cells = thermal.cells().clone();
        let old_states = solid.mechanics().element_states().to_vec();
        let old_velocity = solid.mechanics().velocity_mps().unwrap().to_vec();
        let old_displacement = solid.mechanics().displacement_m().to_vec();
        let old_temperature = solid.mechanics().temperatures_k().to_vec();
        let old_internal = solid.mechanics().internal_energy_j();
        let observer = LateRefusal {
            completed_substeps: 0,
        };
        let failed = solid.advance_with_pressure(
            &thermal,
            &observer,
            2.0,
            &prescribed(0.82),
            &[],
            [0.0; 3],
            &surface,
            1.4e6,
            settings(),
            controls(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(
            failed,
            Err(FiniteThermoplasticError::Observer(
                DepositionError::Cancelled
            ))
        ));
        let short = solid.advance_with_pressure(
            &thermal,
            &(),
            2.0,
            &prescribed(0.82),
            &[],
            [0.0; 3],
            &surface,
            1.4e6,
            settings(),
            controls(),
            thermal_config(),
            FiniteThermoplasticConfig {
                maximum_iterations: 1,
                ..config()
            },
            cx,
        );
        assert!(matches!(
            short,
            Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })
        ));
        assert_eq!(observer.completed_substeps, 0);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_states);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_velocity);
        assert_eq!(solid.mechanics().displacement_m(), old_displacement);
        assert_eq!(solid.mechanics().temperatures_k(), old_temperature);
        assert_eq!(solid.mechanics().internal_energy_j(), old_internal);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);
        let retried = solid
            .advance_with_pressure(
                &thermal,
                &(),
                2.0,
                &prescribed(0.82),
                &[],
                [0.0; 3],
                &surface,
                1.4e6,
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let report = retried.report();
        let energy = &report.material_energy;
        let inertia = report.inertia.as_ref().unwrap();
        assert!(retried.solid().mechanics().element_states()[0].alpha > old_states[0].alpha);
        assert!(inertia.numerical_loss_j > 0.0);
        close(inertia.reference_mass_kg, RHO * V, 1e-18);
        close(
            energy.pressure.work_j,
            1.4e6 * V * (0.9_f64.powi(3) - 0.82_f64.powi(3)),
            1e-12,
        );
        close(
            report.thermal_audit.stored_energy_change_j,
            energy.material.caloric_heat_j,
            2e-9,
        );
        close(
            energy.material.dead_load_endpoint_work_j
                + energy.material.support_endpoint_work_j
                + energy.pressure.work_j,
            energy.material.internal_energy_change_j
                + report.thermal_audit.stored_energy_change_j
                + energy.pressure.external_work_defect_j
                + inertia.kinetic_change_j
                + inertia.numerical_loss_j,
            3e-9,
        );
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_states);
    });
}
