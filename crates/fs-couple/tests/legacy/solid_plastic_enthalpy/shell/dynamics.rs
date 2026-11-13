//! G1 prescribed plastic stretch: native inertia changes work, never caloric heat.

use super::*;
use fs_material::IsotropicElastic;
use fs_material::plastic::section::LayeredPlaneStressSection;
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;
use fs_solid::plastic_shell::{PlasticShellSections, SmallStrainShellProblem};

#[test]
fn g1_dynamic_shell_heats_only_plastic_work_and_retains_history_when_motion_stops() {
    context(|cx| {
        let dt = 2e-6;
        let strain = 0.008;
        let fixture = Fixture::new(100.0);
        let mesh = mesh();
        let thermal = thermal(&fixture, &[0.0, 0.0], 0.0);
        let shell = bind_shell(&fixture, &mesh, &thermal, 2, cx);
        let reference = shell.prepare_reference_inertia(cx).unwrap();
        close(reference.reference_mass_kg(), RHO * VOLUME, 1e-18);
        assert_eq!(
            reference.reference_layer_densities(0),
            Some(&[RHO, RHO][..])
        );
        let prescribed: Vec<_> = (0..18)
            .map(|d| {
                (
                    d,
                    if d % 6 == 0 {
                        strain * mesh.nodes[d / 6][0]
                    } else {
                        0.0
                    },
                )
            })
            .collect();
        let velocity = vec![[0.0; 6]; mesh.nodes.len()];

        // Equal total mass is insufficient: this same geometry has an incorrect
        // density first moment, despite preserving the total reference inventory.
        let d = discretization();
        let misplaced_density_section = LayeredPlaneStressSection::new(
            (0..2)
                .map(|i| PlaneStressLayer {
                    bottom_z_m: -THICKNESS / 2.0 + i as f64 * THICKNESS / 2.0,
                    top_z_m: i as f64 * THICKNESS / 2.0,
                    subdivisions: 2,
                    law: IsotropicElastic::new(E, 0.2, 0.04).unwrap(),
                })
                .collect(),
            d.quadrature,
            d.reduction,
            d.section_budget,
        )
        .unwrap();
        let wrong_problem = SmallStrainShellProblem {
            mesh: &mesh,
            sections: PlasticShellSections::Uniform(&misplaced_density_section),
            coordinate_length_m: d.coordinate_length_m,
            rotation_limit_rad: d.rotation_limit_rad,
            drilling_alpha: d.drilling_alpha,
            prescribed: &prescribed,
            nodal_loads: &[],
            load_path: &[1.0],
            budget: d.shell_budget,
            settings: IncrementSettings::default(),
        };
        // The second wrong inventory preserves EACH layer's mean density,
        // but moves mass between corners and changes the local first moment.
        let corner_density = [vec![
            [0.5 * RHO, RHO, 1.5 * RHO],
            [1.5 * RHO, RHO, 0.5 * RHO],
        ]];
        for density in [
            PlasticShellDensity::Layers(&[0.5 * RHO, 1.5 * RHO]),
            PlasticShellDensity::PerElementCorners(&corner_density),
        ] {
            let wrong_reference = wrong_problem
                .prepare_reference_inertia(density, cx)
                .unwrap();
            close(
                wrong_reference.reference_mass_kg(),
                reference.reference_mass_kg(),
                1e-18,
            );
            assert!(
                shell
                    .advance_backward_euler(
                        &thermal,
                        &(),
                        dt,
                        &prescribed,
                        &[],
                        &velocity,
                        &wrong_reference,
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .is_err()
            );
        }

        let quasistatic = shell
            .advance(
                &thermal,
                &(),
                dt,
                &prescribed,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let dynamic = shell
            .advance_backward_euler(
                &thermal,
                &(),
                dt,
                &prescribed,
                &[],
                &velocity,
                &reference,
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let motion = dynamic.motion().unwrap();
        assert!(quasistatic.motion().is_none());
        assert_eq!(
            dynamic.shell().solution().coordinates,
            quasistatic.shell().solution().coordinates
        );
        assert_eq!(
            dynamic.shell().solution().section_states,
            quasistatic.shell().solution().section_states
        );
        assert!(
            dynamic
                .shell()
                .solution()
                .section_states
                .iter()
                .flat_map(|section| &section.points)
                .any(|point| point.history.alpha > 0.0)
        );

        // Integral N_1^2 dA = A/6 on the original P1 triangle. This is physical
        // translational kinetic energy; affine axial stretch has zero spin gauge.
        let expected_velocity = strain * L / dt;
        let kinetic_j = RHO * VOLUME * expected_velocity.powi(2) / 12.0;
        close(motion.velocity[1][0], expected_velocity, 1e-10);
        close(motion.inertia.physical.kinetic_after_j, kinetic_j, 1e-12);
        close(motion.inertia.physical.numerical_loss_j, kinetic_j, 1e-12);
        close(motion.inertia.drilling_gauge.kinetic_after_j, 0.0, 1e-15);
        close(
            dynamic.report().mechanical.inertial_work_j,
            2.0 * kinetic_j,
            1e-12,
        );
        close(
            dynamic.report().mechanical.external_work_j()
                - quasistatic.report().mechanical.external_work_j(),
            2.0 * kinetic_j,
            1e-11,
        );
        let heat_j = dynamic.report().mechanical.physical_heat_j;
        assert!(heat_j > 0.0);
        assert!(dynamic.report().mechanical.material_numerical_loss_j > 0.0);
        close(
            heat_j,
            quasistatic.report().mechanical.physical_heat_j,
            1e-12,
        );
        close(
            dynamic.report().thermal_audit.stored_energy_change_j,
            heat_j,
            config().energy_tolerance_j,
        );
        for cell in [CellId(100), CellId(101)] {
            close(
                dynamic.thermal().temperature_k(cell).unwrap(),
                quasistatic.thermal().temperature_k(cell).unwrap(),
                1e-10,
            );
            assert!(dynamic.thermal().temperature_k(cell).unwrap() > 300.0);
        }
        assert!(dynamic.report().combined_residual_j.abs() <= config().energy_tolerance_j);

        // Holding the accepted shape arrests its velocity. The backward-Euler
        // loss removes the kinetic energy but must not warm either original cell.
        let held = dynamic
            .shell()
            .advance_backward_euler(
                dynamic.thermal(),
                dynamic.observer(),
                2.0 * dt,
                &prescribed,
                &[],
                &motion.velocity,
                &reference,
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let stopped = held.motion().unwrap();
        for value in stopped.velocity.iter().flatten() {
            close(*value, 0.0, 1e-12);
        }
        close(stopped.inertia.physical.kinetic_before_j, kinetic_j, 1e-12);
        close(stopped.inertia.physical.kinetic_after_j, 0.0, 1e-15);
        close(stopped.inertia.physical.numerical_loss_j, kinetic_j, 1e-12);
        close(held.report().mechanical.inertial_work_j, 0.0, 1e-12);
        close(held.report().mechanical.physical_heat_j, 0.0, 1e-12);
        for cell in [CellId(100), CellId(101)] {
            close(
                held.thermal().temperature_k(cell).unwrap(),
                dynamic.thermal().temperature_k(cell).unwrap(),
                1e-8,
            );
        }
        for (before, after) in dynamic
            .shell()
            .solution()
            .section_states
            .iter()
            .flat_map(|section| &section.points)
            .zip(
                held.shell()
                    .solution()
                    .section_states
                    .iter()
                    .flat_map(|section| &section.points),
            )
        {
            close(after.history.alpha, before.history.alpha, 1e-11);
            for (a, b) in after
                .history
                .plastic_strain
                .iter()
                .zip(before.history.plastic_strain)
            {
                close(*a, b, 1e-11);
            }
        }
        assert!(held.report().combined_residual_j.abs() <= config().energy_tolerance_j);
    });
}
