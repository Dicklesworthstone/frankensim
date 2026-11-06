use super::*;
use fs_material::plastic::energy::EnergyTolerances;
use fs_solid::plastic_shell::energy::ShellEnergyTolerances;

fn gates() -> ShellEnergyTolerances {
    ShellEnergyTolerances {
        material: EnergyTolerances {
            energy_j_m3: 1e-9,
            yield_pa: 1e-8,
        },
        section_work_j_m2: 1e-10,
        shell_work_j: 1e-9,
    }
}

#[test]
fn yielded_cantilever_closes_true_material_heat_and_work_through_unloading() {
    with_cx(|_, cx| {
        let (mesh, fixed, loads) = cantilever(2, 1);
        let section = section(metal());
        let mut previous: Option<SmallStrainShellSolution<J2State>> = None;
        let mut previous_storage = 0.0;
        let mut accumulated_heat = 0.0;
        for factor in [0.25, 0.5, 0.75, 1.0, 0.0] {
            let step = problem(&mesh, &section, &fixed, &loads, &[factor], 0.37)
                .solve_with_energy(previous.as_ref(), gates(), cx)
                .unwrap();
            let energy = &step.energy;
            let mut native_yield_heat = 0.0;
            for e in 0..mesh.tris.len() {
                let area = mesh.facet(e).unwrap().area_m2 / 3.0;
                for q in 0..3 {
                    for (i, point) in section.points().iter().enumerate() {
                        let after = step.solution.section_states[3 * e + q].points[i]
                            .history
                            .alpha;
                        let before = previous.as_ref().map_or(0.0, |old| {
                            old.section_states[3 * e + q].points[i].history.alpha
                        });
                        native_yield_heat += area
                            * point.weight_m
                            * section.layers()[point.layer].law.yield_stress
                            * (after - before);
                    }
                }
            }
            near(energy.physical_dissipation_j, native_yield_heat, 1e-12);
            near(energy.initial.total_j(), previous_storage, 1e-11);
            near(energy.equilibrium_work_residual_j, 0.0, 1e-9);
            near(energy.balance_residual_j, 0.0, 1e-9);
            near(
                energy.dead_load_endpoint_work_j
                    + energy.support_endpoint_work_j
                    + energy.transverse_residual_work_j,
                energy.stored_change_j()
                    + energy.physical_dissipation_j
                    + energy.material_numerical_dissipation_j
                    + energy.drilling_final_j
                    - energy.drilling_initial_j
                    + energy.drilling_numerical_loss_j,
                1e-9,
            );
            accumulated_heat += energy.physical_dissipation_j;
            previous_storage = energy.final_state.total_j();
            previous = Some(step.solution);
        }
        assert!(accumulated_heat > 1e-10);
        let unloaded = previous.unwrap();
        assert!(
            unloaded
                .section_states
                .iter()
                .any(|s| s.points.iter().any(|p| p.history.alpha > 0.0))
        );
        assert!(mean_tip_rotation(&mesh, &unloaded).abs() > 1e-5);
    });
}

#[test]
fn pure_drilling_energy_is_numerical_and_cannot_become_material_heat() {
    with_cx(|_, cx| {
        let mesh = ShellMesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let section = section(metal());
        let prescribed: Vec<_> = (0..18)
            .map(|d| (d, if d == 5 { 0.01 } else { 0.0 }))
            .collect();
        let step = problem(&mesh, &section, &prescribed, &[], &[1.0], 0.37)
            .solve_with_energy(None, gates(), cx)
            .unwrap();
        let energy = step.energy;
        let expected =
            0.5 * DRILLING_ALPHA * (E * H / (2.0 * (1.0 + NU))) * (0.5 / 3.0) * 0.01_f64.powi(2);
        near(energy.drilling_final_j, expected, 1e-15);
        near(energy.drilling_numerical_loss_j, expected, 1e-15);
        near(energy.drilling_work_j, 2.0 * expected, 1e-15);
        near(energy.support_endpoint_work_j, 2.0 * expected, 1e-15);
        near(energy.physical_dissipation_j, 0.0, 0.0);
        near(energy.material_numerical_dissipation_j, 0.0, 0.0);
        near(energy.final_state.total_j(), 0.0, 0.0);
        assert!(
            problem(&mesh, &section, &prescribed, &[], &[0.5, 1.0], 0.37)
                .solve_with_energy(None, gates(), cx)
                .is_err()
        );
    });
}

#[path = "dynamics_energy.rs"]
mod dynamics_energy;
