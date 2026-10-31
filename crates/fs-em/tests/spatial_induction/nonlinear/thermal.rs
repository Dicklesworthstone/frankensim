//! G1/G3/G4: retained tetrahedral H(B,T), conductivity and signed material work.
use super::*;
#[path = "../../../../fs-material/tests/nonlinear_magnetic_thermal/common.rs"]
mod material;
use fs_material::conductor::ElectricalPropertyBasis;

fn conductor(
    at: &fs_material::field_storage::nonlinear::thermal::ResolvedThermalReversibleMagneticMaterial,
) -> ResolvedConductorState {
    ResolvedConductorState::try_from_material_with_basis(
        at.material().clone(),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap()
}

#[test]
fn sourced_thermal_refresh_matches_independent_tetrahedral_work_and_continuation() {
    // Unit cube: |curl(W)|=2, M=3 sigma/40 and source work=J/2 delta(a).
    // On the second reference segment H_ref(2a)=8a-3, so its generalized
    // effort is 16s*a-6s. These scalar equations do not call field assembly.
    for slope in [-0.001, 0.002] {
        for alpha in [0.0, 0.6] {
            let card = material::card(slope);
            let law = material::law(&card);
            let start = law.at_temperature(300.0).unwrap();
            let p = model(
                1,
                SpatialNonlinearMaterial::Reversible(start.magnetic().clone()),
            );
            let initial = mode(&p, 0.8, 13.6);
            let state = history(&initial, &conductor(&start));
            let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
            let velocity: Vec<_> = p
                .positions_m()
                .iter()
                .map(|x| x.map(|v| alpha * v))
                .collect();
            let mut accepted = state.clone();
            let mut a0 = 0.8_f64;
            let mut s0 = 1.0_f64;
            for (index, temperature) in [350.0, 320.0].into_iter().enumerate() {
                let at = law.at_temperature(temperature).unwrap();
                let s = 1.0 + slope * (temperature - 300.0);
                let m = 3.0 * 18.0 / (40.0 * (1.0 + 0.005 * (temperature - 300.0)));
                let dt = 0.02;
                let expected =
                    (m * a0 + dt * (1.0 + 6.0 * s)) / (m + dt * (16.0 * s + 2.0 * alpha * m));
                assert!((0.5..2.0).contains(&expected));
                let rate = (expected - a0) / dt;
                let conductor_rate = rate + 2.0 * alpha * expected;
                let run = || {
                    accepted
                        .step_with_materials_candidate(
                            &[[2.0; 3]; 6],
                            &velocity,
                            dt,
                            vec![SpatialNonlinearMaterial::Reversible(at.magnetic().clone()); 6],
                            vec![Some(conductor(&at)); 6],
                            nonlinear_controls(),
                            &mut || false,
                        )
                        .unwrap()
                };
                let next = run();
                let replay = run();
                assert_eq!(next.state().potential_wb(), replay.state().potential_wb());
                assert_eq!(next.step().heat_j.to_bits(), replay.step().heat_j.to_bits());
                near(next.state().potential_wb()[edge], expected, 2e-9);
                let step = next.step();
                let old_energy = s0 * material::psi(2.0 * a0);
                let new_energy = s * material::psi(2.0 * expected);
                let parameter = (s - s0) * material::psi(2.0 * a0);
                near(step.stored_before_j, old_energy, 2e-9);
                near(step.stored_after_j, new_energy, 2e-9);
                near(step.stored_change_j, new_energy - old_energy, 2e-9);
                near(step.parameter_work_j, parameter, 2e-9);
                near(step.heat_j, m * conductor_rate.powi(2) * dt, 2e-9);
                near(
                    step.mechanical_work_j,
                    -2.0 * alpha * expected * m * conductor_rate * dt,
                    2e-9,
                );
                near(
                    step.numerical_dissipation_j,
                    8.0 * s * (expected - a0).powi(2),
                    2e-9,
                );
                near(step.impressed_work_j, expected - a0, 2e-9);
                near(step.energy_residual_j, 0.0, 2e-9);
                near(
                    step.magnetic_cells.iter().map(|c| c.parameter_work_j).sum(),
                    parameter,
                    2e-9,
                );
                for cell in &step.magnetic_cells {
                    near(
                        cell.free_energy_j[1] - cell.free_energy_j[0],
                        cell.energy_change_j,
                        2e-9,
                    );
                }
                assert_eq!(next.state().tree_edges(), state.tree_edges());
                assert_eq!(
                    next.state().conductors()[0]
                        .as_ref()
                        .unwrap()
                        .source_basis(),
                    ElectricalPropertyBasis::Resistivity
                );
                near(next.state().time_s(), (index + 1) as f64 * dt, 1e-15);
                a0 = expected;
                s0 = s;
                accepted = next.into_parts().0;
            }
            // The borrowed initial owner remains a complete accepted history.
            assert_eq!(state.time_s(), 0.0);
            near(state.potential_wb()[edge], 0.8, 2e-10);
            near(
                state.field(&mut || false).unwrap().stored_energy_j(),
                material::psi(1.6),
                2e-9,
            );
        }
    }
}

#[test]
fn source_changes_missing_conductors_and_late_cancellation_preserve_native_history() {
    let card = material::card(0.002);
    let law = material::law(&card);
    let start = law.at_temperature(300.0).unwrap();
    let p = model(
        1,
        SpatialNonlinearMaterial::Reversible(start.magnetic().clone()),
    );
    let initial = mode(&p, 0.8, 13.6);
    let state = history(&initial, &conductor(&start));
    let at = law.at_temperature(350.0).unwrap();
    let materials = vec![SpatialNonlinearMaterial::Reversible(at.magnetic().clone()); 6];
    let conductors = vec![Some(conductor(&at)); 6];
    let other = material::law(&material::card(0.001))
        .at_temperature(350.0)
        .unwrap();
    assert!(
        state
            .step_with_materials_candidate(
                &[[2.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.02,
                vec![SpatialNonlinearMaterial::Reversible(other.magnetic().clone()); 6],
                vec![Some(conductor(&other)); 6],
                nonlinear_controls(),
                &mut || false,
            )
            .is_err()
    );
    let mut missing = conductors.clone();
    missing[0] = None;
    assert!(
        state
            .step_with_materials_candidate(
                &[[2.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.02,
                materials.clone(),
                missing,
                nonlinear_controls(),
                &mut || false,
            )
            .is_err()
    );
    let mut calls = 0usize;
    state
        .step_with_materials_candidate(
            &[[2.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.02,
            materials.clone(),
            conductors.clone(),
            nonlinear_controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled_calls = 0usize;
    let failure = state.step_with_materials_candidate(
        &[[2.0; 3]; 6],
        &[[0.0; 3]; 8],
        0.02,
        materials,
        conductors,
        nonlinear_controls(),
        &mut || {
            cancelled_calls += 1;
            cancelled_calls >= calls
        },
    );
    assert!(matches!(
        failure,
        Err(SpatialNonlinearInductionError::Spatial(
            SpatialMagneticError::Cancelled
        ))
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.problem().materials(), p.materials());
    near(
        state.field(&mut || false).unwrap().stored_energy_j(),
        material::psi(1.6),
        2e-9,
    );
}
