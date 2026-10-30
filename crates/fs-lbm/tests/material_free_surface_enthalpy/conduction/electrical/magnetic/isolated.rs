//! Native isolated loop references on the actual material/enthalpy/VOF owners.
//! Synthetic source cards are numerical fixtures, not measured liquid data.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::electrical::FreeSurfaceCurrentBoundary;

const ORIGIN: [usize; 3] = [1, 3, 3];
const DX: f64 = 0.001;
const DT: f64 = 0.001;
const SPEED: f64 = 0.02;

fn isolated() -> FreeSurfaceElectricalInput {
    FreeSurfaceElectricalInput {
        origin: ORIGIN,
        dimensions: [2, 2, 1],
        setup: ElectricSetup::default(),
        contacts: vec![],
        maximum_cells: 4,
        maximum_solves: 32,
        power_tolerance_w: 1e-12,
    }
}
fn loop_field(b: f64) -> FreeSurfaceMagneticInput {
    // Bz changes with y, not z, so this original Cartesian fixture is solenoidal.
    // Uniform +y motion has a nonzero circulation of u x B around this square.
    FreeSurfaceMagneticInput {
        magnetic_flux_density_t: vec![[0.0, 0.0, b], [0.0, 0.0, b], [0.0, 0.0, -b], [0.0, 0.0, -b]],
        maximum_velocity_change_m_s: 0.03,
    }
}
fn velocities(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> Vec<[f64; 3]> {
    (0..4)
        .map(|cell| {
            liquid
                .flow()
                .lattice()
                .sample_at(ORIGIN[0] + cell % 2, ORIGIN[1] + cell / 2, ORIGIN[2])
                .unwrap()
                .velocity
        })
        .collect()
}

#[test]
fn isolated_material_loop_brakes_heats_and_preserves_the_gauge_cells_physics() {
    let run = |basis, scale, b| {
        let card = electrical_card(basis, scale, false, false, 500.0);
        with_initial_flow(
            &card,
            basis,
            false,
            Some([0.0, SPEED, 0.0]),
            [0.0; 3],
            |liquid, gate| {
                let mut liquid = liquid
                    .with_isolated_electrical_core(isolated(), 1e-12, gate)
                    .unwrap()
                    .with_imposed_magnetic_field(loop_field(b))
                    .unwrap();
                let initial_velocity = velocities(&liquid);
                let before = liquid.total_enthalpy_j().unwrap();
                let report = liquid
                    .step_with_conduction(&[], &[], controls(), gate)
                    .unwrap();
                let electric = report.electrical.as_ref().unwrap();
                let state = &electric.last_applied;
                let magnetic = state.magnetic.as_ref().unwrap();
                // Four equal edges: R=rho_e*dx/A=1/scale ohm. The two horizontal
                // motional voltages add around the loop: emf=2*u*B*dx, I=emf/(4R).
                // Every voxel owns two half-resistances, INCLUDING gauge cell 0.
                let r = 1.0 / scale;
                let edge_emf = SPEED * b * DX;
                let current = 2.0 * edge_emf / (4.0 * r);
                let power = 4.0 * r * current * current;
                let cell_force_y = -current * DX * b / 2.0;
                assert!(state.solution.electrodes.is_empty());
                assert_eq!(state.solution.delivered_w, 0.0);
                assert_eq!(electric.electrode_energy_j, 0.0);
                assert_eq!(
                    state.current_boundary,
                    FreeSurfaceCurrentBoundary::Isolated {
                        current_tolerance_a: 1e-12
                    }
                );
                assert_eq!(
                    liquid.electrical_current_boundary(),
                    Some(state.current_boundary)
                );
                assert_eq!(magnetic.velocity_m_s, initial_velocity);
                for cell in 0..4 {
                    near(state.solution.joule_w[cell], r * current * current, 2e-10);
                    near(
                        magnetic.lorentz_force_density_n_m3[cell][1],
                        cell_force_y / DX.powi(3),
                        2e-10,
                    );
                    near(
                        electric.lorentz_impulse_n_s[cell][1],
                        cell_force_y * DT,
                        2e-10,
                    );
                    assert_eq!(state.materials[cell].source_basis(), basis);
                    assert_eq!(state.materials[cell].material().properties().len(), 6);
                }
                // Reconstruct all actual edge currents and every cell's KCL,
                // independently of the native gauge and residual accumulation.
                let mut balance = [0.0; 4];
                let phi = &state.solution.potential_v;
                for (a, b, emf) in [
                    (0, 1, edge_emf),
                    (0, 2, 0.0),
                    (1, 3, 0.0),
                    (2, 3, -edge_emf),
                ] {
                    let edge_current = (phi[a] - phi[b] + emf) / r;
                    near(edge_current.abs(), current.abs(), 2e-10);
                    balance[a] += edge_current;
                    balance[b] -= edge_current;
                }
                assert!(balance.iter().all(|residual| residual.abs() < 1e-14));
                near(electric.joule_heat_j, power * DT, 2e-10);
                near(electric.sampled_lorentz_work_j, -power * DT, 2e-10);
                near(report.flow.external_heat_j, power * DT, 2e-10);
                // Subtracting the approximately 29 J body storage to observe
                // 4e-11 J heating needs a floating-point accumulation allowance.
                assert!((liquid.total_enthalpy_j().unwrap() - before - power * DT).abs() < 1e-12);
                assert!(electric.energy_defect_j.abs() < 1e-20);
                assert_eq!(electric.substeps, 4);
                assert_eq!(electric.current_solves, 8);
                let accepted_velocity = velocities(&liquid);
                assert!(accepted_velocity.iter().map(|u| u[1]).sum::<f64>() < 4.0 * SPEED);
                assert!(
                    electric
                        .thermal_endpoint
                        .temperature_k
                        .iter()
                        .zip(&state.temperature_k)
                        .all(|(end, start)| end > start)
                );
                for &p in &state.coordinates {
                    assert_eq!(
                        liquid
                            .flow()
                            .lattice()
                            .force_density_at(p[0], p[1], p[2])
                            .unwrap(),
                        [0.0; 3]
                    );
                }
                let bits = liquid.flow().lattice().grid().state_bits();
                let next = liquid
                    .step_with_conduction(&[], &[], controls(), gate)
                    .unwrap();
                let next_current = next.electrical.as_ref().unwrap();
                assert_eq!(
                    next_current
                        .last_applied
                        .magnetic
                        .as_ref()
                        .unwrap()
                        .velocity_m_s,
                    accepted_velocity
                );
                assert!(next_current.sampled_lorentz_work_j < 0.0);
                assert!(next_current.energy_defect_j.abs() < 1e-20);
                (bits, accepted_velocity, electric.joule_heat_j)
            },
        )
    };
    let reference = run(ElectricalPropertyBasis::Resistivity, 1.0, 10.0);
    let equivalent = run(ElectricalPropertyBasis::Conductivity, 1.0, 10.0);
    assert_eq!(
        reference, equivalent,
        "changing only equivalent source basis must preserve native physics"
    );
    let reversed = run(ElectricalPropertyBasis::Resistivity, 1.0, -10.0);
    assert_eq!(
        reference, reversed,
        "reversing B reverses induced I, so braking and heat retain their sign"
    );
    let conducting = run(ElectricalPropertyBasis::Conductivity, 2.0, 10.0);
    near(conducting.2, 2.0 * reference.2, 2e-10);
    assert_ne!(reference.0, conducting.0);
    assert!(
        conducting.1.iter().map(|u| u[1]).sum::<f64>()
            < reference.1.iter().map(|u| u[1]).sum::<f64>()
    );
}

#[test]
fn isolated_uniform_translation_polarizes_without_inventing_joule_heat_or_drag() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    for magnetic in [false, true] {
        with_initial_flow(
            &card,
            basis,
            false,
            Some([0.0, SPEED, 0.0]),
            [0.0; 3],
            |liquid, gate| {
                let mut liquid = liquid
                    .with_isolated_electrical_core(isolated(), 1e-12, gate)
                    .unwrap();
                if magnetic {
                    liquid = liquid
                        .with_imposed_magnetic_field(FreeSurfaceMagneticInput {
                            magnetic_flux_density_t: vec![[0.0, 0.0, 10.0]; 4],
                            maximum_velocity_change_m_s: 0.03,
                        })
                        .unwrap();
                }
                let report = liquid
                    .step_with_conduction(&[], &[], controls(), gate)
                    .unwrap();
                let electric = report.electrical.as_ref().unwrap();
                assert!(electric.last_applied.solution.electrodes.is_empty());
                assert_eq!(electric.electrode_energy_j, 0.0);
                assert!(electric.joule_heat_j < 1e-24);
                assert!(electric.sampled_lorentz_work_j.abs() < 1e-18);
                assert!(report.flow.external_heat_j.abs() < 1e-24);
                if magnetic {
                    let state = &electric.last_applied;
                    assert!(
                        state
                            .magnetic
                            .as_ref()
                            .unwrap()
                            .lorentz_force_density_n_m3
                            .iter()
                            .flatten()
                            .all(|f| f.abs() < 1e-8)
                    );
                    for (cell, potential) in state.solution.potential_v.iter().enumerate() {
                        assert!((*potential - (cell % 2) as f64 * SPEED * 10.0 * DX).abs() < 1e-12);
                    }
                } else {
                    assert!(electric.last_applied.magnetic.is_none());
                    assert!(electric.lorentz_impulse_n_s.is_empty());
                }
            },
        );
    }
}

#[test]
fn isolated_reference_refusals_and_late_current_budget_preserve_the_native_liquid() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    with_initial_flow(
        &card,
        basis,
        false,
        Some([0.0, SPEED, 0.0]),
        [0.0; 3],
        |liquid, gate| {
            assert!(
                liquid
                    .clone()
                    .with_isolated_electrical_core(current(0.0), 1e-12, gate)
                    .is_err()
            );
            for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
                assert!(
                    liquid
                        .clone()
                        .with_isolated_electrical_core(isolated(), tolerance, gate)
                        .is_err()
                );
            }
            let bits = liquid.flow().lattice().grid().state_bits();
            let energy = liquid.total_enthalpy_j().unwrap();
            let source = liquid.material_at(ORIGIN).unwrap().clone();
            let mut input = isolated();
            input.maximum_solves = 1; // Applied current succeeds; first endpoint refuses.
            let mut liquid = liquid
                .with_isolated_electrical_core(input, 1e-12, gate)
                .unwrap()
                .with_imposed_magnetic_field(loop_field(10.0))
                .unwrap();
            assert!(matches!(
                liquid.step_with_conduction(&[], &[], controls(), gate),
                Err(FreeSurfaceEnthalpyError::Invalid(
                    "electrical current-solve budget exhausted"
                ))
            ));
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
            assert_eq!(liquid.material_at(ORIGIN).unwrap(), &source);
            assert_eq!(liquid.accepted_steps(), 0);
            let mut retry = liquid
                .with_isolated_electrical_core(isolated(), 1e-12, gate)
                .unwrap()
                .with_imposed_magnetic_field(loop_field(10.0))
                .unwrap();
            assert!(
                retry
                    .step_with_conduction(&[], &[], controls(), gate)
                    .unwrap()
                    .electrical
                    .unwrap()
                    .joule_heat_j
                    > 0.0
            );
            assert_eq!(retry.accepted_steps(), 1);
        },
    );
}
