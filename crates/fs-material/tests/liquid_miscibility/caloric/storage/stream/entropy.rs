//! Native liquid mixing entropy uses original Cp/T and accepted flow stages.

use super::*;

// The independently integrated pure Cp/T terms use their original common
// Tdatum=300 K. The pair's entropy is ideal mixing + C*x_A*x_B*ln(T/400).
pub(super) fn expected_s(y: [f64; 2], excess_cp: f64, t: f64) -> f64 {
    let n = [y[0] / 0.018, y[1] / 0.10];
    let total = n[0] + n[1];
    let pure = n[0] * (60.0 * (t / 300.0).ln() + 0.1 * (t - 300.0))
        + n[1] * (80.0 * (t / 300.0).ln() + 0.2 * (t - 300.0));
    let ideal: f64 = n
        .into_iter()
        .filter(|&moles| moles > 0.0)
        .map(|moles| -R * moles * (moles / total).ln())
        .sum();
    pure + ideal + excess_cp * (n[0] * n[1] / total) * (t / 400.0).ln()
}

#[test]
fn native_liquid_mixer_closes_original_entropy_on_be_and_sdirk_stages() {
    for excess_cp in [-40.0, 40.0] {
        let card = card(CompositionBasis::MoleFraction, 1.0);
        let liquid = liquid_with_excess_cp(&card, 100.0, excess_cp);
        let caloric =
            BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0, 0), "pressure")
                .unwrap();
        // Distinct cold brackets still share the same original pure-species
        // entropy reference. A bracket-relative mixing entropy would fail.
        let storage: Vec<_> = [[0.25, 0.75], [1.0, 0.0], [0.0, 1.0]]
            .into_iter()
            .enumerate()
            .map(|(i, fractions)| {
                let mut control = controls();
                control.temperature_bracket[0] = Temperature::new([290.0, 300.0, 320.0][i]);
                caloric
                    .stream_storage_from_mass_fractions(fractions, control, || false)
                    .unwrap()
            })
            .collect();
        let laws = conduction_laws(&storage);
        for (initial_t, inlet_t) in [(380.0, [350.0, 500.0]), (400.0, [400.0; 2])] {
            let inlet_h = [h(&storage[1], inlet_t[0]), h(&storage[2], inlet_t[1])];
            let inlet_s = [
                expected_s([1.0, 0.0], excess_cp, inlet_t[0]),
                expected_s([0.0, 1.0], excess_cp, inlet_t[1]),
            ];
            let flow = mixer_flow(inlet_h);
            let (thermal, conduction) = prepare(&laws, &[initial_t], &[2.0], false);
            let plain =
                ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                    .unwrap()
                    .with_implicit_advection(&thermal, &flow, || Ok(()))
                    .unwrap();
            let audited =
                ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                    .unwrap()
                    .with_implicit_advection(&thermal, &flow, || Ok(()))
                    .unwrap()
                    .with_implicit_advection_entropy(&thermal, || Ok(()))
                    .unwrap();
            let original = thermal.cells().clone();
            let dt = 20.0;
            let initial_h = expected_h([0.25, 0.75], excess_cp, initial_t);
            let incoming_h = 0.25 * inlet_h[0] + 0.75 * inlet_h[1];
            let incoming_s = 0.25 * inlet_s[0] + 0.75 * inlet_s[1];
            for second_order in [false, true] {
                let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
                let stage_h = (2.0 * initial_h + gamma * dt * incoming_h) / (2.0 + gamma * dt);
                let final_h = if second_order {
                    (2.0 * initial_h
                        + (1.0 - gamma) * dt * (incoming_h - stage_h)
                        + gamma * dt * incoming_h)
                        / (2.0 + gamma * dt)
                } else {
                    (2.0 * initial_h + dt * incoming_h) / (2.0 + dt)
                };
                let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
                let advance =
                    |system: &ImplicitConductionChain<'_>,
                     workspace: &mut ImplicitConductionWorkspace| {
                        if second_order {
                            system.advance_sdirk2_to(
                                &thermal,
                                dt,
                                implicit_controls(dt),
                                workspace,
                                || Ok(()),
                            )
                        } else {
                            system.advance_to(
                                &thermal,
                                dt,
                                implicit_controls(dt),
                                workspace,
                                || Ok(()),
                            )
                        }
                    };
                let (next, report) = advance(&audited, &mut workspace).unwrap();
                let (unclassified, unclassified_report) = advance(&plain, &mut workspace).unwrap();
                assert_eq!(next.cells(), unclassified.cells());
                assert_eq!(report.energy, unclassified_report.energy);
                near(next.cells()[&CellId(0)].specific_h_j_kg, final_h, 2e-6);
                let final_t = expected_t([0.25, 0.75], excess_cp, final_h);
                near(next.temperature_k(CellId(0)).unwrap(), final_t, 2e-8);
                if inlet_t == [400.0; 2] {
                    assert!(final_t < 399.0);
                }
                let stages = if second_order {
                    [((1.0 - gamma) * dt, stage_h), (gamma * dt, final_h)]
                } else {
                    [(dt, final_h), (0.0, final_h)]
                };
                let mut boundary = 0.0;
                let mut mixing = 0.0;
                let mut heat_at_cell = 0.0;
                for (weight, energy) in stages {
                    let t = expected_t([0.25, 0.75], excess_cp, energy);
                    let entropy = expected_s([0.25, 0.75], excess_cp, t);
                    let carried_heat = (incoming_h - energy) / t;
                    boundary += weight * (incoming_s - entropy);
                    mixing += weight * (carried_heat - (incoming_s - entropy));
                    heat_at_cell += weight * carried_heat;
                }
                let storage_change = 2.0
                    * (expected_s([0.25, 0.75], excess_cp, final_t)
                        - expected_s([0.25, 0.75], excess_cp, initial_t));
                let account = report.entropy;
                let carried = account.advection.unwrap();
                near(carried.boundary_entropy_input_j_k, boundary, 2e-7);
                near(carried.upwind_mixing_j_k, mixing, 2e-7);
                near(carried.balance_residual_j_k, 0.0, 2e-10);
                near(account.thermal_storage_change_j_k, storage_change, 2e-7);
                near(
                    account.temporal_defect_j_k,
                    storage_change - heat_at_cell,
                    2e-7,
                );
                assert_eq!(account.conduction_generation_j_k, 0.0);
                assert!(!account.has_unclassified_sources);
                assert_eq!(account.unclassified_source_at_cell_j_k, 0.0);
                assert_eq!(report.energy.substeps, 1);
                assert_eq!(report.rejected_intervals, 0);
            }
            assert!(
                audited
                    .advance_to(
                        &thermal,
                        dt,
                        implicit_controls(dt),
                        &mut ImplicitConductionWorkspace::new(1).unwrap(),
                        || Err(DepositionError::Cancelled)
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), &original);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(flow.audit().substeps, 0);
        }
    }
}
