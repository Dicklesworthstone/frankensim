//! Native fixed-R finite-body trajectories on original phase charts.

use super::super::*;
use fs_conduction::lumped::exchange::pair::{LumpedPairExchangeConfig, solve_lumped_pair_exchange};
use fs_material::phase::exchange::pair::FinitePhaseExchangeControls;

fn caloric(id: u8, datum: f64, knots: &[(f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([id; 32]),
        knots
            .iter()
            .map(|&(h, t, liquid)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + h,
                temperature_k: t,
                liquid_mass_fraction: liquid,
                bulk_density_kg_m3: 1000.0 - 200.0 * liquid,
            })
            .collect(),
    )
    .unwrap()
}
fn pair_body(curve: &EquilibriumEnthalpyPhaseCurve, mass: f64) -> LumpedEnthalpyBody<'_> {
    LumpedEnthalpyBody::try_new(
        "finite contact body",
        mass,
        0.01,
        0.0,
        0.0,
        0.001,
        35.0,
        curve,
    )
    .unwrap()
}
fn controls(initial: [f64; 2], duration: f64) -> LumpedPairExchangeConfig {
    LumpedPairExchangeConfig {
        initial_specific_enthalpies_j_kg: initial,
        duration_s: duration,
        energy_tolerance_j: 1e-8,
        integration: FinitePhaseExchangeControls::default(),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_native_direct_pair_uses_actual_mass_resistance_and_latent_residence() {
    let hot = caloric(0x91, 0.0, &[(0.0, 200.0, 0.0), (1200.0, 500.0, 0.0)]);
    let pcm = caloric(
        0x92,
        0.0,
        &[
            (0.0, 300.0, 0.0),
            (100.0, 350.0, 0.0),
            (300.0, 350.0, 1.0),
            (900.0, 500.0, 1.0),
        ],
    );
    let a = pair_body(&hot, 2.0);
    let b = pair_body(&pcm, 1.0);
    let r = contact_resistance();
    with_cx(|cx| {
        for (duration, heat, cold_t, hot_t, liquid, entropy) in [
            (
                1.0,
                153.7117830631862,
                350.0,
                480.7860271171017,
                0.2685589153159312,
                0.1482773780981327,
            ),
            (
                5.0,
                494.8265774837081,
                398.7066443709271,
                438.1466778145365,
                1.0,
                0.3444660844927276,
            ),
        ] {
            let result =
                solve_lumped_pair_exchange(cx, &a, &b, &r, controls([1200.0, 0.0], duration))
                    .unwrap();
            let phase = result.caloric;
            assert_eq!(result.body_identities, [a.identity(), b.identity()]);
            assert_eq!(result.resistance_k_per_w, 1.0);
            assert_eq!(phase.conductance_time_j_k, duration);
            near(phase.heat_first_to_second_j, heat, 2e-10);
            near(phase.states[0].temperature_k(), hot_t, 2e-12);
            near(phase.states[1].temperature_k(), cold_t, 2e-12);
            near(phase.states[1].liquid_mass_fraction(), liquid, 2e-12);
            near(phase.total_entropy_change_j_k, entropy, 3e-13);
            near(phase.specific_heat_j_kg[0], -heat / 2.0, 2e-10);
            near(phase.specific_heat_j_kg[1], heat, 2e-10);
            assert!(phase.energy_residual_j.iter().all(|r| r.abs() <= 1e-8));
            assert!(result.maximum_biot < LUMPED_BIOT_CEILING);
        }
        use fs_conduction::interface::{
            ResistanceUncertainty, SeriesThermalResistance, ThermalResistanceTerm,
        };
        let twice_r = SeriesThermalResistance::new(vec![
            ThermalResistanceTerm::declared(
                "twice the contact resistance",
                2.0,
                ResistanceUncertainty::Unstated,
                "manufactured 2 K/W contact",
            )
            .unwrap(),
        ])
        .unwrap();
        let scaled =
            solve_lumped_pair_exchange(cx, &a, &b, &twice_r, controls([1200.0, 0.0], 10.0))
                .unwrap();
        near(
            scaled.caloric.heat_first_to_second_j,
            494.8265774837081,
            2e-10,
        );
        assert_eq!(scaled.caloric.conductance_time_j_k, 5.0);
        // R is already total K/W; the body's 0.01 m2 is not multiplied into G.
        let stiff =
            solve_lumped_pair_exchange(cx, &a, &b, &r, controls([1200.0, 0.0], 100.0)).unwrap();
        near(stiff.caloric.heat_first_to_second_j, 600.0, 2e-10);
        assert!(stiff.caloric.states[0].temperature_k() >= stiff.caloric.states[1].temperature_k());
    });
}

#[test]
fn g0_equal_and_opposite_state_rounding_cannot_bypass_either_local_energy_gate() {
    let datum = 2.0_f64.powi(52);
    let curve = caloric(0x93, datum, &[(0.0, 300.0, 0.0), (400.0, 500.0, 0.0)]);
    let a = pair_body(&curve, 1.0);
    let b = pair_body(&curve, 1.0);
    let r = contact_resistance();
    with_cx(|cx| {
        let cfg = LumpedPairExchangeConfig {
            energy_tolerance_j: 1.0,
            ..controls([datum + 400.0, datum], 0.003)
        };
        let accepted = solve_lumped_pair_exchange(cx, &a, &b, &r, cfg).unwrap();
        near(
            accepted.caloric.heat_first_to_second_j,
            200.0 * -(-0.003_f64).exp_m1(),
            2e-14,
        );
        assert_eq!(accepted.total_energy_residual_j, 0.0);
        assert!(
            accepted
                .caloric
                .energy_residual_j
                .iter()
                .all(|e| e.abs() > 0.1)
        );
        assert!(
            solve_lumped_pair_exchange(
                cx,
                &a,
                &b,
                &r,
                LumpedPairExchangeConfig {
                    energy_tolerance_j: 1e-6,
                    ..cfg
                }
            )
            .is_err()
        );
        assert!(
            solve_lumped_pair_exchange(
                cx,
                &a,
                &b,
                &r,
                LumpedPairExchangeConfig {
                    duration_s: 1e-10,
                    ..cfg
                }
            )
            .is_err()
        );
    });
}

#[test]
fn g0_native_source_biot_budget_and_cancellation_refusals_publish_no_candidate() {
    let hot = caloric(0x94, 0.0, &[(0.0, 450.0, 0.0), (200.0, 500.0, 0.0)]);
    let pcm = caloric(
        0x95,
        0.0,
        &[
            (0.0, 300.0, 0.0),
            (100.0, 350.0, 0.0),
            (300.0, 350.0, 1.0),
            (900.0, 500.0, 1.0),
        ],
    );
    let a = pair_body(&hot, 2.0);
    let b = pair_body(&pcm, 1.0);
    let r = contact_resistance();
    let high_biot =
        LumpedEnthalpyBody::try_new("spatial body", 2.0, 0.01, 0.0, 0.0, 1.0, 1.0, &hot).unwrap();
    with_cx(|cx| {
        assert!(solve_lumped_pair_exchange(cx, &a, &b, &r, controls([200.0, 0.0], 5.0)).is_err());
        assert!(
            solve_lumped_pair_exchange(cx, &high_biot, &b, &r, controls([200.0, 0.0], 0.1))
                .is_err()
        );
        assert!(solve_lumped_pair_exchange(cx, &a, &b, &r, controls([200.0, 0.0], -1.0)).is_err());
        let zero = solve_lumped_pair_exchange(cx, &a, &b, &r, controls([200.0, 0.0], 0.0)).unwrap();
        assert_eq!(zero.caloric.states, zero.caloric.initial);
        let mut small_budget = controls([200.0, 0.0], 0.1);
        small_budget.integration.maximum_chart_intervals = 1;
        assert!(solve_lumped_pair_exchange(cx, &a, &b, &r, small_budget).is_err());
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_lumped_pair_exchange(cx, &a, &b, &r, controls([200.0, 0.0], 0.1)),
            Err(ConductionError::Cancelled { .. }),
        ));
    });
}
