//! Manufactured original calorics and independent two-body trajectory references.

use super::*;
use crate::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};
use fs_blake3::ContentHash;

fn chart_knots(datum: f64, points: &[(f64, f64, f64)]) -> Vec<EnthalpyPhaseKnot> {
    points
        .iter()
        .map(|&(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: datum + h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 1000.0 - 200.0 * liquid,
        })
        .collect()
}

fn curve(
    id: u8,
    datum: f64,
    points: &[(f64, f64, f64)],
    cp: &[Option<[f64; 2]>],
) -> EquilibriumEnthalpyPhaseCurve {
    let mut law =
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([id; 32]), chart_knots(datum, points))
            .unwrap();
    law.retain_sensible_heat_capacities(cp.to_vec()).unwrap();
    law
}

// An all-solid chart declares its one phase. The two-phase ingress admits only
// a chart that begins fully solid and ends fully liquid.
fn solid_curve(
    id: u8,
    datum: f64,
    points: &[(f64, f64, f64)],
    cp: &[Option<[f64; 2]>],
) -> EquilibriumEnthalpyPhaseCurve {
    let mut law = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([id; 32]),
        SolidLiquidPhase::Solid,
        chart_knots(datum, points),
    )
    .unwrap();
    law.retain_sensible_heat_capacities(cp.to_vec()).unwrap();
    law
}

fn melting(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    curve(
        0xe1,
        datum,
        &[
            (0.0, 300.0, 0.0),
            (100.0, 350.0, 0.0),
            (300.0, 350.0, 1.0),
            (900.0, 500.0, 1.0),
        ],
        &[Some([2.0; 2]), None, Some([4.0; 2])],
    )
}
fn sensible(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    solid_curve(
        0xe2,
        datum,
        &[(0.0, 200.0, 0.0), (1200.0, 500.0, 0.0)],
        &[Some([4.0; 2])],
    )
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64, mass: f64) -> FinitePhaseBody<'_> {
    FinitePhaseBody {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: mass,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_three_stage_melting_and_freezing_follow_independent_finite_body_exponentials() {
    let hot = sensible(0.0);
    let pcm = melting(0.0);
    let start_melt = 1.6 * (200.0_f64 / 137.5).ln();
    let end_melt = start_melt + 8.0 * (137.5_f64 / 112.5).ln();
    for x in [start_melt, 1.0, end_melt, 5.0] {
        let q = if x <= start_melt {
            320.0 * -(-0.625 * x).exp_m1()
        } else if x <= end_melt {
            100.0 + 1100.0 * -(-0.125 * (x - start_melt)).exp_m1()
        } else {
            300.0 + 300.0 * -(-0.375 * (x - end_melt)).exp_m1()
        };
        let cold_t = if q < 100.0 {
            300.0 + q / 2.0
        } else if q <= 300.0 {
            350.0
        } else {
            350.0 + (q - 300.0) / 4.0
        };
        let hot_t = 500.0 - q / 8.0;
        let result = exchange(
            body(&hot, 1200.0, 2.0),
            body(&pcm, 0.0, 1.0),
            x,
            FinitePhaseExchangeControls::default(),
        )
        .unwrap();
        near(result.heat_first_to_second_j, q, 2e-10);
        near(result.states[0].temperature_k(), hot_t, 2e-12);
        near(result.states[1].temperature_k(), cold_t, 2e-12);
        near(
            result.states[1].liquid_mass_fraction(),
            ((q - 100.0) / 200.0).clamp(0.0, 1.0),
            2e-12,
        );
        assert_eq!(result.root_iterations, 0);
        assert_eq!(result.quadrature_panels, 0);
        let cold_entropy = if q <= 100.0 {
            2.0 * (cold_t / 300.0).ln()
        } else if q <= 300.0 {
            2.0 * (350.0_f64 / 300.0).ln() + (q - 100.0) / 350.0
        } else {
            2.0 * (350.0_f64 / 300.0).ln() + 200.0 / 350.0 + 4.0 * (cold_t / 350.0).ln()
        };
        near(
            result.entropy_change_j_k[0],
            8.0 * (hot_t / 500.0).ln(),
            3e-13,
        );
        near(result.entropy_change_j_k[1], cold_entropy, 3e-13);
        assert!(result.total_entropy_change_j_k > 0.0);
        assert!(result.energy_residual_j.iter().all(|e| e.abs() < 2e-10));
    }
    // Cooling liquid transfers to a finite body initially at 250 K. Actual
    // receiver warming controls the freezing residence, not a fixed bath.
    let first = (8.0 / 3.0) * 4.0_f64.ln();
    let frozen = first + 8.0 * 2.0_f64.ln();
    let x = 10.0;
    let q = 600.0 + 40.0 * -(-0.625 * (x - frozen)).exp_m1();
    let result = exchange(
        body(&pcm, 700.0, 1.0),
        body(&hot, 200.0, 2.0),
        x,
        FinitePhaseExchangeControls::default(),
    )
    .unwrap();
    near(result.heat_first_to_second_j, q, 2e-10);
    near(
        result.states[0].temperature_k(),
        300.0 + (700.0 - q) / 2.0,
        2e-12,
    );
    near(result.states[1].temperature_k(), 250.0 + q / 8.0, 2e-12);
    assert_eq!(result.states[0].liquid_mass_fraction(), 0.0);
    assert!(result.states[0].temperature_k() > result.states[1].temperature_k());
}

#[test]
fn g3_relabeling_composition_weak_heat_and_stiff_relaxation_preserve_one_transfer() {
    let hot = sensible(0.0);
    let pcm = melting(0.0);
    let controls = FinitePhaseExchangeControls::default();
    let a = body(&hot, 1200.0, 2.0);
    let b = body(&pcm, 0.0, 1.0);
    let direct = exchange(a, b, 5.0, controls).unwrap();
    let reversed = exchange(b, a, 5.0, controls).unwrap();
    assert_eq!(
        reversed.heat_first_to_second_j,
        -direct.heat_first_to_second_j
    );
    assert_eq!(reversed.states, [direct.states[1], direct.states[0]]);
    let first = exchange(a, b, 1.0, controls).unwrap();
    let second = exchange(
        FinitePhaseBody {
            initial: first.states[0],
            ..a
        },
        FinitePhaseBody {
            initial: first.states[1],
            ..b
        },
        4.0,
        controls,
    )
    .unwrap();
    near(
        first.heat_first_to_second_j + second.heat_first_to_second_j,
        direct.heat_first_to_second_j,
        2e-10,
    );
    near(
        second.states[1].temperature_k(),
        direct.states[1].temperature_k(),
        2e-12,
    );
    let x = 1e-8;
    let weak = exchange(a, b, x, controls).unwrap();
    near(weak.heat_first_to_second_j, 200.0 * x - 62.5 * x * x, 1e-19);
    for x in [100.0, 1e250] {
        let stiff = exchange(a, b, x, controls).unwrap();
        near(stiff.heat_first_to_second_j, 600.0, 2e-10);
        near(stiff.states[0].temperature_k(), 425.0, 2e-12);
        near(stiff.states[1].temperature_k(), 425.0, 2e-12);
        assert!(stiff.states[0].temperature_k() >= stiff.states[1].temperature_k());
    }
    let equilibrium = exchange(
        body(&pcm, 150.0, 1.0),
        body(&pcm, 250.0, 2.0),
        1e250,
        controls,
    )
    .unwrap();
    assert_eq!(equilibrium.states, equilibrium.initial);
    assert_eq!(equilibrium.heat_first_to_second_j, 0.0);
    assert_eq!(equilibrium.entropy_change_j_k, [0.0; 2]);
}

fn nonlinear(
    datum: f64,
    split: bool,
) -> (EquilibriumEnthalpyPhaseCurve, EquilibriumEnthalpyPhaseCurve) {
    if split {
        (
            solid_curve(
                0xf1,
                datum,
                &[(0.0, 400.0, 0.0), (125.0, 450.0, 0.0), (300.0, 500.0, 0.0)],
                &[Some([2.0, 3.0]), Some([3.0, 4.0])],
            ),
            solid_curve(
                0xf2,
                -datum,
                &[(0.0, 300.0, 0.0), (175.0, 350.0, 0.0), (300.0, 400.0, 0.0)],
                &[Some([4.0, 3.0]), Some([3.0, 2.0])],
            ),
        )
    } else {
        (
            solid_curve(
                0xf1,
                datum,
                &[(0.0, 400.0, 0.0), (300.0, 500.0, 0.0)],
                &[Some([2.0, 4.0])],
            ),
            solid_curve(
                0xf2,
                -datum,
                &[(0.0, 300.0, 0.0), (300.0, 400.0, 0.0)],
                &[Some([4.0, 2.0])],
            ),
        )
    }
}

#[test]
fn g1_source_affine_cp_uses_the_original_nonlinear_trajectory_across_source_knots() {
    // Symmetric capacities Cp_hot=2+.02*(T-400), Cp_cold=2-.02*(T-400).
    // With T_hot=400+d, T_cold=400-d, X=ln(100/d)+.01*(100-d).
    // d=50 gives Q=2*50+.01*(100^2-50^2)=175 J, independently of inversion.
    let exposure = 2.0_f64.ln() + 0.5;
    for datum in [0.0, 2.0_f64.powi(35)] {
        for split in [false, true] {
            let (hot, cold) = nonlinear(datum, split);
            let result = exchange(
                body(&hot, datum + 300.0, 1.0),
                body(&cold, -datum, 1.0),
                exposure,
                FinitePhaseExchangeControls::default(),
            )
            .unwrap();
            near(result.heat_first_to_second_j, 175.0, 2e-6);
            near(result.states[0].temperature_k(), 450.0, 2e-6);
            near(result.states[1].temperature_k(), 350.0, 2e-6);
            assert!(result.root_iterations > 0);
            assert!(result.quadrature_panels > 0);
            assert!(result.heat_error_bound_j <= 1e-7 + 1e-10 * result.heat_first_to_second_j);
            assert!(
                result.exposure_residual_j_k.abs() + result.exposure_error_allowance_j_k
                    <= 1e-9 + 1e-10 * exposure
            );
            let entropy =
                -6.0 * (450.0_f64 / 500.0).ln() - 1.0 + 10.0 * (350.0_f64 / 300.0).ln() - 1.0;
            near(result.total_entropy_change_j_k, entropy, 3e-8);
        }
    }
}

#[test]
fn g1_finite_latent_donor_and_nonlinear_receiver_consume_both_original_inventories() {
    let donor = curve(0xe3, 0.0, &[(0.0, 550.0, 0.0), (1e6, 550.0, 1.0)], &[None]);
    let receiver = curve(
        0xe4,
        0.0,
        &[
            (0.0, 300.0, 0.0),
            (150_000.0, 400.0, 0.0),
            (170_000.0, 400.0, 1.0),
            (470_000.0, 500.0, 1.0),
        ],
        &[Some([1000.0, 2000.0]), None, Some([2000.0, 4000.0])],
    );
    // Independent Cp/(550-T) primitives and the one latent residence.
    let x = 3500.0 * (230.0_f64 / 150.0).ln() - 10.0 * 80.0
        + 20_000.0 / 150.0
        + 5000.0 * (150.0_f64 / 90.0).ln()
        - 20.0 * 60.0;
    let result = exchange(
        body(&donor, 800_000.0, 1.0),
        body(&receiver, 22_000.0, 1.0),
        x,
        FinitePhaseExchangeControls::default(),
    )
    .unwrap();
    near(result.heat_first_to_second_j, 304_000.0, 8e-5);
    near(result.states[1].temperature_k(), 460.0, 5e-8);
    assert_eq!(result.states[0].temperature_k(), 550.0);
    near(result.states[0].liquid_mass_fraction(), 0.496, 1e-10);
    near(result.entropy_change_j_k[0], -304_000.0 / 550.0, 2e-7);
    assert!(result.total_entropy_change_j_k > 0.0);
}

#[test]
fn g0_source_work_representation_identity_and_late_cancellation_remain_transactional() {
    let hot = sensible(0.0);
    let pcm = melting(0.0);
    let controls = FinitePhaseExchangeControls::default();
    let a = body(&hot, 1200.0, 2.0);
    let b = body(&pcm, 0.0, 1.0);
    let limited = solid_curve(
        0xe5,
        0.0,
        &[(0.0, 450.0, 0.0), (200.0, 500.0, 0.0)],
        &[Some([4.0; 2])],
    );
    assert!(matches!(
        exchange(body(&limited, 200.0, 2.0), b, 5.0, controls),
        Err(FinitePhaseExchangeError::OutsideSourceSupport),
    ));
    let unsupported_direction = body(&hot, 0.0, 2.0);
    assert_eq!(
        exchange(unsupported_direction, b, 0.0, controls)
            .unwrap()
            .heat_first_to_second_j,
        0.0
    );
    assert!(matches!(
        exchange(
            a,
            b,
            1.0,
            FinitePhaseExchangeControls {
                maximum_chart_intervals: 1,
                ..controls
            }
        ),
        Err(FinitePhaseExchangeError::WorkLimit)
    ));
    let (nonlinear_hot, nonlinear_cold) = nonlinear(0.0, false);
    for (panels, iterations, work_limit) in [(1, 96, true), (131_072, 1, false)] {
        let refused = exchange(
            body(&nonlinear_hot, 300.0, 1.0),
            body(&nonlinear_cold, 0.0, 1.0),
            2.0_f64.ln() + 0.5,
            FinitePhaseExchangeControls {
                maximum_quadrature_panels: panels,
                maximum_root_iterations: iterations,
                ..controls
            },
        );
        if work_limit {
            assert!(matches!(refused, Err(FinitePhaseExchangeError::WorkLimit)));
        } else {
            assert!(matches!(
                refused,
                Err(FinitePhaseExchangeError::NonConvergence)
            ));
        }
    }
    let shifted_hot = sensible(2.0_f64.powi(52));
    let shifted_pcm = melting(2.0_f64.powi(52));
    assert!(matches!(
        exchange(
            body(&shifted_hot, 2.0_f64.powi(52) + 1200.0, 2.0),
            body(&shifted_pcm, 2.0_f64.powi(52), 1.0),
            1e-10,
            controls,
        ),
        Err(FinitePhaseExchangeError::Phase(
            PhaseStateError::LostSpecificEnergyIncrement
        ))
    ));
    assert!(matches!(
        exchange(
            FinitePhaseBody {
                initial: b.initial,
                ..a
            },
            b,
            1.0,
            controls
        ),
        Err(FinitePhaseExchangeError::Phase(
            PhaseStateError::StateCurveMismatch
        ))
    ));
    let before = (hot.clone(), pcm.clone());
    let mut checkpoints = 0;
    assert!(matches!(
        exchange_with_checkpoint(a, b, 1.0, controls, || {
            checkpoints += 1;
            checkpoints == 10
        }),
        Err(FinitePhaseExchangeError::Cancelled)
    ));
    assert_eq!(checkpoints, 10);
    assert_eq!((hot, pcm), before);
}
