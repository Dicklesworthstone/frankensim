//! Independent two-cell heat/species roots and whole-window thermal transport.

use super::adaptive::{BinaryNonisothermalAdvanceError, BinaryNonisothermalTimeControl};
use super::*;
use crate::miscibility::adaptive::BinaryMiscibilityTimeControl;
use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
use fs_qty::AmountConcentration;
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const WORKSPACE: usize = 1024 * 1024;

fn bracket() -> [Temperature; 2] {
    [Temperature::new(250.0), Temperature::new(700.0)]
}

fn prepared(
    chi: f64,
    kappa: f64,
    d: f64,
    cp: [f64; 2],
    area: f64,
) -> (BinaryMiscibilityNetwork, BinaryLiquidCaloricModelV1) {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        chi,
        0.0,
        bracket()[0],
        bracket()[1],
    )
    .unwrap();
    let sources = ["A", "B"]
        .into_iter()
        .zip(cp)
        .map(|(name, capacity)| {
            PureLiquidHeatCapacityV1::new(
                SpeciesId::new(name).unwrap(),
                Temperature::new(300.0),
                bracket()
                    .into_iter()
                    .map(|t| (t, MolarThermalQuantityV1::new(capacity)))
                    .collect(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let caloric =
        BinaryLiquidCaloricModelV1::new(pair.clone(), sources.try_into().unwrap()).unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        kappa,
    )
    .unwrap();
    let network = BinaryMiscibilityNetwork::new(
        model,
        vec![BinaryMiscibilityCell { volume_m3: 0.001 }; 2],
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: area,
            distance_m: 1.0,
            diffusivity_m2_s: d,
        }],
        WORKSPACE,
        || false,
    )
    .unwrap();
    (network, caloric)
}

fn field(
    cp: [f64; 2],
    d: f64,
    k: f64,
    q: f64,
    x: [f64; 2],
    t: [f64; 2],
) -> BinaryNonisothermalField {
    let (network, caloric) = prepared(0.0, 0.0, d, cp, 1.0);
    let state = network
        .state_from_mole_fractions(x.map(|v| [v, 1.0 - v]).to_vec(), || false)
        .unwrap();
    network
        .bind_nonisothermal_field(
            &state,
            caloric,
            t.map(Temperature::new).to_vec(),
            vec![BinaryNonisothermalFace {
                thermal_conductivity_w_m_k: k,
                heat_of_transport_j_per_mol: q,
            }],
            bracket(),
            WORKSPACE,
            || false,
        )
        .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = 0.5 * lo + 0.5 * hi;
        if mid == lo || mid == hi {
            break;
        }
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * lo + 0.5 * hi
}

fn time(duration: f64, temperature_tolerance: f64) -> BinaryNonisothermalTimeControl {
    BinaryNonisothermalTimeControl {
        composition: BinaryMiscibilityTimeControl {
            initial_step_s: duration,
            minimum_step_s: 1e-5,
            maximum_step_s: duration,
            absolute_mole_fraction_tolerance: 1.0,
            relative_mole_fraction_tolerance: 0.0,
            maximum_attempts: 512,
        },
        absolute_temperature_tolerance_k: temperature_tolerance,
        relative_temperature_tolerance: 0.0,
    }
}

fn conserved(initial: &BinaryNonisothermalField, final_state: &BinaryNonisothermalField) {
    for species in 0..2 {
        near(
            final_state
                .component_moles()
                .iter()
                .map(|n| n[species])
                .sum(),
            initial.component_moles().iter().map(|n| n[species]).sum(),
            2e-13,
        );
    }
    for amounts in final_state.component_moles() {
        assert!(amounts.iter().all(|&n| n > 0.0));
        near(amounts[0] + amounts[1], 1.0, 2e-13);
    }
    let mut heat = 0.0;
    for cell in 0..2 {
        let old = initial
            .caloric_model()
            .homogeneous(
                initial.temperatures()[cell],
                initial.component_moles()[cell],
            )
            .unwrap();
        let new = initial
            .caloric_model()
            .homogeneous(
                final_state.temperatures()[cell],
                final_state.component_moles()[cell],
            )
            .unwrap();
        heat += initial
            .caloric_model()
            .homogeneous_transport_change(&old, &new)
            .unwrap()
            .enthalpy_change_j;
    }
    near(heat, 0.0, 1e-8);
}

#[test]
fn g1_fourier_entropy_mobility_matches_independent_implicit_root_without_overshoot() {
    let initial = field([20.0, 40.0], 0.0, 1.0, 0.0, [0.5; 2], [600.0, 300.0]);
    let controls = BinaryNonisothermalSolve::default();
    // Each cell has N=1 and Cp=30. The independent equation is cubic in e:
    // e (600-e/30)(300+e/30) = dt*600*300*(300-2e/30).
    for dt in [0.5, 10_000.0] {
        let expected = root(0.0, 4500.0, |e| {
            e - dt
                * 600.0
                * 300.0
                * ((300.0 - 2.0 * e / 30.0) / ((600.0 - e / 30.0) * (300.0 + e / 30.0)))
        });
        let candidate = initial.step_candidate(dt, controls, || false).unwrap();
        let result = candidate.field();
        let [tl, tr] = [
            result.temperatures()[0].value(),
            result.temperatures()[1].value(),
        ];
        near(candidate.step().transferred_energy_j[0], expected, 2e-8);
        near(tl, 600.0 - expected / 30.0, 1e-9);
        near(tr, 300.0 + expected / 30.0, 1e-9);
        assert!(450.0 <= tl && tl < 600.0 && 300.0 < tr && tr <= 450.0);
        assert_eq!(result.component_moles(), initial.component_moles());
        assert_eq!(candidate.step().transferred_a_mol, [0.0]);
        let entropy = 30.0 * (tl / 600.0).ln() + 30.0 * (tr / 300.0).ln();
        near(candidate.step().entropy_change_j_per_k, entropy, 1e-10);
        assert!(entropy > 0.0);
        conserved(&initial, result);
        assert_eq!(result.elapsed_s(), dt);
    }
}

#[test]
fn g1_unequal_pure_capacities_do_not_invent_thermal_diffusion_without_soret() {
    let initial = field([20.0, 40.0], 1e-3, 0.0, 0.0, [0.25; 2], [300.0, 650.0]);
    let candidate = initial
        .step_candidate(1000.0, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(
        candidate.field().component_moles(),
        initial.component_moles()
    );
    assert_eq!(candidate.field().temperatures(), initial.temperatures());
    assert_eq!(candidate.step().transferred_a_mol, [0.0]);
    assert_eq!(candidate.step().transferred_energy_j, [0.0]);
    assert_eq!(candidate.step().entropy_change_j_per_k, 0.0);
    assert_eq!(candidate.step().energy_residual_j, 0.0);
}

#[test]
fn g1_soret_counterdiffusion_carries_original_enthalpy_and_reciprocal_dufour_heat() {
    let initial = field([20.0, 40.0], 1e-3, 0.0, 4000.0, [0.4; 2], [350.0, 550.0]);
    let dt = 0.2;
    let candidate = initial
        .step_candidate(dt, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let result = candidate.field();
    let q = candidate.step().transferred_a_mol[0];
    assert!(
        q < 0.0,
        "positive heat of transport sends A toward the cold cell"
    );
    for (cell, sign) in [(0, -1.0), (1, 1.0)] {
        near(
            result.component_moles()[cell][0] - initial.component_moles()[cell][0],
            sign * q,
            2e-14,
        );
        near(
            result.component_moles()[cell][1] - initial.component_moles()[cell][1],
            -sign * q,
            2e-14,
        );
    }
    // Independent constant-Cp Gibbs/enthalpy divided difference, datum 300 K.
    let divided_h = |tl: f64, tr: f64| -20.0 * (tl * tr * (tr / tl).ln() / (tr - tl) - 300.0);
    let transported_h = divided_h(350.0, 550.0) + 4000.0;
    near(
        candidate.step().transferred_energy_j[0],
        transported_h * q,
        5e-10,
    );
    let [tl, tr] = [
        result.temperatures()[0].value(),
        result.temperatures()[1].value(),
    ];
    let [left, right] = [result.component_moles()[0], result.component_moles()[1]];
    near(
        (20.0 * left[0] + 40.0 * left[1]) * (tl - 300.0) - 32.0 * 50.0,
        -transported_h * q,
        5e-9,
    );
    near(
        (20.0 * right[0] + 40.0 * right[1]) * (tr - 300.0) - 32.0 * 250.0,
        transported_h * q,
        5e-9,
    );
    let chemical = R * ((right[0] / right[1]).ln() - (left[0] / left[1]).ln());
    let force = -chemical + (transported_h - divided_h(tl, tr)) * (1.0 / tr - 1.0 / tl);
    near(q, dt * (0.4 * 0.6 / R) * force, 2e-12);
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
    conserved(&initial, result);
}

#[test]
fn g1_isothermal_ideal_pair_recovers_an_independent_implicit_log_force_root() {
    let initial = field([20.0; 2], 1e-3, 0.0, 0.0, [0.8, 0.2], [450.0; 2]);
    let dt = 1.7;
    // cD=1, arithmetic face x_A=x_B=1/2, one mole in each cell.
    // Hence q=dt/2 * ln[(.8-q)/(.2+q)] with a unique root in (0,.3).
    let expected = root(0.0, 0.3, |q| q - 0.5 * dt * ((0.8 - q) / (0.2 + q)).ln());
    let candidate = initial
        .step_candidate(dt, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let result = candidate.field();
    near(candidate.step().transferred_a_mol[0], expected, 2e-12);
    near(result.mole_fractions()[0][0], 0.8 - expected, 2e-12);
    near(result.mole_fractions()[1][0], 0.2 + expected, 2e-12);
    for t in result.temperatures() {
        near(t.value(), 450.0, 1e-9);
    }
    assert_eq!(candidate.step().transferred_energy_j, [0.0]);
    let mixing_s = |state: &BinaryNonisothermalField| {
        -R * state
            .component_moles()
            .iter()
            .flat_map(|n| n.iter())
            .map(|&n| n * n.ln())
            .sum::<f64>()
    };
    near(
        candidate.step().entropy_change_j_per_k,
        mixing_s(result) - mixing_s(&initial),
        1e-10,
    );
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
    assert!(result.mole_fractions()[0][0] > 0.5 && result.mole_fractions()[1][0] < 0.5);
    conserved(&initial, result);
}

#[test]
fn g2_adaptive_temperature_error_accepts_two_halves_and_refines_pure_heat() {
    let initial = field([20.0, 40.0], 0.0, 1.0, 0.0, [0.5; 2], [600.0, 300.0]);
    let controls = BinaryNonisothermalSolve::default();
    let first = initial.step_candidate(1.0, controls, || false).unwrap();
    let second = first
        .field()
        .step_candidate(1.0, controls, || false)
        .unwrap();
    let loose = initial
        .advance_candidate(2.0, time(2.0, 1000.0), controls, || false)
        .unwrap();
    assert_eq!(loose.field(), second.field());
    assert_eq!(loose.advance().accepted_intervals, 1);
    assert_eq!(loose.advance().rejected_intervals, 0);
    near(
        loose.advance().transferred_energy_j[0],
        first.step().transferred_energy_j[0] + second.step().transferred_energy_j[0],
        1e-10,
    );
    // Composition controls are unchanged; only the Kelvin allowance tightens.
    let tight = initial
        .advance_candidate(2.0, time(2.0, 0.01), controls, || false)
        .unwrap();
    assert!(tight.advance().accepted_intervals > 1);
    assert!(tight.advance().rejected_intervals > 0);
    assert!(tight.advance().maximum_accepted_error_ratio <= 1.0);
    let exact_hot = 450.0 + 150.0 * (-2.0_f64 * 2.0 / 30.0).exp();
    assert!(
        (tight.field().temperatures()[0].value() - exact_hot).abs()
            < (loose.field().temperatures()[0].value() - exact_hot).abs()
    );
    for candidate in [&loose, &tight] {
        assert_eq!(candidate.field().elapsed_s(), 2.0);
        assert_eq!(
            candidate.field().component_moles(),
            initial.component_moles()
        );
        assert!(candidate.advance().entropy_change_j_per_k > 0.0);
        conserved(&initial, candidate.field());
    }
}

#[test]
fn g4_late_cancellation_and_exhausted_work_leave_the_accepted_field_untouched() {
    let initial = field([20.0, 40.0], 0.0, 1.0, 0.0, [0.5; 2], [600.0, 300.0]);
    let before = initial.clone();
    let controls = BinaryNonisothermalSolve::default();
    let mut polls = 0;
    initial
        .step_candidate(2.0, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 2);
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            initial
                .step_candidate(2.0, controls, || {
                    count += 1;
                    count == stop
                })
                .unwrap_err(),
            BinaryNonisothermalError::Cancelled
        );
        assert_eq!(initial, before);
    }
    let mut bounded = time(0.5, 1000.0);
    bounded.composition.maximum_attempts = 1;
    let error = initial
        .advance_candidate(2.0, bounded, controls, || false)
        .unwrap_err();
    assert!(matches!(
        error,
        BinaryNonisothermalAdvanceError::WorkLimit {
            attempted_intervals: 1,
            accepted_intervals: 1,
            elapsed_s: 0.5,
        }
    ));
    assert_eq!(initial, before);
}

#[test]
fn g3_binding_refuses_unowned_physics_absent_support_and_invalid_transport() {
    let bind = |chi, kappa, x: [f64; 2], k, q, area| {
        let (network, caloric) = prepared(chi, kappa, 0.0, [20.0, 40.0], area);
        let state = network
            .state_from_mole_fractions(x.map(|v| [v, 1.0 - v]).to_vec(), || false)
            .unwrap();
        network.bind_nonisothermal_field(
            &state,
            caloric,
            vec![Temperature::new(600.0), Temperature::new(300.0)],
            vec![BinaryNonisothermalFace {
                thermal_conductivity_w_m_k: k,
                heat_of_transport_j_per_mol: q,
            }],
            bracket(),
            WORKSPACE,
            || false,
        )
    };
    assert!(bind(0.0, 1e-3, [0.5; 2], 1.0, 0.0, 1.0).is_err());
    // The original spatial model is admissible at its 300 K reference:
    // chi=550/300<2. The offered thermal bracket nevertheless reaches
    // chi=550/250=2.2, which this new field binder must refuse.
    let (reference, source) = prepared(0.0, 0.0, 0.0, [20.0, 40.0], 1.0);
    let pair = BinaryRegularSolutionV1::new(
        source.thermodynamics().species().clone(),
        0.0,
        550.0,
        bracket()[0],
        bracket()[1],
    )
    .unwrap();
    let caloric =
        BinaryLiquidCaloricModelV1::new(pair.clone(), source.heat_capacities().clone()).unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        0.0,
    )
    .unwrap();
    let network = BinaryMiscibilityNetwork::new(
        model,
        reference.cells().to_vec(),
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_m2_s: 0.0,
        }],
        WORKSPACE,
        || false,
    )
    .unwrap();
    let state = network
        .state_from_mole_fractions(vec![[0.5; 2]; 2], || false)
        .unwrap();
    assert!(matches!(
        network.bind_nonisothermal_field(
            &state,
            caloric,
            vec![Temperature::new(600.0), Temperature::new(300.0)],
            vec![BinaryNonisothermalFace {
                thermal_conductivity_w_m_k: 1.0,
                heat_of_transport_j_per_mol: 0.0
            }],
            bracket(),
            WORKSPACE,
            || false
        ),
        Err(BinaryNonisothermalError::Caloric(
            BinaryCaloricErrorV1::Invalid(_)
        ))
    ));
    assert!(bind(0.0, 0.0, [0.0, 0.5], 1.0, 0.0, 1.0).is_err());
    assert!(bind(0.0, 0.0, [1.0, 0.5], 1.0, 0.0, 1.0).is_err());
    for (k, q) in [(-1.0, 0.0), (f64::INFINITY, 0.0), (1.0, f64::NAN)] {
        assert!(bind(0.0, 0.0, [0.5; 2], k, q, 1.0).is_err());
    }
    // A representable positive geometry must not silently close a positive
    // heat path whose entropy mobility underflows in the actual step.
    let tiny = bind(0.0, 0.0, [0.5; 2], 1e-100, 0.0, f64::from_bits(1)).unwrap();
    assert!(
        tiny.step_candidate(1.0, BinaryNonisothermalSolve::default(), || false)
            .is_err()
    );
}
