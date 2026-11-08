//! Original-source plug-flow calorics, spatial pinch bounds and finite work.
#![cfg(feature = "material")]

// Reuse source inputs only. Numerical references use independent PR fugacity
// roots, analytic h/s departures, Brent caloric roots and Gauss-Kronrod spatial
// quadrature, rather than the production PH owner or adaptive Simpson rule.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, model, source};
use fs_gas::chemistry::ThermalConductance;
use fs_gas::material::phase_heat_exchanger::{
    MixtureHeatExchangerError as PhysicalError, MixtureHeatExchangerInlet as Inlet,
    spatial::{
        FlowArrangement as Flow, SpatialHeatExchangerControl as Control,
        SpatialHeatExchangerError as E, SpatialHeatExchangerStep as Step,
        exchange_spatial_mixture_heat,
    },
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError, equilibrium::PhaseEnthalpyControls,
};
use fs_qty::{MassFlowRate, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected
    );
}
fn inlet(
    model: &ThermochemicalPengRobinsonModel,
    temperature: f64,
    mass_flow: f64,
    pressure: f64,
) -> Inlet<'_> {
    Inlet {
        model,
        temperature: Temperature::new(temperature),
        pressure: Pressure::new(pressure),
        mass_flow: MassFlowRate::new(mass_flow),
    }
}
fn controls(bracket: [f64; 2], temperature: [f64; 2]) -> Control {
    let phase = PhaseEnthalpyControls::new(temperature.map(Temperature::new));
    Control::new(bracket, phase, phase)
}
fn ordinary_controls() -> Control {
    controls([-10_000.0, 10_000.0], [200.0, 240.0])
}
fn ordinary(
    model: &ThermochemicalPengRobinsonModel,
    flow: Flow,
    control: Control,
    cancelled: impl FnMut() -> bool,
) -> Result<Step, E> {
    exchange_spatial_mixture_heat(
        inlet(model, 210.0, 0.5, 3e6),
        inlet(model, 230.0, 1.0, 3e6),
        ThermalConductance::new(500.0),
        flow,
        control,
        cancelled,
    )
}

#[test]
fn g1_both_flow_arrangements_match_independent_two_phase_spatial_reference() {
    let original = source();
    for (flow, q, temperature, vapor, entropy) in [
        (
            Flow::Counterflow,
            9132.052_736_402_26,
            [212.332_247_327_112_47, 228.852_572_625_836_61],
            [0.123_881_717_329_143_98, 0.388_035_855_485_273_25],
            3.444_378_532_143_218_8,
        ),
        (
            Flow::Parallel,
            9109.826_753_706_286,
            [212.326_452_499_550_07, 228.855_392_002_711_57],
            [0.123_784_432_939_797_8, 0.388_085_040_978_992_27],
            3.436_820_229_515_433_3,
        ),
    ] {
        let result = ordinary(&original, flow, ordinary_controls(), || false).unwrap();
        close(result.balance.duty_to_left.value(), q, 1e-3);
        close(result.balance.entropy_production_w_k, entropy, 1e-6);
        for (side, receipt) in [
            &result.balance.left_equilibrium,
            &result.balance.right_equilibrium,
        ]
        .into_iter()
        .enumerate()
        {
            close(receipt.state().temperature_k(), temperature[side], 3e-7);
            close(
                receipt.state().vapor_mass_fraction().unwrap(),
                vapor[side],
                2e-8,
            );
            assert_eq!(receipt.state().constitutive_law(), original.mechanics());
        }
        assert!(result.balance.energy_residual_w.abs() < 2e-5);
        assert!(
            result
                .balance
                .thermal_residual_w
                .iter()
                .all(|r| r.abs() < 2e-5)
        );
        assert!(result.minimum_temperature_gap_k > 0.0);
        let first = result.profile.first().unwrap();
        let last = result.profile.last().unwrap();
        assert_eq!(first.heat_fraction, 0.0);
        assert_eq!(last.heat_fraction, 1.0);
        assert_eq!(first.conductance_fraction, 0.0);
        assert_eq!(last.conductance_fraction, 1.0);
        assert_eq!(first.left_equilibrium.state().temperature_k(), 210.0);
        assert_eq!(last.left_equilibrium, result.balance.left_equilibrium);
        if flow == Flow::Counterflow {
            assert_eq!(last.right_equilibrium.state().temperature_k(), 230.0);
            assert_eq!(first.right_equilibrium, result.balance.right_equilibrium);
        } else {
            assert_eq!(first.right_equilibrium.state().temperature_k(), 230.0);
            assert_eq!(last.right_equilibrium, result.balance.right_equilibrium);
        }
        assert!(
            result
                .profile
                .windows(2)
                .all(|w| w[0].heat_fraction < w[1].heat_fraction
                    && w[0].conductance_fraction < w[1].conductance_fraction)
        );
        // The exchanged-heat coordinate is not mislabeled as axial/UA position.
        assert!(
            result
                .profile
                .iter()
                .any(|p| (p.heat_fraction - p.conductance_fraction).abs() > 1e-3)
        );
        for station in &result.profile {
            let expected_left = result.balance.duty_to_left.value() * station.heat_fraction;
            let fraction = if flow == Flow::Counterflow {
                1.0 - station.heat_fraction
            } else {
                station.heat_fraction
            };
            let expected_right = -result.balance.duty_to_left.value() * fraction;
            close(
                0.5 * station.left_equilibrium.accepted_increment_j_kg().unwrap(),
                expected_left,
                2e-6,
            );
            close(
                station.right_equilibrium.accepted_increment_j_kg().unwrap(),
                expected_right,
                2e-6,
            );
        }
    }
}

#[test]
fn g1_spatial_enthalpy_heating_crosses_the_original_bubble_boundary() {
    let original = source();
    let control = controls([0.0, 30_000.0], [200.0, 240.0]);
    let result = exchange_spatial_mixture_heat(
        inlet(&original, 205.0, 0.5, 3e6),
        inlet(&original, 235.0, 1.0, 3e6),
        ThermalConductance::new(1000.0),
        Flow::Counterflow,
        control,
        || false,
    )
    .unwrap();
    assert!(
        result.profile[0]
            .left_equilibrium
            .state()
            .phase()
            .two_phase()
            .is_none()
    );
    let left = result.balance.left_equilibrium.state();
    let right = result.balance.right_equilibrium.state();
    close(
        result.balance.duty_to_left.value(),
        25133.875_973_382_834,
        1e-3,
    );
    close(left.temperature_k(), 211.493_305_509_259_13, 3e-7);
    close(right.temperature_k(), 232.093_546_580_616_67, 3e-7);
    close(
        left.vapor_mass_fraction().unwrap(),
        0.109_644_050_201_974_85,
        3e-8,
    );
    close(
        right.vapor_mass_fraction().unwrap(),
        0.447_157_022_983_650_23,
        3e-8,
    );
    close(
        result.balance.entropy_production_w_k,
        13.003_111_853_315_204,
        2e-6,
    );
    assert!(
        result
            .profile
            .iter()
            .any(|p| p.left_equilibrium.state().phase().two_phase().is_some())
    );
    assert!(result.panels > 1);
}

#[test]
fn g1_dilute_gas_recovers_both_analytic_effectiveness_ntu_limits() {
    let original = source();
    let cp = 4.0 * R / MIXTURE_MASS;
    let c_min = 0.5 * cp;
    let ratio = 0.5;
    let ntu = 500.0 / c_min;
    for (flow, reference) in [
        (Flow::Counterflow, 32681.363_552_826_086),
        (Flow::Parallel, 31085.670_282_820_3),
    ] {
        let effectiveness = if flow == Flow::Counterflow {
            let factor = (-ntu * (1.0 - ratio)).exp();
            (1.0 - factor) / (1.0 - ratio * factor)
        } else {
            (1.0 - (-ntu * (1.0 + ratio)).exp()) / (1.0 + ratio)
        };
        let result = exchange_spatial_mixture_heat(
            inlet(&original, 300.0, 0.5, 1.0),
            inlet(&original, 400.0, 1.0, 1.0),
            ThermalConductance::new(500.0),
            flow,
            controls([0.0, 60_000.0], [295.0, 405.0]),
            || false,
        )
        .unwrap();
        close(result.balance.duty_to_left.value(), reference, 1e-4);
        // Finite 1-Pa PR departures shift duty by 0.486/0.552 mW.
        close(
            result.balance.duty_to_left.value(),
            effectiveness * c_min * 100.0,
            1e-3,
        );
        assert!(
            result
                .balance
                .left_equilibrium
                .state()
                .phase()
                .two_phase()
                .is_none()
        );
        assert!(
            result
                .balance
                .right_equilibrium
                .state()
                .phase()
                .two_phase()
                .is_none()
        );
        assert!(result.balance.entropy_production_w_k > 0.0);
    }
}

#[test]
fn g1_counterflow_accepts_outlet_crossing_and_bounds_a_near_pinch_upper_trial() {
    let original = source();
    let cp = 4.0 * R / MIXTURE_MASS;
    let mut control = controls([0.0, 100.0 * cp], [295.0, 405.0]);
    control.absolute_power_tolerance_w = 1e-3;
    control.maximum_refinement_depth = 12;
    control.left_equilibrium.enthalpy.enthalpy_tolerance_j_kg = 1e-7;
    control.right_equilibrium.enthalpy.enthalpy_tolerance_j_kg = 1e-7;
    // The ideal limiting duty gives a tiny POSITIVE finite-PR gap. Its monotone
    // panel gap proof would need millions of intervals. Conductance lower
    // bounds must reject that upper trial early while admitting the real root.
    let solve = |control| {
        exchange_spatial_mixture_heat(
            inlet(&original, 300.0, 1.0, 1.0),
            inlet(&original, 400.0, 1.0, 1.0),
            ThermalConductance::new(100.0 * cp),
            Flow::Counterflow,
            control,
            || false,
        )
    };
    let result = solve(control).unwrap();
    close(
        result.balance.duty_to_left.value(),
        cp * 100.0 * 100.0 / 101.0,
        0.03,
    );
    close(
        result.balance.left_equilibrium.state().temperature_k(),
        399.009_900_990_099,
        3e-5,
    );
    close(
        result.balance.right_equilibrium.state().temperature_k(),
        300.990_099_009_901,
        3e-5,
    );
    assert!(
        result.balance.left_equilibrium.state().temperature_k()
            > result.balance.right_equilibrium.state().temperature_k() + 98.0
    );
    assert!(
        result
            .profile
            .iter()
            .all(|p| p.right_equilibrium.state().temperature_k()
                > p.left_equilibrium.state().temperature_k())
    );
    assert!(result.balance.entropy_production_w_k > 0.0);
    control.maximum_panels = 1;
    assert!(
        matches!(solve(control), Err(E::QuadratureLimit { gap_lower_bound_k, .. }) if gap_lower_bound_k < 0.0)
    );
}

#[test]
fn g3_independent_stream_reference_datums_leave_the_entire_profile_unchanged() {
    let original = source();
    let a = ordinary(&original, Flow::Counterflow, ordinary_controls(), || false).unwrap();
    let shifted_left = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let shifted_right = model([-7e23, 2e22], [-5e18, 3e17], [100.0, 450.0]);
    let b = exchange_spatial_mixture_heat(
        inlet(&shifted_left, 210.0, 0.5, 3e6),
        inlet(&shifted_right, 230.0, 1.0, 3e6),
        ThermalConductance::new(500.0),
        Flow::Counterflow,
        ordinary_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(a.balance.duty_to_left, b.balance.duty_to_left);
    assert_eq!(
        a.balance.enthalpy_flow_change_w,
        b.balance.enthalpy_flow_change_w
    );
    assert_eq!(
        a.balance.entropy_flow_change_w_k,
        b.balance.entropy_flow_change_w_k
    );
    assert_eq!(a.balance.pt_evaluations, b.balance.pt_evaluations);
    assert_eq!(
        a.balance.eos_evaluation_upper_bound,
        b.balance.eos_evaluation_upper_bound
    );
    assert_eq!(
        a.quadrature_error_estimate_w_k,
        b.quadrature_error_estimate_w_k
    );
    assert_eq!(a.profile.len(), b.profile.len());
    for (plain, shifted) in a.profile.iter().zip(&b.profile) {
        assert_eq!(plain.heat_fraction, shifted.heat_fraction);
        assert_eq!(plain.conductance_fraction, shifted.conductance_fraction);
        for (p, s) in [
            (&plain.left_equilibrium, &shifted.left_equilibrium),
            (&plain.right_equilibrium, &shifted.right_equilibrium),
        ] {
            assert_eq!(p.state().temperature_k(), s.state().temperature_k());
            assert_eq!(p.state().phase(), s.state().phase());
        }
    }
    let reverse = exchange_spatial_mixture_heat(
        inlet(&original, 230.0, 1.0, 3e6),
        inlet(&original, 210.0, 0.5, 3e6),
        ThermalConductance::new(500.0),
        Flow::Counterflow,
        ordinary_controls(),
        || false,
    )
    .unwrap();
    close(
        reverse.balance.duty_to_left.value(),
        -a.balance.duty_to_left.value(),
        1e-4,
    );
    close(
        reverse.balance.left_equilibrium.state().temperature_k(),
        a.balance.right_equilibrium.state().temperature_k(),
        3e-7,
    );
    close(
        reverse.balance.entropy_production_w_k,
        a.balance.entropy_production_w_k,
        1e-6,
    );
}

#[test]
fn g0_exact_rest_reuses_both_preparations_with_no_endpoint_pt_work() {
    let original = source();
    let mut control = ordinary_controls();
    control.maximum_pt_evaluations = 2;
    control.maximum_duty_evaluations = 1;
    for flow in [Flow::Parallel, Flow::Counterflow] {
        for (temperature, ua) in [(230.0, 0.0), (210.0, 500.0)] {
            let result = exchange_spatial_mixture_heat(
                inlet(&original, 210.0, 0.5, 3e6),
                inlet(&original, temperature, 1.0, 3e6),
                ThermalConductance::new(ua),
                flow,
                control,
                || false,
            )
            .unwrap();
            assert_eq!(result.balance.duty_to_left.value(), 0.0);
            assert_eq!(result.balance.enthalpy_flow_change_w, [0.0; 2]);
            assert_eq!(result.balance.entropy_flow_change_w_k, [0.0; 2]);
            assert_eq!(result.balance.pt_evaluations, 2);
            assert_eq!(result.balance.left_equilibrium.pt_evaluations(), 0);
            assert_eq!(result.balance.right_equilibrium.pt_evaluations(), 0);
            assert_eq!(result.balance.duty_evaluations, 1);
            assert_eq!(result.caloric_queries, 2);
            assert_eq!(result.panels, 0);
            assert_eq!(result.profile.len(), 2);
            assert_eq!(result.estimated_duty_error_w, 0.0);
        }
    }
}

#[test]
fn g0_total_work_limits_include_preparations_all_trials_and_spatial_samples() {
    let original = source();
    let accepted = ordinary(&original, Flow::Counterflow, ordinary_controls(), || false).unwrap();
    let mut control = ordinary_controls();
    control.maximum_pt_evaluations = accepted.balance.pt_evaluations;
    control.maximum_duty_evaluations = accepted.balance.duty_evaluations;
    assert_eq!(
        ordinary(&original, Flow::Counterflow, control, || false).unwrap(),
        accepted
    );
    control.maximum_pt_evaluations -= 1;
    match ordinary(&original, Flow::Counterflow, control, || false).unwrap_err() {
        E::Exchanger(PhysicalError::Equilibrium {
            completed_pt_evaluations,
            source: ThermochemicalFlashError::WorkLimit { flash_evaluations },
            ..
        }) => {
            assert_eq!(
                completed_pt_evaluations + flash_evaluations,
                control.maximum_pt_evaluations
            );
        }
        other => panic!("unexpected whole-spatial PT refusal: {other:?}"),
    }
    control.maximum_pt_evaluations = accepted.balance.pt_evaluations;
    control.maximum_duty_evaluations -= 1;
    assert!(
        matches!(ordinary(&original, Flow::Counterflow, control, || false),
        Err(E::Exchanger(PhysicalError::WorkLimit { duty_evaluations, .. }))
        if duty_evaluations == control.maximum_duty_evaluations)
    );
    control = ordinary_controls();
    control.left_equilibrium.enthalpy.temperature_bracket =
        [Temperature::new(210.0), Temperature::new(211.0)];
    assert!(matches!(
        ordinary(&original, Flow::Counterflow, control, || false),
        Err(E::Exchanger(PhysicalError::Equilibrium {
            side: 0,
            source: ThermochemicalFlashError::OutsideEnthalpyBracket,
            ..
        }))
    ));
}

#[test]
fn g1_tighter_spatial_accuracy_converges_to_the_independent_continuum_duty() {
    let original = source();
    let mut coarse = controls([0.0, 60_000.0], [295.0, 405.0]);
    coarse.absolute_power_tolerance_w = 1.0;
    coarse.relative_power_tolerance = 0.0;
    coarse.absolute_conductance_tolerance_w_k = 0.1;
    coarse.relative_conductance_tolerance = 0.0;
    let fine = controls([0.0, 60_000.0], [295.0, 405.0]);
    let solve = |control| {
        exchange_spatial_mixture_heat(
            inlet(&original, 300.0, 0.5, 1.0),
            inlet(&original, 400.0, 1.0, 1.0),
            ThermalConductance::new(500.0),
            Flow::Counterflow,
            control,
            || false,
        )
        .unwrap()
    };
    let a = solve(coarse);
    let b = solve(fine);
    let reference = 32681.363_552_826_086;
    assert!(
        (b.balance.duty_to_left.value() - reference).abs()
            < (a.balance.duty_to_left.value() - reference).abs()
    );
    assert!(b.panels > a.panels);
    assert!(b.quadrature_error_estimate_w_k < a.quadrature_error_estimate_w_k);
    close(b.balance.duty_to_left.value(), reference, 1e-4);
}

#[test]
fn g4_final_cancellation_returns_no_partial_profile_and_replay_is_exact() {
    let original = source();
    let mut polls = 0;
    let accepted = ordinary(&original, Flow::Counterflow, ordinary_controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert_eq!(
        ordinary(&original, Flow::Counterflow, ordinary_controls(), || {
            count += 1;
            count == polls
        })
        .unwrap_err(),
        E::Cancelled
    );
    assert_eq!(count, polls);
    assert_eq!(
        ordinary(&original, Flow::Counterflow, ordinary_controls(), || false).unwrap(),
        accepted
    );
}
