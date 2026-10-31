//! Dense gas, sourced finite-wall heat and native inertial compression together.
//! Synthetic sources verify coupling; no empirical dense-fluid calibration claim.

use super::layers::MaterialConductionLayers;
use super::layers::tests::{material, point};
use super::power::PowerSample;
use super::*;
use crate::cavity_phs::piston::nonideal_tests::{controls, drive, parcel};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use std::convert::Infallible;

fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-2.0, 2.0],
        source_tolerance_j: 1.0e-9,
        energy_tolerance_j: 1.0e-6,
        maximum_iterations: 128,
    }
}
fn input() -> PistonInput {
    PistonInput {
        exterior_pressure_pa: 5.2e6,
        external_force_n: 0.0,
        ..drive(0.0)
    }
}
fn body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    enthalpy: f64,
    duration: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(enthalpy).unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: duration,
    }
}
fn layers<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    context: &'a QueryPoint,
) -> MaterialConductionLayers<'a> {
    MaterialConductionLayers {
        material: card,
        context,
        temperature_axis: "T",
        liquid_fraction_axis: None,
        conductivity: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.02,
        gas_distance_m: 1.0e-5,
        body_distance_m: 1.0e-4,
    }
}

#[test]
fn g3_real_gas_compression_and_sourced_wall_heat_close_accepted_state_energy() {
    let (card, key, caloric) = material(Some(0.2), true, false);
    let context = point(300.0);
    let source = layers(&card, &key, &context);
    let mut thermal = body(caloric.curve(), 150_000.0, input().duration_s);
    let h0 = thermal.initial.specific_enthalpy_j_kg();
    let mut piston = GasPiston::new(parcel(), 2.0, 0.01, -0.03).unwrap();
    let initial = piston.clone();
    let mut transferred = 0.0;
    for _ in 0..40 {
        let candidate = piston
            .thermal_candidate(
                input(),
                thermal,
                controls(),
                exchange(),
                |gas, wall| source.resolve(gas, wall),
                || false,
            )
            .unwrap();
        let (next, accepted, report) = candidate.into_parts();
        assert!(report.source_residual_j.abs() <= exchange().source_tolerance_j);
        assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
        assert!(report.piston.gas.work_by_gas_j < 0.0);
        assert_eq!(accepted.external_heat_j, -report.exchanged_heat_j);
        assert_eq!(accepted.report.card_identity(), card.content_hash());
        assert_eq!(accepted.report.query_point().get("lot"), Some(&7.0));
        assert_eq!(
            accepted.report.query_point().get("T"),
            Some(&report.thermal_midpoint.temperature_k())
        );
        let transport = report.piston.midpoint.transport_state().unwrap();
        assert_eq!(
            transport.density_kg_m3(),
            report.piston.midpoint.density_kg_m3()
        );
        assert_eq!(
            transport.specific_heat_j_kg_k(),
            report.piston.midpoint.specific_heat_cp_j_kg_k()
        );
        let kb = accepted.report.property(key.name()).unwrap().value_si();
        let expected_g = source.area_m2
            / (source.gas_distance_m / transport.thermal_conductivity_w_m_k()
                + source.body_distance_m / kb);
        assert!((report.conductance_w_k / expected_g - 1.0).abs() < 1.0e-13);
        transferred += report.exchanged_heat_j;
        thermal.initial = accepted.state;
        piston = next;
    }
    // Independently recover each physical store from the accepted endpoints.
    // No accumulated solver energy-residual field enters this balance.
    let gas_energy = initial.gas().mass_kg()
        * initial
            .gas()
            .model()
            .specific_internal_energy_change_between_densities(
                initial.gas().density_kg_m3(),
                initial.gas().state().temperature_k(),
                piston.gas().density_kg_m3(),
                piston.gas().state().temperature_k(),
            )
            .unwrap();
    let kinetic =
        0.5 * piston.mass_kg() * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
    let wall_energy = thermal.mass_kg * (thermal.initial.specific_enthalpy_j_kg() - h0);
    let exterior =
        input().exterior_pressure_pa * (piston.gas().volume_m3() - initial.gas().volume_m3());
    assert!((gas_energy + kinetic + wall_energy + exterior).abs() < 2.0e-5);
    assert!((wall_energy + transferred).abs() < 1.0e-8);
    // Independent simultaneous solve of PR u(T,rho), momentum, midpoint
    // energy and the WMS/series-conduction law; no native feedback report
    // supplied these references. Tolerances allow accumulated solver budgets.
    assert!((piston.gas().volume_m3() - 0.0009966439609426744).abs() < 1.0e-10);
    assert!((piston.gas().state().temperature_k() - 300.55705895905845).abs() < 1.0e-5);
    assert!((piston.velocity_m_s() - -0.7908538377390196).abs() < 2.0e-5);
    assert!((transferred - 4.8579865501912085).abs() < 2.0e-5);
    assert!(transferred > 0.0 && thermal.initial.temperature_k() < 500.0);
    assert!(piston.gas().volume_m3() < initial.gas().volume_m3());
    assert!(piston.gas().state().temperature_k() > initial.gas().state().temperature_k());
    assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
    assert_eq!(piston.gas().model(), initial.gas().model());
    assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
}

#[test]
fn g3_dense_gas_receives_latent_wall_heat_from_the_actual_fraction_midpoint() {
    let (card, key, _) = material(Some(0.2), true, true);
    let knot = |h, t, x| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: x,
        bulk_density_kg_m3: 4000.0,
    };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            knot(0.0, 200.0, 0.0),
            knot(100_000.0, 350.0, 0.0),
            knot(160_000.0, 350.0, 1.0),
            knot(400_000.0, 800.0, 1.0),
        ],
    )
    .unwrap();
    let context = point(300.0);
    let mut source = layers(&card, &key, &context);
    source.liquid_fraction_axis = Some("melt");
    let input = PistonInput {
        duration_s: 2.0e-4,
        ..input()
    };
    let thermal = body(&curve, 130_000.0, input.duration_s);
    let piston = GasPiston::new(parcel(), 2.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let fixed = PistonSolve {
        volume_bracket_m3: [piston.gas().volume_m3(); 2],
        ..controls()
    };
    let candidate = piston
        .thermal_candidate(
            input,
            thermal,
            fixed,
            exchange(),
            |gas, wall| source.resolve(gas, wall),
            || false,
        )
        .unwrap();
    let report = candidate.step();
    assert!(report.exchanged_heat_j > 0.0);
    assert_eq!(candidate.thermal().state.temperature_k(), 350.0);
    assert!(
        candidate.thermal().state.liquid_mass_fraction() < thermal.initial.liquid_mass_fraction()
    );
    assert_eq!(report.thermal_midpoint.temperature_k(), 350.0);
    assert_eq!(
        candidate.thermal().report.query_point().get("melt"),
        Some(&report.thermal_midpoint.liquid_mass_fraction())
    );
    assert!(
        candidate.piston().gas().state().temperature_k() > piston.gas().state().temperature_k()
    );
    assert_eq!(
        candidate.piston().gas().density_kg_m3(),
        piston.gas().density_kg_m3()
    );
    assert_eq!(report.piston.gas.work_by_gas_j, 0.0);
    assert!(report.total_energy_residual_j.abs() < exchange().energy_tolerance_j);
}

#[test]
fn g3_state_dependent_power_uses_the_real_gas_and_body_energy_midpoints() {
    let (card, key, caloric) = material(Some(0.2), true, false);
    let context = point(300.0);
    let source = layers(&card, &key, &context);
    let thermal = body(caloric.curve(), 150_000.0, input().duration_s);
    let piston = GasPiston::new(parcel(), 2.0, 0.01, -0.03).unwrap();
    let power_controls = HeatExchangeSolve {
        heat_bracket_j: [0.0, 1.0],
        ..exchange()
    };
    let law = |gas: &ResolvedThermochemicalRealGasState, wall: EquilibriumPhaseState| {
        10.0 * (gas.pressure_pa() / 1.0e6) * (wall.temperature_k() - 280.0)
    };
    let candidate = piston
        .powered_thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            power_controls,
            |gas, wall| source.resolve(gas, wall),
            |gas, wall| {
                Ok(PowerSample {
                    power_w: law(gas, wall),
                    source: (gas.temperature_k(), wall.temperature_k()),
                })
            },
            || false,
        )
        .unwrap();
    let coupled = candidate.coupled();
    let report = coupled.step();
    let expected = input().duration_s * law(&report.piston.midpoint, report.thermal_midpoint);
    assert!(
        (candidate.step().applied_energy_j - expected).abs() <= power_controls.source_tolerance_j
    );
    assert_eq!(
        candidate.power().source,
        (
            report.piston.midpoint.temperature_k(),
            report.thermal_midpoint.temperature_k()
        )
    );
    let end = coupled.piston();
    let gas_energy = piston.gas().mass_kg()
        * piston
            .gas()
            .model()
            .specific_internal_energy_change_between_densities(
                piston.gas().density_kg_m3(),
                piston.gas().state().temperature_k(),
                end.gas().density_kg_m3(),
                end.gas().state().temperature_k(),
            )
            .unwrap();
    let wall_energy = thermal.mass_kg
        * (coupled.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let kinetic =
        0.5 * piston.mass_kg() * (end.velocity_m_s().powi(2) - piston.velocity_m_s().powi(2));
    let exterior =
        input().exterior_pressure_pa * (end.gas().volume_m3() - piston.gas().volume_m3());
    assert!(
        (gas_energy + wall_energy + kinetic + exterior - expected).abs()
            < power_controls.energy_tolerance_j
    );
}

#[test]
fn g4_real_gas_isolation_source_refusal_and_final_cancellation_preserve_owners() {
    let (card, key, caloric) = material(Some(0.2), true, false);
    let context = point(300.0);
    let source = layers(&card, &key, &context);
    let thermal = body(caloric.curve(), 150_000.0, input().duration_s);
    let mut piston = GasPiston::new(parcel(), 2.0, 0.01, -0.03).unwrap();
    piston.step(input(), controls(), || false).unwrap();
    let saved = piston.clone();
    let saved_body = thermal.initial;
    let isolated = piston
        .thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |_, _| {
                Ok::<_, Infallible>(ConductanceSample {
                    conductance_w_k: 0.0,
                    source: (),
                })
            },
            || false,
        )
        .unwrap();
    let mut direct = piston.clone();
    let direct_report = direct.step(input(), controls(), || false).unwrap();
    assert_eq!(isolated.piston(), &direct);
    assert_eq!(isolated.step().piston, direct_report);
    assert_eq!(isolated.thermal().state, thermal.initial);
    let wrong_context = context.clone().with("lot", 9.0).unwrap();
    let wrong = layers(&card, &key, &wrong_context);
    assert!(matches!(
        piston.thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |gas, wall| wrong.resolve(gas, wall),
            || false
        ),
        Err(ThermalPistonError::Interface(_))
    ));
    let mut polls = 0;
    piston
        .thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |gas, wall| source.resolve(gas, wall),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        piston.thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |gas, wall| source.resolve(gas, wall),
            || {
                current += 1;
                current == polls
            }
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(piston, saved);
    assert_eq!(thermal.initial, saved_body);
    use fs_ports::PortedModule;
    assert_eq!(
        piston.last_window_duration_s(),
        saved.last_window_duration_s()
    );
    assert_eq!(
        piston.last_window_transfers(),
        saved.last_window_transfers()
    );
}
