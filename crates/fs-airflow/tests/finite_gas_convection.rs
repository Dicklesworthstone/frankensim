//! G0/G1/G3: real material packs, original correlation and finite heat/motion solver.
#[path = "finite_gas_convection/source.rs"]
mod source;
use fs_airflow::graph::thermal::medium::convection::{
    ForcedConvectionSettings, evaluate_forced_from_material, finite_gas::*,
};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve, ThermalPistonError};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Length;

const DT: f64 = 0.00005;
fn surface() -> ForcedGasSurface {
    ForcedGasSurface { area_m2: 0.02, characteristic_length: Length::new(0.1), speed_m_s: 1.0,
        correlation: CorrelationId::FlatPlateLaminarAverage, settings: ForcedConvectionSettings::default() }
}
fn piston(x: f64, t: f64, p: f64) -> GasPiston {
    GasPiston::new(ClosedGasParcel::new(source::model(x), source::point(t, p), 0.001).unwrap(),
        1.0, 0.01, 0.0).unwrap()
}
fn input() -> PistonInput {
    PistonInput { duration_s: DT, heat_j: 0.0, exterior_pressure_pa: 100000.0, external_force_n: 0.0 }
}
fn controls() -> PistonSolve {
    PistonSolve { volume_bracket_m3: [0.0009, 0.0011],
        gas: GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 900.0,
            energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 },
        momentum_tolerance_n_s: 1e-11, energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-0.05, 0.05], source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.001, volume_m3: 1e-6, surface_area_m2: 0.02, duration_s: DT, internal_heat_j: 0.0 }
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a} != {b}"); }

#[test]
fn source_pressure_and_composition_reach_conductance_without_an_air_override() {
    let curve = source::body_curve();
    let wall = body(&curve, 260000.0).initial;
    let low = piston(0.75, 400.0, 100000.0);
    let high = piston(0.75, 400.0, 200000.0);
    let a = surface().resolve(low.gas().state(), wall).unwrap();
    let b = surface().resolve(high.gas().state(), wall).unwrap();
    near(b.conductance_w_k/a.conductance_w_k, 2.0_f64.sqrt(), 1e-10);
    let different = piston(0.25, 400.0, 100000.0);
    let c = surface().resolve(different.gas().state(), wall).unwrap();
    assert!((c.conductance_w_k-a.conductance_w_k).abs() > 1e-3);
    assert_eq!(a.source.gas, *low.gas().state());
    assert_eq!(a.source.gas.source_packs(), low.gas().model().source_packs());
    assert_eq!(a.source.gas.source_packs().len(), 2);
}

#[test]
fn actual_joint_solver_closes_body_gas_motion_energy_and_retains_midpoint_sources() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let thermal = body(&curve, 260000.0); // 500 K, liquid body
    let next = surface().finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).unwrap();
    let report = next.step();
    assert!(report.exchanged_heat_j > 0.0);
    assert!(next.thermal().state.temperature_k() < thermal.initial.temperature_k());
    assert!(next.piston().gas().state().state().temperature > 400.0);
    assert!(next.piston().velocity_m_s() > 0.0);
    let du = old.gas().mass_kg()*(next.piston().gas().state().specific_internal_energy_j_kg()
        - old.gas().state().specific_internal_energy_j_kg());
    let dh = thermal.mass_kg*(next.thermal().state.specific_enthalpy_j_kg()-thermal.initial.specific_enthalpy_j_kg());
    let dk = 0.5*old.mass_kg()*next.piston().velocity_m_s().powi(2);
    let outside = input().exterior_pressure_pa*(next.piston().gas().volume_m3()-old.gas().volume_m3());
    near(du+dh+dk+outside, 0.0, 1e-7);
    near(dh+report.exchanged_heat_j, 0.0, 1e-7);
    let sample = &next.thermal().report;
    assert_eq!(sample.gas, report.piston.midpoint);
    assert_eq!(sample.body, report.thermal_midpoint);
    let evaluated = sample.evaluation.as_ref().unwrap();
    near(report.conductance_w_k, evaluated.coefficient.value.value()*sample.area_m2, 1e-14);
    near(report.exchanged_heat_j, DT*report.conductance_w_k*
        (sample.body.temperature_k()-sample.gas.state().temperature), exchange().source_tolerance_j);
    assert_eq!(old.gas().state().state().temperature, 400.0);
}

#[test]
fn melting_uses_body_enthalpy_and_not_an_artificial_heat_capacity() {
    let curve = source::body_curve();
    let old = piston(0.75, 600.0, 100000.0);
    let thermal = body(&curve, 180000.0); // halfway through latent heat
    let next = surface().finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).unwrap();
    assert!(next.step().exchanged_heat_j < 0.0);
    near(next.thermal().state.temperature_k(), 450.0, 1e-12);
    assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
    near(thermal.mass_kg*60000.0*(next.thermal().state.liquid_mass_fraction()-0.5),
        -next.step().exchanged_heat_j, 1e-7);
    assert_eq!(next.piston().mass_kg(), old.mass_kg()); // no melted-piston inference
}

#[test]
fn zero_area_reuses_the_closed_solver_without_resolving_a_fictitious_correlation() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let thermal = body(&curve, 260000.0);
    let closed = ForcedGasSurface { area_m2: 0.0, speed_m_s: 0.0, ..surface() };
    let a = closed.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).unwrap();
    let b = old.thermal_candidate(input(), thermal, controls(), exchange(),
        |_, _| Ok::<_, std::convert::Infallible>(ConductanceSample {conductance_w_k: 0.0, source: ()}), || false).unwrap();
    assert_eq!(a.piston(), b.piston());
    assert_eq!(a.thermal().state, b.thermal().state);
    assert_eq!(a.step(), b.step());
    assert!(a.thermal().report.evaluation.is_none());
}

#[test]
fn late_cancellation_geometry_and_correlation_refusals_leave_all_inputs_unchanged() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 260000.0);
    for s in [ForcedGasSurface {area_m2: 0.03, ..surface()},
        ForcedGasSurface {speed_m_s: 1e20, ..surface()},
        ForcedGasSurface {characteristic_length: Length::new(0.0), ..surface()}] {
        assert!(s.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).is_err());
        assert_eq!(old, saved);
    }
    let mut total = 0;
    surface().finite_body_candidate(&old, input(), thermal, controls(), exchange(), || { total += 1; false }).unwrap();
    for stop in [1, total/2, total] {
        let mut calls = 0;
        assert!(matches!(surface().finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
            calls += 1; calls == stop
        }), Err(ThermalPistonError::Cancelled)));
        assert_eq!(old, saved);
        assert_eq!(thermal.initial, curve.state_at_specific_enthalpy(260000.0).unwrap());
    }
}

#[test]
fn repeated_heat_exchange_refreshes_state_and_material_substitution_changes_motion() {
    let run = |fraction| {
        let curve = source::body_curve();
        let mut state = piston(fraction, 400.0, 100000.0);
        let initial = state.clone();
        let mut h = 260000.0;
        let mut work = 0.0;
        let mut first_g = 0.0;
        let mut last_g = 0.0;
        for i in 0..32 {
            let next = surface().finite_body_candidate(&state, input(), body(&curve, h),
                controls(), exchange(), || false).unwrap();
            last_g = next.step().conductance_w_k;
            if i == 0 { first_g = last_g; }
            work += next.step().piston.exterior_pressure_work_j;
            let (piston, thermal, _) = next.into_parts();
            state = piston; h = thermal.state.specific_enthalpy_j_kg();
        }
        let du = state.gas().mass_kg()*(state.gas().state().specific_internal_energy_j_kg()
            - initial.gas().state().specific_internal_energy_j_kg());
        near(du+0.001*(h-260000.0)+0.5*state.mass_kg()*state.velocity_m_s().powi(2)+work, 0.0, 2e-6);
        assert!((last_g-first_g).abs() > 1e-9);
        assert_eq!(state.gas().mass_kg(), initial.gas().mass_kg());
        (state.velocity_m_s(), h, state.gas().state().state().sound_speed)
    };
    let a = run(0.75); let b = run(0.25);
    assert!((a.0-b.0).abs() > 1e-8);
    assert!((a.1-b.1).abs() > 1e-3);
    assert!((a.2-b.2).abs() > 100.0);
}

#[test]
fn heat_flow_direction_comes_from_actual_midpoints_not_an_authored_flag() {
    let curve = source::body_curve();
    let old = piston(0.25, 400.0, 100000.0);
    let settings = ForcedConvectionSettings { length_over_hydraulic_diameter: Some(20.0),
        direction: ThermalDirection::CoolingFluid, ..ForcedConvectionSettings::default() };
    let s = ForcedGasSurface { correlation: CorrelationId::DittusBoelter, speed_m_s: 50.0,
        settings, ..surface() };
    let hot = s.resolve(old.gas().state(), body(&curve, 260000.0).initial).unwrap();
    let cold = s.resolve(old.gas().state(), body(&curve, 50000.0).initial).unwrap();
    let expected = evaluate_forced_from_material(s.correlation, old.gas().state(), s.speed_m_s,
        s.characteristic_length, ForcedConvectionSettings {direction: ThermalDirection::HeatingFluid, ..settings}).unwrap();
    near(hot.conductance_w_k, expected.coefficient().value.value()*s.area_m2, 1e-12);
    near(hot.conductance_w_k/cold.conductance_w_k, old.gas().state().state().prandtl.powf(0.1), 1e-10);
}
