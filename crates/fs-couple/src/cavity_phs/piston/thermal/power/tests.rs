//! Actual public coupled owners; linear isolated storage supplies an analytic oracle.
use super::*;
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use std::convert::Infallible;

fn bounds() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(fs_blake3::hash_bytes(b"power feedback caloric fixture"),
        SolidLiquidPhase::Solid, [200.0, 800.0].into_iter().map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: 500.0*(t-200.0), temperature_k: t,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 4000.0,
        }).collect()).unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(100000.0).unwrap(),
        mass_kg: 0.001, volume_m3: 0.001/4000.0, surface_area_m2: 0.02,
        internal_heat_j: 0.0, duration_s: drive(0.0).duration_s }
}
fn link(g: f64) -> impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
    -> Result<ConductanceSample<EquilibriumPhaseState>, Infallible> {
    move |_, state| Ok(ConductanceSample { conductance_w_k: g, source: state })
}
fn constant(watts: f64) -> impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
    -> Result<PowerSample<EquilibriumPhaseState>, Infallible> {
    move |_, state| Ok(PowerSample { power_w: watts, source: state })
}
fn piston() -> GasPiston { GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap() }

#[test]
fn zero_and_constant_power_match_existing_prescribed_heat_owners() {
    let old = piston(); let curve = curve(); let thermal = body(&curve);
    let existing = old.thermal_candidate(drive(0.0), thermal, controls(), bounds(), link(1.0), ||false).unwrap();
    let zero = old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        link(1.0), constant(0.0), ||false).unwrap();
    assert_eq!(zero.coupled(), &existing);
    assert_eq!(zero.step().source_trials, 1);
    assert_eq!(zero.step().exchange_trials, existing.step().exchange_trials);
    assert_eq!(zero.step().material_trials, existing.step().material_trials);
    let powered = old.powered_thermal_candidate(drive(0.02), thermal, controls(), bounds(), bounds(),
        link(1.0), constant(200.0), ||false).unwrap();
    let prescribed = old.thermal_candidate(drive(0.02), UniformEnthalpyStepInput {
        internal_heat_j: 200.0*thermal.duration_s, ..thermal
    }, controls(), bounds(), link(1.0), ||false).unwrap();
    assert!((powered.coupled().thermal().state.temperature_k()-prescribed.thermal().state.temperature_k()).abs()<1e-8);
    assert!((powered.coupled().piston().velocity_m_s()-prescribed.piston().velocity_m_s()).abs()<1e-9);
    assert!(powered.step().total_energy_residual_j.abs()<=bounds().energy_tolerance_j);
    assert_eq!(powered.power().source, powered.coupled().step().thermal_midpoint);
    assert!(powered.step().exchange_trials >= powered.step().source_trials);
}

#[test]
fn implicit_midpoint_power_matches_closed_form_and_falsifies_lagged_heating() {
    let old = piston(); let curve = curve(); let thermal = body(&curve);
    let power = |_: &ResolvedThermochemicalGasState, s: EquilibriumPhaseState| -> Result<_, Infallible> {
        Ok(PowerSample { power_w: 100.0+1000.0*(s.temperature_k()-400.0), source: s })
    };
    let candidate = old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        link(0.0), power, ||false).unwrap();
    // W=dt*(100+B*(Tmid-T0)), Tmid=T0+W/(2*m*cp).
    let dt = thermal.duration_s;
    let expected = dt*100.0/(1.0-dt*1000.0/(2.0*thermal.mass_kg*500.0));
    assert!((candidate.step().applied_energy_j-expected).abs()<1e-9);
    assert!((expected-dt*100.0).abs()>0.04);
    assert!((candidate.coupled().thermal().state.temperature_k()-400.0-expected/0.5).abs()<1e-8);
    assert_eq!(candidate.coupled().piston(), &old);
}

#[test]
fn source_shutoff_and_cooling_preserve_finite_storage_and_actual_gas_motion() {
    let mut piston = piston(); let initial = piston.clone(); let curve = curve();
    let mut thermal = body(&curve); let h0 = thermal.initial.specific_enthalpy_j_kg();
    let mut supplied = 0.0;
    for i in 0..40 {
        let p = if i<20 {200.0} else {0.0};
        let proposal = piston.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
            link(5.0), constant(p), ||false).unwrap();
        supplied += proposal.step().evaluated_energy_j;
        assert!(proposal.step().source_residual_j.abs()<=bounds().source_tolerance_j);
        let (coupled, _, report) = proposal.into_parts();
        assert!(report.total_energy_residual_j.abs()<=bounds().energy_tolerance_j);
        let (next, accepted, _) = coupled.into_parts();
        piston = next; thermal.initial = accepted.state;
    }
    assert!(piston.position_m()>initial.position_m());
    assert!(piston.gas().state().state().sound_speed>initial.gas().state().state().sound_speed);
    assert!(thermal.initial.temperature_k()>400.0);
    let stored = piston.gas().mass_kg()*(piston.gas().state().specific_internal_energy_j_kg()
        -initial.gas().state().specific_internal_energy_j_kg())
        +thermal.mass_kg*(thermal.initial.specific_enthalpy_j_kg()-h0);
    let kinetic = 0.5*piston.mass_kg()*piston.velocity_m_s().powi(2);
    let exterior = 100000.0*(piston.gas().volume_m3()-initial.gas().volume_m3());
    assert!((stored+kinetic+exterior-supplied).abs()<40.0*bounds().energy_tolerance_j);
    let cooling = initial.powered_thermal_candidate(drive(0.0), body(&curve), controls(), bounds(), bounds(),
        link(0.0), constant(-100.0), ||false).unwrap();
    assert!(cooling.coupled().thermal().state.temperature_k()<400.0);
    assert!(cooling.step().evaluated_energy_j<0.0);
}

#[test]
fn source_failures_nonconvergence_and_late_cancellation_return_no_candidate() {
    let old = piston(); let saved = old.clone(); let curve = curve(); let thermal = body(&curve);
    let conductive = |_: &ResolvedThermochemicalGasState, s| -> Result<_, &str> {
        Ok(ConductanceSample { conductance_w_k: 1.0, source: s })
    };
    let failed = |_: &ResolvedThermochemicalGasState, _: EquilibriumPhaseState|
        -> Result<PowerSample<()>, &str> { Err("missing electrical material") };
    assert!(matches!(old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        conductive, failed, ||false), Err(ThermalPistonError::Interface("missing electrical material"))));
    assert!(old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        link(1.0), constant(f64::NAN), ||false).is_err());
    assert!(matches!(old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(),
        HeatExchangeSolve { heat_bracket_j: [-0.001,0.001], ..bounds() },
        link(1.0), constant(100.0), ||false), Err(ThermalPistonError::OutsideBracket {..})));
    let nonlinear = |_: &ResolvedThermochemicalGasState, s: EquilibriumPhaseState| -> Result<_, Infallible> {
        Ok(PowerSample { power_w: 100.0+(s.temperature_k()-400.0).powi(2), source: s })
    };
    assert!(matches!(old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(),
        HeatExchangeSolve { maximum_iterations:1, ..bounds() }, link(0.0), nonlinear, ||false),
        Err(ThermalPistonError::NotConverged {..})));
    let mut total = 0;
    old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        link(1.0), constant(100.0), ||{total+=1;false}).unwrap();
    let mut calls = 0;
    assert!(matches!(old.powered_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        link(1.0), constant(100.0), ||{calls+=1;calls==total}), Err(ThermalPistonError::Cancelled)));
    assert_eq!(old, saved);
}
