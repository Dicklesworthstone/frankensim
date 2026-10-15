use super::*;
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use std::convert::Infallible;

fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-1.0,1.0], source_tolerance_j:1e-9,
        energy_tolerance_j:1e-7, maximum_iterations:128 }
}
fn piston() -> GasPiston { GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap() }
fn source(g: f64) -> impl FnMut(&ResolvedThermochemicalGasState)
    -> Result<ConductanceSample<ResolvedThermochemicalGasState>,Infallible> {
    move |state| Ok(ConductanceSample { conductance_w_k:g, source:state.clone() })
}

#[test]
fn isolated_port_matches_original_motion_and_preserves_source() {
    let old=piston(); let saved=old.clone();
    let next=old.thermal_boundary_candidate(drive(0.02),controls(),exchange(),500.0,source(0.0),||false).unwrap();
    let mut original=old.clone();let report=original.step(drive(0.02),controls(),||false).unwrap();
    assert_eq!(next.piston(),&original);assert_eq!(next.step().piston,report);
    assert_eq!(next.source(),&next.step().piston.midpoint);
    assert_eq!(next.step().exchanged_heat_j,0.0);assert_eq!(next.step().exchange_trials,1);
    assert_eq!(old,saved);
}

#[test]
fn midpoint_boundary_heat_drives_motion_and_reverses_with_cooling() {
    for (temperature,sign) in [(500.0,1.0),(300.0,-1.0)] {
        let old=piston();
        let next=old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),temperature,
            |g| Ok::<_,Infallible>(ConductanceSample { conductance_w_k:10.0*g.state().thermal_conductivity,
                source:g.clone() }),||false).unwrap();
        let report=next.step();
        let expected=drive(0.0).duration_s*report.conductance_w_k
            *(temperature-report.piston.midpoint.state().temperature);
        assert!((expected-report.exchanged_heat_j).abs()<=exchange().source_tolerance_j);
        assert!(sign*report.exchanged_heat_j>0.0);
        assert!(sign*next.piston().velocity_m_s()>0.0);
        assert!(sign*(next.piston().gas().state().state().sound_speed-old.gas().state().state().sound_speed)>0.0);
        assert!(report.total_energy_residual_j.abs()<=exchange().energy_tolerance_j);
        assert!(report.material_trials>=report.piston.material_trials);
        assert_eq!(next.source(),&report.piston.midpoint);
    }
}

#[test]
fn equilibrium_does_not_bypass_source_admission_and_failures_preserve_history() {
    let old=piston();let saved=old.clone();
    let failure=|_:&ResolvedThermochemicalGasState| -> Result<ConductanceSample<()>,&str> { Err("missing source") };
    assert!(matches!(old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),400.0,failure,||false),
        Err(ThermalPistonError::Interface("missing source"))));
    assert!(old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),500.0,source(-1.0),||false).is_err());
    assert!(matches!(old.thermal_boundary_candidate(drive(0.0),controls(),
        HeatExchangeSolve { heat_bracket_j:[-1e-6,1e-6],..exchange() },500.0,source(1.0),||false),
        Err(ThermalPistonError::OutsideBracket { .. })));
    let mut count=0;
    old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),500.0,source(1.0),||{count+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),500.0,source(1.0),||{
        actual+=1;actual==count
    }),Err(ThermalPistonError::Cancelled)));
    for t in [0.0,f64::NAN,f64::INFINITY] {
        assert!(old.thermal_boundary_candidate(drive(0.0),controls(),exchange(),t,source(1.0),||false).is_err());
    }
    assert_eq!(old,saved);
}
