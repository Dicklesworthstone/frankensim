//! Real gas/momentum owner plus source-bound enthalpy/phase charts.
mod real_gas;

use super::*;
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use std::convert::Infallible;

fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn curve(cp: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(fs_blake3::hash_bytes(format!("synthetic cp={cp}").as_bytes()),
        SolidLiquidPhase::Solid, [200.0, 800.0].into_iter().map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: cp*(t-200.0), temperature_k:t,
            liquid_mass_fraction:0.0, bulk_density_kg_m3:4000.0,
        }).collect()).unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial:curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg:0.001, volume_m3:0.001/4000.0, surface_area_m2:0.02,
        internal_heat_j:0.0, duration_s:drive(0.0).duration_s }
}
fn link(g: f64) -> impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
    -> Result<ConductanceSample<EquilibriumPhaseState>, Infallible>
{
    move |_, state| Ok(ConductanceSample { conductance_w_k:g, source:state })
}
fn piston() -> GasPiston { GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap() }

#[test]
fn equilibrium_and_isolation_use_the_same_existing_step() {
    let old=piston(); let curve=curve(500.0);
    let equilibrium=body(&curve,100000.0);
    let candidate=old.thermal_candidate(drive(0.0),equilibrium,controls(),exchange(),link(1.0),||false).unwrap();
    assert_eq!(candidate.piston(),&old);
    assert_eq!(candidate.thermal().state,equilibrium.initial);
    assert_eq!(candidate.step().exchanged_heat_j,0.0);
    assert_eq!(candidate.step().exchange_trials,1);
    let hot=body(&curve,150000.0);
    let candidate=old.thermal_candidate(drive(0.02),hot,controls(),exchange(),link(0.0),||false).unwrap();
    let mut independent=old.clone();
    let independent_step=independent.step(drive(0.02),controls(),||false).unwrap();
    assert_eq!(candidate.piston(),&independent);
    assert_eq!(candidate.step().piston,independent_step);
    assert_eq!(candidate.thermal().state,hot.initial);
    assert_eq!(candidate.thermal().external_heat_j,0.0);
}

fn run(cp: f64, g: f64) -> (GasPiston, f64, f64) {
    let curve=curve(cp); let mut thermal=body(&curve,cp*300.0);
    let mut piston=piston(); let original=piston.clone();
    let h0=thermal.initial.specific_enthalpy_j_kg();
    let mut total_q=0.0;
    for _ in 0..80 {
        let candidate=piston.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(g),||false).unwrap();
        let (next, accepted, report)=candidate.into_parts();
        assert!(report.source_residual_j.abs()<=exchange().source_tolerance_j);
        assert!(report.total_energy_residual_j.abs()<=exchange().energy_tolerance_j);
        assert_eq!(accepted.report,report.thermal_midpoint);
        assert_eq!(accepted.external_heat_j,-report.exchanged_heat_j);
        thermal.initial=accepted.state; piston=next; total_q+=report.exchanged_heat_j;
        let gas_storage=piston.gas().mass_kg()*(piston.gas().state().specific_internal_energy_j_kg()
            -original.gas().state().specific_internal_energy_j_kg());
        let kinetic=0.5*piston.mass_kg()*piston.velocity_m_s().powi(2);
        let wall_storage=thermal.mass_kg*(thermal.initial.specific_enthalpy_j_kg()-h0);
        let exterior=100000.0*(piston.gas().volume_m3()-original.gas().volume_m3());
        assert!((gas_storage+kinetic+wall_storage+exterior).abs()<80.0*exchange().energy_tolerance_j);
        assert_eq!(thermal.initial.phase_curve_identity(),thermal.curve.identity());
    }
    assert!(total_q>0.0);
    assert!(thermal.initial.temperature_k()<500.0);
    assert!(piston.gas().state().state().temperature>400.0);
    assert_eq!(piston.gas().mass_kg(),original.gas().mass_kg());
    (piston,thermal.initial.temperature_k(),total_q)
}
#[test]
fn finite_heat_capacity_and_transport_data_change_the_actual_motion() {
    let (a,ta,qa)=run(200.0,1.0); let (b,tb,qb)=run(2000.0,1.0);
    assert!(tb>ta && qb>qa);
    assert!((a.position_m()-b.position_m()).abs()>1e-8);
    assert!((a.gas().state().state().sound_speed-b.gas().state().state().sound_speed).abs()>1e-3);
    let (c,_,qc)=run(200.0,2.0);
    assert!(qc>qa);
    assert!((a.position_m()-c.position_m()).abs()>1e-8);
}

#[test]
fn latent_heat_changes_phase_without_inventing_a_temperature_rise() {
    let knot=|h,t,x| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,temperature_k:t,
        liquid_mass_fraction:x,bulk_density_kg_m3:4000.0 };
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(fs_blake3::hash_bytes(b"synthetic latent store"),vec![
        knot(0.0,200.0,0.0),knot(100000.0,350.0,0.0),
        knot(160000.0,350.0,1.0),knot(400000.0,800.0,1.0),
    ]).unwrap();
    let thermal=body(&curve,130000.0);
    let old=piston();
    let candidate=old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(1.0),||false).unwrap();
    assert!(candidate.step().exchanged_heat_j<0.0);
    assert_eq!(candidate.thermal().state.temperature_k(),350.0);
    assert!(candidate.thermal().state.liquid_mass_fraction()>thermal.initial.liquid_mass_fraction());
    assert_eq!(candidate.step().thermal_midpoint.temperature_k(),350.0);
    assert!(candidate.piston().gas().state().state().temperature<400.0);
    assert!(candidate.piston().velocity_m_s()<0.0);
    assert!(candidate.step().total_energy_residual_j.abs()<exchange().energy_tolerance_j);
}

#[test]
fn both_prescribed_sources_and_nonlinear_enthalpy_share_one_clock_and_balance() {
    let knot=|h,t| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,temperature_k:t,
        liquid_mass_fraction:0.0,bulk_density_kg_m3:4000.0 };
    let curve=EquilibriumEnthalpyPhaseCurve::try_single_phase(fs_blake3::hash_bytes(b"synthetic variable cp"),
        SolidLiquidPhase::Solid,vec![knot(0.0,200.0),knot(100000.0,400.0),knot(500000.0,800.0)]).unwrap();
    let thermal=UniformEnthalpyStepInput { internal_heat_j:25.0,..body(&curve,90000.0) };
    let candidate=piston().thermal_candidate(drive(0.02),thermal,controls(),exchange(),link(1.0),||false).unwrap();
    let report=candidate.step();
    let averaged=0.5*(thermal.initial.temperature_k()+candidate.thermal().state.temperature_k());
    assert!((report.thermal_midpoint.temperature_k()-averaged).abs()>1.0);
    assert!((report.thermal_energy_change_j+report.exchanged_heat_j-25.0).abs()<1e-7);
    assert!((report.piston.gas.supplied_heat_j-0.02-report.exchanged_heat_j).abs()<1e-12);
    assert!(report.total_energy_residual_j.abs()<exchange().energy_tolerance_j);
}

#[test]
fn source_domain_and_chart_failures_publish_neither_state() {
    let old=piston(); let saved=old.clone(); let curve=curve(500.0);
    let thermal=body(&curve,150000.0); let before=thermal.initial;
    let failure=|_:&ResolvedThermochemicalGasState,_| -> Result<ConductanceSample<()>, &str> { Err("missing interface data") };
    assert!(matches!(old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),failure,||false),
        Err(ThermalPistonError::Interface("missing interface data"))));
    assert!(old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(-1.0),||false).is_err());
    let other=EquilibriumEnthalpyPhaseCurve::try_single_phase(fs_blake3::hash_bytes(b"other source"),
        SolidLiquidPhase::Solid,curve.knots().to_vec()).unwrap();
    assert!(old.thermal_candidate(drive(0.0),UniformEnthalpyStepInput { curve:&other,..thermal },
        controls(),exchange(),link(1.0),||false).is_err());
    assert!(old.thermal_candidate(drive(0.0),UniformEnthalpyStepInput { duration_s:0.1,..thermal },
        controls(),exchange(),link(1.0),||false).is_err());
    assert!(matches!(old.thermal_candidate(drive(0.0),thermal,controls(),
        HeatExchangeSolve { heat_bracket_j:[-1e-6,1e-6],..exchange() },link(1.0),||false),
        Err(ThermalPistonError::OutsideBracket { .. })));
    assert!(old.thermal_candidate(drive(0.0),UniformEnthalpyStepInput { internal_heat_j:1e9,..thermal },
        controls(),exchange(),link(1.0),||false).is_err());
    assert_eq!(old,saved); assert_eq!(thermal.initial,before);
}

#[test]
fn exhausted_work_and_late_cancellation_do_not_return_a_partial_candidate() {
    let old=piston(); let curve=curve(500.0); let thermal=body(&curve,150000.0);
    assert!(matches!(old.thermal_candidate(drive(0.0),thermal,controls(),
        HeatExchangeSolve { maximum_iterations:1,..exchange() },link(1.0),||false),
        Err(ThermalPistonError::NotConverged { .. })));
    let mut total=0;
    old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(1.0),||{total+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(1.0),
        ||{actual+=1;actual==total}),Err(ThermalPistonError::Cancelled)));
    assert!(matches!(old.thermal_candidate(drive(0.0),thermal,controls(),exchange(),link(1.0),||true),
        Err(ThermalPistonError::Cancelled)));
}
