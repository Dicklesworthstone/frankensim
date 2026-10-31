//! G1/G3/G4: source substitutions, original section equations and bounded refusal.
use super::*;
use fs_material::{MaterialError, fluid_transport::FluidTransportState};
use fs_qty::{Area,Length};
fn law(a:f64,b:f64)->MaterialFlowResistance {
    MaterialFlowResistance::try_new(a,b,Area::new(1e-4),Length::new(0.01),(0.0,1e9)).unwrap()
}
fn fluid(rho:f64,mu:f64)->FluidTransportState {
    FluidTransportState::try_new(rho,4000.0,0.5,mu).unwrap()
}
fn budget()->SeriesPressureBudget {SeriesPressureBudget {maximum_sections:8,
    absolute_pressure_tolerance_pa:1e-10,relative_pressure_tolerance:1e-12}}
fn sections()->Vec<PressureSection<FluidTransportState>> {vec![
    PressureSection {resistance:law(1e10,2e8),medium:fluid(1000.0,0.002)},
    PressureSection {resistance:law(3e10,4e8),medium:fluid(500.0,0.001)},
]}
fn near(a:f64,b:f64) {assert!((a-b).abs()<=1e-12*(1.0+a.abs()+b.abs()),"{a:e} != {b:e}");}

#[test]
fn original_series_equations_keep_one_mass_and_distinct_local_volume_and_power() {
    let gate=CancelGate::new_clock_free();let solve=SeriesPressureFlow::solve(sections(),1000.0,budget(),&gate).unwrap();
    // Independently solve 80000*m + 1000000*m^2 = 1000.
    let expected=(-80000.0_f64+(80000.0_f64.powi(2)+4.0e9).sqrt())/(2.0e6);
    near(solve.mass_flow_kg_s(),expected);
    let mut pressure=0.0;let mut power=0.0;
    for (section,state) in solve.sections().iter().zip(solve.evaluations()) {
        assert_eq!(&state.transport,&section.medium);
        assert_eq!(state.mass_flow_kg_per_s,solve.mass_flow_kg_s());
        near(state.volume_flow_m3_s,expected/section.medium.density_kg_m3());
        pressure+=state.pressure_drop_pa;power+=state.dissipated_power_w;
    }
    near(pressure,1000.0);near(power,solve.dissipated_power_w());
    assert!((power-1000.0*solve.evaluations()[0].volume_flow_m3_s).abs()>1e-3);
    near(solve.dpressure_dmass_flow(),80000.0+2.0e6*expected);
    let reverse=SeriesPressureFlow::solve(sections(),-1000.0,budget(),&gate).unwrap();
    assert_eq!(reverse.mass_flow_kg_s(),-solve.mass_flow_kg_s());
    near(reverse.dissipated_power_w(),solve.dissipated_power_w());
    for (a,b) in solve.evaluations().iter().zip(reverse.evaluations()) {
        near(a.pressure_drop_pa,-b.pressure_drop_pa);
        near(a.dissipated_power_w,b.dissipated_power_w);
    }
}

#[test]
fn data_only_viscosity_changes_flow_without_touching_geometry_and_split_sections_agree() {
    let gate=CancelGate::new_clock_free();
    let base=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(4e10,4e8),medium:fluid(1000.0,0.001)}],
        1000.0,budget(),&gate).unwrap();
    let split=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(1e10,1e8),medium:fluid(1000.0,0.001)},
        PressureSection {resistance:law(3e10,3e8),medium:fluid(1000.0,0.001)}],1000.0,budget(),&gate).unwrap();
    near(base.mass_flow_kg_s(),split.mass_flow_kg_s());near(base.dissipated_power_w(),split.dissipated_power_w());
    let viscous=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(4e10,4e8),medium:fluid(1000.0,0.004)}],
        1000.0,budget(),&gate).unwrap();
    assert!(viscous.mass_flow_kg_s()<base.mass_flow_kg_s());
    assert!(viscous.dissipated_power_w()<base.dissipated_power_w());
    assert_eq!(base.sections()[0].resistance,viscous.sections()[0].resistance);
    for (a,b) in [(1e10,0.0),(0.0,1e8)] {
        let rest=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(a,b),medium:fluid(1000.0,0.001)}],
            0.0,budget(),&gate).unwrap();
        assert_eq!(rest.mass_flow_kg_s(),0.0);assert_eq!(rest.dissipated_power_w(),0.0);
    }
}

#[test]
fn one_query_per_original_source_and_thermal_properties_do_not_invent_pressure() {
    struct Source {identity:&'static str,calls:std::cell::Cell<usize>,state:FluidTransportState}
    impl FluidTransport for Source {
        fn transport_state(&self)->core::result::Result<FluidTransportState,MaterialError> {
            self.calls.set(self.calls.get()+1);Ok(self.state)
        }
    }
    let gate=CancelGate::new_clock_free();
    let source=Source {identity:"retained source",calls:std::cell::Cell::new(0),state:fluid(1000.0,0.002)};
    let result=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(1e10,2e8),medium:source}],
        100.0,budget(),&gate).unwrap();
    assert_eq!(result.sections()[0].medium.identity,"retained source");
    assert_eq!(result.sections()[0].medium.calls.get(),1);
    let changed=FluidTransportState::try_new(1000.0,2000.0,2.0,0.002).unwrap();
    let other=SeriesPressureFlow::solve(vec![PressureSection {resistance:law(1e10,2e8),medium:changed}],
        100.0,budget(),&gate).unwrap();
    assert_eq!(result.mass_flow_kg_s(),other.mass_flow_kg_s());
    assert_eq!(result.dissipated_power_w(),other.dissipated_power_w());
}

#[test]
fn inverse_scaling_preserves_extreme_representable_single_term_solutions() {
    for (p,a,b,expected) in [(1.0,1e200,0.0,1e-200),(1.0,0.0,1e-200,1e100),
        (1e200,1e100,0.0,1e100),(1e-200,0.0,1e200,1e-200)] {
        let x=invert(p,a,b).unwrap();assert!((x/expected-1.0).abs()<1e-14);
    }
    for (a,b) in [(1e-200,1e200),(1e200,1e-200)] {
        let x=invert(1.0,a,b).unwrap();near(a*x+(b*x)*x,1.0);
    }
    assert!(invert(1.0,0.0,0.0).is_err());
}

#[test]
fn unsupported_regime_limits_and_cancellation_refuse_without_a_fallback() {
    let gate=CancelGate::new_clock_free();
    let narrow=MaterialFlowResistance::try_new(1e10,0.0,Area::new(1e-4),Length::new(0.01),(1.0,10.0)).unwrap();
    for drive in [0.0,1000.0] {
        assert!(matches!(SeriesPressureFlow::solve(vec![PressureSection {resistance:narrow,medium:fluid(1000.0,0.001)}],
            drive,budget(),&gate),Err(SeriesPressureError::Section {index:0,source:ResistanceError::Reynolds {..}})));
    }
    assert!(SeriesPressureFlow::solve(sections(),f64::NAN,budget(),&gate).is_err());
    assert!(SeriesPressureFlow::solve(sections(),1000.0,SeriesPressureBudget {maximum_sections:1,..budget()},&gate).is_err());
    assert!(SeriesPressureFlow::<FluidTransportState>::solve(vec![],1.0,budget(),&gate).is_err());
    gate.request();
    assert!(matches!(SeriesPressureFlow::solve(sections(),1000.0,budget(),&gate),Err(SeriesPressureError::Cancelled)));
}
