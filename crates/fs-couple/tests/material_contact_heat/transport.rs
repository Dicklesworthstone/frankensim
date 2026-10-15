//! Actual shared thermal owner and complete coupled transaction checks.
use super::*;
use fs_couple::render::plate::impact::ImpactError;
use fs_conduction::lumped::{LumpedThermalEnvironment,LumpedThermalTransport};
use fs_material::phase::{UniformEnthalpyStep,UniformEnthalpyStepInput};
use fs_exec::{Budget,Cx,ExecMode,StreamKey};

fn context<R>(f:impl FnOnce(&CancelGate,&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:0,kernel_id:0x434854,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        f(&gate,&cx)
    })
}
fn environment()->LumpedThermalEnvironment {
    LumpedThermalEnvironment {temperature_k:300.0,radiation_temperature_k:300.0,
        convection_w_per_m2_k:1.0,transport:LumpedThermalTransport::try_declared(100.0,0.0).unwrap(),
        enthalpy_tolerance_j_kg:1e-12,maximum_thermal_residual_j:1e-14}
}
fn adiabatic(input:UniformEnthalpyStepInput<'_>)->UniformEnthalpyStep<()> {
    UniformEnthalpyStep {state:input.curve.state_at_specific_enthalpy(
        input.initial.specific_enthalpy_j_kg()+input.internal_heat_j/input.mass_kg).unwrap(),
        external_heat_j:0.0,energy_residual_tolerance_j:0.0,report:()}
}

#[test]
fn shared_conduction_consumes_actual_contact_heat_and_continues_cooling_after_separation() {
    let c=curve(0.02,1.0);let slower=curve(0.04,1.0);
    context(|gate,cx| {
        let mut a=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
            stores(&c),vec![route("touch")],1e-10).unwrap();
        let mut b=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
            stores(&slower),vec![route("touch")],1e-10).unwrap();
        let env=environment();let initial=a.mechanics().stored_energy_j();let mut peak=300.0_f64;
        let mut boundary=0.0;let mut calls=0;
        for _ in 0..2000 {
            let fa=a.step_with_transport(&[0.0;2],gate,|_,input| {
                calls+=1;env.advance(cx,input).map_err(|e|ImpactError::Owner(e.to_string()))
            }).unwrap();
            b.step_with_transport(&[0.0;2],gate,|_,input|
                env.advance(cx,input).map_err(|e|ImpactError::Owner(e.to_string()))).unwrap();
            assert_eq!(fa.transfers.len(),2);assert_eq!(a.mechanics().state(),b.mechanics().state());
            assert!(fa.impact.combined_energy_residual_j.abs()<1e-9);
            boundary+=fa.boundary_heat_j;peak=peak.max(a.reservoirs()[0].state().temperature_k());
        }
        assert_eq!(calls,4000);assert!(boundary<0.0);
        assert!(peak>a.reservoirs()[0].state().temperature_k());
        assert!(a.reservoirs()[0].state().temperature_k()>b.reservoirs()[0].state().temperature_k());
        near(a.mechanics().stored_energy_j()+a.reservoirs().iter().map(|s|s.storage_change_j()).sum::<f64>()-boundary,initial,1e-9);
        near(a.reservoirs().iter().map(|s|s.boundary_heat_j()).sum::<f64>(),boundary,1e-15);
        // Switching off external exchange does not reset heat to its initial datum.
        let state=a.reservoirs()[0].state();let external=a.reservoirs()[0].boundary_heat_j();
        let f=a.step(&[0.0;2],gate).unwrap();assert_eq!(f.deposited_heat_j,0.0);
        assert_eq!(a.reservoirs()[0].state(),state);assert_eq!(a.reservoirs()[0].boundary_heat_j(),external);
    });
}

#[test]
fn adiabatic_carrier_matches_default_heat_and_keeps_concrete_reports() {
    let c=curve(1.0,1.0);let gate=CancelGate::new_clock_free();
    let make=||HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),stores(&c),vec![route("touch")],1e-12).unwrap();
    let mut a=make();let mut b=make();
    for _ in 0..1400 {
        a.step(&[0.0;2],&gate).unwrap();
        let f=b.step_with_transport(&[0.0;2],&gate,|key,input| {
            let p=adiabatic(input);Ok(UniformEnthalpyStep {state:p.state,external_heat_j:0.0,
                energy_residual_tolerance_j:0.0,report:(key.to_owned(),input.internal_heat_j,input.duration_s)})
        }).unwrap();
        assert_eq!(a.mechanics().state(),b.mechanics().state());
        assert_eq!(f.transfers[0].report.0,"left");assert_eq!(f.transfers[1].report.0,"right");
        assert_eq!(f.transfers[0].report.2,config().dt_s);
        for (x,y) in a.reservoirs().iter().zip(b.reservoirs()) {
            near(x.state().specific_enthalpy_j_kg(),y.state().specific_enthalpy_j_kg(),1e-12);
        }
    }
}

#[test]
fn thermal_source_failure_foreign_state_and_late_cancellation_are_atomic() {
    let c=curve(1.0,1.0);let foreign=curve(2.0,1.0);let gate=CancelGate::new_clock_free();
    let mut s=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),stores(&c),vec![route("touch")],1e-12).unwrap();
    let motion=s.mechanics().state().to_vec();let stores=s.reservoirs().to_vec();
    assert!(s.step_with_transport(&[0.0;2],&gate,|key,input| {
        if key=="right" {return Err(ImpactError::Invalid("injected transport refusal"));}
        Ok(adiabatic(input))
    }).is_err());
    assert!(s.step_with_transport(&[0.0;2],&gate,|_,input| {
        let mut p=adiabatic(input);p.state=foreign.state_at_specific_enthalpy(0.0).unwrap();Ok(p)
    }).is_err());
    assert!(s.step_with_transport(&[0.0;2],&gate,|_,input| {
        let mut p=adiabatic(input);p.state=c.state_at_specific_enthalpy(1e-4).unwrap();Ok(p)
    }).is_err());
    assert!(s.step_with_transport(&[0.0;2],&gate,|_,input| {
        let mut p=adiabatic(input);p.energy_residual_tolerance_j=1.0;Ok(p)
    }).is_err());
    let late=CancelGate::new_clock_free();let mut calls=0;
    assert!(s.step_with_transport(&[0.0;2],&late,|_,input| {
        calls+=1;if calls==2 {late.request();}Ok(adiabatic(input))
    }).is_err());
    assert_eq!(calls,2);assert_eq!(s.mechanics().state(),motion);assert_eq!(s.reservoirs(),stores);
    assert_eq!(s.mechanics().samples(),0);
    // Failed candidate scratch cannot poison a retry of the unchanged state.
    s.step(&[0.0;2],&gate).unwrap();assert_eq!(s.mechanics().samples(),1);
}

#[path="coupled.rs"]
mod coupled;
#[path="refresh.rs"]
mod refresh;
