//! Both collision-heated finite bodies advance through the ORIGINAL heat solver.
use super::*;
use fs_conduction::interface::{ResistanceUncertainty,SeriesThermalResistance,ThermalResistanceTerm};
use fs_conduction::lumped::{LumpedContactStep,LumpedContactStepConfig,LumpedEnthalpyBody,solve_lumped_contact_step};
use fs_couple::render::plate::impact::material::heating::transport::{
    CoupledThermalStep,ReservoirThermalInput,
};

fn resistance()->SeriesThermalResistance {
    SeriesThermalResistance::new(vec![ThermalResistanceTerm::declared("finite link",1000.0,
        ResistanceUncertainty::Unstated,"synthetic contact, no measured authority").unwrap()]).unwrap()
}
fn warm_stores(c:&EquilibriumEnthalpyPhaseCurve)->Vec<CaloricReservoir<'_>> {
    [("left",0.4),("right",0.1)].into_iter().map(|(key,h)|CaloricReservoir::new(key.into(),c,
        c.state_at_specific_enthalpy(h).unwrap(),0.001,1e-6,0.0006).unwrap()).collect()
}
fn body<'a>(name:&str,i:UniformEnthalpyStepInput<'a>)->LumpedEnthalpyBody<'a> {
    LumpedEnthalpyBody::try_new(name,i.mass_kg,i.surface_area_m2,0.0,0.0,
        i.volume_m3/i.surface_area_m2,100.0,i.curve).unwrap()
}
fn exchange(cx:&Cx<'_>,inputs:&[ReservoirThermalInput<'_>])->CoupledThermalStep<LumpedContactStep> {
    assert_eq!(inputs.len(),2);assert_eq!(inputs[0].key,"left");assert_eq!(inputs[1].key,"right");
    let a=inputs[0].thermal;let b=inputs[1].thermal;
    assert_eq!(a.duration_s,b.duration_s);
    let result=solve_lumped_contact_step(cx,&body("left",a),&body("right",b),&resistance(),
        LumpedContactStepConfig {
            specific_enthalpy_a_j_kg:a.initial.specific_enthalpy_j_kg(),
            specific_enthalpy_b_j_kg:b.initial.specific_enthalpy_j_kg(),
            supplied_heat_a_j:a.internal_heat_j,supplied_heat_b_j:b.internal_heat_j,
            duration_s:a.duration_s,energy_tolerance_j:1e-13,
        }).unwrap();
    // Independent constant-Cp identity. Sources participate in the SAME equation.
    let cap=0.001*0.02;
    let delta=a.initial.temperature_k()-b.initial.temperature_k()
        +(a.internal_heat_j-b.internal_heat_j)/cap;
    let q=a.duration_s*delta/1000.0/(1.0+2.0*a.duration_s/(1000.0*cap));
    near(result.heat_a_to_b_j,q,1.1e-13);
    let transfers=[(result.state_a,-result.heat_a_to_b_j),(result.state_b,result.heat_a_to_b_j)]
        .into_iter().map(|(state,external_heat_j)|UniformEnthalpyStep {state,external_heat_j,
            energy_residual_tolerance_j:1e-13,report:result}).collect();
    CoupledThermalStep {transfers,boundary_heat_j:0.0}
}

#[test]
fn finite_partners_exchange_heat_on_one_collision_clock_without_double_debit() {
    let c=curve(0.02,1.0);
    context(|gate,cx| {
        let mut s=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
            warm_stores(&c),vec![route("touch")],1e-10).unwrap();
        let mut plain=collision(Some(0.4),0.0,false).prepare_analytic().unwrap();
        let initial=plain.stored_energy_j();let mut calls=0;let mut peak=0.0_f64;
        let mut heat_after_separation=0.0;
        for tick in 1..=1800 {
            plain.step(&[0.0;2],gate).unwrap();
            let expected=plain.material_contact_observation(0).unwrap();
            let old=s.reservoirs().to_vec();
            let frame=s.step_with_coupled_transport(&[0.0;2],gate,|trial| {
                calls+=1;assert_eq!(trial.sample(),tick);assert_eq!(trial.duration_s(),config().dt_s);
                let observation=trial.contact_observation("touch")?;
                assert_eq!(observation.normal_force_n,expected.normal_force_n);
                peak=peak.max(observation.normal_force_n);
                assert!(trial.contact_observation("missing").is_err());
                for (input,store) in trial.reservoirs().iter().zip(&old) {
                    assert_eq!(input.thermal.initial,store.state());
                }
                Ok(exchange(cx,trial.reservoirs()))
            }).unwrap();
            assert_eq!(s.mechanics().state(),plain.state());assert_eq!(frame.boundary_heat_j,0.0);
            assert_eq!(frame.transfers[0].external_heat_j,-frame.transfers[1].external_heat_j);
            if tick>1500 && expected.normal_force_n==0.0 {heat_after_separation+=frame.transfers[1].external_heat_j;}
            assert!(frame.impact.combined_energy_residual_j.abs()<1e-9);
        }
        assert_eq!(calls,1800);assert!(peak>0.0 && heat_after_separation>0.0);
        assert!(s.reservoirs()[0].state().temperature_k()<320.0);
        assert!(s.reservoirs()[1].state().temperature_k()>305.0);
        near(s.reservoirs().iter().map(|s|s.boundary_heat_j()).sum(),0.0,1e-16);
        near(s.mechanics().stored_energy_j()+s.reservoirs().iter().map(|s|s.storage_change_j()).sum::<f64>(),initial,1e-9);
        // Disabling the physical link resumes from the accepted shared endpoint.
        let before=s.reservoirs().to_vec();s.step(&[0.0;2],gate).unwrap();
        assert_eq!(s.reservoirs(),before);
    });
}

#[test]
fn joint_and_independent_adiabatic_paths_share_exact_acceptance_and_retry() {
    let c=curve(0.02,1.0);let gate=CancelGate::new_clock_free();
    let make=||HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
        warm_stores(&c),vec![route("touch")],1e-10).unwrap();
    let mut a=make();let mut b=make();
    for _ in 0..900 {
        a.step_with_transport(&[0.0;2],&gate,|_,i|Ok(adiabatic(i))).unwrap();
        b.step_with_coupled_transport(&[0.0;2],&gate,|trial|Ok(CoupledThermalStep {
            transfers:trial.reservoirs().iter().map(|i|adiabatic(i.thermal)).collect(),boundary_heat_j:0.0,
        })).unwrap();
    }
    assert_eq!(a.reservoirs(),b.reservoirs());assert_eq!(a.mechanics().state(),b.mechanics().state());
}

#[test]
fn missing_partner_debit_wrong_count_foreign_state_and_late_cancellation_roll_back_every_owner() {
    let c=curve(0.02,1.0);let foreign=curve(0.04,1.0);let gate=CancelGate::new_clock_free();
    let mut s=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
        warm_stores(&c),vec![route("touch")],1e-10).unwrap();
    let saved=s.reservoirs().to_vec();let motion=s.mechanics().state().to_vec();
    for failure in 0..6 {
        let result=s.step_with_coupled_transport(&[0.0;2],&gate,|trial| {
            let mut p=CoupledThermalStep {transfers:trial.reservoirs().iter()
                .map(|i|adiabatic(i.thermal)).collect(),boundary_heat_j:0.0};
            match failure {
                0=>{p.transfers.pop();},
                1=>{p.transfers.push(adiabatic(trial.reservoirs()[0].thermal));},
                2=>{ // Locally exact heating, but no donor or outside source exists.
                    p.transfers[1].external_heat_j=1e-6;
                    p.transfers[1].state=c.state_at_specific_enthalpy(0.101).unwrap();
                },
                3=>{p.transfers[1].state=foreign.state_at_specific_enthalpy(0.1).unwrap();},
                4=>{p.boundary_heat_j=f64::NAN;},
                _=>return Err(ImpactError::Invalid("joint thermal law refused")),
            }
            Ok(p)
        });
        assert!(result.is_err());assert_eq!(s.reservoirs(),saved);
        assert_eq!(s.mechanics().state(),motion);assert_eq!(s.mechanics().samples(),0);
    }
    let late=CancelGate::new_clock_free();let mut calls=0;
    assert!(s.step_with_coupled_transport(&[0.0;2],&late,|trial| {
        calls+=1;let transfers=trial.reservoirs().iter().map(|i|adiabatic(i.thermal)).collect();
        late.request();Ok(CoupledThermalStep {transfers,boundary_heat_j:0.0})
    }).is_err());
    assert_eq!(calls,1);assert_eq!(s.reservoirs(),saved);assert_eq!(s.mechanics().state(),motion);
    s.step(&[0.0;2],&gate).unwrap();assert_eq!(s.mechanics().samples(),1);
}
