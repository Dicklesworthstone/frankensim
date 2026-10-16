use super::*;
use fs_lbm::conjugate::enthalpy::transport::interfaces::implicit::adaptive::*;

fn controls(error:f64) -> AdaptiveJouleConfig {
    AdaptiveJouleConfig {step:config(),maximum_step_s:1.0,minimum_step_s:1e-7,
        max_step_attempts:10_000,cell_error_budget_j:error}
}
fn guess(s:&DepositionSimulation<'_>,_:f64)->Result<Vec<f64>,JouleEnthalpyError> {
    Ok(vec![0.0;s.cells().len()])
}

#[test]
fn g1_error_control_refines_actual_material_current_and_only_debits_the_fine_trajectory() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));let gate=CancelGate::new();
    let request=JouleEnthalpyInput {max_current_solves:100_000,..input(&domain,&drive,&ids,&[])};
    let run=|error|advance_joule_enthalpy_adaptive(&old,1.0,request,&[],&gate,controls(error),guess,
        |_,cell,_|EnthalpyTransport::try_new((),Some(1.0+0.25*cell.specific_h_j_kg),0.0),
        no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
    let (coarse,a)=run(0.01);let (fine,b)=run(0.0025);
    let exact=(0.25_f64.exp()-1.0)/0.25;
    assert!((fine.cells()[&ids[1]].specific_h_j_kg-exact).abs()
        <(coarse.cells()[&ids[1]].specific_h_j_kg-exact).abs());
    assert!(a.rejected_intervals>0 && b.rejected_intervals>0);
    assert!(b.accepted_times_s.len()>a.accepted_times_s.len());
    // Independent scalar backward-Euler recurrence on the ACTUAL accepted grid.
    let (mut h,mut t)=(0.0,0.0);
    for &end in &b.accepted_times_s {
        let dt=end-t;h=(h+dt)/(1.0-0.25*dt);t=end;
    }
    near(fine.cells()[&ids[1]].specific_h_j_kg,h);
    near(b.external_heat_j,h);near(b.electrode_energy_j,h);
    near(b.joule_energy_j,h);near(b.deposited_joule_j,h);
    assert_eq!(fine.time_s(),1.0);
    assert_eq!(b.last_step.material.time_s,1.0);
    assert_eq!(b.last_step.interfaces.time_s,1.0);
    assert!(b.estimated_cell_errors_j.iter().all(|(_,e)|*e<=0.0025*(1.0+1e-10)));
    assert!(b.step_attempts>b.accepted_times_s.len());
    assert!(b.current_solves>b.last_step.current_solves);
    let (repeat,c)=run(0.0025);
    assert_eq!(repeat.cells(),fine.cells());
    assert_eq!(c.accepted_times_s,b.accepted_times_s);
    assert_eq!(c.current_solves,b.current_solves);
    assert_eq!(old.time_s(),0.0);assert_eq!(old.energy_j(),0.0);
}

#[test]
fn g3_latent_plateau_cannot_hide_temporal_error_and_source_off_retains_matter() {
    let mut fractions=Vec::new();
    for l in [100.0,200.0] {
        let law=latent(l);let (old,domain,ids)=initial(&law,&[0.0,40.0,0.0],None,0.0);
        let drive=setup(3,ElectrodeDrive::Current(4.0));let gate=CancelGate::new();
        let request=JouleEnthalpyInput {max_current_solves:100_000,..input(&domain,&drive,&ids,&[])};
        let bulk=|_:usize,cell:&DepositedCell,_:f64| {
            let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
            let rho=1.0+phase.liquid_mass_fraction();EnthalpyTransport::try_new(phase,Some(rho),0.0)
        };
        let (next,r)=advance_joule_enthalpy_adaptive(&old,0.5,request,&[],&gate,controls(0.01),guess,
            bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
        let h=next.cells()[&ids[1]].specific_h_j_kg;
        let exact=40.0+(l+20.0)*(8.0/l).exp_m1();
        assert!((h-exact).abs()<0.012);
        assert_eq!(next.temperature_k(ids[1]).unwrap(),310.0);
        assert!(r.rejected_intervals>0);
        assert!(r.estimated_cell_errors_j.iter().any(|(id,e)|*id==ids[1] && *e>0.0));
        near(h-40.0,r.deposited_joule_j);
        fractions.push(r.last_step.material.cells[1].response.material().liquid_mass_fraction());
        let off=setup(3,ElectrodeDrive::Voltage(0.0));
        let (after,r)=advance_joule_enthalpy_adaptive(&next,1.0,input(&domain,&off,&ids,&[]),&[],&gate,
            controls(0.01),guess,bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
        assert_eq!(after.cells(),next.cells());assert_eq!(r.deposited_joule_j,0.0);
        for (id,c) in old.cells() {assert_eq!(c.mass_kg,after.cells()[id].mass_kg);}
    }
    assert!(fractions[0]>fractions[1]);
}

#[test]
fn g4_failed_and_discarded_nonlinear_work_never_renews_the_current_budget() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));let gate=CancelGate::new();
    let bulk=|_:usize,c:&DepositedCell,_:f64|EnthalpyTransport::try_new((),Some(1.0+0.25*c.specific_h_j_kg),0.0);
    let request=input(&domain,&drive,&ids,&[]);
    let failure=advance_joule_enthalpy_adaptive(&old,1.0,JouleEnthalpyInput {max_current_solves:1,..request},
        &[],&gate,controls(0.001),guess,bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).err().expect("interval must refuse");
    assert_eq!(failure.current_solves,1);assert_eq!(failure.step_attempts,1);
    assert!(matches!(failure.cause,AdaptiveJouleError::Step(ImplicitJouleError::Physical(JouleEnthalpyError::CurrentBudget))));
    let failure=advance_joule_enthalpy_adaptive(&old,1.0,request,&[],&gate,
        AdaptiveJouleConfig {max_step_attempts:3,..controls(1e-12)},guess,
        bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).err().expect("interval must refuse");
    assert_eq!(failure.step_attempts,3);assert!(failure.current_solves>3);
    assert!(matches!(failure.cause,AdaptiveJouleError::AttemptBudget));
    let failure=advance_joule_enthalpy_adaptive(&old,1.0,request,&[],&gate,
        AdaptiveJouleConfig {minimum_step_s:0.75,..controls(1e-12)},guess,
        bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).err().expect("interval must refuse");
    assert!(matches!(failure.cause,AdaptiveJouleError::StepFloor {nonlinear:None,..}));
    assert_eq!(old.time_s(),0.0);assert_eq!(old.energy_j(),0.0);
}

#[test]
fn g4_nonlinear_nonconvergence_can_refine_without_changing_the_material_law() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));let gate=CancelGate::new();
    let mut config=controls(0.01);config.step.max_newton_iterations=1;
    let request=JouleEnthalpyInput {max_current_solves:100_000,..input(&domain,&drive,&ids,&[])};
    let (next,r)=advance_joule_enthalpy_adaptive(&old,0.01,request,&[],&gate,config,guess,
        |_,_,_|EnthalpyTransport::try_new((),Some(1.0),0.0),no_links,no_interfaces,|s| {
            let h=s.cells()[&ids[1]].specific_h_j_kg;Ok([(ids[1],-h*h)].into())
        }).unwrap();
    assert!(r.rejected_intervals>0);assert!(r.step_attempts>r.accepted_times_s.len());
    assert!((next.cells()[&ids[1]].specific_h_j_kg-0.01_f64.tanh()).abs()<1e-7);
    near(next.energy_j(),r.external_heat_j);
    // External cooling is not added to the electrical work a second time.
    near(r.deposited_joule_j,0.01);assert!(r.external_heat_j<r.deposited_joule_j);
}

#[test]
fn g4_late_source_refusal_or_cancellation_returns_no_accepted_prefix() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));
    for cancel in [false,true] {
        let gate=CancelGate::new();
        let failure=advance_joule_enthalpy_adaptive(&old,1.0,input(&domain,&drive,&ids,&[]),&[],&gate,
            AdaptiveJouleConfig {maximum_step_s:0.2,..controls(0.01)},guess,
            |_,_,_|EnthalpyTransport::try_new((),Some(1.0),0.0),no_links,no_interfaces,|s| {
                if s.time_s()>0.3 {
                    if cancel {gate.request();} else {return Err(JouleEnthalpyError::Invalid("unsupported later source"));}
                }
                Ok(BTreeMap::new())
            }).err().expect("interval must refuse");
        assert!(failure.step_attempts>=4);assert!(failure.current_solves>0);
        if !cancel {assert!(matches!(failure.cause,AdaptiveJouleError::Step(
            ImplicitJouleError::Physical(JouleEnthalpyError::Invalid("unsupported later source")))));}
        assert_eq!(old.time_s(),0.0);assert_eq!(old.energy_j(),0.0);
    }
}
