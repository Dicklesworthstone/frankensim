#![cfg(feature = "implicit-enthalpy-coupling")]
//! Real current/contact solves inside the shared endpoint energy residual.
use std::collections::BTreeMap;
use fs_blake3::ContentHash;
use fs_conduction::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::*;
use fs_lbm::conjugate::enthalpy::transport::*;
use fs_lbm::conjugate::enthalpy::transport::interfaces::*;
use fs_lbm::conjugate::enthalpy::transport::interfaces::implicit::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_solver::NewtonKrylovConfig;

fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([21;32]), SolidLiquidPhase::Solid,
        [(-200.0,100.0),(1e6,1_000_300.0)].map(|(h,t)|EnthalpyPhaseKnot {
            specific_enthalpy_j_kg:h, temperature_k:t, liquid_mass_fraction:0.0, bulk_density_kg_m3:1.0,
        }).to_vec()).unwrap()
}
fn latent(l:f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([22;32]),
        [(-100.0,250.0,0.0),(20.0,310.0,0.0),(20.0+l,310.0,1.0),(1020.0+l,560.0,1.0)]
        .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn config() -> ImplicitJouleConfig {
    ImplicitJouleConfig { max_unknowns:8, max_newton_iterations:24, difference_step_j:1e-4,
        cell_energy_tolerance_j:1e-7, total_energy_tolerance_j:1e-7,
        newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-10,
            forcing_minimum:1e-4,linear_restart:8,..NewtonKrylovConfig::default()} }
}
fn setup(n:usize,drive:ElectrodeDrive) -> ElectricSetup {
    ElectricSetup { conductors:vec![], electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive},
        Electrode {lo:[n-1,0,0],hi:[n,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ] }
}
fn initial<'a>(law:&'a EquilibriumEnthalpyPhaseCurve,h:&[f64],g:Option<f64>,bath:f64)
    ->(DepositionSimulation<'a>,VoxelDomain,Vec<CellId>) {
    let n=h.len(); let domain=VoxelDomain::new(n,1,1,1.0).unwrap();
    // Deliberately different from sorted thermal unknown order.
    let ids:Vec<_>=(0..n).map(|c|CellId((n-c) as u64*10)).collect();
    let mut events:Vec<_>=h.iter().enumerate().map(|(c,&h)|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:ids[c],segment_id:c as u64,body_id:1,material_id:MaterialId(1),mass_kg:1.0,specific_h_j_kg:h,
        geometry:CellGeometry {start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],volume_m3:1.0,surface_area_m2:6.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,source_w:0.0,
            reservoir:(c==1 && bath>0.0).then_some(Reservoir {area_m2:1.0,conductance_w_k:bath,temperature_k:300.0})},
    })}).collect();
    if let Some(g)=g {for c in 1..n {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:c as u64,a:ids[c-1],b:ids[c],area_m2:1.0,conductance_w_k:g,
    })});}}
    let sim=DepositionSimulation::new(events,&[(MaterialId(1),law as &dyn ThermalLaw)],DepositionBudget {
        max_cells:8,max_links:8,max_events:16,max_substeps:1000,
    }).unwrap().advance_to(0.0,StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    (sim,domain,ids)
}
fn input<'a>(domain:&'a VoxelDomain,setup:&'a ElectricSetup,ids:&'a [CellId],contacts:&'a [ElectricalContact]) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup,cell_ids:ids,contacts,max_cells:8,max_current_solves:2000,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn no_links<F>(_: &ContactPatch, _: [&EnthalpyTransportSample<F>; 2])
    -> Result<(f64, ()), JouleEnthalpyError> { panic!("no thermal links") }
fn no_interfaces<F>(_: usize, _: &ElectricalContact, _: [&EnthalpyTransportSample<F>; 2])
    -> Result<ElectricalInterfaceResponse<()>, JouleEnthalpyError> { panic!("no electrical interfaces") }
fn near(a:f64,b:f64) {assert!((a-b).abs()<2e-7,"{a} != {b}");}

#[test]
fn g1_closed_transient_resistive_feedback_matches_backward_euler_and_refines() {
    let law=sensible(); let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0)); let gate=CancelGate::new();
    let run=|n:usize| {
        let mut state=old.clone();
        for k in 1..=n {
            let (next,r)=step_joule_enthalpy_implicit(&state,k as f64/n as f64,
                input(&domain,&drive,&ids,&[]),&[],&gate,config(),&[0.0;3],
                |_,cell,_|EnthalpyTransport::try_new((),Some(1.0+0.25*cell.specific_h_j_kg),0.0),
                no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
            assert!(r.nonlinear.converged); assert!(r.current_solves>2);
            assert_eq!(r.endpoint.time_s,next.time_s());
            near(next.energy_j()-state.energy_j(),r.deposited_joule_j);
            assert!(r.maximum_cell_residual_j<config().cell_energy_tolerance_j);
            state=next;
        }
        state.cells()[&ids[1]].specific_h_j_kg
    };
    near(run(1),4.0/3.0); // Explicit one-pass source would produce 1 J/kg.
    let exact=(0.25_f64.exp()-1.0)/0.25;
    let errors=[8,16,32].map(|n|(run(n)-exact).abs());
    for pair in errors.windows(2) {assert!((1.9..2.2).contains(&(pair[0]/pair[1])));}
    assert_eq!(old.time_s(),0.0); assert_eq!(old.energy_j(),0.0);
}

#[test]
fn g1_stiff_reservoir_and_nonlinear_resistance_close_without_explicit_cfl_steps() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,20.0);
    let drive=setup(3,ElectrodeDrive::Voltage(1.0));
    let (next,r)=step_joule_enthalpy_implicit(&old,1.0,input(&domain,&drive,&ids,&[]),&[],
        &CancelGate::new(),config(),&[0.0,0.05,0.0],
        |_,cell,_|EnthalpyTransport::try_new((),Some(1.0+0.2*cell.specific_h_j_kg),0.0),
        no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
    let h=2.0/(21.0*(1.0+(1.0_f64+0.8/21.0).sqrt()));
    near(next.cells()[&ids[1]].specific_h_j_kg,h);
    near(r.thermal.external_heat_j,h);
    near(r.deposited_joule_j-h,20.0*h);
    assert_eq!(r.thermal.substeps,1);
    assert_eq!(old.cells()[&ids[1]].boundary,next.cells()[&ids[1]].boundary);
}

#[test]
fn g3_material_data_changes_implicit_phase_feedback_at_identical_kelvin_and_survives_shutoff() {
    let gate=CancelGate::new();let mut fractions=Vec::new();
    for l in [100.0,200.0] {
        let law=latent(l);let (old,domain,ids)=initial(&law,&[0.0,40.0,0.0],None,0.0);
        let drive=setup(3,ElectrodeDrive::Current(4.0));
        let bulk=|_:usize,cell:&DepositedCell,_:f64| {
            let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
            let r=1.0+phase.liquid_mass_fraction();
            EnthalpyTransport::try_new(phase,Some(r),0.0)
        };
        let (next,r)=step_joule_enthalpy_implicit(&old,0.25,input(&domain,&drive,&ids,&[]),&[],&gate,
            config(),&[0.0;3],bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
        let expected=(40.0+4.0*(1.0-20.0/l))/(1.0-4.0/l);
        near(next.cells()[&ids[1]].specific_h_j_kg,expected);
        near(next.temperature_k(ids[1]).unwrap(),310.0);
        near(r.deposited_joule_j,expected-40.0);
        fractions.push(r.material.cells[1].response.material().liquid_mass_fraction());
        let off=setup(3,ElectrodeDrive::Voltage(0.0));
        let (retained,r)=step_joule_enthalpy_implicit(&next,0.5,input(&domain,&off,&ids,&[]),&[],&gate,
            config(),&[0.0;3],bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new())).unwrap();
        near(retained.energy_j(),next.energy_j());near(r.deposited_joule_j,0.0);
    }
    assert!(fractions[0]>fractions[1]);
}

#[test]
fn g3_current_contacts_transport_and_other_heat_share_the_converged_endpoint() {
    let law=latent(100.0);let (old,domain,ids)=initial(&law,&[0.0,40.0,40.0,0.0],Some(0.5),0.0);
    let drive=setup(4,ElectrodeDrive::Voltage(6.0));let gate=CancelGate::new();
    let contacts=[ElectricalContact {cells:[1,2],resistance_ohm_m2:3.0,heat_fraction_to_first:0.25}];
    let request=input(&domain,&drive,&ids,&contacts);let links=[1,2,3];
    let (next,r)=step_joule_enthalpy_implicit(&old,0.125,request,&links,&gate,config(),&[0.0;4],
        |_,cell,_| {
            let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
            let f=phase.liquid_mass_fraction(); EnthalpyTransport::try_new(phase,Some(1.0+f),1.0+0.1*f)
        },
        |_,s|Ok((0.5*(s[0].response.conductivity_w_m_k()+s[1].response.conductivity_w_m_k()),())),
        |_,_,s|ElectricalInterfaceResponse::try_new((),3.0+s[0].response.material().liquid_mass_fraction()
            +2.0*s[1].response.material().liquid_mass_fraction(),0.25),
        |_|Ok([(ids[1],2.0)].into())).unwrap();
    near(next.energy_j()-old.energy_j(),r.deposited_joule_j+0.25);
    assert_eq!(r.interfaces.time_s,r.material.time_s);
    assert_eq!(r.interfaces.time_s,r.endpoint.time_s);
    let applied:Vec<_>=r.interfaces.contacts.iter().map(|s|s.applied()).collect();
    let fresh=drive.solve_with_cell_contacts(&domain,&r.endpoint.temperature_k,&applied,
        |c,_,_|Ok(r.endpoint.resistivity_ohm_m[c]),&gate).unwrap();
    near(fresh.total_joule_w,r.endpoint.electric.total_joule_w);
    for (id,old_cell) in old.cells() {assert_eq!(next.cells()[id].mass_kg,old_cell.mass_kg);}
    assert_eq!(next.links(),old.links());
    let increments:Vec<_>=old.cells().iter().map(|(id,c)|c.mass_kg*(next.cells()[id].specific_h_j_kg-c.specific_h_j_kg)).collect();
    let stage=BackwardEulerStep::new(&old,0.125,||Ok(())).unwrap();
    let trial=stage.evaluate(&increments,|_| {
        let mut sources:BTreeMap<_,_>=ids.iter().copied().zip(fresh.joule_w.iter().copied()).collect();
        *sources.get_mut(&ids[1]).unwrap()+=2.0;
        Ok::<_,DepositionError>(StepCoefficients {sources_w:sources,
            conductances_w_k:r.material.links.iter().map(|s|(s.patch.link.id,s.conductance_w_k)).collect(),payload:()})
    },||Ok(())).unwrap();
    assert!(trial.maximum_cell_residual_j()<config().cell_energy_tolerance_j);
}

#[test]
fn g4_budget_material_and_nonlinear_refusals_never_advance_thermal_history() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));let gate=CancelGate::new();
    let request=input(&domain,&drive,&ids,&[]);
    let no_links=|_:&ContactPatch,_:[&EnthalpyTransportSample<()>;2]|->Result<(f64,()),JouleEnthalpyError>{panic!("no links")};
    let no_interfaces=|_:usize,_:&ElectricalContact,_:[&EnthalpyTransportSample<()>;2]|->Result<ElectricalInterfaceResponse<()>,JouleEnthalpyError>{panic!("no interfaces")};
    let bulk=|_:usize,_:&DepositedCell,_:f64|EnthalpyTransport::try_new((),Some(1.0),0.0);
    let budget=step_joule_enthalpy_implicit(&old,1.0,JouleEnthalpyInput {max_current_solves:1,..request},&[],&gate,
        config(),&[0.0;3],bulk,no_links,no_interfaces,|_|Ok(BTreeMap::new()));
    assert!(matches!(budget,Err(ImplicitJouleError::Physical(JouleEnthalpyError::CurrentBudget))));
    let nonlinear=step_joule_enthalpy_implicit(&old,1.0,request,&[],&gate,
        ImplicitJouleConfig {max_newton_iterations:1,..config()},&[0.0;3],bulk,no_links,no_interfaces,|s| {
            let h=s.cells()[&ids[1]].specific_h_j_kg; Ok([(ids[1],-h*h)].into())
        });
    assert!(matches!(nonlinear,Err(ImplicitJouleError::NotConverged(_))));
    let unsupported=step_joule_enthalpy_implicit(&old,1.0,request,&[],&gate,
        config(),&[0.0;3],|c,cell,t| {
            if c==1 && cell.specific_h_j_kg>0.01 {Err(JouleEnthalpyError::Invalid("source domain"))} else {bulk(c,cell,t)}
        },no_links,no_interfaces,|_|Ok(BTreeMap::new()));
    assert!(matches!(unsupported,Err(ImplicitJouleError::Physical(JouleEnthalpyError::Invalid("source domain")))));
    let cancelled=CancelGate::new();
    assert!(step_joule_enthalpy_implicit(&old,1.0,request,&[],&cancelled,config(),&[0.0;3],|c,cell,t| {
        cancelled.request();bulk(c,cell,t)
    },no_links,no_interfaces,|_|Ok(BTreeMap::new())).is_err());
    assert_eq!(old.time_s(),0.0);assert_eq!(old.energy_j(),0.0);
}

#[test]
fn g0_loose_newton_stopping_cannot_bypass_actual_local_energy_acceptance() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(1.0));let mut loose=config();
    loose.newton.absolute_tolerance=1e9;
    let result=step_joule_enthalpy_implicit(&old,1.0,input(&domain,&drive,&ids,&[]),&[],&CancelGate::new(),loose,&[0.0;3],
        |_,_,_|EnthalpyTransport::try_new((),Some(1.0),0.0),no_links,
        no_interfaces,|_|Ok(BTreeMap::new()));
    assert!(matches!(result,Err(ImplicitJouleError::Physical(JouleEnthalpyError::Thermal(DepositionError::EnergyBalance(_))))));
    assert_eq!(old.energy_j(),0.0);
}

#[test]
fn g0_insulated_zero_current_terminal_retains_undefined_voltage_and_zero_work() {
    let law=sensible();let (old,domain,ids)=initial(&law,&[0.0;3],None,0.0);
    let drive=setup(3,ElectrodeDrive::Current(0.0));
    let (next,r)=step_joule_enthalpy_implicit(&old,0.125,input(&domain,&drive,&ids,&[]),&[],
        &CancelGate::new(),config(),&[0.0;3],
        |_,_,_|EnthalpyTransport::try_new((),None,0.0),no_links,no_interfaces,
        |_|Ok(BTreeMap::new())).unwrap();
    assert!(r.endpoint.electric.electrodes[0].voltage_v.is_nan());
    assert_eq!(r.endpoint.electric.electrodes[0].current_a,0.0);
    assert_eq!(r.electrode_energy_j,0.0);
    assert_eq!(r.joule_energy_j,0.0);
    assert_eq!(next.cells(),old.cells());
    assert_eq!(next.time_s(),0.125);
}
