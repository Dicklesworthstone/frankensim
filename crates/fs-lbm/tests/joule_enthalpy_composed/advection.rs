//! Actual current solve, original caloric chart, source owner and mass reservoirs.
//! The fixtures are synthetic prescribed through-flow, not a momentum/melt solver.
use super::*;
use fs_conduction::transient::deposition::advection::{Advected, MassBoundary};

fn cooling(rate: f64) -> Advected<Heater> {
    let flow = Advected::new(Heater::default(),vec![],0).unwrap();
    if rate==0.0 { return flow; }
    flow.with_boundaries(vec![
        MassBoundary::Inlet {id:100,cell:CellId(1),kg_s:rate,material_id:MaterialId(1),specific_h_j_kg:0.0},
        MassBoundary::Outlet {id:101,cell:CellId(1),kg_s:rate},
    ],2).unwrap()
}

#[test]
fn flowing_core_keeps_electrode_work_and_reservoir_enthalpy_separate() {
    let law=curve(100.0); let old=initial(&law);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap(); let setup=setup();
    let ids=[CellId(0),CellId(1),CellId(2)]; let flow=cooling(1.0);
    let (next,accepted,r)=advance_joule_enthalpy_with_material_transport_observed(
        &old,1.0,config(),input(&domain,&setup,&ids),&[],&CancelGate::new(),&flow,
        |_,cell,_|resolve(&law,cell),|_,_|Ok((0.0,()))).unwrap();
    // Independent one-cell recurrence with constant 80 W electrical + 20 W heater,
    // cold inlet and local-enthalpy outlet. Use the actual accepted duration grid.
    let nominal: f64=0.125/(1.0+0.125/(0.8/1.0));
    let mut time: f64=0.0; let mut h=0.0;
    while time<1.0 {
        let dt=nominal.min(1.0-time);
        h+=dt*(100.0-h);
        time=if dt==1.0-time {1.0} else {time+dt};
    }
    near(next.cells()[&ids[1]].specific_h_j_kg,h);
    near(r.joule.electrode_energy_j,80.0);
    near(r.joule.deposited_joule_j,80.0);
    near(accepted.inner().joules,20.0);
    assert!(accepted.audit().external_enthalpy_j<0.0);
    near(r.joule.thermal.external_heat_j,80.0+20.0+accepted.audit().external_enthalpy_j);
    near(next.energy_j()-old.energy_j(),r.joule.thermal.external_heat_j);
    near(accepted.audit().inlet_mass_kg,1.0);
    near(accepted.audit().outlet_mass_kg,1.0);
    near(next.mass_kg(),old.mass_kg());
    near(r.endpoint.cells[1].response.material().liquid_mass_fraction(),(h-20.0)/100.0);
    near(next.temperature_k(ids[1]).unwrap(),310.0);
    assert_eq!(accepted.audit().substeps,accepted.inner().steps);
    assert_eq!(r.joule.current_solves,accepted.inner().steps+1);
    assert_eq!(flow.inner(),&Heater::default());
    assert_eq!(flow.audit().substeps,0);
    assert!(next.cells().values().all(|c|c.boundary.source_w==0.0));
}

#[test]
fn data_only_phase_response_and_flow_change_actual_current_at_equal_temperature() {
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap(); let mut setup=setup();
    setup.electrodes[0].drive=ElectrodeDrive::Voltage(20.0);
    let ids=[CellId(0),CellId(1),CellId(2)];
    let mut observed=Vec::new();
    for (latent,rate,contrast) in [(100.0,0.0,1.0),(100.0,2.0,1.0),(200.0,2.0,1.0),(100.0,2.0,2.0)] {
        let law=curve(latent); let old=initial(&law); let flow=cooling(rate);
        let (next,accepted,r)=advance_joule_enthalpy_with_material_transport_observed(
            &old,0.75,config(),input(&domain,&setup,&ids),&[7],&CancelGate::new(),&flow,
            |_,cell,_| {
                let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
                let f=phase.liquid_mass_fraction();
                EnthalpyTransport::try_new(phase,Some(5.0*(1.0+contrast*f)),2.0*(1.0+f))
            },
            |patch,sides| {
                // Original two half-cell series conduction, material-resolved.
                let a=sides[0].response.conductivity_w_m_k();
                let b=sides[1].response.conductivity_w_m_k();
                Ok((patch.link.area_m2/(0.5/a+0.5/b),()))
            }).unwrap();
        let sample=&r.endpoint.cells[1];
        let f=sample.response.material().liquid_mass_fraction();
        let current=r.joule.endpoint.electric.electrodes[0].current_a;
        near(sample.temperature_k,310.0);
        assert!(f>0.0 && f<1.0);
        near(current,20.0/(5.0*(1.0+contrast*f)));
        near(sample.response.conductivity_w_m_k(),2.0*(1.0+f));
        near(sample.cell.specific_h_j_kg,next.cells()[&ids[1]].specific_h_j_kg);
        near(r.endpoint.links[0].conductance_w_k,
            1.0/(0.5/r.endpoint.cells[0].response.conductivity_w_m_k()+0.5/sample.response.conductivity_w_m_k()));
        near(r.joule.thermal.external_heat_j,r.joule.deposited_joule_j+accepted.inner().joules+accepted.audit().external_enthalpy_j);
        near(r.joule.thermal.residual_j,0.0);
        near(r.joule.electrical_defect_j,0.0);
        near(r.joule.transfer_defect_j,0.0);
        assert_eq!(next.links()[&7].conductance_w_k,0.0); // Authored template survives.
        assert_eq!(r.joule.current_solves,r.joule.thermal.substeps+1);
        observed.push((current,f));
    }
    assert!(observed[1].0>observed[0].0); // Through-flow cools; lower r raises I.
    assert!(observed[2].0>observed[1].0); // More latent storage changes phase/r.
    assert!(observed[3].0<observed[1].0); // Only material response data changed.
}

#[test]
fn flow_composition_keeps_current_work_budget_and_endpoint_failure_atomic() {
    let law=curve(100.0); let old=initial(&law);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap(); let setup=setup();
    let ids=[CellId(0),CellId(1),CellId(2)]; let flow=cooling(1.0); let gate=CancelGate::new();
    let mut request=input(&domain,&setup,&ids); request.max_current_solves=1;
    let result=advance_joule_enthalpy_with_material_transport_observed(
        &old,1.0,config(),request,&[],&gate,&flow,
        |_,cell,_|resolve(&law,cell),|_,_|Ok((0.0,())));
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::CurrentBudget))));
    let result=advance_joule_enthalpy_with_material_transport_observed(
        &old,0.0625,config(),input(&domain,&setup,&ids),&[],&gate,&flow,
        |_,cell,_| if cell.specific_h_j_kg>6.0 {Err(JouleEnthalpyError::Invalid("unsupported flow endpoint"))}
            else {resolve(&law,cell)},|_,_|Ok((0.0,())));
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::Invalid(_)))));
    near(old.energy_j(),0.0);
    assert_eq!(old.time_s(),0.0);
    assert_eq!(flow.inner(),&Heater::default());
    assert_eq!(flow.audit().substeps,0);
}
