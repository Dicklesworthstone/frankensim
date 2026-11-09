use super::*;
use super::conduction::PureFluidIsochoricConductionLaw;
use crate::thermal::conduction::ConductivityCurve;
use crate::thermochemical_gas::{ThermochemicalGasModel,
    nonideal::PengRobinsonSpeciesBinding,tests::{species_in_domain,transport}};
use fs_qty::Pressure;
use fs_thermal::enthalpy::{Boundary,CellGeometry,CellId,DepositionBudget,
    DepositionSimulation,Event,MaterialId,StepConfig,TimedEvent};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionLaw;
use fs_thermochem::{Composition,CompositionBasis,SpeciesId,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,mixture::FrozenIdealGasMixtureModelV1};
use fs_thermochem::eos::peng_robinson::saturation::PengRobinsonSaturationOptionsV1;

const MOLAR_MASS:f64=0.016_042_46;
fn model(formation:f64)->ThermochemicalPengRobinsonModel {
    let source=species_in_domain("A",MOLAR_MASS,4.0,0.0,formation,[100.0,250.0]);
    let ideal=FrozenIdealGasMixtureModelV1::new(vec![source],
        Composition::new(CompositionBasis::MoleFraction,vec![1.0]).unwrap()).unwrap();
    let bank=ThermochemicalGasModel::new(ideal,vec![transport("A",1.7e-5,0.025)],None).unwrap();
    ThermochemicalPengRobinsonModel::new(bank,vec![PengRobinsonSpeciesBinding {
        species:SpeciesId::new("A").unwrap(),critical_temperature:Temperature::new(190.564),
        critical_pressure:Pressure::new(4.5992e6),acentric_factor:0.01142}],vec![]).unwrap()
}
fn controls()->PureFluidEquilibriumControls {
    PureFluidEquilibriumControls::new(
        [Temperature::new(140.0),Temperature::new(160.0)],
        [Pressure::new(1e5),Pressure::new(4.4e6)])
}
fn density(model:&ThermochemicalPengRobinsonModel)->Density {
    let sat=model.mechanics().pure_saturation_at_temperature(Temperature::new(150.0),
        controls().pressure_bracket,PengRobinsonSaturationOptionsV1 {
            log_fugacity_tolerance:1e-12,relative_pressure_tolerance:1e-12,..Default::default()
        }).unwrap();
    Density::new(sat.at_vapor_fraction(0.35).unwrap().mass_density_kg_m3())
}
fn storage(model:&ThermochemicalPengRobinsonModel)->PureFluidIsochoricThermalStorage<'_> {
    model.prepare_isochoric_thermal_storage(density(model),controls(),||false).unwrap()
}
fn close(a:f64,b:f64,tolerance:f64) {
    assert!((a-b).abs()<=tolerance,"{a:.17e} != {b:.17e} (allowance {tolerance})");
}
fn boundary(power:f64)->Boundary {
    Boundary{ambient_k:150.0,htc_w_m2_k:0.0,reservoir:None,source_w:power}
}
fn cell(law:&PureFluidIsochoricThermalStorage<'_>,u:f64,power:f64)->DepositedCell {
    let mass=0.2;
    DepositedCell{id:CellId(1),segment_id:1,body_id:1,material_id:MaterialId(1),
        mass_kg:mass,specific_h_j_kg:u,geometry:CellGeometry {
            start_m:[0.0;3],end_m:[1.0,0.0,0.0],
            volume_m3:mass/law.density().value(),surface_area_m2:1.0},
        boundary:boundary(power)}
}

#[test]
fn g1_original_uv_storage_has_positive_moving_phase_capacity_and_fixed_density() {
    let model=model(0.0);let original=model.clone();let law=storage(&model);
    assert_eq!(law.storage_basis(),ThermalStorageBasis::InternalEnergy);
    let [lo,hi]=law.specific_internal_energy_domain_j_kg();
    assert_eq!(lo,0.0);assert!(hi>0.0);
    assert_eq!(law.preparation_temperature_evaluations(),2);
    assert!(law.preparation_eos_evaluations()<=controls().maximum_eos_evaluations);
    close(law.maximum_temperature_slope(),MOLAR_MASS/(3.0*R),2e-16);
    let mut previous=None;
    for f in [0.0,0.1,0.5,0.9,1.0] {
        let u=f*hi;let result=law.resolve_specific_internal_energy(u,||false).unwrap();
        let material=result.material();
        assert_eq!(material.density_kg_m3(),law.density().value());
        assert_eq!(law.bulk_density_kg_m3(u).unwrap(),Some(law.density().value()));
        assert!(material.vapor_mass_fraction()>0.0&&material.vapor_mass_fraction()<1.0);
        assert_eq!(material.reference(),&model.source().thermodynamics()
            .evaluate(Temperature::new(material.temperature_k()),Pressure::new(material.pressure_pa())).unwrap());
        assert!(result.temperature_slope_kg_k_j()>0.0);
        assert!(result.temperature_slope_kg_k_j()<=law.maximum_temperature_slope());
        assert!(result.pressure_slope_pa_kg_j()>0.0);
        assert!(result.solution().residual_j_kg().abs()<=controls().energy_tolerance_j_kg);
        close(result.solution().accepted_increment_j_kg().unwrap(),u,controls().energy_tolerance_j_kg);
        if let Some((t,p))=previous {
            assert!(material.temperature_k()>t);assert!(material.pressure_pa()>p);
        }
        previous=Some((material.temperature_k(),material.pressure_pa()));
    }
    assert_eq!(law.resolve_specific_internal_energy(0.0,||false).unwrap().solution().eos_evaluations(),0);
    assert!(law.advected_composition().is_none());
    assert_eq!(model,original);
}

#[test]
fn g1_original_thermal_timeline_heats_and_cools_a_sealed_two_phase_cell() {
    let model=model(0.0);let law=storage(&model);
    let hi=law.specific_internal_energy_domain_j_kg()[1];
    let u0=0.2*hi;let u1=0.8*hi;let heat=0.2*(u1-u0);
    let born=cell(&law,u0,heat);
    let timeline=DepositionSimulation::new(vec![
        TimedEvent{time_s:0.0,event:Event::Birth(born.clone())},
        TimedEvent{time_s:1.0,event:Event::Boundary(CellId(1),boundary(-heat))}],
        &[(MaterialId(1),&law as &dyn ThermalLaw)],
        DepositionBudget{max_cells:1,max_links:0,max_events:2,max_substeps:128}).unwrap();
    let step=StepConfig{max_dt_s:0.125,energy_tolerance_j:1e-5};
    let (initial,_)=timeline.advance_to(0.0,step,||Ok(())).unwrap();
    let t0=initial.temperature_k(CellId(1)).unwrap();
    let (hot,audit)=initial.advance_to(1.0,step,||Ok(())).unwrap();
    assert!(hot.temperature_k(CellId(1)).unwrap()>t0);
    close(audit.stored_energy_change_j,heat,2e-8);
    close(audit.external_heat_j,heat,2e-8);
    assert_eq!(hot.cells()[&CellId(1)].geometry,born.geometry);
    assert_eq!(hot.mass_kg(),born.mass_kg);
    let first=law.resolve_specific_internal_energy(u0,||false).unwrap();
    let last=law.resolve_specific_internal_energy(hot.cells()[&CellId(1)].specific_h_j_kg,||false).unwrap();
    assert!(last.material().pressure_pa()>first.material().pressure_pa());
    assert_ne!(last.material().vapor_mass_fraction(),first.material().vapor_mass_fraction());
    let (cooled,audit)=hot.advance_to(2.0,step,||Ok(())).unwrap();
    close(cooled.temperature_k(CellId(1)).unwrap(),t0,2e-8);
    close(audit.stored_energy_change_j,-heat,2e-8);
    assert_eq!(initial.time_s(),0.0);
    let mut wrong=born;wrong.geometry.volume_m3*=1.001;
    assert!(law.validate_cell(&wrong).is_err());
}

#[test]
fn g1_implicit_heat_uses_equilibrium_cv_and_zero_rigid_volume_derivative() {
    let model=model(0.0);let law=storage(&model);
    let conduction=PureFluidIsochoricConductionLaw::try_new(&law,
        ConductivityCurve::new(vec![[140.0,0.1],[160.0,0.12]]).unwrap()).unwrap();
    let hi=law.specific_internal_energy_domain_j_kg()[1];let u=0.5*hi;
    let state=law.resolve_specific_internal_energy(u,||false).unwrap();
    let thermal=conduction.state_at_specific_internal_energy(u).unwrap();
    close(thermal.temperature_slope_kg_k_j,MOLAR_MASS/state.response().molar_cv_j_mol_k,1e-18);
    let step=1.0;
    let minus=law.resolve_specific_internal_energy(u-step,||false).unwrap();
    let plus=law.resolve_specific_internal_energy(u+step,||false).unwrap();
    close((plus.material().temperature_k()-minus.material().temperature_k())/(2.0*step),
        thermal.temperature_slope_kg_k_j,1e-9);
    close((plus.material().pressure_pa()-minus.material().pressure_pa())/(2.0*step),
        state.pressure_slope_pa_kg_j(),1e-4);
    close(thermal.conductivity_slope_kg_w_m_k_j,0.001*thermal.temperature_slope_kg_k_j,1e-18);
    assert_eq!(conduction.specific_volume_derivative_at_specific_enthalpy(u).unwrap(),Some(0.0));
    let observed=conduction.time_observables(u).unwrap();
    assert_eq!(observed.temperature_k,state.material().temperature_k());
    assert_eq!(observed.phase_mass_fraction,Some(state.material().vapor_mass_fraction()));
    assert!(conduction.admit_harmonic_operating_point(u).is_ok());
    let a=0.2*hi;let b=0.8*hi;
    let entropy=law.thermal_entropy_change_j_kg_k(a,b,||false).unwrap();
    let source=law.specific_source_entropy_change_j_kg_k(a,b,||false).unwrap();
    close(entropy,source,3e-8);
    let ta=law.temperature_k(a).unwrap();let tb=law.temperature_k(b).unwrap();
    assert!(entropy/(b-a)>=1.0/tb&&entropy/(b-a)<=1.0/ta);
}

#[test]
fn g3_large_formation_reference_does_not_erase_weak_rigid_cell_heating() {
    let ordinary=model(0.0);let shifted=model(1e14);
    let a=storage(&ordinary);let b=storage(&shifted);
    assert_eq!(a.specific_internal_energy_domain_j_kg(),b.specific_internal_energy_domain_j_kg());
    let u=0.6*a.specific_internal_energy_domain_j_kg()[1];
    for target in [u,u+0.001] {
        let x=a.resolve_specific_internal_energy(target,||false).unwrap();
        let y=b.resolve_specific_internal_energy(target,||false).unwrap();
        assert_eq!(x.material().temperature_k(),y.material().temperature_k());
        assert_eq!(x.material().pressure_pa(),y.material().pressure_pa());
        assert_eq!(x.material().vapor_mass_fraction(),y.material().vapor_mass_fraction());
        assert_eq!(x.response(),y.response());
        assert_eq!(x.solution().residual_j_kg(),y.solution().residual_j_kg());
        assert!(y.material().specific_internal_energy_j_kg()>1e13);
    }
    let x=b.resolve_specific_internal_energy(u,||false).unwrap();
    let y=b.resolve_specific_internal_energy(u+0.001,||false).unwrap();
    assert!(y.material().temperature_k()>x.material().temperature_k());
    let ds=a.specific_source_entropy_change_j_kg_k(u,u+0.001,||false).unwrap();
    assert_eq!(ds,b.specific_source_entropy_change_j_kg_k(u,u+0.001,||false).unwrap());
    close(ds,0.001/x.material().temperature_k(),3e-8);
}

#[test]
fn g4_phase_departure_source_work_and_cancellation_refuse_without_state_change() {
    let model=model(0.0);let original=model.clone();let law=storage(&model);
    let hi=law.specific_internal_energy_domain_j_kg()[1];
    for u in [-f64::MIN_POSITIVE,hi.next_up(),f64::NAN,f64::INFINITY] {
        assert!(law.resolve_specific_internal_energy(u,||false).is_err());
    }
    let sat=model.mechanics().pure_saturation_at_temperature(Temperature::new(150.0),
        controls().pressure_bracket,controls().saturation_options).unwrap();
    let endpoint=Density::new(sat.liquid().mass_density_kg_per_m3());
    assert!(model.prepare_isochoric_thermal_storage(endpoint,controls(),||false).is_err());
    let mut too_small=controls();too_small.maximum_eos_evaluations=4;
    assert!(model.prepare_isochoric_thermal_storage(law.density(),too_small,||false).is_err());
    let mut unsupported=controls();unsupported.temperature_bracket[0]=Temperature::new(99.0);
    assert!(model.prepare_isochoric_thermal_storage(law.density(),unsupported,||false).is_err());
    let u=0.4*hi;let mut count=0;
    let expected=law.resolve_specific_internal_energy(u,||{count+=1;false}).unwrap();
    let mut calls=0;
    assert!(matches!(law.resolve_specific_internal_energy(u,||{calls+=1;calls==count}),Err(E::Cancelled)));
    assert_eq!(count,calls);
    assert_eq!(law.resolve_specific_internal_energy(u,||false).unwrap(),expected);
    assert_eq!(model,original);
}
