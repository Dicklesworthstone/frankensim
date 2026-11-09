use super::*;
use super::conduction::PureFluidIsobaricConductionLaw;
use crate::thermal::conduction::ConductivityCurve;
use crate::thermochemical_gas::{ThermochemicalGasModel,
    nonideal::PengRobinsonSpeciesBinding,
    tests::{point, species_with_energy_reference, transport}};
use fs_evidence::ValidityDomain;
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionSimulation, Event, MaterialId, StepConfig, TimedEvent};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionLaw;
use fs_thermochem::{Composition, CompositionBasis, Nasa9RegionV1, Nasa9StandardStateModelV1,
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1};

const MOLAR_MASS:f64=0.016_042_46;
fn model(formation:f64)->ThermochemicalPengRobinsonModel {
    let old=species_with_energy_reference("A",MOLAR_MASS,4.0,0.0,formation);
    let mut card=old.regions()[0].card().clone();
    card.validity=ValidityDomain::unconstrained().with("T",100.0,250.0);
    let source=Nasa9StandardStateModelV1::new(old.species().clone(),old.molar_mass(),
        old.convention().clone(),vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap();
    let ideal=FrozenIdealGasMixtureModelV1::new(vec![source],
        Composition::new(CompositionBasis::MoleFraction,vec![1.0]).unwrap()).unwrap();
    let bank=ThermochemicalGasModel::new(ideal,vec![transport("A",1.7e-5,0.025)],None).unwrap();
    ThermochemicalPengRobinsonModel::new(bank,vec![PengRobinsonSpeciesBinding {
        species:SpeciesId::new("A").unwrap(),critical_temperature:Temperature::new(190.564),
        critical_pressure:Pressure::new(4.5992e6),acentric_factor:0.01142}],vec![]).unwrap()
}
fn storage(model:&ThermochemicalPengRobinsonModel)->PureFluidIsobaricThermalStorage<'_> {
    model.prepare_isobaric_thermal_storage(Temperature::new(150.0),
        [Pressure::new(1e5),Pressure::new(4.4e6)],
        [Temperature::new(120.0),Temperature::new(180.0)],
        PureFluidIsobaricSolveV1::default(),PengRobinsonSaturationOptionsV1 {
            log_fugacity_tolerance:1e-12,relative_pressure_tolerance:1e-12,..Default::default()},||false).unwrap()
}
fn close(a:f64,b:f64,tolerance:f64) {
    assert!((a-b).abs()<=tolerance,"{a:.17e} != {b:.17e}");
}
fn boundary(power:f64)->Boundary {
    Boundary {ambient_k:150.0,htc_w_m2_k:0.0,reservoir:None,source_w:power}
}

#[test]
fn g1_source_bound_latent_storage_preserves_material_roots_receipts_and_volume() {
    let model=model(-1e6); let before=model.clone(); let law=storage(&model);
    assert_eq!(law.storage_basis(),ThermalStorageBasis::Enthalpy);
    let l=law.coexistence_specific_enthalpy_domain_j_kg()[1];
    close(l,6621.875_622_767_658/MOLAR_MASS,3e-6);
    assert!(law.maximum_temperature_slope()>=MOLAR_MASS/(3.0*R));
    assert!(law.maximum_temperature_slope()<=MOLAR_MASS/(3.0*R)*(1.0+1e-12));
    let sat=law.original_calorics().saturation();
    for q in [0.0,0.2,0.5,0.8,1.0] {
        let h=q*l;
        let result=law.resolve_specific_enthalpy(h,||false).unwrap();
        let material=result.material();
        assert_eq!(material.temperature_k(),150.0);
        assert_eq!(material.pressure_pa(),law.pressure().value());
        close(material.vapor_mass_fraction(),q,3e-16);
        assert_eq!(material.reference(),&model.source().thermodynamics()
            .evaluate(Temperature::new(150.0),law.pressure()).unwrap());
        let expected=(1.0-q)*sat.liquid().molar_volume_m3_per_mol()
            +q*sat.vapor().molar_volume_m3_per_mol();
        close(MOLAR_MASS/material.density_kg_m3(),expected,4e-19);
        assert_eq!(law.bulk_density_kg_m3(h).unwrap(),Some(material.density_kg_m3()));
        assert_eq!(result.caloric().state.temperature_slope_kg_k_j,0.0);
        assert!(result.caloric().residual_j_kg.abs()<=1e-6);
    }
    // Existing gas transport claims do not cover these fluid temperatures and
    // must not be repurposed as a liquid/two-phase transport model.
    assert!(model.source().evaluate(&point(150.0,law.pressure().value())).is_err());
    assert!(law.advected_composition().is_none());
    assert_eq!(model,before);
}

#[test]
fn g1_existing_thermal_timeline_crosses_boiling_and_condensation_conservatively() {
    let model=model(0.0); let law=storage(&model);
    let latent=law.coexistence_specific_enthalpy_domain_j_kg()[1];
    let mass=0.2; let h0=-0.1*latent; let h1=1.1*latent;
    let heat=mass*(h1-h0);
    let initial_temperature=law.temperature_k(h0).unwrap();
    let initial_volume=mass/law.bulk_density_kg_m3(h0).unwrap().unwrap();
    let born=DepositedCell {id:CellId(1),segment_id:1,body_id:1,material_id:MaterialId(1),
        mass_kg:mass,specific_h_j_kg:h0,geometry:CellGeometry {
            start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:initial_volume,surface_area_m2:1.0},
        boundary:boundary(heat)};
    let timeline=DepositionSimulation::new(vec![
        TimedEvent{time_s:0.0,event:Event::Birth(born.clone())},
        TimedEvent{time_s:1.0,event:Event::Boundary(CellId(1),boundary(-heat))}],
        &[(MaterialId(1),&law as &dyn ThermalLaw)],
        DepositionBudget{max_cells:1,max_links:0,max_events:2,max_substeps:128}).unwrap();
    let step=StepConfig{max_dt_s:0.125,energy_tolerance_j:1e-5};
    let (initial,_)=timeline.advance_to(0.0,step,||Ok(())).unwrap();
    let (half,audit)=initial.advance_to(0.5,step,||Ok(())).unwrap();
    assert_eq!(half.temperature_k(CellId(1)).unwrap(),150.0);
    close(audit.stored_energy_change_j,0.5*heat,1e-8);
    let half_state=law.resolve_specific_enthalpy(half.cells()[&CellId(1)].specific_h_j_kg,||false).unwrap();
    close(half_state.material().vapor_mass_fraction(),0.5,2e-15);
    assert!(mass/half_state.material().density_kg_m3()>initial_volume);
    assert_eq!(half.cells()[&CellId(1)].geometry,born.geometry); // reference mesh only
    let (hot,audit)=initial.advance_to(1.0,step,||Ok(())).unwrap();
    assert!(hot.temperature_k(CellId(1)).unwrap()>150.0);
    close(audit.stored_energy_change_j,heat,2e-8);
    close(audit.external_heat_j,heat,2e-8);
    let (cooled,audit)=hot.advance_to(2.0,step,||Ok(())).unwrap();
    close(cooled.temperature_k(CellId(1)).unwrap(),initial_temperature,2e-8);
    close(audit.stored_energy_change_j,-heat,2e-8);
    assert_eq!(cooled.mass_kg(),mass);
    assert_eq!(initial.time_s(),0.0);
}

#[test]
fn g1_implicit_response_retains_zero_thermal_tangent_and_nonzero_latent_expansion() {
    let model=model(0.0);let law=storage(&model);
    let conductor=PureFluidIsobaricConductionLaw::try_new(&law,
        ConductivityCurve::new(vec![[120.0,0.1],[180.0,0.22]]).unwrap()).unwrap();
    let latent=law.coexistence_specific_enthalpy_domain_j_kg()[1];
    for q in [0.0,0.3,0.7,1.0] {
        let state=conductor.state_at_specific_enthalpy(q*latent).unwrap();
        assert_eq!(state.temperature_k,150.0);
        assert_eq!(state.temperature_slope_kg_k_j,0.0);
        assert_eq!(state.conductivity_slope_kg_w_m_k_j,0.0);
        close(state.conductivity_w_m_k,0.16,1e-15);
        let dv=conductor.specific_volume_derivative_at_specific_enthalpy(q*latent).unwrap().unwrap();
        assert!(dv>0.0);
        let observables=conductor.time_observables(q*latent).unwrap();
        close(observables.phase_mass_fraction.unwrap(),q,3e-16);
        assert_eq!(observables.temperature_k,150.0);
    }
    // Latent endpoints retain valid one-sided finite-step tangents above, but
    // crossing either corner has no two-sided harmonic derivative.
    for h in [0.0,latent] {
        assert!(conductor.admit_harmonic_operating_point(h).is_err());
    }
    assert!(conductor.admit_harmonic_operating_point(0.5*latent).is_ok());
    assert!(conductor.admit_harmonic_operating_point(f64::NAN).is_err());
    let dh=0.6*latent;
    close(law.thermal_entropy_change_j_kg_k(0.2*latent,0.8*latent,||false).unwrap(),dh/150.0,2e-12);
    let source=law.specific_source_entropy_change_j_kg_k(0.2*latent,0.8*latent,||false).unwrap();
    close(dh-150.0*source,0.6*law.latent_closure_residual_j_kg(),3e-10);
    let sensible=law.caloric_state_at_temperature(Temperature::new(170.0),||false).unwrap();
    let response=conductor.state_at_specific_enthalpy(sensible.specific_enthalpy_j_kg).unwrap();
    assert!(conductor.admit_harmonic_operating_point(sensible.specific_enthalpy_j_kg).is_ok());
    close(response.conductivity_slope_kg_w_m_k_j,0.002*response.temperature_slope_kg_k_j,1e-18);
    assert!(conductor.specific_volume_derivative_at_specific_enthalpy(sensible.specific_enthalpy_j_kg)
        .unwrap().unwrap()>0.0);
}

#[test]
fn g3_large_formation_datums_leave_the_storage_coordinate_and_heat_unchanged() {
    let ordinary=model(0.0);let shifted=model(1e14);
    let a=storage(&ordinary);let b=storage(&shifted);
    assert_eq!(a.specific_enthalpy_domain_j_kg(),b.specific_enthalpy_domain_j_kg());
    let latent=a.coexistence_specific_enthalpy_domain_j_kg()[1];
    for h in [-0.05*latent,0.0,0.5*latent,latent,1.05*latent] {
        assert_eq!(a.caloric_state_at_specific_enthalpy(h,||false).unwrap(),
            b.caloric_state_at_specific_enthalpy(h,||false).unwrap());
        let x=b.resolve_specific_enthalpy(h,||false).unwrap();
        assert!(x.material().specific_enthalpy_j_kg()>1e13);
    }
}

#[test]
fn g4_source_domain_undetermined_amount_and_late_cancellation_refuse_atomically() {
    let model=model(0.0);let original=model.clone();let law=storage(&model);
    assert!(matches!(law.caloric_state_at_temperature(Temperature::new(150.0),||false),
        Err(E::Caloric(PureFluidIsobaricErrorV1::CoexistenceAmountRequired))));
    let [lo,hi]=law.specific_enthalpy_domain_j_kg();
    for h in [lo.next_down(),hi.next_up(),f64::NAN,f64::INFINITY] {
        assert!(law.resolve_specific_enthalpy(h,||false).is_err());
    }
    assert!(model.prepare_isobaric_thermal_storage(Temperature::new(150.0),
        [Pressure::new(1e5),Pressure::new(4.4e6)],
        [Temperature::new(99.0),Temperature::new(180.0)],Default::default(),Default::default(),||false).is_err());
    let h=0.4*law.coexistence_specific_enthalpy_domain_j_kg()[1];
    let mut count=0;
    let expected=law.resolve_specific_enthalpy(h,||{count+=1;false}).unwrap();
    let mut calls=0;
    assert!(matches!(law.resolve_specific_enthalpy(h,||{calls+=1;calls==count}),Err(E::Cancelled)));
    assert_eq!(count,calls);
    assert_eq!(law.resolve_specific_enthalpy(h,||false).unwrap(),expected);
    assert_eq!(model,original);
}
