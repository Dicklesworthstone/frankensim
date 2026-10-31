//! Real gas/EOS/caloric/source integration over the existing synthetic fixtures.
use super::*;
use super::super::reservoir::{IsothermalGasMembrane,IsothermalGasMembraneSolve,
    IsothermalGasMembraneError as GasMembraneError};
use super::super::slab::finite::FiniteMembraneDiffusionWorkspace;
use crate::thermochemical_gas::{ThermochemicalGasModel,tests::{species_with_energy_reference,
    transport_in_domain,point as gas_point,model as mixed_model}};
use crate::thermochemical_gas::energy::{GasEnergySolve,parcel::ClosedGasParcel};
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::GasSpeciesBoundarySolve;
use fs_thermochem::{Composition,CompositionBasis,mixture::FrozenIdealGasMixtureModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const M:f64=0.028;
fn solve()->IsothermalGasMembraneSolve {
    IsothermalGasMembraneSolve {membrane:controls(),gas:[GasSpeciesBoundarySolve {
        gas:GasEnergySolve {minimum_temperature_k:200.0,maximum_temperature_k:1000.0,
            energy_tolerance_j_kg:1e-7,maximum_iterations:128},
        mass_tolerance_kg:1e-17,energy_tolerance_j:1e-9};2],mass_tolerance_kg:1e-17}
}
fn gas_model(reference:f64,pressure_domain:[f64;2])->ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A",M,3.5,0.001,reference)],
        Composition::new(CompositionBasis::MoleFraction,vec![1.0]).unwrap()).unwrap(),
        vec![transport_in_domain("A",1.7e-5,0.025,[200.0,1000.0],pressure_domain)],None).unwrap()
}
fn gas(model:&ThermochemicalGasModel,pressure:f64,volume:f64,lot:f64)->ClosedGasParcel {
    ClosedGasParcel::new(model.clone(),gas_point(300.0,pressure).with("local-lot",lot).unwrap(),volume).unwrap()
}
fn resolved<'a>(card:&'a MaterialCard,q:&'a QueryPoint,loading:f64,pressures:[f64;2])->ResolvedMembraneDiffusion<'a> {
    ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(card,q,8)],&[],&[loading;8],pressures,||false).unwrap()
}
fn owner<'a>(card:&'a MaterialCard,q:&'a QueryPoint,model:&ThermochemicalGasModel,
    loading:f64,pressures:[f64;2])->IsothermalGasMembrane<'a>
{
    resolved(card,q,loading,pressures).with_finite_gas_reservoirs(
        [gas(model,pressures[0],1e-7,7.0),gas(model,pressures[1],3e-7,9.0)],solve(),||false).unwrap()
}
fn near(a:f64,b:f64,tolerance:f64) {assert!((a-b).abs()<=tolerance,"{a} != {b}");}

#[test]
fn actual_gases_deplete_into_storage_and_keep_original_eos_transport_and_sources() {
    let card=material("dense A",[0.0,2e5]);let q=point(50000.0);
    let model=gas_model(0.0,[1000.0,1e7]);let state=owner(&card,&q,&model,0.0,[1e5,10000.0]);
    let saved_gases=state.gases().clone();let saved_film=state.membrane().clone();
    let mut work=FiniteMembraneDiffusionWorkspace::new(8).unwrap();
    let next=state.step_candidate(1.0,&mut work,||false).unwrap();let report=next.step();
    assert!(report.sorbed_mass_change_kg>0.0);
    assert!(report.actual_gas_mass_changes_kg[0]<0.0);
    assert!(report.actual_gas_mass_changes_kg[1]<0.0);
    assert!(report.mass_residual_kg.abs()<solve().mass_tolerance_kg);
    assert!(report.total_mass_residual_kg.abs()<solve().mass_tolerance_kg);
    for i in 0..2 {
        let accepted=&next.state().gases()[i];let old=&saved_gases[i];
        assert_eq!(accepted.state().state().temperature,300.0);
        assert_eq!(accepted.volume_m3(),old.volume_m3());
        assert_eq!(accepted.point().axes().get("local-lot"),old.point().axes().get("local-lot"));
        assert_eq!(accepted.model().transport_materials(),old.model().transport_materials());
        assert_eq!(accepted.model().thermodynamics().components(),old.model().thermodynamics().components());
        assert_eq!(accepted.model().source_packs(),old.model().source_packs());
        assert_eq!(&accepted.model().evaluate(accepted.point()).unwrap(),accepted.state());
        near(accepted.state().state().density*accepted.volume_m3(),accepted.mass_kg(),1e-17);
        assert_eq!(accepted.state().state().pressure,next.state().finite_state().reservoir_pressures_pa()[i]);
        near(report.stream_enthalpy_into_gas_j[i]+report.thermostat_heat_into_gas_j[i],report.gas[i].energy_into_gas_j,1e-12);
        assert_eq!(report.gas[i].energy_residual_j,0.0);
        assert!(report.pressure_mismatch_mass_kg[i].abs()<=solve().gas[i].mass_tolerance_kg);
    }
    for receipt in &next.state().source_support().layers[0] {assert_eq!(receipt.card_identity(),card.content_hash());}
    assert_eq!(state.gases(),&saved_gases);assert_eq!(state.membrane(),&saved_film);
}

#[test]
fn huge_common_formation_reference_cancels_before_isothermal_reconstruction() {
    let card=material("dense A",[0.0,2e5]);let q=point(50000.0);let mut ordinary=None;
    for reference in [0.0,1e24] {
        let model=gas_model(reference,[1000.0,1e7]);let state=owner(&card,&q,&model,0.0,[1e5,10000.0]);
        let mut work=FiniteMembraneDiffusionWorkspace::new(8).unwrap();
        let next=state.step_candidate(5.0,&mut work,||false).unwrap();
        let values=(next.state().membrane().partial_pressures_pa().to_vec(),
            next.state().gases().each_ref().map(|g|g.mass_kg()),
            next.state().gases().each_ref().map(|g|g.state().state().pressure),next.step().thermostat_heat_into_gas_j);
        if let Some(expected)=&ordinary {assert_eq!(&values,expected);}else{ordinary=Some(values);}
        for gas in next.state().gases() {assert_eq!(gas.state().state().temperature,300.0);}
        assert_eq!(next.step().gas_energy_residual_j,0.0);
        assert!(next.step().mass_residual_kg.abs()<solve().mass_tolerance_kg);
        if reference>0.0 {
            assert!(next.step().stream_enthalpy_into_gas_j[0].abs()
                >next.step().thermostat_heat_into_gas_j[0].abs()*1e12);
        }
    }
}

#[test]
fn actual_gas_equilibrium_includes_sorbed_capacity_and_closes_mass_after_reseeding() {
    let card=material("dense A",[0.0,2e5]);let q=point(50000.0);let model=gas_model(0.0,[1000.0,1e7]);
    let mut state=owner(&card,&q,&model,0.0,[1e5,10000.0]);let initial=state.total_mass_kg().unwrap();
    let capacity=(state.gases()[0].volume_m3()+state.gases()[1].volume_m3())/(R*300.0)+0.01*2e-5*1e-4;
    let expected=(initial/M)/capacity;let mut converged=false;
    let mut work=FiniteMembraneDiffusionWorkspace::new(8).unwrap();
    for _ in 0..64 {
        let next=state.step_candidate(100.0,&mut work,||false).unwrap();
        assert!(next.step().mass_residual_kg.abs()<solve().mass_tolerance_kg);
        assert!(next.step().total_mass_residual_kg.abs()<solve().mass_tolerance_kg);
        state=next.into_parts().0;
        near(state.total_mass_kg().unwrap(),initial,1e-16);
        converged=state.gases().iter().map(|gas|gas.state().state().pressure)
            .chain(state.membrane().partial_pressures_pa().iter().copied())
            .all(|p|(p-expected).abs()<=1e-5*expected);
        if converged {break;}
    }
    assert!(converged,"finite actual gas reservoirs did not approach total-capacity equilibrium");
}

#[test]
fn unloading_and_reverse_permeation_use_signed_original_gas_streams() {
    let card=material("dense A",[0.0,2e5]);let q=point(50000.0);let model=gas_model(0.0,[1000.0,1e7]);
    let loaded=owner(&card,&q,&model,150000.0,[10000.0,20000.0]);
    let mut work=FiniteMembraneDiffusionWorkspace::new(8).unwrap();
    let release=loaded.step_candidate(1.0,&mut work,||false).unwrap();
    assert!(release.step().actual_gas_mass_changes_kg.into_iter().all(|dm|dm>0.0));
    assert!(release.step().thermostat_heat_into_gas_j.into_iter().all(|heat|heat<0.0));
    assert!(release.step().sorbed_mass_change_kg<0.0);
    let reversed=owner(&card,&q,&model,50000.0,[10000.0,100000.0]);
    let reverse=reversed.step_candidate(1.0,&mut work,||false).unwrap();
    assert!(reverse.step().actual_gas_mass_changes_kg[0]>0.0);
    assert!(reverse.step().actual_gas_mass_changes_kg[1]<0.0);
    assert!(reverse.step().mass_residual_kg.abs()<solve().mass_tolerance_kg);
}

#[test]
fn source_temperature_mixture_and_cancellation_refusals_preserve_accepted_owners() {
    let card=material("dense A",[0.0,2e5]);let q=point(50000.0);let model=gas_model(0.0,[1000.0,1e7]);
    let hot=ClosedGasParcel::new(model.clone(),gas_point(400.0,1e5),1e-7).unwrap();
    assert!(resolved(&card,&q,0.0,[1e5;2]).with_finite_gas_reservoirs(
        [gas(&model,1e5,1e-7,7.0),hot],solve(),||false).is_err());
    let mixed=ClosedGasParcel::new(mixed_model(vec![0.5,0.5],0.0),gas_point(300.0,1e5),1e-7).unwrap();
    assert!(resolved(&card,&q,0.0,[1e5;2]).with_finite_gas_reservoirs(
        [gas(&model,1e5,1e-7,7.0),mixed],solve(),||false).is_err());
    let changed=gas_model(1e6,[1000.0,1e7]);
    assert!(resolved(&card,&q,0.0,[1e5;2]).with_finite_gas_reservoirs(
        [gas(&model,1e5,1e-7,7.0),gas(&changed,1e5,1e-7,9.0)],solve(),||false).is_err());
    let bounded=gas_model(0.0,[90000.0,200000.0]);
    let state=resolved(&card,&q,0.0,[1e5;2]).with_finite_gas_reservoirs(
        [gas(&bounded,1e5,1e-8,7.0),gas(&bounded,1e5,1e-8,9.0)],solve(),||false).unwrap();
    let saved=state.gases().clone();let film=state.membrane().clone();
    let mut work=FiniteMembraneDiffusionWorkspace::new(8).unwrap();
    assert!(state.step_candidate(10000.0,&mut work,||false).is_err());
    assert_eq!(state.gases(),&saved);assert_eq!(state.membrane(),&film);
    let mut polls=0;
    assert!(matches!(state.step_candidate(1.0,&mut work,||{polls+=1;polls==20}),Err(GasMembraneError::Cancelled)));
    assert_eq!(state.gases(),&saved);assert_eq!(state.membrane(),&film);
}
