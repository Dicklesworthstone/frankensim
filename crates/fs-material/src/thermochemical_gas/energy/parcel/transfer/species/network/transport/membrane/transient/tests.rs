//! Original-source admission, local domains and immutable candidate boundaries.
//! Synthetic values test the implementation, not a named polymer's validity.
use super::*;
use crate::state_point::MaterialPropertySelection;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,UncertaintyModel,InterfaceSystemCard,SurfaceSpec,SystemContext};
use fs_qty::{Dims,Pressure,QuantitySpec,semantic::{QuantityKind,SemanticType,ValueForm}};

fn typed(kind:QuantityKind)->QuantitySpec {QuantitySpec::semantic(SemanticType::new(kind,ValueForm::Static))}
fn point(pressure:f64)->QueryPoint {
    QueryPoint::new().with_quantity("temperature",typed(QuantityKind::AbsoluteTemperature),300.0).unwrap()
        .with_quantity("pressure",typed(QuantityKind::Pressure),pressure).unwrap()
}
fn claim(name:&str,value:f64,dims:Dims,pressure:[f64;2])->PropertyClaim {
    PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
        validity:ValidityDomain::unconstrained()
            .with_quantity("temperature",typed(QuantityKind::AbsoluteTemperature),200.0,600.0)
            .with_quantity("pressure",typed(QuantityKind::Pressure),pressure[0],pressure[1]),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"synthetic dense membrane diffusion".into(),
            license:"CC0-1.0".into(),artifact:None}}
}
fn material(name:&str,pressure:[f64;2])->MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(claim(MEMBRANE_DIFFUSIVITY_PROPERTY,1e-10,MEMBRANE_DIFFUSIVITY_DIMS,pressure)).unwrap();
    claims.insert_claim(claim(MOLAR_HENRY_SOLUBILITY_PROPERTY,2e-5,MOLAR_HENRY_SOLUBILITY_DIMS,pressure)).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:name.into(),phase:"dense solid".into(),
        process:"analytic test".into(),revision:0},claims,vec![]).unwrap()
}
fn layer<'a>(material:&'a MaterialCard,context:&'a QueryPoint,cells:usize)->MembraneDiffusionLayerSource<'a> {
    MembraneDiffusionLayerSource {thickness_m:1e-4,cells,material:MembranePermeabilityMaterial {
        species:SpeciesId::new("A").unwrap(),material,context,selection:MaterialPropertySelection::SingleClaimOnly}}
}
fn interface(left:&MaterialCard,right:&MaterialCard,value:f64,pressure:[f64;2])->InterfaceSystemCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(claim(MOLAR_INTERFACE_PERMEANCE_PROPERTY,value,MOLAR_PERMEANCE_DIMS,pressure)).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {material:left.id().clone(),texture_frame:"left dense face".into()},
        SurfaceSpec {material:right.id().clone(),texture_frame:"right dense face".into()},
        SystemContext {medium:"dense-layer contact".into(),third_body:None,
            environment:"one pure permeant".into(),history:"unchanged".into()},claims,vec![]).unwrap()
}
fn bound_interface<'a>(interface:&'a InterfaceSystemCard,context:&'a QueryPoint)->MembraneInterfaceMaterial<'a> {
    MembraneInterfaceMaterial {species:SpeciesId::new("A").unwrap(),interface,context,
        selection:MaterialPropertySelection::SingleClaimOnly}
}
fn controls()->MembraneDiffusionSolve {
    MembraneDiffusionSolve {absolute_amount_tolerance_mol:1e-14,relative_amount_tolerance:1e-9}
}

#[test]
fn breakthrough_retains_both_original_coefficient_sources_and_actual_pressure_support() {
    let card=material("layer A",[0.0,1e5]);let context=point(50000.0);
    let state=ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(&card,&context,8)],&[],
        &[0.0;8],[1e5,0.0],||false).unwrap();
    assert_eq!(state.species(),&SpeciesId::new("A").unwrap());
    assert_eq!(state.temperature_k(),300.0);
    for (source,p) in state.source_support().layers[0].iter().zip([0.0,1e5]) {
        assert_eq!(source.card_identity(),card.content_hash());
        assert_eq!(source.state_coordinate("pressure"),Some(p));
        assert_eq!(source.property(MEMBRANE_DIFFUSIVITY_PROPERTY).unwrap().value_si(),1e-10);
        assert_eq!(source.property(MOLAR_HENRY_SOLUBILITY_PROPERTY).unwrap().value_si(),2e-5);
    }
    let mut work=MembraneDiffusionWorkspace::new(8).unwrap();
    let next=state.step_candidate([1e5,0.0],1.0,controls(),&mut work,||false).unwrap();
    assert!(next.step().stored_after_mol>0.0);
    assert!(next.step().left_reservoir_to_membrane_mol>next.step().membrane_to_right_reservoir_mol);
    assert_eq!(state.state().stored_amount_mol().unwrap(),0.0);
    for source in &next.state().source_support().layers[0] {
        assert_eq!(source.card_identity(),card.content_hash());
        for name in [MEMBRANE_DIFFUSIVITY_PROPERTY,MOLAR_HENRY_SOLUBILITY_PROPERTY] {
            assert_eq!(source.property(name).unwrap().answer().receipt.selected,
                state.source_support().layers[0][0].property(name).unwrap().answer().receipt.selected);
        }
    }
}

#[test]
fn unsupported_boundary_or_candidate_refuses_even_when_initial_cells_and_anchor_are_supported() {
    let context=point(50000.0);let narrow=material("narrow",[0.0,90000.0]);
    assert!(ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(&narrow,&context,4)],&[],
        &[0.0;4],[1e5,0.0],||false).is_err());
    let card=material("bounded",[0.0,1e5]);
    let state=ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(&card,&context,4)],&[],
        &[50000.0;4],[50000.0;2],||false).unwrap();
    let saved=state.state().clone();let mut work=MembraneDiffusionWorkspace::new(4).unwrap();
    assert!(state.step_candidate([2e5,0.0],1.0,controls(),&mut work,||false).is_err());
    assert_eq!(state.state(),&saved);
    assert!(matches!(state.step_candidate([1e5,0.0],1.0,controls(),&mut work,||true),Err(Error::Cancelled)));
    assert_eq!(state.state(),&saved);
}

#[test]
fn each_layer_and_interface_uses_its_actual_local_domain_and_ordered_identity() {
    let a=material("left",[49000.0,1e5]);let b=material("right",[0.0,51000.0]);
    let (qa,qb,qi)=(point(75000.0),point(25000.0),point(50000.0));
    let layers=[layer(&a,&qa,1),layer(&b,&qb,1)];
    let contact=interface(&a,&b,1e-9,[49000.0,51000.0]);
    let state=ResolvedMembraneDiffusion::new(300.0,0.01,&layers,
        &[Some(bound_interface(&contact,&qi))],&[75000.0,25000.0],[1e5,0.0],||false).unwrap();
    let source=state.source_support().interfaces[0].as_ref().unwrap();
    for side in source {
        assert_eq!(side.card_identity(),contact.content_hash());
        assert_eq!(side.surface_materials(),&[a.id().clone(),b.id().clone()]);
    }
    assert!(source[0].state_coordinate("pressure").unwrap()>49000.0);
    assert!(source[1].state_coordinate("pressure").unwrap()<51000.0);
    let mut work=MembraneDiffusionWorkspace::new(2).unwrap();
    state.step_candidate([1e5,0.0],1.0,controls(),&mut work,||false).unwrap();
    let wrong_order=interface(&b,&a,1e-9,[0.0,1e5]);
    assert!(ResolvedMembraneDiffusion::new(300.0,0.01,&layers,
        &[Some(bound_interface(&wrong_order,&qi))],&[75000.0,25000.0],[1e5,0.0],||false).is_err());
    let unsupported=interface(&a,&b,1e-9,[49900.0,50100.0]);
    assert!(ResolvedMembraneDiffusion::new(300.0,0.01,&layers,
        &[Some(bound_interface(&unsupported,&qi))],&[75000.0,25000.0],[1e5,0.0],||false).is_err());
}

#[test]
fn two_explicit_seals_keep_a_supported_middle_inventory_exactly() {
    let (a,b,c)=(material("a",[0.0,1e5]),material("b",[0.0,1e5]),material("c",[0.0,1e5]));
    let q=point(50000.0);let (ab,bc)=(interface(&a,&b,0.0,[0.0,1e5]),interface(&b,&c,0.0,[0.0,1e5]));
    let state=ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(&a,&q,1),layer(&b,&q,1),layer(&c,&q,1)],
        &[Some(bound_interface(&ab,&q)),Some(bound_interface(&bc,&q))],
        &[0.0,12345.0,0.0],[1e5,0.0],||false).unwrap();
    let mut work=MembraneDiffusionWorkspace::new(3).unwrap();
    let next=state.step_candidate([1e5,0.0],100.0,controls(),&mut work,||false).unwrap();
    assert_eq!(next.state().state().partial_pressures_pa()[1],12345.0);
    assert_eq!(&next.step().face_molar_fluxes_mol_m2_s[1..3],&[0.0;2]);
}

#[test]
fn absent_wrong_basis_and_pressure_dependent_sorption_are_not_inferred() {
    let base=material("specified",[0.0,1e5]);let q=point(50000.0);
    for defect in 0..3 {
        let mut claims=ClaimSet::new();
        claims.insert_claim(base.claims_for(MEMBRANE_DIFFUSIVITY_PROPERTY)[0].1.clone()).unwrap();
        if defect!=0 {
            let mut s=base.claims_for(MOLAR_HENRY_SOLUBILITY_PROPERTY)[0].1.clone();
            if defect==1 {
                s.key=PropertyKey::new(MOLAR_HENRY_SOLUBILITY_PROPERTY,MOLAR_PERMEANCE_DIMS);
                s.value=PropertyValue::Scalar {value:2e-5,dims:MOLAR_PERMEANCE_DIMS};
            }else{
                s.value=PropertyValue::Curve {abscissa:"pressure".into(),abscissa_dims:Pressure::DIMS,
                    knots:vec![(0.0,2e-5),(1e5,2e-5)],dims:MOLAR_HENRY_SOLUBILITY_DIMS};
                s.interpolation=InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(s).unwrap();
        }
        let card=MaterialCard::assemble(base.id().clone(),claims,vec![]).unwrap();
        assert!(ResolvedMembraneDiffusion::new(300.0,0.01,&[layer(&card,&q,1)],&[],
            &[0.0],[1e5,0.0],||false).is_err());
    }
}
