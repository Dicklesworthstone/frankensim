//! Same-specimen caloric/transport selection, not composition-derived Cp.
use super::*;
use fs_material::sheet::thermal::*;
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_PROPERTY as CP, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY as K, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind,
    SemanticType, ValueForm}};

fn thermal_card(interpolation:InterpolationPolicy) -> MaterialCard {
    let mut claims=card(&entries()[..10]).claims().clone();
    for (name,dims,values) in [(CP,SPECIFIC_HEAT_CAPACITY_DIMS,vec![(290.0,1000.0),(300.0,1200.0),(310.0,1300.0)]),
        (K,THERMAL_CONDUCTIVITY_DIMS,vec![(290.0,0.04),(310.0,0.08)])] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Dims([0,0,0,1,0,0]),knots:values,dims},
            interpolation,validity:ValidityDomain::unconstrained().with("T",290.0,310.0).with("specimen",7.0,7.0),
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"synthetic thermal specimen".into(),license:"test-only".into(),artifact:None}
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),phase:"solid".into(),
        process:"same specimen".into(),revision:0},claims,vec![]).unwrap()
}
fn bundle(c:&MaterialCard,specimen:f64) -> ResolvedSheet {
    let mut requirements=Vec::new();
    for (key,_,dims) in entries().into_iter().take(10) {
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let q=QueryPoint::new().with("T",300.0).unwrap().with("specimen",specimen).unwrap();
    ResolvedSheet::try_from_material(resolve_material_state_point(c,&q,&requirements,Pick::SingleClaimOnly).unwrap()).unwrap()
}
fn claim(c:&MaterialCard,key:&str,dims:Dims)->fs_matdb::ClaimId {
    let requirement=ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::StrictlyPositive).unwrap();
    resolve_material_state_point(c,&QueryPoint::new().with("T",300.0).unwrap().with("specimen",7.0).unwrap(),
        &[requirement],Pick::SingleClaimOnly).unwrap().property(key).unwrap().answer().receipt.selected
}
#[test]
fn integrates_selected_cp_and_retains_complete_sheet_and_sample_context() {
    let c=thermal_card(InterpolationPolicy::LinearInside);let original=bundle(&c,7.0);
    let caloric=SheetCaloric::try_new(original.clone(),&c,claim(&c,CP,SPECIFIC_HEAT_CAPACITY_DIMS),
        [290.0,310.0],-30.0,0.01).unwrap();
    assert_eq!(caloric.sheet(),&original);
    assert_eq!(caloric.heat_capacity_samples().len(),3);
    for s in caloric.heat_capacity_samples() {
        assert_eq!(s.card_identity(),c.content_hash());assert_eq!(s.query_point().iter().find(|(name,_)| name=="specimen").map(|(_,v)| *v),Some(7.0));
    }
    // Exact source integral: 10*(1000+1200)/2 + 10*(1200+1300)/2.
    let a=caloric.state_at_temperature(290.0).unwrap();let b=caloric.state_at_temperature(310.0).unwrap();
    close(b.specific_enthalpy_j_kg()-a.specific_enthalpy_j_kg(),23500.0);
    assert_eq!(b.liquid_mass_fraction(),0.0);
    assert!(caloric.state_at_temperature(310.1).is_err());
    let k=conductivity_at_temperature(&original,&c,claim(&c,K,THERMAL_CONDUCTIVITY_DIMS),305.0).unwrap();
    close(k.property(K).unwrap().value_si(),0.07);
    assert_eq!(k.query_point().iter().find(|(name,_)| name=="specimen").map(|(_,v)| *v),Some(7.0));
}
#[test]
fn exact_knots_do_not_authorize_transient_interpolation_or_other_specimens() {
    let exact=thermal_card(InterpolationPolicy::TabulatedOnly);let cp=claim(&exact,CP,SPECIFIC_HEAT_CAPACITY_DIMS);
    assert!(SheetCaloric::try_new(bundle(&exact,7.0),&exact,cp,[290.0,310.0],0.0,0.01).is_err());
    let c=thermal_card(InterpolationPolicy::LinearInside);let cp=claim(&c,CP,SPECIFIC_HEAT_CAPACITY_DIMS);
    assert!(SheetCaloric::try_new(bundle(&c,8.0),&c,cp,[290.0,310.0],0.0,0.01).is_err());
    assert!(SheetCaloric::try_new(bundle(&c,7.0),&c,cp,[289.0,310.0],0.0,0.01).is_err());
    assert!(SheetCaloric::try_new(bundle(&c,7.0),&exact,cp,[290.0,310.0],0.0,0.01).is_err());
}
#[test]
fn conductivity_is_not_a_substitute_for_caloric_storage() {
    let c=thermal_card(InterpolationPolicy::LinearInside);let k=claim(&c,K,THERMAL_CONDUCTIVITY_DIMS);
    let sheet=bundle(&c,7.0);
    assert!(conductivity_at_temperature(&sheet,&c,k,300.0).is_ok());
    assert!(SheetCaloric::try_new(sheet.clone(),&c,k,[290.0,310.0],0.0,0.01).is_err());
    assert!(conductivity_at_temperature(&sheet,&c,k,320.0).is_err());
}

fn conditioned_thermal_card(quantity:Option<QuantitySpec>,range:[f64;2],
    interpolation:InterpolationPolicy,gap:bool)->MaterialCard {
    let mut claims=card(&entries()[..10]).claims().clone();
    let (temperatures,capacities,triangles)=if gap {
        (vec![290.0,295.0,305.0,310.0],vec![1000.0,1100.0,1250.0,1300.0],
            vec![[0,1,4],[1,5,4],[2,3,6],[3,7,6]])
    } else {
        (vec![290.0,300.0,310.0],vec![1000.0,1200.0,1300.0],
            vec![[0,1,3],[1,4,3],[1,2,4],[2,5,4]])
    };
    for (name,dims) in [(CP,SPECIFIC_HEAT_CAPACITY_DIMS),(K,THERMAL_CONDUCTIVITY_DIMS)] {
        let mut samples=Vec::new();
        for (row,&condition) in range.iter().enumerate() {
            for (i,&temperature) in temperatures.iter().enumerate() {
                let value=if name==CP {capacities[i]+1000.0*row as f64}
                    else {0.04+0.002*(temperature-290.0)+0.04*row as f64};
                samples.push(([temperature,condition],value));
            }
        }
        let validity=ValidityDomain::unconstrained().with("T",290.0,310.0);
        let validity=match quantity {
            Some(quantity)=>validity.with_quantity("condition",quantity,range[0],range[1]),
            None=>validity.with("condition",range[0],range[1]),
        };
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Surface2D {axes:["T".into(),"condition".into()],
                axis_dims:[fs_qty::Temperature::DIMS,quantity.map_or(Dims::NONE,|q|q.dims())],
                samples,triangles:triangles.clone(),dims},
            interpolation,validity,uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured conditioned-sheet thermal surface".into(),
                license:"test-only; not measured data".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic conditioned sheet".into(),
        phase:"solid".into(),process:"fixed environment".into(),revision:0},claims,vec![]).unwrap()
}
fn conditioned_point(quantity:Option<QuantitySpec>,condition:f64,t:f64)->QueryPoint {
    let point=QueryPoint::new().with("T",t).unwrap();
    match quantity {
        Some(quantity)=>point.with_quantity("condition",quantity,condition).unwrap(),
        None=>point.with("condition",condition).unwrap(),
    }
}
fn conditioned_sheet(c:&MaterialCard,point:&QueryPoint)->ResolvedSheet {
    resolve_sheet(c,point,SheetCapabilities {membrane:true,bending:true},Pick::SingleClaimOnly).unwrap()
}
fn static_quantity(kind:QuantityKind)->QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind,ValueForm::Static))
}

#[test]
fn g1_held_environment_reaches_sheet_enthalpy_and_conductivity_through_original_surface() {
    for (quantity,range) in [
        (static_quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),[0.0,0.2]),
        (static_quantity(QuantityKind::Pressure),[1e5,2e5]),
    ] {
        let c=conditioned_thermal_card(Some(quantity),range,InterpolationPolicy::LinearInside,false);
        let cp=c.claims_for(CP)[0].0;
        let k=c.claims_for(K)[0].0;
        for fraction in [0.0,0.5,1.0] {
            let condition=range[0]+fraction*(range[1]-range[0]);
            let point=conditioned_point(Some(quantity),condition,300.0);
            let sheet=conditioned_sheet(&c,&point);
            let caloric=SheetCaloric::try_new(sheet.clone(),&c,cp,[290.0,310.0],-30.0,0.01).unwrap();
            assert_eq!(caloric.sheet(),&sheet);
            for sample in caloric.heat_capacity_samples() {
                assert_eq!(sample.card_identity(),c.content_hash());
                assert_eq!(sample.query().axes()["condition"],condition);
                assert_eq!(sample.axis_quantities()["condition"],quantity);
                assert_eq!(sample.property(CP).unwrap().answer().receipt.selected,cp);
            }
            if fraction==0.5 {
                let temperatures=caloric.heat_capacity_samples().iter()
                    .map(|s|s.query().axes()["T"]).collect::<Vec<_>>();
                assert_eq!(temperatures,vec![290.0,295.0,300.0,305.0,310.0]);
            }
            let a=caloric.state_at_temperature(290.0).unwrap();
            let b=caloric.state_at_temperature(310.0).unwrap();
            // The original Cp integral is 23500 J/kg plus a held offset of
            // 1000*fraction J/(kg K) over the same twenty-kelvin interval.
            close(b.specific_enthalpy_j_kg()-a.specific_enthalpy_j_kg(),23500.0+20000.0*fraction);
            assert_eq!(b.liquid_mass_fraction(),0.0);
            let conductivity=conductivity_at_temperature(&sheet,&c,k,305.0).unwrap();
            close(conductivity.property(K).unwrap().value_si(),0.07+0.04*fraction);
            assert_eq!(conductivity.query(),conditioned_point(Some(quantity),condition,305.0));
            assert_eq!(conductivity.property(K).unwrap().answer().receipt.selected,k);
        }
    }
}

#[test]
fn g0_sheet_thermal_surfaces_require_continuous_pinned_support_and_owned_conditions() {
    let moisture=static_quantity(QuantityKind::Moisture(MoistureBasis::DryMass));
    for (interpolation,gap) in [(InterpolationPolicy::TabulatedOnly,false),
        (InterpolationPolicy::LinearInside,true)] {
        let c=conditioned_thermal_card(Some(moisture),[0.0,0.2],interpolation,gap);
        let point=conditioned_point(Some(moisture),0.0,300.0);
        let sheet=conditioned_sheet(&c,&point);
        let cp=c.claims_for(CP)[0].0;
        // Both endpoints are genuine supported samples; they do not license
        // interpolation through exact-only data or a hole between source cells.
        for t in [290.0,310.0] {
            resolve_material_state_point(&c,&conditioned_point(Some(moisture),0.0,t),
                &[ScalarPropertyRequirement::try_new(CP,SPECIFIC_HEAT_CAPACITY_DIMS,
                    ScalarAdmissibility::StrictlyPositive).unwrap()],
                Pick::PinnedByProperty(vec![(CP.into(),cp)])).unwrap();
        }
        assert!(matches!(SheetCaloric::try_new(sheet,&c,cp,[290.0,310.0],0.0,0.01),
            Err(SheetThermalError::Support(_))));
    }
    for quantity in [None,Some(static_quantity(QuantityKind::Stress)),
        Some(static_quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)))] {
        let c=conditioned_thermal_card(quantity,[0.0,1.0],InterpolationPolicy::LinearInside,false);
        let sheet=conditioned_sheet(&c,&conditioned_point(quantity,0.5,300.0));
        assert!(SheetCaloric::try_new(sheet.clone(),&c,c.claims_for(CP)[0].0,
            [290.0,310.0],0.0,0.01).is_err());
        assert!(conductivity_at_temperature(&sheet,&c,c.claims_for(K)[0].0,305.0).is_err());
    }
    let c=conditioned_thermal_card(Some(moisture),[0.0,0.2],InterpolationPolicy::LinearInside,false);
    let wrong=conditioned_sheet(&c,&conditioned_point(Some(static_quantity(
        QuantityKind::Moisture(MoistureBasis::WetMass))),0.1,300.0));
    let missing=conditioned_sheet(&c,&QueryPoint::new().with("T",300.0).unwrap());
    for sheet in [wrong,missing] {
        assert!(SheetCaloric::try_new(sheet,&c,c.claims_for(CP)[0].0,[290.0,310.0],0.0,0.01).is_err());
    }
    // An independently named temperature law must not become constant merely
    // because the consumer evolves T while that other coordinate is held.
    let absolute=static_quantity(QuantityKind::AbsoluteTemperature);
    let mut claims=card(&entries()[..10]).claims().clone();
    for (name,dims,value) in [(CP,SPECIFIC_HEAT_CAPACITY_DIMS,1000.0),
        (K,THERMAL_CONDUCTIVITY_DIMS,0.04)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"source_temperature".into(),
                abscissa_dims:fs_qty::Temperature::DIMS,
                knots:vec![(290.0,value),(310.0,2.0*value)],dims},
            interpolation:InterpolationPolicy::LinearInside,
            validity:ValidityDomain::unconstrained().with("T",290.0,310.0)
                .with_quantity("source_temperature",absolute,290.0,310.0),
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured foreign-temperature law".into(),
                license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    let c=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),
        phase:"solid".into(),process:"foreign temperature".into(),revision:0},claims,vec![]).unwrap();
    let point=QueryPoint::new().with("T",300.0).unwrap()
        .with_quantity("source_temperature",absolute,300.0).unwrap();
    let sheet=conditioned_sheet(&c,&point);
    assert!(SheetCaloric::try_new(sheet.clone(),&c,c.claims_for(CP)[0].0,
        [290.0,310.0],0.0,0.01).is_err());
    assert!(conductivity_at_temperature(&sheet,&c,c.claims_for(K)[0].0,305.0).is_err());
}
