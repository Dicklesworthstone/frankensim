//! G1/G3/G4: authored elastic source paths, not measured material calibration.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::elastic::thermal::AffineIsotropicThermoelastic;
use fs_material::state_point::{MaterialPropertySelection,ResolvedMaterialStatePoint,
    ScalarAdmissibility,ScalarPropertyRequirement,resolve_material_state_point,
    YOUNG_MODULUS_PROPERTY,POISSON_RATIO_PROPERTY,LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY};
use fs_qty::{Dims,Pressure,Temperature,QuantitySpec,semantic::{QuantityKind,SemanticType,ValueForm}};

fn quantity()->QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static))
}
fn point(t:f64)->QueryPoint {QueryPoint::new().with_quantity("T",quantity(),t).unwrap()}
fn source_card(slope:f64,alpha:f64,nu_slope:f64,alpha_slope:f64,nonaffine:bool)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,dims,value,derivative) in [
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,1000.0,slope),
        (POISSON_RATIO_PROPERTY,Dims([0;6]),0.2,nu_slope),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,Dims([0,0,0,-1,0,0]),alpha,alpha_slope),
    ] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*derivative),
                    (300.0,value+if nonaffine && name==YOUNG_MODULUS_PROPERTY {1.0}else{0.0}),
                    (360.0,value+60.0*derivative)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",quantity(),280.0,360.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic affine elastic fixture".into(),
                license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"authored thermoelastic solid".into(),
        phase:"solid".into(),process:"numerical source".into(),revision:0},claims,vec![]).unwrap()
}
fn sources(card:&MaterialCard,point:&QueryPoint)->(ResolvedMaterialStatePoint,ResolvedMaterialStatePoint) {
    let requirements=[
        ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY,Pressure::DIMS,ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new(POISSON_RATIO_PROPERTY,Dims([0;6]),ScalarAdmissibility::Finite).unwrap(),
    ];
    let alpha=ScalarPropertyRequirement::try_new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        Dims([0,0,0,-1,0,0]),ScalarAdmissibility::Finite).unwrap();
    (resolve_material_state_point(card,point,&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap(),
        resolve_material_state_point(card,point,&[alpha],MaterialPropertySelection::SingleClaimOnly).unwrap())
}
fn law(card:&MaterialCard)->AffineIsotropicThermoelastic {
    let (elastic,alpha)=sources(card,&point(300.0));
    AffineIsotropicThermoelastic::new(card,&elastic,&alpha,[280.0,360.0],0.05).unwrap()
}
fn close(a:f64,b:f64,tolerance:f64) {assert!((a-b).abs()<=tolerance,"{a} != {b}");}

// Independent isotropic invariants; shear coordinates here are Mandel.
fn expected(t:f64,strain:[f64;6],slope:f64,alpha:f64,volume:f64)->[f64;4] {
    let youngs=1000.0+slope*(t-300.0); let free=alpha*(t-300.0);
    let mut e=strain; for v in &mut e[..3] {*v-=free;}
    let trace=e[0]+e[1]+e[2];
    let dot=e.iter().map(|v|v*v).sum::<f64>();
    let lambda=youngs*0.2/(1.2*0.6); let mu=youngs/(2.0*1.2);
    let psi=volume*(0.5*lambda*trace*trace+mu*dot);
    let stress_trace=youngs*trace/0.6;
    let entropy=volume*alpha*stress_trace-(slope/youngs)*psi;
    [psi,entropy,psi+t*entropy,
        t*volume*(2.0*slope*alpha*trace/0.6-3.0*youngs*alpha*alpha/0.6)]
}

#[test]
fn g1_thermoelastic_potentials_differentiate_at_fixed_total_mandel_strain() {
    let strain=[0.013,-0.004,0.002,0.007,-0.003,0.005];
    for slope in [-2.0,0.0,2.0] {
        for alpha in [-0.0002,0.0,0.0001] {
            let card=source_card(slope,alpha,0.0,0.0,false); let law=law(&card);
            let t=321.0; let volume=0.7;
            let at=law.at_temperature(t).unwrap(); let actual=at.state(strain,volume).unwrap();
            let values=[actual.free_energy_j,actual.entropy_j_k,actual.internal_energy_j,actual.heat_capacity_j_k];
            for (a,b) in values.into_iter().zip(expected(t,strain,slope,alpha,volume)) {close(a,b,2e-13);}
            let shear=law.at_temperature(300.0).unwrap().state([0.0,0.0,0.0,0.01,0.0,0.0],volume).unwrap();
            let mu=1000.0/(2.0*1.2);
            close(shear.stress_mandel_pa[3],2.0*mu*0.01,2e-14);
            close(shear.free_energy_j,volume*mu*0.01*0.01,2e-14);
            let low=law.at_temperature(t-0.001).unwrap().state(strain,volume).unwrap();
            let high=law.at_temperature(t+0.001).unwrap().state(strain,volume).unwrap();
            close(-(high.free_energy_j-low.free_energy_j)/0.002,actual.entropy_j_k,2e-10);
            close((high.internal_energy_j-low.internal_energy_j)/0.002,actual.heat_capacity_j_k,2e-10);
            let free=at.free_strain().unwrap();
            let relaxed=at.state([free,free,free,0.0,0.0,0.0],volume).unwrap();
            assert_eq!(relaxed.stress_mandel_pa,[0.0;6]);
            assert_eq!(relaxed.free_energy_j,0.0);
            assert_eq!(relaxed.entropy_j_k,0.0);
            assert_eq!(relaxed.internal_energy_j,0.0);
            if slope==0.0 && alpha==0.0 {
                assert_eq!(actual.entropy_j_k,0.0);
                assert_eq!(actual.heat_capacity_j_k,0.0);
                assert_eq!(actual.free_energy_j,actual.internal_energy_j);
            }
            let bound=law.minimum_heat_capacity_j_k(volume).unwrap();
            for t in [280.0,300.0,360.0] {for sign in [-1.0,1.0] {
                let extreme=[0.02*sign;6];
                assert!(law.at_temperature(t).unwrap().state(extreme,volume).unwrap().heat_capacity_j_k>=bound);
            }}
            let (actual_elastic,actual_alpha)=sources(&card,&point(t));
            law.at_materials(&actual_elastic,&actual_alpha).unwrap();
            assert_eq!(at.elastic.card_identity(),card.content_hash());
        }
    }
}

#[test]
fn g4_original_claims_and_supported_thermal_potential_remain_mandatory() {
    let card=source_card(2.0,0.0001,0.0,0.0,false); let law=law(&card);
    assert!(law.at_temperature(360.1).is_err());
    assert!(law.at_temperature(320.0).unwrap().state([0.051;6],1.0).is_err());
    for invalid in [
        source_card(2.0,0.0001,0.01,0.0,false),
        source_card(2.0,0.0001,-0.03,0.0,false),
        source_card(2.0,0.0001,0.0,0.0,true),
    ] {
        let (elastic,alpha)=sources(&invalid,&point(300.0));
        assert!(AffineIsotropicThermoelastic::new(&invalid,&elastic,&alpha,[280.0,360.0],0.05).is_err());
    }
    let other=source_card(3.0,0.0001,0.0,0.0,false);
    let (elastic,alpha)=sources(&other,&point(300.0));
    assert!(law.at_materials(&elastic,&alpha).is_err());

    // A pressure-valued semantic claim has the same SI dimensions as E,
    // but cannot replace the native static elastic coefficient schema.
    let semantic_key=PropertyKey::with_quantity(YOUNG_MODULUS_PROPERTY,
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure,ValueForm::Static)));
    let mut claims=ClaimSet::new();
    for (key,value) in [(semantic_key.clone(),1000.0),
        (PropertyKey::new(POISSON_RATIO_PROPERTY,Dims([0;6])),0.2),
        (PropertyKey::new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,Dims([0,0,0,-1,0,0])),0.0001)] {
        let dims=key.quantity().dims();
        claims.insert_claim(PropertyClaim {key,value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",quantity(),280.0,360.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"wrong semantic kind fixture".into(),
                license:"test-only".into(),artifact:None}}).unwrap();
    }
    let semantic_card=MaterialCard::assemble(MaterialStateId {chemistry:"semantic fixture".into(),
        phase:"solid".into(),process:"numerical".into(),revision:0},claims,vec![]).unwrap();
    let elastic=resolve_material_state_point(&semantic_card,&point(300.0),&[
        ScalarPropertyRequirement::try_with_key(&semantic_key,ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new(POISSON_RATIO_PROPERTY,Dims([0;6]),ScalarAdmissibility::Finite).unwrap()],
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let expansion=resolve_material_state_point(&semantic_card,&point(300.0),&[
        ScalarPropertyRequirement::try_new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0,0,0,-1,0,0]),ScalarAdmissibility::Finite).unwrap()],
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(AffineIsotropicThermoelastic::new(&semantic_card,&elastic,&expansion,[280.0,360.0],0.05).is_err());
}


// Independent tensor-invariant energy, evaluated directly from authored source
// functions. Thermal derivatives below use a five-point difference of THIS
// energy, not the runtime derivative formula or resolved material helpers.
fn varying_energy(t:f64,strain:[f64;6],young_slope:f64,nu_slope:f64,
    alpha:f64,alpha_slope:f64,volume:f64)->f64 {
    let dt=t-300.0;
    let young=1000.0+young_slope*dt;
    let nu=0.2+nu_slope*dt;
    let free=alpha*dt+0.5*alpha_slope*dt*dt;
    let mut e=strain;
    for v in &mut e[..3] {*v-=free;}
    let trace=e[..3].iter().sum::<f64>();
    let lambda=young*nu/((1.0+nu)*(1.0-2.0*nu));
    let mu=young/(2.0*(1.0+nu));
    volume*(0.5*lambda*trace*trace+mu*e.iter().map(|v|v*v).sum::<f64>())
}

#[test]
fn g1_varying_poisson_and_instantaneous_expansion_share_native_stress_and_thermal_derivatives() {
    let strain=[0.013,-0.004,0.002,0.007,-0.003,0.005];
    let volume=0.7;
    for (young_slope,nu_slope,alpha,alpha_slope) in [
        (-2.0,-0.0007,0.0001,-0.000002),
        (2.0,0.0008,-0.0001,0.000003),
        (0.0,0.001,0.0,0.0),
        (0.0,0.0,0.0,0.000002),
    ] {
        let card=source_card(young_slope,alpha,nu_slope,alpha_slope,false);
        let law=law(&card);
        for t in [281.0,300.0,321.0,359.0] {
            let at=law.at_temperature(t).unwrap();
            let actual=at.state(strain,volume).unwrap();
            let dt=t-300.0;
            let free=alpha*dt+0.5*alpha_slope*dt*dt;
            close(at.free_strain().unwrap(),free,2e-17);
            let energy=|temperature|varying_energy(temperature,strain,young_slope,nu_slope,
                alpha,alpha_slope,volume);
            let h=0.01;
            let samples=[energy(t-2.0*h),energy(t-h),energy(t),energy(t+h),energy(t+2.0*h)];
            let entropy=-(samples[0]-8.0*samples[1]+8.0*samples[3]-samples[4])/(12.0*h);
            let capacity=-t*(-samples[0]+16.0*samples[1]-30.0*samples[2]
                +16.0*samples[3]-samples[4])/(12.0*h*h);
            close(actual.free_energy_j,samples[2],2e-13);
            close(actual.entropy_j_k,entropy,2e-11);
            close(actual.internal_energy_j,samples[2]+t*entropy,1e-8);
            close(actual.heat_capacity_j_k,capacity,2e-6);
            let low=law.at_temperature(t-0.001).unwrap().state(strain,volume).unwrap();
            let high=law.at_temperature(t+0.001).unwrap().state(strain,volume).unwrap();
            close((high.internal_energy_j-low.internal_energy_j)/0.002,
                actual.heat_capacity_j_k,2e-9);
            let young=1000.0+young_slope*dt;
            let nu=0.2+nu_slope*dt;
            let mu=young/(2.0*(1.0+nu));
            let lambda=young*nu/((1.0+nu)*(1.0-2.0*nu));
            let trace=strain[..3].iter().sum::<f64>()-3.0*free;
            for k in 0..6 {
                let expected=if k<3 {2.0*mu*(strain[k]-free)+lambda*trace}
                    else {2.0*mu*strain[k]};
                close(actual.stress_mandel_pa[k],expected,2e-13);
            }
            let (elastic,expansion)=sources(&card,&point(t));
            assert_eq!(law.at_materials(&elastic,&expansion).unwrap(),at);
            let relaxed=at.state([free,free,free,0.0,0.0,0.0],volume).unwrap();
            close(relaxed.free_energy_j,0.0,1e-27);
            close(relaxed.entropy_j_k,0.0,1e-16);
        }
    }
}

#[test]
fn g3_capacity_floor_includes_poisson_curvature_and_all_mandel_shears() {
    for (young_slope,nu_slope,alpha,alpha_slope) in [
        (2.0,0.001,0.0001,0.000002),
        (-2.0,-0.001,-0.0001,-0.000002),
        (0.0,0.001,0.0,0.0),
        (0.0,0.0,0.0,0.000002),
    ] {
        let card=source_card(young_slope,alpha,nu_slope,alpha_slope,false);
        let law=law(&card);
        let bound=law.minimum_heat_capacity_j_k(0.7).unwrap();
        assert!(bound<0.0);
        let mut saw_negative=false;
        for step in 0..=40 {
            let at=law.at_temperature(280.0+2.0*f64::from(step)).unwrap();
            for corner in 0..64 {
                let strain=core::array::from_fn(|k|if corner&(1<<k)==0 {-0.03}else{0.03});
                let actual=at.state(strain,0.7).unwrap();
                assert!(actual.heat_capacity_j_k>=bound,
                    "capacity {} below full-envelope bound {bound}",actual.heat_capacity_j_k);
                saw_negative|=actual.heat_capacity_j_k<0.0;
            }
        }
        assert!(saw_negative);
    }
    let card=source_card(0.0,0.0,0.0,0.0,false);
    assert_eq!(law(&card).minimum_heat_capacity_j_k(0.7).unwrap(),0.0);
}

#[path = "thermoelastic_source/acoustic.rs"]
mod acoustic;
