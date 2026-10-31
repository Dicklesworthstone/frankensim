//! Independent synthetic constitutive fixtures; no measured-material claim.
use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use crate::state_point::{ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Temperature, semantic::{QuantityKind, SemanticType, ValueForm}};

const E: [f64;3] = [1200.0,900.0,650.0];
const ET: [f64;3] = [-0.7,0.4,-0.2];
const NU: [f64;3] = [0.16,0.12,0.2];
const NUT: [f64;3] = [0.0006,-0.0004,0.0005];
const G: [f64;3] = [380.0,260.0,310.0];
const GT: [f64;3] = [-0.12,0.08,0.05];
const ALPHA: [f64;3] = [0.0001,-0.00005,0.0002];
const ALPHAT: [f64;3] = [0.000001,-0.0000002,-0.0000008];
const QS: [[f64;3];3] = [[0.8,-0.6,0.0],[0.6,0.8,0.0],[0.0,0.0,1.0]];
const QA: [[f64;3];3] = [[0.8,0.0,0.6],[0.0,1.0,0.0],[-0.6,0.0,0.8]];

fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static))
}
fn point(t:f64) -> QueryPoint { QueryPoint::new().with_quantity("T",quantity(),t).unwrap() }
fn card(nonaffine:bool) -> MaterialCard {
    let mut claims=ClaimSet::new();
    let mut insert=|name:&str,dims:Dims,base:f64,slope:f64| {
        let middle=base+if nonaffine && name==ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0] {1.0} else {0.0};
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,base-20.0*slope),(300.0,middle),(360.0,base+60.0*slope)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",quantity(),280.0,360.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"authored orthotropic thermal fixture".into(),
                license:"test-only".into(),artifact:None}}).unwrap();
    };
    for i in 0..3 {
        insert(ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i],Pressure::DIMS,E[i],ET[i]);
        insert(ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i],Dims::NONE,NU[i],NUT[i]);
        insert(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i],Pressure::DIMS,G[i],GT[i]);
        insert(ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i],Dims([0,0,0,-1,0,0]),ALPHA[i],ALPHAT[i]);
    }
    insert(YOUNG_MODULUS_PROPERTY,Pressure::DIMS,E[0],ET[0]);
    insert(POISSON_RATIO_PROPERTY,Dims::NONE,NU[0],NUT[0]);
    insert(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,Dims([0,0,0,-1,0,0]),ALPHA[0],ALPHAT[0]);
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic anisotropic solid".into(),phase:"solid".into(),
        process:"thermal derivative test".into(),revision:0},claims,vec![]).unwrap()
}
fn sources(card:&MaterialCard,t:f64,elastic_isotropic:bool,expansion_isotropic:bool)
    -> (ResolvedMaterialStatePoint,ResolvedMaterialStatePoint) {
    let names:Vec<_>=if elastic_isotropic {vec![(YOUNG_MODULUS_PROPERTY,Pressure::DIMS),(POISSON_RATIO_PROPERTY,Dims::NONE)]}
        else {ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.map(|n|(n,Pressure::DIMS)).into_iter()
            .chain(ORTHOTROPIC_POISSON_RATIO_PROPERTIES.map(|n|(n,Dims::NONE)))
            .chain(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.map(|n|(n,Pressure::DIMS))).collect()};
    let requirements:Vec<_>=names.into_iter().map(|(n,d)|ScalarPropertyRequirement::try_new(n,d,ScalarAdmissibility::Finite).unwrap()).collect();
    let expansion_names=if expansion_isotropic {vec![LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY]}
        else {ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES.to_vec()};
    let expansion:Vec<_>=expansion_names.into_iter().map(|n|ScalarPropertyRequirement::try_new(n,Dims([0,0,0,-1,0,0]),ScalarAdmissibility::Finite).unwrap()).collect();
    (resolve_material_state_point(card,&point(t),&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap(),
        resolve_material_state_point(card,&point(t),&expansion,MaterialPropertySelection::SingleClaimOnly).unwrap())
}
fn law(card:&MaterialCard,qs:[[f64;3];3],qa:[[f64;3];3]) -> AffineThermoelastic {
    let (e,a)=sources(card,300.0,false,false);
    AffineThermoelastic::new(card,&e,&a,
        ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:qs},
        ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:qa},[280.0,360.0],0.1).unwrap()
}
fn close(a:f64,b:f64,tolerance:f64) { assert!((a-b).abs()<=tolerance,"{a:.16e} != {b:.16e}"); }
fn tensor(v:[f64;6]) -> [[f64;3];3] {
    let s=core::f64::consts::FRAC_1_SQRT_2;
    [[v[0],s*v[3],s*v[5]],[s*v[3],v[1],s*v[4]],[s*v[5],s*v[4],v[2]]]
}
fn mandel(a:[[f64;3];3]) -> [f64;6] {
    let s=core::f64::consts::SQRT_2;[a[0][0],a[1][1],a[2][2],s*a[0][1],s*a[1][2],s*a[2][0]]
}
fn rotate(a:[[f64;3];3],q:[[f64;3];3]) -> [[f64;3];3] {
    core::array::from_fn(|i|core::array::from_fn(|j|
        (0..3).map(|p|(0..3).map(|r|q[i][p]*a[p][r]*q[j][r]).sum::<f64>()).sum()))
}
fn transpose(q:[[f64;3];3]) -> [[f64;3];3] {core::array::from_fn(|i|core::array::from_fn(|j|q[j][i]))}
fn product(a:[[f64;3];3],b:[[f64;3];3]) -> [[f64;3];3] {
    core::array::from_fn(|i|core::array::from_fn(|j|(0..3).map(|k|a[i][k]*b[k][j]).sum()))
}

// Independent physical-tensor oracle: rotate strain into principal axes,
// solve the 3x3 compliance by elimination (not the production cofactors),
// and contract physical stress/strain with the tensor shear multiplicity.
fn oracle(t:f64,total:[f64;6],qs:[[f64;3];3],qa:[[f64;3];3],volume:f64) -> (f64,[f64;6]) {
    let dt=t-300.0;
    let e=core::array::from_fn::<_,3,_>(|i|E[i]+ET[i]*dt);
    let nu=core::array::from_fn::<_,3,_>(|i|NU[i]+NUT[i]*dt);
    let g=core::array::from_fn::<_,3,_>(|i|G[i]+GT[i]*dt);
    let mut thermal=[[0.0;3];3];
    for i in 0..3 {thermal[i][i]=ALPHA[i]*dt+0.5*ALPHAT[i]*dt*dt;}
    let thermal=rotate(thermal,qa); let total=tensor(total);
    let elastic=core::array::from_fn(|i|core::array::from_fn(|j|total[i][j]-thermal[i][j]));
    let principal=rotate(elastic,transpose(qs));
    let mut system=[
        [1.0/e[0],-nu[0]/e[0],-nu[1]/e[0],principal[0][0]],
        [-nu[0]/e[0],1.0/e[1],-nu[2]/e[1],principal[1][1]],
        [-nu[1]/e[0],-nu[2]/e[1],1.0/e[2],principal[2][2]],
    ];
    for k in 0..3 {
        let pivot=system[k][k];for j in k..4 {system[k][j]/=pivot;}
        for i in 0..3 {if i!=k {let factor=system[i][k];for j in k..4 {system[i][j]-=factor*system[k][j];}}}
    }
    let s12=2.0*g[0]*principal[0][1];let s23=2.0*g[1]*principal[1][2];let s31=2.0*g[2]*principal[2][0];
    let stress=[[system[0][3],s12,s31],[s12,system[1][3],s23],[s31,s23,system[2][3]]];
    let energy=0.5*volume*(0..3).map(|i|(0..3).map(|j|principal[i][j]*stress[i][j]).sum::<f64>()).sum::<f64>();
    (energy,mandel(rotate(stress,qs)))
}

#[test]
fn anisotropic_source_potential_differentiates_original_compliance_and_integrated_expansion() {
    let card=card(false);let law=law(&card,QS,QA);
    let strain=[0.013,-0.004,0.002,0.007,-0.003,0.005];let volume=0.7;
    for t in [281.0,300.0,321.0,359.0] {
        let at=law.at_temperature(t).unwrap();let state=at.state(strain,volume).unwrap();
        let (psi,stress)=oracle(t,strain,QS,QA,volume);
        close(state.free_energy_j,psi,2e-13);
        for (a,b) in state.stress_mandel_pa.into_iter().zip(stress) {close(a,b,2e-12);}
        let h=0.02;
        let samples=[-2.0,-1.0,0.0,1.0,2.0].map(|k|oracle(t+k*h,strain,QS,QA,volume).0);
        let entropy=-(samples[0]-8.0*samples[1]+8.0*samples[3]-samples[4])/(12.0*h);
        let capacity=-t*(-samples[0]+16.0*samples[1]-30.0*samples[2]+16.0*samples[3]-samples[4])/(12.0*h*h);
        close(state.entropy_j_k,entropy,2e-10);close(state.heat_capacity_j_k,capacity,2e-6);
        close(state.internal_energy_j,psi+t*entropy,1e-7);
        let low=oracle(t-0.001,strain,QS,QA,volume).1;let high=oracle(t+0.001,strain,QS,QA,volume).1;
        for k in 0..6 {close(state.stress_temperature_derivative_mandel_pa_k[k],(high[k]-low[k])/0.002,2e-9);}
        let relaxed=at.state(at.free_strain_mandel(),volume).unwrap();
        assert_eq!(relaxed.stress_mandel_pa,[0.0;6]);assert_eq!(relaxed.free_energy_j,0.0);
        assert_eq!(relaxed.entropy_j_k,0.0);assert!(at.free_strain().is_err());
        let (elastic,expansion)=sources(&card,t,false,false);
        assert_eq!(law.at_materials(&elastic,&expansion).unwrap(),at);
    }
}

#[test]
fn rotating_geometry_and_both_independent_material_frames_preserves_thermal_work_and_acoustics() {
    let card=card(false);let first=law(&card,QS,QA);let q=transpose(QS);
    let rotated=law(&card,product(q,QS),product(q,QA));
    let strain=[0.013,-0.004,0.002,0.007,-0.003,0.005];
    let transformed=mandel(rotate(tensor(strain),q));
    let a=first.at_temperature(321.0).unwrap();let b=rotated.at_temperature(321.0).unwrap();
    let sa=a.state(strain,0.7).unwrap();let sb=b.state(transformed,0.7).unwrap();
    for (a,b) in [sa.free_energy_j,sa.entropy_j_k,sa.internal_energy_j,sa.heat_capacity_j_k].into_iter()
        .zip([sb.free_energy_j,sb.entropy_j_k,sb.internal_energy_j,sb.heat_capacity_j_k]) {close(a,b,2e-12);}
    for (actual,expected) in sb.stress_mandel_pa.into_iter().zip(mandel(rotate(tensor(sa.stress_mandel_pa),q))) {close(actual,expected,5e-12);}
    let input=|total|super::super::acoustic::AdiabaticThermoelasticInput {total_strain_mandel:total,
        reference_density_kg_m3:2.0,caloric_heat_capacity_j_m3_k:1000.0,
        mechanical_material_identity:ContentHash([1;32]),caloric_state_identity:ContentHash([2;32])};
    let aa=a.adiabatic_acoustics(input(strain)).unwrap();let bb=b.adiabatic_acoustics(input(transformed)).unwrap();
    let expected=rotate_mandel(*aa.acoustic().stiffness_mandel_pa(),q);
    for (actual,expected) in bb.acoustic().stiffness_mandel_pa().iter().flatten().zip(expected.into_iter().flatten()) {close(*actual,expected,2e-11);}
    assert_eq!(aa.isothermal_stiffness_mandel_pa(),a.stiffness_mandel_pa());
}

#[test]
fn full_interval_capacity_bound_and_source_refusals_do_not_depend_on_endpoint_samples() {
    let card=card(false);let law=law(&card,QS,QA);let bound=law.minimum_heat_capacity_j_k(0.7).unwrap();
    let mut seed=41u64;
    for i in 0..129 {
        let t=280.0+80.0*(i as f64/128.0);let at=law.at_temperature(t).unwrap();
        let strain=core::array::from_fn(|_| {seed=seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            0.08*((seed>>11) as f64/(1u64<<53) as f64-0.5)});
        assert!(at.state(strain,0.7).unwrap().heat_capacity_j_k>=bound);
    }
    assert!(law.at_temperature(360.01).is_err());
    assert!(law.at_temperature(321.0).unwrap().state([0.2;6],1.0).is_err());
    let bad=self::card(true);let (e,a)=sources(&bad,300.0,false,false);
    assert!(AffineThermoelastic::new(&bad,&e,&a,ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:QS},
        ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:QA},[280.0,360.0],0.1).is_err());
    assert!(law.at_materials(&e,&a).is_err());
    let (e,a)=sources(&card,300.0,false,false);
    let mut reflection=IDENTITY;reflection[0][0]=-1.0;
    assert!(AffineThermoelastic::new(&card,&e,&a,ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:reflection},
        ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:QA},[280.0,360.0],0.1).is_err());
}

#[test]
fn elastic_and_expansion_symmetries_remain_independent_and_isotropic_arithmetic_is_preserved() {
    let card=card(false);
    for elastic_isotropic in [false,true] {for expansion_isotropic in [false,true] {
        let (e,a)=sources(&card,300.0,elastic_isotropic,expansion_isotropic);
        let law=AffineThermoelastic::new(&card,&e,&a,
            if elastic_isotropic {ThermoelasticSymmetry::Isotropic} else {ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:QS}},
            if expansion_isotropic {ThermoelasticSymmetry::Isotropic} else {ThermoelasticSymmetry::OrientedOrthotropic {principal_to_target:QA}},
            [280.0,360.0],0.1).unwrap();
        let at=law.at_temperature(321.0).unwrap();
        assert_eq!(at.state(at.free_strain_mandel(),1.0).unwrap().free_energy_j,0.0);
        if elastic_isotropic && expansion_isotropic {
            let original=AffineIsotropicThermoelastic::new(&card,&e,&a,[280.0,360.0],0.1).unwrap();
            assert_eq!(at,original.at_temperature(321.0).unwrap());
            assert_eq!(law.minimum_heat_capacity_j_k(0.7).unwrap(),original.minimum_heat_capacity_j_k(0.7).unwrap());
        }
    }}
}
