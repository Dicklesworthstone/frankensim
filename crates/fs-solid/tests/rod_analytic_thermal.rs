//! G1/G4: analytical rod derivatives respect narrow source strain domains.
//! Synthetic source data exercise the original material and equilibrium owners.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::material_field::supports::RodSupport;

fn source(alphas:&[f64])->(ConcentricShearSection,ConcentricThermoelasticSection) {
    let reference=QueryPoint::new().with("T",300.0).unwrap();
    let current=QueryPoint::new().with("T",301.0).unwrap();
    let mut layers=Vec::new();let mut paths=Vec::new();
    for (i,&alpha) in alphas.iter().enumerate() {
        let mut claims=ClaimSet::new();
        for (name,value,dims) in [("density",1000.0,Density::DIMS),
            ("young_modulus",1e6,Pressure::DIMS),("shear_modulus",4e5,Pressure::DIMS),
            ("linear_thermal_expansion_coefficient",alpha,INVERSE_TEMPERATURE_DIMS)] {
            claims.insert_claim(PropertyClaim {
                key:PropertyKey::new(name,dims),value:PropertyValue::Scalar{value,dims},
                interpolation:InterpolationPolicy::ConstantWithinValidity,
                validity:ValidityDomain::unconstrained().with("T",300.0,400.0),
                uncertainty:UncertaintyModel::Unstated,observations:vec![],
                provenance:Provenance{source:"authored narrow-strain analytical rod".into(),license:"test-only".into(),artifact:None},
            }).unwrap();
        }
        let card=MaterialCard::assemble(MaterialStateId{chemistry:format!("synthetic layer {i}"),
            phase:"solid".into(),process:"test".into(),revision:0},claims,vec![]).unwrap();
        let requirements=[("density",Density::DIMS),("young_modulus",Pressure::DIMS),
            ("shear_modulus",Pressure::DIMS)].map(|(name,dims)|
                ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive).unwrap());
        let material=resolve_material_state_point(&card,&current,&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap();
        paths.push(integrate_isotropic_thermal_expansion(&card,&reference,&current,
            MaterialPropertySelection::SingleClaimOnly).unwrap());
        layers.push(ConcentricMaterialLayer{outer_radius_m:0.01*(i+1) as f64/alphas.len() as f64,material});
    }
    let elastic=ConcentricElasticSection::try_new(0.0,layers).unwrap();
    let thermal=ConcentricThermoelasticSection::try_new(elastic.clone(),paths).unwrap();
    (elastic.with_shear_modulus("shear_modulus",0.8).unwrap(),thermal)
}

#[test]
fn g1_free_expansion_inside_a_domain_far_narrower_than_old_derivative_probes() {
    let gate=CancelGate::new();
    let (section,thermal)=source(&[1e-8]);
    let mut rod=SpatialMaterialRod::try_straight(0.1,vec![section]).unwrap()
        .with_thermal_sections(vec![thermal],1e-7,&gate).unwrap();
    let mass=rod.reference_mass().total_mass_kg();
    rod.solve_static(&TipLoad::default(),1,1e-10,&gate).unwrap();
    assert!((rod.positions_m()[1][0]-0.1*(1.0+1e-8)).abs()<5e-14);
    assert!(rod.segment_resultants(0).unwrap().0[0].abs()<1e-10);
    assert_eq!(rod.reference_mass().total_mass_kg(),mass);
    let before=rod.positions_m().to_vec();
    let energy=rod.energy_j().unwrap();
    assert!(rod.solve_static(&TipLoad{force:[1.0,0.0,0.0],moment:[0.0;3]},1,1e-10,&gate).is_err());
    assert_eq!(rod.positions_m(),before.as_slice());assert_eq!(rod.energy_j().unwrap(),energy);
    gate.request();
    assert!(matches!(rod.solve_static(&TipLoad::default(),1,1e-10,&gate),Err(SpatialRodError::Cancelled)));
    assert_eq!(rod.positions_m(),before.as_slice());
}

#[test]
fn g1_clamped_thermal_reactions_require_no_out_of_domain_probe() {
    let gate=CancelGate::new();
    let (section,thermal)=source(&[1e-8]);
    let ea=section.elastic().axial_stiffness_n();
    let mut rod=SpatialMaterialRod::try_straight(0.1,vec![section]).unwrap()
        .with_thermal_sections(vec![thermal],1.1e-8,&gate).unwrap()
        .with_supports(&[RodSupport{node:1,position_m:[Some(0.1),Some(0.0),Some(0.0)],frame:Some(So3::identity())}],&gate).unwrap();
    assert_eq!(rod.solve_static(&TipLoad::default(),1,1e-10,&gate).unwrap(),vec![0.0]);
    let root=rod.root_support_reaction(&gate).unwrap();
    let tip=rod.support_reactions(&TipLoad::default(),&gate).unwrap();
    assert!((root.force_n[0].unwrap()-ea*1e-8).abs()<1e-15);
    assert!((tip.nodes[0].force_n[0].unwrap()+ea*1e-8).abs()<1e-15);
    assert_eq!(tip.free_residual_norm,0.0);
    assert!(rod.energy_j().unwrap()>0.0);
    assert_eq!(rod.positions_m()[1],[0.1,0.0,0.0]);
}

#[test]
fn g1_opposing_layer_stresses_survive_analytical_thermal_equilibrium() {
    let gate=CancelGate::new();
    let (section,thermal)=source(&[1e-8,2e-8]);
    let mut rod=SpatialMaterialRod::try_straight(0.1,vec![section]).unwrap()
        .with_thermal_sections(vec![thermal],1e-7,&gate).unwrap();
    rod.solve_static(&TipLoad::default(),1,1e-10,&gate).unwrap();
    let response=rod.thermal_axial_response(0).unwrap().unwrap();
    assert!((rod.segment_strains(0).unwrap().0[0]-1.75e-8).abs()<4e-13);
    assert!(response.force_n().abs()<1e-10);
    assert!((response.layers()[0].stress_pa()-0.0075).abs()<4e-7);
    assert!((response.layers()[1].stress_pa()+0.0025).abs()<4e-7);
    assert!(rod.energy_j().unwrap()>0.0);
    assert_eq!(rod.thermal_sections().unwrap()[0].section(),rod.sections()[0].elastic());
}
