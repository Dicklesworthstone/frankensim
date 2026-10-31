//! G1/G3: actual shared heating and acoustic pressure, with data-only expansion.
use fs_couple::material_assembly::{AcousticMaterialBindings,
    concentric::thermal::heating::compile_heated_concentric_material_assembly as compile};
use fs_material::conductor::{OhmicDrive,ELECTRICAL_RESISTIVITY_DIMS,ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::phase::{HeatCapacitySourceQuery,ResolvedHeatCapacityCurve,SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer,thermal::AxialConstraint};
use fs_material::section_temperature::{ConcentricSectionTemperatureLaw,
    heating::{HeatedConcentricSection,SectionThermalInterface}};
use fs_material::state_point::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{MaterialCard,MaterialStateId,ClaimSet,PropertyKey,PropertyClaim,PropertyValue,
    InterpolationPolicy,Provenance,UncertaintyModel,QueryPoint,SelectionPolicy};
use fs_qty::{Density,Pressure,Temperature};
use fs_scenario::{AcousticAssembly,AmbientGas,Listener,Pluck,PrestressedString};
use fs_conduction::transient::deposition::{Boundary,DepositionBudget,StepConfig,CellId};
const LENGTH:f64=0.65;
fn point(t:f64)->QueryPoint {QueryPoint::new().with("T",t).unwrap()}
fn card(alpha:f64)->(MaterialCard,Vec<ScalarPropertyRequirement>) {
    let mut claims=ClaimSet::new();let mut req=Vec::new();
    for (name,dims,lo,hi) in [(DENSITY_PROPERTY,Density::DIMS,8000.0,8000.0),
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,2e11,2e11),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS,500.0,500.0),
        (THERMAL_CONDUCTIVITY_PROPERTY,THERMAL_CONDUCTIVITY_DIMS,20.0,20.0),
        (ELECTRICAL_RESISTIVITY_PROPERTY,ELECTRICAL_RESISTIVITY_DIMS,1e-7,2e-7),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS,alpha,alpha)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.0,lo),(400.0,hi)],dims},
            validity:ValidityDomain::unconstrained().with("T",300.0,400.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic constrained section".into(),license:"test-only".into(),artifact:None},
        }).unwrap();
        req.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    (MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),process:"test".into(),revision:0},claims,vec![]).unwrap(),req)
}
fn caloric(card:&MaterialCard)->ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
        lower:&point(300.0),upper:&point(400.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0,0.001).unwrap()
}
fn body<'a>(card:&'a MaterialCard,req:&[ScalarPropertyRequirement],chart:&'a ResolvedHeatCapacityCurve)->HeatedConcentricSection<'a> {
    let material=resolve_material_state_point(card,&point(300.0),req,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section=ConcentricElasticSection::try_new(0.0,[0.00025,0.0005].map(|outer_radius_m|
        ConcentricMaterialLayer {outer_radius_m,material:material.clone()}).to_vec()).unwrap();
    let law=ConcentricSectionTemperatureLaw::try_new(section,&[card,card],&["T","T"]).unwrap();
    let boundary=Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0};
    HeatedConcentricSection::try_new(law,&[chart,chart],LENGTH,&[boundary;2],vec![SectionThermalInterface {
        area_m2:0.0005,resistance_geometry_per_m:[1.0;2],contact_resistance_k_w:0.0,
    }],DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:10000}).unwrap()
}
fn template()->AcousticAssembly {
    AcousticAssembly {ambient:AmbientGas::sea_level(),string:Some(PrestressedString {
        length_m:LENGTH,tension_n:10.0,lin_density_kg_m:1.0,axial_stiffness_n:0.0,width_m:0.01,
        n_modes:4,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,
        kelvin_voigt_bending:None,relaxation_bending:None,polarization_detune:0.0,moving_end:false,
        tension_loss_factor:None,
        bending_loss_factor:None,
    }),duct:None,pluck:Some(Pluck {station_frac:0.25,height_m:1e-5}),bow:None,blow:None,
        reed:None,soundboard:None,body_modes:vec![],plate:None,cavity:None,obstacles:vec![],
        contact_texture:None,listener:Listener {distance_m:1.0},
        // Midpoint string band: a pole rate near 1100 rad/s needs dt <= 3.7e-5 s.
        sample_rate_hz:32000,duration_s:0.008}
}
fn config()->StepConfig {StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7}}

#[test]
fn expansion_data_alone_changes_fixed_end_tension_and_actual_pressure() {
    let run=|alpha| {
        let (card,req)=card(alpha);let chart=caloric(&card);
        let initial=body(&card,&req,&chart).with_axial_constraint(AxialConstraint::FixedStrain(1e-4),0.01).unwrap();
        let (heated,report)=initial.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
        let template=template();let saved=template.clone();let bindings=AcousticMaterialBindings::default();
        let cold=compile(&template,&initial,&bindings).unwrap().realize().unwrap();
        let compiled=compile(&template,&heated,&bindings).unwrap();
        assert!(std::ptr::eq(compiled.source(),&heated));
        let tension=compiled.compiled().compiled().assembly().string.as_ref().unwrap().tension_n;
        assert_eq!(tension,heated.axial_state().unwrap().response().force_n());
        assert_eq!(initial.section().axial_stiffness_n(),heated.section().axial_stiffness_n());
        let hot=compiled.realize().unwrap();assert_eq!(cold.gas,hot.gas);assert_eq!(template,saved);
        assert_eq!(hot.pressure_pa.len(),256);assert!(hot.pressure_pa.iter().any(|p|*p!=0.0));
        assert!(hot.pressure_pa.iter().all(|p|p.is_finite()));
        (cold.pressure_pa,hot.pressure_pa,tension,heated.thermal().temperature_k(CellId(0)).unwrap(),report.electrical_energy_j)
    };
    let zero=run(0.0);let expanding=run(2e-5);
    assert_eq!(zero.0,zero.1);assert_eq!(zero.0,expanding.0);assert_ne!(expanding.0,expanding.1);
    assert!(expanding.2<zero.2);assert_eq!(zero.3,expanding.3);assert_eq!(zero.4,expanding.4);
}

#[test]
fn stale_geometry_missing_constraints_and_moving_supports_refuse() {
    let (card,req)=card(2e-5);let chart=caloric(&card);let original=body(&card,&req,&chart);
    let bindings=AcousticMaterialBindings::default();let mut template=template();
    assert!(compile(&template,&original,&bindings).is_err());
    let force=original.clone().with_axial_constraint(AxialConstraint::FixedForce(10.0),0.01).unwrap();
    assert!(compile(&template,&force,&bindings).is_err());
    let fixed=original.with_axial_constraint(AxialConstraint::FixedStrain(1e-4),0.01).unwrap();
    template.string.as_mut().unwrap().length_m*=2.0;
    assert!(compile(&template,&fixed,&bindings).is_err());
    template.string.as_mut().unwrap().length_m=LENGTH;template.string.as_mut().unwrap().moving_end=true;
    assert!(compile(&template,&fixed,&bindings).is_err());
    template.string.as_mut().unwrap().moving_end=false;template.string.as_mut().unwrap().damping_ratio=0.01;
    assert!(compile(&template,&fixed,&bindings).is_err());
    assert_eq!(fixed.thermal().time_s(),0.0);
}

#[test]
fn real_thermal_compression_is_not_clamped_into_a_taut_string() {
    let (card,req)=card(2e-5);let chart=caloric(&card);
    let fixed=body(&card,&req,&chart).with_axial_constraint(AxialConstraint::FixedStrain(1e-7),0.01).unwrap();
    let (heated,_)=fixed.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
    assert!(heated.axial_state().unwrap().response().force_n()<0.0);
    assert!(compile(&template(),&heated,&AcousticMaterialBindings::default()).is_err());
    assert!(fixed.axial_state().unwrap().response().force_n()>0.0);
    assert_eq!(fixed.thermal().time_s(),0.0);
}
