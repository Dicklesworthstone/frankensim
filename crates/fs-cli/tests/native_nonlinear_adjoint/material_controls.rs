//! Change actual normalized material packs and re-run the native pipeline.
use super::*;
use fs_project::{ConductionRadiation, RadiatingSurface, ThermalBoundaryCondition, spec::dims};
use fs_qty::QtyAny;

fn card(scale: f64, nonlinear: bool) -> (Vec<u8>,String,String) {
    let source=fs_blake3::hash_bytes(format!("synthetic conductivity scale={scale}, nonlinear={nonlinear}; fixed emissivity=0.8").as_bytes());
    let provenance=|| Provenance {source:"synthetic material-control regression, not experimental data".into(),
        license:"CC-BY-4.0; redistribution permitted".into(),artifact:Some(source)};
    let mut claims=ClaimSet::new();
    let observation=claims.register_observation(ObservationDataset {specimen:"material-control fixture".into(),
        method:"declared conductivity multiplier; fixed emissivity".into(),artifact:source,
        caveats:"numerical verification only".into(),provenance:provenance()}).unwrap();
    let kd=fs_conduction::CONDUCTIVITY_DIMS;
    let conductivity=if nonlinear { PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:dims::TEMPERATURE,
        knots:vec![(250.0,scale),(450.0,11.0*scale)],dims:kd} }
        else {PropertyValue::Scalar {value:6.0*scale,dims:kd}};
    for (name,dims,value,interpolation) in [
        ("thermal-conductivity",kd,conductivity,if nonlinear {InterpolationPolicy::LinearInside}
            else {InterpolationPolicy::ConstantWithinValidity}),
        (fs_conduction::SURFACE_EMISSIVITY_PROPERTY,fs_qty::Dims::NONE,
            PropertyValue::Scalar {value:0.8,dims:fs_qty::Dims::NONE},InterpolationPolicy::ConstantWithinValidity),
    ] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value,
            validity:fs_evidence::ValidityDomain::unconstrained().with("T",250.0,450.0),
            uncertainty:UncertaintyModel::Unstated,interpolation,
            observations:vec![observation],provenance:provenance()}).unwrap();
    }
    let pack=NormalizedMaterialCardPack::new(MaterialStateId {chemistry:"conductivity-control".into(),
        phase:"solid".into(),process:"synthetic".into(),revision:0},
        NormalizedPack::new("conductivity-control","synthetic-v1",source,
            "CC-BY-4.0; redistribution permitted",claims,Vec::new(),Vec::new()).unwrap()).unwrap();
    let bytes=pack.to_bytes();
    let admitted=CardPackSet::admit(vec![RawCardPack {kind:CardPackKind::Material,
        source:"conductivity-control-fixture".into(),bytes:bytes.clone(),expect:None}]).unwrap();
    let selected=&admitted.materials()[0];
    (bytes,selected.card().to_hex(),selected.identity().to_string())
}

fn material(fixture: &Fixture, project: &mut fs_project::ProjectSpec, scale: f64, nonlinear: bool) {
    let (bytes,identity,state)=card(scale,nonlinear);
    // Only this fixture's private input is rewritten; retained child runs keep
    // their own exact original packs and the real material resolver still runs.
    std::fs::write(fixture.dir.join("nonlinear.fsmcdpk"),bytes).unwrap();
    let binding=&mut project.materials.as_mut().unwrap()[0];
    binding.card=identity.clone();binding.state=state;binding.claim=None;
    binding.temp_lo.value=250.0;binding.temp_hi.value=450.0;
    binding.conductivity_tolerance=None;
    if let Some(radiation)=project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().radiation.as_mut() {
        for surface in &mut radiation.surfaces {surface.card=identity.clone();surface.claim=None;}
    }
}

#[test]
fn native_conductivity_multiplier_matches_full_material_and_cooling_resolves() {
    for nonlinear in [false,true] { for air in [false,true] { for radiation in [false,true] {
        let fixture=Fixture::new();
        let mut base=fixture.project.clone();
        base.power.as_mut().unwrap()[0].watts.value=500.0;
        let setup=base.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        if air { setup.boundaries[0].condition=ThermalBoundaryCondition::AirflowConvection {
            branch:"air".into(),order:0,inlet_temperature:QtyAny::new(300.0,dims::TEMPERATURE),
            hydraulic_diameter:QtyAny::new(0.02,dims::LENGTH),flow_area:QtyAny::new(0.004,dims::AREA),
            channel_length:QtyAny::new(0.3,dims::LENGTH),correlation:"convection.gnielinski".into(),
        }; }
        if radiation { setup.radiation=Some(ConductionRadiation {surfaces:vec![RadiatingSurface {
            name:"gray-wall".into(),target:"air".into(),card:String::new(),claim:None,
            query_temperature:QtyAny::new(300.0,dims::TEMPERATURE),
            reservoir_temperature:QtyAny::new(285.0,dims::TEMPERATURE)}],max_iterations:128,
            temperature_tolerance:QtyAny::new(1e-11,dims::TEMPERATURE),
            heat_tolerance:QtyAny::new(1e-8,dims::POWER),relaxation:0.5}); }
        material(&fixture,&mut base,1.0,nonlinear);
        let (plain,plain_field)=fixture.solve(&base,0);
        let mut requested=base.clone();
        requested.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
            name:"temperature-max-adjoint".into(),kind:"report".into(),region:None});
        let (receipt,field)=fixture.solve(&requested,1);
        assert_eq!(plain_field.get("temperature"),field.get("temperature"));
        for key in ["energy","conjugate","radiation"] {assert_eq!(plain.get(key),receipt.get(key),"{key}");}
        let report=receipt.get("nominal_adjoint").unwrap();
        let rows:Vec<_>=report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|r| r.str_field("target")==Some("conductivity-multiplier")).collect();
        assert_eq!(rows.len(),1);
        assert_eq!(rows[0].str_field("entity"),Some("air"));
        assert_eq!(rows[0].str_field("parameter_unit"),Some("1"));
        assert_eq!(rows[0].f64_field("reference_value"),Some(1.0));
        assert_eq!(rows[0].str_field("material_card"),Some(base.materials.as_ref().unwrap()[0].card.as_str()));
        let actual=rows[0].f64_field("derivative").unwrap();
        let vertex=report.f64_field("selected_vertex").unwrap() as usize;
        let epsilon=1e-3;
        let mut values=Vec::new();
        for (side,sign) in [-1.0,1.0].into_iter().enumerate() {
            let mut shifted=base.clone();material(&fixture,&mut shifted,1.0+sign*epsilon,nonlinear);
            let (_,field)=fixture.solve(&shifted,2+side);
            values.push(field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap());
        }
        let expected=(values[1]-values[0])/(2.0*epsilon);
        assert!(expected.abs()>1e-5,"the fixture must resolve material sensitivity");
        assert!((actual-expected).abs()<2e-3*expected.abs().max(1e-4),
            "nonlinear={nonlinear}, air={air}, radiation={radiation}: {actual:e} vs {expected:e}");
    } } }
}
