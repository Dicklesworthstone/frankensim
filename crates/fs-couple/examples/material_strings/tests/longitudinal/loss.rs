//! G1/G3/G5: material-only axial loss reaches the original playable network.
use super::*;
use fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY as AXIAL_ETA;

fn source_pack(name: &str, eta: f64, band: f64) -> NormalizedPack {
    let mut claims = pack(name,7800.0,200_000.0,None).claims().clone();
    let provenance = Provenance { source: "synthetic axial viscosity; not measured wire data".into(),
        license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([23;32])) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: name.into(), method: "manufactured uniaxial loss".into(),
        artifact: fs_blake3::ContentHash([23;32]), caveats: "not a measured material".into(),
        provenance: provenance.clone(),
    }).unwrap();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AXIAL_ETA,DynViscosity::DIMS),
        value: PropertyValue::Scalar { value: eta, dims: DynViscosity::DIMS },
        validity: ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,band),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation], provenance,
    }).unwrap();
    let packed = NormalizedPack::new(name,"synthetic-axial-test-v1",fs_blake3::ContentHash([23;32]),
        "synthetic redistribution permitted",claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&packed.to_bytes()).unwrap()
}
fn selected_text() -> String {
    text().replace(RECORD,&format!("{RECORD}\nlongitudinal-loss kelvin-voigt"))
}
fn build(cx: &Cx<'_>, gate: &CancelGate, input: &str, sleeve_eta: f64, band: f64) -> Result<Ready> {
    input::Spec::parse(input.as_bytes())?.compile(cx,gate,|name| Ok(match name {
        "wire.fsmatpk" => source_pack(name,1e5,band),
        "sleeve.fsmatpk" => source_pack(name,sleeve_eta,band),
        _ => pack(name,7800.0,200_000.0,None), // Body does NOT need axial loss data.
    }))
}

#[test]
fn sleeve_viscosity_alone_changes_axial_decay_and_actual_coupled_pressure() {
    scoped(|cx,gate| {
        let input = selected_text();
        let mut first = build(cx,gate,&input,1e5,200_000.0).unwrap();
        let mut changed = build(cx,gate,&input,4e5,200_000.0).unwrap();
        assert_eq!(first.spec.bytes,changed.spec.bytes);
        assert_eq!(first.catalog["wire.fsmatpk"].pack.content_hash(),changed.catalog["wire.fsmatpk"].pack.content_hash());
        assert_ne!(first.catalog["sleeve.fsmatpk"].pack.content_hash(),changed.catalog["sleeve.fsmatpk"].pack.content_hash());
        let a = first.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0].clone();
        let b = changed.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0].clone();
        assert_eq!(a.model().axial_omegas(),b.model().axial_omegas());
        assert_eq!(a.solid_linear_density_kg_m(),b.solid_linear_density_kg_m());
        assert_eq!(a.axial_rigidity_n(),b.axial_rigidity_n());
        let loss = b.axial_material_loss().unwrap();
        assert_eq!(loss.viscosity_property(),AXIAL_ETA);
        assert_eq!(loss.elastic().layers().len(),2);
        let g = loss.elastic().layer_integrals();
        let expected = 1e5*g[0].area_m2()+4e5*g[1].area_m2();
        assert!((loss.viscous_rigidity_n_s()-expected).abs() < 1e-13*expected);
        let ai = a.axial_component();
        assert_eq!(ai,b.axial_component());
        for (i,model) in first.performance.bank().system().components().iter().enumerate() {
            if i != ai { assert_eq!(model.modes(),changed.performance.bank().system().components()[i].modes()); }
        }
        let modes = changed.performance.bank().system().components()[ai].modes();
        for mode in modes {
            let expected = loss.longitudinal_damping_ratio(mode.angular_frequency_rad_s).unwrap();
            assert_eq!(mode.damping_ratio,expected);
            assert!(mode.damping_ratio>0.0);
        }
        let x = pressure(&mut first,37);
        let y = pressure(&mut changed,128);
        assert!(x.iter().chain(&y).all(|p| p.is_finite()));
        assert!(x.iter().zip(&y).any(|(x,y)|(x-y).abs()>1e-18));
        assert_eq!(y,pressure(&mut build(cx,gate,&input,4e5,200_000.0).unwrap(),1));
        let frame = changed.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn explicit_zero_preserves_elastic_pressure_but_missing_and_out_of_band_loss_refuse() {
    scoped(|cx,gate| {
        let input = selected_text();
        let zero = |name: &str| Ok(if name=="body.fsmatpk" { pack(name,7800.0,200_000.0,None) }
            else { source_pack(name,0.0,200_000.0) });
        let mut damped = input::Spec::parse(input.as_bytes()).unwrap().compile(cx,gate,zero).unwrap();
        let mut elastic = input::Spec::parse(text().as_bytes()).unwrap().compile(cx,gate,zero).unwrap();
        assert!(damped.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0]
            .axial_material_loss().is_some());
        assert!(elastic.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0]
            .axial_material_loss().is_none());
        assert_eq!(pressure(&mut damped,37),pressure(&mut elastic,128));
        // Existing elastic source has a wide band; only axial loss is too narrow.
        assert!(build(cx,gate,&input,4e5,50_000.0).is_err());
        assert!(build(cx,gate,&input,-1.0,200_000.0).is_err());
        let missing_sleeve = input::Spec::parse(input.as_bytes()).unwrap().compile(cx,gate,|name| {
            Ok(if name=="wire.fsmatpk" { source_pack(name,1e5,200_000.0) }
                else { pack(name,7800.0,200_000.0,None) })
        });
        assert!(missing_sleeve.is_err());
    });
}

#[test]
fn source_selection_is_explicit_and_cannot_be_ignored_without_axial_coordinates() {
    for record in ["longitudinal-loss fitted-q", "longitudinal-loss kelvin-voigt extra",
        "longitudinal-loss kelvin-voigt\nlongitudinal-loss elastic"] {
        assert!(input::Spec::parse(text().replace(RECORD,&format!("{RECORD}\n{record}")).as_bytes()).is_err());
    }
    assert!(input::Spec::parse(text().replace(RECORD,"longitudinal-loss kelvin-voigt").as_bytes()).is_err());
    scoped(|cx,gate| {
        let mut s=input::Spec::parse(selected_text().as_bytes()).unwrap();
        s.longitudinal_waves=None;
        let mut calls=0;
        assert!(s.compile(cx,gate,|name| {calls+=1; Ok(pack(name,7800.0,200_000.0,None))}).is_err());
        assert_eq!(calls,0,"refuse ignored source physics before pack I/O");
    });
}

#[test]
fn finite_strain_selection_changes_material_reaction_and_shared_pressure_without_double_loss() {
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::ModalPotential;
    scoped(|cx,gate| {
        let input = selected_text().replace("longitudinal-loss kelvin-voigt",
            "longitudinal-loss finite-strain-kelvin-voigt");
        let mut a = build(cx,gate,&input,1e5,200_000.0).unwrap();
        let mut b = build(cx,gate,&input,4e5,200_000.0).unwrap();
        let law = b.performance.bank().source().bank().longitudinal_waves().unwrap().strings()[0].clone();
        assert!(law.has_finite_strain_viscosity());
        assert!(b.performance.bank().system().components()[law.axial_component()].modes()
            .iter().all(|mode|mode.damping_ratio==0.0));
        // Pure transverse instantaneous motion stretches material even with
        // every axial velocity zero. Linear u_xt loss would be exactly zero.
        let mut q = vec![0.0;law.coordinate_count()];
        let mut v = q.clone();
        let mut reaction = q.clone();
        q[0]=1e-6;v[0]=0.01;
        law.dissipative_reaction(&q,&v,&mut reaction).unwrap();
        assert!(reaction[0]<0.0);
        assert!(-reaction.iter().zip(&v).map(|(f,v)|f*v).sum::<f64>()>0.0);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        let first=pressure(&mut a,37);
        let changed=pressure(&mut b,128);
        assert!(first.iter().zip(&changed).any(|(x,y)|(x-y).abs()>1e-18));
        assert_eq!(changed,pressure(&mut build(cx,gate,&input,4e5,200_000.0).unwrap(),1));
        let frame=b.performance.bank().system().last_frame().unwrap();
        assert!(frame.nonlinear.unwrap().nonlinear_dissipation_j>=0.0);
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        assert!(build(cx,gate,&input,4e5,50_000.0).is_err());
        assert!(build(cx,gate,&input,-1.0,200_000.0).is_err());
        assert!(input::Spec::parse(input.as_bytes()).unwrap().compile(cx,gate,|name|
            Ok(pack(name,7800.0,200_000.0,None))).is_err());
    });
}

#[test]
fn zero_finite_strain_viscosity_preserves_the_existing_elastic_run() {
    scoped(|cx,gate| {
        let input = selected_text().replace("longitudinal-loss kelvin-voigt",
            "longitudinal-loss finite-strain-kelvin-voigt");
        let source = |name:&str| Ok(if name=="body.fsmatpk" { pack(name,7800.0,200_000.0,None) }
            else { source_pack(name,0.0,200_000.0) });
        let mut loss = input::Spec::parse(input.as_bytes()).unwrap().compile(cx,gate,source).unwrap();
        let mut elastic = input::Spec::parse(text().as_bytes()).unwrap().compile(cx,gate,source).unwrap();
        assert_eq!(pressure(&mut loss,37),pressure(&mut elastic,128));
    });
}
