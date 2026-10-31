//! Original specimen spectra share the performer and PCM with layered/body memory.
//! Manufactured sources provide integration evidence, not measured instrument data.
use super::*;
use fs_couple::string_specimen::{EQUILIBRIUM_YOUNG_MODULUS_PROPERTY as EQUILIBRIUM,
    ResolvedStringSpecimen};
use fs_scenario::acoustic::BendingRelaxationProperties;

const LOSS: &str = "bending-loss prony 1 relaxing_bending_modulus bending_relaxation_time";

fn source_pack(name: &str, delta: f64, time: f64, band: f64, omit: Option<&str>) -> NormalizedPack {
    let exclusive = matches!(name, "wire.fsmatpk" | "relaxing-wire.fsmatpk" | "layered-sleeve.fsmatpk");
    let base = pack(name,7800.0,200_000.0,exclusive.then_some(ETA));
    let mut claims = base.claims().clone();
    for (key,dims,value) in [(EQUILIBRIUM,Pressure::DIMS,200e9),
        ("relaxing_bending_modulus",Pressure::DIMS,delta),
        ("bending_relaxation_time",Dims([0,0,1,0,0,0]),time)] {
        if omit == Some(key) { continue; }
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(key,dims), value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,band),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity, observations:vec![],
            provenance:Provenance {source:"synthetic string memory integration".into(),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
    }
    NormalizedPack::new(base.pack_id(),"synthetic-string-memory-v1",base.source_artifact(),
        base.redistribution_terms(),claims,vec![],vec![]).unwrap()
}

fn selected() -> input::Spec {
    let text = include_str!("../pluck.material-strings")
        .replace("48000 4096 0.01 128","48000 512 0.0001 128")
        .replace("cavity neck",&format!("{LOSS}\ncavity neck"));
    input::Spec::parse(text.as_bytes()).unwrap()
}

fn specimen(source: &MaterialStringSource) -> &ResolvedStringSpecimen {
    match source {
        MaterialStringSource::Homogeneous(s) => s,
        MaterialStringSource::WithTorsion(s) => specimen(s.bending()),
        MaterialStringSource::WithPolarizations(s) => specimen(s.source()),
        _ => panic!("fixture lost its original homogeneous source"),
    }
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, delta: f64) -> Ready {
    spec.compile(cx,gate,|name| Ok(source_pack(name,
        if name=="body.fsmatpk" {20e9} else {delta},0.001,200_000.0,None))).unwrap()
}

#[test]
fn source_spectrum_changes_shared_pressure_and_pcm_with_original_mode_geometry() {
    scoped(|cx,gate| {
        let mut a = ready(cx,gate,selected(),20e9);
        let mut b = ready(cx,gate,selected(),160e9);
        let source = a.performance.bank().source().bank();
        let wire = specimen(&source.strings()[0].source);
        assert!(wire.material().property(ETA).is_none(),"unselected viscosity was still required");
        let law = wire.string().relaxation_bending.unwrap();
        assert_eq!(law.material_state_identity,Some(wire.material().identity()));
        assert_eq!(law.source_properties.unwrap(),vec![BendingRelaxationProperties {
            modulus:"relaxing_bending_modulus".into(),relaxation_time:"bending_relaxation_time".into(),
        }]);
        assert_eq!(a.performance.bank().system().material_memory().len(),4);
        assert!(a.performance.bank().system().material_memory().iter().any(|z| *z != 0.0),
            "the shared static preload did not initialize relaxed string history");
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),
            b.catalog["body.fsmatpk"].pack.content_hash());
        for (left,right) in a.performance.bank().system().components().iter()
            .zip(b.performance.bank().system().components()) {
            assert_eq!(left.modes(),right.modes(),"memory selection changed the equilibrium pencil");
        }
        let p = pressure(&mut a,37);
        let q = pressure(&mut b,128);
        assert!(p.iter().any(|v| v.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(x,y)| (x-y).abs() > 1e-12));
        assert_eq!(p,pressure(&mut ready(cx,gate,selected(),20e9),1));
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p,a.spec.rate,a.spec.full_scale).unwrap().0;
        let changed = fs_couple::pcm_wav::encode_pcm16_wav(&q,b.spec.rate,b.spec.full_scale).unwrap().0;
        assert_ne!(expected[44..],changed[44..],"material changes never reached PCM samples");
        let (wav,report) = ready(cx,gate,selected(),20e9).write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),expected);
        assert_eq!(report.samples,512);
        for image in [&a,&b] {
            let frame=image.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        }
    });
}

#[test]
fn hereditary_bending_and_causal_fluid_keep_the_independent_tension_loss() {
    scoped(|cx, gate| {
        let mut s = selected();
        s.samples = 128;
        s.strings[0].bending_modes = 2;
        s.cylinder_gas_loading = Some(CylinderFlowResolution {
            minimum_wavenumber_radius: 1e-3, maximum_wavenumber_radius: 30.0, intervals: 2,
        });
        let mut a = ready(cx, gate, s.clone(), 80e9);
        let row = input::tension_loss_rows().unwrap().into_iter()
            .find(|row| row.id == "lampis-2024-dominant-cello-g").unwrap();
        s.strings[0].tension_loss = Some(row.clone());
        let mut b = ready(cx, gate, s, 80e9);
        let branches_a = a.performance.bank().system().material_relaxation_branches();
        let branches_b = b.performance.bank().system().material_relaxation_branches();
        assert_eq!(branches_a.len(), branches_b.len());
        for (a, b) in branches_a.iter().zip(branches_b) {
            assert_eq!(a.projection, b.projection);
            assert_eq!(a.stiffness, b.stiffness);
            assert_eq!(a.relaxation_time_s, b.relaxation_time_s);
        }
        assert_eq!(a.performance.bank().system().material_memory(),
            b.performance.bank().system().material_memory());
        let (plain, loss) = (&a.performance.bank().cards()[0], &b.performance.bank().cards()[0]);
        assert_eq!(loss.viscous_bending_n_m2_s, 0.0, "Prony storage acquired a second bending loss");
        for k in 1..=loss.mode_count {
            let wave = k as f64 * core::f64::consts::PI / loss.length_m;
            let omega = loss.mode_omega_rad_s(k);
            let expected = 0.5 * row.eta_f * loss.tension_n * wave * wave
                / (loss.linear_density_kg_m * omega * omega);
            assert!(((loss.zetas[k-1] - plain.zetas[k-1]) / expected - 1.0).abs() < 1e-10);
        }
        assert!(pressure(&mut a, 128).iter().zip(pressure(&mut b, 37))
            .any(|(p, q)| (p-q).abs() > 1e-15));
        let frame = b.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn string_and_body_memories_share_gas_torsion_axial_and_multibow_performer() {
    scoped(|cx,gate| {
        let mut s = selected();
        s.samples=96;
        s.strings[0].polarizations=Some(
            fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections::try_new(
                [1.0,0.0],[0.8,0.6]).unwrap());
        s.cylinder_gas_loading=Some(CylinderFlowResolution {
            minimum_wavenumber_radius:1e-3,maximum_wavenumber_radius:30.0,intervals:2,
        });
        s.longitudinal_waves=Some((LongitudinalStringResolution {
            axial_modes:vec![1],maximum_slope:0.2,maximum_strain:0.01,
        },NonlinearModalConfig {max_iterations:64,maximum_generalized_force:1e5,
            force_absolute_tolerance:1e-11,force_relative_tolerance:1e-10}));
        s.forces=None;
        s.schedule=GestureSchedule::try_new(1000,vec![GestureTrack {
            id:"bow".into(),target:GestureTarget::BowStroke {string:0},
            initial:GestureValue::Bow {velocity_m_per_s:0.01,normal_force_n:0.0001,station:0.2},
            events:vec![],
        }]).unwrap();
        let text=include_str!("../orthotropic.body").replace("elastic orthotropic12","elastic isotropic");
        let profile=body::Profile::parse(format!("{text}loss regional-maxwell relaxed 8 1000000 1 relaxing_bending_modulus bending_relaxation_time\n").as_bytes()).unwrap();
        let build = || s.clone().compile_with_body(cx,gate,Some(profile.clone()),|name|
            Ok(source_pack(name,if name=="body.fsmatpk" {20e9} else {80e9},0.001,200_000.0,None))).unwrap();
        let mut a=build(); let mut b=build();
        let source=a.performance.bank().source().bank();
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.longitudinal_waves().is_some());
        assert!(source.material_cavity().is_some());
        let body_states=source.material_chart_relaxation().unwrap().coupling_rows_per_s().len();
        assert_eq!(body_states,2);
        assert_eq!(a.performance.bank().system().material_memory().len(),8+body_states);
        let options=source.options();
        let n=options.network.max_modes;
        let links=options.network.max_connections;
        let reserved=8+8;
        assert!(options.network.max_setup_terms >= n*(links+1).pow(2)+(2*n+reserved).pow(2)+2*n*reserved);
        let p=pressure(&mut a,1); let q=pressure(&mut b,37);
        assert_eq!(p,q);
        assert!(p.iter().any(|v| *v != 0.0));
        assert_eq!(a.performance.bank().system().material_memory(),b.performance.bank().system().material_memory());
        assert!(a.performance.bank().system().material_memory()[..8].iter().any(|z| *z != 0.0));
        assert!(a.performance.bank().system().material_memory()[8..].iter().any(|z| *z != 0.0));
    });
}

#[test]
fn exclusive_and_additional_string_spectra_share_one_body_history_and_budget() {
    scoped(|cx,gate| {
        let mut s=selected();s.samples=128;
        s.strings[0].bending_modes=2;s.strings[0].twist_modes=1;
        let mut other=s.strings[0].clone();
        other.pack="parallel.fsmatpk".into();
        other.bending_loss=input::BendingLoss::KelvinVoigt;
        other.prestress=StringPrestress::FixedTension(22.0);
        s.strings.push(other);
        s.bending_memory=Some(bending_memory::Input::new(4,vec![
            bending_memory::Input::uniform_selection(1,
                vec![("relaxing_bending_modulus".into(),"bending_relaxation_time".into())],
                &s.strings).unwrap(),
        ],&s.strings).unwrap());
        let text=include_str!("../orthotropic.body")
            .replace("elastic orthotropic12","elastic isotropic");
        let profile=body::Profile::parse(format!("{text}loss regional-maxwell relaxed 4 1000000 1 relaxing_bending_modulus bending_relaxation_time\n").as_bytes()).unwrap();
        let build=|tau|s.clone().compile_with_body(cx,gate,Some(profile.clone()),|name|
            Ok(source_pack(name,if name=="body.fsmatpk" {20e9} else {80e9},
                if name=="parallel.fsmatpk" {tau} else {0.001},200_000.0,None))).unwrap();
        let mut a=build(0.002);let mut b=build(0.008);
        let source=a.performance.bank().source().bank();
        assert_eq!(source.string_bending_relaxation().unwrap().channels().len(),2);
        assert_eq!(source.material_chart_relaxation().unwrap().coupling_rows_per_s().len(),2);
        assert_eq!(a.performance.bank().system().material_memory().len(),6);
        assert!(specimen(&source.strings()[0].source).string().relaxation_bending.is_some());
        assert!(specimen(&source.strings()[1].source).string().kelvin_voigt_bending.is_some());
        let options=source.options();let n=options.network.max_modes;let links=options.network.max_connections;
        let reserved=2+4+4;
        assert!(options.network.max_setup_terms >= n*(links+1).pow(2)+(2*n+reserved).pow(2)+2*n*reserved);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),b.catalog["body.fsmatpk"].pack.content_hash());
        let p=pressure(&mut a,1);let q=pressure(&mut b,37);
        assert!(p.iter().zip(&q).any(|(x,y)|(x-y).abs()>1e-15));
        assert_eq!(p,pressure(&mut build(0.002),37));
        let frame=a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn selected_sources_unsupported_models_and_duplicate_laws_refuse_before_they_can_be_ignored() {
    let text=String::from_utf8(selected().bytes).unwrap();
    for bad in ["bending-loss prony 0","bending-loss prony 9",
        "bending-loss prony 1 duplicate duplicate",
        "bending-loss prony 2 e tau e tau","bending-loss kelvin-voigt ignored"] {
        assert!(input::Spec::parse(text.replace(LOSS,bad).as_bytes()).is_err());
    }
    let duplicate="bending-memory 4 1\nstring-memory 0 1 relaxing_bending_modulus bending_relaxation_time\n";
    assert!(input::Spec::parse(text.replace("cavity ",&format!("{duplicate}cavity ")).as_bytes()).is_err());
    assert!(input::Spec::parse(include_bytes!("../relaxing-strings.material-strings")).is_ok());
    scoped(|cx,gate| {
        for missing in [EQUILIBRIUM,"relaxing_bending_modulus","bending_relaxation_time"] {
            assert!(selected().compile(cx,gate,|name|
                Ok(source_pack(name,80e9,0.001,200_000.0,Some(missing)))).is_err());
        }
        for (delta,time,band) in [(-1.0,0.001,200_000.0),(80e9,0.0,200_000.0),
            (80e9,1e-9,200_000.0),(80e9,0.001,1000.0)] {
            assert!(selected().compile(cx,gate,|name|Ok(source_pack(name,delta,time,band,None))).is_err());
        }
        let mut unsupported=selected();
        unsupported.prepared_linear=Some(PreparedModalBudget {
            max_modes:64,max_dense_bytes:16_777_216,max_multiply_terms:16_777_216,
        });
        let plain=spec();
        let mut twice=selected();
        twice.bending_memory=Some(bending_memory::Input::new(4,vec![
            bending_memory::Input::uniform_selection(0,
                vec![("relaxing_bending_modulus".into(),"bending_relaxation_time".into())],
                &plain.strings).unwrap(),
        ],&plain.strings).unwrap());
        let mut too_many=selected();
        let mut other=too_many.strings[0].clone();other.bending_loss=input::BendingLoss::KelvinVoigt;
        too_many.strings.push(other);
        too_many.bending_memory=Some(bending_memory::Input::new(253,vec![
            bending_memory::Input::uniform_selection(1,
                vec![("relaxing_bending_modulus".into(),"bending_relaxation_time".into())],
                &too_many.strings).unwrap(),
        ],&too_many.strings).unwrap());
        for rejected in [unsupported,twice,too_many] {
            let mut calls=0;
            assert!(rejected.compile(cx,gate,|name| {
                calls+=1; Ok(source_pack(name,80e9,0.001,200_000.0,None))
            }).is_err());
            assert_eq!(calls,0);
        }
    });
}


// Native selectors retain the same manufactured material owners used by the
// piano loss-angle and circular thermoelastic paths. They do not select a
// viscosity, replace the body, or imply measured instrument qualification.
const INTRINSIC: &str = "intrinsic_bending_loss";
const EXPANSION: &str = "linear_thermal_expansion_coefficient";
const CAPACITY: &str = "specific_heat_capacity";
const CONDUCTIVITY: &str = "thermal_conductivity";

fn modal_spec(selection: &str) -> input::Spec {
    let text = String::from_utf8(selected().bytes).unwrap().replace(LOSS, selection);
    input::Spec::parse(text.as_bytes()).unwrap()
}

fn modal_pack(name: &str, tangent: f64, tangent_dims: Dims, band: f64,
    omit: Option<&str>) -> NormalizedPack
{
    modal_pack_with_temperature(name, tangent, tangent_dims, band, omit, true)
}

fn modal_pack_with_temperature(name: &str, tangent: f64, tangent_dims: Dims,
    band: f64, omit: Option<&str>, typed_temperature: bool) -> NormalizedPack
{
    if name != "wire.fsmatpk" {
        return pack(name, 7800.0, 200_000.0, None);
    }
    let base = pack(name, 7800.0, 200_000.0, Some(ETA));
    let kelvin = fs_qty::QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
        fs_qty::semantic::QuantityKind::AbsoluteTemperature, fs_qty::semantic::ValueForm::Static));
    let validity = |maximum| {
        let domain = if typed_temperature {
            ValidityDomain::unconstrained().with_quantity("T", kelvin, 280.0, 320.0)
        } else { ValidityDomain::unconstrained().with("T", 280.0, 320.0) };
        domain.with("omega", 0.0, maximum)
    };
    // These manufactured sources declare their own temperature semantics.
    // Every original wire claim and observation remains present; the native
    // material query must preserve this declaration rather than invent it.
    let mut claims = ClaimSet::new();
    for id in base.claims().observation_ids() {
        assert_eq!(claims.register_observation(
            base.claims().observation(id).unwrap().clone()).unwrap(), id);
    }
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = validity(200_000.0);
        claims.insert_claim(claim).unwrap();
    }
    for (key, dims, value, maximum) in [
        (INTRINSIC, tangent_dims, tangent, band),
        (EXPANSION, fs_material::state_point::INVERSE_TEMPERATURE_DIMS, 12e-6, 200_000.0),
        (CAPACITY, fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS, 450.0, 150_000.0),
        (CONDUCTIVITY, fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS, 50.0, 200_000.0),
    ] {
        if omit == Some(key) { continue; }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: validity(maximum),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "manufactured native intrinsic-loss integration".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    let source = NormalizedPack::new(base.pack_id(), "synthetic-intrinsic-loss-v1",
        base.source_artifact(), base.redistribution_terms(), claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&source.to_bytes()).unwrap()
}

fn modal_source(source: &MaterialStringSource)
    -> &fs_material::section::modal_loss::ConcentricModalLossSection
{
    match source {
        MaterialStringSource::WithPolarizations(p) => modal_source(p.source()),
        MaterialStringSource::WithTorsion(t) => modal_source(t.bending()),
        MaterialStringSource::ConcentricModalLoss { section, .. } => section,
        _ => panic!("native selection lost its original intrinsic bending law"),
    }
}

#[test]
fn intrinsic_loss_selection_changes_decay_and_pressure_without_viscosity_or_retuning() {
    scoped(|cx, gate| {
        let selection = format!("bending-loss loss-angle {INTRINSIC}");
        let s = modal_spec(&selection);
        assert!(matches!(&s.strings[0].bending_loss,
            input::BendingLoss::LossAngle { property } if property == INTRINSIC));
        let build = |tangent| s.clone().compile(cx, gate, |name|
            Ok(modal_pack(name, tangent, Dims::NONE, 100_000.0, None))).unwrap();
        let mut a = build(0.002);
        let mut b = build(0.02);
        let law = modal_source(&a.performance.bank().source().bank().strings()[0].source);
        assert_eq!(law.loss_angle().unwrap().loss_property(), INTRINSIC);
        assert!(law.thermoelastic().is_none());
        assert_eq!(law.omega_band_rad_s(), (0.0, 100_000.0));
        let material = &law.elastic().layers()[0].material;
        assert!(material.property(ETA).is_none(), "unselected bending viscosity was still required");
        assert_eq!(material.property(INTRINSIC).unwrap().value_si(), 0.002);
        let (left, right) = (&a.performance.bank().cards()[0], &b.performance.bank().cards()[0]);
        assert_eq!(left.viscous_bending_n_m2_s, 0.0);
        assert_eq!(right.viscous_bending_n_m2_s, 0.0);
        assert_eq!(left.tension_n, right.tension_n);
        assert_eq!(left.linear_density_kg_m, right.linear_density_kg_m);
        assert_eq!(left.bending_stiffness_n_m2, right.bending_stiffness_n_m2);
        for k in 1..=left.mode_count {
            assert_eq!(left.mode_omega_rad_s(k), right.mode_omega_rad_s(k));
            assert!(right.zetas[k - 1] > left.zetas[k - 1]);
        }
        assert_eq!(a.spec.bytes, b.spec.bytes);
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),
            b.catalog["body.fsmatpk"].pack.content_hash());
        let p = pressure(&mut a, 37);
        let q = pressure(&mut b, 128);
        assert!(p.iter().any(|p| p.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(p, q)| (p - q).abs() > 1e-12));
        assert_eq!(p, pressure(&mut build(0.002), 1));
        for ready in [&a, &b] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn thermal_and_excess_selectors_keep_original_sources_bands_and_loaded_inertia() {
    scoped(|cx, gate| {
        let mut only = modal_spec("bending-loss thermoelastic");
        let mut excess = modal_spec(&format!("bending-loss thermoelastic excess {INTRINSIC}"));
        assert!(matches!(&only.strings[0].bending_loss,
            input::BendingLoss::Thermoelastic { excess_property: None }));
        assert!(matches!(&excess.strings[0].bending_loss,
            input::BendingLoss::Thermoelastic { excess_property: Some(property) }
                if property == INTRINSIC));
        let resolution = CylinderFlowResolution {
            minimum_wavenumber_radius: 1e-3, maximum_wavenumber_radius: 30.0, intervals: 2,
        };
        only.cylinder_gas_loading = Some(resolution);
        excess.cylinder_gas_loading = Some(resolution);
        let build = |s: input::Spec| s.compile(cx, gate, |name|
            Ok(modal_pack(name, 0.003, Dims::NONE, 100_000.0, None))).unwrap();
        let a = build(only);
        let b = build(excess);
        let thermal_only = modal_source(&a.performance.bank().source().bank().strings()[0].source);
        let combined = modal_source(&b.performance.bank().source().bank().strings()[0].source);
        assert!(thermal_only.loss_angle().is_none());
        assert_eq!(thermal_only.omega_band_rad_s(), (0.0, 150_000.0));
        assert_eq!(combined.omega_band_rad_s(), (0.0, 100_000.0));
        assert!(thermal_only.hysteretic_stiffness_n_m2(120_000.0).is_ok());
        assert!(combined.hysteretic_stiffness_n_m2(120_000.0).is_err(),
            "selected excess disappeared beyond its source support");
        let thermal = combined.thermoelastic().unwrap();
        let loss = combined.loss_angle().unwrap();
        assert_eq!(loss.elastic(), thermal.elastic());
        assert_eq!(loss.loss_property(), INTRINSIC);
        assert_eq!(thermal.law().alpha_t, 12e-6);
        assert_eq!(thermal.law().cp, 450.0);
        assert_eq!(thermal.law().conductivity, 50.0);
        assert_eq!(thermal.law().t0, 293.15);
        let state = &combined.elastic().layers()[0].material;
        assert!(state.property(ETA).is_none());
        let kelvin = fs_qty::QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
            fs_qty::semantic::QuantityKind::AbsoluteTemperature, fs_qty::semantic::ValueForm::Static));
        assert_eq!(state.axis_quantities().get("T"), Some(&kelvin));
        let (plain, total) = (&a.performance.bank().cards()[0], &b.performance.bank().cards()[0]);
        assert!(total.linear_density_kg_m > combined.elastic().line_density_kg_m());
        assert_eq!(plain.linear_density_kg_m, total.linear_density_kg_m);
        assert_eq!(total.viscous_bending_n_m2_s, 0.0);
        for k in 1..=total.mode_count {
            let omega = total.mode_omega_rad_s(k);
            assert_eq!(omega, plain.mode_omega_rad_s(k));
            let wave = k as f64 * core::f64::consts::PI / total.length_m;
            let expected = 0.5 * loss.hysteretic_stiffness_n_m2() * wave.powi(4)
                / total.linear_density_kg_m / omega / omega;
            assert!(((total.zetas[k - 1] - plain.zetas[k - 1]) / expected - 1.0).abs() < 1e-10);
        }
    });
}

#[test]
fn missing_invalid_or_out_of_band_intrinsic_sources_refuse_without_fallback() {
    scoped(|cx, gate| {
        let bare = modal_spec("bending-loss thermoelastic").compile(cx, gate, |name|
            Ok(modal_pack_with_temperature(name, 0.003, Dims::NONE, 100_000.0, None, false)))
            .err().expect("bare source temperature must not be relabelled as absolute kelvin");
        assert!(bare.to_string().contains("absolute-temperature"), "{bare}");
        for (selection, missing) in [
            (format!("bending-loss loss-angle {INTRINSIC}"), INTRINSIC),
            ("bending-loss thermoelastic".into(), EXPANSION),
            ("bending-loss thermoelastic".into(), CAPACITY),
            ("bending-loss thermoelastic".into(), CONDUCTIVITY),
            (format!("bending-loss thermoelastic excess {INTRINSIC}"), INTRINSIC),
        ] {
            assert!(modal_spec(&selection).compile(cx, gate, |name|
                Ok(modal_pack(name, 0.003, Dims::NONE, 100_000.0, Some(missing)))).is_err());
        }
        for selection in [
            format!("bending-loss loss-angle {INTRINSIC}"),
            format!("bending-loss thermoelastic excess {INTRINSIC}"),
        ] {
            for (tangent, dims, band) in [
                (-0.001, Dims::NONE, 100_000.0),
                (0.003, DynViscosity::DIMS, 100_000.0),
                (0.003, Dims::NONE, 1000.0),
            ] {
                assert!(modal_spec(&selection).compile(cx, gate, |name|
                    Ok(modal_pack(name, tangent, dims, band, None))).is_err());
            }
        }
    });
}

#[test]
fn malformed_native_loss_and_unsupported_winding_selections_refuse() {
    let text = String::from_utf8(selected().bytes).unwrap();
    for invalid in ["bending-loss loss-angle", "bending-loss loss-angle phi ignored",
        "bending-loss thermoelastic ignored", "bending-loss thermoelastic excess",
        "bending-loss thermoelastic excess phi ignored"]
    {
        assert!(input::Spec::parse(text.replace(LOSS, invalid).as_bytes()).is_err());
    }
    let memory = "bending-memory 4 1\nstring-memory 0 1 relaxing_bending_modulus bending_relaxation_time\n";
    for selection in [
        format!("bending-loss loss-angle {INTRINSIC}"),
        "bending-loss thermoelastic".into(),
        format!("bending-loss thermoelastic excess {INTRINSIC}"),
    ] {
        let duplicate = text.replace(LOSS, &selection)
            .replace("cavity ", &format!("{memory}cavity "));
        assert!(input::Spec::parse(duplicate.as_bytes()).is_err(),
            "another bending spectrum must not replace or duplicate selected modal loss");
    }
    scoped(|cx, gate| {
        let mut wound = modal_spec("bending-loss kelvin-voigt");
        wound.strings[0].twist_modes = 0;
        wound.strings[0].section = input::WireSection::CloseWound {
            core_radius_m: 0.0002, pack: "wrap.fsmatpk".into(),
            wire_radius_m: 0.00005, pitch_m: 0.000101,
        };
        assert!(wound.clone().compile(cx, gate, |name|
            Ok(pack(name, 7800.0, 200_000.0, None))).is_ok(),
            "comparison winding geometry must be independently admitted");
        for selection in [
            input::BendingLoss::LossAngle { property: INTRINSIC.into() },
            input::BendingLoss::Thermoelastic { excess_property: None },
            input::BendingLoss::Thermoelastic { excess_property: Some(INTRINSIC.into()) },
        ] {
            let mut selected = wound.clone();
            selected.strings[0].bending_loss = selection;
            let error = selected.compile(cx, gate, |name|
                Ok(modal_pack(name, 0.003, Dims::NONE, 100_000.0, None)))
                .err().expect("a close-wound core needs its own retained loss model");
            assert!(error.to_string().contains(
                "modal bending loss requires an original homogeneous or bonded material section"),
                "{error}");
        }
    });
}

#[path = "relaxation/layered.rs"]
mod layered;
