//! Synthetic source data changes the ordinary coupled performer and PCM.
use super::*;

const LOSS: &str = "loss regional-maxwell relaxed 32 1000000 1 maxwell_modulus maxwell_time";
fn body_pack(name: &str, modulus: f64, time: f64, band: f64, omit: Option<&str>) -> NormalizedPack {
    // One missing tensor-viscosity coefficient deliberately proves that the
    // Maxwell selector does not keep the old loss path in parallel.
    let base = pack(name, 0.0, Some(VISCOSITY[0]));
    let mut claims = base.claims().clone();
    for (key,dims,value) in [("maxwell_modulus",Pressure::DIMS,modulus),
        ("maxwell_time",Dims([0,0,1,0,0,0]),time)] {
        if omit == Some(key) { continue; }
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(key,dims), value:PropertyValue::Scalar { value,dims },
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,band),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity, observations:vec![],
            provenance:Provenance { source:"synthetic body material memory".into(),license:"CC0-1.0".into(),artifact:None },
        }).unwrap();
    }
    NormalizedPack::new(base.pack_id(),"synthetic-memory-v1",base.source_artifact(),
        base.redistribution_terms(),claims,vec![],vec![]).unwrap()
}
fn memory_profile(regional: bool) -> Profile {
    let text = if regional {
        include_str!("../regional.body").replace("region right @body","region right right.fsmatpk")
    } else { include_str!("../orthotropic.body").into() };
    Profile::parse(format!("{text}{LOSS}\n").as_bytes()).unwrap()
}
fn memory_ready(cx: &Cx<'_>, gate: &CancelGate, s: input::Spec, profile: Profile,
    right_modulus: f64) -> Ready
{
    s.compile_with_body(cx,gate,Some(profile),|name| Ok(match name {
        "body.fsmatpk" => body_pack(name,4e9,0.002,200_000.0,None),
        "right.fsmatpk" => body_pack("body.fsmatpk",right_modulus,0.0005,200_000.0,None),
        _ => pack(name,1.0,None),
    })).unwrap()
}

#[test]
fn regional_material_memory_changes_bridge_reactions_and_streamed_pressure() {
    scoped(|cx,gate| {
        let build = |modulus| memory_ready(cx,gate,spec(),memory_profile(true),modulus);
        let mut a=build(2e9); let mut changed=build(7e9);
        let a_source=a.performance.bank().source().bank();
        let b_source=changed.performance.bank().source().bank();
        assert!(a_source.material_chart_viscosity().is_none());
        let memory=a_source.material_chart_relaxation().unwrap();
        assert_eq!(memory.source().material_spectra().len(),2);
        assert_eq!(memory.coupling_rows_per_s().len(),6);
        assert_eq!(a.performance.bank().system().material_memory().len(),6);
        assert!(memory.coupling_rows_per_s().iter().any(|row|
            row.iter().filter(|c| c.abs()>1e-12).count()>1));
        assert_eq!(a_source.material_chart().unwrap().specimen.regions()[0].input.material.identity(),
            b_source.material_chart().unwrap().specimen.regions()[0].input.material.identity());
        assert_ne!(memory.source().material_spectra()[1].terms,
            b_source.material_chart_relaxation().unwrap().source().material_spectra()[1].terms);
        assert_eq!(a.performance.bank().system().components()[1].modes(),
            changed.performance.bank().system().components()[1].modes());
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),
            changed.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),
            changed.catalog["body.fsmatpk"].pack.content_hash());
        let p=pressure(&mut a,37);let q=pressure(&mut changed,128);
        assert!(p.iter().any(|v|v.abs()>1e-12));
        assert!(p.iter().zip(&q).any(|(a,b)|(a-b).abs()>1e-12));
        assert_ne!(a.performance.bank().system().components()[0].states(),
            changed.performance.bank().system().components()[0].states(),
            "material history must react through the bridge onto the original string");
        let mut replay=build(2e9);assert_eq!(p,pressure(&mut replay,1));
        let expected=fs_couple::pcm_wav::encode_pcm16_wav(&p,a.spec.rate,a.spec.full_scale).unwrap().0;
        let altered=fs_couple::pcm_wav::encode_pcm16_wav(&q,changed.spec.rate,changed.spec.full_scale).unwrap().0;
        assert_ne!(expected[44..],altered[44..],"source substitution must reach actual PCM samples");
        let (wav,summary)=build(2e9).write(std::io::Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),expected);assert_eq!(summary.samples,256);
        for image in [&a,&changed,&replay] {
            let frame=image.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        }
    });
}

#[test]
fn explicit_zero_spectrum_keeps_elastic_playback_and_retains_original_sources() {
    scoped(|cx,gate| {
        let mut bare=compile(cx,gate,Some(profile()),0.0);
        let mut zero=spec().compile_with_body(cx,gate,Some(memory_profile(false)),|name|
            Ok(if name=="body.fsmatpk" {body_pack(name,0.0,0.002,200_000.0,None)}
                else {pack(name,0.0,None)})).unwrap();
        let source=zero.performance.bank().source().bank();
        let memory=source.material_chart_relaxation().unwrap();
        assert!(memory.coupling_rows_per_s().is_empty());
        assert_eq!(memory.source().material_spectra()[0].terms,vec![(0.0,0.002)]);
        assert!(source.material_chart_viscosity().is_none());
        assert_eq!(pressure(&mut bare,128),pressure(&mut zero,37));
    });
}

#[test]
fn body_memory_is_bound_after_original_bonded_torsional_gas_and_axial_coordinates() {
    scoped(|cx,gate| {
        let text=include_str!("../longitudinal-cylinder-bow.material-strings")
            .replace("48000 4096 0.01 128","48000 64 0.01 64")
            .replace("prepared-linear 64 67108864 200000000000\n","");
        let build=||memory_ready(cx,gate,input::Spec::parse(text.as_bytes()).unwrap(),memory_profile(false),2e9);
        let mut a=build();let mut b=build();
        let source=a.performance.bank().source().bank();
        assert!(source.longitudinal_waves().is_some());
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        assert_eq!(source.material_chart_relaxation().unwrap().mode_count(),2);
        assert_eq!(a.performance.bank().system().components().iter().map(|c|c.modes().len())
            .collect::<Vec<_>>(),vec![3,2,1,12,2]);
        let p=pressure(&mut a,1);assert_eq!(p,pressure(&mut b,37));
        assert!(p.iter().any(|v|v.abs()>1e-15));
        let frame=a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn incomplete_unresolved_and_over_budget_material_history_refuses() {
    scoped(|cx,gate| {
        for (modulus,time,band,omit) in [(4e9,0.002,200_000.0,Some("maxwell_time")),
            (-1.0,0.002,200_000.0,None),(4e9,0.0,200_000.0,None),
            (4e9,0.002,2000.0,None)] {
            assert!(spec().compile_with_body(cx,gate,Some(memory_profile(false)),|name|
                Ok(if name=="body.fsmatpk" {body_pack(name,modulus,time,band,omit)}
                    else {pack(name,1.0,None)})).is_err());
        }
        let tiny=Profile::parse(format!("{}{}\n",include_str!("../orthotropic.body"),
            LOSS.replace("relaxed 32","relaxed 1")).as_bytes()).unwrap();
        assert!(spec().compile_with_body(cx,gate,Some(tiny),|name|Ok(if name=="body.fsmatpk" {
            body_pack(name,4e9,0.002,200_000.0,None)} else {pack(name,1.0,None)})).is_err());
        let prepared=input::Spec::parse(include_bytes!("../longitudinal-cylinder-bow.material-strings")).unwrap();
        let mut reads=0;
        assert!(prepared.compile_with_body(cx,gate,Some(memory_profile(false)),|name| {
            reads+=1;Ok(pack(name,1.0,None))
        }).is_err());
        assert_eq!(reads,0,"unsupported prepared owner must refuse before material I/O");
    });
}

#[test]
fn optional_loss_grammar_is_bounded_and_never_ignores_extra_fields() {
    let body=include_str!("../orthotropic.body");
    assert!(Profile::parse(format!("{body}loss tensor\n").as_bytes()).is_ok());
    for loss in ["loss wood","loss tensor ignored",
        "loss regional-maxwell relaxed 257 1000000 1 e tau",
        "loss regional-maxwell relaxed 32 0 1 e tau",
        "loss regional-maxwell relaxed 32 1000000 999999999 e tau",
        "loss regional-maxwell guessed 32 1000000 1 e tau",
        "loss regional-maxwell relaxed 32 1000000 1 e",
        "loss regional-maxwell relaxed 32 1000000 1 e tau ignored",
        "loss regional-maxwell relaxed 32 1000000 2 e tau e tau"] {
        assert!(Profile::parse(format!("{body}{loss}\n").as_bytes()).is_err(),"accepted {loss}");
    }
    assert!(Profile::parse(format!("{body}{LOSS}\n{LOSS}\n").as_bytes()).is_err());
}
