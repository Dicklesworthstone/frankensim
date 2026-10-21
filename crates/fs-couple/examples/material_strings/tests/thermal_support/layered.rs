//! G1/G3/G5: each bonded layer's expansion changes the same coupled sound.
use super::*;
use fs_material::section::ConcentricKelvinVoigtSection;

const NAMES: [&str; 3] = ["wire.fsmatpk", "inner.fsmatpk", "outer.fsmatpk"];
const RADII: [f64; 4] = [0.0, 0.0002, 0.00025, 0.0003];
const YOUNG: [f64; 3] = [200e9, 70e9, 80e9];
const DENSITY: [f64; 3] = [7800.0, 2700.0, 19000.0];
const SUPPORT: &str = "thermal-supports 1\nthermal-support 0 290\n";

fn layered_input() -> String {
    include_str!("../../layered-thermal-bow.material-strings")
        .replace("48000 4096", "48000 256")
        .replace("layer wire.fsmatpk 0.00025", "layer inner.fsmatpk 0.00025")
        .replace("layer wire.fsmatpk 0.0003", "layer outer.fsmatpk 0.0003")
}

fn layer_pack(index: usize, factor: f64) -> NormalizedPack {
    let base = pack(NAMES[index], DENSITY[index], 200_000.0, Some("young_modulus"));
    let mut claims = base.claims().clone();
    let observations = claims.claims_for("density")[0].1.observations.clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("young_modulus", Pressure::DIMS),
        value: PropertyValue::Scalar { value: YOUNG[index], dims: Pressure::DIMS },
        validity: ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations,
        provenance: Provenance { source: "manufactured layered elasticity".into(),
            license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([9;32])) },
    }).unwrap();
    let base = NormalizedPack::new(NAMES[index], "synthetic-layer-v1", fs_blake3::ContentHash([9;32]),
        "CC0-1.0", claims, vec![], vec![]).unwrap();
    with_expansion(base, factor)
}

fn compile_layered(cx: &Cx<'_>, gate: &CancelGate, text: &str, factors: [f64; 3]) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx, gate, |name| Ok(match NAMES.iter().position(|n| *n == name) {
        Some(i) => layer_pack(i, factors[i]),
        None => pack(name, 7800.0, 200_000.0, None),
    }))
}

fn section(source: &MaterialStringSource) -> &ConcentricKelvinVoigtSection {
    match source {
        MaterialStringSource::Concentric { section, .. } => section,
        MaterialStringSource::WithTorsion(t) => section(t.bending()),
        MaterialStringSource::WithPolarizations(p) => section(p.source()),
        _ => panic!("original bonded source must survive lowering"),
    }
}

fn expected(factors: [f64; 3], integral: f64) -> (f64, f64) {
    let strain = (0.5 - 0.49975) / 0.49975;
    let (mut tension, mut mu, mut ei) = (0.0, 0.0, 0.0);
    for i in 0..3 {
        let area = std::f64::consts::PI * (RADII[i+1].powi(2)-RADII[i].powi(2));
        let inertia = std::f64::consts::PI/4.0 * (RADII[i+1].powi(4)-RADII[i].powi(4));
        tension += YOUNG[i]*area*(strain-integral*factors[i]);
        mu += DENSITY[i]*area;
        ei += YOUNG[i]*inertia;
    }
    let k = std::f64::consts::PI/0.5;
    (tension, ((tension*k*k+ei*k.powi(4))/mu).sqrt()/std::f64::consts::TAU)
}

#[test]
fn independent_layer_expansion_changes_tension_frequency_and_pressure_not_other_physics() {
    scoped(|cx, gate| {
        let text = layered_input();
        let factors = [1.0, 0.5, 1.5];
        let changed = [1.0, 0.5, 2.0];
        let mut cold = compile_layered(cx,gate,&text.replace("state 300","state 290"),factors).unwrap();
        let mut hot = compile_layered(cx,gate,&text,factors).unwrap();
        let mut swapped = compile_layered(cx,gate,&text,changed).unwrap();
        assert_eq!(hot.spec.bytes,swapped.spec.bytes);
        assert_eq!(hot.spec.strings[0].thermal_reference_kelvin,Some(290.0));
        assert!(matches!(hot.spec.strings[0].prestress,StringPrestress::FixedExtension { .. }));
        assert_eq!(hot.catalog.len(),4);
        assert_ne!(hot.catalog[NAMES[2]].pack.content_hash(),swapped.catalog[NAMES[2]].pack.content_hash());
        for name in &NAMES[..2] {
            assert_eq!(hot.catalog[*name].pack.content_hash(),swapped.catalog[*name].pack.content_hash());
        }
        for (ready,factors,integral) in [(&cold,factors,0.0),(&hot,factors,0.0002),(&swapped,changed,0.0002)] {
            let bank = ready.performance.bank();
            let (tension,frequency) = expected(factors,integral);
            assert!((bank.cards()[0].tension_n-tension).abs()<1e-10);
            assert!((bank.cards()[0].fundamental_hz()-frequency).abs()<1e-10);
            let source = section(&bank.source().bank().strings()[0].source);
            for (layer,name) in source.elastic().layers().iter().zip(NAMES) {
                assert_eq!(layer.material.card_identity(),ready.catalog[name].card.content_hash());
                assert_eq!(layer.material.query_point(), &[("T".into(),ready.spec.temperature),("omega".into(),1000.0)]);
            }
            assert_eq!(bank.system().components()[0].modes().len(),11);
            assert_eq!(bank.system().components()[1].modes(),hot.performance.bank().system().components()[1].modes());
            assert_eq!(bank.source().bank().gas(),hot.performance.bank().source().bank().gas());
            assert!(bank.source().bank().material_cavity().is_some());
        }
        assert!(swapped.performance.bank().cards()[0].fundamental_hz()<hot.performance.bank().cards()[0].fundamental_hz());
        let a=pressure(&mut cold,1); let b=pressure(&mut hot,37); let c=pressure(&mut swapped,128);
        assert!(a.iter().chain(&b).chain(&c).all(|p|p.is_finite()));
        assert!(a.iter().zip(&b).any(|(a,b)|(a-b).abs()>1e-12));
        assert!(b.iter().zip(&c).any(|(b,c)|(b-c).abs()>1e-12));
        assert_eq!(b,pressure(&mut compile_layered(cx,gate,&text,factors).unwrap(),128));
        let frame=hot.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn zero_expansion_keeps_original_extension_bits_but_still_requires_every_source() {
    scoped(|cx,gate| {
        for (text,factors) in [(layered_input().replace("state 300","state 290"),[1.0,0.5,1.5]),
            (layered_input(),[0.0;3])] {
            let mut thermal=compile_layered(cx,gate,&text,factors).unwrap();
            let mut ordinary=compile_layered(cx,gate,&text.replace(SUPPORT,""),factors).unwrap();
            assert_eq!(thermal.performance.bank().cards()[0].tension_n,ordinary.performance.bank().cards()[0].tension_n);
            assert_eq!(pressure(&mut thermal,37),pressure(&mut ordinary,128));
        }
        for missing in NAMES {
            let load = |name: &str| Ok(if name == missing {
                pack(name,7800.0,200_000.0,None)
            } else if let Some(i)=NAMES.iter().position(|n|*n==name) {
                layer_pack(i,1.0)
            } else {pack(name,7800.0,200_000.0,None)});
            let cold=layered_input().replace("state 300","state 290");
            assert!(input::Spec::parse(cold.as_bytes()).unwrap().compile(cx,gate,load).is_err(),"missing {missing}");
            assert!(input::Spec::parse(cold.replace(SUPPORT,"").as_bytes()).unwrap().compile(cx,gate,load).is_ok());
        }
    });
}

#[test]
fn one_material_partitioned_into_layers_matches_the_homogeneous_sound() {
    scoped(|cx,gate| {
        let text=layered_input().replace("inner.fsmatpk","wire.fsmatpk").replace("outer.fsmatpk","wire.fsmatpk");
        let mut calls=0;
        let mut bonded=input::Spec::parse(text.as_bytes()).unwrap().compile(cx,gate,|name| {
            calls+=1;
            Ok(if name=="wire.fsmatpk" {wire(1.0)} else {pack(name,7800.0,200_000.0,None)})
        }).unwrap();
        let mut uniform=ready(cx,gate,&input(),1.0).unwrap();
        assert_eq!(calls,2,"one wire artifact and one body artifact, not one load per layer");
        assert_eq!(section(&bonded.performance.bank().source().bank().strings()[0].source).elastic().layers().len(),3);
        assert!((bonded.performance.bank().cards()[0].tension_n-uniform.performance.bank().cards()[0].tension_n).abs()<1e-10);
        assert!((bonded.performance.bank().cards()[0].fundamental_hz()-uniform.performance.bank().cards()[0].fundamental_hz()).abs()<1e-10);
        let a=pressure(&mut bonded,37); let b=pressure(&mut uniform,128);
        let energy=b.iter().map(|p|p*p).sum::<f64>();
        let error=a.iter().zip(&b).map(|(a,b)|(a-b)*(a-b)).sum::<f64>();
        assert!(energy>0.0 && error<=1e-16*energy,"partitioning a uniform source must not change its sound");
    });
}

#[test]
fn each_layer_and_the_taut_resultant_keep_the_original_thermal_limits() {
    scoped(|cx,gate| {
        let text=layered_input();
        // Positive resultant, but the contracting outside layer exceeds the
        // elastic strain limit. Checking only total force would accept this.
        assert!(expected([1.0,0.5,-5.0],0.0002).0>0.0);
        assert!(compile_layered(cx,gate,&text,[1.0,0.5,-5.0]).is_err());
        assert!(compile_layered(cx,gate,&text,[3.0;3]).is_err(),"thermal slack is not a retuning request");
        assert!(compile_layered(cx,gate,&text.replace("thermal-support 0 290","thermal-support 0 270"),[1.0;3]).is_err());
        assert!(compile_layered(cx,gate,&text.replace("0.49975 0.001","0.49975 2"),[0.0;3]).is_err());
        let duplicate=text.replace(SUPPORT,"thermal-supports 2\nthermal-support 0 290\nthermal-support 0 295\n");
        assert!(input::Spec::parse(duplicate.as_bytes()).is_err());
        assert!(input::Spec::parse(text.replace("extension 0.49975 0.001","tension 20").as_bytes()).is_err());
    });
}
