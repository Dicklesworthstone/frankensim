//! G1/G3: sourced thermal eigenstrain reaches the same bowed string/body sound.
use super::*;
use fs_couple::string_specimen::ResolvedStringSpecimen;
use fs_material::state_point::{LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
    INVERSE_TEMPERATURE_DIMS};
use fs_qty::Temperature;

fn input() -> String {
    include_str!("../thermal-bow.material-strings").replace("48000 4096", "48000 256")
}

fn wire(factor: f64) -> NormalizedPack {
    let mut claims = pack("wire.fsmatpk", 7800.0, 200_000.0, None).claims().clone();
    let provenance = Provenance { source: "synthetic piecewise thermal expansion".into(),
        license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([8;32])) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic wire".into(), method: "manufactured alpha(T)".into(),
        artifact: fs_blake3::ContentHash([8;32]), caveats: "not measured material data".into(),
        provenance: provenance.clone(),
    }).unwrap();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ALPHA, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: [(280.0,1e-5),(290.0,1e-5),(295.0,3e-5),(300.0,1e-5),(320.0,1e-5)]
                .map(|(t,alpha)|(t,factor*alpha)).to_vec(), dims: INVERSE_TEMPERATURE_DIMS },
        validity: ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![observation], provenance,
    }).unwrap();
    let pack = NormalizedPack::new("thermal-wire", "synthetic-test-v1", fs_blake3::ContentHash([8;32]),
        "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, text: &str, factor: f64) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx, gate, |name| Ok(if name == "wire.fsmatpk" {
        wire(factor)
    } else { pack(name, 7800.0, 200_000.0, None) }))
}

fn specimen(source: &MaterialStringSource) -> &ResolvedStringSpecimen {
    match source {
        MaterialStringSource::Homogeneous(s) => s,
        MaterialStringSource::WithTorsion(t) => specimen(t.bending()),
        MaterialStringSource::WithPolarizations(p) => specimen(p.source()),
        _ => panic!("thermal fixture must retain its homogeneous source"),
    }
}

#[test]
fn integrated_expansion_changes_tension_modes_and_bowed_pressure_with_original_sources() {
    scoped(|cx, gate| {
        let mut cold = ready(cx, gate, &input().replace("state 300", "state 290"), 1.0).unwrap();
        let mut hot = ready(cx, gate, &input(), 1.0).unwrap();
        let mut swapped = ready(cx, gate, &input(), 2.0).unwrap();
        assert_eq!(hot.spec.bytes, swapped.spec.bytes, "only the material alpha(T) changes");
        assert_ne!(hot.catalog["wire.fsmatpk"].pack.content_hash(), swapped.catalog["wire.fsmatpk"].pack.content_hash());
        for (performance, free_strain) in [(&cold,0.0),(&hot,0.0002),(&swapped,0.0004)] {
            let bank = performance.performance.bank();
            let source = specimen(&bank.source().bank().strings()[0].source);
            let expansion = source.thermal_expansion().unwrap();
            // Integrate both triangular halves of the crossed 295 K knot.
            // Endpoint alpha*dT would give only 0.0001 for the hot specimen.
            assert!((expansion.free_linear_strain() - free_strain).abs() < 1e-17);
            assert_eq!(expansion.current().resolved().card_identity(), source.material().card_identity());
            assert_eq!(expansion.current().resolved().query_point(), source.material().query_point());
            let area = std::f64::consts::PI * 0.0003_f64.powi(2);
            let tension = 200e9 * area * ((0.5 - 0.49975) / 0.49975 - free_strain);
            assert!((bank.cards()[0].tension_n - tension).abs() < 1e-10);
            let k = std::f64::consts::PI / 0.5;
            let ei = 200e9 * area * 0.0003_f64.powi(2) / 4.0;
            let frequency = ((tension*k*k + ei*k.powi(4)) / (7800.0*area)).sqrt()
                / (2.0*std::f64::consts::PI);
            assert!((bank.cards()[0].fundamental_hz() - frequency).abs() < 1e-10);
            assert_eq!(source.string().length_m, 0.5);
            assert_eq!(source.radius_m(), 0.0003);
            assert_eq!(bank.system().components()[0].modes().len(), 11, "two bending planes and sourced twist retained");
            assert_eq!(bank.system().components()[1].modes(), cold.performance.bank().system().components()[1].modes());
            assert_eq!(bank.source().bank().gas(), cold.performance.bank().source().bank().gas());
        }
        assert!(hot.performance.bank().cards()[0].fundamental_hz() < cold.performance.bank().cards()[0].fundamental_hz());
        assert!(swapped.performance.bank().cards()[0].fundamental_hz() < hot.performance.bank().cards()[0].fundamental_hz());
        let a = pressure(&mut cold, 1); let b = pressure(&mut hot, 37); let c = pressure(&mut swapped, 128);
        assert!(a.iter().zip(&b).any(|(a,b)|(a-b).abs() > 1e-12));
        assert!(b.iter().zip(c).any(|(b,c)|(b-c).abs() > 1e-12));
        assert_eq!(b, pressure(&mut ready(cx, gate, &input(), 1.0).unwrap(), 128));
        let frame = hot.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn thermal_support_is_explicit_and_refuses_missing_paths_slack_or_unsupported_sections() {
    let block = "thermal-supports 1\nthermal-support 0 290\n";
    for text in [input().replace("thermal-supports 1", "thermal-supports 2"),
        input().replace("thermal-support 0 290", "thermal-support 1 290"),
        input().replace("thermal-support 0 290", "thermal-support 0 0"),
        input().replace("extension 0.49975 0.001", "tension 20"),
        input().replace("section homogeneous", "section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003")] {
        assert!(input::Spec::parse(text.as_bytes()).is_err());
    }
    scoped(|cx, gate| {
        let cold = input().replace("state 300", "state 290");
        let mut thermal = ready(cx, gate, &cold, 1.0).unwrap();
        let mut plain = ready(cx, gate, &cold.replace(block, ""), 1.0).unwrap();
        assert_eq!(pressure(&mut thermal, 37), pressure(&mut plain, 128), "zero thermal path preserves mechanics");
        let invalid_limit = cold.replace("extension 0.49975 0.001", "extension 0.49975 2");
        assert!(ready(cx,gate,&invalid_limit,1.0).is_err(), "zero thermal path must retain linear-strain admission");
        assert!(ready(cx,gate,&invalid_limit.replace(block,""),1.0).is_err());
        let no_expansion = |name: &str| Ok(pack(name,7800.0,200_000.0,None));
        assert!(input::Spec::parse(input().as_bytes()).unwrap().compile(cx,gate,no_expansion).is_err());
        assert!(input::Spec::parse(input().replace(block, "").as_bytes()).unwrap()
            .compile(cx,gate,no_expansion).is_ok(), "ordinary extension must not require alpha");
        assert!(ready(cx,gate,&input().replace("thermal-support 0 290", "thermal-support 0 270"),1.0).is_err());
        assert!(ready(cx,gate,&input(),3.0).is_err(), "thermal slack must not be clamped to a tuned tension");
    });
}
