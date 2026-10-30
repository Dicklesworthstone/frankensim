//! Measured-section declarations enter the original native compiler and pressure path.
//! The test packs are manufactured adapters with published nominal numbers;
//! the documented user workflow compiles the original cited TSV, without a new importer.
use super::*;
use fs_couple::material_assembly::measured_string::{
    MeasuredStringProperties, ResolvedMeasuredStringSpecimen,
};

fn text() -> String {
    include_str!("../measured-cello-g.material-strings")
        .replace("48000 4096 0.01 128", "48000 128 0.01 128")
}
fn properties() -> MeasuredStringProperties {
    MeasuredStringProperties {
        length: "playing_length".into(),
        diameter: "diameter_dominant".into(),
        tension: "tension_dominant".into(),
        linear_density: "linear_density_dominant".into(),
        bending_stiffness: "bending_stiffness_dominant".into(),
    }
}
fn measured_pack(mass_factor: f64) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let artifact = fs_blake3::ContentHash([91; 32]);
    let provenance = Provenance {
        source: "manufactured measured-section frontend fixture; nominal Dominant values from the Lampis source pack".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured whole-string adapter fixture".into(),
            method: "explicit scalar section coefficients, not a bulk material reconstruction"
                .into(),
            artifact,
            caveats: "not a measurement, cello body, or reference-atmosphere calibration".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (requirement, value) in properties().requirements().unwrap().into_iter().zip([
        0.7,
        0.00115,
        116.3,
        0.0062 * mass_factor,
        0.000303,
    ]) {
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(requirement.name(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("temperature_known", 0.0, 0.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    // Even if a source contains mean Q, selecting it as EI cannot manufacture a loss law.
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("mean_quality_factor_first_30_modes_dominant", Dims::NONE),
            value: PropertyValue::Scalar {
                value: 1064.0,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("temperature_known", 0.0, 0.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let pack = NormalizedPack::new(
        "measured-frontend-fixture",
        "synthetic-test-v1",
        artifact,
        "manufactured test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, mass_factor: f64) -> Result<Ready> {
    spec.compile(cx, gate, |name| {
        Ok(if name == "cello-g-whole.fsmatpk" {
            measured_pack(mass_factor)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        })
    })
}
fn specimen(source: &MaterialStringSource) -> &ResolvedMeasuredStringSpecimen {
    match source {
        MaterialStringSource::Measured(s) => s,
        MaterialStringSource::WithPolarizations(p) => specimen(p.source()),
        _ => panic!("the original measured whole-string source must survive"),
    }
}

#[test]
fn measured_section_reaches_shared_pressure_without_fabricated_bulk_or_loss_properties() {
    scoped(|cx, gate| {
        let build = |factor| {
            ready(
                cx,
                gate,
                input::Spec::parse(text().as_bytes()).unwrap(),
                factor,
            )
            .unwrap()
        };
        let mut original = build(1.0);
        let mut changed = build(2.0);
        let bank = original.performance.bank();
        let source = specimen(&bank.source().bank().strings()[0].source);
        assert_eq!(source.properties(), &properties());
        assert_eq!(source.measurements().properties().len(), 5);
        for absent in [
            "density",
            "young_modulus",
            "shear_modulus",
            ETA,
            SHEAR_ETA,
            "mean_quality_factor_first_30_modes_dominant",
        ] {
            assert!(source.measurements().property(absent).is_none(), "{absent}");
        }
        assert_eq!(source.string().axial_stiffness_n, 0.0);
        assert_eq!(source.string().damping_ratio, 0.0);
        assert!(source.string().kelvin_voigt_bending.is_none());
        assert_eq!(bank.cards()[0].linear_density_kg_m, 0.0062);
        assert_eq!(bank.cards()[0].bending_stiffness_n_m2, 0.000303);
        assert_eq!(bank.cards()[0].tension_n, 116.3);
        assert_eq!(bank.system().components()[0].modes().len(), 8);
        // The real-source example must not add fluid inertia to apparent mass
        // inferred from acoustic frequencies with unknown reference atmosphere.
        assert!(bank.source().bank().cylinder_gas_loading().is_none());
        assert!(
            bank.cards()[0].zetas.iter().all(|zeta| *zeta > 0.0),
            "the declared lab gas still supplies the existing sampled resistance"
        );
        assert_eq!(original.spec.bytes, changed.spec.bytes);
        assert_eq!(
            bank.system().components()[1].modes(),
            changed.performance.bank().system().components()[1].modes()
        );
        assert!(
            changed.performance.bank().cards()[0].fundamental_hz()
                < bank.cards()[0].fundamental_hz()
        );
        let p = pressure(&mut original, 37);
        let q = pressure(&mut changed, 1);
        assert!(p.iter().any(|value| value.abs() > 1e-15));
        assert!(p.iter().zip(q).any(|(a, b)| (a - b).abs() > 1e-15));
        assert_eq!(p, pressure(&mut build(1.0), 128));
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, 48000, 0.01).unwrap();
        let (wav, summary) = build(1.0)
            .write(Cursor::new(Vec::<u8>::new()), gate)
            .unwrap();
        assert_eq!(wav.into_inner(), expected.0);
        assert_eq!(summary.samples, 128);
        let frame = original.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn unknown_measurement_temperature_is_independent_of_body_and_declared_lab_gas() {
    scoped(|cx, gate| {
        let a = ready(
            cx,
            gate,
            input::Spec::parse(text().as_bytes()).unwrap(),
            1.0,
        )
        .unwrap();
        let changed = text()
            .replace("state 293.15 1000", "state 310 1000")
            .replace("gas mixture 293.15 101325", "gas mixture 300 60000");
        let b = ready(
            cx,
            gate,
            input::Spec::parse(changed.as_bytes()).unwrap(),
            1.0,
        )
        .unwrap();
        let source = |ready: &Ready| {
            specimen(&ready.performance.bank().source().bank().strings()[0].source)
                .measurements()
                .identity()
        };
        assert_eq!(source(&a), source(&b));
        for ready in [&a, &b] {
            let bank = ready.performance.bank().source().bank();
            assert_eq!(
                specimen(&bank.strings()[0].source)
                    .measurements()
                    .query_point(),
                &[("temperature_known".into(), 0.0)]
            );
        }
        assert_ne!(
            a.performance
                .bank()
                .source()
                .bank()
                .material_plate()
                .unwrap()
                .specimen
                .material()
                .query_point(),
            b.performance
                .bank()
                .source()
                .bank()
                .material_plate()
                .unwrap()
                .specimen
                .material()
                .query_point()
        );
        assert_ne!(
            a.performance.bank().source().bank().gas(),
            b.performance.bank().source().bank().gas()
        );
        for replacement in ["0", "1 T 293.15", "1 temperature_known 1"] {
            let invalid = text().replace("1 temperature_known 0", replacement);
            let spec = input::Spec::parse(invalid.as_bytes()).unwrap();
            assert!(
                ready(cx, gate, spec, 1.0).is_err(),
                "accepted {replacement}"
            );
        }
    });
}

#[test]
fn measured_setup_and_explicit_elastic_scope_are_not_silently_replaced() {
    let original = text();
    for (from, to) in [
        (" no-intrinsic-loss", ""),
        ("no-intrinsic-loss", "infer-loss-from-Q"),
        (
            "playing_length diameter_dominant",
            "diameter_dominant diameter_dominant",
        ),
        (
            "1 temperature_known 0",
            "2 temperature_known 0 temperature_known 0",
        ),
        ("cavity none", "bending-loss kelvin-voigt\ncavity none"),
        (
            "cavity none",
            "bending-loss prony 1 modulus time\ncavity none",
        ),
        (
            "cavity none",
            "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2\nlongitudinal-memory 8 equilibrium 1 delta tau\ncavity none",
        ),
    ] {
        assert!(
            input::Spec::parse(original.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    scoped(|cx, gate| {
        for (from, to, field) in [
            ("0.7 0.000575", "0.71 0.000575", "length"),
            ("0.7 0.000575", "0.7 0.000576", "radius"),
            ("tension 116.3", "tension 117", "tension"),
        ] {
            let spec = input::Spec::parse(original.replace(from, to).as_bytes()).unwrap();
            let error = ready(cx, gate, spec, 1.0).err().unwrap().to_string();
            assert!(
                error.contains(field) && error.contains("measured setup"),
                "{error}"
            );
        }
        let invalid = original.replace(
            "bending_stiffness_dominant",
            "mean_quality_factor_first_30_modes_dominant",
        );
        assert!(
            ready(
                cx,
                gate,
                input::Spec::parse(invalid.as_bytes()).unwrap(),
                1.0
            )
            .is_err()
        );
        for (from, to) in [
            ("tension 116.3", "extension 0.699 0.01"),
            ("116.3 4 0 0.8", "116.3 4 1 0.8"),
            (
                "cavity none",
                "tension-loss lampis-2024-dominant-cello-g\ncavity none",
            ),
            (
                "cavity none",
                "geometric-stretching 0.2 64 100000 1e-11 1e-10\ncavity none",
            ),
            (
                "cavity none",
                "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2\ncavity none",
            ),
        ] {
            let spec = input::Spec::parse(original.replace(from, to).as_bytes()).unwrap();
            let error = ready(cx, gate, spec, 1.0).err().unwrap().to_string();
            assert!(error.contains("measured whole strings"), "{error}");
        }
    });
}
