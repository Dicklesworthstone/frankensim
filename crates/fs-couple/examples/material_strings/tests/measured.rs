//! Measured-section declarations enter the original native compiler and pressure path.
//! The test packs are manufactured adapters with published nominal numbers;
//! the documented user workflow compiles the original cited TSV, without a new importer.
use super::*;
use fs_couple::material_assembly::measured_string::{
    MeasuredStringLossProperties, MeasuredStringProperties, ResolvedMeasuredStringSpecimen,
};

const PREPARED: &str = "prepared-linear 64 67108864 200000000000";
const FRICTION: &str = "loss_coefficient_friction_fixture";
const BENDING: &str = "loss_coefficient_bending_fixture";
const AIR: &str = "loss_coefficient_air_fixture";
/// Eq. (8)'s friction and bending terms of Woodhouse 2004, p. 956, as zeta = eta / 2,
/// for this fixture's T, mu and EI with eta_F 2e-5 and eta_B 0.02, modes 1 to 4 of
/// each span. Worked from the printed form by a separate script, not by this crate.
const LAW_0_7_M: [f64; 4] = [
    1.052_421_479_942_407_79e-5,
    1.209_652_915_817_314_69e-5,
    1.471_595_347_701_406_16e-5,
    1.838_084_015_348_785_03e-5,
];
const LAW_0_875_M: [f64; 4] = [
    1.033_550_380_950_752_54e-5,
    1.134_188_004_107_344_04e-5,
    1.301_872_324_000_466_94e-5,
    1.536_535_810_214_646_96e-5,
];

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
    // Selected only by `loss-factors`. The air row's unit keeps it out of both roles.
    for (name, value, dims) in [
        (FRICTION, 2e-5, Dims::NONE),
        (BENDING, 0.02, Dims::NONE),
        (AIR, 1.2, Dims([0, 0, -1, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
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
fn a_speaking_length_source_binds_nut_to_bridge_and_a_test_span_source_never_does() {
    // A 0.875 m span with the bridge at station 0.8 puts nut to bridge at the source's measured 0.7 m.
    let longer = text().replace("0.7 0.000575", "0.875 0.000575");
    let speaking = |card: &str| card.replace("no-intrinsic-loss", "no-intrinsic-loss speaking-length");
    let parse = |card: &str| input::Spec::parse(card.as_bytes()).unwrap();
    scoped(|cx, gate| {
        // Planted negative: an unmarked (test-span) source refuses the longer span, although its bridge station
        // reproduces the measured length from the nut.
        let error = ready(cx, gate, parse(&longer), 1.0).err().unwrap().to_string();
        assert!(error.contains("length [m]") && error.contains("measured setup"), "{error}");
        // Marked speaking-length, the same source binds nut to bridge and keeps its coefficients on the whole span.
        let mut bound = ready(cx, gate, parse(&speaking(&longer)), 1.0).unwrap();
        {
            let source = specimen(&bound.performance.bank().source().bank().strings()[0].source);
            assert_eq!(source.string().length_m, 0.875);
            assert_eq!(source.speaking_length_m(), Some(0.7));
            assert_eq!(source.properties(), &properties());
            let card = &bound.performance.bank().cards()[0];
            assert_eq!(card.linear_density_kg_m, 0.0062);
            assert_eq!(card.bending_stiffness_n_m2, 0.000303);
            assert_eq!(card.tension_n, 116.3);
        }
        let p = pressure(&mut bound, 64);
        assert!(p.iter().all(|value| value.is_finite()) && p.iter().any(|value| value.abs() > 1e-15));
        // A speaking-length source whose span times bridge station misses the measured length refuses, including
        // the original 0.7 m span (0.56 m nut to bridge).
        for card in [speaking(&longer).replace("116.3 4 0 0.8 ", "116.3 4 0 0.81 "), speaking(&text())] {
            let error = ready(cx, gate, parse(&card), 1.0).err().unwrap().to_string();
            assert!(error.contains("speaking length") && error.contains("measured setup"), "{error}");
        }
        // test-span spelled out is the default meaning: the same string and the same pressure.
        let explicit = text().replace("no-intrinsic-loss", "no-intrinsic-loss test-span");
        let a = pressure(&mut ready(cx, gate, parse(&text()), 1.0).unwrap(), 64);
        let b = pressure(&mut ready(cx, gate, parse(&explicit), 1.0).unwrap(), 64);
        assert_eq!(a, b);
    });
    for token in ["scale", "speaking-length speaking-length", "test-span extra"] {
        let invalid = text().replace("no-intrinsic-loss", &format!("no-intrinsic-loss {token}"));
        assert!(input::Spec::parse(invalid.as_bytes()).is_err(), "accepted {token}");
    }
}

fn with_loss(card: &str) -> String {
    card.replace("no-intrinsic-loss", &format!("loss-factors {FRICTION} {BENDING}"))
}

#[test]
fn loss_factors_add_the_printed_law_to_every_card_and_leave_the_computed_drag_alone() {
    let parse = |card: &str| input::Spec::parse(card.as_bytes()).unwrap();
    let speaking = text()
        .replace("0.7 0.000575", "0.875 0.000575")
        .replace("no-intrinsic-loss", "no-intrinsic-loss speaking-length");
    scoped(|cx, gate| {
        for (card, law) in [(text(), LAW_0_7_M), (speaking, LAW_0_875_M)] {
            let plain = ready(cx, gate, parse(&card), 1.0).unwrap();
            let mut bound = ready(cx, gate, parse(&with_loss(&card)), 1.0).unwrap();
            let (a, b) = (plain.performance.bank(), bound.performance.bank());
            assert_eq!(a.cards().len(), b.cards().len());
            for (a, b) in a.cards().iter().zip(b.cards()) {
                assert!(!a.zetas.is_empty() && a.zetas.len() == b.zetas.len());
                // Polarizations repeat the same four modes.
                assert_eq!(a.zetas.len() % law.len(), 0);
                for (index, (a, b)) in a.zetas.iter().zip(&b.zetas).enumerate() {
                    let expected = law[index % law.len()];
                    assert!(*a > 0.0, "the computed drag is still there");
                    let added = b - a;
                    assert!(
                        (added / expected - 1.0).abs() < 1e-9,
                        "mode index {index}: added {added:e} vs {expected:e}"
                    );
                }
            }
            let source = specimen(&b.source().bank().strings()[0].source);
            assert_eq!(
                source.loss_properties(),
                Some(&MeasuredStringLossProperties { friction: FRICTION.into(), bending: BENDING.into() })
            );
            assert_eq!(source.measurements().properties().len(), 7);
            // Both rows stay on the receipts; eta_A is never read.
            assert!(source.measurements().property(FRICTION).is_some());
            assert!(source.measurements().property(BENDING).is_some());
            assert!(source.measurements().property(AIR).is_none());
            let p = pressure(&mut bound, 64);
            assert!(p.iter().all(|value| value.is_finite()) && p.iter().any(|value| value.abs() > 1e-15));
        }
        // Planted negative: the card without the selection reaches no rows and adds nothing.
        let plain = ready(cx, gate, parse(&text()), 1.0).unwrap();
        let source = specimen(&plain.performance.bank().source().bank().strings()[0].source);
        assert!(source.loss_properties().is_none());
        assert_eq!(source.string().tension_loss_factor, None);
        assert_eq!(source.string().bending_loss_factor, None);
        for absent in [FRICTION, BENDING, AIR] {
            assert!(source.measurements().property(absent).is_none(), "{absent}");
        }
    });
}

#[test]
fn a_loss_selection_refuses_every_row_it_cannot_own() {
    let original = text();
    for selection in [
        "loss-factors".to_owned(),
        format!("loss-factors {FRICTION}"),
        format!("loss-factors {FRICTION} {FRICTION}"),
        format!("loss-factors {FRICTION} {BENDING} extra"),
        format!("loss-factor {FRICTION} {BENDING}"),
    ] {
        let invalid = original.replace("no-intrinsic-loss", &selection);
        assert!(input::Spec::parse(invalid.as_bytes()).is_err(), "accepted {selection}");
    }
    // A loss factor does not open the door to the other bending laws.
    let invalid = with_loss(&original).replace(PREPARED, &format!("bending-loss kelvin-voigt\n{PREPARED}"));
    let error = input::Spec::parse(invalid.as_bytes()).err().unwrap().to_string();
    assert!(error.contains("measured whole strings do not supply an intrinsic bending-loss law"), "{error}");
    scoped(|cx, gate| {
        for invalid in [
            original.replace("no-intrinsic-loss", &format!("loss-factors {AIR} {BENDING}")),
            original.replace("no-intrinsic-loss", &format!("loss-factors {FRICTION} {AIR}")),
            original.replace("no-intrinsic-loss", &format!("loss-factors not_in_this_source {BENDING}")),
            original.replace(
                "no-intrinsic-loss",
                &format!("loss-factors {FRICTION} mean_quality_factor_first_30_modes_absent"),
            ),
        ] {
            let spec = input::Spec::parse(invalid.as_bytes()).unwrap();
            assert!(ready(cx, gate, spec, 1.0).is_err(), "accepted {invalid}");
        }
        let doubled = with_loss(&original)
            .replace(PREPARED, &format!("tension-loss lampis-2024-dominant-cello-g\n{PREPARED}"));
        let spec = input::Spec::parse(doubled.as_bytes()).unwrap();
        let error = ready(cx, gate, spec, 1.0).err().unwrap().to_string();
        assert!(error.contains("measured whole strings"), "{error}");
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
    ] {
        assert!(
            input::Spec::parse(original.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    for record in [
        "bending-loss kelvin-voigt",
        "bending-loss prony 1 modulus time",
    ] {
        let invalid = original.replace(PREPARED, &format!("{record}\n{PREPARED}"));
        let error = input::Spec::parse(invalid.as_bytes())
            .err()
            .unwrap()
            .to_string();
        assert!(
            error.contains("measured whole strings do not supply an intrinsic bending-loss law"),
            "{error}"
        );
    }
    // Unsupported history must reach source admission, not a prepared-owner refusal.
    let unprepared = original.replace(&format!("{PREPARED}\n"), "");
    let invalid = unprepared.replace("cavity none",
        "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2\nlongitudinal-memory 8 equilibrium 1 delta tau\ncavity none");
    let error = input::Spec::parse(invalid.as_bytes())
        .err()
        .unwrap()
        .to_string();
    assert!(
        error.contains("measured whole strings do not supply axial material history"),
        "{error}"
    );
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
        for invalid in [
            original.replace("tension 116.3", "extension 0.699 0.01"),
            original.replace("116.3 4 0 0.8", "116.3 4 1 0.8"),
            original.replace(
                PREPARED,
                &format!("tension-loss lampis-2024-dominant-cello-g\n{PREPARED}"),
            ),
            unprepared.replace(
                "cavity none",
                "geometric-stretching 0.2 64 100000 1e-11 1e-10\ncavity none",
            ),
            unprepared.replace(
                "cavity none",
                "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2\ncavity none",
            ),
        ] {
            let spec = input::Spec::parse(invalid.as_bytes()).unwrap();
            let error = ready(cx, gate, spec, 1.0).err().unwrap().to_string();
            assert!(error.contains("measured whole strings"), "{error}");
        }
    });
}
