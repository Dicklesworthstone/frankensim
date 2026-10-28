//! G0/G1/G3: source-bound constant-mass density uses finite free expansion.
//! Authored fixtures test algebra, not experimental material calibration.
#![deny(unsafe_code)]

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::thermal::transport::{
    ThermalTransportError,
    expansion::{
        FreeExpansionBasis, FreeExpansionSourceQuery, FreeExpansionStrainSources,
        PreparedFreeExpansionDensity,
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const FRAME: ContentHash = ContentHash([7; 32]);

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", absolute(), t)
        .unwrap()
}
fn domain(low: f64, high: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("temperature", absolute(), low, high)
}
fn keys() -> [PropertyKey; 4] {
    [
        PropertyKey::new("density", RHO),
        PropertyKey::new("engineering-free-strain-x", Dims::NONE),
        PropertyKey::new("engineering-free-strain-y", Dims::NONE),
        PropertyKey::new("engineering-free-strain-z", Dims::NONE),
    ]
}
fn curve(key: &PropertyKey, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        key: key.clone(),
        validity: domain(knots[0].0, knots[knots.len() - 1].0),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: key.quantity().dims(),
        },
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.02,
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "authored expansion regression fixture; not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn density(t: f64, value: f64) -> PropertyClaim {
    let mut claim = curve(&keys()[0], vec![(t, value), (t + 1.0, value)]);
    claim.value = PropertyValue::Scalar { value, dims: RHO };
    claim.validity = domain(t, t);
    claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claim
}
fn material(set: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "declared expansion test model".into(),
            phase: "solid".into(),
            process: "explicit reference density and free-expansion assumptions".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap()
}
fn card(claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    material(set)
}
fn query<'a>(
    keys: &'a [PropertyKey; 4],
    lower: &'a QueryPoint,
    upper: &'a QueryPoint,
    reference: &'a QueryPoint,
) -> FreeExpansionSourceQuery<'a> {
    FreeExpansionSourceQuery {
        density: &keys[0],
        density_selection: POLICY,
        reference,
        strains: FreeExpansionStrainSources::Isotropic {
            strain: &keys[1],
            selection: POLICY,
        },
        temperature_axis: "temperature",
        lower,
        upper,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_g3_isotropic_normalization_preserves_reference_residual_and_original_evidence() {
    let keys = keys();
    let card = card([
        density(310.0, 998.25),
        curve(&keys[1], vec![(280.0, -0.04), (360.0, 0.12)]),
    ]);
    let (lower, upper, reference) = (point(280.0), point(360.0), point(310.0));
    let q = query(&keys, &lower, &upper, &reference);
    let prepared = PreparedFreeExpansionDensity::try_new(&card, q).unwrap();
    assert_eq!(
        prepared.temperatures_k().collect::<Vec<_>>(),
        [280.0, 310.0, 360.0]
    );
    assert_eq!(prepared.material_card_identity(), card.content_hash());
    assert_eq!(prepared.query_endpoints(), [&lower, &upper]);
    assert_eq!(prepared.reference_query_point(), &reference);
    assert_eq!(prepared.temperature_axis(), "temperature");
    assert_eq!(prepared.basis(), FreeExpansionBasis::Isotropic);
    let reference_answer = card
        .claims()
        .query_typed(&keys[0], &reference, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(prepared.reference_density_answer(), &reference_answer);
    let original = card
        .claims()
        .query_line_typed(&keys[1], "temperature", &lower, &upper, POLICY)
        .unwrap();
    assert_eq!(prepared.source_lines(), &[original]);
    for (_, answer) in &prepared.source_lines()[0].knots {
        card.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::RelativeHalfWidth {
                fraction: 0.02,
                confidence: 0.95
            }
        );
    }
    card.claims()
        .verify_receipt(&prepared.reference_density_answer().receipt)
        .unwrap();
    let at_reference = prepared.at_temperature(310.0).unwrap();
    assert_eq!(at_reference.density_kg_m3().to_bits(), 998.25_f64.to_bits());
    assert_eq!(at_reference.principal_stretches(), [1.0; 3]);
    assert_eq!(at_reference.relative_volume(), 1.0);
    // The reference is inside an original source segment and its epsilon is
    // 0.02, not zero. Replacing the original datum would change the model.
    for t in [280.0, 285.0, 295.0, 310.0, 320.0, 345.0, 360.0] {
        let stretch = (1.0 - 0.04 + 0.002 * (t - 280.0)) / 1.02;
        let sample = prepared.at_temperature(t).unwrap();
        close(sample.relative_volume(), stretch.powi(3));
        close(sample.density_kg_m3(), 998.25 / stretch.powi(3));
        assert_eq!(sample.path_identity(), prepared.identity());
        assert_eq!(sample.temperature_k(), t);
    }
    let mid = prepared.at_temperature(295.0).unwrap();
    let low = prepared.at_temperature(280.0).unwrap();
    let high = prepared.at_temperature(310.0).unwrap();
    assert!(
        (mid.relative_volume() - 0.5 * (low.relative_volume() + high.relative_volume())).abs()
            > 0.001
    );
    assert!((mid.density_kg_m3() - 0.5 * (low.density_kg_m3() + high.density_kg_m3())).abs() > 1.0);
    assert_eq!(
        PreparedFreeExpansionDensity::try_new(&card, q)
            .unwrap()
            .identity(),
        prepared.identity()
    );
    let density_pin = ClaimSelection::Pinned(reference_answer.receipt.selected);
    let strain_pin = ClaimSelection::Pinned(prepared.source_lines()[0].knots[0].1.receipt.selected);
    let pinned = PreparedFreeExpansionDensity::try_new(
        &card,
        FreeExpansionSourceQuery {
            density_selection: density_pin,
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &keys[1],
                selection: strain_pin,
            },
            ..q
        },
    )
    .unwrap();
    assert_ne!(pinned.identity(), prepared.identity());
    assert_eq!(
        pinned.at_temperature(310.0).unwrap().density_kg_m3(),
        998.25
    );
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &card,
            FreeExpansionSourceQuery {
                strains: FreeExpansionStrainSources::Isotropic {
                    strain: &keys[1],
                    selection: density_pin
                },
                ..q
            }
        )
        .is_err()
    );
}

#[test]
fn g1_principal_volume_uses_union_inputs_and_remains_safe_under_axis_permutations() {
    let keys = keys();
    let card = card([
        density(320.0, 1000.0),
        curve(&keys[1], vec![(280.0, -0.1), (300.0, 0.0), (360.0, 0.3)]),
        curve(&keys[2], vec![(280.0, -0.2), (340.0, 0.1), (360.0, 0.5)]),
        curve(&keys[3], vec![(280.0, -0.05), (360.0, 0.1)]),
    ]);
    let (lower, upper, reference) = (point(280.0), point(360.0), point(320.0));
    let q = FreeExpansionSourceQuery {
        strains: FreeExpansionStrainSources::FixedPrincipal {
            frame: FRAME,
            strains: [&keys[1], &keys[2], &keys[3]],
            selections: [POLICY; 3],
        },
        ..query(&keys, &lower, &upper, &reference)
    };
    let prepared = PreparedFreeExpansionDensity::try_new(&card, q).unwrap();
    assert_eq!(
        prepared.basis(),
        FreeExpansionBasis::FixedPrincipal { frame: FRAME }
    );
    assert_eq!(
        prepared.temperatures_k().collect::<Vec<_>>(),
        [280.0, 300.0, 320.0, 340.0, 360.0]
    );
    let source = |i, t| {
        card.claims()
            .query_typed(&keys[i], &point(t), SelectionPolicy::SingleClaimOnly)
            .unwrap()
            .evidence
            .value
            .value
    };
    for t in [280.0, 295.0, 310.0, 320.0, 330.0, 350.0, 360.0] {
        let stretches: [f64; 3] =
            std::array::from_fn(|i| (1.0 + source(i + 1, t)) / (1.0 + source(i + 1, 320.0)));
        let volume = stretches.into_iter().product::<f64>();
        let sample = prepared.at_temperature(t).unwrap();
        for (actual, expected) in sample.principal_stretches().into_iter().zip(stretches) {
            close(actual, expected);
        }
        close(sample.relative_volume(), volume);
        close(sample.density_kg_m3(), 1000.0 / volume);
    }
    assert!(prepared.at_temperature(280.0).unwrap().density_kg_m3() > 1000.0);
    let hot = prepared.at_temperature(360.0).unwrap();
    assert!(hot.density_kg_m3() < 1000.0);
    let trace_approximation = 1.0
        + hot
            .principal_stretches()
            .into_iter()
            .map(|s| s - 1.0)
            .sum::<f64>();
    assert!((hot.relative_volume() - trace_approximation).abs() > 0.05);
    let other_frame = PreparedFreeExpansionDensity::try_new(
        &card,
        FreeExpansionSourceQuery {
            strains: FreeExpansionStrainSources::FixedPrincipal {
                frame: ContentHash([8; 32]),
                strains: [&keys[1], &keys[2], &keys[3]],
                selections: [POLICY; 3],
            },
            ..q
        },
    )
    .unwrap();
    assert_ne!(prepared.identity(), other_frame.identity());
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &card,
            FreeExpansionSourceQuery {
                strains: FreeExpansionStrainSources::FixedPrincipal {
                    frame: ContentHash([0; 32]),
                    strains: [&keys[1], &keys[2], &keys[3]],
                    selections: [POLICY; 3]
                },
                ..q
            }
        )
        .is_err()
    );

    // All original lengths are positive and representable. At 360 K, normalized
    // stretches [1e200,1e200,1e-300] have finite J despite naive ordered overflow.
    let extreme = self::card([
        density(280.0, 1.0),
        curve(&keys[1], vec![(280.0, 0.0), (360.0, 1e200)]),
        curve(&keys[2], vec![(280.0, 0.0), (360.0, 1e200)]),
        curve(&keys[3], vec![(280.0, 1e300), (360.0, 0.0)]),
    ]);
    for order in [[1, 2, 3], [3, 1, 2]] {
        let path = PreparedFreeExpansionDensity::try_new(
            &extreme,
            FreeExpansionSourceQuery {
                reference: &lower,
                strains: FreeExpansionStrainSources::FixedPrincipal {
                    frame: FRAME,
                    strains: order.map(|i| &keys[i]),
                    selections: [POLICY; 3],
                },
                ..q
            },
        )
        .unwrap();
        let before = path.clone();
        close(path.at_temperature(360.0).unwrap().relative_volume(), 1e100);
        close(path.at_temperature(360.0).unwrap().density_kg_m3(), 1e-100);
        // Here the actual interior J is unrepresentable; the prepared model is
        // still usable afterward and all failures are transactional.
        assert!(matches!(
            path.at_temperature(320.0),
            Err(ThermalTransportError::Invalid(_))
        ));
        for invalid in [279.0, 361.0, f64::NAN, f64::INFINITY] {
            assert!(path.at_temperature(invalid).is_err());
        }
        assert_eq!(path.at_temperature(280.0).unwrap().density_kg_m3(), 1.0);
        assert_eq!(path, before);
    }
}

#[test]
fn g0_unsupported_source_spans_and_reference_context_cannot_be_promoted_to_expansion() {
    let keys = keys();
    let (lower, upper, reference) = (point(280.0), point(360.0), point(310.0));
    let q = query(&keys, &lower, &upper, &reference);
    let strain = || curve(&keys[1], vec![(280.0, -0.04), (360.0, 0.12)]);
    let supported = card([density(310.0, 998.25), strain()]);
    let prepare = |claims| PreparedFreeExpansionDensity::try_new(&card(claims), q);
    let mut discrete = strain();
    discrete.interpolation = InterpolationPolicy::TabulatedOnly;
    assert!(matches!(
        prepare(vec![density(310.0, 998.25), discrete]),
        Err(ThermalTransportError::Source { .. })
    ));
    let gapped = vec![
        density(310.0, 998.25),
        curve(&keys[1], vec![(280.0, -0.04), (300.0, 0.0)]),
        curve(&keys[1], vec![(340.0, 0.08), (360.0, 0.12)]),
    ];
    assert!(matches!(
        prepare(gapped),
        Err(ThermalTransportError::Source { .. })
    ));
    let conflict = curve(&keys[1], vec![(280.0, -0.03), (360.0, 0.13)]);
    assert!(prepare(vec![density(310.0, 998.25), strain(), conflict]).is_err());
    let mut unknown = strain();
    unknown.value = PropertyValue::Scalar {
        value: 0.0,
        dims: Dims::NONE,
    };
    unknown.validity = ValidityDomain::unconstrained();
    unknown.interpolation = InterpolationPolicy::ConstantWithinValidity;
    assert!(matches!(
        prepare(vec![density(310.0, 998.25), unknown]),
        Err(ThermalTransportError::Invalid(_))
    ));
    for value in [-1.0, -1.1] {
        assert!(
            prepare(vec![
                density(310.0, 998.25),
                curve(&keys[1], vec![(280.0, value), (360.0, 0.12)])
            ])
            .is_err()
        );
    }
    for value in [0.0, -1.0] {
        assert!(prepare(vec![density(310.0, value), strain()]).is_err());
    }
    let inverse_temperature =
        PropertyKey::new("linear-expansion-coefficient", Dims([0, 0, 0, -1, 0, 0]));
    for wrong in [&keys[0], &inverse_temperature] {
        assert!(
            PreparedFreeExpansionDensity::try_new(
                &supported,
                FreeExpansionSourceQuery {
                    strains: FreeExpansionStrainSources::Isotropic {
                        strain: wrong,
                        selection: POLICY
                    },
                    ..q
                }
            )
            .is_err()
        );
    }
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            FreeExpansionSourceQuery {
                density: &keys[1],
                ..q
            }
        )
        .is_err()
    );
    for t in [279.0, 361.0, 300.0] {
        // 300 K is inside the strain span but the density does not support it.
        assert!(
            PreparedFreeExpansionDensity::try_new(
                &supported,
                FreeExpansionSourceQuery {
                    reference: &point(t),
                    ..q
                }
            )
            .is_err()
        );
    }
    for (lo, hi) in [(&upper, &lower), (&lower, &lower)] {
        assert!(
            PreparedFreeExpansionDensity::try_new(
                &supported,
                FreeExpansionSourceQuery {
                    lower: lo,
                    upper: hi,
                    ..q
                }
            )
            .is_err()
        );
    }
    let difference = |t| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::TemperatureDifference,
                    ValueForm::Static,
                )),
                t,
            )
            .unwrap()
    };
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            query(
                &keys,
                &difference(280.0),
                &difference(360.0),
                &difference(310.0)
            )
        )
        .is_err()
    );
    let untyped = |t| QueryPoint::new().with("temperature", t).unwrap();
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            query(&keys, &untyped(280.0), &untyped(360.0), &untyped(310.0))
        )
        .is_err()
    );
    let conditioned = |t, value| {
        point(t)
            .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    };
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            query(
                &keys,
                &conditioned(280.0, 0.0),
                &conditioned(360.0, 0.0),
                &conditioned(310.0, 1.0)
            )
        )
        .is_err()
    );
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            query(
                &keys,
                &conditioned(280.0, 0.0),
                &conditioned(360.0, 1.0),
                &conditioned(310.0, 0.0)
            )
        )
        .is_err()
    );
    let dynamic = |t| {
        point(t)
            .with_quantity(
                "test-rate",
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                1.0,
            )
            .unwrap()
    };
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &supported,
            query(&keys, &dynamic(280.0), &dynamic(360.0), &dynamic(310.0))
        )
        .is_err()
    );
}

#[test]
fn g1_g3_compiled_nist_relative_length_uses_an_explicitly_assumed_reference_density() {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-expansion-density-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace.join("data/matdb/seed-v1/stainless-304-nist-cryogenic/manifest.tsv"))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace)
        .current_dir(workspace)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert!(pack.claims().claims_for("density").is_empty());
    let mut set = pack.claims().clone();
    // This is a caller's engineering assumption for an integration fixture,
    // not a NIST density observation or a matched specimen claim. The original
    // seed pack remains unmodified, and density is supported at 293 K only.
    let density_key = PropertyKey::new("assumed-reference-density", RHO);
    let mut assumption = density(293.0, 8000.0);
    assumption.key = density_key.clone();
    assumption.provenance.source = "AUTHORED ASSUMPTION: 8000 kg/m3 at 293 K for a free-expansion integration fixture; not NIST data or matched specimens".into();
    assumption.uncertainty = UncertaintyModel::Unstated;
    set.insert_claim(assumption).unwrap();
    let card = material(set);
    let context = |t| {
        point(t)
            .with_quantity(
                "source-pressure-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap()
    };
    let strain_key = PropertyKey::new("linear-expansion-relative-to-293k", Dims::NONE);
    let (lower, upper, reference) = (context(77.0), context(300.0), context(293.0));
    let q = FreeExpansionSourceQuery {
        density: &density_key,
        density_selection: POLICY,
        reference: &reference,
        strains: FreeExpansionStrainSources::Isotropic {
            strain: &strain_key,
            selection: POLICY,
        },
        temperature_axis: "temperature",
        lower: &lower,
        upper: &upper,
    };
    let prepared = PreparedFreeExpansionDensity::try_new(&card, q).unwrap();
    assert_eq!(prepared.source_lines()[0].knots.len(), 11);
    assert_eq!(
        prepared.at_temperature(293.0).unwrap().density_kg_m3(),
        8000.0
    );
    let residual = prepared.source_lines()[0]
        .knots
        .iter()
        .find(|(t, _)| *t == 293.0)
        .unwrap()
        .1
        .evidence
        .value
        .value;
    assert_eq!(residual, 7.4646191727e-7);
    assert_ne!(residual, 0.0);
    for t in [77.0, 87.0, 150.0, 205.0, 293.0, 299.0, 300.0] {
        let original = card
            .claims()
            .query_typed(&strain_key, &context(t), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        let stretch = (1.0 + original.evidence.value.value) / (1.0 + residual);
        close(
            prepared.at_temperature(t).unwrap().density_kg_m3(),
            8000.0 / stretch.powi(3),
        );
    }
    for (_, answer) in &prepared.source_lines()[0].knots {
        assert!(answer.receipt.observation_backed);
        assert!(!answer.receipt.source_hashes.is_empty());
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        card.claims().verify_receipt(&answer.receipt).unwrap();
    }
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &card,
            FreeExpansionSourceQuery {
                lower: &context(76.0),
                ..q
            }
        )
        .is_err()
    );
    let wrong_context = reference
        .clone()
        .with("source-pressure-known", 1.0)
        .unwrap();
    assert!(
        PreparedFreeExpansionDensity::try_new(
            &card,
            FreeExpansionSourceQuery {
                reference: &wrong_context,
                ..q
            }
        )
        .is_err()
    );
}
