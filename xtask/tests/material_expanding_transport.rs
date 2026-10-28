//! G0/G1/G3: compose original thermal sources with derived expansion density.
//! Numerical fixtures and declared reference assumptions are not measured data.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::thermal::transport::{
    ThermalTransportError,
    expanding::{ExpandingThermalSourceQuery, PreparedExpandingThermalTransport},
    expansion::{
        FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
}
fn keys() -> [PropertyKey; 4] {
    [
        PropertyKey::new("assumed-reference-density", RHO),
        PropertyKey::new("engineering-free-strain", Dims::NONE),
        PropertyKey::new("specific-heat-capacity", CP),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            semantic(QuantityKind::ThermalConductivity),
        ),
    ]
}
fn curve(key: &PropertyKey, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        key: key.clone(),
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            knots[0].0,
            knots[knots.len() - 1].0,
        ),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: key.quantity().dims(),
        },
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.01,
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "authored expanding-transport regression; not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn density(t: f64, value: f64) -> PropertyClaim {
    let mut claim = curve(&keys()[0], vec![(t, value), (t + 1.0, value)]);
    claim.value = PropertyValue::Scalar { value, dims: RHO };
    claim.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        t,
        t,
    );
    claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claim
}
fn material(set: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "declared expanding thermal test model".into(),
            phase: "solid".into(),
            process: "explicit density-reference and free-expansion assumptions".into(),
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
fn claims() -> Vec<PropertyClaim> {
    let k = keys();
    vec![
        density(310.0, 1000.0),
        curve(&k[1], vec![(280.0, -0.1), (310.0, 0.05), (360.0, 0.3)]),
        curve(&k[2], vec![(280.0, 400.0), (320.0, 800.0), (360.0, 900.0)]),
        curve(&k[3], vec![(280.0, 1.0), (340.0, 3.0), (360.0, 6.0)]),
    ]
}
fn expansion(
    card: &MaterialCard,
    k: &[PropertyKey; 4],
    lower: &QueryPoint,
    upper: &QueryPoint,
    reference: &QueryPoint,
) -> PreparedFreeExpansionDensity {
    PreparedFreeExpansionDensity::try_new(
        card,
        FreeExpansionSourceQuery {
            density: &k[0],
            density_selection: POLICY,
            reference,
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &k[1],
                selection: POLICY,
            },
            temperature_axis: "temperature",
            lower,
            upper,
        },
    )
    .unwrap()
}
fn query<'a>(
    k: &'a [PropertyKey; 4],
    lower: &'a QueryPoint,
    upper: &'a QueryPoint,
) -> ExpandingThermalSourceQuery<'a> {
    ExpandingThermalSourceQuery {
        specific_heat: &k[2],
        conductivity: &k[3],
        temperature_axis: "temperature",
        lower,
        upper,
        selections: [POLICY; 2],
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_g3_original_density_and_thermal_sources_compose_without_interpolating_products() {
    let k = keys();
    let card = card(claims());
    let (lower, upper, reference) = (point(280.0), point(360.0), point(310.0));
    let density = expansion(&card, &k, &lower, &upper, &reference);
    let q = query(&k, &lower, &upper);
    let prepared = PreparedExpandingThermalTransport::try_new(&card, density.clone(), q).unwrap();
    assert_eq!(prepared.density_model(), &density);
    assert_eq!(prepared.material_card_identity(), card.content_hash());
    assert_eq!(prepared.query_endpoints(), [&lower, &upper]);
    assert_eq!(prepared.temperature_axis(), "temperature");
    assert_eq!(prepared.temperature_span_k(), [280.0, 360.0]);
    assert_eq!(
        prepared.temperatures_k().collect::<Vec<_>>(),
        [280.0, 320.0, 340.0, 360.0]
    );
    assert_eq!(
        density.temperatures_k().collect::<Vec<_>>(),
        [280.0, 310.0, 360.0]
    );
    for (key, line) in k[2..].iter().zip(prepared.source_lines()) {
        assert_eq!(
            *line,
            card.claims()
                .query_line_typed(key, "temperature", &lower, &upper, POLICY)
                .unwrap()
        );
        for (_, answer) in &line.knots {
            card.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.95
                }
            );
        }
    }
    for t in [280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 345.0, 360.0] {
        let source = |index| {
            card.claims()
                .query_typed(&k[index], &point(t), SelectionPolicy::SingleClaimOnly)
                .unwrap()
                .evidence
                .value
                .value
        };
        let rho = 1000.0 / ((1.0 + source(1)) / 1.05).powi(3);
        let (cp, conductivity) = (source(2), source(3));
        let sample = prepared.at_temperature(t).unwrap();
        close(sample.density_kg_m3(), rho);
        close(sample.specific_heat_j_kg_k(), cp);
        close(sample.conductivity_w_m_k(), conductivity);
        close(sample.volumetric_heat_capacity_j_m3_k(), rho * cp);
        close(sample.diffusivity_m2_s(), conductivity / (rho * cp));
        close(
            sample.effusivity_j_m2_k_sqrt_s(),
            (conductivity * rho * cp).sqrt(),
        );
        assert_eq!(sample.path_identity(), prepared.identity());
        assert_eq!(sample.temperature_k(), t);
    }
    assert_eq!(
        prepared.at_temperature(310.0).unwrap().density_kg_m3(),
        1000.0
    );
    let (a, m, b) = (
        prepared.at_temperature(320.0).unwrap(),
        prepared.at_temperature(330.0).unwrap(),
        prepared.at_temperature(340.0).unwrap(),
    );
    assert!((m.density_kg_m3() - 0.5 * (a.density_kg_m3() + b.density_kg_m3())).abs() > 1.0);
    assert!(
        (m.volumetric_heat_capacity_j_m3_k()
            - 0.5 * (a.volumetric_heat_capacity_j_m3_k() + b.volumetric_heat_capacity_j_m3_k()))
        .abs()
            > 100.0
    );
    let selections = prepared
        .source_lines()
        .each_ref()
        .map(|line| ClaimSelection::Pinned(line.knots[0].1.receipt.selected));
    let pinned = PreparedExpandingThermalTransport::try_new(
        &card,
        density.clone(),
        ExpandingThermalSourceQuery { selections, ..q },
    )
    .unwrap();
    assert_ne!(pinned.identity(), prepared.identity());
    assert_eq!(pinned.density_model(), &density);
    assert_eq!(
        pinned.at_temperature(330.0).unwrap().density_kg_m3(),
        m.density_kg_m3()
    );
    assert_eq!(
        PreparedExpandingThermalTransport::try_new(&card, density, q)
            .unwrap()
            .identity(),
        prepared.identity()
    );
}

#[test]
fn g0_composition_requires_identical_identity_context_complete_support_and_finite_results() {
    let k = keys();
    let card = card(claims());
    let (lower, upper, reference) = (point(280.0), point(360.0), point(310.0));
    let density_model = expansion(&card, &k, &lower, &upper, &reference);
    let q = query(&k, &lower, &upper);
    let other_card = MaterialCard::supersede(&card, card.claims().clone(), Vec::new()).unwrap();
    assert!(
        PreparedExpandingThermalTransport::try_new(&other_card, density_model.clone(), q).is_err()
    );
    let changed = lower
        .clone()
        .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    let dimensional = lower
        .clone()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            280.0,
        )
        .unwrap();
    for endpoint in [point(290.0), changed, dimensional] {
        assert!(
            PreparedExpandingThermalTransport::try_new(
                &card,
                density_model.clone(),
                ExpandingThermalSourceQuery {
                    lower: &endpoint,
                    ..q
                }
            )
            .is_err()
        );
    }
    assert!(
        PreparedExpandingThermalTransport::try_new(
            &card,
            density_model.clone(),
            ExpandingThermalSourceQuery {
                temperature_axis: "T",
                ..q
            }
        )
        .is_err()
    );
    assert!(
        PreparedExpandingThermalTransport::try_new(
            &card,
            density_model.clone(),
            ExpandingThermalSourceQuery {
                specific_heat: &k[0],
                ..q
            }
        )
        .is_err()
    );
    assert!(
        PreparedExpandingThermalTransport::try_new(
            &card,
            density_model.clone(),
            ExpandingThermalSourceQuery {
                conductivity: &k[2],
                ..q
            }
        )
        .is_err()
    );
    let cp_pin = card.claims().claims_for(k[2].name())[0].0;
    assert!(
        PreparedExpandingThermalTransport::try_new(
            &card,
            density_model.clone(),
            ExpandingThermalSourceQuery {
                selections: [POLICY, ClaimSelection::Pinned(cp_pin)],
                ..q
            }
        )
        .is_err()
    );
    let mut bad_sources = Vec::new();
    let mut discrete = claims();
    discrete[2].interpolation = InterpolationPolicy::TabulatedOnly;
    bad_sources.push(discrete);
    let mut short = claims();
    short[2] = curve(&k[2], vec![(280.0, 400.0), (320.0, 800.0)]);
    bad_sources.push(short);
    let mut unbounded = claims();
    unbounded[2].value = PropertyValue::Scalar {
        value: 500.0,
        dims: CP,
    };
    unbounded[2].validity = ValidityDomain::unconstrained();
    unbounded[2].interpolation = InterpolationPolicy::ConstantWithinValidity;
    bad_sources.push(unbounded);
    let mut nonpositive = claims();
    nonpositive[2] = curve(&k[2], vec![(280.0, 400.0), (320.0, 0.0), (360.0, 900.0)]);
    bad_sources.push(nonpositive);
    for bad in bad_sources {
        let unsupported = self::card(bad);
        let density = expansion(&unsupported, &k, &lower, &upper, &reference);
        assert!(PreparedExpandingThermalTransport::try_new(&unsupported, density, q).is_err());
    }
    // Cp/k alone have only endpoint knots. A density-only knot at 310 K makes
    // rho*Cp overflow, so checking only the thermal grid would miss this refusal.
    let middle_overflow = self::card([
        density(310.0, 1e300),
        curve(&k[1], vec![(280.0, 1e100), (310.0, 0.0), (360.0, 1e100)]),
        curve(&k[2], vec![(280.0, 1e100), (360.0, 1e100)]),
        curve(&k[3], vec![(280.0, 1.0), (360.0, 1.0)]),
    ]);
    let density = expansion(&middle_overflow, &k, &lower, &upper, &reference);
    assert!(matches!(
        PreparedExpandingThermalTransport::try_new(&middle_overflow, density, q),
        Err(ThermalTransportError::Invalid(_))
    ));

    // All input knots have representable outputs; an interior product need
    // not. Synthetic extreme kelvin values exercise the arithmetic boundary.
    let (lo, hi) = (point(1e-200), point(1.0));
    let extreme = self::card([
        self::density(1e-200, 1e200),
        curve(&k[1], vec![(1e-200, 0.0), (1.0, 1e100)]),
        curve(&k[2], vec![(1e-200, 1e-100), (1.0, 1e300)]),
        curve(&k[3], vec![(1e-200, 1.0), (1.0, 1.0)]),
    ]);
    let density = expansion(&extreme, &k, &lo, &hi, &lo);
    let prepared =
        PreparedExpandingThermalTransport::try_new(&extreme, density, query(&k, &lo, &hi)).unwrap();
    let before = prepared.clone();
    assert!(matches!(
        prepared.at_temperature(1e-100),
        Err(ThermalTransportError::Invalid(_))
    ));
    for t in [0.0, 2.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.at_temperature(t).is_err());
    }
    close(
        prepared
            .at_temperature(1e-200)
            .unwrap()
            .volumetric_heat_capacity_j_m3_k(),
        1e100,
    );
    close(
        prepared
            .at_temperature(1.0)
            .unwrap()
            .volumetric_heat_capacity_j_m3_k(),
        1e200,
    );
    assert_eq!(prepared, before);
}

#[test]
fn g1_g3_compiled_nist_sources_supply_transport_with_an_explicit_reference_density_assumption() {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-expanding-transport-{}-{nonce}.fsmatpk",
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
    // Caller-authored assumption, not a NIST density value or matched specimen.
    // The original pack remains unchanged and the added density is point-only.
    let mut assumption = density(293.0, 8000.0);
    assumption.provenance.source = "AUTHORED ASSUMPTION: 8000 kg/m3 at 293 K for expanding transport; not NIST data or matched specimens".into();
    assumption.uncertainty = UncertaintyModel::Unstated;
    set.insert_claim(assumption).unwrap();
    let card = material(set);
    let mut k = keys();
    k[1] = PropertyKey::new("linear-expansion-relative-to-293k", Dims::NONE);
    let context = |t| {
        point(t)
            .with_quantity(
                "source-pressure-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap()
    };
    let (lower, upper, reference) = (context(77.0), context(300.0), context(293.0));
    let density = expansion(&card, &k, &lower, &upper, &reference);
    let prepared = PreparedExpandingThermalTransport::try_new(
        &card,
        density.clone(),
        query(&k, &lower, &upper),
    )
    .unwrap();
    assert_eq!(prepared.density_model(), &density);
    // Transport exists throughout the sourced interval without pretending that
    // the 293 K reference is a measured density at other temperatures.
    assert!(
        card.claims()
            .query_typed(&k[0], &context(77.0), SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    assert_eq!(
        prepared.at_temperature(293.0).unwrap().density_kg_m3(),
        8000.0
    );
    let residual = 7.4646191727e-7;
    for t in [77.0, 87.0, 150.0, 205.0, 275.0, 293.0, 299.0, 300.0] {
        let source = |index| {
            card.claims()
                .query_typed(&k[index], &context(t), SelectionPolicy::SingleClaimOnly)
                .unwrap()
                .evidence
                .value
                .value
        };
        let rho = 8000.0 / ((1.0 + source(1)) / (1.0 + residual)).powi(3);
        let (cp, conductivity) = (source(2), source(3));
        let sample = prepared.at_temperature(t).unwrap();
        close(sample.density_kg_m3(), rho);
        close(sample.specific_heat_j_kg_k(), cp);
        close(sample.conductivity_w_m_k(), conductivity);
        close(sample.volumetric_heat_capacity_j_m3_k(), rho * cp);
        close(sample.diffusivity_m2_s(), conductivity / (rho * cp));
        close(
            sample.effusivity_j_m2_k_sqrt_s(),
            (conductivity * rho * cp).sqrt(),
        );
    }
    for line in prepared
        .source_lines()
        .iter()
        .chain(prepared.density_model().source_lines())
    {
        assert_eq!(line.knots.len(), 11);
        for (_, answer) in &line.knots {
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    for t in [76.0, 301.0] {
        assert!(prepared.at_temperature(t).is_err());
    }
}
