//! Primary-source composite curves through canonical queries and finite expansion.
//! Test-only missing-axis and unit-density laws are explicit caller assumptions.
#![deny(unsafe_code)]

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::thermal::transport::{
    expanding::{ExpandingThermalSourceQuery, PreparedExpandingThermalTransport},
    expansion::{
        FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
    },
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command, sync::OnceLock};

const K: &str = "thermal-conductivity";
const NORMAL: &str = "linear-expansion-normal-relative-to-293k";
const WARP: &str = "linear-expansion-warp-relative-to-293k";
const TEST: &str = "linear-expansion-test-direction-relative-to-293k";
const GRAIN: &str = "linear-expansion-grain-relative-to-293k";
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const FLAGS: &[(&str, f64)] = &[
    ("engineering-reference-model", 1.0),
    ("source-pressure-known", 0.0),
    ("source-specimens-matched", 0.0),
    ("source-cure-known", 0.0),
    ("source-fiber-volume-fraction-known", 0.0),
    ("source-resin-formulation-known", 0.0),
];

struct Fit {
    property: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    knots: usize,
    context: &'static [(&'static str, f64)],
}
struct Reference {
    slug: &'static str,
    claims: usize,
    context: &'static [(&'static str, f64)],
    fits: &'static [Fit],
}
const fn strain(
    property: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    knots: usize,
) -> Fit {
    Fit {
        property,
        range,
        coefficients,
        knots,
        context: &[],
    }
}
const fn conductivity(
    range: [f64; 2],
    coefficients: &'static [f64],
    knots: usize,
    context: &'static [(&'static str, f64)],
) -> Fit {
    Fit {
        property: K,
        range,
        coefficients,
        knots,
        context,
    }
}
// Independent transcription of the published equations, not sampled ordinates.
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "g10cr-nist-cryogenic",
        claims: 5,
        context: &[],
        fits: &[
            strain(
                NORMAL,
                [10.0, 300.0],
                &[-719.8, 0.4455, 0.007505, -2.219e-6, 0.0],
                150,
            ),
            strain(
                WARP,
                [12.0, 300.0],
                &[-246.9, 0.2064, 0.003072, -3.226e-6, 0.0],
                149,
            ),
        ],
    },
    Reference {
        slug: "glass-fabric-polyester-nist-thermal-reference",
        claims: 5,
        context: &[],
        fits: &[
            conductivity(
                [38.0, 300.0],
                &[
                    689.532, -2543.63, 3967.067, -3400.366, 1731.725, -524.309, 87.4249, -6.19597,
                    0.0,
                ],
                136,
                &[
                    ("conductivity-direction", 2.0),
                    ("source-conductivity-medium", 1.0),
                ],
            ),
            conductivity(
                [80.0, 300.0],
                &[
                    -2141.58, 4639.74, -3249.405, 51.72425, 1101.977, -613.8563, 141.1432,
                    -12.3133, 0.0,
                ],
                114,
                &[
                    ("conductivity-direction", 2.0),
                    ("source-conductivity-medium", 2.0),
                ],
            ),
            conductivity(
                [84.0, 300.0],
                &[
                    2909.905, -8616.64, 10542.69, -6832.068, 2475.328, -475.7, 37.9003, 0.0, 0.0,
                ],
                112,
                &[
                    ("conductivity-direction", 1.0),
                    ("source-conductivity-medium", 2.0),
                ],
            ),
            strain(
                NORMAL,
                [80.0, 285.0],
                &[-717.9, -3.157, 0.05251, -0.0001947, 0.0000002752],
                105,
            ),
            strain(
                WARP,
                [80.0, 280.0],
                &[-308.97, 1.0245, -0.0029503, 0.000018323, -0.000000027013],
                102,
            ),
        ],
    },
    Reference {
        slug: "glass-mat-epoxy-chopped-nist-expansion-reference",
        claims: 1,
        context: &[("source-test-direction-known", 0.0)],
        fits: &[strain(
            TEST,
            [77.0, 293.0],
            &[-658.98, 4.7697, -0.029638, 0.00011501, -0.00000014763],
            112,
        )],
    },
    Reference {
        slug: "glass-mat-epoxy-continuous-nist-expansion-reference",
        claims: 1,
        context: &[("source-test-direction-known", 0.0)],
        fits: &[strain(
            TEST,
            [77.0, 293.0],
            &[-439.2, 1.525, -0.002384, 0.000008665, -0.000000002857],
            112,
        )],
    },
    Reference {
        slug: "beechwood-phenolic-crosslaminate-nist-thermal-reference",
        claims: 3,
        context: &[("source-cross-laminate-0-90", 1.0)],
        fits: &[
            conductivity(
                [92.0, 300.0],
                &[
                    -1375.11, 3740.69, -4238.465, 2559.333, -868.6067, 157.1018, -11.82957, 0.0,
                    0.0,
                ],
                108,
                &[("conductivity-direction", 2.0)],
            ),
            conductivity(
                [92.0, 300.0],
                &[
                    1035.33, -2191.85, 1470.505, 39.845, -541.9035, 289.844, -65.2253, 5.59956, 0.0,
                ],
                108,
                &[("conductivity-direction", 1.0)],
            ),
            strain(
                GRAIN,
                [77.0, 293.0],
                &[-944.081, 2.909782, 0.002971, -0.00001, 0.0000000123],
                112,
            ),
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(name: &str) -> PropertyKey {
    if name == K {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, Dims::NONE)
    }
}
fn base_at(reference: &Reference, temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap();
    for (name, value) in FLAGS.iter().chain(reference.context.iter()) {
        point = point
            .with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap();
    }
    if reference.slug != "g10cr-nist-cryogenic" {
        point = point
            .with_quantity(
                format!("source-family-{}", reference.slug),
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    point
}
fn at(reference: &Reference, fit: &Fit, temperature: f64) -> QueryPoint {
    let mut point = base_at(reference, temperature);
    for (name, value) in fit.context {
        point = point
            .with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap();
    }
    point
}
fn pack(index: usize) -> &'static NormalizedPack {
    static PACKS: OnceLock<Vec<NormalizedPack>> = OnceLock::new();
    &PACKS.get_or_init(|| {
        REFERENCES
            .iter()
            .map(|reference| {
                let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
                let stamp = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos();
                let output = std::env::temp_dir().join(format!(
                    "nist-composites-{}-{stamp}-{}.fsmatpk",
                    std::process::id(),
                    reference.slug
                ));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(root.join(format!(
                        "data/matdb/seed-v1/{}/manifest.tsv",
                        reference.slug
                    )))
                    .arg("--out")
                    .arg(&output)
                    .current_dir(root)
                    .env("CARGO_WORKSPACE_DIR", root)
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{}: {}",
                    reference.slug,
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(output).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.claims().claims_ordered().count(), reference.claims);
                pack
            })
            .collect()
    })[index]
}
fn original(fit: &Fit, temperature: f64) -> f64 {
    let x = if fit.property == K {
        temperature.log10()
    } else {
        temperature
    };
    let y = fit
        .coefficients
        .iter()
        .rev()
        .fold(0.0, |sum, coefficient| sum * x + coefficient);
    if fit.property == K {
        10.0_f64.powf(y)
    } else {
        y * 1e-5
    }
}
fn value(pack: &NormalizedPack, reference: &Reference, fit: &Fit, temperature: f64) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &key(fit.property),
            &at(reference, fit, temperature),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    answer.evidence.value.value
}
fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn composite_curves_preserve_equations_units_domains_and_fit_residuals() {
    let mut knots = 0;
    for (index, reference) in REFERENCES.iter().enumerate() {
        let pack = pack(index);
        for fit in reference.fits {
            let line = pack
                .claims()
                .query_line_typed(
                    &key(fit.property),
                    "temperature",
                    &at(reference, fit, fit.range[0]),
                    &at(reference, fit, fit.range[1]),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(line.knots.len(), fit.knots);
            knots += line.knots.len();
            for (temperature, answer) in &line.knots {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                // Loose only for cancellation in independent f64 log-polynomials.
                close(
                    answer.evidence.value.value,
                    original(fit, *temperature),
                    5e-8,
                    1e-15,
                );
            }
            for pair in line.knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let temperature = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let actual = value(pack, reference, fit, temperature);
                    let expected = original(fit, temperature);
                    if fit.property == K {
                        close(actual, expected, 0.0004, 0.0);
                    } else {
                        close(actual, expected, 0.0, 2e-7);
                    }
                }
            }
            if fit.property != K && fit.range[1] >= 293.0 {
                let residual = value(pack, reference, fit, 293.0);
                assert_ne!(
                    residual, 0.0,
                    "rounded source datum must not be forced to zero"
                );
                close(residual, original(fit, 293.0), 0.0, 1e-15);
            }
        }
    }
    assert_eq!(knots, 1420);
}

#[test]
fn context_missing_directions_and_nominal_datum_do_not_extend_support() {
    for (index, reference) in REFERENCES.iter().enumerate() {
        let pack = pack(index);
        for fit in reference.fits {
            let property = key(fit.property);
            for temperature in [fit.range[0] - 0.01, fit.range[1] + 0.01] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &property,
                            &at(reference, fit, temperature),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            let point = at(reference, fit, 150.0);
            for (name, current) in FLAGS.iter().chain(reference.context.iter()) {
                let changed = point.clone().with(*name, 1.0 - current).unwrap();
                assert!(
                    pack.claims()
                        .query_typed(&property, &changed, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{} {name}",
                    reference.slug
                );
            }
            for (name, _) in fit.context {
                let changed = point.clone().with(*name, 99.0).unwrap();
                assert!(
                    pack.claims()
                        .query_typed(&property, &changed, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
            assert!(
                pack.claims()
                    .query_typed(
                        &property,
                        &QueryPoint::new(),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        let point = base_at(reference, 250.0);
        for missing in [
            PropertyKey::new("density", RHO),
            key("linear-expansion-fill-relative-to-293k"),
            PropertyKey::new(
                "linear-thermal-expansion-coefficient",
                Dims([0, 0, 0, -1, 0, 0]),
            ),
        ] {
            assert!(
                pack.claims()
                    .query_typed(&missing, &point, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        if index != 0 {
            assert!(
                pack.claims()
                    .query_typed(
                        &PropertyKey::new("specific-heat-capacity", CP),
                        &point,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
    let polyester = &REFERENCES[1];
    for fit in polyester.fits.iter().filter(|fit| fit.property != K) {
        // 293 K defines the original denominator; it is outside these equations.
        assert!(
            pack(1)
                .claims()
                .query_typed(
                    &key(fit.property),
                    &at(polyester, fit, 293.0),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
    let unsupported_helium_normal = base_at(polyester, 150.0)
        .with_quantity(
            "conductivity-direction",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "source-conductivity-medium",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    assert!(
        pack(1)
            .claims()
            .query_typed(
                &key(K),
                &unsupported_helium_normal,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

fn authored_scalar(name: &str, dims: Dims, value: f64, range: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        validity: ValidityDomain::unconstrained().with_quantity("temperature", semantic(QuantityKind::AbsoluteTemperature), range[0], range[1]),
        value: PropertyValue::Scalar { value, dims },
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance { source: "Explicit regression-only normalization/absent-direction model; not measured NIST material data".into(), license: "CC0-1.0".into(), artifact: None },
    }
}

#[test]
fn explicit_missing_axis_models_drive_bounded_principal_expansion_and_g10_transport() {
    for (index, reference) in REFERENCES.iter().enumerate() {
        let pack = pack(index);
        let fits: Vec<_> = reference
            .fits
            .iter()
            .filter(|fit| fit.property != K)
            .collect();
        let first = key(fits[0].property);
        let missing = key("authored-missing-axis-zero-strain");
        let second = fits
            .get(1)
            .map_or_else(|| missing.clone(), |fit| key(fit.property));
        let rho = PropertyKey::new("authored-unit-reference-density", RHO);
        let point = |temperature| {
            // Fixed normal conduction selection is caller-supplied for G-10 only.
            let p = base_at(reference, temperature);
            if index == 0 {
                p.with_quantity(
                    "conductivity-direction",
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                )
                .unwrap()
            } else {
                p
            }
        };
        let (lower, upper, datum) = (point(100.0), point(280.0), point(250.0));
        let query = FreeExpansionSourceQuery {
            density: &rho,
            density_selection: SINGLE,
            reference: &datum,
            strains: FreeExpansionStrainSources::FixedPrincipal {
                frame: ContentHash([11; 32]),
                strains: [&first, &second, &missing],
                selections: [SINGLE; 3],
            },
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
        };
        let state = MaterialStateId { chemistry: reference.slug.into(), phase: "solid".into(), process: "source test directions embedded in a caller-declared fixed orthogonal frame; absent directions explicitly held constant".into(), revision: 0 };
        let mut claims = pack.claims().clone();
        let source_only =
            MaterialCard::assemble(state.clone(), claims.clone(), Vec::new()).unwrap();
        assert!(PreparedFreeExpansionDensity::try_new(&source_only, query).is_err());
        // Unit reference density makes the numerical result a density ratio.
        claims
            .insert_claim(authored_scalar(rho.name(), RHO, 1.0, [250.0, 250.0]))
            .unwrap();
        let incomplete = MaterialCard::assemble(state.clone(), claims.clone(), Vec::new()).unwrap();
        assert!(
            PreparedFreeExpansionDensity::try_new(&incomplete, query).is_err(),
            "missing principal directions must never be inferred from warp or a scalar test"
        );
        // This zero law is a hypothetical closure, not a G-10 fill measurement,
        // transverse isotropy, or a catalog claim about an unreported test axis.
        claims
            .insert_claim(authored_scalar(
                missing.name(),
                Dims::NONE,
                0.0,
                [100.0, 280.0],
            ))
            .unwrap();
        let card = MaterialCard::assemble(state, claims, Vec::new()).unwrap();
        let density = PreparedFreeExpansionDensity::try_new(&card, query).unwrap();
        assert!(density.at_temperature(99.99).is_err());
        assert!(density.at_temperature(280.01).is_err());
        for temperature in [100.0, 150.0, 200.0, 250.0, 280.0] {
            let mut expected = [1.0; 3];
            for (axis, fit) in fits.iter().enumerate() {
                expected[axis] = (1.0 + original(fit, temperature)) / (1.0 + original(fit, 250.0));
            }
            let sample = density.at_temperature(temperature).unwrap();
            for (actual, expected) in sample.principal_stretches().into_iter().zip(expected) {
                close(actual, expected, 8e-14, 0.0);
            }
            let volume: f64 = expected.into_iter().product();
            close(sample.relative_volume(), volume, 8e-14, 0.0);
            close(sample.density_kg_m3(), 1.0 / volume, 8e-14, 0.0);
        }
        let cp = PropertyKey::new("specific-heat-capacity", CP);
        let k = key(K);
        let transport = PreparedExpandingThermalTransport::try_new(
            &card,
            density,
            ExpandingThermalSourceQuery {
                specific_heat: &cp,
                conductivity: &k,
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selections: [SINGLE; 2],
            },
        );
        if index == 0 {
            let transport = transport.unwrap();
            for temperature in [100.0, 150.0, 250.0, 280.0] {
                let sample = transport.at_temperature(temperature).unwrap();
                let k = card
                    .claims()
                    .query_typed(&k, &point(temperature), SelectionPolicy::SingleClaimOnly)
                    .unwrap()
                    .evidence
                    .value
                    .value;
                close(
                    sample.diffusivity_m2_s(),
                    k / (sample.density_kg_m3() * sample.specific_heat_j_kg_k()),
                    8e-14,
                    0.0,
                );
            }
        } else {
            assert!(
                transport.is_err(),
                "unreported heat capacity must not be synthesized"
            );
        }
    }
}
