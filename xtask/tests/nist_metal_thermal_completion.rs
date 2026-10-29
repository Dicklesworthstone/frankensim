//! G0/G1/G3: NIST metal fits through the source compiler and thermal consumers.
//! Numerical fit/interpolation checks do not qualify a manufactured specimen.
#![deny(unsafe_code)]

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
use std::{path::Path, process::Command};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const E: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fit {
    property: &'static str,
    range: [f64; 2],
    log: bool,
    scale: f64,
    coefficients: &'static [f64],
    knots: usize,
}
struct Reference {
    slug: &'static str,
    family: &'static str,
    temper_known: f64,
    fits: &'static [Fit],
}
// Independent transcription of the source equations, not copied ordinates.
const AL_CP: &[f64] = &[
    46.6467, -314.292, 866.662, -1298.3, 1162.27, -637.795, 210.351, -38.3094, 2.96344,
];
const AL_STRAIN: &[f64] = &[-412.77, -0.30389, 0.0087696, -9.9821e-6, 0.0];
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "aluminum-3003-f-nist-thermal-reference",
        family: "aluminum-3003-f",
        temper_known: 1.0,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                range: [4.0, 300.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    0.63736, -1.1437, 7.4624, -12.6905, 11.9165, -6.18721, 1.63939, -0.172667, 0.0,
                ],
                knots: 140,
            },
            Fit {
                property: "specific-heat-capacity",
                range: [4.0, 300.0],
                log: true,
                scale: 1.0,
                coefficients: AL_CP,
                knots: 140,
            },
            Fit {
                property: "linear-expansion-relative-to-293k",
                range: [20.0, 300.0],
                log: false,
                scale: 1e-5,
                coefficients: AL_STRAIN,
                knots: 62,
            },
        ],
    },
    Reference {
        slug: "aluminum-5083-o-nist-thermal-reference",
        family: "aluminum-5083-o",
        temper_known: 1.0,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                range: [4.0, 300.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    -0.90933, 5.751, -11.112, 13.612, -9.3977, 3.6873, -0.77295, 0.067336, 0.0,
                ],
                knots: 140,
            },
            Fit {
                property: "specific-heat-capacity",
                range: [4.0, 300.0],
                log: true,
                scale: 1.0,
                coefficients: AL_CP,
                knots: 140,
            },
            Fit {
                property: "linear-expansion-relative-to-293k",
                range: [20.0, 300.0],
                log: false,
                scale: 1e-5,
                coefficients: AL_STRAIN,
                knots: 62,
            },
            Fit {
                property: "young-modulus",
                range: [2.0, 295.0],
                log: false,
                scale: 1e9,
                coefficients: &[
                    80.83212,
                    0.01061708,
                    -3.016100e-4,
                    7.561340e-7,
                    -6.994800e-10,
                ],
                knots: 65,
            },
        ],
    },
    Reference {
        slug: "invar-fe36ni-nist-thermal-reference",
        family: "invar-fe36ni",
        temper_known: 0.0,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                range: [4.0, 300.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    -2.7064, 8.5191, -15.923, 18.276, -11.9116, 4.40318, -0.86018, 0.068508, 0.0,
                ],
                knots: 140,
            },
            Fit {
                property: "specific-heat-capacity",
                range: [4.0, 20.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    28.08, -228.23, 777.587, -1448.423, 1596.567, -1040.294, 371.2125, -56.004, 0.0,
                ],
                knots: 124,
            },
            Fit {
                property: "linear-expansion-relative-to-293k",
                range: [90.0, 300.0],
                log: false,
                scale: 1e-5,
                coefficients: &[-52.65, 0.1009, 8.395e-4, -1.973e-6, 8.794e-11],
                knots: 47,
            },
            Fit {
                property: "young-modulus",
                range: [3.0, 298.0],
                log: false,
                scale: 1e9,
                coefficients: &[141.565, 0.0254435, -0.00100842, 6.72797e-6, -1.08230e-8],
                knots: 66,
            },
        ],
    },
    Reference {
        slug: "platinum-nist-thermal-reference",
        family: "platinum",
        temper_known: 0.0,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                range: [4.0, 295.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    -7.33450054,
                    80.8550484,
                    -268.441084,
                    481.629105,
                    -503.890454,
                    314.812622,
                    -115.699394,
                    23.0957119,
                    -1.93361717,
                ],
                knots: 138,
            },
            Fit {
                property: "specific-heat-capacity",
                range: [1.0, 295.0],
                log: true,
                scale: 1.0,
                coefficients: &[
                    -1.6135538,
                    0.95823584,
                    1.4317770,
                    -3.5963989,
                    5.1299735,
                    -2.4186452,
                    -0.12560841,
                    0.34342394,
                    -0.06198179,
                ],
                knots: 139,
            },
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(name: &str) -> PropertyKey {
    match name {
        "thermal-conductivity" => {
            PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
        }
        "specific-heat-capacity" => PropertyKey::new(name, CP),
        "young-modulus" => PropertyKey::new(name, E),
        "linear-expansion-relative-to-293k" => PropertyKey::new(name, Dims::NONE),
        _ => panic!("unexpected property"),
    }
}
fn at(reference: &Reference, temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-processing-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
        ("source-product-form-known", 0.0),
        ("source-temper-known", reference.temper_known),
    ] {
        point = point
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
        .with_quantity(
            format!("source-family-nist-{}", reference.family),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
}
fn compile(reference: &Reference) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-metals-{}-{stamp}-{}.fsmatpk",
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
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), reference.fits.len());
    pack
}
fn original(fit: &Fit, temperature: f64) -> f64 {
    let x = if fit.log {
        temperature.log10()
    } else {
        temperature
    };
    let sum = fit
        .coefficients
        .iter()
        .rev()
        .fold(0.0, |value, coefficient| value * x + coefficient);
    if fit.log {
        10.0_f64.powf(sum) * fit.scale
    } else {
        sum * fit.scale
    }
}
fn value(pack: &NormalizedPack, reference: &Reference, property: &str, temperature: f64) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &key(property),
            &at(reference, temperature),
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
fn compiled_nist_curves_preserve_equations_units_and_original_domains() {
    let mut total_knots = 0;
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let line = pack
                .claims()
                .query_line_typed(
                    &key(fit.property),
                    "temperature",
                    &at(reference, fit.range[0]),
                    &at(reference, fit.range[1]),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(line.knots.len(), fit.knots);
            total_knots += line.knots.len();
            let mut observed: f64 = 0.0;
            for (temperature, answer) in &line.knots {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                // f64 polynomial cancellation is bounded here independently of
                // the 50-digit source-table evaluations and NIST fit errors.
                close(
                    answer.evidence.value.value,
                    original(fit, *temperature),
                    2e-8,
                    1e-15,
                );
            }
            for pair in line.knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let temperature = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let actual = value(&pack, reference, fit.property, temperature);
                    let expected = original(fit, temperature);
                    let difference = if fit.property == "linear-expansion-relative-to-293k" {
                        (actual - expected).abs()
                    } else {
                        (actual / expected - 1.0).abs()
                    };
                    observed = observed.max(difference);
                }
            }
            let sampled_limit = if fit.property == "linear-expansion-relative-to-293k" {
                5.1e-7
            } else {
                0.0016
            };
            assert!(
                observed < sampled_limit,
                "{} {}: {observed}",
                reference.slug,
                fit.property
            );
        }
    }
    assert_eq!(total_knots, 1403);
}

#[test]
fn source_context_branch_joins_missing_density_and_warm_invar_cp_refuse() {
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let property = key(fit.property);
            for t in [fit.range[0] - 0.01, fit.range[1] + 0.01] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &property,
                            &at(reference, t),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            let midpoint = (fit.range[0] + fit.range[1]) * 0.5;
            for flag in [
                "engineering-reference-model",
                "source-processing-known",
                "source-pressure-known",
                "source-specimens-matched",
                "source-orientation-known",
                "source-product-form-known",
                "source-temper-known",
            ] {
                let point = at(reference, midpoint);
                let changed = point.clone().with(flag, 1.0 - point.axes()[flag]).unwrap();
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
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new("density", RHO),
                    &at(reference, 293.0),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(
                        "linear-thermal-expansion-coefficient",
                        Dims([0, 0, 0, -1, 0, 0])
                    ),
                    &at(reference, 293.0),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        if reference.family == "invar-fe36ni" {
            for t in [20.01, 27.0, 293.0] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &key("specific-heat-capacity"),
                            &at(reference, t),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            assert!(
                pack.claims()
                    .query_line_typed(
                        &key("linear-expansion-relative-to-293k"),
                        "temperature",
                        &at(reference, 79.0),
                        &at(reference, 90.0),
                        SINGLE
                    )
                    .is_err()
            );
        } else if reference.temper_known == 1.0 {
            assert!(
                pack.claims()
                    .query_line_typed(
                        &key("linear-expansion-relative-to-293k"),
                        "temperature",
                        &at(reference, 17.0),
                        &at(reference, 20.0),
                        SINGLE
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn sourced_finite_strain_drives_existing_expansion_and_transport_owners() {
    for reference in &REFERENCES[..3] {
        let pack = compile(reference);
        let reference_point = at(reference, 293.0);
        let mut claims = pack.claims().clone();
        // An explicitly authored normalization model, NOT a measured density.
        // Keeping rho_ref=1 makes the output a visible density ratio; nothing
        // adds a density observation to the actual catalog pack.
        let density_key = PropertyKey::new("authored-unit-reference-density", RHO);
        claims.insert_claim(PropertyClaim {
            key: density_key.clone(),
            validity: ValidityDomain::unconstrained().with_quantity("temperature", semantic(QuantityKind::AbsoluteTemperature), 293.0, 293.0),
            value: PropertyValue::Scalar { value: 1.0, dims: RHO },
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance { source: "Authored unit-density normalization for this regression; not NIST material data".into(), license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: reference.slug.into(), phase: "solid".into(),
            process: "source reference plus explicit isotropic free-expansion and unit-density assumptions".into(), revision: 0,
        }, claims, Vec::new()).unwrap();
        let strain_key = key("linear-expansion-relative-to-293k");
        let lower = at(reference, 100.0);
        let upper = at(reference, 295.0);
        let density = PreparedFreeExpansionDensity::try_new(
            &card,
            FreeExpansionSourceQuery {
                density: &density_key,
                density_selection: SINGLE,
                reference: &reference_point,
                strains: FreeExpansionStrainSources::Isotropic {
                    strain: &strain_key,
                    selection: SINGLE,
                },
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
            },
        )
        .unwrap();
        let reference_strain = value(&pack, reference, strain_key.name(), 293.0);
        // The rounded published fit does not land at zero at its nominal datum.
        assert!(reference_strain != 0.0);
        for temperature in [100.0, 150.0, 200.0, 250.0, 293.0, 295.0] {
            let strain = value(&pack, reference, strain_key.name(), temperature);
            let expected = ((1.0 + reference_strain) / (1.0 + strain)).powi(3);
            close(
                density.at_temperature(temperature).unwrap().density_kg_m3(),
                expected,
                8e-14,
                0.0,
            );
        }
        let cp = key("specific-heat-capacity");
        let conductivity = key("thermal-conductivity");
        let result = PreparedExpandingThermalTransport::try_new(
            &card,
            density,
            ExpandingThermalSourceQuery {
                specific_heat: &cp,
                conductivity: &conductivity,
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selections: [SINGLE; 2],
            },
        );
        if reference.family == "invar-fe36ni" {
            assert!(result.is_err(), "warm Invar Cp must not be synthesized");
        } else {
            let transport = result.unwrap();
            for temperature in [100.0, 160.0, 200.0, 293.0, 295.0] {
                let sample = transport.at_temperature(temperature).unwrap();
                let capacity = value(&pack, reference, cp.name(), temperature);
                let k = value(&pack, reference, conductivity.name(), temperature);
                close(sample.specific_heat_j_kg_k(), capacity, 8e-14, 0.0);
                close(
                    sample.diffusivity_m2_s(),
                    k / (capacity * sample.density_kg_m3()),
                    8e-14,
                    0.0,
                );
            }
        }
    }
}
