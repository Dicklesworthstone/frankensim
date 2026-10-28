//! Source compiler and bounded thermal-query checks for NIST polymer references.
//! The interpolation comparisons are numerical observations, not confidence bounds
//! on measurements, physical qualification, or a density/phase-transition model.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, InterpolationPolicy, NormalizedPack, PropertyKey, PropertyValue, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fit {
    property: &'static str,
    dims: Dims,
    range: (f64, f64),
    coefficients: [f64; 9],
    knots: usize,
}

struct Reference {
    slug: &'static str,
    family: &'static str,
    nominal_density: Option<f64>,
    helium: bool,
    fits: &'static [Fit],
}

// Independently transcribed a..i coefficients from the cited NIST pages.
// These evaluate the original polynomial, never the source-table ordinates.
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "polyamide-nylon-nist-thermal-reference",
        family: "nylon",
        nominal_density: None,
        helium: false,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                dims: K,
                range: (4.0, 300.0),
                knots: 91,
                coefficients: [
                    -2.6135, 2.3239, -4.7586, 7.1602, -4.9155, 1.6324, -0.2507, 0.0131, 0.0,
                ],
            },
            Fit {
                property: "specific-heat-capacity",
                dims: CP,
                range: (4.0, 300.0),
                knots: 91,
                coefficients: [
                    -5.2929, 25.301, -54.874, 71.061, -52.236, 21.648, -4.7317, 0.42518, 0.0,
                ],
            },
        ],
    },
    Reference {
        slug: "polyimide-kapton-nist-thermal-reference",
        family: "kapton",
        nominal_density: None,
        helium: false,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                dims: K,
                range: (4.0, 300.0),
                knots: 91,
                coefficients: [
                    5.73101, -39.5199, 79.9313, -83.8572, 50.9157, -17.9835, 3.42413, -0.27133, 0.0,
                ],
            },
            Fit {
                property: "specific-heat-capacity",
                dims: CP,
                range: (4.0, 300.0),
                knots: 91,
                coefficients: [
                    -1.3684, 0.65892, 2.8719, 0.42651, -3.0088, 1.9558, -0.51998, 0.051574, 0.0,
                ],
            },
        ],
    },
    Reference {
        slug: "pet-mylar-nist-thermal-reference",
        family: "pet-mylar",
        nominal_density: None,
        helium: false,
        fits: &[Fit {
            property: "thermal-conductivity",
            dims: K,
            range: (2.0, 80.0),
            knots: 87,
            coefficients: [
                -1.37737, -3.40668, 20.5842, -53.1244, 73.2476, -57.6546, 26.1192, -6.34790,
                0.640331,
            ],
        }],
    },
    Reference {
        slug: "polyurethane-helium-49p02-nist-thermal-reference",
        family: "polyurethane-helium",
        nominal_density: Some(49.02),
        helium: true,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                dims: K,
                range: (30.0, 300.0),
                knots: 87,
                coefficients: [
                    -33.898, 117.81, -178.376, 142.038, -63.034, 14.958, -1.5468, 0.020625, 0.0,
                ],
            },
            Fit {
                property: "specific-heat-capacity",
                dims: CP,
                range: (23.0, 300.0),
                knots: 88,
                coefficients: [
                    89.69, -269.32, 333.276, -214.635, 76.2052, -14.1137, 1.061, 0.0, 0.0,
                ],
            },
        ],
    },
    Reference {
        slug: "polystyrene-99p96-nist-thermal-reference",
        family: "polystyrene",
        nominal_density: Some(99.96),
        helium: false,
        fits: &[
            Fit {
                property: "thermal-conductivity",
                dims: K,
                range: (4.0, 300.0),
                knots: 91,
                coefficients: [
                    7.39582, -59.6737, 160.58, -240.33, 218.817, -124.155, 42.9088, -8.26683,
                    0.68082,
                ],
            },
            Fit {
                property: "specific-heat-capacity",
                dims: CP,
                range: (100.0, 300.0),
                knots: 84,
                coefficients: [
                    -5911.474, 14991.16, -15513.753, 8232.5666, -2237.3586, 225.423, 20.33505,
                    -4.6169, 0.0,
                ],
            },
        ],
    },
];

fn compile(reference: &Reference) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "frankensim-polymer-temperature-{}-{nanos}-{}.fsmatpk",
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
    assert_eq!(pack.pack_id(), reference.slug);
    assert_eq!(
        pack.claims().claims_ordered().count(),
        reference.fits.len() + usize::from(reference.nominal_density.is_some())
    );
    pack
}

fn context(reference: &Reference) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-grade-known", 0.0),
        ("source-processing-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
    ] {
        point = point
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point = point
        .with_quantity(
            format!("source-family-nist-{}", reference.family),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    if let Some(density) = reference.nominal_density {
        point = point
            .with_quantity(
                "source-reference-density",
                QuantitySpec::dimensional(RHO),
                density,
            )
            .unwrap();
        let (name, value) = if reference.helium {
            ("source-cell-gas-helium", 1.0)
        } else {
            ("source-cell-gas-known", 0.0)
        };
        point = point
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn at(reference: &Reference, kelvin: f64) -> QueryPoint {
    context(reference)
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if name == "thermal-conductivity" {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn polynomial(fit: &Fit, temperature: f64) -> f64 {
    let log_t = temperature.log10();
    10_f64.powf(
        fit.coefficients
            .iter()
            .rev()
            .fold(0.0, |sum, coefficient| sum * log_t + coefficient),
    )
}

fn query(pack: &NormalizedPack, fit: &Fit, point: &QueryPoint) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &key(fit.property, fit.dims),
            point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(answer.evidence.value.dims, fit.dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn refuses(pack: &NormalizedPack, property: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(
                &key(property, dims),
                point,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err(),
        "{} {property} unexpectedly accepts {point:?}",
        pack.pack_id()
    );
}

#[test]
fn compiled_curves_follow_independent_nist_fits_and_bound_sampled_linearization_error() {
    let mut total_knots = 0;
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let claims = pack.claims().claims_for(fit.property);
            assert_eq!(claims.len(), 1);
            let claim = claims[0].1;
            assert_eq!(claim.interpolation, InterpolationPolicy::LinearInside);
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            let PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                knots,
                dims,
            } = &claim.value
            else {
                panic!("expected temperature curve")
            };
            assert_eq!(abscissa, "temperature");
            assert_eq!(*abscissa_dims, TEMP);
            assert_eq!(*dims, fit.dims);
            assert_eq!(knots.len(), fit.knots);
            assert_eq!(
                (knots.first().unwrap().0, knots.last().unwrap().0),
                fit.range
            );
            for &(temperature, expected) in knots {
                assert!((expected / polynomial(fit, temperature) - 1.0).abs() < 1e-8);
                assert_eq!(query(&pack, fit, &at(reference, temperature)), expected);
            }
            let mut observed_max = 0.0_f64;
            for pair in knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let temperature = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let actual = query(&pack, fit, &at(reference, temperature));
                    observed_max =
                        observed_max.max((actual / polynomial(fit, temperature) - 1.0).abs());
                }
            }
            // A sampled 0.25% regression budget for our interpolation, separate
            // from NIST's reported 1–4% fit error and any measurement uncertainty.
            assert!(
                observed_max < 0.0025,
                "{} {}: {observed_max}",
                reference.slug,
                fit.property
            );
            println!(
                "{} {} sampled interpolation difference {:.8}%",
                reference.slug,
                fit.property,
                observed_max * 100.0
            );
            total_knots += knots.len();
        }
    }
    assert_eq!(total_knots, 801);
}

#[test]
fn source_limits_reference_opt_in_gas_and_nominal_density_do_not_leak() {
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let midpoint = (fit.range.0 + fit.range.1) / 2.0;
            for temperature in [fit.range.0 - 0.01, fit.range.1 + 0.01] {
                refuses(&pack, fit.property, fit.dims, &at(reference, temperature));
            }
            refuses(&pack, fit.property, fit.dims, &QueryPoint::new());
            for (flag, value) in [
                ("engineering-reference-model", 0.0),
                ("source-specimens-matched", 1.0),
                ("source-grade-known", 1.0),
                ("source-processing-known", 1.0),
                ("source-moisture-condition-known", 1.0),
                ("source-orientation-known", 1.0),
            ] {
                refuses(
                    &pack,
                    fit.property,
                    fit.dims,
                    &at(reference, midpoint).with(flag, value).unwrap(),
                );
            }
            let wrong_temperature = context(reference)
                .with_quantity("temperature", QuantitySpec::dimensional(TEMP), midpoint)
                .unwrap();
            refuses(&pack, fit.property, fit.dims, &wrong_temperature);
            if let Some(density) = reference.nominal_density {
                refuses(
                    &pack,
                    fit.property,
                    fit.dims,
                    &at(reference, midpoint)
                        .with("source-reference-density", density + 1.0)
                        .unwrap(),
                );
                let (gas, wrong) = if reference.helium {
                    ("source-cell-gas-helium", 0.0)
                } else {
                    ("source-cell-gas-known", 1.0)
                };
                refuses(
                    &pack,
                    fit.property,
                    fit.dims,
                    &at(reference, midpoint).with(gas, wrong).unwrap(),
                );
            }
            let line = pack
                .claims()
                .query_line_typed(
                    &key(fit.property, fit.dims),
                    "temperature",
                    &at(reference, fit.range.0),
                    &at(reference, fit.range.1),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(line.knots.len(), fit.knots);
            assert!(
                pack.claims()
                    .query_line_typed(
                        &key(fit.property, fit.dims),
                        "temperature",
                        &at(reference, fit.range.0 - 0.01),
                        &at(reference, fit.range.1),
                        SINGLE
                    )
                    .is_err()
            );
        }
        // Source column density labels must never become a thermal solver's rho.
        refuses(&pack, "density", RHO, &at(reference, 200.0));
        if let Some(density) = reference.nominal_density {
            let point = context(reference)
                .with_quantity(
                    "source-density-temperature-known",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                )
                .unwrap();
            let answer = pack
                .claims()
                .query_typed(
                    &key("nominal-source-density", RHO),
                    &point,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert_eq!(answer.evidence.value.value, density);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            refuses(
                &pack,
                "nominal-source-density",
                RHO,
                &point.with("source-density-temperature-known", 1.0).unwrap(),
            );
        }
        if reference.family == "pet-mylar" {
            refuses(&pack, "specific-heat-capacity", CP, &at(reference, 77.0));
            refuses(&pack, "thermal-conductivity", K, &at(reference, 293.15));
        }
    }
}

#[test]
fn heat_capacity_lines_support_temperature_integrals_without_freezing_the_property() {
    for reference in REFERENCES {
        let Some(fit) = reference
            .fits
            .iter()
            .find(|fit| fit.property == "specific-heat-capacity")
        else {
            continue;
        };
        let pack = compile(reference);
        let line = pack
            .claims()
            .query_line_typed(
                &key(fit.property, CP),
                "temperature",
                &at(reference, 100.0),
                &at(reference, 300.0),
                SINGLE,
            )
            .unwrap();
        // Exact integral of the admitted piecewise-affine Cp, with each original
        // source knot retained by the canonical query path. Units are J/kg.
        let integral: f64 = line
            .knots
            .windows(2)
            .map(|pair| {
                (pair[1].0 - pair[0].0)
                    * (pair[0].1.evidence.value.value + pair[1].1.evidence.value.value)
                    * 0.5
            })
            .sum();
        // Independent Simpson quadrature of the published NIST polynomial.
        let intervals = 4096;
        let step = 200.0 / f64::from(intervals);
        let mut sum = polynomial(fit, 100.0) + polynomial(fit, 300.0);
        for i in 1..intervals {
            sum +=
                if i % 2 == 0 { 2.0 } else { 4.0 } * polynomial(fit, 100.0 + step * f64::from(i));
        }
        let source_integral = sum * step / 3.0;
        let observed_difference = (integral / source_integral - 1.0).abs();
        assert!(
            observed_difference < 0.0005,
            "{}: {observed_difference}",
            reference.slug
        );
        let frozen_at_100k = query(&pack, fit, &at(reference, 100.0)) * 200.0;
        assert!((frozen_at_100k / integral - 1.0).abs() > 0.15);
        println!(
            "{} integral Cp(100..300 K)={integral:.8} J/kg; observed polynomial-integral difference {:.8}%",
            reference.slug,
            observed_difference * 100.0
        );
    }
}
