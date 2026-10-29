//! G0/G1/G3: primary elemental/mineral fits reach shared thermal geometry and flux.
#![deny(unsafe_code)]

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::expansion::{
        FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
    },
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const K: &str = "thermal-conductivity";
const CP: &str = "specific-heat-capacity";
const STRAIN: &str = "linear-expansion-relative-to-293k";
const ALPHA: &str = "linear-thermal-expansion-coefficient";
const A: &str = "linear-expansion-a-axis-relative-to-293k";
const C: &str = "linear-expansion-c-axis-relative-to-293k";
const POLY: &str = "linear-expansion-polycrystalline-relative-to-293k";

struct Fit {
    name: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    logarithmic: bool,
    scale: f64,
    knots: usize,
    orientation: f64,
}
struct Reference {
    slug: &'static str,
    family: &'static str,
    fits: &'static [Fit],
}

// Independently transcribed primary coefficients and bounds; stored ordinates
// were evaluated separately with 60-digit mpmath, including the Si erf blend.
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "molybdenum-nist-thermal-reference",
        family: "molybdenum",
        fits: &[
            Fit {
                name: K,
                range: [4.0, 300.0],
                logarithmic: true,
                scale: 1.0,
                knots: 218,
                orientation: 0.0,
                coefficients: &[
                    10.78259, -72.13065, 228.57351, -384.50447, 381.43825, -228.83783, 81.26658,
                    -15.69097, 1.26814,
                ],
            },
            Fit {
                name: STRAIN,
                range: [20.0, 350.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 173,
                orientation: 0.0,
                coefficients: &[
                    -90.912613,
                    -0.127173,
                    0.00266801,
                    -0.0000050432,
                    0.0000000035183,
                ],
            },
        ],
    },
    Reference {
        slug: "beryllium-nist-thermal-reference",
        family: "beryllium",
        fits: &[
            Fit {
                name: CP,
                range: [14.0, 284.0],
                logarithmic: true,
                scale: 1.0,
                knots: 211,
                orientation: 0.0,
                coefficients: &[
                    -526.84477, 2755.4105, -6209.8985, 7859.2257, -6106.2095, 2982.9958,
                    -894.99967, 150.85256, -10.943615,
                ],
            },
            Fit {
                name: A,
                range: [75.0, 500.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 230,
                orientation: 1.0,
                coefficients: &[
                    -134.173,
                    -0.44051,
                    0.0034063,
                    -0.00000065593,
                    -0.0000000018112,
                ],
            },
            Fit {
                name: C,
                range: [75.0, 500.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 230,
                orientation: 1.0,
                coefficients: &[-86.081, -0.26372, 0.0013623, 0.0000032709, -0.0000000048719],
            },
            Fit {
                name: POLY,
                range: [75.0, 500.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 230,
                orientation: 0.0,
                coefficients: &[
                    -120.298,
                    -0.34687,
                    0.0025691,
                    0.00000089848,
                    -0.0000000029422,
                ],
            },
        ],
    },
    Reference {
        slug: "sapphire-nist-thermal-expansion-reference",
        family: "sapphire",
        fits: &[
            Fit {
                name: ALPHA,
                range: [5.0, 80.0],
                logarithmic: true,
                scale: 1e-6,
                knots: 204,
                orientation: 0.0,
                coefficients: &[
                    10.97236, -97.23540, 240.2436, -294.9933, 195.9244, -66.89247, 9.19921, 0.0,
                    0.0,
                ],
            },
            Fit {
                name: STRAIN,
                range: [20.0, 293.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 155,
                orientation: 0.0,
                coefficients: &[
                    -78.850,
                    -0.022346,
                    0.00010185,
                    0.0000055594,
                    -0.0000000085422,
                ],
            },
        ],
    },
    Reference {
        slug: "silicon-nist-thermal-expansion-reference",
        family: "silicon",
        fits: &[Fit {
            name: ALPHA,
            range: [13.0, 600.0],
            logarithmic: false,
            scale: 1e-8,
            knots: 1177,
            orientation: 1.0,
            coefficients: &[
                1.00500e-5,
                -5.99688e-6,
                1.25574e-6,
                -1.12086e-7,
                3.63225e-9,
                2.67708e-2,
                -1.22829e-4,
                1.62544e-18,
                4.72374e2,
                -3.58796e4,
                -1.24191e7,
                1.25972e9,
            ],
        }],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(name: &str) -> PropertyKey {
    match name {
        K => PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity)),
        CP => PropertyKey::new(name, Dims([2, 0, -2, -1, 0, 0])),
        ALPHA => PropertyKey::new(name, Dims([0, 0, 0, -1, 0, 0])),
        _ => PropertyKey::new(name, Dims::NONE),
    }
}
fn point(reference: &Reference, fit: &Fit, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-processing-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", fit.orientation),
        ("source-product-form-known", 0.0),
    ] {
        p = p
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if fit.name == ALPHA {
        p = p
            .with_quantity(
                "source-coefficient-length-convention-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
    }
    p.with_quantity(
        format!("source-family-nist-{}", reference.family),
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
    )
    .unwrap()
    .with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        t,
    )
    .unwrap()
}
fn compile(reference: &Reference) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-element-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        reference.family
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
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), reference.fits.len());
    pack
}
fn original(reference: &Reference, fit: &Fit, t: f64) -> f64 {
    if reference.family == "silicon" {
        let p = fit.coefficients;
        // This is the original NIST equation image te.png, including the nested
        // w15 blend, half powers, ninth power, both switches and source1e-8 scale.
        let w15 = 0.5 * (1.0 + fs_math::det::erf(t - 15.0));
        let w52 = 0.5 * (1.0 + fs_math::det::erf(0.2 * (t - 52.0)));
        let w200 = 0.5 * (1.0 + fs_math::det::erf(0.1 * (t - 200.0)));
        let low = 4.8e-5 * t.powi(3)
            + (p[0] * t.powi(5)
                + p[1] * t.powf(5.5)
                + p[2] * t.powi(6)
                + p[3] * t.powf(6.5)
                + p[4] * t.powi(7))
                * w15;
        let mid = -47.6
            + p[5] * (t - 76.0).powi(2)
            + p[6] * (t - 76.0).powi(3)
            + p[7] * (t - 76.0).powi(9);
        let high = p[8] + p[9] / t + p[10] / t.powi(2) + p[11] / t.powi(3);
        return 1e-8 * (low * (1.0 - w52) + mid * w52 * (1.0 - w200) + high * w200);
    }
    let x = if fit.logarithmic { t.log10() } else { t };
    let y = fit.coefficients.iter().rev().fold(0.0, |y, a| y * x + a);
    if fit.logarithmic {
        10.0_f64.powf(y) * fit.scale
    } else {
        y * fit.scale
    }
}
fn value(pack: &NormalizedPack, reference: &Reference, fit: &Fit, t: f64) -> f64 {
    let a = pack
        .claims()
        .query_typed(
            &key(fit.name),
            &point(reference, fit, t),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    pack.claims().verify_receipt(&a.receipt).unwrap();
    assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
    a.evidence.value.value
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs() + absolute,
        "{a:.16e} != {b:.16e}"
    );
}
fn card(reference: &Reference, pack: &NormalizedPack) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: reference.slug.into(),
            phase: "solid".into(),
            process:
                "original source with caller-declared geometry interpretation; no density claim"
                    .into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )
    .unwrap()
}

#[test]
fn g1_g3_original_equations_units_receipts_and_signed_coefficients_are_retained() {
    let mut count = 0;
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let line = pack
                .claims()
                .query_line_typed(
                    &key(fit.name),
                    "temperature",
                    &point(reference, fit, fit.range[0]),
                    &point(reference, fit, fit.range[1]),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(line.knots.len(), fit.knots);
            count += line.knots.len();
            for (t, answer) in &line.knots {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
                // Independent f64 log-polynomial evaluation has cancellation
                // error; that arithmetic allowance is not source uncertainty.
                close(
                    answer.evidence.value.value,
                    original(reference, fit, *t),
                    2e-7,
                    5e-18,
                );
            }
            for pair in line.knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let (relative, absolute) = if reference.family == "silicon" {
                        (0.0, 2.3e-11)
                    } else if fit.name == STRAIN || [A, C, POLY].contains(&fit.name) {
                        (0.0, 3.3e-8)
                    } else {
                        (7e-4, 1e-16)
                    };
                    close(
                        value(&pack, reference, fit, t),
                        original(reference, fit, t),
                        relative,
                        absolute,
                    );
                }
            }
        }
        if reference.family == "silicon" {
            let f = &reference.fits[0];
            assert!(value(&pack, reference, f, 13.0) > 0.0);
            assert!(value(&pack, reference, f, 76.0) < 0.0);
            assert!(value(&pack, reference, f, 300.0) > 0.0);
        }
    }
    assert_eq!(count, 2828);
}

#[test]
fn g0_unknown_conditions_coefficient_conventions_and_absent_laws_do_not_disappear() {
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let t = 0.5 * (fit.range[0] + fit.range[1]);
            for axis in [
                "engineering-reference-model",
                "source-processing-known",
                "source-pressure-known",
                "source-specimens-matched",
                "source-product-form-known",
                "source-orientation-known",
            ] {
                let p = point(reference, fit, t);
                assert!(
                    pack.claims()
                        .query_typed(
                            &key(fit.name),
                            &p.clone().with(axis, 1.0 - p.axes()[axis]).unwrap(),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            for t in [fit.range[0].next_down(), fit.range[1].next_up()] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &key(fit.name),
                            &point(reference, fit, t),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            assert!(
                pack.claims()
                    .query_typed(
                        &key(fit.name),
                        &QueryPoint::new(),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            let wrong = point(reference, fit, t)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    t,
                )
                .unwrap();
            assert!(
                pack.claims()
                    .query_typed(&key(fit.name), &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
            if fit.name == ALPHA {
                let known = point(reference, fit, t)
                    .with("source-coefficient-length-convention-known", 1.0)
                    .unwrap();
                assert!(
                    pack.claims()
                        .query_typed(&key(ALPHA), &known, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
                let lower = point(reference, fit, fit.range[0]);
                let upper = point(reference, fit, fit.range[1]);
                // An inverse-K coefficient cannot silently become a finite strain.
                assert!(
                    PreparedFreeExpansion::try_new(
                        &card(reference, &pack),
                        FreeExpansionGeometryQuery {
                            reference: &lower,
                            strains: FreeExpansionStrainSources::Isotropic {
                                strain: &key(ALPHA),
                                selection: SINGLE
                            },
                            temperature_axis: "temperature",
                            lower: &lower,
                            upper: &upper
                        }
                    )
                    .is_err()
                );
            } else if fit.name != K && fit.name != CP {
                assert_ne!(
                    value(&pack, reference, fit, 293.0),
                    0.0,
                    "retain the rounded source datum residual"
                );
            }
        }
        assert!(pack.claims().claims_for("density").is_empty());
        if reference.family != "molybdenum" {
            assert!(pack.claims().claims_for(K).is_empty());
        }
        if reference.family != "beryllium" {
            assert!(pack.claims().claims_for(CP).is_empty());
        }
    }
}

#[test]
fn g1_g3_real_sources_drive_density_independent_geometry_and_shared_conduction() {
    for (index, fit_index) in [(0, 1), (2, 1)] {
        let reference = &REFERENCES[index];
        let fit = &reference.fits[fit_index];
        let pack = compile(reference);
        let material = card(reference, &pack);
        let lower = point(reference, fit, 100.0);
        let upper = point(reference, fit, 250.0);
        // Isotropy is an explicit geometry-fixture assumption, not qualification
        // of NIST's unreported direction. No dummy density is inserted.
        let geometry = PreparedFreeExpansion::try_new(
            &material,
            FreeExpansionGeometryQuery {
                reference: &upper,
                strains: FreeExpansionStrainSources::Isotropic {
                    strain: &key(STRAIN),
                    selection: SINGLE,
                },
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
            },
        )
        .unwrap();
        for t in [100.0, 150.0, 200.0, 250.0] {
            let ratio =
                (1.0 + original(reference, fit, t)) / (1.0 + original(reference, fit, 250.0));
            let sample = geometry.at_temperature(t).unwrap();
            for stretch in sample.principal_stretches() {
                close(stretch, ratio, 5e-14, 0.0);
            }
            close(sample.relative_volume(), ratio.powi(3), 5e-14, 0.0);
        }
        assert!(geometry.at_temperature(99.0).is_err());
        for line in geometry.source_lines() {
            for (_, answer) in &line.knots {
                material.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
        assert!(material.claims().claims_for("density").is_empty());
    }

    let reference = &REFERENCES[1];
    let pack = compile(reference);
    let material = card(reference, &pack);
    let a = &reference.fits[1];
    let c = &reference.fits[2];
    let lower = point(reference, a, 100.0);
    let upper = point(reference, a, 250.0);
    let a_key = key(A);
    let c_key = key(C);
    // The fixture declares mutually orthogonal basal a/a and axial c directions.
    // This uses the two source fits without substituting the separate polycrystal.
    let geometry = PreparedFreeExpansion::try_new(
        &material,
        FreeExpansionGeometryQuery {
            reference: &upper,
            strains: FreeExpansionStrainSources::FixedPrincipal {
                frame: ContentHash([9; 32]),
                strains: [&a_key, &a_key, &c_key],
                selections: [SINGLE; 3],
            },
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
        },
    )
    .unwrap();
    let sample = geometry.at_temperature(100.0).unwrap();
    let ra = (1.0 + original(reference, a, 100.0)) / (1.0 + original(reference, a, 250.0));
    let rc = (1.0 + original(reference, c, 100.0)) / (1.0 + original(reference, c, 250.0));
    for (actual, expected) in sample.principal_stretches().into_iter().zip([ra, ra, rc]) {
        close(actual, expected, 5e-14, 0.0);
    }
    close(sample.relative_volume(), ra * ra * rc, 5e-14, 0.0);
    assert_ne!(ra, rc);

    let reference = &REFERENCES[0];
    let fit = &reference.fits[0];
    let pack = compile(reference);
    let line = pack
        .claims()
        .query_line_typed(
            &key(K),
            "temperature",
            &point(reference, fit, 100.0),
            &point(reference, fit, 250.0),
            SINGLE,
        )
        .unwrap();
    let curve = ConductivityCurve::new(
        line.knots
            .iter()
            .map(|(t, a)| [*t, a.evidence.value.value])
            .collect(),
    )
    .unwrap();
    let model = TwoLayerConduction::new([&curve; 2], [0.01, 0.03], 0.002, 0.0).unwrap();
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    };
    let result = model.evaluate([250.0, 100.0], config, || Ok(())).unwrap();
    // Independent Simpson integral of the original NIST conductivity equation.
    let n = 2048;
    let h = 150.0 / f64::from(n);
    let mut integral = original(reference, fit, 100.0) + original(reference, fit, 250.0);
    for i in 1..n {
        integral +=
            if i % 2 == 0 { 2.0 } else { 4.0 } * original(reference, fit, 100.0 + f64::from(i) * h);
    }
    close(
        result.power_a_to_b_w,
        0.002 / 0.04 * integral * h / 3.0,
        1.5e-4,
        1e-8,
    );
    let reverse = model.evaluate([100.0, 250.0], config, || Ok(())).unwrap();
    close(reverse.power_a_to_b_w, -result.power_a_to_b_w, 5e-13, 1e-9);
    assert!(result.layer_heat_residual_w.abs() < 1e-8);
}
