//! G0/G1/G3: primary alloy references reach heating, flux and finite geometry.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::{
        expansion::{
            FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
        },
        storage::{HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage},
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
const E: &str = "young-modulus";
const STRAIN: &str = "linear-expansion-relative-to-293k";

#[derive(Clone, Copy)]
enum Form {
    Log10,
    Polynomial,
    Constant,
}
struct Fit {
    property: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    form: Form,
    scale: f64,
    values: usize,
    cp_branch: u8,
}
struct Alloy {
    slug: &'static str,
    family: &'static str,
    fits: &'static [Fit],
}
const fn log_fit(
    property: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    values: usize,
    cp_branch: u8,
) -> Fit {
    Fit {
        property,
        range,
        coefficients,
        form: Form::Log10,
        scale: 1.0,
        values,
        cp_branch,
    }
}
const fn poly_fit(
    property: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    scale: f64,
    values: usize,
) -> Fit {
    Fit {
        property,
        range,
        coefficients,
        form: Form::Polynomial,
        scale,
        values,
        cp_branch: 0,
    }
}
const fn plateau(range: [f64; 2], coefficients: &'static [f64]) -> Fit {
    Fit {
        property: STRAIN,
        range,
        coefficients,
        form: Form::Constant,
        scale: 1e-5,
        values: 1,
        cp_branch: 0,
    }
}

// Independent transcription of the original NIST tables. Stored ordinates
// were evaluated separately using 70-digit decimal arithmetic. Identical
// published coefficients do not establish independent experimental replication.
const STAINLESS_STRAIN: &[f64] = &[-295.54, -0.39811, 0.0092683, -0.000020261, 0.000000017127];
const NICKEL_CP: &[f64] = &[
    15503.108, -37280.377, 26788.417, 7010.0877, -22731.651, 15386.526, -5175.7968, 896.97274,
    -64.055866,
];
const NICKEL_STRAIN: &[f64] = &[-210.4, -0.05699, 0.005072, -0.00001381, 0.00000001897];
const ALLOYS: &[Alloy] = &[
    Alloy {
        slug: "stainless-304l-nist-thermal-reference",
        family: "stainless-304l",
        fits: &[
            log_fit(
                K,
                [4.0, 300.0],
                &[
                    -1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199,
                ],
                188,
                0,
            ),
            log_fit(
                CP,
                [4.0, 20.0],
                &[
                    -351.51, 3123.695, -12017.28, 26143.99, -35176.33, 29981.75, -15812.78,
                    4719.64, -610.515,
                ],
                166,
                0,
            ),
            poly_fit(
                STRAIN,
                [23.000000000000004, 300.0],
                STAINLESS_STRAIN,
                1e-5,
                182,
            ),
            plateau([4.0, 22.999999999999996], &[-300.04]),
        ],
    },
    Alloy {
        slug: "stainless-310-nist-thermal-reference",
        family: "stainless-310",
        fits: &[
            log_fit(
                K,
                [4.0, 300.0],
                &[
                    -0.81907, -2.1967, 9.1059, -13.078, 10.853, -5.1269, 1.2583, -0.12395, 0.0,
                ],
                188,
                0,
            ),
            log_fit(
                CP,
                [4.0, 47.0],
                &[
                    20.694, -171.007, 600.6256, -1162.748, 1361.931, -986.2934, 430.093, -102.85,
                    10.275,
                ],
                200,
                1,
            ),
            log_fit(
                CP,
                [47.0, 300.0],
                &[
                    -2755.63, 9704.831, -14618.36, 12202.74, -6092.339, 1818.555, -300.458,
                    21.1942, 0.0,
                ],
                176,
                2,
            ),
            poly_fit(
                E,
                [8.0, 290.0],
                &[
                    206.6588,
                    0.006375129,
                    -0.0004345200,
                    0.000001129930,
                    -0.000000001191600,
                ],
                1e9,
                85,
            ),
            poly_fit(
                STRAIN,
                [23.000000000000004, 300.0],
                STAINLESS_STRAIN,
                1e-5,
                182,
            ),
            plateau([4.0, 22.999999999999996], &[-300.04]),
        ],
    },
    Alloy {
        slug: "titanium-15-3-3-3-nist-thermal-reference",
        family: "titanium-15-3-3-3",
        fits: &[log_fit(
            K,
            [1.4, 300.0],
            &[
                -2.398794842,
                8.970743802,
                -29.19286973,
                54.87139779,
                -59.67137228,
                38.89321714,
                -14.94175848,
                3.111616089,
                -0.270452768,
            ],
            205,
            0,
        )],
    },
    Alloy {
        slug: "beryllium-copper-nist-thermal-reference",
        family: "beryllium-copper",
        fits: &[
            log_fit(
                K,
                [2.0, 80.0],
                &[
                    -0.50015, 1.93190, -1.69540, 0.71218, 1.27880, -1.61450, 0.68722, -0.10501, 0.0,
                ],
                177,
                0,
            ),
            poly_fit(
                STRAIN,
                [24.000000000000004, 300.0],
                &[-311.50, -0.44498, 0.010133, -0.000024718, 0.000000026277],
                1e-5,
                182,
            ),
            plateau([4.0, 23.999999999999996], &[-316.68]),
        ],
    },
    Alloy {
        slug: "nickel-steel-fe2p25ni-nist-thermal-reference",
        family: "nickel-steel-fe2p25ni",
        fits: &[
            log_fit(
                K,
                [4.0, 300.0],
                &[
                    -3.5785, 20.804, -55.826, 87.812, -83.5016, 49.08, -17.4348, 3.4277, -0.28622,
                ],
                188,
                0,
            ),
            log_fit(CP, [60.0, 300.0], NICKEL_CP, 174, 0),
            poly_fit(
                E,
                [0.0, 297.0],
                &[
                    217.656,
                    0.0111923,
                    -0.000354650,
                    0.000000831168,
                    -0.000000000703680,
                ],
                1e9,
                94,
            ),
            poly_fit(STRAIN, [6.000000000000001, 265.0], NICKEL_STRAIN, 1e-5, 151),
            plateau([4.0, 5.999999999999999], &[-210.56]),
        ],
    },
    Alloy {
        slug: "nickel-steel-fe5ni-nist-thermal-reference",
        family: "nickel-steel-fe5ni",
        fits: &[
            log_fit(
                K,
                [4.0, 300.0],
                &[
                    -0.57245, -0.33593, 3.9379, -5.4589, 4.2337, -1.8653, 0.43223, -0.041207, 0.0,
                ],
                188,
                0,
            ),
            log_fit(CP, [60.0, 300.0], NICKEL_CP, 174, 0),
            poly_fit(
                E,
                [0.0, 297.0],
                &[
                    209.182,
                    0.0166974,
                    -0.000206410,
                    -0.000000182360,
                    0.000000000927168,
                ],
                1e9,
                94,
            ),
            poly_fit(STRAIN, [6.000000000000001, 265.0], NICKEL_STRAIN, 1e-5, 151),
            plateau([4.0, 5.999999999999999], &[-210.56]),
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(property: &str) -> PropertyKey {
    match property {
        K => PropertyKey::with_quantity(K, semantic(QuantityKind::ThermalConductivity)),
        CP => PropertyKey::new(CP, Dims([2, 0, -2, -1, 0, 0])),
        E => PropertyKey::new(E, Dims([-1, 1, -2, 0, 0, 0])),
        STRAIN => PropertyKey::new(STRAIN, Dims::NONE),
        _ => unreachable!(),
    }
}
fn point(alloy: &Alloy, fit: &Fit, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (axis, value) in [
        ("engineering-reference-model", 1.0),
        ("source-processing-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
        ("source-product-form-known", 0.0),
        ("source-temper-known", 0.0),
    ] {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if fit.property == CP {
        p = p
            .with_quantity(
                "source-caloric-constraint-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
    }
    if fit.cp_branch != 0 {
        p = p
            .with_quantity(
                "source-specific-heat-fit-branch",
                QuantitySpec::dimensional(Dims::NONE),
                f64::from(fit.cp_branch),
            )
            .unwrap();
    }
    p.with_quantity(
        format!("source-family-nist-{}", alloy.family),
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
fn compile(alloy: &Alloy) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-alloy-completion-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        alloy.slug
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", alloy.slug)))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        alloy.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), alloy.fits.len());
    pack
}
fn card(alloy: &Alloy, pack: &NormalizedPack) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId {
        chemistry: alloy.slug.into(), phase: "solid".into(),
        process: "NIST reference; authored sensible-heating, mass and direction; no density or common-coupon claim".into(),
        revision: 0,
    }, pack.claims().clone(), vec![]).unwrap()
}
fn original(fit: &Fit, t: f64) -> f64 {
    let x = match fit.form {
        Form::Log10 => t.log10(),
        Form::Polynomial => t,
        Form::Constant => 0.0,
    };
    let y = fit.coefficients.iter().rev().fold(0.0, |y, a| y * x + a);
    fit.scale
        * if matches!(fit.form, Form::Log10) {
            10.0_f64.powf(y)
        } else {
            y
        }
}
fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn original_integral(fit: &Fit, lo: f64, hi: f64, entropy: bool) -> f64 {
    let value = |t| original(fit, t) / if entropy { t } else { 1.0 };
    let n = 2048;
    let step = (hi - lo) / f64::from(n);
    let mut sum = value(lo) + value(hi);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(i) * step);
    }
    sum * step / 3.0
}

#[test]
fn g1_g3_primary_equations_units_piecewise_accuracy_and_receipts_survive_compilation() {
    let mut total_values = 0;
    let mut total_claims = 0;
    for alloy in ALLOYS {
        let pack = compile(alloy);
        total_claims += alloy.fits.len();
        for fit in alloy.fits {
            let property = key(fit.property);
            let line = pack
                .claims()
                .query_line_typed(
                    &property,
                    "temperature",
                    &point(alloy, fit, fit.range[0]),
                    &point(alloy, fit, fit.range[1]),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(
                line.knots.len(),
                if matches!(fit.form, Form::Constant) {
                    2
                } else {
                    fit.values
                }
            );
            total_values += fit.values;
            for (t, answer) in &line.knots {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
                // High-order Cp polynomials have substantial f64 cancellation.
                // This arithmetic allowance is not source uncertainty.
                close(answer.evidence.value.value, original(fit, *t), 2e-6, 2e-15);
            }
            for pair in line.knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let answer = pack
                        .claims()
                        .query_typed(
                            &property,
                            &point(alloy, fit, t),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .unwrap();
                    pack.claims().verify_receipt(&answer.receipt).unwrap();
                    if fit.property == STRAIN {
                        close(answer.evidence.value.value, original(fit, t), 0.0, 4.1e-8);
                    } else {
                        close(
                            answer.evidence.value.value,
                            original(fit, t),
                            1.55e-4,
                            1e-10,
                        );
                    }
                }
            }
        }
    }
    // Finite interpolation checks are neither continuous error certificates
    // nor measured uncertainties. Five plateaus each store one value.
    assert_eq!(total_claims, 24);
    assert_eq!(total_values, 3150);
}

#[test]
fn g0_data_support_unknown_conditions_and_source_branches_cannot_be_broadened() {
    for alloy in ALLOYS {
        let pack = compile(alloy);
        for fit in alloy.fits {
            let property = key(fit.property);
            let mid = 0.5 * (fit.range[0] + fit.range[1]);
            let p = point(alloy, fit, mid);
            for axis in [
                "engineering-reference-model",
                "source-processing-known",
                "source-pressure-known",
                "source-specimens-matched",
                "source-orientation-known",
                "source-product-form-known",
                "source-temper-known",
                &format!("source-family-nist-{}", alloy.family),
            ] {
                let wrong = p.clone().with(axis, 1.0 - p.axes()[axis]).unwrap();
                assert!(
                    pack.claims()
                        .query_typed(&property, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
            if fit.property == CP {
                let known_constraint = p
                    .clone()
                    .with("source-caloric-constraint-known", 1.0)
                    .unwrap();
                assert!(
                    pack.claims()
                        .query_typed(
                            &property,
                            &known_constraint,
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            for t in [fit.range[0].next_down(), fit.range[1].next_up()] {
                if t < 0.0 {
                    continue;
                }
                assert!(
                    pack.claims()
                        .query_typed(
                            &property,
                            &point(alloy, fit, t),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .is_err(),
                    "{} {} unexpectedly supports {t}",
                    alloy.slug,
                    fit.property
                );
            }
            let wrong_kind = p
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    mid,
                )
                .unwrap();
            assert!(
                pack.claims()
                    .query_typed(&property, &wrong_kind, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
            let no_context = QueryPoint::new()
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::AbsoluteTemperature),
                    mid,
                )
                .unwrap();
            assert!(
                pack.claims()
                    .query_typed(&property, &no_context, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        for missing in [
            "density",
            "volumetric-heat-capacity",
            "poisson-ratio",
            "yield-strength",
            "linear-thermal-expansion-coefficient",
        ] {
            assert!(pack.claims().claims_for(missing).is_empty());
        }
        if let Some(fit) = alloy
            .fits
            .iter()
            .find(|f| f.property == STRAIN && matches!(f.form, Form::Polynomial))
        {
            let join = fit.range[0].next_down();
            assert!(
                pack.claims()
                    .query_typed(
                        &key(STRAIN),
                        &point(alloy, fit, join),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            assert!(
                pack.claims()
                    .query_line_typed(
                        &key(STRAIN),
                        "temperature",
                        &point(alloy, fit, 4.5),
                        &point(alloy, fit, 100.0),
                        SINGLE,
                    )
                    .is_err()
            );
        }
        if ["titanium-15-3-3-3", "beryllium-copper"].contains(&alloy.family) {
            assert!(pack.claims().claims_for(CP).is_empty());
        }
        if alloy.family == "stainless-304l" {
            assert!(pack.claims().claims_for(E).is_empty());
            assert!(
                pack.claims()
                    .query_typed(
                        &key(CP),
                        &point(alloy, &alloy.fits[1], 25.0),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
    let alloy = &ALLOYS[1];
    let pack = compile(alloy);
    let material = card(alloy, &pack);
    let low_fit = &alloy.fits[1];
    let high_fit = &alloy.fits[2];
    let cp = key(CP);
    // Both rounded 310 fits cover 47 K but do not agree exactly. Selection is
    // explicit; neither averaging nor silent stitching is allowed.
    let low_at_join = pack
        .claims()
        .query_typed(
            &cp,
            &point(alloy, low_fit, 47.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let high_at_join = pack
        .claims()
        .query_typed(
            &cp,
            &point(alloy, high_fit, 47.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!((low_at_join.evidence.value.value - high_at_join.evidence.value.value).abs() > 0.01);
    let no_branch = point(alloy, &alloy.fits[0], 47.0)
        .with_quantity(
            "source-caloric-constraint-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    assert!(
        pack.claims()
            .query_typed(&cp, &no_branch, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    let lower = point(alloy, low_fit, 4.0);
    for upper in [point(alloy, low_fit, 300.0), point(alloy, high_fit, 300.0)] {
        assert!(
            PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
                id: 1,
                card: &material,
                capacity: &cp,
                amount: HeatStorageAmount::MassKg(0.1),
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SINGLE,
            }])
            .is_err()
        );
    }
}

#[test]
fn g1_g3_source_alloys_change_mass_heating_conduction_and_density_free_geometry() {
    let config = LayerConductionConfig {
        max_iterations: 256,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-11,
    };
    let mut powers = Vec::new();
    for alloy in ALLOYS {
        let pack = compile(alloy);
        let material = card(alloy, &pack);
        let fit = &alloy.fits[0];
        let line = pack
            .claims()
            .query_line_typed(
                &key(K),
                "temperature",
                &point(alloy, fit, 60.0),
                &point(alloy, fit, 75.0),
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
        // Fixed-direction, homogeneous layers are authored fixture choices.
        // The source's unspecified direction is not measured isotropy.
        let model = TwoLayerConduction::new([&curve, &curve], [0.004, 0.006], 0.02, 0.0).unwrap();
        let forward = model.evaluate([75.0, 60.0], config, || Ok(())).unwrap();
        close(
            forward.power_a_to_b_w,
            2.0 * original_integral(fit, 60.0, 75.0, false),
            1.55e-4,
            1e-8,
        );
        let reverse = model.evaluate([60.0, 75.0], config, || Ok(())).unwrap();
        close(reverse.power_a_to_b_w, -forward.power_a_to_b_w, 1e-10, 1e-8);
        assert!(forward.layer_heat_residual_w.abs() < 1e-7);
        powers.push(forward.power_a_to_b_w);

        for fit in alloy.fits.iter().filter(|f| f.property == CP) {
            let [lo, hi] = if fit.range[1] == 20.0 {
                [8.0, 15.0]
            } else if fit.cp_branch == 1 {
                [10.0, 40.0]
            } else {
                [100.0, 200.0]
            };
            let lower = point(alloy, fit, lo);
            let upper = point(alloy, fit, hi);
            let cp = key(CP);
            // NIST calls this "specific heat" without identifying Cp versus Cv.
            // Using it for sensible storage is authored here, not reclassified
            // as a measured constant-pressure or constant-volume experiment.
            let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
                id: 1,
                card: &material,
                capacity: &cp,
                amount: HeatStorageAmount::MassKg(0.1),
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SINGLE,
            }])
            .unwrap();
            let energy = storage.energy_change_j(lo, hi).unwrap();
            let cp_line = pack
                .claims()
                .query_line_typed(&cp, "temperature", &lower, &upper, SINGLE)
                .unwrap();
            // Independent trapezoid sum, separate from the production
            // piecewise-affine primitive and inversion implementation.
            let exact_piecewise: f64 = cp_line
                .knots
                .windows(2)
                .map(|p| {
                    0.05 * (p[1].0 - p[0].0)
                        * (p[0].1.evidence.value.value + p[1].1.evidence.value.value)
                })
                .sum();
            close(energy, exact_piecewise, 1e-12, 1e-10);
            close(
                energy,
                0.1 * original_integral(fit, lo, hi, false),
                1.55e-4,
                1e-8,
            );
            close(
                storage.sensible_entropy_change_j_k(lo, hi).unwrap(),
                0.1 * original_integral(fit, lo, hi, true),
                1.55e-4,
                1e-9,
            );
            let target = lo + 0.7 * (hi - lo);
            let delivered = storage.energy_change_j(lo, target).unwrap();
            close(
                storage
                    .temperature_after_energy_change(lo, delivered)
                    .unwrap(),
                target,
                1e-11,
                1e-9,
            );
            assert!(
                storage
                    .temperature_after_energy_change(lo, energy * 1.01)
                    .is_err()
            );
        }

        if let Some(fit) = alloy
            .fits
            .iter()
            .find(|f| f.property == STRAIN && matches!(f.form, Form::Polynomial))
        {
            let strain = key(STRAIN);
            let lower = point(alloy, fit, 100.0);
            let upper = point(alloy, fit, 250.0);
            let reference = point(alloy, fit, 200.0);
            // Isotropy is only an explicit geometry-fixture assumption. L293
            // in the source formula does not require a query at 293 K.
            let geometry = PreparedFreeExpansion::try_new(
                &material,
                FreeExpansionGeometryQuery {
                    reference: &reference,
                    strains: FreeExpansionStrainSources::Isotropic {
                        strain: &strain,
                        selection: SINGLE,
                    },
                    temperature_axis: "temperature",
                    lower: &lower,
                    upper: &upper,
                },
            )
            .unwrap();
            for t in [100.0, 150.0, 200.0, 250.0] {
                let ratio = (1.0 + original(fit, t)) / (1.0 + original(fit, 200.0));
                let sample = geometry.at_temperature(t).unwrap();
                for stretch in sample.principal_stretches() {
                    close(stretch, ratio, 1e-13, 0.0);
                }
                close(sample.relative_volume(), ratio.powi(3), 3e-13, 0.0);
            }
            if alloy.family.starts_with("nickel-steel") {
                let unsupported_reference = point(alloy, fit, 293.0);
                assert!(
                    PreparedFreeExpansion::try_new(
                        &material,
                        FreeExpansionGeometryQuery {
                            reference: &unsupported_reference,
                            strains: FreeExpansionStrainSources::Isotropic {
                                strain: &strain,
                                selection: SINGLE
                            },
                            temperature_axis: "temperature",
                            lower: &lower,
                            upper: &upper,
                        }
                    )
                    .is_err()
                );
            }
        }
        assert!(material.claims().claims_for("density").is_empty());
    }
    // Source transport substitution is useful without claiming a controlled
    // study of grade, temper, direction, or common specimens.
    assert!(powers[3] > 4.0 * powers[2]);
    assert!((powers[4] - powers[5]).abs() > 1.0);
}
