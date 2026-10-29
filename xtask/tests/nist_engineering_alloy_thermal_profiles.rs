//! G0/G1/G3: engineering-alloy source fits reach shared conduction and expansion.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::{
        expanding::{ExpandingThermalSourceQuery, PreparedExpandingThermalTransport},
        expansion::{
            FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
        },
    },
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const STRAIN: &str = "linear-expansion-relative-to-293k";
const K: &str = "thermal-conductivity";
const CP: &str = "specific-heat-capacity";
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);

struct Fit {
    name: &'static str,
    range: [f64; 2],
    coefficients: &'static [f64],
    logarithmic: bool,
    scale: f64,
    knots: usize,
}
struct Reference {
    slug: &'static str,
    family: &'static str,
    temper: f64,
    // Source T_low and its separate constant finite-strain branch.
    plateau: Option<[f64; 2]>,
    fits: &'static [Fit],
}

// Independently transcribed published equations, not stored knot ordinates.
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "aluminum-6063-t5-nist-thermal-reference",
        family: "aluminum-6063-t5",
        temper: 1.0,
        plateau: None,
        fits: &[Fit {
            name: K,
            range: [4.0, 295.0],
            logarithmic: true,
            scale: 1.0,
            knots: 184,
            coefficients: &[
                22.401433, -141.13433, 394.95461, -601.15377, 547.83202, -305.99691, 102.38656,
                -18.810237, 1.4576882,
            ],
        }],
    },
    Reference {
        slug: "titanium-ti6al4v-nist-thermal-reference",
        family: "titanium-ti6al4v",
        temper: 0.0,
        plateau: Some([24.0, -0.0017361]),
        fits: &[
            Fit {
                name: K,
                range: [23.0, 300.0],
                logarithmic: true,
                scale: 1.0,
                knots: 180,
                coefficients: &[
                    -5107.8774, 19240.422, -30789.064, 27134.756, -14226.379, 4438.2154,
                    -763.07767, 55.796592, 0.0,
                ],
            },
            Fit {
                name: STRAIN,
                range: [24.000000000000004, 300.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 147,
                coefficients: &[-171.1, -0.2140, 0.004807, -0.000007111, 0.0],
            },
        ],
    },
    Reference {
        slug: "inconel-718-nist-thermal-reference",
        family: "inconel-718",
        temper: 0.0,
        plateau: Some([20.0, -0.0023887]),
        fits: &[
            Fit {
                name: K,
                range: [6.0, 275.0],
                logarithmic: true,
                scale: 1.0,
                knots: 179,
                coefficients: &[
                    -8.28921, 39.4470, -83.4353, 98.1690, -67.2088, 26.7082, -5.7205, 0.51115, 0.0,
                ],
            },
            Fit {
                name: STRAIN,
                range: [20.000000000000004, 300.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 150,
                coefficients: &[-236.8, -0.2120, 0.005497, -0.000006882, 0.0],
            },
        ],
    },
    Reference {
        slug: "nickel-steel-fe9ni-nist-thermal-reference",
        family: "nickel-steel-fe9ni",
        temper: 0.0,
        plateau: Some([6.0, -0.0021056]),
        fits: &[
            Fit {
                name: K,
                range: [4.0, 300.0],
                logarithmic: true,
                scale: 1.0,
                knots: 186,
                coefficients: &[
                    -0.0712785, -3.48735, 10.6547, -12.9153, 8.89066, -3.51482, 0.743643,
                    -0.0657884, 0.0,
                ],
            },
            Fit {
                name: CP,
                range: [60.0, 300.0],
                logarithmic: true,
                scale: 1.0,
                knots: 174,
                coefficients: &[
                    15503.108, -37280.377, 26788.417, 7010.0877, -22731.651, 15386.526, -5175.7968,
                    896.97274, -64.055866,
                ],
            },
            Fit {
                name: "young-modulus",
                range: [1.0, 298.0],
                logarithmic: false,
                scale: 1e9,
                knots: 69,
                coefficients: &[
                    205.335,
                    0.0174835,
                    -0.000365760,
                    0.000000871545,
                    -0.000000000778130,
                ],
            },
            Fit {
                name: STRAIN,
                range: [6.000000000000001, 265.0],
                logarithmic: false,
                scale: 1e-5,
                knots: 135,
                coefficients: &[-210.4, -0.05699, 0.005072, -0.00001381, 0.00000001897],
            },
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(name: &str) -> PropertyKey {
    match name {
        K => PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity)),
        CP => PropertyKey::new(name, Dims([2, 0, -2, -1, 0, 0])),
        STRAIN => PropertyKey::new(name, Dims::NONE),
        "young-modulus" => PropertyKey::new(name, Dims([-1, 1, -2, 0, 0, 0])),
        _ => panic!("unexpected property"),
    }
}
fn point(reference: &Reference, temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-processing-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
        ("source-product-form-known", 0.0),
        ("source-temper-known", reference.temper),
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
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-alloy-{}-{nonce}-{}.fsmatpk",
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
    assert_eq!(
        pack.claims().claims_ordered().count(),
        reference.fits.len() + usize::from(reference.plateau.is_some())
    );
    pack
}
fn original(fit: &Fit, t: f64) -> f64 {
    let x = if fit.logarithmic { t.log10() } else { t };
    let y = fit.coefficients.iter().rev().fold(0.0, |y, a| y * x + a);
    if fit.logarithmic {
        10.0_f64.powf(y) * fit.scale
    } else {
        y * fit.scale
    }
}
fn value(pack: &NormalizedPack, reference: &Reference, name: &str, t: f64) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &key(name),
            &point(reference, t),
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
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_g3_all_equation_knots_units_receipts_and_interpolation_are_preserved() {
    let mut knots = 0;
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let line = pack
                .claims()
                .query_line_typed(
                    &key(fit.name),
                    "temperature",
                    &point(reference, fit.range[0]),
                    &point(reference, fit.range[1]),
                    SINGLE,
                )
                .unwrap();
            assert_eq!(line.knots.len(), fit.knots);
            knots += line.knots.len();
            for (t, answer) in &line.knots {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                // Stored tables used 60-digit arithmetic. This independent f64
                // polynomial has cancellation error, distinct from source error.
                close(answer.evidence.value.value, original(fit, *t), 2e-7, 1e-15);
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
            }
            for pair in line.knots.windows(2) {
                for fraction in [0.25, 0.5, 0.75] {
                    let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                    let actual = value(&pack, reference, fit.name, t);
                    let expected = original(fit, t);
                    if fit.name == STRAIN {
                        close(actual, expected, 0.0, 5.2e-8);
                    } else {
                        close(actual, expected, 5.5e-4, 1e-15);
                    }
                }
            }
        }
    }
    assert_eq!(knots, 1404);
}

#[test]
fn g0_conditions_missing_properties_and_undefined_expansion_joins_refuse() {
    for reference in REFERENCES {
        let pack = compile(reference);
        for fit in reference.fits {
            let t = 0.5 * (fit.range[0] + fit.range[1]);
            for flag in [
                "engineering-reference-model",
                "source-processing-known",
                "source-pressure-known",
                "source-specimens-matched",
                "source-orientation-known",
                "source-product-form-known",
                "source-temper-known",
            ] {
                let p = point(reference, t);
                let changed = p.clone().with(flag, 1.0 - p.axes()[flag]).unwrap();
                assert!(
                    pack.claims()
                        .query_typed(&key(fit.name), &changed, SelectionPolicy::SingleClaimOnly)
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
            // The upper endpoint never borrows an equation-only or data-only tail.
            assert!(
                pack.claims()
                    .query_typed(
                        &key(fit.name),
                        &point(reference, fit.range[1].next_up()),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            if fit.name != STRAIN {
                assert!(
                    pack.claims()
                        .query_typed(
                            &key(fit.name),
                            &point(reference, fit.range[0].next_down()),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
        }
        let density = PropertyKey::new("density", RHO);
        assert!(
            pack.claims()
                .query_typed(
                    &density,
                    &point(reference, 250.0),
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
                    &point(reference, 250.0),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        if let Some([threshold, plateau]) = reference.plateau {
            assert_eq!(
                value(&pack, reference, STRAIN, threshold.next_down()),
                plateau
            );
            assert_eq!(value(&pack, reference, STRAIN, 4.0), plateau);
            assert!(
                pack.claims()
                    .query_typed(
                        &key(STRAIN),
                        &point(reference, threshold),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            let fit = reference.fits.iter().find(|f| f.name == STRAIN).unwrap();
            let upper = value(&pack, reference, STRAIN, threshold.next_up());
            close(upper, original(fit, threshold.next_up()), 2e-15, 1e-18);
            assert_ne!(
                upper, plateau,
                "source's rounded branch mismatch must survive"
            );
            assert!(
                pack.claims()
                    .query_line_typed(
                        &key(STRAIN),
                        "temperature",
                        &point(reference, threshold - 1.0),
                        &point(reference, threshold + 1.0),
                        SINGLE
                    )
                    .is_err()
            );
            if reference.family == "nickel-steel-fe9ni" {
                assert!(
                    pack.claims()
                        .query_typed(
                            &key(STRAIN),
                            &point(reference, 293.0),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
                assert!(
                    pack.claims()
                        .query_typed(
                            &key(CP),
                            &point(reference, 55.0),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            } else {
                assert_ne!(value(&pack, reference, STRAIN, 293.0), 0.0);
            }
        }
        if reference.family != "nickel-steel-fe9ni" {
            assert!(
                pack.claims()
                    .query_typed(
                        &key(CP),
                        &point(reference, 250.0),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}

fn card(reference: &Reference, pack: &NormalizedPack, density_model: bool) -> MaterialCard {
    let mut claims = pack.claims().clone();
    if density_model {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new("authored-unit-reference-density", RHO),
            value: PropertyValue::Scalar { value: 1.0, dims: RHO },
            validity: ValidityDomain::unconstrained().with_quantity("temperature",
                semantic(QuantityKind::AbsoluteTemperature), 250.0, 250.0),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "Authored unit-density test normalization at 250 K; not a NIST density measurement".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: reference.slug.into(), phase: "solid".into(),
        process: "source reference with explicit test-only isotropy and density normalization when requested".into(), revision: 0 },
        claims, vec![]).unwrap()
}

#[test]
fn g1_g3_sourced_alloys_drive_shared_conduction_and_bounded_free_expansion() {
    let mut powers = Vec::new();
    for reference in REFERENCES {
        let pack = compile(reference);
        let line = pack
            .claims()
            .query_line_typed(
                &key(K),
                "temperature",
                &point(reference, 100.0),
                &point(reference, 250.0),
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
        // Explicit planar, fixed-solid, scalar/aligned conductivity geometry.
        // Neither a density nor a caloric law is needed by this steady owner.
        let layer = TwoLayerConduction::new([&curve; 2], [0.01, 0.03], 0.002, 0.0).unwrap();
        let config = LayerConductionConfig {
            max_iterations: 128,
            absolute_heat_tolerance_w: 1e-9,
            relative_heat_tolerance: 1e-12,
        };
        let result = layer.evaluate([250.0, 100.0], config, || Ok(())).unwrap();
        let integral: f64 = line
            .knots
            .windows(2)
            .map(|p| {
                0.5 * (p[0].1.evidence.value.value + p[1].1.evidence.value.value)
                    * (p[1].0 - p[0].0)
            })
            .sum();
        close(result.power_a_to_b_w, 0.002 / 0.04 * integral, 2e-12, 1e-8);
        let fit = &reference.fits[0];
        // Independent Simpson integration of the ORIGINAL analytic fit checks
        // the actual material heat flow without using the stored ordinates.
        let n = 2048;
        let h = 150.0 / f64::from(n);
        let mut simpson = original(fit, 100.0) + original(fit, 250.0);
        for i in 1..n {
            simpson += if i % 2 == 0 { 2.0 } else { 4.0 } * original(fit, 100.0 + f64::from(i) * h);
        }
        close(
            result.power_a_to_b_w,
            0.002 / 0.04 * simpson * h / 3.0,
            5.5e-4,
            1e-8,
        );
        assert!(result.layer_heat_residual_w.abs() <= 1e-8);
        let reverse = layer.evaluate([100.0, 250.0], config, || Ok(())).unwrap();
        close(reverse.power_a_to_b_w, -result.power_a_to_b_w, 2e-12, 1e-8);
        let contact = TwoLayerConduction::new([&curve; 2], [0.01, 0.03], 0.002, 0.001)
            .unwrap()
            .evaluate([250.0, 100.0], config, || Ok(()))
            .unwrap();
        assert!(contact.power_a_to_b_w > 0.0 && contact.power_a_to_b_w < result.power_a_to_b_w);
        close(
            contact.interface_temperature_k[0] - contact.interface_temperature_k[1],
            0.001 * contact.power_a_to_b_w / 0.002,
            1e-12,
            1e-10,
        );
        assert!(layer.evaluate([250.0, 99.0], config, || Ok(())).is_err());
        powers.push(result.power_a_to_b_w);

        if reference.plateau.is_none() {
            continue;
        }
        let density_key = PropertyKey::new("authored-unit-reference-density", RHO);
        let strain = key(STRAIN);
        let lower = point(reference, 100.0);
        let upper = point(reference, 250.0);
        let query = FreeExpansionSourceQuery {
            density: &density_key,
            density_selection: SINGLE,
            reference: &upper,
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &strain,
                selection: SINGLE,
            },
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
        };
        assert!(
            PreparedFreeExpansionDensity::try_new(&card(reference, &pack, false), query).is_err()
        );
        let modeled = card(reference, &pack, true);
        let density = PreparedFreeExpansionDensity::try_new(&modeled, query).unwrap();
        let datum = value(&pack, reference, STRAIN, 250.0);
        for t in [100.0, 150.0, 200.0, 250.0] {
            let ratio = (1.0 + datum) / (1.0 + value(&pack, reference, STRAIN, t));
            close(
                density.at_temperature(t).unwrap().density_kg_m3(),
                ratio * ratio * ratio,
                5e-14,
                0.0,
            );
        }
        let transport = PreparedExpandingThermalTransport::try_new(
            &modeled,
            density,
            ExpandingThermalSourceQuery {
                specific_heat: &key(CP),
                conductivity: &key(K),
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selections: [SINGLE; 2],
            },
        );
        if reference.family == "nickel-steel-fe9ni" {
            let transport = transport.unwrap();
            for t in [100.0, 150.0, 200.0, 250.0] {
                let sample = transport.at_temperature(t).unwrap();
                close(
                    sample.specific_heat_j_kg_k(),
                    value(&pack, reference, CP, t),
                    5e-14,
                    0.0,
                );
                close(
                    sample.diffusivity_m2_s(),
                    value(&pack, reference, K, t)
                        / (sample.density_kg_m3() * value(&pack, reference, CP, t)),
                    5e-14,
                    0.0,
                );
            }
        } else {
            assert!(
                transport.is_err(),
                "an unmeasured heat capacity must not appear"
            );
        }
    }
    assert!(powers[0] > powers[3] && powers[3] > powers[2] && powers[2] > powers[1]);
}
