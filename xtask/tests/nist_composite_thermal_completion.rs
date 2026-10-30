//! Primary composite/nonmetal fits reach bounded conduction and mass storage.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::storage::{HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage},
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Reference {
    slug: &'static str,
    family: &'static str,
    range: [f64; 2],
    knots: usize,
    coefficients: &'static [f64],
    extra: &'static [(&'static str, f64)],
    density_lb_ft3: Option<f64>,
    capacity: bool,
}

// Literal primary-page coefficients, independent of the 70-digit decimal
// evaluator that generated the stored ordinates. Six coefficients select the
// published Kevlar erf blend; nine select the published log10 polynomial.
const REFERENCES: &[Reference] = &[
    Reference {
        slug: "kevlar49-fiber-nist-thermal-reference",
        family: "kevlar49-fiber",
        range: [0.1, 291.0],
        knots: 392,
        coefficients: &[-2.4219, 1.986637, 1.257441, 0.961209, -9.6106, 0.777857],
        extra: &[
            ("source-fiber-volume-fraction-known", 0.0),
            ("source-fiber-sizing-known", 0.0),
        ],
        density_lb_ft3: None,
        capacity: false,
    },
    Reference {
        slug: "kevlar49-composite-nist-thermal-reference",
        family: "kevlar49-composite",
        range: [6.0, 302.0],
        knots: 128,
        coefficients: &[-2.65, 1.986637, 1.24851, 0.57, -8.0, 0.777857],
        extra: &[
            ("source-fiber-volume-fraction-known", 0.0),
            ("source-fiber-sizing-known", 0.0),
            ("source-cure-known", 0.0),
            ("source-resin-formulation-known", 0.0),
            ("source-layup-known", 0.0),
        ],
        density_lb_ft3: None,
        capacity: false,
    },
    Reference {
        slug: "balsa-6lb-ft3-nist-thermal-reference",
        family: "balsa",
        range: [88.0, 300.0],
        knots: 26,
        coefficients: &[
            4172.447, -11309.97, 12745.09, -7647.584, 2577.309, -462.538, 34.5351, 0.0, 0.0,
        ],
        extra: &[
            ("source-species-known", 0.0),
            ("source-moisture-basis-known", 0.0),
            ("source-density-column-lb-ft3", 6.0),
        ],
        density_lb_ft3: Some(6.0),
        capacity: false,
    },
    Reference {
        slug: "balsa-11lb-ft3-nist-thermal-reference",
        family: "balsa",
        range: [88.0, 300.0],
        knots: 42,
        coefficients: &[
            4815.4, -12969.63, 14520.76, -8654.164, 2895.712, -515.7272, 38.19218, 0.0, 0.0,
        ],
        extra: &[
            ("source-species-known", 0.0),
            ("source-moisture-basis-known", 0.0),
            ("source-density-column-lb-ft3", 11.0),
        ],
        density_lb_ft3: Some(11.0),
        capacity: false,
    },
    Reference {
        slug: "apiezon-n-nist-thermal-reference",
        family: "apiezon-n",
        range: [2.0, 200.0],
        knots: 319,
        coefficients: &[
            -1.61975, 3.10923, -0.712719, 4.93675, -9.37993, 7.58304, -3.11048, 0.628309,
            -0.0482634,
        ],
        extra: &[
            ("source-formulation-known", 0.0),
            ("source-lot-known", 0.0),
            ("source-caloric-constraint-known", 0.0),
        ],
        density_lb_ft3: None,
        capacity: true,
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn key(reference: &Reference) -> PropertyKey {
    if reference.capacity {
        PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0]))
    } else {
        PropertyKey::with_quantity(
            "thermal-conductivity",
            semantic(QuantityKind::ThermalConductivity),
        )
    }
}

fn point(reference: &Reference, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (axis, value) in [
        ("engineering-reference-model", 1.0),
        ("source-processing-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-orientation-known", 0.0),
        ("source-specimens-matched", 0.0),
    ]
    .iter()
    .chain(reference.extra)
    {
        p = p
            .with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), *value)
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
        "nist-composite-{}-{nonce}-{}.fsmatpk",
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
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(
        pack.claims().claims_ordered().count(),
        1 + usize::from(reference.density_lb_ft3.is_some())
    );
    pack
}

fn original(reference: &Reference, t: f64) -> f64 {
    let x = t.log10();
    let exponent = match reference.coefficients {
        [a, b, c, d, e, f] => {
            let weight = 0.5 * (1.0 + fs_math::det::erf(2.0 * (x - c)));
            (a + b * x) * (1.0 - weight) + (d + e * (-x / f).exp()) * weight
        }
        polynomial => polynomial.iter().rev().fold(0.0, |y, a| y * x + a),
    };
    10.0_f64.powf(exponent)
}

fn original_integral(reference: &Reference, lo: f64, hi: f64, entropy: bool) -> f64 {
    // Composite Simpson quadrature of the original equation, independent of
    // the authored samples and the production piecewise-linear integrator.
    let value = |t| original(reference, t) / if entropy { t } else { 1.0 };
    let n = 2048;
    let step = (hi - lo) / f64::from(n);
    let mut total = value(lo) + value(hi);
    for i in 1..n {
        total += if i % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(i) * step);
    }
    total * step / 3.0
}

fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs() + absolute,
        "{a:.16e} != {b:.16e}"
    );
}

fn curve(pack: &NormalizedPack, reference: &Reference) -> ConductivityCurve {
    assert!(!reference.capacity);
    let line = pack
        .claims()
        .query_line_typed(
            &key(reference),
            "temperature",
            &point(reference, reference.range[0]),
            &point(reference, reference.range[1]),
            SINGLE,
        )
        .unwrap();
    ConductivityCurve::new(
        line.knots
            .iter()
            .map(|(t, a)| [*t, a.evidence.value.value])
            .collect(),
    )
    .unwrap()
}

#[test]
fn g1_g3_primary_equations_units_receipts_and_density_labels_survive_compilation() {
    let mut numerical_values = 0;
    for reference in REFERENCES {
        let pack = compile(reference);
        let line = pack
            .claims()
            .query_line_typed(
                &key(reference),
                "temperature",
                &point(reference, reference.range[0]),
                &point(reference, reference.range[1]),
                SINGLE,
            )
            .unwrap();
        assert_eq!(line.knots.len(), reference.knots);
        numerical_values += line.knots.len();
        for (t, answer) in &line.knots {
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            // Allows cancellation in independent f64 polynomial evaluation;
            // this numerical tolerance is not experimental uncertainty.
            close(
                answer.evidence.value.value,
                original(reference, *t),
                2e-7,
                1e-16,
            );
        }
        for pair in line.knots.windows(2) {
            for fraction in [0.25, 0.5, 0.75] {
                let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                let answer = pack
                    .claims()
                    .query_typed(
                        &key(reference),
                        &point(reference, t),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap();
                close(
                    answer.evidence.value.value,
                    original(reference, t),
                    8.1e-5,
                    1e-16,
                );
            }
        }
        if let Some(lb_ft3) = reference.density_lb_ft3 {
            let density_key = PropertyKey::new("nominal-source-density", Dims([-3, 1, 0, 0, 0, 0]));
            let p = point(reference, 150.0)
                .with_quantity(
                    "source-density-temperature-known",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                )
                .unwrap();
            let answer = pack
                .claims()
                .query_typed(&density_key, &p, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(
                answer.evidence.value.value,
                lb_ft3 * 0.45359237 / 0.3048_f64.powi(3),
                2e-15,
                0.0,
            );
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(
                pack.claims()
                    .query_typed(
                        &density_key,
                        &p.with("source-density-temperature-known", 1.0).unwrap(),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            numerical_values += 1;
        }
    }
    assert_eq!(numerical_values, 909);
}

#[test]
fn g0_cure_orientation_caloric_context_and_unsupported_bulk_laws_are_not_inferred() {
    for reference in REFERENCES {
        let pack = compile(reference);
        let p = point(reference, 150.0);
        for (axis, value) in p
            .axes()
            .iter()
            .filter(|(axis, _)| axis.as_str() != "temperature")
        {
            let wrong = p
                .clone()
                .with(axis, if *value == 0.0 { 1.0 } else { 0.0 })
                .unwrap();
            assert!(
                pack.claims()
                    .query_typed(&key(reference), &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "{} accepted {axis}",
                reference.slug
            );
        }
        for t in [reference.range[0].next_down(), reference.range[1].next_up()] {
            assert!(
                pack.claims()
                    .query_typed(
                        &key(reference),
                        &point(reference, t),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        assert!(
            pack.claims()
                .query_typed(
                    &key(reference),
                    &QueryPoint::new(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        let wrong_temperature = p
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                150.0,
            )
            .unwrap();
        assert!(
            pack.claims()
                .query_typed(
                    &key(reference),
                    &wrong_temperature,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        for absent in [
            "density",
            "volumetric-heat-capacity",
            "linear-thermal-expansion-coefficient",
            "young-modulus",
            "thermal-contact-resistance",
            if reference.capacity {
                "thermal-conductivity"
            } else {
                "specific-heat-capacity"
            },
        ] {
            assert!(
                pack.claims().claims_for(absent).is_empty(),
                "{} supplied {absent}",
                reference.slug
            );
        }
    }
    let fiber = compile(&REFERENCES[0]);
    assert!(
        fiber
            .claims()
            .query_typed(
                &key(&REFERENCES[0]),
                &point(&REFERENCES[1], 150.0),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let balsa = compile(&REFERENCES[2]);
    assert!(
        balsa
            .claims()
            .query_typed(
                &key(&REFERENCES[2]),
                &point(&REFERENCES[3], 150.0),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn g1_g3_distinct_fiber_composite_and_wood_transport_and_grease_mass_heating() {
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    };
    let curves: Vec<_> = REFERENCES[..4]
        .iter()
        .map(|r| curve(&compile(r), r))
        .collect();
    let mut powers = Vec::new();
    for (reference, k) in REFERENCES[..4].iter().zip(&curves) {
        // This fixture explicitly authors scalar fixed-direction layers; it
        // does not assert source isotropy or a known balsa grain direction.
        let model = TwoLayerConduction::new([k, k], [0.03, 0.07], 0.04, 0.0).unwrap();
        let solution = model.evaluate([200.0, 100.0], config, || Ok(())).unwrap();
        close(
            solution.power_a_to_b_w,
            0.04 / 0.1 * original_integral(reference, 100.0, 200.0, false),
            8.1e-5,
            1e-10,
        );
        close(
            solution.interface_temperature_k[0],
            solution.interface_temperature_k[1],
            0.0,
            1e-9,
        );
        let reverse = model.evaluate([100.0, 200.0], config, || Ok(())).unwrap();
        close(
            reverse.power_a_to_b_w,
            -solution.power_a_to_b_w,
            1e-12,
            1e-9,
        );
        assert!(
            model
                .evaluate([reference.range[1].next_up(), 100.0], config, || Ok(()))
                .is_err()
        );
        powers.push(solution.power_a_to_b_w);
    }
    assert!((powers[0] - powers[1]).abs() > 1.0);
    assert!((powers[2] - powers[3]).abs() > 0.1);
    let mixed =
        TwoLayerConduction::new([&curves[0], &curves[1]], [0.03, 0.07], 0.04, 0.12).unwrap();
    let solution = mixed.evaluate([200.0, 100.0], config, || Ok(())).unwrap();
    let [ta, tb] = solution.interface_temperature_k;
    close(
        solution.power_a_to_b_w,
        0.04 / 0.03 * original_integral(&REFERENCES[0], ta, 200.0, false),
        8.1e-5,
        1e-9,
    );
    close(
        solution.power_a_to_b_w,
        0.04 / 0.07 * original_integral(&REFERENCES[1], 100.0, tb, false),
        8.1e-5,
        1e-9,
    );
    close(ta - tb, solution.power_a_to_b_w * 0.12 / 0.04, 1e-11, 1e-9);
    assert!(solution.layer_heat_residual_w.abs() < 1e-9);

    let grease = &REFERENCES[4];
    let pack = compile(grease);
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: grease.slug.into(),
            phase: "source-phase-unspecified".into(),
            process: "NIST specific heat; supplied mass; authored sensible-storage interpretation"
                .into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )
    .unwrap();
    let capacity = key(grease);
    let lower = point(grease, 2.0);
    let upper = point(grease, 200.0);
    let query = |mass| HeatStorageComponentQuery {
        id: 1,
        card: &card,
        capacity: &capacity,
        amount: HeatStorageAmount::MassKg(mass),
        temperature_axis: "temperature",
        lower: &lower,
        upper: &upper,
        selection: SINGLE,
    };
    // The source does not identify Cp versus Cv. This explicitly authored
    // sensible-storage model retains that unknown; it does not qualify Cp.
    let storage = PreparedHeatStorage::try_new(&[query(0.01)]).unwrap();
    assert_eq!(storage.temperature_span_k(), [2.0, 200.0]);
    close(
        storage.heat_capacity_j_k(100.0).unwrap(),
        0.01 * original(grease, 100.0),
        2e-7,
        1e-10,
    );
    let energy = storage.energy_change_j(100.0, 150.0).unwrap();
    close(
        energy,
        0.01 * original_integral(grease, 100.0, 150.0, false),
        8.1e-5,
        1e-8,
    );
    close(
        storage.energy_change_j(150.0, 100.0).unwrap(),
        -energy,
        1e-13,
        1e-9,
    );
    close(
        storage
            .temperature_after_energy_change(100.0, energy)
            .unwrap(),
        150.0,
        0.0,
        2e-9,
    );
    close(
        storage.sensible_entropy_change_j_k(100.0, 150.0).unwrap(),
        0.01 * original_integral(grease, 100.0, 150.0, true),
        8.1e-5,
        1e-9,
    );
    assert!(storage.heat_capacity_j_k(1.0).is_err());
    assert!(storage.temperature_after_energy_change(150.0, 1e9).is_err());
    let doubled = PreparedHeatStorage::try_new(&[query(0.02)]).unwrap();
    let final_temperature = doubled
        .temperature_after_energy_change(100.0, energy)
        .unwrap();
    assert!(final_temperature > 100.0 && final_temperature < 150.0);
    assert_ne!(storage.identity(), doubled.identity());
}
