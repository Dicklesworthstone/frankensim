//! Primary NIST foam columns remain separate and drive bounded insulation flux.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, TwoLayerConduction,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Foam {
    slug: &'static str,
    family: &'static str,
    gas_axis: &'static str,
    gas_value: f64,
    density_label: f64,
    lower_k: f64,
    knots: usize,
    coefficients: &'static [f64],
}

// Literal primary-page coefficients transcribed independently of the 70-digit
// decimal evaluator used to prepare the committed ordinates. Density is solely
// the conductivity-column label; no density(T) or Cp has been attached to it.
const FOAMS: &[Foam] = &[
    Foam {
        slug: "polyurethane-freon-31p88-nist-thermal-reference",
        family: "polyurethane-freon",
        gas_axis: "source-cell-gas-freon",
        gas_value: 1.0,
        density_label: 31.88,
        lower_k: 76.0,
        knots: 77,
        coefficients: &[
            -3218.679, 9201.61, -10956.66, 6950.102, -2476.94, 470.284, -37.1669, 0.0, 0.0,
        ],
    },
    Foam {
        slug: "polyurethane-co2-32p04-nist-thermal-reference",
        family: "polyurethane-co2",
        gas_axis: "source-cell-gas-co2",
        gas_value: 1.0,
        density_label: 32.04,
        lower_k: 100.0,
        knots: 35,
        coefficients: &[
            3788.48, -7642.66, 4592.448, 778.8423, -2214.434, 1090.293, -235.6349, 19.66088, 0.0,
        ],
    },
    Foam {
        slug: "polyurethane-freon-64p08-nist-thermal-reference",
        family: "polyurethane-freon",
        gas_axis: "source-cell-gas-freon",
        gas_value: 1.0,
        density_label: 64.08,
        lower_k: 88.0,
        knots: 69,
        coefficients: &[
            789.79, -2347.94, 3024.61, -2206.78, 989.238, -273.18, 43.065, -2.9863, 0.0,
        ],
    },
    Foam {
        slug: "polystyrene-freon-31p88-nist-thermal-reference",
        family: "polystyrene-freon",
        gas_axis: "source-cell-gas-freon",
        gas_value: 1.0,
        density_label: 31.88,
        lower_k: 90.0,
        knots: 61,
        coefficients: &[
            -1557.5, 3984.7, -3940.245, 1649.668, 12.80097, -294.2616, 119.5898, -20.6301, 1.36067,
        ],
    },
    Foam {
        slug: "polystyrene-32p04-nist-thermal-reference",
        family: "polystyrene",
        gas_axis: "source-cell-gas-known",
        gas_value: 0.0,
        density_label: 32.04,
        lower_k: 33.0,
        knots: 52,
        coefficients: &[
            -1145.45, 4086.02, -6234.293, 5260.106, -2649.914, 797.0372, -132.5244, 9.3968, 0.0,
        ],
    },
    Foam {
        slug: "polystyrene-49p98-nist-thermal-reference",
        family: "polystyrene",
        gas_axis: "source-cell-gas-known",
        gas_value: 0.0,
        density_label: 49.98,
        lower_k: 7.0,
        knots: 100,
        coefficients: &[
            -1.5194, -4.6449, 11.643, -15.969, 12.722, -5.821, 1.4174, -0.14128, 0.0,
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal-conductivity",
        semantic(QuantityKind::ThermalConductivity),
    )
}

fn point(foam: &Foam, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-grade-known", 0.0),
        ("source-processing-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
        ("source-cell-gas-composition-known", 0.0),
        (foam.gas_axis, foam.gas_value),
    ] {
        p = p
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p.with_quantity(
        format!("source-family-nist-{}", foam.family),
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
    )
    .unwrap()
    .with_quantity(
        "source-reference-density",
        QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
        foam.density_label,
    )
    .unwrap()
    .with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        t,
    )
    .unwrap()
}

fn compile(foam: &Foam) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-foam-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        foam.slug
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", foam.slug)))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        foam.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), 2);
    pack
}

fn original(foam: &Foam, t: f64) -> f64 {
    let x = t.log10();
    let exponent = foam.coefficients.iter().rev().fold(0.0, |y, a| y * x + a);
    10.0_f64.powf(exponent)
}

fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs() + absolute,
        "{a:.16e} != {b:.16e}"
    );
}

fn original_integral(foam: &Foam, lo: f64, hi: f64) -> f64 {
    // Independent composite Simpson quadrature of the source equation, not
    // the production chart integrator or its stored ordinate samples.
    let n = 2048;
    let step = (hi - lo) / f64::from(n);
    let mut total = original(foam, lo) + original(foam, hi);
    for i in 1..n {
        total += if i % 2 == 0 { 2.0 } else { 4.0 } * original(foam, lo + f64::from(i) * step);
    }
    total * step / 3.0
}

fn curve(pack: &NormalizedPack, foam: &Foam) -> ConductivityCurve {
    let line = pack
        .claims()
        .query_line_typed(
            &key(),
            "temperature",
            &point(foam, foam.lower_k),
            &point(foam, 300.0),
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
fn g1_g3_primary_coefficients_units_and_column_metadata_survive_compilation() {
    let mut numerical_values = 0;
    for foam in FOAMS {
        let pack = compile(foam);
        let line = pack
            .claims()
            .query_line_typed(
                &key(),
                "temperature",
                &point(foam, foam.lower_k),
                &point(foam, 300.0),
                SINGLE,
            )
            .unwrap();
        assert_eq!(line.knots.len(), foam.knots);
        numerical_values += line.knots.len() + 1;
        for (t, answer) in &line.knots {
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            // The tolerance accounts for cancellation in this independent f64
            // evaluation of large rounded coefficients; it is not uncertainty.
            close(answer.evidence.value.value, original(foam, *t), 2e-7, 1e-16);
        }
        for pair in line.knots.windows(2) {
            for fraction in [0.25, 0.5, 0.75] {
                let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                let a = pack
                    .claims()
                    .query_typed(&key(), &point(foam, t), SelectionPolicy::SingleClaimOnly)
                    .unwrap();
                close(a.evidence.value.value, original(foam, t), 8.1e-5, 1e-16);
            }
        }
        let density_key = PropertyKey::new("nominal-source-density", Dims([-3, 1, 0, 0, 0, 0]));
        let density_point = point(foam, 150.0)
            .with_quantity(
                "source-density-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        let density = pack
            .claims()
            .query_typed(
                &density_key,
                &density_point,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(density.evidence.value.value, foam.density_label);
        pack.claims().verify_receipt(&density.receipt).unwrap();
        assert!(
            pack.claims()
                .query_typed(
                    &density_key,
                    &density_point
                        .with("source-density-temperature-known", 1.0)
                        .unwrap(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
    assert_eq!(numerical_values, 400);
}

#[test]
fn g0_gas_density_conditions_and_missing_bulk_laws_are_not_inferred() {
    for foam in FOAMS {
        let pack = compile(foam);
        let p = point(foam, 150.0);
        for axis in [
            "engineering-reference-model",
            "source-grade-known",
            "source-processing-known",
            "source-moisture-condition-known",
            "source-pressure-known",
            "source-specimens-matched",
            "source-orientation-known",
            "source-cell-gas-composition-known",
            foam.gas_axis,
        ] {
            let wrong = p.clone().with(axis, 1.0 - p.axes()[axis]).unwrap();
            assert!(
                pack.claims()
                    .query_typed(&key(), &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        let wrong_density = p
            .clone()
            .with("source-reference-density", foam.density_label + 0.01)
            .unwrap();
        assert!(
            pack.claims()
                .query_typed(&key(), &wrong_density, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        for t in [foam.lower_k.next_down(), 300.0_f64.next_up()] {
            assert!(
                pack.claims()
                    .query_typed(&key(), &point(foam, t), SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        assert!(
            pack.claims()
                .query_typed(&key(), &QueryPoint::new(), SelectionPolicy::SingleClaimOnly)
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
                .query_typed(&key(), &wrong_temperature, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        for missing in [
            "density",
            "specific-heat-capacity",
            "volumetric-heat-capacity",
            "linear-thermal-expansion-coefficient",
            "young-modulus",
            "porosity",
        ] {
            assert!(pack.claims().claims_for(missing).is_empty());
        }
    }
}

#[test]
fn g1_g3_real_foam_curves_change_shared_insulation_conduction() {
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    };
    let curves: Vec<_> = FOAMS.iter().map(|f| curve(&compile(f), f)).collect();
    let mut powers = Vec::new();
    for (foam, k) in FOAMS.iter().zip(&curves) {
        // A scalar fixed-direction homogeneous-layer interpretation is an
        // authored test-model choice, not evidence of isotropic foam response.
        let model = TwoLayerConduction::new([k, k], [0.025, 0.04], 0.2, 0.0).unwrap();
        let solution = model.evaluate([300.0, 125.0], config, || Ok(())).unwrap();
        close(
            solution.power_a_to_b_w,
            0.2 / 0.065 * original_integral(foam, 125.0, 300.0),
            8.1e-5,
            1e-10,
        );
        assert!(solution.layer_heat_residual_w.abs() < 1e-9);
        let reverse = model.evaluate([125.0, 300.0], config, || Ok(())).unwrap();
        close(
            reverse.power_a_to_b_w,
            -solution.power_a_to_b_w,
            1e-12,
            1e-10,
        );
        assert!(model.evaluate([300.01, 125.0], config, || Ok(())).is_err());
        powers.push(solution.power_a_to_b_w);
    }
    // The columns represent distinct source cases, not a controlled study
    // isolating gas or density. Their observed transport difference is retained.
    assert!((powers[0] - powers[1]).abs() > 1.0);
    assert!((powers[4] - powers[5]).abs() > 1.0);

    let mixed =
        TwoLayerConduction::new([&curves[0], &curves[1]], [0.025, 0.04], 0.2, 0.12).unwrap();
    let solution = mixed.evaluate([300.0, 125.0], config, || Ok(())).unwrap();
    let [ta, tb] = solution.interface_temperature_k;
    close(
        solution.power_a_to_b_w,
        0.2 / 0.025 * original_integral(&FOAMS[0], ta, 300.0),
        8.1e-5,
        1e-9,
    );
    close(
        solution.power_a_to_b_w,
        0.2 / 0.04 * original_integral(&FOAMS[1], 125.0, tb),
        8.1e-5,
        1e-9,
    );
    close(ta - tb, solution.power_a_to_b_w * 0.12 / 0.2, 1e-11, 1e-10);
    assert!(ta > tb);
    assert!(solution.layer_heat_residual_w.abs() < 1e-9);
}
