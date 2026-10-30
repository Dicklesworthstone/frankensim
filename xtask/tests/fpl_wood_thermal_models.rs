//! USDA general-wood correlations reach fixed-moisture physical storage.
//! Model evaluations are not independent experiments or species qualification.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, PropertyValue,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command, sync::OnceLock};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const MOISTURE: &str = "moisture-content-dry-basis";
const GRAVITY: &str = "specific-gravity-oven-dry-mass-current-moisture-volume";
const SLUGS: [&str; 3] = [
    "wood-dry-heat-capacity-fpl-gtr282-model",
    "wood-moist-heat-capacity-fpl-gtr282-model",
    "wood-cross-grain-conductivity-fpl-gtr282-model",
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn key(model: usize) -> PropertyKey {
    if model == 2 {
        PropertyKey::with_quantity(
            "thermal-conductivity",
            semantic(QuantityKind::ThermalConductivity),
        )
    } else {
        PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0]))
    }
}

fn point(model: usize, t: f64, m: f64, g: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    let mut context = vec![
        ("engineering-reference-model", 1.0),
        ("source-family-fpl-wood-thermal-model", 1.0),
        ("source-species-qualified", 0.0),
        ("source-pressure-known", 0.0),
        ("source-moisture-basis-dry-mass", 1.0),
        (MOISTURE, m),
    ];
    if model == 2 {
        context.extend([
            (GRAVITY, g),
            ("source-heat-flow-cross-grain", 1.0),
            (
                "specific-gravity-basis-oven-dry-mass-current-moisture-volume",
                1.0,
            ),
            ("source-temperature-nominal-24-c", 1.0),
        ]);
    } else {
        context.extend([
            ("source-below-fiber-saturation", 1.0),
            ("fixed-moisture-state", 1.0),
            (
                "capacity-mass-basis-total-moist-mass",
                if model == 1 { 1.0 } else { 0.0 },
            ),
        ]);
    }
    for (axis, value) in context {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p.with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        t,
    )
    .unwrap()
}

fn packs() -> &'static [NormalizedPack; 3] {
    static PACKS: OnceLock<[NormalizedPack; 3]> = OnceLock::new();
    PACKS.get_or_init(|| {
        SLUGS.map(|slug| {
            let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
            let nonce = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let output = std::env::temp_dir().join(format!(
                "fpl-wood-{}-{nonce}-{slug}.fsmatpk",
                std::process::id()
            ));
            let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                .args(["matdb-pack", "--manifest"])
                .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
                .arg("--out")
                .arg(&output)
                .current_dir(root)
                .env("CARGO_WORKSPACE_DIR", root)
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{slug}: {}",
                String::from_utf8_lossy(&run.stderr)
            );
            let bytes = std::fs::read(output).unwrap();
            let raw = NormalizedPack::from_bytes(&bytes).unwrap();
            let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
            assert_eq!(pack.to_bytes(), bytes);
            assert_eq!(pack.claims().claims_ordered().count(), 1);
            assert_eq!(pack.claims().observation_ids().count(), 3);
            pack
        })
    })
}

// Literal SI equations from FPL-GTR-282 Chapter 4, Equations 4-17..4-20.
// Independent f64 evaluation of the published formulas, not the decimal
// generator or a fit to the stored samples. x is percent of oven-dry mass.
fn original_cp(t: f64, moisture_dry_mass_ratio: f64) -> f64 {
    let x = 100.0 * moisture_dry_mass_ratio;
    let dry_kj_kg_k = 0.1031 + 0.003867 * t;
    let water_kj_kg_k = 4.18;
    let bond_adjustment = x * (-0.06191 + 2.36e-4 * t - 1.33e-4 * x);
    1000.0 * ((dry_kj_kg_k + water_kj_kg_k * x / 100.0) / (1.0 + x / 100.0) + bond_adjustment)
}

fn original_k(g_current_moisture_volume: f64, moisture_dry_mass_ratio: f64) -> f64 {
    g_current_moisture_volume * (0.1941 + 0.004064 * 100.0 * moisture_dry_mass_ratio) + 0.01864
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:.16e} != {expected:.16e}"
    );
}

fn value(model: usize, p: &QueryPoint) -> f64 {
    let pack = &packs()[model];
    let answer = pack
        .claims()
        .query_typed(&key(model), p, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    answer.evidence.value.value
}

#[test]
fn g1_original_handbook_equations_and_authored_surface_approximations() {
    let mut responses = 0;
    for (model, pack) in packs().iter().enumerate() {
        let (_, claim) = pack.claims().claims_ordered().next().unwrap();
        // One source correlation plus context and sampling provenance. None of
        // the hundreds of model evaluations becomes an independent experiment.
        assert_eq!(claim.observations.len(), 3);
        match &claim.value {
            PropertyValue::Curve { knots, .. } => {
                assert_eq!(model, 0);
                assert_eq!(knots.len(), 2);
                responses += knots.len();
                for t in [280.0, 291.37, 300.0, 371.6, 420.0] {
                    close(
                        value(model, &point(model, t, 0.0, 0.0)),
                        original_cp(t, 0.0),
                        2e-15,
                        1e-12,
                    );
                }
            }
            PropertyValue::Surface2D {
                axes,
                samples,
                triangles,
                ..
            } => {
                assert_eq!(axes[0], if model == 1 { "temperature" } else { GRAVITY });
                assert_eq!(axes[1], MOISTURE);
                assert_eq!(
                    (samples.len(), triangles.len()),
                    if model == 1 { (328, 560) } else { (325, 576) }
                );
                responses += samples.len();
                let equation = |x, m| {
                    if model == 1 {
                        original_cp(x, m)
                    } else {
                        original_k(x, m)
                    }
                };
                for &([x, m], y) in samples {
                    close(y, equation(x, m), 2e-15, 1e-12);
                }
                for (index, triangle) in triangles.iter().enumerate() {
                    let vertices = triangle.map(|i| samples[i as usize]);
                    for weights in [
                        [1.0 / 3.0; 3],
                        [0.5, 0.25, 0.25],
                        [0.25, 0.5, 0.25],
                        [0.25, 0.25, 0.5],
                        [0.5, 0.5, 0.0],
                        [0.5, 0.0, 0.5],
                        [0.0, 0.5, 0.5],
                    ] {
                        let sum = |component: usize| -> f64 {
                            (0..3)
                                .map(|i| {
                                    weights[i]
                                        * if component < 2 {
                                            vertices[i].0[component]
                                        } else {
                                            vertices[i].1
                                        }
                                })
                                .sum()
                        };
                        let [x, m, interpolated] = [sum(0), sum(1), sum(2)];
                        close(interpolated, equation(x, m), 6e-4, 1e-12);
                        // Exercise native point selection/receipts throughout
                        // the mesh; check every triangle's numerical reduction.
                        if index % 47 == 0 && weights[0] == 0.5 && weights[1] == 0.25 {
                            let p = if model == 1 {
                                point(model, x, m, 0.0)
                            } else {
                                point(model, 297.15, m, x)
                            };
                            close(value(model, &p), interpolated, 3e-14, 1e-12);
                        }
                    }
                }
            }
            _ => panic!("wood correlations require their authored curve or surface"),
        }
    }
    assert_eq!(responses, 655);
    // The independently printed Table 4-8 values at 300 and 340 K are rounded
    // to 0.1 kJ/(kg K), and are derived examples rather than new experiments.
    for (t, expected) in [(300.0, [1.3, 1.4, 1.7, 1.9]), (340.0, [1.4, 1.6, 1.9, 2.2])] {
        for (m, expected) in [0.0, 0.05, 0.12, 0.20].into_iter().zip(expected) {
            close(
                (original_cp(t, m) / 100.0).round() / 10.0,
                expected,
                0.0,
                1e-14,
            );
        }
    }
}

#[test]
fn g0_basis_domain_temperature_and_fixed_moisture_are_required() {
    for model in 0..3 {
        let p = point(
            model,
            if model == 2 { 297.15 } else { 300.0 },
            if model == 0 { 0.0 } else { 0.12 },
            0.51,
        );
        let claims = packs()[model].claims();
        let refuses = |wrong: &QueryPoint| {
            assert!(
                claims
                    .query_typed(&key(model), wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "model {model} admitted {wrong:?}"
            )
        };
        for (axis, v) in p
            .axes()
            .iter()
            .filter(|(axis, _)| !["temperature", MOISTURE, GRAVITY].contains(&axis.as_str()))
        {
            refuses(
                &p.clone()
                    .with(axis, if *v == 0.0 { 1.0 } else { 0.0 })
                    .unwrap(),
            );
        }
        refuses(
            &p.clone()
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    300.0,
                )
                .unwrap(),
        );
        // A Celsius numerical argument is not a Kelvin coordinate; percent 12
        // is not the canonical dry-basis ratio 0.12. Neither is silently fixed.
        refuses(&p.clone().with("temperature", 20.0).unwrap());
        refuses(&p.clone().with(MOISTURE, 12.0).unwrap());
        refuses(&p.clone().with(MOISTURE, -0.001).unwrap());
        refuses(&QueryPoint::new().with("temperature", 300.0).unwrap());
        let bounds: [f64; 2] = if model == 2 {
            [297.15, 297.15]
        } else {
            [280.0, 420.0]
        };
        for t in [bounds[0].next_down(), bounds[1].next_up()] {
            refuses(&p.clone().with("temperature", t).unwrap());
        }
        let max_m: f64 = [0.0, 0.20, 0.24][model];
        refuses(&p.clone().with(MOISTURE, max_m.next_up()).unwrap());
        if model == 2 {
            for g in [0.30, 0.31_f64.next_down(), 0.78_f64.next_up()] {
                refuses(&p.clone().with(GRAVITY, g).unwrap());
            }
        }
        for absent in [
            "density",
            "volumetric-heat-capacity",
            "heat-of-sorption",
            "latent-heat",
            "moisture-diffusivity",
        ] {
            assert!(
                claims
                    .query(absent, &p, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
    }
    assert!(
        packs()[1]
            .claims()
            .query_line_typed(
                &key(1),
                "temperature",
                &point(1, 280.0, 0.05, 0.0),
                &point(1, 420.0, 0.12, 0.0),
                SINGLE
            )
            .is_err()
    );
    // The conductivity nominal point cannot become a temperature curve.
    assert!(
        packs()[2]
            .claims()
            .query_line_typed(
                &key(2),
                "temperature",
                &point(2, 297.15, 0.12, 0.51),
                &point(2, 300.0, 0.12, 0.51),
                SINGLE
            )
            .is_err()
    );
    let card = card(1);
    let capacity = key(1);
    let lower = point(1, 280.0, 0.12, 0.0);
    let upper = point(1, 420.0, 0.12, 0.0);
    assert!(
        PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
            id: 1,
            card: &card,
            capacity: &capacity,
            amount: HeatStorageAmount::ReferenceVolumeM3(0.01),
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selection: SINGLE,
        }])
        .is_err()
    );
}

fn card(model: usize) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: SLUGS[model].into(),
            phase: "wood-general-empirical-model".into(),
            process: "authored fixed-moisture sensible storage; no drying or species qualification"
                .into(),
            revision: 0,
        },
        packs()[model].claims().clone(),
        vec![],
    )
    .unwrap()
}

fn storage(model: usize, m: f64, total_mass_kg: f64) -> PreparedHeatStorage {
    let card = card(model);
    let capacity = key(model);
    let lower = point(model, 280.0, m, 0.0);
    let upper = point(model, 420.0, m, 0.0);
    PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 1,
        card: &card,
        capacity: &capacity,
        amount: HeatStorageAmount::MassKg(total_mass_kg),
        temperature_axis: "temperature",
        lower: &lower,
        upper: &upper,
        selection: SINGLE,
    }])
    .unwrap()
}

fn analytic_storage(m: f64, mass: f64, lo: f64, hi: f64, entropy: bool) -> f64 {
    // Algebraic integration of the literal source equation at fixed moisture,
    // independent of the production piecewise-affine integral and inverse.
    let x = 100.0 * m;
    let a = 1000.0 * ((0.1031 + 4.18 * m) / (1.0 + m) + x * (-0.06191 - 1.33e-4 * x));
    let b = 1000.0 * (0.003867 / (1.0 + m) + x * 2.36e-4);
    mass * if entropy {
        a * (hi / lo).ln() + b * (hi - lo)
    } else {
        a * (hi - lo) + 0.5 * b * (hi * hi - lo * lo)
    }
}

#[test]
fn g1_g3_dry_and_retained_moisture_heating_entropy_and_inverse() {
    let dry_mass_kg = 2.0;
    let dry = storage(0, 0.0, dry_mass_kg);
    let fixed_heat = dry.energy_change_j(290.0, 350.0).unwrap();
    for (model, m) in [(0, 0.0), (1, 0.0), (1, 0.05), (1, 0.12), (1, 0.20)] {
        // Same dry wood quantity; moisture adds retained water and changes Cp.
        // The production mass is total moist mass, not its oven-dry reference.
        let mass = dry_mass_kg * (1.0 + m);
        let body = storage(model, m, mass);
        assert_eq!(body.temperature_span_k(), [280.0, 420.0]);
        close(
            body.heat_capacity_j_k(319.75).unwrap(),
            mass * original_cp(319.75, m),
            3e-14,
            1e-10,
        );
        for [lo, hi] in [[290.0, 350.0], [283.5, 410.25], [410.25, 283.5]] {
            let energy = analytic_storage(m, mass, lo, hi, false);
            close(body.energy_change_j(lo, hi).unwrap(), energy, 3e-14, 1e-8);
            close(
                body.sensible_entropy_change_j_k(lo, hi).unwrap(),
                analytic_storage(m, mass, lo, hi, true),
                3e-13,
                1e-9,
            );
            close(
                body.temperature_after_energy_change(lo, energy).unwrap(),
                hi,
                0.0,
                2e-10,
            );
        }
        let final_t = body
            .temperature_after_energy_change(290.0, fixed_heat)
            .unwrap();
        if m > 0.0 {
            assert!(final_t > 290.0 && final_t < 350.0);
        } else {
            close(final_t, 350.0, 0.0, 2e-10);
        }
        assert!(body.heat_capacity_j_k(280.0_f64.next_down()).is_err());
        assert!(body.temperature_after_energy_change(419.0, 1e9).is_err());
    }
    // A moisture between stored rows follows the declared triangulated model.
    // Its storage integral and entropy track the original equation within the
    // independently checked numerical approximation, without a new evaluator.
    let m = 0.123;
    let mass = dry_mass_kg * (1.0 + m);
    let body = storage(1, m, mass);
    let energy = body.energy_change_j(291.3, 401.7).unwrap();
    close(
        energy,
        analytic_storage(m, mass, 291.3, 401.7, false),
        6e-4,
        1e-8,
    );
    close(
        body.sensible_entropy_change_j_k(291.3, 401.7).unwrap(),
        analytic_storage(m, mass, 291.3, 401.7, true),
        6e-4,
        1e-9,
    );
    close(
        body.temperature_after_energy_change(291.3, energy).unwrap(),
        401.7,
        0.0,
        2e-9,
    );
    assert_ne!(body.identity(), storage(1, m, 2.0 * mass).identity());
}
