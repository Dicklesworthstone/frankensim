//! G0/G1/G3: cured epoxy calorimetry and qualified producer adhesive data.
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
const EA: &str = "henkel-loctite-ea9394-cured-2015";
const ARALDITE: &str = "huntsman-araldite-2011-cured-reference";
const CP: &str = "specific-heat-capacity";
const K: &str = "thermal-conductivity";
const NOMINAL_K: &str = "nominal-source-thermal-conductivity";
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

struct Epoxy {
    slug: &'static str,
    grade: &'static str,
    ratio: f64,
    nominal_density: f64,
    knots: usize,
    coefficients: [f64; 5],
    figure_ordinate_at_5k: [f64; 2],
}
// Independently transcribed from the rendered original Table 2. Its PDF text
// encoding corrupts some digits. Fig. 1 uses J/(g K^4) for C/T^3, whereas the
// unlabeled polynomial output requires the explicit J/(kg K) interpretation.
const EPOXIES: &[Epoxy] = &[
    Epoxy {
        slug: "stycast-2850ft-23lv-rostem-2021-caloric",
        grade: "stycast-2850ft-23lv",
        ratio: 100.0 / 7.5,
        nominal_density: 2400.0,
        knots: 252,
        coefficients: [-2.11494050, 2.51709056, 2.59374465, -2.99675188, 0.86125486],
        figure_ordinate_at_5k: [0.9e-5, 1.1e-5],
    },
    Epoxy {
        slug: "stycast-1266-rostem-2021-caloric",
        grade: "stycast-1266",
        ratio: 100.0 / 28.0,
        nominal_density: 1100.0,
        knots: 256,
        coefficients: [-1.72272035, 2.74684260, 2.20865182, -2.62360390, 0.73084529],
        figure_ordinate_at_5k: [2.7e-5, 3.1e-5],
    },
    Epoxy {
        slug: "epotek-377-rostem-2021-caloric",
        grade: "epotek-377",
        ratio: 1.0,
        nominal_density: 1200.0,
        knots: 245,
        coefficients: [-1.73020906, 3.65343010, 0.30871953, -1.22590643, 0.37382964],
        figure_ordinate_at_5k: [3.3e-5, 3.7e-5],
    },
    Epoxy {
        slug: "epotek-377h-rostem-2021-caloric",
        grade: "epotek-377h",
        ratio: 1.0,
        nominal_density: 1300.0,
        knots: 245,
        coefficients: [-1.77036189, 3.21084905, 1.13701708, -1.91294933, 0.58275509],
        figure_ordinate_at_5k: [2.4e-5, 2.8e-5],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn quantity(p: QueryPoint, axis: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(axis, QuantitySpec::dimensional(dims), value)
        .unwrap()
}
fn flag(p: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    quantity(p, axis, Dims::NONE, value)
}
fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values
        .iter()
        .fold(QueryPoint::new(), |p, (name, v)| flag(p, name, *v))
}
fn temperature(p: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    p.with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), value)
        .unwrap()
}
fn key(name: &str) -> PropertyKey {
    if name == K || name == NOMINAL_K {
        return PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity));
    }
    PropertyKey::new(
        name,
        match name {
            CP => Dims([2, 0, -2, -1, 0, 0]),
            "nominal-source-density" | "density" => RHO,
            "linear-thermal-expansion-coefficient" => Dims([0, 0, 0, -1, 0, 0]),
            "tensile-elongation-at-break" => Dims::NONE,
            _ => PA,
        },
    )
}
fn epoxy_base(e: &Epoxy) -> QueryPoint {
    flag(
        flags(&[
            ("resin-to-hardener-mass-ratio", e.ratio),
            ("source-cure-history-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-specimens-matched", 0.0),
        ]),
        &format!("source-grade-{}", e.grade),
        1.0,
    )
}
fn epoxy_point(e: &Epoxy, t: f64) -> QueryPoint {
    let mut p = epoxy_base(e);
    for (axis, value) in [
        ("engineering-reference-model", 1.0),
        ("source-measurement-vacuum", 1.0),
        ("source-pressure-known", 0.0),
        ("source-caloric-constraint-known", 0.0),
        ("source-fit-unit-labelled", 0.0),
        ("source-fit-unit-interpretation-j-per-kg-k", 1.0),
    ] {
        p = flag(p, axis, value);
    }
    temperature(p, "temperature", t)
}
fn ea_base(cure_known: f64) -> QueryPoint {
    flags(&[
        ("source-grade-ea9394", 1.0),
        ("resin-to-hardener-mass-ratio", 100.0 / 17.0),
        ("source-moisture-conditioning-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-test-direction-known", 0.0),
        ("source-cure-history-known", cure_known),
    ])
}
fn ea_at(t: f64, cure_known: f64) -> QueryPoint {
    temperature(
        flag(ea_base(cure_known), "source-test-temperature-known", 1.0),
        "temperature",
        t,
    )
}
fn ea_nominal_k() -> QueryPoint {
    flag(
        flag(ea_base(0.0), "source-test-temperature-known", 0.0),
        "source-test-method-known",
        0.0,
    )
}
fn ea_tension() -> QueryPoint {
    quantity(
        flag(ea_at(298.15, 0.0), "source-test-rate-known", 0.0),
        "casting-thickness",
        LENGTH,
        0.00318,
    )
}
fn ea_shear(moisture: f64) -> QueryPoint {
    flag(
        flag(
            flag(ea_at(298.15, 0.0), "source-test-rate-known", 0.0),
            "source-test-method-known",
            0.0,
        ),
        "source-moisture-label",
        moisture,
    )
}
fn ea_compression() -> QueryPoint {
    let mut p = flag(ea_at(298.15, 0.0), "source-test-rate-known", 0.0);
    for (axis, value) in [
        ("specimen-width", 0.0127),
        ("specimen-length", 0.0254),
        ("specimen-height", 0.0127),
    ] {
        p = quantity(p, axis, LENGTH, value);
    }
    p
}
fn ea_expansion(t: f64) -> QueryPoint {
    flag(
        flag(
            ea_at(t, 0.0),
            "source-coefficient-length-convention-known",
            0.0,
        ),
        "source-thermal-scan-conditions-known",
        0.0,
    )
}
fn ea_joint(t: f64) -> QueryPoint {
    let mut p = temperature(ea_at(t, 1.0), "cure-temperature", 298.15);
    p = quantity(p, "cure-age", TIME, 432000.0);
    for (axis, value) in [
        ("source-joint-2024-t3-bare-paa", 1.0),
        ("source-bondline-thickness-known", 0.0),
        ("source-overlap-dimensions-known", 0.0),
        ("source-test-rate-known", 0.0),
    ] {
        p = flag(p, axis, value);
    }
    p
}
fn araldite(cure_temperature: f64, age: f64) -> QueryPoint {
    let p = flags(&[
        ("source-grade-araldite-2011", 1.0),
        ("source-mix-ratio-known", 0.0),
        ("source-mix-ratio-unit-consistent", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-test-direction-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-cure-history-known", 1.0),
        ("source-test-temperature-known", 1.0),
    ]);
    quantity(
        temperature(
            temperature(p, "temperature", 296.15),
            "cure-temperature",
            cure_temperature,
        ),
        "cure-age",
        TIME,
        age,
    )
}
fn araldite_flexure() -> QueryPoint {
    flag(
        flag(
            araldite(313.15, 57600.0),
            "source-flexural-geometry-known",
            0.0,
        ),
        "source-test-rate-known",
        0.0,
    )
}
fn compile(slug: &str, expected_claims: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "adhesive-resin-{}-{nonce}-{slug}.fsmatpk",
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
    assert_eq!(pack.claims().claims_ordered().count(), expected_claims);
    pack
}
fn card(slug: &str, pack: &NormalizedPack) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId {
        chemistry: slug.into(), phase: "cured-solid".into(),
        process: "Original source population; sensible caloric constraint explicitly authored by fixture".into(), revision: 0,
    }, pack.claims().clone(), vec![]).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn sample(pack: &NormalizedPack, property: &str, p: &QueryPoint, expected: f64) -> f64 {
    let a = pack
        .claims()
        .query_typed(&key(property), p, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(a.evidence.value.value, expected, 1e-11, 1e-13);
    assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
    pack.claims().verify_receipt(&a.receipt).unwrap();
    a.evidence.value.value
}
fn refuses(pack: &NormalizedPack, property: &str, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(property), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} {property} accepted {p:?}",
        pack.pack_id()
    );
}
fn original(e: &Epoxy, t: f64) -> f64 {
    let x = t.log10();
    10.0_f64.powf(e.coefficients.iter().rev().fold(0.0, |y, a| y * x + a))
}
fn integral(e: &Epoxy, lo: f64, hi: f64, entropy: bool) -> f64 {
    let value = |t| original(e, t) / if entropy { t } else { 1.0 };
    let n = 2048;
    let step = (hi - lo) / f64::from(n);
    let mut sum = value(lo) + value(hi);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(i) * step);
    }
    sum * step / 3.0
}

#[test]
fn g1_g3_original_cured_epoxy_fits_and_producer_properties_retain_units_and_receipts() {
    let mut values = 20; // EA9394 and Araldite scalar claims.
    for e in EPOXIES {
        let pack = compile(e.slug, 2);
        let line = pack
            .claims()
            .query_line_typed(
                &key(CP),
                "temperature",
                &epoxy_point(e, 2.0),
                &epoxy_point(e, 20.0),
                SINGLE,
            )
            .unwrap();
        assert_eq!(line.knots.len(), e.knots);
        values += e.knots + 1;
        for (t, a) in &line.knots {
            close(a.evidence.value.value, original(e, *t), 1e-11, 1e-13);
            assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
            pack.claims().verify_receipt(&a.receipt).unwrap();
        }
        for pair in line.knots.windows(2) {
            for f in [0.25, 0.5, 0.75] {
                let t = pair[0].0 + f * (pair[1].0 - pair[0].0);
                let a = pack
                    .claims()
                    .query_typed(
                        &key(CP),
                        &epoxy_point(e, t),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap();
                close(a.evidence.value.value, original(e, t), 8.1e-5, 1e-13);
                pack.claims().verify_receipt(&a.receipt).unwrap();
            }
        }
        let cp5 = sample(&pack, CP, &epoxy_point(e, 5.0), original(e, 5.0));
        let ordinate = cp5 / 1000.0 / 5.0_f64.powi(3);
        // Broad visual ordinate intervals from Fig. 1, not additional measured
        // points. They independently catch a 1000-fold unit misinterpretation.
        assert!((e.figure_ordinate_at_5k[0]..=e.figure_ordinate_at_5k[1]).contains(&ordinate));
        let density_context = flag(
            flag(epoxy_base(e), "source-density-temperature-known", 0.0),
            "source-density-cured-state-known",
            0.0,
        );
        sample(
            &pack,
            "nominal-source-density",
            &density_context,
            e.nominal_density,
        );
    }
    assert_eq!(values, 1022);
    let ea = compile(EA, 17);
    sample(&ea, NOMINAL_K, &ea_nominal_k(), 0.331);
    for (property, value) in [
        ("tensile-modulus", 4237e6),
        ("tensile-strength", 46e6),
        ("tensile-elongation-at-break", 0.0166),
    ] {
        sample(&ea, property, &ea_tension(), value);
    }
    sample(&ea, "shear-modulus", &ea_shear(1.0), 1461e6);
    sample(&ea, "shear-modulus", &ea_shear(2.0), 1027e6);
    sample(&ea, "compressive-strength", &ea_compression(), 68.9e6);
    for (t, alpha) in [(313.15, 55.6e-6), (373.15, 80.6e-6)] {
        sample(
            &ea,
            "linear-thermal-expansion-coefficient",
            &ea_expansion(t),
            alpha,
        );
    }
    for (kelvin, mpa) in [
        (218.15, 22.7),
        (298.15, 28.9),
        (355.15, 20.7),
        (366.15, 20.0),
        (394.15, 15.8),
        (422.15, 11.0),
        (450.15, 8.3),
        (477.15, 4.1),
    ] {
        sample(&ea, "lap-shear-strength", &ea_joint(kelvin), mpa * 1e6);
    }
    let ar = compile(ARALDITE, 3);
    sample(&ar, K, &araldite(373.15, 1200.0), 0.22);
    sample(&ar, "flexural-strength", &araldite_flexure(), 60.4e6);
    sample(&ar, "flexural-modulus", &araldite_flexure(), 1904.1e6);
}

#[test]
fn g0_cure_geometry_temperature_and_missing_thermal_laws_cannot_be_silently_joined() {
    for e in EPOXIES {
        let pack = compile(e.slug, 2);
        let p = epoxy_point(e, 10.0);
        for axis in [
            "engineering-reference-model",
            "source-measurement-vacuum",
            "source-pressure-known",
            "source-caloric-constraint-known",
            "source-fit-unit-labelled",
            "source-fit-unit-interpretation-j-per-kg-k",
            "source-cure-history-known",
            "source-moisture-condition-known",
            "source-specimens-matched",
            &format!("source-grade-{}", e.grade),
        ] {
            refuses(
                &pack,
                CP,
                &p.clone().with(axis, 1.0 - p.axes()[axis]).unwrap(),
            );
        }
        refuses(
            &pack,
            CP,
            &p.clone()
                .with("resin-to-hardener-mass-ratio", e.ratio + 0.01)
                .unwrap(),
        );
        for t in [2.0_f64.next_down(), 20.0_f64.next_up()] {
            refuses(&pack, CP, &epoxy_point(e, t));
        }
        refuses(
            &pack,
            CP,
            &p.with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                10.0,
            )
            .unwrap(),
        );
        refuses(
            &pack,
            CP,
            &temperature(QueryPoint::new(), "temperature", 10.0),
        );
        for missing in ["density", "volumetric-heat-capacity", K, "young-modulus"] {
            assert!(pack.claims().claims_for(missing).is_empty());
        }
    }
    let ea = compile(EA, 17);
    for missing in [
        K,
        CP,
        "density",
        "young-modulus",
        "poisson-ratio",
        "shear-strength",
        "linear-expansion-relative-to-293k",
    ] {
        assert!(ea.claims().claims_for(missing).is_empty());
    }
    refuses(&ea, K, &ea_nominal_k());
    refuses(
        &ea,
        "tensile-modulus",
        &ea_tension().with("casting-thickness", 0.002).unwrap(),
    );
    refuses(
        &ea,
        "tensile-modulus",
        &ea_tension().with("source-cure-history-known", 1.0).unwrap(),
    );
    refuses(
        &ea,
        "tensile-modulus",
        &ea_tension().with("temperature", 299.15).unwrap(),
    );
    refuses(&ea, "shear-modulus", &ea_shear(0.0));
    refuses(
        &ea,
        "compressive-strength",
        &ea_compression().with("specimen-height", 0.01).unwrap(),
    );
    refuses(
        &ea,
        "lap-shear-strength",
        &ea_joint(298.15).with("cure-age", 604800.0).unwrap(),
    );
    refuses(
        &ea,
        "lap-shear-strength",
        &ea_joint(298.15)
            .with("source-joint-2024-t3-bare-paa", 0.0)
            .unwrap(),
    );
    refuses(
        &ea,
        "lap-shear-strength",
        &ea_joint(298.15)
            .with("source-bondline-thickness-known", 1.0)
            .unwrap(),
    );
    refuses(&ea, "lap-shear-strength", &ea_joint(308.15));
    assert!(
        ea.claims()
            .query_line_typed(
                &key("lap-shear-strength"),
                "temperature",
                &ea_joint(298.15),
                &ea_joint(355.15),
                SINGLE
            )
            .is_err()
    );
    assert!(
        ea.claims()
            .query_line_typed(
                &key("linear-thermal-expansion-coefficient"),
                "temperature",
                &ea_expansion(313.15),
                &ea_expansion(373.15),
                SINGLE
            )
            .is_err()
    );
    let ar = compile(ARALDITE, 3);
    for missing in [
        CP,
        "density",
        "young-modulus",
        "linear-thermal-expansion-coefficient",
    ] {
        assert!(ar.claims().claims_for(missing).is_empty());
    }
    let k_context = araldite(373.15, 1200.0);
    for t in [296.15_f64.next_down(), 296.15_f64.next_up()] {
        refuses(&ar, K, &k_context.clone().with("temperature", t).unwrap());
    }
    refuses(&ar, K, &araldite(313.15, 57600.0));
    refuses(
        &ar,
        K,
        &k_context
            .clone()
            .with("source-mix-ratio-unit-consistent", 1.0)
            .unwrap(),
    );
    let wrong_flexure = araldite_flexure()
        .with("cure-temperature", 373.15)
        .unwrap()
        .with("cure-age", 1200.0)
        .unwrap();
    refuses(&ar, "flexural-modulus", &wrong_flexure);
    assert!(
        ar.claims()
            .query_line_typed(
                &key(K),
                "temperature",
                &k_context.clone().with("temperature", 295.15).unwrap(),
                &k_context.with("temperature", 297.15).unwrap(),
                SINGLE
            )
            .is_err()
    );
}

#[test]
fn g1_g3_real_epoxy_heating_and_explicit_adhesive_conduction_approximations() {
    let packs: Vec<_> = EPOXIES.iter().map(|e| compile(e.slug, 2)).collect();
    let cards: Vec<_> = EPOXIES
        .iter()
        .zip(&packs)
        .map(|(e, p)| card(e.slug, p))
        .collect();
    let bounds: Vec<_> = EPOXIES
        .iter()
        .map(|e| [epoxy_point(e, 2.0), epoxy_point(e, 20.0)])
        .collect();
    let cp = key(CP);
    let query = |i: usize, mass: f64| HeatStorageComponentQuery {
        id: (i + 1) as u64,
        card: &cards[i],
        capacity: &cp,
        amount: HeatStorageAmount::MassKg(mass),
        temperature_axis: "temperature",
        lower: &bounds[i][0],
        upper: &bounds[i][1],
        selection: SINGLE,
    };
    let mut energies = Vec::new();
    for (i, e) in EPOXIES.iter().enumerate() {
        // The fixture explicitly chooses a sensible caloric interpretation and
        // supplied mass. Nominal source-density metadata is never consumed.
        let storage = PreparedHeatStorage::try_new(&[query(i, 0.01)]).unwrap();
        let energy = storage.energy_change_j(3.0, 18.0).unwrap();
        close(energy, 0.01 * integral(e, 3.0, 18.0, false), 8.1e-5, 1e-10);
        close(
            storage.sensible_entropy_change_j_k(3.0, 18.0).unwrap(),
            0.01 * integral(e, 3.0, 18.0, true),
            8.1e-5,
            1e-11,
        );
        let delivered = storage.energy_change_j(3.0, 13.5).unwrap();
        close(
            storage
                .temperature_after_energy_change(3.0, delivered)
                .unwrap(),
            13.5,
            1e-11,
            1e-10,
        );
        assert!(
            storage
                .temperature_after_energy_change(
                    3.0,
                    1.01 * storage.energy_change_j(3.0, 20.0).unwrap()
                )
                .is_err()
        );
        let line = packs[i]
            .claims()
            .query_line_typed(
                &cp,
                "temperature",
                &epoxy_point(e, 3.0),
                &epoxy_point(e, 18.0),
                SINGLE,
            )
            .unwrap();
        let trapezoids: f64 = line
            .knots
            .windows(2)
            .map(|p| {
                0.005
                    * (p[1].0 - p[0].0)
                    * (p[0].1.evidence.value.value + p[1].1.evidence.value.value)
            })
            .sum();
        close(energy, trapezoids, 1e-12, 1e-11);
        energies.push(energy);
    }
    assert!(energies[1] > 2.0 * energies[0]);
    assert!(energies[2] > energies[3]);
    // Co-temperature assembly of distinct cured solids, not a chemical mixing
    // rule: the paper explicitly observes nonlinear filler-dependent capacities.
    let assembly = PreparedHeatStorage::try_new(&[query(0, 0.01), query(2, 0.02)]).unwrap();
    close(
        assembly.energy_change_j(3.0, 18.0).unwrap(),
        energies[0] + 2.0 * energies[2],
        1e-12,
        1e-10,
    );

    let ea = compile(EA, 17);
    let ar = compile(ARALDITE, 3);
    let k_ea = sample(&ea, NOMINAL_K, &ea_nominal_k(), 0.331);
    let k_ar = sample(&ar, K, &araldite(373.15, 1200.0), 0.22);
    // An explicitly authored two-kelvin constant-k approximation, deliberately
    // outside source claims: EA9394 has unspecified test temperature and Araldite
    // supplies only one point. These charts are not written back to either card
    // and do not establish qualified temperature-dependent producer behavior.
    let curves =
        [k_ea, k_ar].map(|k| ConductivityCurve::new(vec![[295.15, k], [297.15, k]]).unwrap());
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    };
    let lengths = [0.0002, 0.0003];
    let contact = 0.0001;
    let area = 0.01;
    let model = TwoLayerConduction::new([&curves[0], &curves[1]], lengths, area, contact).unwrap();
    let response = model.evaluate([297.15, 295.15], config, || Ok(())).unwrap();
    let expected = 2.0 * area / (lengths[0] / k_ea + contact + lengths[1] / k_ar);
    close(response.power_a_to_b_w, expected, 1e-11, 1e-9);
    close(
        response.interface_temperature_k[0] - response.interface_temperature_k[1],
        expected * contact / area,
        1e-10,
        1e-9,
    );
    let reverse = model.evaluate([295.15, 297.15], config, || Ok(())).unwrap();
    close(reverse.power_a_to_b_w, -expected, 1e-11, 1e-9);
    let substitute =
        TwoLayerConduction::new([&curves[1], &curves[1]], lengths, area, contact).unwrap();
    assert!(
        substitute
            .evaluate([297.15, 295.15], config, || Ok(()))
            .unwrap()
            .power_a_to_b_w
            < response.power_a_to_b_w
    );
}
