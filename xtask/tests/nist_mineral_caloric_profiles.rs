//! G1/G3: first solid-phase NIST Shomate branches as qualified mass-based Cp.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CP: &str = "specific-heat-capacity";

struct Mineral {
    slug: &'static str,
    phase: &'static str,
    molecular_weight: f64,
    upper_k: f64,
    knots: usize,
    coefficients: [f64; 5],
    cp_300_molar: f64,
    table_rounding: f64,
}

// Independently transcribed first solid columns and molecular weights from
// the primary WebBook pages cited in each pack. The 300 K values are rounded
// source-generated table entries, not new measured observations. MgO's source
// says solid; it does not qualify a particular periclase producer specimen.
const MINERALS: &[Mineral] = &[
    Mineral {
        slug: "quartz-alpha-nist-shomate-caloric",
        phase: "quartz-alpha",
        molecular_weight: 60.0843,
        upper_k: 847.0,
        knots: 179,
        coefficients: [-6.076591, 251.6755, -324.7964, 168.5604, 0.002548],
        cp_300_molar: 44.77,
        table_rounding: 0.005,
    },
    Mineral {
        slug: "alumina-alpha-nist-shomate-caloric",
        phase: "alumina-alpha",
        molecular_weight: 101.9613,
        upper_k: 2327.0,
        knots: 138,
        coefficients: [102.4290, 38.74980, -15.91090, 2.628181, -3.007551],
        cp_300_molar: 79.28,
        table_rounding: 0.005,
    },
    Mineral {
        slug: "magnesium-oxide-solid-nist-shomate-caloric",
        phase: "magnesium-oxide-solid",
        molecular_weight: 40.3044,
        upper_k: 3105.0,
        knots: 115,
        coefficients: [47.25995, 5.681621, -0.872665, 0.104300, -1.053955],
        cp_300_molar: 37.18,
        table_rounding: 0.005,
    },
    Mineral {
        slug: "rutile-nist-shomate-caloric",
        phase: "rutile",
        molecular_weight: 79.866,
        upper_k: 2000.0,
        knots: 137,
        coefficients: [67.29830, 18.70940, -11.57900, 2.449561, -1.485471],
        cp_300_molar: 55.43,
        table_rounding: 0.005,
    },
    Mineral {
        slug: "hematite-nist-shomate-caloric",
        phase: "hematite",
        molecular_weight: 159.688,
        upper_k: 950.0,
        knots: 74,
        coefficients: [93.43834, 108.3577, -50.86447, 25.58683, -1.611330],
        cp_300_molar: 104.2,
        table_rounding: 0.05,
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key() -> PropertyKey {
    PropertyKey::new(CP, Dims([2, 0, -2, -1, 0, 0]))
}
fn point(m: &Mineral, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
        .with_quantity(
            "source-standard-pressure",
            semantic(QuantityKind::Pressure),
            100000.0,
        )
        .unwrap();
    for (axis, value) in [
        ("source-standard-state", 1.0),
        ("source-fit-branch", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (&format!("source-phase-{}", m.phase), 1.0),
    ] {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p
}
fn compile(m: &Mineral) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-mineral-caloric-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        m.slug
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", m.slug)))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        m.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), 1);
    assert_eq!(pack.claims().observation_ids().count(), 4);
    pack
}
fn card(m: &Mineral, pack: &NormalizedPack) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId {
        chemistry: m.slug.into(), phase: m.phase.into(),
        process: "Selected JANAF solid reference at one bar; fixed mass, no phase crossing or specimen qualification".into(),
        revision: 0,
    }, pack.claims().clone(), vec![]).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn original_cp(m: &Mineral, temperature: f64) -> f64 {
    let [a, b, c, d, e] = m.coefficients;
    let t = temperature / 1000.0;
    (a + b * t + c * t * t + d * t * t * t + e / (t * t)) / (0.001 * m.molecular_weight)
}
// Exact analytic Shomate primitives, independently of the stored linear knots.
// F-H and G cancel in differences. This does not consume formation enthalpy,
// impose H(298.15)=0 on a rounded fit, or turn a difference into absolute entropy.
fn enthalpy_primitive(m: &Mineral, temperature: f64) -> f64 {
    let [a, b, c, d, e] = m.coefficients;
    let t = temperature / 1000.0;
    1000.0 * (a * t + b * t * t / 2.0 + c * t.powi(3) / 3.0 + d * t.powi(4) / 4.0 - e / t)
        / (0.001 * m.molecular_weight)
}
fn entropy_primitive(m: &Mineral, temperature: f64) -> f64 {
    let [a, b, c, d, e] = m.coefficients;
    let t = temperature / 1000.0;
    (a * t.ln() + b * t + c * t * t / 2.0 + d * t.powi(3) / 3.0 - e / (2.0 * t * t))
        / (0.001 * m.molecular_weight)
}
fn curvature_relative_bound(m: &Mineral, lo: f64, hi: f64) -> f64 {
    let [a, b, c, d, e] = m.coefficients;
    let (tl, th) = (lo / 1000.0, hi / 1000.0);
    let lower = a
        + (b * tl).min(b * th)
        + (c * tl * tl).min(c * th * th)
        + (d * tl.powi(3)).min(d * th.powi(3))
        + (e / tl.powi(2)).min(e / th.powi(2));
    assert!(lower > 0.0);
    let curvature =
        2.0 * c.abs() / 1e6 + 6.0 * d.abs() * hi / 1e9 + 6.0 * e.abs() * 1e6 / lo.powi(4);
    (hi - lo).powi(2) * curvature / (8.0 * lower)
}

#[test]
fn g1_original_molar_equations_mass_units_and_source_receipts_compile() {
    let cp = key();
    let mut values = 0;
    for m in MINERALS {
        let pack = compile(m);
        let line = pack
            .claims()
            .query_line_typed(
                &cp,
                "temperature",
                &point(m, 298.0),
                &point(m, m.upper_k),
                SINGLE,
            )
            .unwrap();
        assert_eq!(line.knots.len(), m.knots);
        values += line.knots.len();
        for (t, a) in &line.knots {
            close(a.evidence.value.value, original_cp(m, *t), 1e-12, 1e-10);
            assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(a.receipt.source_hashes.len(), 5); // Claim plus all four source observations.
            pack.claims().verify_receipt(&a.receipt).unwrap();
        }
        let a = pack
            .claims()
            .query_typed(&cp, &point(m, 300.0), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        close(
            a.evidence.value.value * 0.001 * m.molecular_weight,
            m.cp_300_molar,
            0.0,
            m.table_rounding,
        );
        for w in line.knots.windows(2) {
            let bound = curvature_relative_bound(m, w[0].0, w[1].0);
            assert!(bound <= 5.0001e-5);
            for f in [0.25, 0.5, 0.75] {
                let t = w[0].0 + f * (w[1].0 - w[0].0);
                let a = pack
                    .claims()
                    .query_typed(&cp, &point(m, t), SelectionPolicy::SingleClaimOnly)
                    .unwrap();
                close(
                    a.evidence.value.value,
                    original_cp(m, t),
                    bound + 1e-12,
                    1e-10,
                );
                pack.claims().verify_receipt(&a.receipt).unwrap();
            }
        }
    }
    assert_eq!(values, 643);
}

#[test]
fn g1_phase_branch_standard_state_and_temperature_boundaries_refuse() {
    for m in MINERALS {
        let pack = compile(m);
        let cp = key();
        let valid = point(m, 300.0);
        for (axis, value) in [
            ("temperature", 298.0_f64.next_down()),
            ("temperature", m.upper_k.next_up()),
            ("source-standard-pressure", 101325.0),
            ("source-standard-state", 0.0),
            ("source-fit-branch", 2.0),
            ("engineering-reference-model", 0.0),
            ("source-specimen-qualified", 1.0),
            (&format!("source-phase-{}", m.phase), 0.0),
        ] {
            let unsupported = valid.clone().with(axis, value).unwrap();
            assert!(
                pack.claims()
                    .query_typed(&cp, &unsupported, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "{} accepted {axis}={value}",
                m.slug
            );
        }
        assert!(
            pack.claims()
                .query_typed(&cp, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        let wrong_units = PropertyKey::new(CP, Dims([2, 0, -2, 0, 0, 0]));
        assert!(
            pack.claims()
                .query_typed(&wrong_units, &valid, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        for absent in [
            "density",
            "thermal-conductivity",
            "latent-heat",
            "specific-heat-capacity-cv",
        ] {
            assert!(pack.claims().claims_for(absent).is_empty());
        }
        let c = card(m, &pack);
        let lower = point(m, 298.0);
        let upper = point(m, m.upper_k.next_up());
        assert!(
            PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
                id: 1,
                card: &c,
                capacity: &cp,
                amount: HeatStorageAmount::MassKg(0.1),
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SINGLE,
            }])
            .is_err()
        );
        // Inclusive upper endpoint is the selected first-branch limit, never
        // permission to select the next polymorph/fit at that same temperature.
        let endpoint = pack
            .claims()
            .query_typed(&cp, &point(m, m.upper_k), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        close(
            endpoint.evidence.value.value,
            original_cp(m, m.upper_k),
            1e-12,
            1e-10,
        );
        assert!(
            pack.claims()
                .query_typed(
                    &cp,
                    &point(m, m.upper_k).with("source-fit-branch", 2.0).unwrap(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn g3_real_mineral_heating_energy_inverse_entropy_and_explicit_composites() {
    let packs: Vec<_> = MINERALS.iter().map(compile).collect();
    let cards: Vec<_> = MINERALS
        .iter()
        .zip(&packs)
        .map(|(m, p)| card(m, p))
        .collect();
    let bounds: Vec<_> = MINERALS
        .iter()
        .map(|m| [point(m, 298.0), point(m, m.upper_k)])
        .collect();
    let cp = key();
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
    let mut entropies = Vec::new();
    let mut heated = Vec::new();
    for (i, m) in MINERALS.iter().enumerate() {
        // Declared fixture: fixed 0.1 kg, one bar standard-state Cp, uniform T,
        // no loss or phase crossing, and no non-pV work. Heat equals delta H.
        let storage = PreparedHeatStorage::try_new(&[query(i, 0.1)]).unwrap();
        let expected = 0.1 * (enthalpy_primitive(m, 700.0) - enthalpy_primitive(m, 300.0));
        let energy = storage.energy_change_j(300.0, 700.0).unwrap();
        let entropy = storage.sensible_entropy_change_j_k(300.0, 700.0).unwrap();
        close(energy, expected, 5.1e-5, 1e-9);
        close(
            entropy,
            0.1 * (entropy_primitive(m, 700.0) - entropy_primitive(m, 300.0)),
            5.1e-5,
            1e-10,
        );
        close(
            storage
                .temperature_after_energy_change(300.0, expected)
                .unwrap(),
            700.0,
            0.0,
            0.03,
        );
        close(
            storage
                .temperature_after_energy_change(300.0, energy)
                .unwrap(),
            700.0,
            1e-12,
            1e-9,
        );
        close(
            storage.energy_change_j(700.0, 300.0).unwrap(),
            -energy,
            1e-12,
            1e-9,
        );
        close(
            storage.sensible_entropy_change_j_k(700.0, 300.0).unwrap(),
            -entropy,
            1e-12,
            1e-10,
        );
        let line = packs[i]
            .claims()
            .query_line_typed(
                &cp,
                "temperature",
                &point(m, 300.0),
                &point(m, 700.0),
                SINGLE,
            )
            .unwrap();
        let trapezoids: f64 = line
            .knots
            .windows(2)
            .map(|w| {
                0.05 * (w[1].0 - w[0].0)
                    * (w[0].1.evidence.value.value + w[1].1.evidence.value.value)
            })
            .sum();
        close(energy, trapezoids, 1e-12, 1e-9);
        let temperature = storage
            .temperature_after_energy_change(300.0, 100.0 * 100.0)
            .unwrap();
        let frozen_cp_temperature = 300.0 + 10000.0 / (0.1 * original_cp(m, 300.0));
        assert!(frozen_cp_temperature - temperature > 1.0);
        let doubled_mass = PreparedHeatStorage::try_new(&[query(i, 0.2)]).unwrap();
        close(
            doubled_mass
                .temperature_after_energy_change(300.0, 10000.0)
                .unwrap(),
            storage
                .temperature_after_energy_change(300.0, 5000.0)
                .unwrap(),
            1e-12,
            1e-9,
        );
        assert!(
            storage
                .temperature_after_energy_change(
                    300.0,
                    1.01 * storage.energy_change_j(300.0, m.upper_k).unwrap()
                )
                .is_err()
        );
        println!(
            "{}: 0.1 kg, 300->700 K: {energy:.6} J, {entropy:.9} J/K; 10 kJ from 300 K -> {temperature:.9} K",
            m.slug
        );
        energies.push(energy);
        entropies.push(entropy);
        heated.push(temperature);
    }
    assert!(heated[4] > heated[2] + 30.0); // Source substitution changes the actual thermal response.
    // An explicitly co-temperature assembly of quartz, alpha alumina and
    // rutile at fixed masses. No reacted mixture, effective k, or compatible
    // sintered specimen is inferred from these distinct reference curves.
    let assembly =
        PreparedHeatStorage::try_new(&[query(0, 0.1), query(1, 0.2), query(3, 0.3)]).unwrap();
    assert_eq!(assembly.components().len(), 3);
    close(
        assembly.energy_change_j(300.0, 700.0).unwrap(),
        energies[0] + 2.0 * energies[1] + 3.0 * energies[3],
        1e-12,
        1e-8,
    );
    close(
        assembly.sensible_entropy_change_j_k(300.0, 700.0).unwrap(),
        entropies[0] + 2.0 * entropies[1] + 3.0 * entropies[3],
        1e-12,
        1e-9,
    );
    let heat = assembly.energy_change_j(300.0, 600.0).unwrap();
    close(
        assembly
            .temperature_after_energy_change(300.0, heat)
            .unwrap(),
        600.0,
        1e-12,
        1e-9,
    );
    assert!(
        assembly
            .energy_change_j(300.0, 847.0_f64.next_up())
            .is_err()
    );
}
