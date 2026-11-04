//! The derived loss rows of the Lampis Dominant seed pack (MagentaBadger 144, 2026-11-03) reach the engine as zeta = 1 / (2 Q) on every string mode.
//! The numbers are the seed pack's own, read from its `properties.tsv` at compile time; the two cards differ only in their loss selection,
//! so the card without the rows is the control.
use super::*;
use fs_couple::material_assembly::measured_string::MeasuredStringProperties;

const SEED: &str =
    include_str!("../../../../../data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/properties.tsv");
const FRICTION: &str = "loss_factor_tension_dominant";
const BENDING: &str = "loss_factor_bending_dominant";
const MEAN_Q: &str = "mean_quality_factor_first_30_modes_dominant";
const PREPARED: &str = "prepared-linear 64 67108864 200000000000";
/// 64 bending modes on two polarizations plus the plate need a larger network allowance than the 4-mode fixture.
const PREPARED_64: &str = "prepared-linear 512 4294967296 1000000000000000";
const PROPERTY_NAMES: [&str; 5] = [
    "playing_length",
    "diameter_dominant",
    "tension_dominant",
    "linear_density_dominant",
    "bending_stiffness_dominant",
];

fn seed(name: &str) -> f64 {
    SEED.lines()
        .find_map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            (c.len() > 4 && c[0] == "scalar" && c[3] == name).then(|| c[4].parse::<f64>().unwrap())
        })
        .unwrap_or_else(|| panic!("the seed pack has no scalar {name}"))
}

/// The Dominant's five section coefficients and its two loss rows, exactly as the seed pack prints them.
fn seed_pack() -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let artifact = fs_blake3::ContentHash([92; 32]);
    let provenance = Provenance {
        source: "numbers read from data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/properties.tsv".into(),
        license: "CC-BY-4.0".into(),
        artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "the seed pack's Dominant coefficients, re-entered as a manufactured fixture".into(),
            method: "explicit scalar section coefficients and the derived constant-loss rows".into(),
            artifact,
            caveats: "not a measurement; the compiler path of the seed pack is xtask matdb-pack".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let properties = MeasuredStringProperties {
        length: PROPERTY_NAMES[0].into(),
        diameter: PROPERTY_NAMES[1].into(),
        tension: PROPERTY_NAMES[2].into(),
        linear_density: PROPERTY_NAMES[3].into(),
        bending_stiffness: PROPERTY_NAMES[4].into(),
    };
    let mut rows: Vec<(&str, Dims)> = Vec::new();
    for (requirement, name) in properties.requirements().unwrap().into_iter().zip(PROPERTY_NAMES) {
        rows.push((name, requirement.dims()));
    }
    rows.push((FRICTION, Dims::NONE));
    rows.push((BENDING, Dims::NONE));
    for (name, dims) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value: seed(name), dims },
                validity: ValidityDomain::unconstrained().with("temperature_known", 0.0, 0.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let pack = NormalizedPack::new(
        "seed-pack-numbers-fixture",
        "synthetic-test-v1",
        artifact,
        "manufactured test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

/// A fixture with 64 bending modes (the posted cello's count) and a short run.
fn card(fixture: &str) -> String {
    let card = fixture
        .replace("48000 4096 0.01 128", "48000 128 0.01 128")
        .replace("116.3 4 0 0.8", "116.3 64 0 0.8")
        .replace(PREPARED, PREPARED_64);
    assert!(card.contains(PREPARED_64) && card.contains("116.3 64 0 0.8"), "the fixture's rows moved");
    card
}

#[test]
fn the_dominant_loss_rows_give_every_string_mode_q_1064_through_the_engines_own_law() {
    let plain_text = include_str!("../measured-cello-g.material-strings");
    let loss_text = include_str!("../measured-cello-g-loss.material-strings");
    // The two fixtures differ only in the loss selection.
    assert_eq!(
        loss_text.replace(&format!("loss-factors {FRICTION} {BENDING}"), "no-intrinsic-loss"),
        plain_text
    );
    // The rows are the constant-loss reading of the printed mean Q.
    let q = seed(MEAN_Q);
    assert_eq!(q, 1064.0);
    let eta = 1.0 / q;
    for row in [FRICTION, BENDING] {
        assert!((seed(row) / eta - 1.0).abs() < 1e-15, "{row} is not 1/Q");
    }
    scoped(|cx, gate| {
        let compile = |fixture: &str| {
            input::Spec::parse(card(fixture).as_bytes())
                .unwrap()
                .compile(cx, gate, |name| {
                    Ok(if name == "cello-g-whole.fsmatpk" { seed_pack() } else { pack(name, 7800.0, 200_000.0, None) })
                })
                .unwrap()
        };
        let plain = compile(plain_text);
        let bound = compile(loss_text);
        let (a, b) = (plain.performance.bank(), bound.performance.bank());
        assert_eq!(a.cards().len(), b.cards().len());
        let mut checked = 0;
        for (a, b) in a.cards().iter().zip(b.cards()) {
            assert!(a.zetas.len() >= 64 && a.zetas.len() == b.zetas.len());
            // Polarizations repeat the same 64 modes.
            assert_eq!(a.zetas.len() % 64, 0);
            for (index, (a, b)) in a.zetas.iter().zip(&b.zetas).enumerate() {
                // The computed air drag is in both; the loss rows add exactly 1 / (2 Q) to each mode: Q_k = Q for k = 1..64.
                let added = b - a;
                assert!(
                    (added / (0.5 * eta) - 1.0).abs() < 1e-9,
                    "mode index {index}: the rows added {added:e}, 1/(2 Q) is {:e}",
                    0.5 * eta
                );
                checked += 1;
            }
        }
        assert!(checked >= 64, "checked {checked} modes");
    });
}
