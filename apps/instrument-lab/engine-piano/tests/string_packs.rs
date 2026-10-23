//! Real packs through the committed constructions (bead pexfp.369, gap 1). The compiled seed packs in
//! tests/fixtures/packs are byte-for-byte the offline compiler's output for data/matdb/seed-v1
//! (xtask/tests/piano_string_packs.rs guards that), applied by the committed 88-key constructions to
//! the embedded Model D scale. Every expected number is computed here from the packs' own stored claims
//! and the construction rows; nothing is tuned to agree.
use fs_live_piano::{Config, Description, PianoAudio, kernel::SCALE_CSV,
    preparation::{self, materials::{StringMaterials, operating_point}},
    refmods::{bank_string_losses, geometry::{self, Course}, linear::Bank}};
use fs_matdb::{NormalizedMaterialCardPack, PropertyClaim, PropertyValue};
use fs_material::{gas::{GasSpec, GasState}, visco::circular_ted::circular_relaxation_shape};
use std::{collections::{BTreeMap, BTreeSet}, f64::consts::{PI, TAU}, sync::Arc};

const GOLD: &[u8] = include_bytes!("fixtures/packs/gold-pure-nbs-c447.fsmatpk");
const STEEL: &[u8] = include_bytes!("fixtures/packs/steel-music-wire-string-heptonstall-ligo-p1400229.fsmatpk");
const QQ: &[u8] = include_bytes!("fixtures/packs/steel-music-wire-min-strength-qq-w-470b.fsmatpk");
const COPPER: &[u8] = include_bytes!("fixtures/packs/copper-pure-nbs-c447-string-wrap.fsmatpk");
const SILVER: &[u8] = include_bytes!("fixtures/packs/silver-pure-nbs-c447.fsmatpk");
const SILVER_UTS: &[u8] = include_bytes!("fixtures/packs/silver-drawn-hard-uts-nbs-c101.fsmatpk");
const GOLD_CONCERT: &str = include_str!("../data/constructions/gold-concert-pitch-rt0425.csv");
const GOLD_HALF: &str = include_str!("../data/constructions/gold-half-strength-rt0425.csv");
const STEEL_WOUND: &str = include_str!("../data/constructions/steel-copper-wound-rt0425.csv");
const SILVER_CONCERT: &str = include_str!("../data/constructions/silver-concert-pitch-rt0425.csv");
/// The constructions' measured-excess property (the steel pack's in-band loss tangent).
const LOSS: &str = "bending_loss_tangent";
const T0: f64 = 293.15;
/// The wound bass of every 88-key construction here is back-solved (gap 3): no cited per-key geometry.
const GEOMETRY: &str = "\"geometryCaveat\":\"core and wrap diameters of keys 21-40 are back-solved, not measured or \
    cited; every tension, stress, limit, pitch and loss computed from them inherits that\"";

fn materials(packs: &[&[u8]], construction: &str) -> Result<StringMaterials, String> {
    StringMaterials::from_packs(packs, construction, &operating_point("T", T0, "frequency", 1000.0)?, LOSS)
}
fn scale() -> Vec<Course> { geometry::read_scale(SCALE_CSV).unwrap() }
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs(), "{a:e} != {b:e} (relative {tolerance:e})");
}

/// One pack's single claim for `property`.
fn claim(pack: &NormalizedMaterialCardPack, property: &str) -> PropertyClaim {
    let claims = pack.card().claims().claims_for(property);
    assert_eq!(claims.len(), 1, "{property}");
    claims[0].1.clone()
}
/// The stored SI value of a scalar claim.
fn scalar(pack: &[u8], property: &str) -> f64 {
    match claim(&NormalizedMaterialCardPack::from_bytes(pack).unwrap(), property).value {
        PropertyValue::Scalar { value, .. } => value,
        other => panic!("{property}: not a scalar: {other:?}"),
    }
}

/// Per key, the outermost core radius [m] and whether the key is wound.
fn rows(construction: &str) -> (BTreeMap<u8, f64>, BTreeSet<u8>) {
    let (mut cores, mut wound) = (BTreeMap::new(), BTreeSet::new());
    for line in construction.lines().map(|line| line.split('#').next().unwrap_or("").trim()) {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 6 { continue; }
        let Ok(key) = fields[0].parse::<u8>() else { continue };
        match fields[1] {
            "core" => { cores.insert(key, fields[2].parse::<f64>().unwrap()); }
            "wrap" => { wound.insert(key); }
            _ => {}
        }
    }
    (cores, wound)
}

/// Numeric field of one flat JSON object; `None` for null.
fn field(object: &str, name: &str) -> Option<f64> {
    let start = object.find(&format!("\"{name}\":")).unwrap_or_else(|| panic!("{name} in {object}")) + name.len() + 3;
    let text = &object[start..][..object[start..].find([',', '}', ']']).unwrap()];
    if text == "null" { None } else { Some(text.parse().unwrap()) }
}
/// Each overstressed key's object in a structured refusal, in key order.
fn overstressed(refusal: &str) -> Vec<String> {
    let list = &refusal[refusal.find("\"overstressed\":[").expect("overstressed list") + 16..];
    list.split("{\"key\":").skip(1).map(|object| format!("{{\"key\":{object}")).collect()
}
/// A section's total line density [kg/m] from the materials receipt.
fn line_density(json: &str, key: u8) -> f64 {
    let section = &json[json.find(&format!("{{\"key\":{key},\"stateIdentities\"")).expect("section")..];
    field(section, "lineDensityKgM").unwrap()
}

/// Check a structured refusal that names every key: each key's limit is `limit_pa`, it is above its
/// limit, and its highest pitch is the pinned stiff-string fundamental at the core's capacity
/// `limit pi r^2`, with EI = E pi r^4 / 4 and the string's line density (rho pi r^2 for a plain core,
/// computed here; the receipt's for a wound one). Returns the lowest and highest f_max.
fn check_refusal(error: &str, materials: &StringMaterials, construction: &str, limit_pa: f64, core: &[u8])
    -> (f64, f64) {
    let refusal = &error[error.find("{\"refusal\"").unwrap_or_else(|| panic!("structured refusal: {error}"))..];
    assert!(refusal.contains("\"keys\":88") && refusal.contains(GEOMETRY), "{refusal}");
    assert_eq!(materials.back_solved_keys(), Some((21, 40)));
    let objects = overstressed(refusal);
    assert_eq!(objects.len(), 88, "{refusal}");
    let (cores, wound) = rows(construction);
    let (rho, young) = (scalar(core, "density"), scalar(core, "young_modulus"));
    let json = materials.json();
    let (mut lowest, mut highest) = (f64::INFINITY, 0.0_f64);
    for (object, course) in objects.iter().zip(scale()) {
        assert_eq!(field(object, "key"), Some(f64::from(course.midi)));
        assert_eq!(field(object, "limitPa"), Some(limit_pa), "{object}");
        assert!(field(object, "utilization").unwrap() > 1.0, "{object}");
        let r = cores[&course.midi];
        let mu = if wound.contains(&course.midi) { line_density(&json, course.midi) } else { rho * PI * r * r };
        if !wound.contains(&course.midi) { near(line_density(&json, course.midi), mu, 1e-12); }
        let k = PI / course.length_m;
        let f_max = ((limit_pa * PI * r * r * k * k + young * PI * r.powi(4) / 4.0 * k.powi(4)) / mu).sqrt() / TAU;
        let printed = field(object, "fMaxHz").unwrap();
        near(printed, f_max, 1e-12);
        lowest = lowest.min(printed);
        highest = highest.max(printed);
    }
    (lowest, highest)
}

// G1 (DONE-WHEN): the real gold pack at concert pitch refuses with one structured refusal that lists
// all 88 keys. Each key carries the 0.05 in. cold-rolled wire's strength, taken at every core diameter
// by the construction's declared extrapolation and stated as such, and its highest attainable pitch.
#[test]
fn g1_real_gold_at_concert_pitch_refuses_all_88_keys_naming_each_highest_pitch() {
    let gold = materials(&[GOLD, GOLD], GOLD_CONCERT).unwrap();
    assert!(gold.extrapolate_diameter() && !gold.diameter_free());
    let error = gold.apply(&scale()).unwrap_err();
    let property = "ultimate_tensile_strength_cold_rolled_wire";
    assert!(error.contains(&format!("\"limitProperty\":\"{property}\"")), "{error}");
    let caveat = format!("strength extrapolated; '{property}' from gold-pure-nbs-c447 holds only for 1.2700 mm wire \
        and is taken at every core diameter, diameter effect unbounded");
    assert!(error.contains(&format!("\"limitCaveats\":[\"{caveat}\"]")), "{error}");
    let wire = claim(&NormalizedMaterialCardPack::from_bytes(GOLD).unwrap(), property);
    let (lo, hi) = wire.validity.bounds()["wire_diameter"];
    assert_eq!(lo, hi, "one specimen diameter");
    for object in overstressed(&error) {
        assert_eq!(field(&object, "limitAtDiameterM"), Some(lo), "{object}");
        assert!(object.contains("\"limitExtrapolated\":true"), "{object}");
    }
    let (lowest, highest) = check_refusal(&error, &gold, GOLD_CONCERT, scalar(GOLD, property), GOLD);
    println!("gold at concert pitch: 88/88 keys refuse; highest attainable pitch {lowest:.2}-{highest:.2} Hz");
}

// G1 (DONE-WHEN): at half the wire's strength the real gold construction admits every key, and A0, C4
// and C7 play. Each tension is half of the core's capacity, so B = EI pi^2/(T L^2) rises against the
// steel construction of the same geometry by (E_gold/E_steel)(T_steel/T_gold) > 1, and each partial's
// strain-rate loss is the gold core's computed thermoelastic loss alone (gold's anelastic excess has no
// source), against steel's thermoelastic loss plus its measured in-band excess.
#[test]
fn g1_real_gold_at_half_strength_admits_a0_c4_c7_and_moves_b_and_decay_as_computed() {
    let gold = Arc::new(materials(&[GOLD, GOLD], GOLD_HALF).unwrap());
    let steel = Arc::new(materials(&[STEEL, COPPER, QQ], STEEL_WOUND).unwrap());
    let mut gold_desc = Description::embedded().unwrap();
    preparation::prepare_with_materials(&mut gold_desc, None, None, Some(&gold)).unwrap();
    let mut steel_desc = Description::embedded().unwrap();
    preparation::prepare_with_materials(&mut steel_desc, None, None, Some(&steel)).unwrap();
    let report = gold.strength_report(&gold_desc.courses).unwrap();
    assert!(report.keys.iter().all(|k| k.utilization < 1.0));
    assert!(report.refusal_json().contains(GEOMETRY), "the gold-wound caveat");

    let uts = scalar(GOLD, "ultimate_tensile_strength_cold_rolled_wire");
    let (cores, wound) = rows(GOLD_HALF);
    assert_eq!(rows(STEEL_WOUND).0, cores, "one geometry");
    // Thermoelastic constants of each core: Delta = E alpha^2 T0 / (rho c), base time rho c r^2 / k.
    let constants = |pack: &[u8]| {
        let (rho, young) = (scalar(pack, "density"), scalar(pack, "young_modulus"));
        let (alpha, c, k) = (scalar(pack, "linear_thermal_expansion_coefficient"),
            scalar(pack, "specific_heat_capacity"), scalar(pack, "thermal_conductivity"));
        (young, young * alpha * alpha * T0 / (rho * c), rho * c / k)
    };
    let (gold_e, gold_delta, gold_base) = constants(GOLD);
    let (steel_e, steel_delta, steel_base) = constants(STEEL);
    let excess = claim(&NormalizedMaterialCardPack::from_bytes(STEEL).unwrap(), LOSS);
    let (band_lo, band_hi) = excess.validity.bounds()["frequency"];
    let PropertyValue::Scalar { value: tan_delta, .. } = excess.value else { panic!("excess loss tangent") };

    let banks = [&gold_desc, &steel_desc].map(|desc|
        Bank::new_with_string_damping(&desc.courses, &desc.board, 48_000, 0.45 * 48_000.0, 8, true, None, false)
            .unwrap());
    let losses = banks.each_ref().map(|bank| bank_string_losses(bank).0.to_vec());
    let embedded = scale();
    for key in [21_u8, 60, 96] {
        let i = embedded.iter().position(|c| c.midi == key).unwrap();
        let (g, s) = (&gold_desc.courses[i], &steel_desc.courses[i]);
        let r = cores[&key];
        let ei = |e: f64| e * PI * r.powi(4) / 4.0;
        near(g.tension_n, 0.5 * uts * PI * r * r, 1e-12);
        near(g.flexural_rigidity_nm2, ei(gold_e), 1e-12);
        near(s.flexural_rigidity_nm2, ei(steel_e), 1e-12);
        // Steel keeps the installed fundamental: T = mu (2 L f)^2 - EI (pi / L)^2.
        let f1 = embedded[i].partial_hz(1, embedded[i].tension_n);
        near(s.tension_n, s.linear_density_kg_m * (2.0 * s.length_m * f1).powi(2)
            - ei(steel_e) * (PI / s.length_m).powi(2), 1e-9);
        if !wound.contains(&key) {
            near(g.linear_density_kg_m, scalar(GOLD, "density") * PI * r * r, 1e-12);
        }
        let b = |c: &Course| c.flexural_rigidity_nm2 * PI * PI / (c.tension_n * c.length_m * c.length_m);
        let predicted = (gold_e / steel_e) * (s.tension_n / g.tension_n);
        near(b(g) / b(s), predicted, 1e-12);
        assert!(predicted > 1.0, "key {key}: B ratio {predicted}");
        println!("key {key}: gold f1 {:.2} Hz (steel {:.2}); B {:.3e} vs {:.3e}, x{predicted:.3}",
            g.partial_hz(1, g.tension_n), f1, b(g), b(s));
        // Each partial's strain-rate loss against its prediction, gold and steel, on every unison
        // member at its own detuned tension.
        for (bank_index, (course, delta, base, young, measured)) in
            [(g, gold_delta, gold_base, gold_e, false), (s, steel_delta, steel_base, steel_e, true)]
            .into_iter().enumerate() {
            let members = banks[bank_index].strings.iter()
                .filter(|string| string.course == i && string.polarization == 0 && !string.duplex).collect::<Vec<_>>();
            assert_eq!(members.len() as f64, course.unison as f64, "key {key}");
            for string in members {
                let cents = (string.member as f64 - 0.5 * (course.unison - 1) as f64) * course.detune_cents;
                let tension = course.tension_at_cents(cents).unwrap();
                for (n, index) in string.modes.clone().enumerate().take(4) {
                    let k = (n + 1) as f64 * PI / course.length_m;
                    let hz = course.partial_hz(n + 1, tension);
                    let (omega, mu) = (TAU * hz, course.linear_density_kg_m);
                    let excess = if measured && (band_lo..=band_hi).contains(&hz) { tan_delta } else { 0.0 };
                    let rate = ei(young) * (delta * circular_relaxation_shape(omega, base * r * r) + excess)
                        / (mu * omega) * k.powi(4);
                    near(losses[bank_index][index].strain_rate, rate, 1e-9);
                }
            }
        }
        let first = |bank_index: usize| {
            let string = banks[bank_index].strings.iter().find(|string| string.course == i).unwrap();
            losses[bank_index][string.modes.start].strain_rate
        };
        println!("key {key}: partial-1 strain-rate loss gold {:.3e} 1/s, steel {:.3e} 1/s", first(0), first(1));
    }
    // A0, C4 and C7 play from the gold strings.
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), T0, 101_325.0).unwrap();
    let mut audio = PianoAudio::new_with_materials(&gold_desc, Config::live(48_000), 1, gas, gold.clone()).unwrap();
    for key in [21, 60, 96] {
        audio.note_on(key, 1.0).unwrap();
        let mut pcm = vec![0.0_f32; 4096];
        audio.render(&mut pcm).unwrap();
        assert!(pcm.iter().all(|x| x.is_finite()) && pcm.iter().any(|x| *x != 0.0), "key {key}");
    }
}

// G1 (ruling 6326): drawn silver's one sourced strength states no wire diameter. The real silver
// construction declares diameter-free, so every receipt says the strength is not diameter-resolved,
// and all 88 keys refuse at concert pitch with their highest pitch. With no sourced yield, a
// yield-fraction construction on silver refuses by name; without the flag the claim refuses.
#[test]
fn g1_real_silver_refuses_every_key_on_a_strength_that_is_not_diameter_resolved() {
    let packs: [&[u8]; 3] = [SILVER, SILVER, SILVER_UTS];
    let silver = materials(&packs, SILVER_CONCERT).unwrap();
    assert!(silver.diameter_free() && !silver.extrapolate_diameter());
    assert_eq!(silver.strength_pack(), Some(2));
    let property = "ultimate_tensile_strength_drawn_hard";
    let error = silver.apply(&scale()).unwrap_err();
    let caveat = format!("strength not diameter-resolved; '{property}' from silver-drawn-hard-uts-nbs-c101 states no \
        wire diameter, diameter effect unbounded");
    assert!(error.contains(&format!("\"limitCaveats\":[\"{caveat}\"]")), "{error}");
    for object in overstressed(&error) {
        assert_eq!(field(&object, "limitAtDiameterM"), None, "{object}");
    }
    let (lowest, highest) = check_refusal(&error, &silver, SILVER_CONCERT, scalar(SILVER_UTS, property), SILVER);
    println!("silver at concert pitch: 88/88 keys refuse; highest attainable pitch {lowest:.2}-{highest:.2} Hz");

    let half = SILVER_CONCERT.replace("\ntension,target-pitch\n", "\ntension,yield-fraction,0.5\n");
    assert_ne!(half, SILVER_CONCERT);
    let error = materials(&packs, &half).unwrap_err();
    assert!(error.contains("yield-fraction") && error.contains("silver-drawn-hard-uts-nbs-c101")
        && error.contains("nothing back-solves a yield"), "{error}");
    let unflagged = SILVER_CONCERT.replace(",diameter-free,", ",");
    assert_ne!(unflagged, SILVER_CONCERT);
    let error = materials(&packs, &unflagged).unwrap_err();
    assert!(error.contains("state no wire diameter") && error.contains("diameter-free"), "{error}");
}

// G1 (ruling 6339): every core of the steel construction lies inside QQ-W-470B's transcribed rows, so
// each takes the minimum strength at its own diameter by the declared linear rule, with no
// extrapolation, and the concert scale holds. Planted negative (DONE-WHEN) on the real data: a core
// drawn below the transcribed rows refuses with its numbers.
#[test]
fn g1_real_steel_cores_take_qq_w_470b_at_their_own_diameters_and_one_outside_refuses() {
    let packs: [&[u8]; 3] = [STEEL, COPPER, QQ];
    let steel = materials(&packs, STEEL_WOUND).unwrap();
    let knots = match claim(&NormalizedMaterialCardPack::from_bytes(QQ).unwrap(), "minimum_tensile_strength").value {
        PropertyValue::Curve { abscissa, knots, .. } => { assert_eq!(abscissa, "wire_diameter"); knots }
        other => panic!("not a curve: {other:?}"),
    };
    let (lo, hi) = (knots[0].0, knots[knots.len() - 1].0);
    let (cores, _) = rows(STEEL_WOUND);
    assert_eq!(cores.len(), 88);
    let (mut smallest, mut largest) = ((0_u8, f64::INFINITY), (0_u8, 0.0_f64));
    for (&key, &r) in &cores {
        let d = 2.0 * r;
        if d < smallest.1 { smallest = (key, d); }
        if d > largest.1 { largest = (key, d); }
        assert!((lo..=hi).contains(&d), "key {key}: {d} m outside {lo}-{hi} m");
        let receipt = &steel.limit_receipts()[&key][0];
        assert_eq!((receipt.pack, receipt.diameter_m, receipt.evaluated_at_m), (2, d, Some(d)));
        assert!(!receipt.extrapolated && receipt.caveat.is_none());
        let upper = knots.iter().position(|knot| knot.0 >= d).unwrap();
        if knots[upper].0 == d {
            assert_eq!(receipt.rows, vec![knots[upper]]);
            assert_eq!(receipt.limit_pa, knots[upper].1);
        } else {
            let (a, b) = (knots[upper - 1], knots[upper]);
            assert_eq!(receipt.rows, vec![a, b]);
            near(receipt.limit_pa, a.1 + (b.1 - a.1) * (d - a.0) / (b.0 - a.0), 1e-12);
        }
    }
    println!("steel cores {:.5} mm (key {}) to {:.5} mm (key {}), inside QQ-W-470B rows {:.3}-{:.3} mm",
        smallest.1 * 1e3, smallest.0, largest.1 * 1e3, largest.0, lo * 1e3, hi * 1e3);
    let courses = steel.apply(&scale()).unwrap();
    let report = steel.strength_report(&courses).unwrap();
    assert!(report.limit_caveats.is_empty() && report.refusal_json().contains(GEOMETRY));
    assert!(steel.json().contains("\"geometry\":{\"backSolvedKeys\":[21,40],"));
    let worst = report.keys.iter().max_by(|a, b| a.utilization.total_cmp(&b.utilization)).unwrap();
    assert!(worst.utilization < 1.0);
    println!("steel at concert pitch: highest utilization {:.3} of the QQ-W-470B minimum, key {}",
        worst.utilization, worst.key);

    let planted = STEEL_WOUND.replace("\n108,core,0.000367,,,0", "\n108,core,0.0002,,,0");
    assert_ne!(planted, STEEL_WOUND);
    let error = materials(&packs, &planted).unwrap_err();
    assert!(error.contains("key 108: its 0.4000 mm core lies outside the 0.5080-2.4130 mm wire diameters")
        && error.contains("missing data"), "{error}");
}
