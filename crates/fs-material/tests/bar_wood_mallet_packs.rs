//! Physicality gates for the tuned-bar wood and mallet matdb packs (lane P9, beads
//! `frankensim-irp-v5-instrument-realism-pexfp.267` and `.158`, plan items M8-W2c and M8-W2e). Every pack holds numerical facts
//! with attribution; no table or figure layout is copied.
//!
//! The gates encode what cannot be wrong whichever specimen was measured: positivity, bulk wood below the cell-wall density, each
//! paper's own identities (c^2 rho = E and A rho = c for the Honduran rosewood plates, tan delta = E''/E' and the printed shear-over-fibre
//! ratio for padauk within their printed rounding, Eq. 3 inverted for the boxwood mallet giving a reduced mass no larger than the mallet),
//! and one cross-source check: eta/E of Doutaut and Chaigne 1994 reproduces the eta = 2.35e-7 s that the repo's bar loss row takes from
//! Chaigne and Doutaut 1997. Schema hygiene: every scalar has an uncertainty and a validity row.

use std::collections::BTreeMap;

/// Bulk wood is a porous solid of cell wall; its density stays below the cell-wall substance (about 1500 kg/m^3).
const CELL_WALL_DENSITY: f64 = 1500.0;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn read(path: &str) -> String {
    std::fs::read_to_string(repo_root().join(path)).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// Scalars by property name, and the absolute one-sigma band of the scalars that have one.
struct Pack {
    value: BTreeMap<String, f64>,
    band: BTreeMap<String, f64>,
}

impl Pack {
    fn get(&self, property: &str) -> f64 {
        *self
            .value
            .get(property)
            .unwrap_or_else(|| panic!("missing {property}"))
    }
    fn sd(&self, property: &str) -> f64 {
        *self
            .band
            .get(property)
            .unwrap_or_else(|| panic!("missing band of {property}"))
    }
}

fn load_pack(dir: &str, cited: &str) -> Pack {
    let manifest = read(&format!("data/matdb/seed-v1/{dir}/manifest.tsv"));
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(manifest.contains("Published-Scientific-Facts-Attribution-Retained"));
    assert!(manifest.contains(cited), "{dir}: citation lacks {cited}");
    let text = read(&format!("data/matdb/seed-v1/{dir}/properties.tsv"));
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let (mut value, mut band, mut name_of) = (BTreeMap::new(), BTreeMap::new(), BTreeMap::new());
    let (mut ids, mut unc, mut val) = (Vec::new(), Vec::new(), Vec::new());
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let v: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    v.is_finite() && v > 0.0,
                    "{dir}: {} must be positive",
                    cols[3]
                );
                value.insert(cols[3].to_string(), v);
                name_of.insert(cols[1].to_string(), cols[3].to_string());
                ids.push(cols[1].to_string());
            }
            "uncertainty" => {
                if cols[2] == "absolute" {
                    let sd: f64 = cols[3].parse().expect("band");
                    band.insert(name_of[cols[1]].clone(), sd);
                }
                unc.push(cols[1].to_string());
            }
            "validity" => val.push(cols[1].to_string()),
            "observation" => {}
            other => panic!("unknown row kind {other}"),
        }
    }
    for id in &ids {
        assert!(
            unc.contains(id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(val.contains(id), "{dir}: scalar {id} lacks a validity row");
    }
    Pack { value, band }
}

/// Half a unit in the last printed digit of a value printed with `decimals` decimals.
fn half_digit(decimals: i32) -> f64 {
    0.5 * 10f64.powi(-decimals)
}

#[test]
fn honduran_rosewood_plates_pass_physicality_gates() {
    let p = load_pack(
        "rosewood-honduran-plates-danihelova-app2015",
        "Acta Physica Polonica A 127(1):106-109",
    );
    let (rho, e, c, a) = (
        p.get("density"),
        p.get("young_modulus_longitudinal") * 1e6,
        p.get("longitudinal_wave_speed"),
        p.get("acoustical_constant_sqrt_e_over_rho_cubed"),
    );
    assert!(rho < CELL_WALL_DENSITY);
    // The table's own identities (Eq. 2.2): c^2 rho = E within E's printed SD, and A rho = c within c's printed SD.
    let e_from_c = c * c * rho;
    assert!(
        (e_from_c - e).abs() <= p.sd("young_modulus_longitudinal") * 1e6,
        "c^2 rho = {e_from_c} against E = {e}"
    );
    assert!(
        (a * rho - c).abs() <= p.sd("longitudinal_wave_speed"),
        "A rho = {} against c = {c}",
        a * rho
    );
    for name in [
        "density",
        "young_modulus_longitudinal",
        "longitudinal_wave_speed",
        "acoustical_constant_sqrt_e_over_rho_cubed",
    ] {
        assert!(
            p.sd(name) < 0.1 * p.get(name),
            "{name}: specimen spread above 10 percent"
        );
    }
}

#[test]
fn boxwood_mallet_and_rosewood_bar_pass_physicality_gates() {
    let p = load_pack(
        "mallet-boxwood-rosewood-bar-doutaut-chaigne-jp4-1994",
        "doi 10.1051/jp4:19945131",
    );
    let card = read(
        "data/matdb/seed-v1/contact-boxwood-mallet-on-rosewood-bar-doutaut-chaigne1994/contact.tsv",
    );
    assert!(card.starts_with("frankensim.matdb-contact.v1"));
    // A card row by its record kind and, for validity rows, its axis.
    let row = |kind: &str, axis: Option<&str>| -> Vec<String> {
        card.lines()
            .map(|l| l.split('\t').map(str::to_string).collect::<Vec<_>>())
            .find(|c| c[0] == kind && axis.is_none_or(|a| c.get(1).is_some_and(|x| x == a)))
            .unwrap_or_else(|| panic!("card row {kind} {axis:?}"))
    };
    let law = row("law", None);
    let (k, alpha, chi): (f64, f64, f64) = (
        law[3].parse().unwrap(),
        law[5].parse().unwrap(),
        law[7].parse().unwrap(),
    );
    let force_hi: f64 = row("validity", Some("force_n"))[3].parse().unwrap();
    let speed_hi: f64 = row("validity", Some("velocity_m_s"))[3].parse().unwrap();
    // The paper's law is elastic: no contact loss was identified.
    assert!(k > 0.0 && chi.abs() < f64::EPSILON);
    assert!(
        (alpha - 1.5).abs() < f64::EPSILON,
        "Eq. 1 is Hertz, exponent 1.5"
    );

    let (m, tau, v) = (
        p.get("mallet_equivalent_mass"),
        p.get("measured_contact_duration"),
        p.get("initial_mallet_velocity"),
    );
    // Eq. 3 inverted: tau = 3.21 (mu^2/(K^2 V))^(1/5) gives the reduced mass of bar and mallet, which can be no larger than the
    // mallet's equivalent mass and, for a bar heavier than the mallet, no smaller than half of it.
    let mu = (tau.powi(5) * k * k * v / 3.21f64.powi(5)).sqrt();
    assert!(
        mu <= m && mu > 0.5 * m,
        "reduced mass {mu} kg against the mallet's {m} kg"
    );
    // The hand adds mass: the equivalent mass is above the head's.
    assert!(p.get("equivalent_mass_over_head_mass_roughly") > 1.0);
    // The card covers the identified impact: Hertz peak force at 2.1 m/s with mu at most the equivalent mass.
    let peak = k.powf(0.4) * (1.25 * m * v * v).powf(0.6);
    assert!(
        peak <= force_hi && v <= speed_hi,
        "peak {peak} N at {v} m/s outside the card"
    );

    let (e, rho, eta) = (
        p.get("young_modulus_longitudinal") * 1e6,
        p.get("density"),
        p.get("viscoelastic_constant_eta"),
    );
    assert!(rho < CELL_WALL_DENSITY);
    // Cross-source: eta/E of the 1994 Table I against the 1997 constant in the repo's bar loss row, within the printed precision of
    // the three numbers (eta 5.0e3, E 2.13e10, eta/E 2.35e-7: half a unit in the last digit of each, relative).
    let loss = read("data/instruments/descriptions/bars/rosewood-loss-chaigne-doutaut1997.tsv");
    let eta_1997: f64 = loss
        .lines()
        .find(|l| l.starts_with("scalar\tviscoelastic-constant"))
        .and_then(|l| l.split('\t').nth(4))
        .expect("1997 eta row")
        .parse()
        .unwrap();
    let rel = (eta / e - eta_1997).abs() / eta_1997;
    assert!(
        rel <= 0.05 / 5.0 + 0.005 / 2.13 + 0.005 / 2.35,
        "eta/E = {} s against {eta_1997} s",
        eta / e
    );
}

#[test]
fn padauk_heartwood_passes_physicality_gates() {
    let p = load_pack(
        "padauk-heartwood-bremaud-wst2010",
        "Wood Science and Technology 44(3):355-367",
    );
    // tan delta = E''/E' for the apparent and the corrected moduli, within the printed rounding of all three numbers.
    let tan = p.get("loss_tangent_along_fibre");
    let tan_band = (tan - half_digit(4), tan + half_digit(4));
    for (loss, storage) in [
        (
            "specific_loss_modulus_along_fibre_per_specific_gravity_apparent",
            "specific_storage_modulus_along_fibre_per_specific_gravity_apparent",
        ),
        (
            "specific_loss_modulus_along_fibre_per_specific_gravity_corrected_for_extractives",
            "specific_storage_modulus_along_fibre_per_specific_gravity_corrected_for_extractives",
        ),
    ] {
        let (l, s) = (p.get(loss), p.get(storage));
        // Loss moduli are printed to 1 MPa, storage moduli to 0.1 GPa.
        let (lo, hi) = ((l - 0.5) / (s + 50.0), (l + 0.5) / (s - 50.0));
        assert!(
            lo <= tan_band.1 && hi >= tan_band.0,
            "{loss}/{storage} in [{lo}, {hi}] misses tan delta {tan}"
        );
    }
    // The printed shear-over-fibre ratio agrees with the two printed loss tangents.
    let shear = p.get("loss_tangent_shear_fibre_tangential");
    let ratio = p.get("loss_tangent_shear_over_loss_tangent_fibre");
    let (lo, hi) = (
        (shear - half_digit(4)) / tan_band.1,
        (shear + half_digit(4)) / tan_band.0,
    );
    assert!(
        (lo - half_digit(2)..=hi + half_digit(2)).contains(&ratio),
        "ratio {ratio} outside [{lo}, {hi}]"
    );
    // Removing the extractives' mass from the specific gravity raises both specific moduli.
    assert!(
        p.get(
            "specific_storage_modulus_along_fibre_per_specific_gravity_corrected_for_extractives"
        ) > p.get("specific_storage_modulus_along_fibre_per_specific_gravity_apparent")
    );
    assert!(
        p.get("specific_loss_modulus_along_fibre_per_specific_gravity_corrected_for_extractives")
            > p.get("specific_loss_modulus_along_fibre_per_specific_gravity_apparent")
    );
    // Wood is stiffest along the fibre: a' (E3/G32 - 2 nu32) and b' (E3/E2) above 1; Poisson ratio inside (0, 0.5).
    assert!(
        p.get("anisotropy_a_prime_e3_over_g32_minus_2nu32") > 1.0
            && p.get("anisotropy_b_prime_e3_over_e2") > 1.0
    );
    assert!(p.get("poisson_ratio_32_about") < 0.5);
    assert!(p.get("fit_rmse_loss_tangent") < tan);
    assert!(
        (p.get("specific_gravity") + p.get("specific_gravity_printed_half_width")) * 1000.0
            < CELL_WALL_DENSITY
    );
}

#[test]
fn orff_xylophone_rosewood_and_madagascar_dalbergia_pass_physicality_gates() {
    let p = load_pack(
        "rosewood-orff-xylophone-bars-rucz-jasa2021",
        "J. Acoust. Soc. Am. 149(5):3200-3212",
    );
    assert!(p.get("density") + p.sd("density") < CELL_WALL_DENSITY);
    assert!(p.sd("density") < 0.1 * p.get("density"));

    let d = load_pack(
        "dalbergia-sp-madagascar-brancheriau2006",
        "doi 10.1007/s10086-005-0755-2",
    );
    assert!(d.get("density") < CELL_WALL_DENSITY);
    let (alpha, tan) = (
        d.get("first_flexural_mode_temporal_damping_alpha1"),
        d.get("loss_tangent_first_mode_eq3_two_alpha_over_pi_f"),
    );
    assert!(tan < 0.1);
    assert!(d.get("first_mode_temporal_damping_determination_accuracy") < 0.1);
    assert!(d.get("loss_tangent_determination_accuracy") < 0.1);
    // Eq. 3 (tan delta = 2 alpha1/(pi f1)) implies an f1 for this specimen; it must lie within three SD of the batch's printed
    // first-mode frequency (J. Wood Sci. p.2: mean 1008 Hz, SD 96 Hz over the 50 specimens of the same size).
    let f1 = 2.0 * alpha / (std::f64::consts::PI * tan);
    assert!(
        (1008.0 - 3.0 * 96.0..=1008.0 + 3.0 * 96.0).contains(&f1),
        "implied f1 {f1} Hz"
    );
}
