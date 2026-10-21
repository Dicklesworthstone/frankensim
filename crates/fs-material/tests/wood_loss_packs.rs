//! Physicality gates for the directional wood-loss matdb packs (bead `frankensim-irp-v5-instrument-realism-pexfp.193`, plan item
//! M8-W2c): western redcedar and Norway spruce (Yano, Mukudai and Onishi 1990), Sitka spruce boards (Ono 1993, Ono 1996, Obataya et
//! al. 2001), Norway spruce plates and wedges (Viala, Placet and Cogan 2020), spruce and sycamore plates (Laemmlein et al. 2020;
//! Carlier et al. 2018), wavy sycamore (Alkadri et al. 2018), the graph-read axial damping of Spycher et al. 2008 Fig. 4, the Picea
//! rows of Bremaud et al. 2012, the compiled Engelmann row of Bremaud 2012 (JASA), seasoned spruce (Obataya et al. 2020) and the
//! moisture sensitivity of Bremaud and Gril 2021.
//!
//! The gates encode what the sources print or what cannot be wrong for air-dry wood: every loss factor lies in the small-damping
//! regime where tan delta, Q^-1 and twice the modal damping coincide; within one source and one specimen set the radial and the LR
//! shear loss exceed the longitudinal loss (each ratio is checked per source, never merged); Yano's printed acoustic conversion
//! efficiency sqrt(E/gamma^3)/tan delta and Ono's printed Q_R/Q_L follow from the transcribed means; the Spycher digitization
//! reproduces every claim from its pixel reads and keeps the paper's text values outside the figure's reading bound (the documented
//! text-versus-figure conflict); no pack encodes a frequency law (scalars only); moisture raises radial damping faster than
//! longitudinal damping. Schema hygiene: licence, DOI, facts-only terms, one uncertainty row and at least one validity row per scalar.

use std::collections::BTreeMap;
use std::f64::consts::PI;

const FACTS: &str = "Published-Scientific-Facts-Attribution-Retained";

const PACKS: [&str; 21] = [
    "cedar-western-red-loss-yano1990",
    "spruce-norway-loss-yano1990",
    "spruce-sitka-loss-ono1993",
    "spruce-sitka-loss-ono1996",
    "maple-sugar-loss-ono1996",
    "spruce-sitka-loss-obataya2001",
    "spruce-norway-loss-viala2020",
    "spruce-unstated-species-loss-lammlein2020",
    "maple-sycamore-loss-lammlein2020",
    "spruce-unstated-species-loss-carlier2018",
    "maple-unstated-species-loss-carlier2018",
    "maple-sycamore-wavy-loss-alkadri2018",
    "spruce-norway-loss-frequency-spycher2008",
    "maple-sycamore-loss-frequency-spycher2008",
    "spruce-norway-loss-bremaud2012",
    "spruce-sitka-loss-bremaud2012",
    "spruce-engelmann-loss-bremaud2012jasa",
    "spruce-sitka-loss-obataya2020",
    "spruce-norway-loss-obataya2020",
    "spruce-norway-loss-moisture-rate-bremaudgril2021",
    "maple-sycamore-loss-moisture-rate-bremaudgril2021",
];

struct Pack {
    scalars: BTreeMap<String, f64>,
    validity: BTreeMap<String, BTreeMap<String, (f64, f64)>>,
    observations: String,
}

impl Pack {
    fn get(&self, id: &str) -> f64 {
        *self
            .scalars
            .get(id)
            .unwrap_or_else(|| panic!("missing claim {id}"))
    }
    fn band(&self, id: &str, axis: &str) -> (f64, f64) {
        self.validity[id][axis]
    }
}

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn read(dir: &str, file: &str) -> String {
    std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/{file}")))
        .unwrap_or_else(|e| panic!("{dir}/{file}: {e}"))
}

/// Parse a pack after asserting the manifest licence, the DOI, the facts-only terms and schema hygiene.
fn load_pack(dir: &str) -> Pack {
    let manifest = read(dir, "manifest.tsv");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains(&format!("\npack_id\t{dir}\n")),
        "{dir}: pack id"
    );
    assert!(
        manifest.contains(&format!("\nlicense\t{FACTS}\n")),
        "{dir}: licence"
    );
    assert!(manifest.contains("doi:10."), "{dir}: DOI required");
    assert!(
        manifest.contains("facts"),
        "{dir}: facts-only terms must say the pack keeps facts"
    );
    let text = read(dir, "properties.tsv");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = Pack {
        scalars: BTreeMap::new(),
        validity: BTreeMap::new(),
        observations: String::new(),
    };
    let mut unc = Vec::new();
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                assert_eq!(cols.len(), 7, "{dir}: scalar row width");
                let value: f64 = cols[4].parse().expect("scalar value");
                assert!(
                    value.is_finite() && value > 0.0,
                    "{dir}: {} must be positive",
                    cols[1]
                );
                assert!(
                    pack.scalars.insert(cols[1].to_string(), value).is_none(),
                    "{dir}: duplicate claim {}",
                    cols[1]
                );
            }
            "uncertainty" => {
                assert_eq!(cols.len(), 7, "{dir}: uncertainty row width");
                unc.push(cols[1].to_string());
            }
            "validity" => {
                assert_eq!(cols.len(), 6, "{dir}: validity row width");
                let (lo, hi): (f64, f64) = (cols[3].parse().unwrap(), cols[4].parse().unwrap());
                assert!(lo <= hi, "{dir}: {} {} bounds", cols[1], cols[2]);
                pack.validity
                    .entry(cols[1].to_string())
                    .or_default()
                    .insert(cols[2].to_string(), (lo, hi));
            }
            "observation" => {
                assert_eq!(cols.len(), 5, "{dir}: observation row width");
                pack.observations.push_str(line);
            }
            // Different specimens at different frequencies are not a frequency law: these packs hold scalars only.
            other => panic!("{dir}: row kind {other} is not allowed in a wood-loss pack"),
        }
    }
    for id in pack.scalars.keys() {
        assert_eq!(
            unc.iter().filter(|u| *u == id).count(),
            1,
            "{dir}: scalar {id} needs exactly one uncertainty row"
        );
        assert!(
            pack.validity.contains_key(id),
            "{dir}: scalar {id} lacks a validity row"
        );
    }
    pack
}

/// Same-specimen anisotropy of one source: the larger loss must exceed the smaller one.
fn exceeds(dir: &str, p: &Pack, larger: &str, smaller: &str) -> f64 {
    let ratio = p.get(larger) / p.get(smaller);
    assert!(ratio > 1.0, "{dir}: {larger}/{smaller} = {ratio}");
    ratio
}

#[test]
fn every_wood_loss_pack_is_schema_clean_and_in_the_small_damping_regime() {
    for dir in PACKS {
        let p = load_pack(dir);
        for (id, v) in &p.scalars {
            if id.starts_with("loss-factor") && !id.ends_with("moisture-rate") {
                // tan delta = Q^-1 = 2 zeta = Lambda/pi only for loss << 0.1 (Bremaud et al. 2012, Sect. 2.4.1).
                assert!(*v < 0.1, "{dir}: {id} = {v}");
            }
            if id.starts_with("log-decrement") {
                assert!(*v / PI < 0.1, "{dir}: {id} = {v}");
            }
        }
    }
}

#[test]
fn radial_and_shear_loss_exceed_longitudinal_loss_within_each_source() {
    let mut ratios = Vec::new();
    for dir in [
        "cedar-western-red-loss-yano1990",
        "spruce-norway-loss-yano1990",
    ] {
        let p = load_pack(dir);
        ratios.push(exceeds(
            dir,
            &p,
            "loss-factor-radial",
            "loss-factor-longitudinal",
        ));
    }
    let p = load_pack("spruce-sitka-loss-ono1993");
    for b in ["1", "2", "3"] {
        ratios.push(exceeds(
            "ono1993",
            &p,
            &format!("loss-factor-radial-board{b}"),
            &format!("loss-factor-longitudinal-board{b}"),
        ));
    }
    // Ono 1996 prints flexural L and R and torsional LR loss of one board: L < LR < R for each board, spruce and maple alike.
    let p = load_pack("spruce-sitka-loss-ono1996");
    for b in ["1", "2"] {
        exceeds(
            "ono1996",
            &p,
            &format!("loss-factor-shear-lr-board{b}"),
            &format!("loss-factor-longitudinal-board{b}"),
        );
        ratios.push(exceeds(
            "ono1996",
            &p,
            &format!("loss-factor-radial-board{b}"),
            &format!("loss-factor-shear-lr-board{b}"),
        ));
    }
    let p = load_pack("maple-sugar-loss-ono1996");
    exceeds(
        "ono1996-maple",
        &p,
        "loss-factor-shear-lr",
        "loss-factor-longitudinal",
    );
    exceeds(
        "ono1996-maple",
        &p,
        "loss-factor-radial",
        "loss-factor-shear-lr",
    );
    let p = load_pack("spruce-sitka-loss-obataya2001");
    ratios.push(exceeds(
        "obataya2001",
        &p,
        "loss-factor-radial",
        "loss-factor-longitudinal",
    ));
    for dir in [
        "spruce-norway-loss-viala2020",
        "spruce-unstated-species-loss-carlier2018",
    ] {
        let p = load_pack(dir);
        exceeds(dir, &p, "loss-factor-shear-lr", "loss-factor-longitudinal");
        ratios.push(exceeds(
            dir,
            &p,
            "loss-factor-radial",
            "loss-factor-shear-lr",
        ));
    }
    let p = load_pack("maple-unstated-species-loss-carlier2018");
    exceeds(
        "carlier-maple",
        &p,
        "loss-factor-radial",
        "loss-factor-longitudinal",
    );
    exceeds(
        "carlier-maple",
        &p,
        "loss-factor-shear-lr",
        "loss-factor-longitudinal",
    );
    let p = load_pack("maple-sycamore-loss-lammlein2020");
    exceeds(
        "laemmlein-maple",
        &p,
        "loss-factor-shear-lr",
        "loss-factor-longitudinal",
    );
    // Alkadri prints plate means only: every plate's radial loss exceeds every plate's longitudinal loss.
    let p = load_pack("maple-sycamore-wavy-loss-alkadri2018");
    exceeds(
        "alkadri",
        &p,
        "loss-factor-longitudinal-high",
        "loss-factor-longitudinal-low",
    );
    exceeds(
        "alkadri",
        &p,
        "loss-factor-radial-high",
        "loss-factor-radial-low",
    );
    exceeds(
        "alkadri",
        &p,
        "loss-factor-radial-low",
        "loss-factor-longitudinal-high",
    );
    // Each ratio stays with its source; none of them is a property of wood in general.
    assert_eq!(ratios.len(), 10);
}

#[test]
fn laemmlein_spruce_radial_loss_is_the_flagged_outlier() {
    // The authors write that the radial value was probably estimated higher because of a close second mode. Its R/L ratio sits
    // above the largest softwood axial-to-radial damping ratio Bremaud, Gril and Thibaut 2011 (Table 1) collected, 4.19 over 14
    // softwoods; the pack keeps the number and the caveat together.
    let p = load_pack("spruce-unstated-species-loss-lammlein2020");
    let ratio = p.get("loss-factor-radial") / p.get("loss-factor-longitudinal");
    assert!(ratio > 4.19, "R/L {ratio}");
    assert!(p.observations.contains("probably estimated higher"));
    assert!(p.observations.contains("not a pure shear loss"));
}

#[test]
fn yano_printed_acoustic_conversion_efficiency_follows_from_the_transcribed_means() {
    // Table I p.1208 prints sqrt(E/gamma^3)/tan delta (cgs) per direction: 2.69e8 and 3.65e7 for Thuja plicata, 1.77e8 and 1.83e7 for
    // Picea abies. The printed column is a mean of per-specimen ratios, the transcription holds the means of E/gamma, gamma and
    // tan delta, so the two agree within a few percent; reading the R loss into the L row misses by more than a factor of two.
    let efficiency = |p: &Pack, d: &str, loss: &str| {
        let specific = p.get(&format!("specific-dynamic-modulus-{d}")) * 1.0e10; // GPa per unit SG -> dyn/cm2
        let gamma = p.get(&format!("specific-gravity-{d}-specimens"));
        (specific / (gamma * gamma)).sqrt() / p.get(loss)
    };
    for (dir, l_printed, r_printed) in [
        ("cedar-western-red-loss-yano1990", 2.69e8, 3.65e7),
        ("spruce-norway-loss-yano1990", 1.77e8, 1.83e7),
    ] {
        let p = load_pack(dir);
        let l = efficiency(&p, "longitudinal", "loss-factor-longitudinal");
        let r = efficiency(&p, "radial", "loss-factor-radial");
        assert!((l / l_printed - 1.0).abs() < 0.03, "{dir}: L {l}");
        assert!((r / r_printed - 1.0).abs() < 0.03, "{dir}: R {r}");
        let swapped = efficiency(&p, "longitudinal", "loss-factor-radial");
        assert!(l_printed / swapped > 2.0, "{dir}: swapped {swapped}");
        let (mc, _) = p.band("loss-factor-longitudinal", "moisture_content");
        assert!(
            (8.0..11.0).contains(&mc),
            "{dir}: printed MC about 8.3 or 10.2 percent"
        );
    }
}

#[test]
fn ono_printed_radial_to_longitudinal_ratios_follow_from_the_transcribed_rows() {
    // Table 1 p.400 prints Q_R^-1/Q_L^-1 for each long-board average. With Q printed to two significant digits the ratio carries a
    // print resolution of about 0.02.
    let p = load_pack("spruce-sitka-loss-ono1993");
    for (b, printed) in [("1", 2.25), ("2", 2.03), ("3", 2.12)] {
        let ratio = p.get(&format!("loss-factor-radial-board{b}"))
            / p.get(&format!("loss-factor-longitudinal-board{b}"));
        assert!((ratio - printed).abs() <= 0.02, "board {b}: {ratio}");
        let (fl, _) = p.band(&format!("loss-factor-longitudinal-board{b}"), "frequency");
        let (fr, _) = p.band(&format!("loss-factor-radial-board{b}"), "frequency");
        // A 2 mm spruce square bends far more easily across the grain: the R mode sits well below the L mode.
        assert!(fr < 0.5 * fl, "board {b}: f_R {fr} f_L {fl}");
    }
}

fn digitized(sp: &str, n: &str) -> (f64, f64) {
    // Points read off Fig. 4 for the conflict check below: (claim value, frequency centre of the validity band).
    let dir = format!("{sp}-loss-frequency-spycher2008");
    let p = load_pack(&dir);
    let id = format!("log-decrement-longitudinal-point{n}");
    let (lo, hi) = p.band(&id, "frequency");
    (p.get(&id), (lo * hi).sqrt())
}

#[test]
fn spycher_fig4_digitization_reproduces_every_claim() {
    for dir in [
        "spruce-norway-loss-frequency-spycher2008",
        "maple-sycamore-loss-frequency-spycher2008",
    ] {
        let p = load_pack(dir);
        assert!(
            p.observations
                .contains("GRAPH-READ in properties.tsv; digitization error bounded only")
        );
        assert!(p.observations.contains("NOT a frequency law"));
        let table = read(dir, "digitization.tsv");
        let lines: Vec<&str> = table.lines().collect();
        assert_eq!(lines[0], "frankensim.figure-extraction.v1");
        let header = lines
            .iter()
            .position(|l| l.starts_with("claim\t"))
            .expect("header");
        let mut reads = 0;
        for line in &lines[header + 1..] {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 13);
            let n = |i: usize| f[i].parse::<f64>().unwrap();
            let (x0, x1, rx, xr) = (n(3), n(4), n(5), n(6));
            let (y0, y1, ry, yr) = (n(7), n(8), n(9), n(10));
            let (nominal, half) = (n(11), n(12));
            assert!(xr >= 0.5 && yr >= 0.5, "{dir}: reading radius");
            // Log axes: x spans 100 Hz to 10000 Hz, y spans 0.10 down to 0.01.
            let freq = |x: f64| 10f64.powf(2.0 + 2.0 * (x - x0) / (x1 - x0));
            let delta = |y: f64| 10f64.powf(-1.0 - (y - y0) / (y1 - y0));
            for y in [ry - yr, ry, ry + yr] {
                assert!(
                    (delta(y) - nominal).abs() <= half + 1e-9,
                    "{dir} {}: pass {} delta {}",
                    f[0],
                    f[1],
                    delta(y)
                );
            }
            assert_eq!(p.get(f[0]), nominal, "{dir}: claim equals the nominal read");
            let (flo, fhi) = p.band(f[0], "frequency");
            for x in [rx - xr, rx + xr] {
                assert!(
                    (flo..=fhi).contains(&freq(x)),
                    "{dir} {}: f {}",
                    f[0],
                    freq(x)
                );
            }
            reads += 1;
        }
        assert_eq!(reads, 9, "{dir}: three points x three passes");
        // The source states delta rose with resonance frequency (p.336): ordered bands, rising delta.
        let mut prev = (0.0, 0.0);
        for k in ["1", "2", "3"] {
            let id = format!("log-decrement-longitudinal-point{k}");
            let (lo, hi) = p.band(&id, "frequency");
            assert!(lo > prev.0 && p.get(&id) > prev.1, "{dir}: point {k}");
            prev = (hi, p.get(&id));
        }
    }
    // The documented conflict: the text values (p.335) sit outside the reading bound of the figure point at that frequency.
    for (sp, text_delta, text_hz) in [
        ("spruce-norway", 0.032, 860.0),
        ("maple-sycamore", 0.043, 540.0),
    ] {
        let (fig, hz) = digitized(sp, "1");
        assert!(
            (hz / text_hz - 1.0).abs() < 0.05,
            "{sp}: same frequency {hz}"
        );
        assert!(
            text_delta / fig > 1.2,
            "{sp}: text {text_delta} figure {fig}"
        );
    }
}

#[test]
fn bremaud_picea_rows_follow_the_standard_damping_relation() {
    let abies = load_pack("spruce-norway-loss-bremaud2012");
    let sitka = load_pack("spruce-sitka-loss-bremaud2012");
    for (dir, p) in [("abies", &abies), ("sitka", &sitka)] {
        let derived =
            p.get("specific-dynamic-modulus-longitudinal") * p.get("specific-gravity") * 1000.0;
        assert!(
            (p.get("dynamic-young-modulus-longitudinal") / derived - 1.0).abs() < 1.0e-9,
            "{dir}: E' = (E'/gamma) x gamma"
        );
    }
    // Table III: the stiffer Sitka rows carry the lower damping, as the tan delta - E'/gamma relation the paper tests says.
    assert!(
        sitka.get("specific-dynamic-modulus-longitudinal")
            > abies.get("specific-dynamic-modulus-longitudinal")
    );
    assert!(sitka.get("loss-factor-longitudinal") < abies.get("loss-factor-longitudinal"));
}

#[test]
fn seasoned_spruce_rows_carry_their_printed_moisture() {
    for (dir, mc) in [
        ("spruce-sitka-loss-obataya2020", 11.7),
        ("spruce-norway-loss-obataya2020", 10.9),
    ] {
        let p = load_pack(dir);
        assert_eq!(
            p.band("loss-factor-longitudinal", "moisture_content"),
            (mc, mc)
        );
        // Table II prints Q_L^-1 between 5.9e-3 and 7.7e-3 over all 16 lots.
        assert!((5.9e-3..=7.7e-3).contains(&p.get("loss-factor-longitudinal")));
    }
}

#[test]
fn moisture_raises_radial_damping_faster_than_longitudinal_damping() {
    for dir in [
        "spruce-norway-loss-moisture-rate-bremaudgril2021",
        "maple-sycamore-loss-moisture-rate-bremaudgril2021",
    ] {
        let p = load_pack(dir);
        let (l, r) = (
            p.get("loss-factor-longitudinal-moisture-rate"),
            p.get("loss-factor-radial-moisture-rate"),
        );
        assert!(r > 2.0 * l && l > 0.0, "{dir}: R {r} L {l}");
        let (lo, hi) = p.band("loss-factor-radial-moisture-rate", "moisture_content");
        // Above the damping minimum near 5.6 percent EMC only.
        assert!(lo >= 5.5 && hi <= 17.4, "{dir}: {lo}..{hi}");
    }
}
