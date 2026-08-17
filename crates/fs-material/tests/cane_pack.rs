//! g2 physicality gates for the cane (Arundo donax) matdb packs
//! (music bead `frankensim-music-v8-root-3ez8g.3.6`) — the reed card's
//! material population, licensing-first (three CC-BY primaries).
//!
//! The gates encode what CANNOT be wrong regardless of specimen:
//! positivity, the fiber anisotropy direction (E_L >> E_T), the
//! moisture law (wet is softer than dry on BOTH axes), cross-source
//! range OVERLAP where two sources measure the same zone, and the
//! blank-vs-heel spread being a DISJOINT-zones fact (cortex gradient,
//! recorded never averaged), plus schema hygiene (every scalar carries
//! an uncertainty row and at least one validity row).

use std::collections::BTreeMap;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

struct Pack {
    scalars: BTreeMap<String, f64>,
    uncertainty_rows: Vec<String>,
    validity_rows: Vec<String>,
    observation_notes: String,
}

fn load_pack(dir: &str) -> Pack {
    let root = repo_root();
    let manifest =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/manifest.tsv")))
            .expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v1"));
    assert!(
        manifest.contains("CC-BY-4.0"),
        "{dir}: the cane packs are licensing-first CC-BY"
    );
    assert!(manifest.contains("https://doi.org/"), "{dir}: DOI required");
    let text =
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{dir}/properties.tsv")))
            .expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1"));
    let mut pack = Pack {
        scalars: BTreeMap::new(),
        uncertainty_rows: Vec::new(),
        validity_rows: Vec::new(),
        observation_notes: String::new(),
    };
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "scalar" => {
                let value: f64 = cols[4].parse().expect("scalar value");
                pack.scalars.insert(cols[3].to_string(), value);
            }
            "uncertainty" => pack.uncertainty_rows.push(cols[1].to_string()),
            "validity" => pack.validity_rows.push(cols[1].to_string()),
            "observation" => pack.observation_notes.push_str(cols[4]),
            other => panic!("unknown row kind {other}"),
        }
    }
    // Schema hygiene: every scalar row-id has an uncertainty row and at
    // least one validity row.
    let scalar_ids: Vec<String> = text
        .lines()
        .filter(|l| l.starts_with("scalar\t"))
        .map(|l| l.split('\t').nth(1).expect("id").to_string())
        .collect();
    for id in &scalar_ids {
        assert!(
            pack.uncertainty_rows.iter().any(|u| u == id),
            "{dir}: scalar {id} lacks an uncertainty row"
        );
        assert!(
            pack.validity_rows.iter().any(|v| v == id),
            "{dir}: scalar {id} lacks a validity row"
        );
    }
    for (name, &value) in &pack.scalars {
        assert!(
            value.is_finite() && value > 0.0,
            "{dir}: {name} must be positive and finite"
        );
    }
    pack
}

#[test]
fn g2_cane_packs_pass_physicality_and_cross_source_gates() {
    let blank = load_pack("cane-arundo-reedblank-mdpi-ma18122759");
    let heel = load_pack("cane-arundo-reedheel-mdpi-ma13204566");
    let damping = load_pack("cane-arundo-damping-scielo-mr20170795");

    // Fiber anisotropy: E_L / E_T in the strongly-orthotropic band the
    // source itself states (10.0 +- 2.1 dry).
    let el_dry = heel.scalars["young_modulus_longitudinal"];
    let et_dry = heel.scalars["young_modulus_transverse"];
    let ratio_dry = el_dry / et_dry;
    assert!(
        (5.0..20.0).contains(&ratio_dry),
        "dry anisotropy ratio {ratio_dry:.1} outside the plausible band"
    );

    // The moisture law: water softens BOTH axes, transverse more.
    let el_wet = heel.scalars["young_modulus_longitudinal_water_soaked"];
    let et_wet = heel.scalars["young_modulus_transverse_water_soaked"];
    assert!(el_wet < el_dry, "wet must be softer along the fiber");
    assert!(et_wet < et_dry, "wet must be softer across the fiber");
    assert!(
        el_wet / et_wet > ratio_dry,
        "water widens the anisotropy (the source's own finding)"
    );

    // Cross-source agreement on the SAME zone: the heel's static E_L
    // band (5.0 +- 0.7 GPa) must OVERLAP the whole-wall DMA storage
    // range (5250..6250 MPa) — two independent labs, two methods.
    let heel_lo = el_dry - 700.0;
    let heel_hi = el_dry + 700.0;
    let dma_lo = damping.scalars["storage_modulus_printed_range_low"];
    let dma_hi = damping.scalars["storage_modulus_printed_range_high"];
    assert!(
        heel_hi >= dma_lo && dma_hi >= heel_lo,
        "the heel static band [{heel_lo},{heel_hi}] must overlap the DMA range [{dma_lo},{dma_hi}]"
    );

    // The blank-vs-heel spread is a DISJOINT-ZONES fact: cortex-weighted
    // blanks are far stiffer, and both observations must NAME the
    // radial gradient so nobody averages across zones.
    let blank_lo = blank.scalars["young_modulus_longitudinal_printed_range_low"];
    assert!(
        blank_lo > heel_hi,
        "blank (cortex) stiffness must sit above the heel band, disjoint"
    );
    for (pack, name) in [(&blank, "blank"), (&heel, "heel")] {
        assert!(
            pack.observation_notes.contains("gradient")
                || pack.observation_notes.contains("cortex"),
            "{name} observation must name the radial gradient / cortex zoning"
        );
    }

    // Density range sanity (dry blank material, cortex-weighted).
    let rho_lo = blank.scalars["density_printed_range_low"];
    let rho_hi = blank.scalars["density_printed_range_high"];
    assert!(rho_lo < rho_hi && (400.0..1200.0).contains(&rho_lo) && rho_hi < 1200.0);

    // Damping rows: the printed log-decrement range, low < high, and
    // the derived loss factor delta/pi lands in the woody 0.01..0.05
    // class (the conversion is OURS and the pack notes say so).
    let d_lo = damping.scalars["logarithmic_decrement_printed_range_low"];
    let d_hi = damping.scalars["logarithmic_decrement_printed_range_high"];
    assert!(d_lo < d_hi);
    let eta_lo = d_lo / core::f64::consts::PI;
    let eta_hi = d_hi / core::f64::consts::PI;
    assert!(
        eta_lo > 0.01 && eta_hi < 0.05,
        "derived loss factor [{eta_lo:.3},{eta_hi:.3}] outside the woody class"
    );
    assert!(
        damping.observation_notes.contains("delta/pi")
            || damping.observation_notes.contains("eta = delta/pi"),
        "the conversion ownership must be stated in the pack"
    );

    // The named absences stay named: no shear row anywhere (the hunt
    // found NO licensable shear number — asserting its absence keeps a
    // future guessed row from sneaking in without a source).
    for (pack, name) in [(&blank, "blank"), (&heel, "heel"), (&damping, "damping")] {
        assert!(
            !pack.scalars.keys().any(|k| k.contains("shear")),
            "{name}: a shear row appeared without a licensable source"
        );
    }
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"g2-cane-packs\",\"verdict\":\"pass\",\
         \"heel_el_dry_mpa\":{el_dry},\"anisotropy_dry\":{ratio_dry:.1},\
         \"blank_el_range_mpa\":[{blank_lo},{}],\"eta_derived\":[{eta_lo:.3},{eta_hi:.3}]}}",
        blank.scalars["young_modulus_longitudinal_printed_range_high"]
    );
}
