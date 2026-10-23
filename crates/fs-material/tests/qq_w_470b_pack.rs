//! Gates for the music-wire minimum tensile strength by diameter of Federal Specification QQ-W-470B (1965, a US
//! Government work; bead `frankensim-irp-v5-instrument-realism-pexfp.369`, U2 follow-up: per-diameter strength).
//!
//! The gates encode what cannot be wrong: the pack reproduces the printed Table II rows (diameter in inches, minimum
//! tensile strength in thousand p.s.i.) at the stated conversion factors; the minimum falls strictly as the diameter
//! grows, as drawing work-hardens the thinner wire; the curve is declared linear and its validity is exactly its knot
//! span, so nothing outside the transcribed rows is claimed; every value lies in the cold-drawn high-carbon band and
//! above the stress a steel string carries at concert pitch at the treble end of the scale.

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

const MPA_PER_KSI: f64 = 6.894757293168;
const MM_PER_IN: f64 = 25.4;

/// QQ-W-470B Table II as printed: diameter of wire [in], minimum tensile strength [thousand p.s.i.].
const TABLE_II: [(f64, f64); 25] = [
    (0.020, 350.0),
    (0.022, 345.0),
    (0.024, 341.0),
    (0.026, 337.0),
    (0.028, 333.0),
    (0.030, 330.0),
    (0.032, 327.0),
    (0.034, 324.0),
    (0.036, 321.0),
    (0.038, 318.0),
    (0.040, 315.0),
    (0.042, 313.0),
    (0.045, 309.0),
    (0.048, 306.0),
    (0.051, 303.0),
    (0.055, 300.0),
    (0.059, 296.0),
    (0.063, 293.0),
    (0.067, 290.0),
    (0.072, 287.0),
    (0.076, 284.0),
    (0.080, 282.0),
    (0.085, 279.0),
    (0.090, 276.0),
    (0.095, 274.0),
];

#[test]
fn qq_w_470b_pack_reproduces_table_ii_and_passes_physicality_gates() {
    let dir = repo_root().join("data/matdb/seed-v1/steel-music-wire-min-strength-qq-w-470b");
    let manifest = std::fs::read_to_string(dir.join("manifest.tsv")).expect("manifest");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"));
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    assert!(manifest.contains("QQ-W-470B") && manifest.contains("Table II"));
    assert!(manifest.contains("6696ea854ad811c5295e9bab33159ba19663c83799ac7e0592a82b0ad2e1c623"));
    assert!(manifest.contains("\naxis\tT\tT\tabsolute-temperature\n"));
    assert!(manifest.contains("\naxis\twire_diameter\twire_diameter\tdimensional\n"));
    let text = std::fs::read_to_string(dir.join("properties.tsv")).expect("properties");
    assert!(text.starts_with("frankensim.matdb-source.v1\n"));
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    let kinds: Vec<&str> = rows.iter().map(|row| row[0]).collect();
    assert_eq!(
        kinds,
        [
            "observation",
            "curve",
            "uncertainty",
            "validity",
            "validity"
        ]
    );
    let curve = &rows[1];
    assert_eq!(
        &curve[1..7],
        [
            "min-tensile-strength",
            "qqw470b-table2",
            "minimum_tensile_strength",
            "wire_diameter",
            "mm",
            "MPa"
        ]
    );
    assert_eq!(
        curve[8], "linear",
        "the interpolation rule between rows is declared"
    );
    let knots: Vec<(f64, f64)> = curve[7]
        .split(',')
        .map(|knot| {
            let (x, y) = knot.split_once(':').expect("x:y knot");
            (x.parse().expect("diameter"), y.parse().expect("strength"))
        })
        .collect();
    assert_eq!(knots.len(), TABLE_II.len());
    for (&(mm, mpa), &(inch, ksi)) in knots.iter().zip(&TABLE_II) {
        assert!(
            (mm - inch * MM_PER_IN).abs() < 1e-9,
            "{inch} in is {} mm, pack says {mm}",
            inch * MM_PER_IN
        );
        assert!(
            (mpa / MPA_PER_KSI / ksi - 1.0).abs() < 5e-8,
            "{inch} in: {mpa} MPa is not {ksi} ksi"
        );
    }
    for pair in knots.windows(2) {
        assert!(
            pair[1].0 > pair[0].0 && pair[1].1 < pair[0].1,
            "strength must fall as diameter grows: {pair:?}"
        );
    }
    assert_eq!(
        &rows[2][1..],
        ["min-tensile-strength", "unstated", "-", "-", "-", "-"]
    );
    let (lo, hi) = (knots[0].0, knots[knots.len() - 1].0);
    assert_eq!(&rows[3][1..3], ["min-tensile-strength", "wire_diameter"]);
    assert_eq!(rows[3][5], "mm");
    let (vlo, vhi): (f64, f64) = (
        rows[3][3].parse().expect("lower"),
        rows[3][4].parse().expect("upper"),
    );
    assert!(
        (vlo - lo).abs() < 1e-12 && (vhi - hi).abs() < 1e-12,
        "validity is exactly the knot span"
    );
    assert_eq!(
        &rows[4][1..],
        ["min-tensile-strength", "T", "293.15", "293.15", "K"]
    );
    // A plain string at fixed pitch carries 4 L^2 f^2 rho whatever its diameter: at C8 on the Model D (RT-0425: 0.052 m,
    // about 4186 Hz) steel of 7800 kg/m^3 needs about 1.5 GPa. Every guaranteed minimum exceeds it, and none leaves the
    // band of cold-drawn high-carbon steel wire.
    let concert_c8 = 4.0 * 0.052_f64.powi(2) * 4186.0_f64.powi(2) * 7800.0;
    for &(mm, mpa) in &knots {
        let pa = mpa * 1e6;
        assert!(pa > concert_c8 && pa < 3.5e9, "{mm} mm: {pa} Pa");
    }
    // The Model D's plain cores (RT-0425, MIDI 41-108: 0.734-1.117 mm) lie inside the transcribed span.
    assert!(lo < 0.734 && 1.117 < hi);
}
