//! Apply a string construction (frankensim-piano-materials-v2) and compiled
//! canonical material packs (`xtask matdb-pack --material-state ...`) to the
//! embedded Model D scale's courses for the keys the construction names (native).
//!
//! Prints one JSON object: either the structured refusal (for example every
//! key a gold string cannot hold at concert pitch), or, per key, the line
//! density, tension and inharmonicity coefficient B = EI pi^2 / (T L^2) of
//! the embedded scale beside those of the material strings, with relative
//! differences and the maxima over all keys. Nothing is tuned to agree.
//!
//! usage: string_materials_check <construction.csv> <bending-loss property> <pack>...
//! The operating point is T = 293.15 K on axis `T`, 1000 Hz on axis `frequency`.
use fs_live_piano::{kernel::SCALE_CSV, preparation::materials::{StringMaterials, operating_point},
    refmods::geometry::{self, Course}};
use std::collections::BTreeSet;

fn inharmonicity(c: &Course) -> f64 {
    c.flexural_rigidity_nm2 * (std::f64::consts::PI / c.length_m).powi(2) / c.tension_n
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert!(args.len() >= 3, "usage: string_materials_check <construction.csv> <loss property> <pack>...");
    let construction = std::fs::read_to_string(&args[0]).expect("construction file");
    let packs: Vec<Vec<u8>> = args[2..].iter().map(|path| std::fs::read(path).expect("pack file")).collect();
    let keys: BTreeSet<u8> = construction.lines()
        .filter_map(|line| line.split('#').next())
        .filter(|line| line.contains(",core,") || line.contains(",wrap,"))
        .map(|line| line.split(',').next().unwrap().trim().parse().expect("row key"))
        .collect();
    let scale: Vec<Course> = geometry::read_scale(SCALE_CSV).expect("embedded scale")
        .into_iter().filter(|c| keys.contains(&c.midi)).collect();
    let point = operating_point("T", 293.15, "frequency", 1000.0).expect("operating point");
    let slices: Vec<&[u8]> = packs.iter().map(Vec::as_slice).collect();
    let materials = match StringMaterials::from_packs(&slices, &construction, &point, &args[1]) {
        Ok(materials) => materials,
        Err(e) => { println!("{{\"constructionRefused\":{:?}}}", e); std::process::exit(2); }
    };
    let bound = match materials.apply(&scale) {
        Ok(bound) => bound,
        Err(e) => { println!("{{\"applyRefused\":{:?}}}", e); std::process::exit(3); }
    };
    let report = materials.strength_report(&bound).expect("strength report");
    let (mut worst_mu, mut worst_t, mut worst_b, mut worst_u) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let rows = scale.iter().zip(&bound).zip(&report.keys).map(|((old, new), stress)| {
        let rel = |a: f64, b: f64| b / a - 1.0;
        let (mu, t, b) = (rel(old.linear_density_kg_m, new.linear_density_kg_m),
            rel(old.tension_n, new.tension_n), rel(inharmonicity(old), inharmonicity(new)));
        worst_mu = worst_mu.max(mu.abs()); worst_t = worst_t.max(t.abs()); worst_b = worst_b.max(b.abs());
        worst_u = worst_u.max(stress.utilization);
        format!(concat!("{{\"key\":{},\"muScaleKgM\":{:e},\"muMaterialKgM\":{:e},\"muRel\":{:e},",
            "\"tensionScaleN\":{:e},\"tensionMaterialN\":{:e},\"tensionRel\":{:e},",
            "\"bScale\":{:e},\"bMaterial\":{:e},\"bRel\":{:e},\"utilization\":{:e},\"fMaxHz\":{:e}}}"),
            old.midi, old.linear_density_kg_m, new.linear_density_kg_m, mu, old.tension_n, new.tension_n, t,
            inharmonicity(old), inharmonicity(new), b, stress.utilization, stress.f_max_hz)
    }).collect::<Vec<_>>().join(",\n");
    println!("{{\"keys\":{},\"maxAbsMuRel\":{:e},\"maxAbsTensionRel\":{:e},\"maxAbsBRel\":{:e},\"maxUtilization\":{:e},\n\"rows\":[\n{}]}}",
        bound.len(), worst_mu, worst_t, worst_b, worst_u, rows);
}
