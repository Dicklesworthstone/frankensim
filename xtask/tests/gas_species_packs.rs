//! Helium, neon and sulfur hexafluoride arrive as source data (IRP L9
//! dispatch: the gas half of pexfp.362).
//!
//! Sources: NASA/TP-2002-211556 (McBride, Zehe and Gordon 2002; a U.S.
//! Government work) Appendix B rows and Appendix D coefficients; Hurly and
//! Mehl 2007, J. Res. NIST 112:75, Table 5 (helium zero-density viscosity and
//! thermal conductivity, facts with attribution); Quinones-Cisneros, Huber and
//! Deiters 2012 and Assael et al. 2012 (SF6 dilute-gas correlations, facts with
//! attribution). Oracles here are test-only; production resolves the compiled
//! packs.

use std::path::{Path, PathBuf};
use std::process::Command;

use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack,
    NormalizedSpeciesPack, PropertyValue, QueryPoint,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::packs::standard_state_from_packs;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

fn workspace_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent")
        .join(relative)
}

static NEXT_PACK: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

fn compile_bytes(slug: &str) -> Vec<u8> {
    let directory = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("gas-species-packs-{}", std::process::id()));
    std::fs::create_dir_all(&directory).expect("pack directory");
    // Tests run in parallel threads of one process: every compile gets its own file.
    let sequence = NEXT_PACK.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let output = directory.join(format!("{slug}-{sequence}.pack"));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-pack")
        .arg("--manifest")
        .arg(workspace_path(&format!(
            "data/matdb/seed-v1/{slug}/manifest.tsv"
        )))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace_path(""))
        .output()
        .expect("run xtask matdb-pack");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    std::fs::read(&output).expect("compiled pack")
}

struct Source {
    species: NormalizedSpeciesPack,
    nasa9: NormalizedModelPack,
    transport: Option<NormalizedMaterialCardPack>,
}

fn source(slug: &str, transport: Option<&str>) -> Source {
    let bytes = compile_bytes(slug);
    let species = NormalizedSpeciesPack::from_bytes(&bytes).unwrap();
    let species =
        NormalizedSpeciesPack::from_bytes_verified(species.content_hash(), &bytes).unwrap();
    let bytes = compile_bytes(&format!("{slug}-nasa9-mcbride2002"));
    let nasa9 = NormalizedModelPack::from_bytes(&bytes).unwrap();
    let nasa9 = NormalizedModelPack::from_bytes_verified(nasa9.content_hash(), &bytes).unwrap();
    let transport = transport.map(|transport_slug| {
        let bytes = compile_bytes(transport_slug);
        let claims = NormalizedPack::from_bytes(&bytes).unwrap();
        let claims = NormalizedPack::from_bytes_verified(claims.content_hash(), &bytes).unwrap();
        NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: species.association().species().as_str().into(),
                phase: "dilute gas".into(),
                process: "explicit primary correlation source association".into(),
                revision: 0,
            },
            claims,
        )
        .unwrap()
    });
    Source {
        species,
        nasa9,
        transport,
    }
}

/// Appendix B prints cp, H(298.15) and S(298.15) at 1 bar to three decimals.
/// The band is half the printed unit plus 1e-5 relative: NASA's
/// R = 8.314510 J/(mol K) (TP p. 4) differs from fs-thermochem's CODATA value
/// by 5.7e-6 relative.
#[test]
fn he_ne_sf6_nasa9_packs_reproduce_their_printed_reference_rows() {
    // (slug, id, Appendix B molar mass g/mol, cp J/(mol K), dH_f kJ/mol, S J/(mol K))
    let rows = [
        ("helium", "He", 4.00260, 20.786, 0.0, 126.154),
        ("neon", "Ne", 20.17970, 20.786, 0.0, 146.330),
        (
            "sulfur-hexafluoride",
            "SF6",
            146.05542,
            97.069,
            -1219.400,
            291.678,
        ),
    ];
    for (slug, id, mass, cp, dh, s) in rows {
        let source = source(slug, None);
        let association = source.species.association();
        assert_eq!(association.species().as_str(), id);
        assert!(
            (association.molar_mass() / (mass * 1e-3) - 1.0).abs() <= 1e-15,
            "{id} molar mass"
        );
        let model = standard_state_from_packs(&source.species, &source.nasa9).unwrap();
        let at = |t: f64| model.evaluate(Temperature::new(t)).unwrap().properties();
        let p = at(298.15);
        let band = |printed: f64, half_unit: f64| half_unit + 1e-5 * printed.abs();
        let (cp_m, h_m, s_m) = (p.cp().value(), p.h().value() * 1e-3, p.s().value());
        assert!(
            (cp_m - cp).abs() <= band(cp, 0.0005),
            "{id} cp {cp_m} vs {cp}"
        );
        assert!((h_m - dh).abs() <= band(dh, 0.0005), "{id} H {h_m} vs {dh}");
        assert!((s_m - s).abs() <= band(s, 0.0005), "{id} S {s_m} vs {s}");
        // The two coefficient regions meet at 1000 K.
        let (below, above) = (at(1000.0 - 1e-9), at(1000.0 + 1e-9));
        assert!(
            (below.cp().value() / above.cp().value() - 1.0).abs() < 1e-6,
            "{id} cp seam"
        );
        println!(
            "{{\"suite\":\"gas-species-packs\",\"species\":\"{id}\",\"cp_298\":{cp_m},\"printed_cp\":{cp},\"H_298_kJ\":{h_m},\"printed_H\":{dh},\"S_298\":{s_m},\"printed_S\":{s},\"verdict\":\"matched\"}}"
        );
    }
}

/// Hurly and Mehl 2007 Table 5 rows (T K, viscosity uPa s, conductivity
/// mW/(m K)), transcribed in full: the source prescribes cubic-spline
/// interpolation between them.
const HM2007: [(f64, f64, f64); 79] = [
    (1.0, 0.3292, 2.632),
    (1.2, 0.3405, 2.720),
    (1.4, 0.3583, 2.845),
    (1.6, 0.3844, 3.033),
    (1.8, 0.4183, 3.283),
    (2.0, 0.4586, 3.588),
    (2.25, 0.5161, 4.030),
    (2.5, 0.5791, 4.519),
    (2.75, 0.6457, 5.037),
    (3.0, 0.7141, 5.569),
    (3.5, 0.8511, 6.637),
    (4.0, 0.9834, 7.666),
    (4.5, 1.1078, 8.636),
    (5.0, 1.2239, 9.542),
    (6.0, 1.4339, 11.184),
    (7.0, 1.6209, 12.650),
    (8.0, 1.7919, 13.992),
    (9.0, 1.9514, 15.244),
    (10.0, 2.1023, 16.427),
    (11.0, 2.2463, 17.556),
    (12.0, 2.3846, 18.641),
    (14.0, 2.6472, 20.699),
    (16.0, 2.8947, 22.638),
    (18.0, 3.1299, 24.481),
    (20.0, 3.3550, 26.244),
    (22.0, 3.5715, 27.939),
    (23.0, 3.6768, 28.764),
    (24.0, 3.7804, 29.575),
    (25.0, 3.8824, 30.374),
    (26.0, 3.9828, 31.160),
    (28.0, 4.1795, 32.700),
    (30.0, 4.3709, 34.199),
    (35.0, 4.8301, 37.793),
    (40.0, 5.2659, 41.204),
    (45.0, 5.6828, 44.466),
    (50.0, 6.0837, 47.603),
    (60.0, 6.8465, 53.570),
    (70.0, 7.5673, 59.208),
    (80.0, 8.2549, 64.585),
    (90.0, 8.9149, 69.746),
    (100.0, 9.5519, 74.726),
    (120.0, 10.7689, 84.240),
    (140.0, 11.9245, 93.272),
    (160.0, 13.0310, 101.919),
    (180.0, 14.0968, 110.248),
    (200.0, 15.1284, 118.308),
    (225.0, 16.3769, 128.062),
    (250.0, 17.5862, 137.510),
    (273.16, 18.6765, 146.027),
    (275.0, 18.7620, 146.695),
    (298.15, 19.8245, 154.994),
    (300.0, 19.9084, 155.649),
    (325.0, 21.0288, 164.400),
    (350.0, 22.1260, 172.970),
    (375.0, 23.2024, 181.375),
    (400.0, 24.2598, 189.633),
    (450.0, 26.3241, 205.753),
    (500.0, 28.3298, 221.413),
    (600.0, 32.1959, 251.596),
    (700.0, 35.9049, 280.550),
    (800.0, 39.4880, 308.517),
    (900.0, 42.9671, 335.670),
    (1000.0, 46.3584, 362.135),
    (1200.0, 52.924, 413.367),
    (1400.0, 59.256, 462.771),
    (1600.0, 65.402, 510.71),
    (1800.0, 71.395, 557.46),
    (2000.0, 77.259, 603.20),
    (2500.0, 91.470, 714.02),
    (3000.0, 105.185, 820.95),
    (3500.0, 118.527, 924.97),
    (4000.0, 131.578, 1026.70),
    (4500.0, 144.394, 1126.60),
    (5000.0, 157.018, 1224.99),
    (6000.0, 181.810, 1418.18),
    (7000.0, 206.135, 1607.72),
    (8000.0, 230.116, 1794.54),
    (9000.0, 253.835, 1979.31),
    (10000.0, 277.355, 2162.52),
];

/// Natural cubic spline second derivatives (Thomas algorithm), the same
/// arithmetic as the pack generator.
fn spline(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    let n = xs.len();
    let h: Vec<f64> = xs.windows(2).map(|w| w[1] - w[0]).collect();
    let (mut sub, mut diag, mut sup, mut rhs) =
        (vec![0.0; n], vec![1.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 1..n - 1 {
        sub[i] = h[i - 1];
        diag[i] = 2.0 * (h[i - 1] + h[i]);
        sup[i] = h[i];
        rhs[i] = 6.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
    }
    for i in 1..n {
        let w = sub[i] / diag[i - 1];
        diag[i] -= w * sup[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    let mut m = vec![0.0; n];
    m[n - 1] = rhs[n - 1] / diag[n - 1];
    for i in (0..n - 1).rev() {
        m[i] = (rhs[i] - sup[i] * m[i + 1]) / diag[i];
    }
    m
}

fn spline_at(xs: &[f64], ys: &[f64], m: &[f64], t: f64) -> f64 {
    let mut k = 0;
    while k + 2 < xs.len() && xs[k + 1] < t {
        k += 1;
    }
    let h = xs[k + 1] - xs[k];
    let a = (xs[k + 1] - t) / h;
    let b = (t - xs[k]) / h;
    a * ys[k]
        + b * ys[k + 1]
        + ((a * a * a - a) * m[k] + (b * b * b - b) * m[k + 1]) * (h * h) / 6.0
}

fn query(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            p,
        )
        .unwrap()
}

/// Shared representation check of one sampled transport curve against its
/// source oracle: every stored knot to 1e-12, a 64/128/256 linear ladder that
/// converges at second order to <= 1e-4, and the actual resolver between all
/// stored knots to 1e-4.
fn representation_ladder(
    card: &NormalizedMaterialCardPack,
    species: &str,
    name: &str,
    low: f64,
    oracle: &dyn Fn(f64) -> f64,
) {
    let claim = card.card().claims().claims_for(name)[0].1;
    let PropertyValue::Curve { knots, .. } = &claim.value else {
        panic!("transport must be an explicit curve")
    };
    for &(t, value) in knots {
        assert!(
            (value / oracle(t) - 1.0).abs() <= 1e-12,
            "{species} {name} knot {t}: {value}"
        );
    }
    let mut errors = Vec::new();
    for intervals in [64_u32, 128, 256] {
        let step = (1000.0 - low) / f64::from(intervals);
        let mut error: f64 = 0.0;
        for cell in 0..intervals {
            let t0 = low + f64::from(cell) * step;
            let (a, b) = (oracle(t0), oracle(t0 + step));
            for fraction in [0.25, 0.5, 0.75] {
                let t = t0 + fraction * step;
                error = error.max(((a * (1.0 - fraction) + b * fraction) / oracle(t) - 1.0).abs());
            }
        }
        errors.push(error);
    }
    let order = (errors[1] / errors[2]).ln() / 2.0_f64.ln();
    assert!(
        errors[2] < errors[1] && errors[1] < errors[0],
        "{species} {name}: {errors:?}"
    );
    assert!(
        order >= 1.8,
        "{species} {name}: interpolation order {order}"
    );
    assert!(
        errors[2] <= 1e-4,
        "{species} {name}: final interpolation error {}",
        errors[2]
    );
    let requirement = ScalarPropertyRequirement::try_new(
        name,
        claim.key.dims(),
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    let mut worst: f64 = 0.0;
    for pair in knots.windows(2) {
        let t = 0.5 * (pair[0].0 + pair[1].0);
        let bundle = resolve_material_state_point(
            card.card(),
            &query(t, 101_325.0),
            std::slice::from_ref(&requirement),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let rel = (bundle.property(name).unwrap().value_si() / oracle(t) - 1.0).abs();
        worst = worst.max(rel);
        assert!(rel <= 1e-4, "{species} {name} resolver at {t}: {rel}");
    }
    println!(
        "{{\"suite\":\"gas-species-packs\",\"species\":\"{species}\",\"property\":\"{name}\",\"rungs\":[64,128,256],\"relative_errors\":{errors:?},\"order\":{order},\"resolver_worst\":{worst:e},\"tolerance\":0.0001,\"verdict\":\"matched\"}}"
    );
}

/// The helium transport pack is the source's own interpolation, sampled: the
/// spline through Table 5 reproduces every printed row inside 200..1000 K, and
/// the stored curve represents that spline (see `representation_ladder`).
#[test]
fn helium_transport_pack_is_the_sources_spline_sampled() {
    let he = source("helium", Some("helium-transport-hurly-mehl2007"));
    let xs: Vec<f64> = HM2007.iter().map(|r| r.0).collect();
    for (column, name, scale) in [
        (1, "dynamic_viscosity", 1e-6),
        (2, "thermal_conductivity", 1e-3),
    ] {
        let ys: Vec<f64> = HM2007
            .iter()
            .map(|r| if column == 1 { r.1 } else { r.2 })
            .collect();
        let m = spline(&xs, &ys);
        let oracle = |t: f64| spline_at(&xs, &ys, &m, t) * scale;
        for row in HM2007.iter().filter(|r| (200.0..=1000.0).contains(&r.0)) {
            let printed = if column == 1 { row.1 } else { row.2 } * scale;
            assert!(
                (oracle(row.0) / printed - 1.0).abs() <= 1e-12,
                "{name} row {}",
                row.0
            );
        }
        representation_ladder(he.transport.as_ref().unwrap(), "He", name, 200.0, &oracle);
    }
}

/// SF6 zero-density viscosity, Quinones-Cisneros, Huber and Deiters 2012
/// Eq. (6) and Table 2 (mPa s; signs read on the page image), Tc = 318.7232 K.
fn sf6_viscosity(t: f64) -> f64 {
    let tr = t / 318.7232;
    let half = tr.sqrt();
    let quarter = half.sqrt();
    (0.118561 - 0.378103 * quarter + 0.416428 * half - 0.165295 * (half * quarter) + 0.0245381 * tr)
        * 1e-3
}

/// SF6 zero-density thermal conductivity, Assael et al. 2012 Eq. (6) (mW/(m K)).
fn sf6_conductivity(t: f64) -> f64 {
    let num = 1461860.0 - 18539.4 * t + 77.7891 * t * t + 0.0241059 * t * t * t;
    let den = 29661.7 + 505.67 * t + t * t;
    num / den * 1e-3
}

/// The SF6 transport pack is the two printed correlations, sampled over
/// 225..1000 K (inside both sources' stated ranges); the 298.15 K viscosity is
/// logged for reference.
#[test]
fn sf6_transport_pack_is_the_printed_correlations_sampled() {
    let sf6 = source(
        "sulfur-hexafluoride",
        Some("sulfur-hexafluoride-transport-quinones2012-assael2012"),
    );
    let card = sf6.transport.as_ref().unwrap();
    representation_ladder(card, "SF6", "dynamic_viscosity", 225.0, &sf6_viscosity);
    representation_ladder(
        card,
        "SF6",
        "thermal_conductivity",
        225.0,
        &sf6_conductivity,
    );
    let eta = sf6_viscosity(298.15);
    println!(
        "{{\"suite\":\"gas-species-packs\",\"species\":\"SF6\",\"viscosity_298_15\":{eta:e}}}"
    );
}
