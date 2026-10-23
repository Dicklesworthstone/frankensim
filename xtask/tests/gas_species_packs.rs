//! Helium, neon and sulfur hexafluoride arrive as source data, and changing only
//! a wind performance's gas changes its sound the way the physics says
//! (IRP L9 dispatch: pexfp.18 and the gas half of pexfp.362).
//!
//! Sources: NASA/TP-2002-211556 (McBride, Zehe and Gordon 2002; a U.S.
//! Government work) Appendix B rows and Appendix D coefficients; Hurly and
//! Mehl 2007, J. Res. NIST 112:75, Table 5 (helium zero-density viscosity and
//! thermal conductivity, facts with attribution); Kestin et al. 1984 (neon
//! dilute-gas correlation, facts with attribution). Oracles here are test-only;
//! production resolves the compiled packs.

use std::path::{Path, PathBuf};
use std::process::Command;

use fs_couple::pcm_wav::observation::{ObservationPolicy, PressureRenderer};
use fs_couple::render::schedule::reed::ReedPerformance;
use fs_duct::{Duct, HoleState, LossModel, Segment, Termination, input_impedance};
use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack,
    NormalizedSpeciesPack, PropertyValue, QueryPoint,
};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::packs::{SpeciesPackRef, standard_state_from_packs};
use fs_math::c64::C64;
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

fn binding(source: &Source) -> SpeciesPackRef<'_> {
    SpeciesPackRef {
        species: &source.species,
        nasa9: Some(&source.nasa9),
        transport: Some((
            source
                .transport
                .as_ref()
                .expect("a gas record needs transport"),
            MaterialPropertySelection::SingleClaimOnly,
        )),
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

/// Neon zero-density viscosity and thermal conductivity, as Kestin et al. 1984
/// (JPCRD 13:229) program them. Every value below was read on the NIST page
/// images, because the text layer is a poor OCR:
/// - Table A1: k, N_A and R;
/// - Table A2: M = 20.179 g/mol;
/// - Table A3: sigma = 0.2755 nm, eps/k = 42.00 K;
/// - Table A6: rho* = 0.0784, V0* = 11.09e5;
/// - Eq. B7b for 1.2 <= T* <= 10, and Eq. B7c above;
/// - Eqs. C3, C8, C10, C12 and C13.
///
/// The source's tables are "not designed for linear interpolation", so the
/// correlation itself is the reference.
fn kestin_neon(t: f64) -> [f64; 2] {
    const K_B: f64 = 1.380662e-23;
    const N_A: f64 = 6.022045e23;
    const R: f64 = 8.31441;
    let (molar_mass, sigma, eps_k, rho, v0) = (20.179e-3, 0.2755e-9, 42.00, 0.0784, 11.09e5);
    let ts = t / eps_k;
    let x = ts.ln();
    let (ln_omega, slope) = if ts <= 10.0 {
        let c = [0.46641, -0.56991, 0.19591, -0.03879, 0.00259];
        (
            c[0] + c[1] * x + c[2] * x.powi(2) + c[3] * x.powi(3) + c[4] * x.powi(4),
            c[1] + 2.0 * c[2] * x + 3.0 * c[3] * x.powi(2) + 4.0 * c[4] * x.powi(3),
        )
    } else {
        let a10 = (v0 / 10.0_f64).ln();
        let q = (a10 * rho).powi(-2);
        let a2 = -33.0838 + q * (20.0862 + 72.1059 / a10 + (8.27648 / a10).powi(2));
        let a3 = 101.571 - q * (56.4472 + 286.393 / a10 + (17.7610 / a10).powi(2));
        let a4 = -87.7036 + q * (46.3130 + 277.146 / a10 + (19.0573 / a10).powi(2));
        let alpha = v0.ln() - x;
        let bracket = 1.04 + a2 * x.powi(-2) + a3 * x.powi(-3) + a4 * x.powi(-4);
        let dbracket = -2.0 * a2 * x.powi(-3) - 3.0 * a3 * x.powi(-4) - 4.0 * a4 * x.powi(-5);
        (
            2.0 * rho.ln() + 2.0 * alpha.ln() + bracket.ln(),
            -2.0 / alpha + dbracket / bracket,
        )
    };
    let e_star = 1.0 + 0.25 * slope;
    let f_eta = 1.0 + 3.0 / 196.0 * (8.0 * e_star - 7.0).powi(2);
    let f_lambda = 1.0 + 1.0 / 42.0 * (8.0 * e_star - 7.0).powi(2);
    let mass = molar_mass / N_A;
    let eta = 5.0 / 16.0 * (mass * K_B * t / core::f64::consts::PI).sqrt() * f_eta
        / (sigma * sigma * ln_omega.exp());
    [eta, 15.0 / 4.0 * eta * R / molar_mass * f_lambda / f_eta]
}

/// The neon transport pack is Kestin's correlation, sampled.
///
/// The programmed correlation is checked against the printed Table 2 (p. 238)
/// at every row inside 200..1000 K:
/// - viscosity to half the last printed digit (0.005 uPa s);
/// - conductivity within the source's stated +-0.7 %. It reproduces the table
///   only to about 3e-4 relative, a systematic and unresolved offset, logged
///   rather than hidden.
///
/// The stored curves represent the correlation (see `representation_ladder`).
#[test]
fn neon_transport_pack_is_the_kestin_correlation_sampled() {
    const TABLE2: [(f64, f64, f64); 19] = [
        (200.0, 24.29, 37.63),
        (250.0, 28.36, 43.96),
        (300.0, 32.10, 49.77),
        (273.15, 30.13, 46.71),
        (293.15, 31.60, 49.00),
        (313.15, 33.04, 51.22),
        (333.15, 34.43, 53.39),
        (353.15, 35.80, 55.51),
        (373.15, 37.13, 57.58),
        (423.15, 40.35, 62.57),
        (473.15, 43.42, 67.34),
        (523.15, 46.38, 71.93),
        (573.15, 49.25, 76.38),
        (623.15, 52.03, 80.69),
        (673.15, 54.74, 84.90),
        (723.15, 57.39, 89.00),
        (773.15, 59.99, 93.02),
        (873.15, 65.02, 100.83),
        (973.15, 69.89, 108.36),
    ];
    let mut worst = (0.0_f64, 0.0_f64);
    for (t, eta, lambda) in TABLE2 {
        let [e, l] = kestin_neon(t);
        let (de, dl) = ((e * 1e6 - eta).abs(), (l * 1e3 / lambda - 1.0).abs());
        assert!(
            de <= 0.005 + 1e-9,
            "Ne viscosity at {t} K: {} vs printed {eta}",
            e * 1e6
        );
        assert!(
            dl <= 0.007,
            "Ne conductivity at {t} K: {} vs printed {lambda}",
            l * 1e3
        );
        worst = (worst.0.max(de), worst.1.max(dl));
    }
    println!(
        "{{\"suite\":\"gas-species-packs\",\"species\":\"Ne\",\"table2_rows\":19,\"worst_viscosity_abs_uPa_s\":{},\"worst_conductivity_rel\":{:e},\"bands\":[0.005,0.007]}}",
        worst.0, worst.1
    );
    let ne = source("neon", Some("neon-transport-kestin1984"));
    let card = ne.transport.as_ref().unwrap();
    representation_ladder(card, "Ne", "dynamic_viscosity", 200.0, &|t| {
        kestin_neon(t)[0]
    });
    representation_ladder(card, "Ne", "thermal_conductivity", 200.0, &|t| {
        kestin_neon(t)[1]
    });
}

const D3: &str = "data/instruments/descriptions/clarinet/performances/BbClar.mf.D3.performance";
const AIR: &str = "ambient thermochemical 293.15 101325 mole 3 N2 0.7812 O2 0.2095 Ar 0.0093";
const HELIOX: &str = "ambient thermochemical 293.15 101325 mole 2 He 0.8 O2 0.2";

/// The performance's own bore, read with the same grammar the reed loader uses.
fn bore(text: &str) -> (Duct, Termination) {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let termination = loop {
        let line = lines.next().expect("termination row");
        if let Some(kind) = line.strip_prefix("termination ") {
            break match kind {
                "flanged" => Termination::FlangedOpen,
                "unflanged" => Termination::UnflangedOpen,
                "closed" => Termination::Closed,
                "ideal-open" => Termination::IdealOpen,
                other => panic!("termination {other}"),
            };
        }
    };
    let count: usize = lines
        .next()
        .unwrap()
        .strip_prefix("segments ")
        .unwrap()
        .parse()
        .unwrap();
    let segments = (0..count)
        .map(|_| {
            let words: Vec<&str> = lines.next().unwrap().split_ascii_whitespace().collect();
            let v = |i: usize| words[i].parse::<f64>().unwrap();
            match words[0] {
                "cylinder" => Segment::Cylinder {
                    radius: v(1),
                    length: v(2),
                },
                "cone" => Segment::Cone {
                    inlet_radius: v(1),
                    outlet_radius: v(2),
                    length: v(3),
                },
                "hole" => Segment::ToneHole {
                    hole_radius: v(1),
                    chimney_height: v(2),
                    bore_radius: v(3),
                    state: match v(4) {
                        s if s == 0.0 => HoleState::Closed,
                        s if s == 1.0 => HoleState::Open,
                        s => HoleState::Vent(s),
                    },
                },
                other => panic!("segment {other}"),
            }
        })
        .collect();
    (Duct { segments }, termination)
}

/// The lowest input-impedance maximum of the bore in this gas, located on a
/// 0.25 Hz scan and refined by golden-section search to 1e-6 Hz.
fn first_impedance_peak(duct: &Duct, termination: Termination, gas: &GasState) -> f64 {
    let z = |f: f64| {
        input_impedance(
            duct,
            gas,
            core::f64::consts::TAU * f,
            LossModel::Bessel,
            termination,
        )
        .unwrap()
        .impedance
        .abs()
    };
    let scan: Vec<(f64, f64)> = (0..4000)
        .map(|i| 50.0 + 0.25 * f64::from(i))
        .map(|f| (f, z(f)))
        .collect();
    let global = scan.iter().map(|p| p.1).fold(0.0, f64::max);
    let i = (1..scan.len() - 1)
        .find(|&i| {
            scan[i].1 > scan[i - 1].1 && scan[i].1 >= scan[i + 1].1 && scan[i].1 >= 0.3 * global
        })
        .expect("an input-impedance maximum");
    let (mut a, mut b) = (scan[i - 1].0, scan[i + 1].0);
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    while b - a > 1e-6 {
        let (c, d) = (b - g * (b - a), a + g * (b - a));
        if z(c) > z(d) { b = d } else { a = c }
    }
    0.5 * (a + b)
}

/// The strongest local maxima of a note's Hann-windowed spectrum over
/// `window` seconds, from a 1 Hz Goertzel scan over 50..4000 Hz.
fn spectral_peaks(pressure: &[f64], rate: f64, window: (f64, f64)) -> Vec<(f64, f64)> {
    let x = &pressure[(window.0 * rate) as usize..(window.1 * rate) as usize];
    let n = x.len() as f64;
    let windowed: Vec<f64> = x
        .iter()
        .enumerate()
        .map(|(i, v)| v * (0.5 - 0.5 * (core::f64::consts::TAU * i as f64 / n).cos()))
        .collect();
    let magnitude: Vec<(f64, f64)> = (50..=4000)
        .map(|f| {
            let coefficient = 2.0 * (core::f64::consts::TAU * f64::from(f) / rate).cos();
            let (mut s1, mut s2) = (0.0, 0.0);
            for v in &windowed {
                let s0 = v + coefficient * s1 - s2;
                s2 = s1;
                s1 = s0;
            }
            (
                f64::from(f),
                (s1 * s1 + s2 * s2 - coefficient * s1 * s2).max(0.0).sqrt(),
            )
        })
        .collect();
    let top = magnitude.iter().map(|m| m.1).fold(0.0, f64::max);
    let mut peaks: Vec<(f64, f64)> = magnitude
        .windows(3)
        .filter(|w| w[1].1 > w[0].1 && w[1].1 >= w[2].1 && w[1].1 >= 0.05 * top)
        .map(|w| (w[1].0, w[1].1 / top))
        .collect();
    peaks.sort_by(|a, b| b.1.total_cmp(&a.1));
    peaks.truncate(8);
    peaks
}

/// Local maxima of |Z_in| on a 1 Hz scan over 50..4000 Hz, with heights
/// relative to the largest.
fn impedance_maxima(duct: &Duct, termination: Termination, gas: &GasState) -> Vec<(f64, f64)> {
    let z: Vec<(f64, f64)> = (50..=4000)
        .map(|f| {
            let omega = core::f64::consts::TAU * f64::from(f);
            let response =
                input_impedance(duct, gas, omega, LossModel::Bessel, termination).unwrap();
            (f64::from(f), response.impedance.abs())
        })
        .collect();
    let top = z.iter().map(|m| m.1).fold(0.0, f64::max);
    z.windows(3)
        .filter(|w| w[1].1 > w[0].1 && w[1].1 >= w[2].1)
        .map(|w| (w[1].0, w[1].1 / top))
        .collect()
}

/// The played frequency over `window` seconds: the strongest spectral line,
/// refined on a 0.01 Hz Goertzel scan within 1 Hz of the coarse peak and a
/// parabola through the log magnitudes.
fn played_frequency(pressure: &[f64], rate: f64, window: (f64, f64)) -> f64 {
    let coarse = spectral_peaks(pressure, rate, window)[0].0;
    let x = &pressure[(window.0 * rate) as usize..(window.1 * rate) as usize];
    let n = x.len() as f64;
    let windowed: Vec<f64> = x
        .iter()
        .enumerate()
        .map(|(i, v)| v * (0.5 - 0.5 * (core::f64::consts::TAU * i as f64 / n).cos()))
        .collect();
    let magnitude = |f: f64| {
        let coefficient = 2.0 * (core::f64::consts::TAU * f / rate).cos();
        let (mut s1, mut s2) = (0.0, 0.0);
        for v in &windowed {
            let s0 = v + coefficient * s1 - s2;
            s2 = s1;
            s1 = s0;
        }
        (s1 * s1 + s2 * s2 - coefficient * s1 * s2).max(1e-300).ln()
    };
    let scan: Vec<(f64, f64)> = (0..=200)
        .map(|i| coarse - 1.0 + 0.01 * f64::from(i))
        .map(|f| (f, magnitude(f)))
        .collect();
    let i = (1..scan.len() - 1)
        .max_by(|&a, &b| scan[a].1.total_cmp(&scan[b].1))
        .unwrap();
    let (y0, y1, y2) = (scan[i - 1].1, scan[i].1, scan[i + 1].1);
    scan[i].0 + 0.01 * 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2)
}

/// The performance's reed row: rest opening, width, closing pressure, mass,
/// stiffness and damping ratio.
struct Reed {
    rest_opening: f64,
    width: f64,
    closing: f64,
    mass: f64,
    stiffness: f64,
    damping_ratio: f64,
}

fn reed(text: &str) -> Reed {
    let v: Vec<f64> = text
        .lines()
        .find_map(|line| line.strip_prefix("reed "))
        .expect("reed row")
        .split_ascii_whitespace()
        .map(|w| w.parse().unwrap())
        .collect();
    Reed {
        rest_opening: v[0],
        width: v[1],
        closing: v[2],
        mass: v[3],
        stiffness: v[4],
        damping_ratio: v[5],
    }
}

/// Small-signal flow admittance of the reed seen by the bore, dU/dp at the
/// static operating point of mouth pressure `p` (fs-duct's e^{-iwt}):
/// `Y_r = U0 (H0/(P_M h0)) H(w) - U0/(2p) + i w (S_r^2/k) H(w)` with
/// Bernoulli flow `U0 = w h0 sqrt(2 p/rho)`, static opening
/// `h0 = H0 (1 - p/P_M)`, effective area `S_r = k H0 / P_M`, and the reed's
/// dynamic gain `H = 1/(1 - (w/w_r)^2 - i q w/w_r)`, `q = 2 zeta`. The first
/// two terms are the modulated and the direct Bernoulli flow; the last is the
/// flow the moving reed sweeps.
fn reed_admittance(reed: &Reed, rho: f64, p: f64, omega: f64) -> C64 {
    let omega_r = (reed.stiffness / reed.mass).sqrt();
    let x = omega / omega_r;
    let gain = C64::new(1.0 - x * x, -2.0 * reed.damping_ratio * x).recip();
    let h0 = reed.rest_opening * (1.0 - p / reed.closing);
    let u0 = reed.width * h0 * (2.0 * p / rho).sqrt();
    let area = reed.stiffness * reed.rest_opening / reed.closing;
    let modulated = gain.scale(u0 * reed.rest_opening / (reed.closing * h0));
    let swept = gain * C64::new(0.0, -omega * area * area / reed.stiffness);
    modulated - C64::new(u0 / (2.0 * p), 0.0) - swept
}

/// One gas's bore input admittance, memoized by frequency: it does not depend
/// on the mouth pressure, so a threshold bisection reuses every evaluation.
struct BoreAdmittance<'a> {
    duct: &'a Duct,
    termination: Termination,
    gas: &'a GasState,
    cache: std::cell::RefCell<std::collections::HashMap<u64, C64>>,
}

impl BoreAdmittance<'_> {
    fn at(&self, hz: f64) -> C64 {
        *self
            .cache
            .borrow_mut()
            .entry(hz.to_bits())
            .or_insert_with(|| {
                let omega = core::f64::consts::TAU * hz;
                input_impedance(
                    self.duct,
                    self.gas,
                    omega,
                    LossModel::Bessel,
                    self.termination,
                )
                .unwrap()
                .impedance
                .recip()
            })
    }
}

/// Barkhausen point near one impedance maximum: the frequency in `window`
/// where the imaginary parts of the reed and bore admittances balance, and the
/// net conductance `Re(Y_r - Y_b)` there (positive grows).
fn barkhausen(
    reed: &Reed,
    bore: &BoreAdmittance<'_>,
    p: f64,
    window: (f64, f64),
    center: f64,
) -> Option<(f64, f64)> {
    let f = |hz: f64| {
        reed_admittance(reed, bore.gas.density, p, core::f64::consts::TAU * hz) - bore.at(hz)
    };
    let steps = 400;
    let mut best: Option<(f64, f64)> = None;
    let mut previous = (window.0, f(window.0));
    for i in 1..=steps {
        let hz = window.0 + (window.1 - window.0) * f64::from(i) / f64::from(steps);
        let value = f(hz);
        if previous.1.im.signum() != value.im.signum() {
            let (mut a, mut b) = (previous.0, hz);
            for _ in 0..60 {
                let mid = 0.5 * (a + b);
                if f(mid).im.signum() == f(a).im.signum() {
                    a = mid
                } else {
                    b = mid
                }
            }
            let root = 0.5 * (a + b);
            if best.is_none_or(|(r, _)| (root - center).abs() < (r - center).abs()) {
                best = Some((root, f(root).re));
            }
        }
        previous = (hz, value);
    }
    best
}

/// Linear analysis of one gas, independent of any render. For every
/// input-impedance maximum:
/// - its threshold: the lowest mouth pressure at which its Barkhausen point
///   turns unstable (bisected over 0.05..0.99 of the closing pressure), with
///   the frequency there;
/// - its Barkhausen frequency and small-signal growth rate at the held
///   pressure, `sigma = Re F Im(F') / |F'|^2` from the first-order root of
///   `F(w) = Y_r - Y_b` (e^{-iwt}: positive grows).
struct Analysis {
    maxima: Vec<(f64, f64)>,
    thresholds: Vec<Option<(f64, f64)>>,
    held: Vec<Option<(f64, f64)>>,
}

impl Analysis {
    /// The first maximum to destabilise as the pressure rises.
    fn lowest_threshold(&self) -> usize {
        (0..self.maxima.len())
            .filter(|&n| self.thresholds[n].is_some())
            .min_by(|&a, &b| {
                self.thresholds[a]
                    .unwrap()
                    .0
                    .total_cmp(&self.thresholds[b].unwrap().0)
            })
            .expect("an unstable impedance maximum")
    }

    /// The maximum that grows fastest at the held pressure.
    fn fastest_growth(&self) -> usize {
        (0..self.maxima.len())
            .filter(|&n| {
                self.thresholds[n].is_some() && self.held[n].is_some_and(|(_, sigma)| sigma > 0.0)
            })
            .max_by(|&a, &b| self.held[a].unwrap().1.total_cmp(&self.held[b].unwrap().1))
            .expect("a growing impedance maximum at the held pressure")
    }

    /// Band of maximum n: max(10 cents, twice its Barkhausen drift between
    /// threshold and held pressure).
    fn band_cents(&self, n: usize) -> f64 {
        let drift = 1200.0 * (self.held[n].unwrap().0 / self.thresholds[n].unwrap().1).log2();
        (2.0 * drift.abs()).max(10.0)
    }

    /// Index of the maximum nearest a played frequency.
    fn nearest(&self, f: f64) -> usize {
        (0..self.maxima.len())
            .min_by(|&a, &b| {
                (self.maxima[a].0 / f)
                    .ln()
                    .abs()
                    .total_cmp(&(self.maxima[b].0 / f).ln().abs())
            })
            .unwrap()
    }

    /// Whether `played` is maximum `n` within its band of the predicted
    /// held-pressure frequency.
    fn plays(&self, n: usize, played: f64) -> bool {
        self.nearest(played) == n
            && (1200.0 * (played / self.held[n].unwrap().0).log2()).abs() <= self.band_cents(n)
    }
}

/// Frequency window around impedance maximum n: 0.4 of the way to each neighbour.
fn mode_window(maxima: &[(f64, f64)], n: usize) -> (f64, f64) {
    let f_n = maxima[n].0;
    let below = if n == 0 { 50.0 } else { maxima[n - 1].0 };
    let above = maxima.get(n + 1).map_or(4000.0, |m| m.0);
    (f_n - 0.4 * (f_n - below), f_n + 0.4 * (above - f_n))
}

/// Barkhausen frequency and small-signal growth rate [1/s] of maximum n at
/// mouth pressure p: `sigma = Re F Im(F') / |F'|^2` from the first-order root
/// of `F(w) = Y_r - Y_b` (e^{-iwt}: positive grows).
fn growth_at(
    reed: &Reed,
    bore: &BoreAdmittance<'_>,
    maxima: &[(f64, f64)],
    n: usize,
    p: f64,
) -> Option<(f64, f64)> {
    barkhausen(reed, bore, p, mode_window(maxima, n), maxima[n].0).map(|(f, re)| {
        let w = core::f64::consts::TAU * f;
        let dw = 1e-6 * w;
        let at = |omega: f64| {
            reed_admittance(reed, bore.gas.density, p, omega)
                - bore.at(omega / core::f64::consts::TAU)
        };
        let slope = (at(w + dw) - at(w - dw)).scale(0.5 / dw);
        (f, re * slope.im / slope.norm_sq())
    })
}

fn analyse(
    reed: &Reed,
    duct: &Duct,
    termination: Termination,
    gas: &GasState,
    held: f64,
) -> Analysis {
    let maxima = impedance_maxima(duct, termination, gas);
    let bore = BoreAdmittance {
        duct,
        termination,
        gas,
        cache: std::cell::RefCell::new(std::collections::HashMap::new()),
    };
    let window_of = |n: usize| mode_window(&maxima, n);
    let mut thresholds = Vec::new();
    for (n, &(f_n, _)) in maxima.iter().enumerate() {
        let window = window_of(n);
        let grows = |p: f64| barkhausen(reed, &bore, p, window, f_n).is_some_and(|(_, g)| g > 0.0);
        let (mut lo, mut hi) = (0.05 * reed.closing, 0.99 * reed.closing);
        if !grows(hi) || grows(lo) {
            thresholds.push(None);
            continue;
        }
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if grows(mid) { hi = mid } else { lo = mid }
        }
        let p_th = 0.5 * (lo + hi);
        thresholds.push(Some((
            p_th,
            barkhausen(reed, &bore, p_th, window, f_n).unwrap().0,
        )));
    }
    let held_points = (0..maxima.len())
        .map(|n| growth_at(reed, &bore, &maxima, n, held))
        .collect();
    Analysis {
        maxima,
        thresholds,
        held: held_points,
    }
}

/// Render a performance through the source registry.
fn render(text: &str, registry: &ThermochemicalSources<'_>) -> (Vec<f64>, GasState, f64) {
    let block = 113;
    let mut performance = ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(),
        block,
        ObservationPolicy::Diagnostic,
        Some(registry),
    )
    .unwrap();
    let gas = *performance.ambient_gas();
    let rate = f64::from(performance.info().sample_rate_hz);
    let mut pressure = vec![0.0; performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) {
        performance.block(chunk).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    (pressure, gas, rate)
}

/// Heliox, by data alone, gated by a two-arm regime prediction (orchestrator
/// ruling 6187, declared on pexfp.362 before either arm was run). The clarinet
/// D3 performance (all holes closed) is rendered in CO2-free dry air and in
/// 80/20 He-O2 at the same temperature and pressure; within each arm only the
/// ambient record differs between the gases.
///
/// Independently of the render, a linear reed-bore analysis computes, for
/// every bore input-impedance maximum in each gas:
/// - its threshold;
/// - its Barkhausen frequency and growth rate at the held 4045 Pa.
///
/// It uses Bernoulli flow, the reed's dynamic gain and its swept flow (from the
/// reed row), against fs-duct's bore admittance.
///
/// The arms:
/// - (a) **Step start** (the committed attack: 0 -> 4045 Pa in 20 ms). The
///   fastest small-signal growth at 4045 Pa must win.
/// - (b) **Quasi-static ramp.** 3000 Pa (below every threshold) from 0.15 s,
///   then 3000 -> 4045 Pa over 2.0 s from 0.5 s, held to 4.3 s. The first
///   maximum to destabilise must establish and persist.
///
/// Band, declared before any prediction was computed: twice the maximum's
/// Barkhausen drift between threshold and held pressure, floor 10 cents.
///
/// - (b') **Heliox hold** (ruling 6205). 3500 Pa, between heliox's n=0 and
///   n=4 thresholds, held for 1.5 ln(1e4)/sigma_0, with sigma_0 the linear
///   growth of n=0 at 3500 Pa. Then 3500 -> 4045 Pa over 2 s, held 1.6 s. n=0
///   must establish during the hold and persist after the ramp. The fixture
///   declares compile_limits 400000 16384 for the longer schedule.
///
/// Gated, per ruling 6205: both gases on (a), and air on (b).
/// Printed as FAILED, ungated rows, open findings on pexfp.362:
/// - heliox (b): n=0 never establishes on the 2 s ramp; n=4 wins;
/// - heliox (b'), declared gated and FAILED on its first run: n=0 establishes
///   during the hold (268.29 Hz), then n=4 takes over after the ramp
///   (2076.5 Hz). The model therefore claims a heliox squeak at 4045 Pa
///   whatever the attack;
/// - the lowest-threshold criterion on the heliox step start: n=0 at 3471.8 Pa
///   against n=4 at 3533.1 Pa.
///
/// Also checked:
/// - Planted negative: the heliox step-start render must FAIL against air's
///   prediction.
/// - The cents below the nearest maximum are logged.
#[test]
fn heliox_regime_follows_the_attack_as_the_linear_reed_bore_analysis_predicts() {
    let sources = [
        source("nitrogen", Some("nitrogen-transport-lemmon2004")),
        source("oxygen", Some("oxygen-transport-lemmon2004")),
        source("argon", Some("argon-transport-lemmon2004")),
        source("helium", Some("helium-transport-hurly-mehl2007")),
    ];
    let bindings: Vec<SpeciesPackRef<'_>> = sources.iter().map(binding).collect();
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let original = std::fs::read_to_string(workspace_path(D3)).unwrap();
    let held = 4.045e3;
    let step_events = "events\t2\nevent\t1.5e-1\t2e-2\tpressure-pa\t4.045e3\nevent\t2.45e0\t6e-2\tpressure-pa\t0e0";
    assert_eq!(original.matches("\nambient 293.15 101325 0\n").count(), 1);
    assert_eq!(
        original.matches(step_events).count(),
        1,
        "the committed step-start schedule"
    );
    assert_eq!(original.matches("\naudio 48000 127200 20\n").count(), 1);
    let ramp_events = format!(
        "events\t3\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}",
        0.15, 0.02, 3000.0, 0.5, 2.0, held, 4.3, 0.06, 0.0
    );
    let ramp = original
        .replace("\naudio 48000 127200 20\n", "\naudio 48000 216000 20\n")
        .replace(step_events, &ramp_events);
    let mut step_results = Vec::new();
    for (label, ambient) in [("air", AIR), ("heliox", HELIOX)] {
        let with_gas =
            |text: &str| text.replace("\nambient 293.15 101325 0\n", &format!("\n{ambient}\n"));
        let step_text = with_gas(&original);
        let (step, gas, rate) = render(&step_text, &registry);
        let (ramped, _, _) = render(&with_gas(&ramp), &registry);
        let (duct, termination) = bore(&step_text);
        let analysis = analyse(&reed(&step_text), &duct, termination, &gas, held);
        let played_step = played_frequency(&step, rate, (1.2, 2.4));
        let mut hold_arm = None;
        if label == "heliox" {
            let first = analysis.lowest_threshold();
            let bore_admittance = BoreAdmittance {
                duct: &duct,
                termination,
                gas: &gas,
                cache: std::cell::RefCell::new(std::collections::HashMap::new()),
            };
            let hold_pa = 3500.0;
            let sigma_0 = growth_at(
                &reed(&step_text),
                &bore_admittance,
                &analysis.maxima,
                first,
                hold_pa,
            )
            .unwrap()
            .1;
            let sigma_fast = growth_at(
                &reed(&step_text),
                &bore_admittance,
                &analysis.maxima,
                analysis.fastest_growth(),
                hold_pa,
            )
            .map_or(f64::NAN, |g| g.1);
            assert!(sigma_0 > 0.0, "3500 Pa must lie above n=0's threshold");
            let hold_s = 1.5 * 1e4_f64.ln() / sigma_0;
            let ramp_at = 0.17 + hold_s;
            let total = ramp_at + 2.0 + 1.6;
            let samples = (total * 48_000.0).ceil() as u64;
            let events = format!(
                "events\t3\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}\nevent\t{:e}\t{:e}\tpressure-pa\t{:e}",
                0.15,
                0.02,
                hold_pa,
                ramp_at,
                2.0,
                held,
                total - 0.2,
                0.06,
                0.0
            );
            assert_eq!(
                original.matches("\ncompile_limits 100000 4096\n").count(),
                1
            );
            let hold_text = with_gas(&original)
                .replace(
                    "\naudio 48000 127200 20\n",
                    &format!("\naudio 48000 {samples} 20\n"),
                )
                .replace(
                    "\ncompile_limits 100000 4096\n",
                    "\ncompile_limits 400000 16384\n",
                )
                .replace(step_events, &events);
            let (held_render, _, _) = render(&hold_text, &registry);
            let played_hold = played_frequency(&held_render, rate, (ramp_at - 1.0, ramp_at));
            let played_after = played_frequency(&held_render, rate, (ramp_at + 2.6, ramp_at + 3.4));
            let established = analysis.nearest(played_hold) == first;
            let persists = analysis.plays(first, played_after);
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_hold_arm\",\"hold_pa\":{hold_pa},\"sigma_0_per_s\":{sigma_0},\"sigma_fastest_at_hold_per_s\":{sigma_fast},\"hold_s\":{hold_s},\"ramp_at_s\":{ramp_at},\"played_end_of_hold_Hz\":{played_hold},\"maximum_end_of_hold\":{},\"played_after_ramp_Hz\":{played_after},\"maximum_after_ramp\":{},\"predicted\":{first},\"established\":{established},\"persists\":{persists}}}",
                analysis.nearest(played_hold),
                analysis.nearest(played_after)
            );
            hold_arm = Some(established && persists);
        }
        let played_ramp = played_frequency(&ramped, rate, (3.3, 4.3));
        // Time course of the ramp (logged): RMS and the two strongest lines per
        // 0.25 s, with the held mouth pressure at each window's centre.
        for k in 0..15 {
            let (a, b) = (0.55 + 0.25 * f64::from(k), 0.8 + 0.25 * f64::from(k));
            let x = &ramped[(a * rate) as usize..(b * rate) as usize];
            let rms = (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt();
            let lines: Vec<String> = spectral_peaks(&ramped, rate, (a, b))
                .iter()
                .take(2)
                .map(|(f, m)| format!("[{f:.0},{m:.3}]"))
                .collect();
            let centre = 0.5 * (a + b);
            let mouth = if centre < 0.5 {
                3000.0
            } else {
                (3000.0 + (held - 3000.0) * (centre - 0.5) / 2.0).min(held)
            };
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_ramp_course\",\"gas\":\"{label}\",\"t_s\":[{a:.2},{b:.2}],\"mouth_pa\":{mouth:.0},\"rms_pa\":{rms:e},\"lines_Hz_rel\":[{}]}}",
                lines.join(",")
            );
        }
        let (fastest, first) = (analysis.fastest_growth(), analysis.lowest_threshold());
        let rows: Vec<String> = (0..analysis.maxima.len())
            .map(|n| {
                let t = analysis.thresholds[n]
                    .map_or_else(|| "null,null".into(), |(p, f)| format!("{p:.1},{f:.3}"));
                let h = analysis.held[n]
                    .map_or_else(|| "null,null".into(), |(f, s)| format!("{f:.3},{s:.3}"));
                format!("[{:.1},{t},{h}]", analysis.maxima[n].0)
            })
            .collect();
        let cents_below = |f: f64| 1200.0 * (f / analysis.maxima[analysis.nearest(f)].0).log2();
        println!(
            "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_two_arm\",\"gas\":\"{label}\",\"rho\":{},\"c\":{},\"maxima_Hz_threshold_Pa_Hz_held_Hz_growth_per_s\":[{}],\"step\":{{\"played_Hz\":{played_step},\"played_maximum\":{},\"cents_below_nearest_maximum\":{},\"predicted_fastest\":{fastest},\"predicted_Hz\":{},\"band_cents\":{},\"pass\":{}}},\"ramp\":{{\"played_Hz\":{played_ramp},\"played_maximum\":{},\"cents_below_nearest_maximum\":{},\"predicted_first\":{first},\"predicted_Hz\":{},\"band_cents\":{},\"pass\":{}}},\"lowest_threshold_on_step_ungated\":{{\"predicted\":{first},\"pass\":{}}}}}",
            gas.density,
            gas.sound_speed,
            rows.join(","),
            analysis.nearest(played_step),
            cents_below(played_step),
            analysis.held[fastest].unwrap().0,
            analysis.band_cents(fastest),
            analysis.plays(fastest, played_step),
            analysis.nearest(played_ramp),
            cents_below(played_ramp),
            analysis.held[first].unwrap().0,
            analysis.band_cents(first),
            analysis.plays(first, played_ramp),
            analysis.plays(first, played_step)
        );
        step_results.push((label, analysis, played_step, played_ramp, hold_arm));
    }
    let air = &step_results[0];
    let heliox = &step_results[1];
    println!(
        "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_ungated\",\"row\":\"heliox lowest threshold on the step start\",\"pass\":{},\"status\":\"FAILED before the two-arm gate, open finding on pexfp.362\"}}",
        heliox.1.plays(heliox.1.lowest_threshold(), heliox.2)
    );
    assert!(
        !air.1.plays(air.1.fastest_growth(), heliox.2),
        "planted negative: heliox's step start must not match air's prediction"
    );
    for (label, analysis, played_step, played_ramp, hold_arm) in &step_results {
        assert!(
            analysis.plays(analysis.fastest_growth(), *played_step),
            "{label} step start: played {played_step} Hz, not the fastest-growing maximum"
        );
        let ramp_pass = analysis.plays(analysis.lowest_threshold(), *played_ramp);
        if *label == "air" {
            assert!(
                ramp_pass,
                "air quasi-static ramp: played {played_ramp} Hz, not n=0"
            );
        } else {
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_ungated\",\"row\":\"heliox quasi-static ramp (b)\",\"pass\":{ramp_pass},\"status\":\"{}\"}}",
                if ramp_pass {
                    "passed"
                } else {
                    "FAILED, open finding on pexfp.362"
                }
            );
        }
        if let Some(hold_pass) = hold_arm {
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"heliox_regime_ungated\",\"row\":\"heliox hold arm (b')\",\"pass\":{hold_pass},\"status\":\"{}\"}}",
                if *hold_pass {
                    "passed"
                } else {
                    "FAILED, open finding on pexfp.362: n=0 established then n=4 took over"
                }
            );
        }
    }
}

/// Whether a held mouth pressure sustains oscillation. The attack leaves an
/// oscillation that either decays exponentially or grows to a limit cycle,
/// which can saturate before the early window. So it is called decaying when
/// the late RMS of the radiated pressure (2.0..2.4 s) falls below 0.95 of the
/// early RMS (0.6..1.0 s), a growth rate below ln(0.95)/1.4 s = -0.037/s.
/// The growth rate's measured slope near the sea-level threshold, about
/// 0.025 /s per Pa (-0.48/s at 3085 Pa, +1.42/s at 3162 Pa), puts that
/// ambiguous band within about 1.5 Pa.
fn held_pressure_sustains(original: &str, ambient: &str, pressure_pa: f64) -> (bool, f64, f64) {
    let event = "event\t1.5e-1\t2e-2\tpressure-pa\t4.045e3";
    assert_eq!(original.matches(event).count(), 1);
    let text = original
        .replace("\nambient 293.15 101325 0\n", &format!("\n{ambient}\n"))
        .replace(
            event,
            &format!("event\t1.5e-1\t2e-2\tpressure-pa\t{pressure_pa:e}"),
        );
    let block = 113;
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let rate = f64::from(performance.info().sample_rate_hz);
    let mut pressure = vec![0.0; performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) {
        performance.block(chunk).unwrap();
    }
    let rms = |a: f64, b: f64| {
        let x = &pressure[(a * rate) as usize..(b * rate) as usize];
        (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt()
    };
    let (early, late) = (rms(0.6, 1.0), rms(2.0, 2.4));
    (late >= 0.95 * early, early, late)
}

/// Quasi-static Bernoulli threshold of a beating reed on a bore: the reed's
/// negative flow conductance -dU/dp, with U = w H0 (1 - p/P_M) sqrt(2 p/rho),
/// equals the bore's input admittance 1/Z_peak at its first resonance. With
/// gamma = p/P_M and A = Z_peak w H0 sqrt(2/(rho P_M)) this is
/// (3 gamma - 1)/(2 sqrt gamma) = 1/A, so sqrt(gamma) = (2/A + sqrt(4/A^2 + 12))/6.
fn bernoulli_threshold_pa(text: &str, gas: &GasState) -> (f64, f64) {
    let reed: Vec<f64> = text
        .lines()
        .find_map(|line| line.strip_prefix("reed "))
        .expect("reed row")
        .split_ascii_whitespace()
        .map(|w| w.parse().unwrap())
        .collect();
    let (rest_opening, width, closing) = (reed[0], reed[1], reed[2]);
    let (duct, termination) = bore(text);
    let f1 = first_impedance_peak(&duct, termination, gas);
    let z_peak = input_impedance(
        &duct,
        gas,
        core::f64::consts::TAU * f1,
        LossModel::Bessel,
        termination,
    )
    .unwrap()
    .impedance
    .abs();
    let a = z_peak * width * rest_opening * (2.0 / (gas.density * closing)).sqrt();
    let s = (2.0 / a + (4.0 / (a * a) + 12.0).sqrt()) / 6.0;
    (s * s * closing, z_peak)
}

/// Altitude, by data alone: 0.7 atm at the same temperature (legacy dry-air
/// record) against 1 atm. The sound speed is unchanged, but density falls and
/// viscothermal loss rises, so Bernoulli flow and the bore's resonance height
/// both change. The oscillation threshold of the D3 performance is bisected in
/// each gas from held mouth pressures (eight steps, about 0.4 % resolution),
/// changing only the ambient record and the held pressure. The prediction is
/// the quasi-static Bernoulli threshold from the reed row and fs-duct's
/// resonance height in each gas. The band, 2 % on the threshold ratio, was
/// declared before measuring: reed inertia, swept flow and tone-hole details
/// shift the absolute threshold but do not depend on density, so they largely
/// cancel in the ratio.
#[test]
fn altitude_raises_the_oscillation_threshold_as_bernoulli_predicts() {
    let original = std::fs::read_to_string(workspace_path(D3)).unwrap();
    let mut results = Vec::new();
    for (label, ambient) in [
        ("sea-level", "ambient 293.15 101325 0"),
        ("0.7-atm", "ambient 293.15 70927.5 0"),
    ] {
        let text = original.replace("\nambient 293.15 101325 0\n", &format!("\n{ambient}\n"));
        let gas = *ReedPerformance::from_bytes(text.as_bytes(), 113)
            .unwrap()
            .ambient_gas();
        let (predicted, z_peak) = bernoulli_threshold_pa(&text, &gas);
        let closing = 7076.1822222222;
        let (mut lo, mut hi) = (0.25 * closing, 0.6 * closing);
        for (bound, expected) in [(lo, false), (hi, true)] {
            let (sustains, early, late) = held_pressure_sustains(&original, ambient, bound);
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"altitude_threshold\",\"gas\":\"{label}\",\"held_pa\":{bound},\"early_rms_pa\":{early:e},\"late_rms_pa\":{late:e},\"sustains\":{sustains}}}"
            );
            assert_eq!(sustains, expected, "{label}: bracket end {bound} Pa");
        }
        for _ in 0..8 {
            let mid = 0.5 * (lo + hi);
            let (sustains, early, late) = held_pressure_sustains(&original, ambient, mid);
            println!(
                "{{\"suite\":\"gas-species-packs\",\"test\":\"altitude_threshold\",\"gas\":\"{label}\",\"held_pa\":{mid},\"early_rms_pa\":{early:e},\"late_rms_pa\":{late:e},\"sustains\":{sustains}}}"
            );
            if sustains { hi = mid } else { lo = mid }
        }
        let measured = 0.5 * (lo + hi);
        println!(
            "{{\"suite\":\"gas-species-packs\",\"test\":\"altitude_threshold\",\"gas\":\"{label}\",\"rho\":{},\"c\":{},\"z_peak\":{z_peak:e},\"bernoulli_threshold_pa\":{predicted},\"measured_threshold_pa\":{measured},\"bracket_pa\":{}}}",
            gas.density,
            gas.sound_speed,
            hi - lo
        );
        results.push((predicted, measured));
    }
    let predicted = results[1].0 / results[0].0;
    let measured = results[1].1 / results[0].1;
    let deviation = measured / predicted - 1.0;
    println!(
        "{{\"suite\":\"gas-species-packs\",\"test\":\"altitude_threshold\",\"predicted_ratio\":{predicted},\"measured_ratio\":{measured},\"deviation\":{deviation:e},\"band\":0.02}}"
    );
    assert!(
        measured > 1.0,
        "thinner air must raise the threshold: ratio {measured}"
    );
    assert!(
        deviation.abs() <= 0.02,
        "threshold ratio {measured} vs Bernoulli {predicted}"
    );
}
