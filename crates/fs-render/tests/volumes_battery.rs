//! fs-render volumes battery (bead qfx.3): Woodcock unbiasedness
//! against closed forms and majorant looseness, phase-function moment
//! gates, spectral emission vs the analytic slab, the ZERO-COPY live
//! LBM binding, and bitwise replay. Statistical gates are 3σ with the
//! empirical σ printed — seeds are fixed, so failures are real.

use fs_lbm::freesurface::{ContactModel, dam_break};
use fs_render::volumes::{
    MajorantGrid, Ray, VolumeGrid, beer_lambert, hg_sample_cos, pixel_stream, planck,
    rayleigh_sample_cos, render_transmittance, woodcock_emission, woodcock_transmittance,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn stream(tile: u32) -> fs_rand::Stream {
    pixel_stream(0x1001_2026_0708_0300, tile)
}

/// vol-001: homogeneous slab — the Woodcock mean matches
/// exp(−σL) within 3σ_stat for several optical depths, and the
/// Beer–Lambert fast path is exact.
#[test]
fn vol_001_homogeneous_slab() {
    let n = 40_000u32;
    for (case, sigma, l) in [(0u32, 0.5f64, 1.0f64), (1, 1.5, 1.0), (2, 3.0, 0.7)] {
        let data = vec![sigma; 8];
        let grid = VolumeGrid::new([2, 2, 2], &data, [0.0; 3], [0.5, 0.5, l / 2.0]);
        let mut s = stream(100 + case);
        let mut acc = 0.0;
        let mut acc2 = 0.0;
        let m_at = |_p: [f64; 3]| sigma;
        for _ in 0..n {
            let ray = Ray {
                origin: [0.25, 0.25, -1.0],
                dir: [0.0, 0.0, 1.0],
                t0: 1.0,
                t1: 1.0 + l,
            };
            let (tr, _) = woodcock_transmittance(&grid, &m_at, sigma, ray, &mut s);
            acc += tr;
            acc2 += tr * tr;
        }
        let mean = acc / f64::from(n);
        let var = (acc2 / f64::from(n) - mean * mean).max(0.0);
        let se = (var / f64::from(n)).sqrt();
        let want = beer_lambert(sigma, l);
        verdict(
            &format!("vol-001-slab-{case}"),
            (mean - want).abs() < 3.0 * se.max(1e-6),
            &format!(
                "sigma*L={:.2}: mean {mean:.4} vs exp {want:.4} (3se {:.4})",
                sigma * l,
                3.0 * se
            ),
        );
    }
}

/// vol-002: unbiasedness under majorant LOOSENESS — a heterogeneous
/// blob renders the same mean transmittance with a tight and a 3×
/// loose bound (only the null-collision ledger moves), and both match
/// a fine deterministic quadrature reference.
#[test]
fn vol_002_majorant_unbiasedness() {
    let dims = [16usize, 16, 16];
    let mut data = vec![0.0f64; 16 * 16 * 16];
    for k in 0..16 {
        for j in 0..16 {
            for i in 0..16 {
                let d = [(i, 8.0f64), (j, 8.0), (k, 8.0)]
                    .iter()
                    .map(|&(a, c)| (a as f64 + 0.5 - c).powi(2))
                    .sum::<f64>();
                data[i + 16 * (j + 16 * k)] = 2.5 * (-d / 18.0).exp();
            }
        }
    }
    let grid = VolumeGrid::new(dims, &data, [0.0; 3], [1.0 / 16.0; 3]);
    let majorant = MajorantGrid::build(&grid, 4);
    let origin = [0.53f64, 0.47, -0.5];
    let dir = [0.0f64, 0.0, 1.0];
    let (t0, t1) = (0.5, 1.5);
    // Deterministic reference: fine midpoint quadrature of ∫σ.
    let steps = 40_000usize;
    let mut tau = 0.0;
    for q in 0..steps {
        let t = t0 + (t1 - t0) * ((q as f64 + 0.5) / steps as f64);
        let p = [origin[0], origin[1], dir[2].mul_add(t, origin[2])];
        tau += grid.sigma_at(p) * (t1 - t0) / steps as f64;
    }
    let reference = (-tau).exp();
    let bound = grid.global_majorant();
    let n = 60_000u32;
    let run = |factor: f64, tile: u32| -> (f64, f64, u64) {
        let mut s = stream(tile);
        let mut acc = 0.0;
        let mut acc2 = 0.0;
        let mut nulls = 0u64;
        let m_at = |p: [f64; 3]| majorant.at(&grid, p);
        for _ in 0..n {
            let (tr, nu) = woodcock_transmittance(
                &grid,
                &m_at,
                bound * factor,
                Ray {
                    origin,
                    dir,
                    t0,
                    t1,
                },
                &mut s,
            );
            acc += tr;
            acc2 += tr * tr;
            nulls += u64::from(nu);
        }
        let mean = acc / f64::from(n);
        let var = (acc2 / f64::from(n) - mean * mean).max(0.0);
        (mean, (var / f64::from(n)).sqrt(), nulls)
    };
    let (m_tight, se_t, nulls_tight) = run(1.0, 200);
    let (m_loose, se_l, nulls_loose) = run(3.0, 201);
    verdict(
        "vol-002-unbiased-tight",
        (m_tight - reference).abs() < 3.0 * se_t.max(1e-6),
        &format!(
            "tight bound: {m_tight:.4} vs quadrature {reference:.4} (3se {:.4})",
            3.0 * se_t
        ),
    );
    verdict(
        "vol-002-unbiased-loose",
        (m_loose - reference).abs() < 3.0 * se_l.max(1e-6),
        &format!(
            "3x-loose bound: {m_loose:.4} vs quadrature {reference:.4} (3se {:.4})",
            3.0 * se_l
        ),
    );
    verdict(
        "vol-002-null-ledger",
        nulls_loose > 2 * nulls_tight,
        &format!(
            "LEDGER null collisions: tight {nulls_tight} vs loose {nulls_loose} (cost of looseness, not bias)"
        ),
    );
}

/// vol-003: phase functions — HG mean cosine equals g, Rayleigh has
/// zero mean and E[cos²] = 2/5 (exact Cardano inversion under test).
#[test]
fn vol_003_phase_functions() {
    let n = 80_000u32;
    let g = 0.6;
    let mut s = stream(300);
    let (mut m1, mut m2) = (0.0f64, 0.0f64);
    for _ in 0..n {
        let c = hg_sample_cos(g, s.next_f64());
        m1 += c;
        m2 += c * c;
    }
    let mean = m1 / f64::from(n);
    let var = (m2 / f64::from(n) - mean * mean).max(0.0);
    let se = (var / f64::from(n)).sqrt();
    verdict(
        "vol-003-hg-mean-cosine",
        (mean - g).abs() < 3.0 * se,
        &format!("HG g=0.6: E[cos] {mean:.4} (3se {:.4})", 3.0 * se),
    );
    let mut s = stream(301);
    let (mut r1, mut r2, mut r4) = (0.0f64, 0.0f64, 0.0f64);
    for _ in 0..n {
        let c = rayleigh_sample_cos(s.next_f64());
        r1 += c;
        r2 += c * c;
        r4 += c * c * c * c;
    }
    let mean1 = r1 / f64::from(n);
    let mean2 = r2 / f64::from(n);
    let var2 = (r4 / f64::from(n) - mean2 * mean2).max(0.0);
    let se2 = (var2 / f64::from(n)).sqrt();
    verdict(
        "vol-003-rayleigh-moments",
        mean1.abs() < 3.0 * (mean2 / f64::from(n)).sqrt().max(1e-4)
            && (mean2 - 0.4).abs() < 3.0 * se2,
        &format!(
            "Rayleigh: E[cos] {mean1:.4}, E[cos2] {mean2:.4} vs 0.4 (3se {:.4})",
            3.0 * se2
        ),
    );
}

/// vol-004: emissive hot slab — the collision estimator matches
/// B_λ(T)·(1 − exp(−σL)) per hero wavelength, and Planck weights are
/// physically ordered (hotter is brighter; hotter is bluer).
#[test]
fn vol_004_spectral_emission() {
    let sigma = 1.2f64;
    let l = 1.0f64;
    let t_k = 4500.0;
    let data = vec![sigma; 8];
    let grid = VolumeGrid::new([2, 2, 2], &data, [0.0; 3], [0.5, 0.5, 0.5]);
    let n = 60_000u32;
    for (case, lambda) in [(0u32, 480.0f64), (1, 600.0), (2, 720.0)] {
        let b = planck(lambda, t_k);
        let src = move |_p: [f64; 3]| b;
        let mut s = stream(400 + case);
        let mut acc = 0.0;
        for _ in 0..n {
            let ray = Ray {
                origin: [0.25, 0.25, -1.0],
                dir: [0.0, 0.0, 1.0],
                t0: 1.0,
                t1: 1.0 + l,
            };
            acc += woodcock_emission(&grid, sigma, &src, ray, &mut s);
        }
        let mean = acc / f64::from(n);
        let want = b * (1.0 - beer_lambert(sigma, l));
        let p_hit = 1.0 - beer_lambert(sigma, l);
        let rel_se = ((1.0 - p_hit) / (p_hit * f64::from(n))).sqrt();
        verdict(
            &format!("vol-004-emission-{case}"),
            (mean / want - 1.0).abs() < 3.0 * rel_se,
            &format!(
                "lambda {lambda} nm: mean/analytic = {:.4} (3 rel-se {:.4})",
                mean / want,
                3.0 * rel_se
            ),
        );
    }
    let hotter_brighter = planck(600.0, 5000.0) > planck(600.0, 4000.0);
    let hotter_bluer = planck(480.0, 6000.0) / planck(720.0, 6000.0)
        > planck(480.0, 3000.0) / planck(720.0, 3000.0);
    verdict(
        "vol-004-planck-ordering",
        hotter_brighter && hotter_bluer,
        "hotter is brighter at fixed lambda; hotter is bluer across lambdas",
    );
}

/// vol-005: the LIVE LBM binding — a dam-break mass field renders
/// through a ZERO-COPY borrow (the API takes `&[f64]`; the buffer is
/// the simulation's own), the image is bitwise replayable, and the
/// free surface is visible as a transmittance gradient (gas above
/// bright, fluid below dark).
#[test]
fn vol_005_lbm_binding() {
    let mut sim = dam_break(40, 24, 10, 1e-4, 0.0, ContactModel::Neutral);
    for _ in 0..120 {
        sim.step();
    }
    // Bind the simulation's OWN mass buffer — no copy, no transform.
    let grid = VolumeGrid::new([40, 24, 1], &sim.mass, [0.0, 0.0, 0.0], [1.0, 1.0, 24.0]);
    let majorant = MajorantGrid::build(&grid, 8);
    let img1 = render_transmittance(&grid, &majorant, 24, 24, 0xBEEF);
    let img2 = render_transmittance(&grid, &majorant, 24, 24, 0xBEEF);
    let bitwise = img1
        .iter()
        .zip(&img2)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    // Compare INTERIOR regions over the dam's footprint (the tank's
    // outermost pixel rows/columns sample the mass-0 WALL cells and
    // are transparent by construction — the first gate draft included
    // them and measured 0.875 where the fluid reads ~0): left-side
    // columns, low rows = fluid (opaque), same columns high rows =
    // gas (transparent).
    let region_mean = |img: &[f64], rows: std::ops::Range<usize>| -> f64 {
        let mut acc = 0.0;
        let mut cnt = 0.0;
        for py in rows {
            for px in 2..8 {
                acc += img[py * 24 + px];
                cnt += 1.0;
            }
        }
        acc / cnt
    };
    let bottom = region_mean(&img1, 1..3);
    let top = region_mean(&img1, 20..22);
    let spread = img1
        .iter()
        .fold((1.0f64, 0.0f64), |(lo, hi), &v| (lo.min(v), hi.max(v)));
    verdict(
        "vol-005-lbm-zero-copy-render",
        bitwise && top > bottom + 0.3 && spread.1 > spread.0,
        &format!(
            "bitwise replay; free surface visible: top transmittance {top:.3} vs bottom {bottom:.3}; range [{:.3}, {:.3}]",
            spread.0, spread.1
        ),
    );
}

/// vol-006: per-pixel stream isolation — any pixel recomputed
/// standalone equals its value in the full render (tile order CANNOT
/// matter because no state crosses pixels).
#[test]
fn vol_006_pixel_isolation() {
    let data: Vec<f64> = (0..64).map(|i| 0.4 + 0.05 * f64::from(i % 8)).collect();
    let grid = VolumeGrid::new([4, 4, 4], &data, [0.0; 3], [0.25; 3]);
    let majorant = MajorantGrid::build(&grid, 2);
    let img = render_transmittance(&grid, &majorant, 8, 64, 0xF00D);
    let (lo, hi) = grid.bounds();
    let bound = grid.global_majorant();
    let mut worst = 0.0f64;
    for &(px, py) in &[(0usize, 0usize), (3, 5), (7, 7), (5, 2)] {
        let pixel = u32::try_from(py * 8 + px).expect("small");
        let mut s = pixel_stream(0xF00D, pixel);
        let x = lo[0] + (hi[0] - lo[0]) * ((px as f64 + 0.5) / 8.0);
        let y = lo[1] + (hi[1] - lo[1]) * ((py as f64 + 0.5) / 8.0);
        let mut acc = 0.0;
        let m_at = |p: [f64; 3]| majorant.at(&grid, p);
        for _ in 0..64 {
            let ray = Ray {
                origin: [x, y, hi[2] + 1.0],
                dir: [0.0, 0.0, -1.0],
                t0: 1.0,
                t1: 1.0 + (hi[2] - lo[2]),
            };
            let (tr, _) = woodcock_transmittance(&grid, &m_at, bound, ray, &mut s);
            acc += tr;
        }
        worst = worst.max((acc / 64.0 - img[py * 8 + px]).abs());
    }
    verdict(
        "vol-006-pixel-isolation",
        worst == 0.0,
        &format!("standalone pixel recompute deviation {worst:.2e} (bitwise zero required)"),
    );
}
