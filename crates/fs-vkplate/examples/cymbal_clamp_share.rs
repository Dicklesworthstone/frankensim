//! Felt-clamp share of each stand-mounted cymbal mode (.372, rulings 6411 and 6423).
//!
//! usage: cymbal_clamp_share [splash] [crash] [ride]
//!
//! For each cymbal, K is the linear stiffness about the rest shape (bending, felt support and
//! the von Karman curvature Hessian), assembled column by column from the body's own operator
//! (`Body::rest_stiffness_apply`), with the lumped mass M. After an asymmetry gate at rounding
//! level (1e-14; the linear part is exactly symmetric and the curvature Hessian's central
//! difference at 1e-12 m sits at its ~2e-15 roundoff floor), the modes
//! of M^-1/2 K M^-1/2 (fs-la `householder_eigh`) give, for mass-normalised phi_a,
//! f_support,a = ks h^2 sum_annulus phi_ia^2 / omega_a^2, the share of the mode's strain energy
//! in the felt springs. The clamp loss is eta_clamp,a = omega_a (R/K)_felt f_support,a, with
//! (R/K) = gamma tau0 / (1 - gamma) from the piano hammer-felt regressions (Stulov's
//! slow-loading form) over keys 1..88: flags composition-proxy and strain-regime-proxy. It is
//! compared with the bell-bronze material loss (van Houten et al. 1997, four bells in vacuum):
//! the inverse-variance mean and the [min, max] bracket, flags composition-proxy and
//! extrapolate-frequency. A mode is clamp-dominated if its clamp bracket lies above the bronze
//! bracket, bronze-dominated if below, undecided otherwise; modes with fewer than 6 FD nodes per
//! bending wavelength are listed apart (their shares are grid numbers).
use fs_la::eigen::householder_eigh;
use fs_vkplate::kit::{self, BRONZE, CymbalSpec};
use std::time::Instant;

const AUDIBLE_HZ: f64 = 20_000.0;

/// Scalar rows of a material-tsv-v1 seed, by property name: (value, unit).
fn seed_scalars(slug: &str) -> Vec<(String, f64, String)> {
    let path = format!(
        "{}/../../data/matdb/seed-v1/{slug}/properties.tsv",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines()
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f.first() == Some(&"scalar")).then(|| {
                (
                    f[3].to_owned(),
                    f[4].parse::<f64>().expect("scalar value"),
                    f[5].to_owned(),
                )
            })
        })
        .collect()
}

fn scalar(rows: &[(String, f64, String)], name: &str) -> (f64, String) {
    let (_, v, u) = rows
        .iter()
        .find(|(n, _, _)| n == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    (*v, u.clone())
}

/// Bronze loss: the inverse-variance mean of bells A-D and the [min, max] of the four claims.
fn bronze_loss() -> (f64, f64, f64) {
    let rows = seed_scalars("bronze-bell-eijsbouts-vacuum-loss-vanhouten1997");
    let (mut num, mut den, mut lo, mut hi) = (0.0, 0.0, f64::INFINITY, 0.0_f64);
    for bell in ["a", "b", "c", "d"] {
        let (eta, _) = scalar(&rows, &format!("material_loss_factor_bell_{bell}"));
        let (sd, _) = scalar(
            &rows,
            &format!("material_loss_factor_bell_{bell}_printed_standard_deviation"),
        );
        num += eta / (sd * sd);
        den += 1.0 / (sd * sd);
        lo = lo.min(eta);
        hi = hi.max(eta);
    }
    (num / den, lo, hi)
}

/// (R/K) = gamma tau0 / (1 - gamma) [s] of the hammer-felt regressions at keys 1..88.
fn felt_relaxation_times() -> Vec<f64> {
    let rows = seed_scalars("hammer-felt-abel-stulov2005");
    let micro = |name: &str| {
        let (v, u) = scalar(&rows, name);
        assert!(u == "\u{b5}s" || u == "\u{3bc}s", "{name} unit {u}");
        v * 1e-6
    };
    let (g0, _) = scalar(&rows, "hammer_hereditary_amplitude_regression_constant");
    let (g1, _) = scalar(
        &rows,
        "hammer_hereditary_amplitude_regression_linear_coefficient",
    );
    let t0 = micro("hammer_relaxation_time_regression_constant");
    let t1 = micro("hammer_relaxation_time_regression_linear_coefficient_magnitude_negative");
    let t2 = micro("hammer_relaxation_time_regression_quadratic_coefficient");
    (1..=88)
        .map(|n| {
            let n = f64::from(n);
            let gamma = g0 + g1 * n;
            let tau = t0 - t1 * n + t2 * n * n;
            gamma * tau / (1.0 - gamma)
        })
        .collect()
}

/// Eigenpairs of the symmetric n x n `a` through the four D2 parity sectors of the centred square
/// grid. The cymbal's thickness, rest shape and felt annulus depend on r only and its stencils are
/// reflection-symmetric, so `a` commutes with x -> -x and y -> -y and the sectors decouple; each is
/// small enough for fs-la's dense householder_eigh. Returns eigenvalues ascending and eigenvectors
/// as columns of a row-major n x n matrix (householder_eigh's layout). The largest cross-sector
/// coupling is printed; where the full matrix is admissible the sector eigenvalues are checked
/// against the full solve.
fn sector_eigh(
    a: &[f64],
    n: usize,
    xy: &[[f64; 2]],
    h: f64,
    name: &str,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let cell = |v: f64| (v / h).round() as i64;
    let mut at = std::collections::BTreeMap::new();
    for (i, p) in xy.iter().enumerate() {
        if (p[0] - cell(p[0]) as f64 * h).abs() > 1e-9 * h
            || (p[1] - cell(p[1]) as f64 * h).abs() > 1e-9 * h
        {
            return Err(format!("{name}: node {i} is off the centred grid"));
        }
        at.insert((cell(p[0]), cell(p[1])), i);
    }
    // Orthonormal symmetry-adapted vectors, sparse: (node, coefficient) lists per sector.
    let mut sectors: [Vec<Vec<(usize, f64)>>; 4] = Default::default();
    for (&(i, j), _) in at.iter().filter(|((i, j), _)| *i >= 0 && *j >= 0) {
        for (s, (sx, sy)) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]
            .into_iter()
            .enumerate()
        {
            let mut coef = std::collections::BTreeMap::<usize, f64>::new();
            for (fx, fy) in [(1i64, 1i64), (-1, 1), (1, -1), (-1, -1)] {
                let node = *at.get(&(fx * i, fy * j)).ok_or_else(|| {
                    format!("{name}: grid not reflection-symmetric at ({i}, {j})")
                })?;
                let c = (if fx < 0 { sx } else { 1.0 }) * (if fy < 0 { sy } else { 1.0 });
                *coef.entry(node).or_default() += c;
            }
            let norm = coef.values().map(|c| c * c).sum::<f64>().sqrt();
            if norm > 1e-12 {
                sectors[s].push(
                    coef.into_iter()
                        .filter(|(_, c)| c.abs() > 0.0)
                        .map(|(k, c)| (k, c / norm))
                        .collect(),
                );
            }
        }
    }
    let total: usize = sectors.iter().map(Vec::len).sum();
    if total != n {
        return Err(format!(
            "{name}: symmetry-adapted basis has {total} vectors for {n} nodes"
        ));
    }
    let project = |u: &[(usize, f64)], v: &[(usize, f64)]| -> f64 {
        u.iter()
            .map(|&(i, ci)| v.iter().map(|&(j, cj)| ci * cj * a[i * n + j]).sum::<f64>())
            .sum()
    };
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    // Cross-sector coupling over a deterministic sample of basis pairs.
    let mut coupling = 0.0_f64;
    for s in 0..4 {
        for t in (s + 1)..4 {
            for p in (0..sectors[s].len()).step_by(7) {
                for q in (0..sectors[t].len()).step_by(7) {
                    coupling = coupling.max(project(&sectors[s][p], &sectors[t][q]).abs() / scale);
                }
            }
        }
    }
    let mut pairs: Vec<(f64, Vec<f64>)> = Vec::with_capacity(n);
    for basis in &sectors {
        let m = basis.len();
        if m == 0 {
            continue;
        }
        let mut block = vec![0.0; m * m];
        for p in 0..m {
            for q in p..m {
                let v = project(&basis[p], &basis[q]);
                block[p * m + q] = v;
                block[q * m + p] = v;
            }
        }
        let (lam, y) = householder_eigh(&block, m);
        for (c, &l) in lam.iter().enumerate() {
            let mut full = vec![0.0; n];
            for p in 0..m {
                let yp = y[p * m + c];
                for &(i, coef) in &basis[p] {
                    full[i] += coef * yp;
                }
            }
            pairs.push((l, full));
        }
    }
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
    let lambda: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let mut vectors = vec![0.0; n * n];
    for (c, (_, v)) in pairs.iter().enumerate() {
        for i in 0..n {
            vectors[i * n + c] = v[i];
        }
    }
    let full_check = if fs_la::eigen::admit_jacobi_eigh(n).is_ok() {
        let (full, _) = householder_eigh(a, n);
        let worst = full
            .iter()
            .zip(&lambda)
            .map(|(f, s)| (f - s).abs() / f.abs().max(f64::MIN_POSITIVE))
            .fold(0.0_f64, f64::max);
        format!("{worst:e}")
    } else {
        "\"inadmissible\"".into()
    };
    println!(
        "{{\"cymbal\":\"{name}\",\"sectors\":[{},{},{},{}],\"cross_sector_coupling_rel\":{coupling:e},\"sector_vs_full_max_rel_eigenvalue_diff\":{full_check}}}",
        sectors[0].len(),
        sectors[1].len(),
        sectors[2].len(),
        sectors[3].len()
    );
    Ok((lambda, vectors))
}

fn analyse(spec: &CymbalSpec, bronze: (f64, f64, f64), tau: (f64, f64, f64)) -> Result<(), String> {
    let started = Instant::now();
    let mut piece = kit::cymbal(spec, &[(0.005, 1.0, 0.5, 0.005)])?;
    let bi = piece
        .instrument
        .bodies
        .iter()
        .position(|b| b.name == spec.name)
        .ok_or("cymbal body")?;
    let body = &mut piece.instrument.bodies[bi];
    let (n, h) = (body.n, body.h);
    let mass = body.mass.clone();
    let (r0, r1, ks, _cs) = kit::felt_support(spec);
    let annulus: Vec<usize> = (0..n)
        .filter(|&i| {
            let r = (body.xy[i][0].powi(2) + body.xy[i][1].powi(2)).sqrt();
            r >= r0 && r <= r1
        })
        .collect();
    // K column by column from the body's own operator.
    let mut k = vec![0.0; n * n];
    let mut unit = vec![0.0; n];
    for j in 0..n {
        unit[j] = 1.0;
        let col = body.rest_stiffness_apply(&unit)?;
        unit[j] = 0.0;
        for i in 0..n {
            k[i * n + j] = col[i];
        }
    }
    let assembled_s = started.elapsed().as_secs_f64();
    // Diagnose: the linear part alone, and the curvature Hessian B = K - K_lin.
    let mut k_lin = vec![0.0; n * n];
    for j in 0..n {
        unit[j] = 1.0;
        let col = body.linear_stiffness_apply(&unit)?;
        unit[j] = 0.0;
        for i in 0..n {
            k_lin[i * n + j] = col[i];
        }
    }
    let asym = |a: &dyn Fn(usize, usize) -> f64| -> (f64, f64, usize, usize) {
        let (mut d2, mut n2, mut worst) = (0.0_f64, 0.0_f64, (0.0_f64, 0usize, 0usize));
        for i in 0..n {
            for j in 0..n {
                let d = a(i, j) - a(j, i);
                d2 += d * d;
                n2 += a(i, j) * a(i, j);
                if d.abs() > worst.0 {
                    worst = (d.abs(), i, j);
                }
            }
        }
        ((d2 / n2).sqrt(), worst.0, worst.1, worst.2)
    };
    let lin = asym(&|i, j| k_lin[i * n + j]);
    let hess = asym(&|i, j| k[i * n + j] - k_lin[i * n + j]);
    let r_of = |i: usize| (body.xy[i][0].powi(2) + body.xy[i][1].powi(2)).sqrt();
    println!(
        "{{\"cymbal\":\"{}\",\"diagnosis\":true,\"linear_asymmetry\":{:e},\"linear_worst_abs\":{:e},\"linear_worst_radii_m\":[{:.4},{:.4}],\"hessian_asymmetry\":{:e},\"hessian_worst_abs\":{:e},\"hessian_worst_radii_m\":[{:.4},{:.4}],\"bell_radius_m\":{},\"bell_blend_end_m\":{}}}",
        spec.name,
        lin.0,
        lin.1,
        r_of(lin.2),
        r_of(lin.3),
        hess.0,
        hess.1,
        r_of(hess.2),
        r_of(hess.3),
        spec.bell_radius,
        spec.bell_radius + 0.02
    );
    drop(k_lin);
    let (mut diff2, mut norm2, mut worst) = (0.0_f64, 0.0_f64, (0.0_f64, 0usize, 0usize));
    for i in 0..n {
        for j in 0..n {
            let d = k[i * n + j] - k[j * n + i];
            diff2 += d * d;
            norm2 += k[i * n + j] * k[i * n + j];
            if d.abs() > worst.0 {
                worst = (d.abs(), i, j);
            }
        }
    }
    let asymmetry = (diff2 / norm2).sqrt();
    let radius = |i: usize| (body.xy[i][0].powi(2) + body.xy[i][1].powi(2)).sqrt();
    println!(
        "{{\"cymbal\":\"{}\",\"nodes\":{n},\"grid_m\":{h},\"annulus_nodes\":{},\"ks_n_m3\":{ks:e},\"asymmetry_rel_frobenius\":{asymmetry:e},\"worst_pair\":[{},{}],\"worst_pair_radii_m\":[{:.4},{:.4}],\"worst_pair_abs\":{:e},\"assembly_s\":{assembled_s:.1}}}",
        spec.name,
        annulus.len(),
        worst.1,
        worst.2,
        radius(worst.1),
        radius(worst.2),
        worst.0
    );
    if asymmetry > 1e-14 {
        return Err(format!(
            "{}: K is not symmetric to rounding (relative {asymmetry:e}; worst pair nodes {} and {} at r {:.4} and {:.4} m): stop and report",
            spec.name,
            worst.1,
            worst.2,
            radius(worst.1),
            radius(worst.2)
        ));
    }
    // The symmetric form; the gate above bounds what the symmetrisation removes.
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = 0.5 * (k[i * n + j] + k[j * n + i]) / (mass[i] * mass[j]).sqrt();
        }
    }
    drop(k);
    let (lambda, vectors) = sector_eigh(&a, n, &body.xy, h, spec.name)?;
    let solved_s = started.elapsed().as_secs_f64();
    let (e, nu, rho) = BRONZE;
    let t_rim = spec.bow_thickness;
    let d_rim = e * t_rim.powi(3) / (12.0 * (1.0 - nu * nu));
    let (mut counts, mut coarse) = ([0usize; 3], 0usize);
    let mut dominated = Vec::new();
    for (mode, &lam) in lambda.iter().enumerate() {
        if lam <= 0.0 {
            return Err(format!(
                "{}: mode {mode} has a non-positive eigenvalue {lam:e}",
                spec.name
            ));
        }
        let omega = lam.sqrt();
        let f_hz = omega / (2.0 * std::f64::consts::PI);
        if f_hz > AUDIBLE_HZ {
            continue;
        }
        let share: f64 = annulus
            .iter()
            .map(|&i| {
                let phi = vectors[i * n + mode] / mass[i].sqrt();
                ks * h * h * phi * phi
            })
            .sum::<f64>()
            / lam;
        let k_b = (rho * t_rim * omega * omega / d_rim).powf(0.25);
        let nodes_per_wavelength = 2.0 * std::f64::consts::PI / k_b / h;
        let clamp = [
            omega * tau.0 * share,
            omega * tau.1 * share,
            omega * tau.2 * share,
        ];
        let verdict = if clamp[0] > bronze.2 {
            0
        } else if clamp[2] < bronze.1 {
            1
        } else {
            2
        };
        let label = ["clamp-dominated", "bronze-dominated", "undecided"][verdict];
        let resolved = nodes_per_wavelength >= 6.0;
        if resolved {
            counts[verdict] += 1;
            if verdict != 1 {
                dominated.push(format!("{f_hz:.1}"));
            }
        } else {
            coarse += 1;
        }
        println!(
            "{{\"cymbal\":\"{}\",\"mode\":{mode},\"f_hz\":{f_hz:.2},\"f_support\":{share:e},\"eta_clamp_min_median_max\":[{:e},{:e},{:e}],\"eta_bronze_mean_min_max\":[{:e},{:e},{:e}],\"nodes_per_bending_wavelength\":{nodes_per_wavelength:.2},\"resolved\":{resolved},\"verdict\":\"{label}\"}}",
            spec.name, clamp[0], clamp[1], clamp[2], bronze.0, bronze.1, bronze.2
        );
    }
    println!(
        "{{\"cymbal\":\"{}\",\"summary\":true,\"resolved_clamp_dominated\":{},\"resolved_bronze_dominated\":{},\"resolved_undecided\":{},\"unresolved_below_6_per_wavelength\":{coarse},\"default_may_switch\":{},\"clamp_or_undecided_hz\":[{}],\"flags\":[\"composition-proxy (felt: piano hammer wool to washer felt; bronze: Eijsbouts bell bronze to B20)\",\"strain-regime-proxy (felt)\",\"extrapolate-frequency (bronze above 3 kHz)\",\"ks declared estimate\"],\"assembly_s\":{assembled_s:.1},\"total_s\":{solved_s:.1},\"process_peak_rss_mb\":{}}}",
        spec.name,
        counts[0],
        counts[1],
        counts[2],
        counts[0] + counts[2] == 0,
        dominated.join(","),
        peak_rss_mb()
    );
    Ok(())
}

/// The process's peak resident set so far [MB] (Linux VmHWM), or "null" where /proc is absent.
/// The cymbals run smallest first and free their dense arrays between, so after each cymbal it
/// is that cymbal's peak.
fn peak_rss_mb() -> String {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find_map(|l| l.strip_prefix("VmHWM:"))
                .and_then(|v| v.trim().trim_end_matches("kB").trim().parse::<f64>().ok())
        })
        .map_or_else(|| "null".into(), |kb| format!("{:.0}", kb / 1024.0))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let names: Vec<&str> = if args.is_empty() {
        vec!["splash", "crash", "ride"]
    } else {
        args.iter().map(String::as_str).collect()
    };
    let bronze = bronze_loss();
    let mut taus = felt_relaxation_times();
    taus.sort_by(f64::total_cmp);
    let tau = (taus[0], taus[taus.len() / 2], taus[taus.len() - 1]);
    println!(
        "{{\"bronze_inverse_variance_mean\":{:e},\"bronze_min\":{:e},\"bronze_max\":{:e},\"felt_r_over_k_s_min_median_max\":[{:e},{:e},{:e}]}}",
        bronze.0, bronze.1, bronze.2, tau.0, tau.1, tau.2
    );
    let mut failed = false;
    for name in names {
        let spec = match name {
            "splash" => kit::splash(),
            "crash" => kit::crash(),
            "ride" => kit::ride(),
            other => {
                eprintln!("unknown cymbal {other}");
                std::process::exit(2);
            }
        };
        if let Err(e) = analyse(&spec, bronze, tau) {
            eprintln!("refused: {e}");
            failed = true;
        }
    }
    if failed {
        std::process::exit(1);
    }
}
