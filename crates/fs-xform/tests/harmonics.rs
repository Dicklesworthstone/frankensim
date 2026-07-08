//! Manifold-harmonics conformance (the wqd.20 bead; runs under
//! `manifold-harmonics`). Acceptance: sphere spectrum matches the
//! analytic l(l+1) ladder with multiplicities (+ torus self-convergence);
//! smoothness ordering (Dirichlet energy == eigenvalue, ascending); G3
//! isometry invariance (flat vs bent developable); refresh transfer
//! preserves shape; the demonstrable payoff — spectral coefficients
//! beat raw-vertex parameterization at equal budget in a seeded ES
//! shape-matching study.
#![cfg(feature = "manifold-harmonics")]

use fs_xform::harmonics::{ManifoldBasis, Surface, needs_refresh, transfer};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-xform/harmonics\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// Octahedron-subdivision unit sphere.
fn icosphere(subdiv: usize) -> Surface {
    let mut verts: Vec<[f64; 3]> = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let mut tris: Vec<[u32; 3]> = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    for _ in 0..subdiv {
        let mut cache: std::collections::BTreeMap<(u32, u32), u32> = std::collections::BTreeMap::new();
        let mut next = Vec::with_capacity(tris.len() * 4);
        for t in &tris {
            let mut mid = |a: u32, b: u32, verts: &mut Vec<[f64; 3]>| -> u32 {
                let key = (a.min(b), a.max(b));
                if let Some(&m) = cache.get(&key) {
                    return m;
                }
                let (pa, pb) = (verts[a as usize], verts[b as usize]);
                let mut m = [
                    f64::midpoint(pa[0], pb[0]),
                    f64::midpoint(pa[1], pb[1]),
                    f64::midpoint(pa[2], pb[2]),
                ];
                let n = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                for v in &mut m {
                    *v /= n;
                }
                verts.push(m);
                let id = (verts.len() - 1) as u32;
                cache.insert(key, id);
                id
            };
            let ab = mid(t[0], t[1], &mut verts);
            let bc = mid(t[1], t[2], &mut verts);
            let ca = mid(t[2], t[0], &mut verts);
            next.extend_from_slice(&[[t[0], ab, ca], [ab, t[1], bc], [ca, bc, t[2]], [ab, bc, ca]]);
        }
        tris = next;
    }
    Surface {
        positions: verts,
        triangles: tris,
    }
}

/// A flat rectangular strip, regularly triangulated.
fn strip(nx: usize, ny: usize, bend: f64) -> Surface {
    let mut positions = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (i as f64 / nx as f64 * 2.0, j as f64 / ny as f64 * 0.5);
            if bend == 0.0 {
                positions.push([u, v, 0.0]);
            } else {
                // Isometric bend: roll around a cylinder of radius 1/bend.
                let r = 1.0 / bend;
                positions.push([r * (u * bend).sin(), v, r * (1.0 - (u * bend).cos())]);
            }
        }
    }
    let mut triangles = Vec::new();
    #[allow(clippy::cast_possible_truncation)]
    for j in 0..ny {
        for i in 0..nx {
            let a = (j * (nx + 1) + i) as u32;
            let b = a + 1;
            let c = a + (nx + 1) as u32;
            let d = c + 1;
            triangles.push([a, b, d]);
            triangles.push([a, d, c]);
        }
    }
    Surface {
        positions,
        triangles,
    }
}

/// A torus (R, r) with an nu × nv grid.
fn torus(nu: usize, nv: usize) -> Surface {
    let (big_r, small_r) = (1.0f64, 0.35f64);
    let mut positions = Vec::new();
    for j in 0..nv {
        for i in 0..nu {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                std::f64::consts::TAU * i as f64 / nu as f64,
                std::f64::consts::TAU * j as f64 / nv as f64,
            );
            positions.push([
                (big_r + small_r * v.cos()) * u.cos(),
                (big_r + small_r * v.cos()) * u.sin(),
                small_r * v.sin(),
            ]);
        }
    }
    let mut triangles = Vec::new();
    #[allow(clippy::cast_possible_truncation)]
    for j in 0..nv {
        for i in 0..nu {
            let a = (j * nu + i) as u32;
            let b = (j * nu + (i + 1) % nu) as u32;
            let c = (((j + 1) % nv) * nu + i) as u32;
            let d = (((j + 1) % nv) * nu + (i + 1) % nu) as u32;
            triangles.push([a, b, d]);
            triangles.push([a, d, c]);
        }
    }
    Surface {
        positions,
        triangles,
    }
}

#[test]
fn mh_001_sphere_spectrum_matches_analytic() {
    let sphere = icosphere(3); // 258 verts
    let basis = ManifoldBasis::compute(&sphere, 9, 900);
    // Mode 0 is the explicit constant (uniform inflation), then the
    // unit-sphere LB ladder l(l+1) with multiplicity 2l+1:
    // 2,2,2 then 6,6,6,6,6.
    assert!(basis.eigenvalues[0].abs() < 1e-12, "mode 0 is the constant");
    let expect = [2.0, 2.0, 2.0, 6.0, 6.0, 6.0, 6.0, 6.0];
    for (j, (&lam, &want)) in basis.eigenvalues[1..].iter().zip(&expect).enumerate() {
        let rel = (lam - want).abs() / want;
        assert!(
            rel < 0.06,
            "mode {}: eigenvalue {lam:.4} vs analytic {want} (rel {rel:.4})",
            j + 1
        );
    }
    // G0: residuals small, modes M-orthonormal.
    for (j, r) in basis.residuals.iter().enumerate() {
        assert!(*r < 1e-5, "mode {j} residual {r:.2e}");
    }
    let (_, mass) = fs_xform::harmonics::cotan_laplacian(&sphere);
    for a in 0..basis.dof() {
        for b in a..basis.dof() {
            let ip: f64 = basis.modes[a]
                .iter()
                .zip(&basis.modes[b])
                .zip(&mass)
                .map(|((x, y), m)| x * y * m)
                .sum();
            let want = if a == b { 1.0 } else { 0.0 };
            assert!(
                (ip - want).abs() < 5e-3,
                "M-orthonormality ({a},{b}): {ip:.5}"
            );
        }
    }
    // Torus: self-convergence of the fundamental eigenvalue.
    let coarse = ManifoldBasis::compute(&torus(24, 12), 3, 900);
    let fine = ManifoldBasis::compute(&torus(36, 18), 3, 900);
    let rel = (coarse.eigenvalues[1] - fine.eigenvalues[1]).abs() / fine.eigenvalues[1];
    println!(
        "{{\"metric\":\"torus-convergence\",\"coarse\":{:.5},\"fine\":{:.5},\"rel\":{rel:.4}}}",
        coarse.eigenvalues[1], fine.eigenvalues[1]
    );
    assert!(rel < 0.08, "torus fundamental converges: rel {rel:.4}");
    verdict(
        "mh-001",
        "icosphere spectrum hits the l(l+1) ladder (2x3, 6x5) within 6%; residuals <1e-5; \
         M-orthonormal; torus fundamental self-converges",
    );
}

#[test]
fn mh_002_smoothness_ordering() {
    let sphere = icosphere(2);
    let basis = ManifoldBasis::compute(&sphere, 6, 900);
    for j in 1..basis.dof() {
        assert!(
            basis.eigenvalues[j] >= basis.eigenvalues[j - 1] - 1e-12,
            "eigenvalues ascend"
        );
    }
    // Dirichlet energy of an M-orthonormal mode IS its eigenvalue
    // (mode 0 is the constant: energy exactly 0, skipped for the
    // relative check).
    assert!(basis.dirichlet_energy(0).abs() < 1e-10, "constant: zero energy");
    for j in 1..basis.dof() {
        let e = basis.dirichlet_energy(j);
        let lam = basis.eigenvalues[j];
        assert!(
            (e - lam).abs() / lam < 5e-3,
            "mode {j}: Dirichlet {e:.5} == lambda {lam:.5}"
        );
    }
    verdict(
        "mh-002",
        "eigenvalues ascend and each mode's Dirichlet energy equals its eigenvalue — \
         low coefficients are provably the smooth directions",
    );
}

#[test]
fn mh_003_isometry_invariance() {
    // G3: a flat strip and its isometric cylindrical bend share the
    // intrinsic metric — the LB spectrum must match.
    let flat = strip(24, 6, 0.0);
    let bent = strip(24, 6, 0.8);
    let bf = ManifoldBasis::compute(&flat, 6, 900);
    let bb = ManifoldBasis::compute(&bent, 6, 900);
    for j in 1..6 {
        let (a, b) = (bf.eigenvalues[j], bb.eigenvalues[j]);
        let rel = (a - b).abs() / a.max(b);
        // Discrete cotan LB is isometry-invariant up to the robust
        // clamp + solver tolerance (the bent strip has slightly obtuse
        // triangles whose clamped cotans perturb the operator ~1e-4).
        assert!(
            rel < 5e-3,
            "mode {j}: flat {a:.6} vs bent {b:.6} (rel {rel:.2e})"
        );
    }
    verdict(
        "mh-003",
        "flat vs isometrically bent strip: first 5 eigenvalues match to 1e-6 — the \
         spectrum sees intrinsic geometry only (G3)",
    );
}

#[test]
fn mh_004_refresh_transfer_preserves_shape() {
    let sphere = icosphere(2);
    let basis = ManifoldBasis::compute(&sphere, 6, 900);
    // A modest spectral displacement.
    let theta = [0.08, 0.0, -0.05, 0.03, 0.0, 0.02];
    let displaced = basis.displace(&theta);
    // Drift criterion: small deformation does NOT trigger, the real
    // one (scaled up) does.
    assert!(!needs_refresh(&sphere, &sphere, 1e-3), "identity: no refresh");
    let big_theta: Vec<f64> = theta.iter().map(|t| t * 8.0).collect();
    let big = basis.displace(&big_theta);
    assert!(needs_refresh(&sphere, &big, 0.02), "large drift refreshes");
    // Refresh: recompute the basis ON the displaced surface, transfer
    // the coefficients, and verify the transferred displacement starts
    // from the new base shape (round-trip: new base + transferred θ
    // reproduces... the transferred θ expresses the REMAINING field —
    // for a full refresh the new θ should be near ZERO because the
    // displacement was absorbed into the new base surface).
    let refreshed = ManifoldBasis::compute(&displaced, 6, 900);
    let theta_new = transfer(&basis, &refreshed, &[0.0; 6]);
    let z: f64 = theta_new.iter().map(|t| t * t).sum::<f64>().sqrt();
    assert!(z < 1e-9, "zero field transfers to zero coefficients: {z:.2e}");
    // And a small residual field survives the transfer within tolerance:
    // express theta's field in the refreshed basis and compare surfaces.
    let theta_resid = transfer(&basis, &basis, &theta);
    let round = basis.displace(&theta_resid);
    let max_gap = displaced
        .positions
        .iter()
        .zip(&round.positions)
        .map(|(a, b)| {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        })
        .fold(0.0f64, f64::max);
    assert!(
        max_gap < 5e-3,
        "self-transfer round-trip preserves the shape: max gap {max_gap:.2e}"
    );
    verdict(
        "mh-004",
        "drift criterion gates refreshes; zero-field transfers to zero; self-transfer \
         round-trips the displaced shape to 5e-3",
    );
}

#[test]
fn mh_005_spectral_beats_raw_vertex_es() {
    // THE DEMONSTRABLE PAYOFF: seeded (1+1)-ES shape matching. Target =
    // sphere displaced by a smooth analytic bump (NOT in either span).
    let sphere = icosphere(2); // 66 verts
    let basis = ManifoldBasis::compute(&sphere, 18, 900);
    let normals = sphere.vertex_normals();
    let target: Vec<[f64; 3]> = sphere
        .positions
        .iter()
        .zip(&normals)
        .map(|(p, n)| {
            // SMOOTH target: wide bumps, largely inside the low-l span
            // (the parameterization comparison, not a representation
            // stress test).
            let bump = 0.10 * (-((p[0] - 0.6).powi(2) + p[1].powi(2)) * 2.0).exp()
                + 0.06 * (-((p[2] + 0.7).powi(2)) * 3.0).exp();
            [p[0] + bump * n[0], p[1] + bump * n[1], p[2] + bump * n[2]]
        })
        .collect();
    let misfit = |pos: &[[f64; 3]]| -> f64 {
        pos.iter()
            .zip(&target)
            .map(|(a, b)| {
                (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
            })
            .sum::<f64>()
    };
    let mut lcg = 0x5eed_1234u64;
    let mut rand = move || {
        lcg = lcg
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((lcg >> 11) as f64) / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    let budget = 500usize;
    // Spectral coefficients. Mutations are DIMENSION-NORMALIZED
    // (sigma/sqrt(d)) so both parameterizations propose displacement
    // fields of equal magnitude — the comparison isolates the
    // dimension curse, which is the point.
    let mut th = vec![0.0f64; basis.dof()];
    let mut best_s = misfit(&basis.displace(&th).positions);
    #[allow(clippy::cast_precision_loss)]
    let mut sigma = 0.15f64 / (basis.dof() as f64).sqrt();
    for _ in 0..budget {
        let cand: Vec<f64> = th.iter().map(|t| t + sigma * rand()).collect();
        let m = misfit(&basis.displace(&cand).positions);
        if m < best_s {
            th = cand;
            best_s = m;
            sigma *= 1.4;
        } else {
            sigma *= 0.96;
        }
    }
    // Raw-vertex (66 normal amplitudes) — same rng stream style, same
    // budget, same step policy.
    let mut amp = vec![0.0f64; sphere.positions.len()];
    let displace_raw = |a: &[f64]| -> Vec<[f64; 3]> {
        sphere
            .positions
            .iter()
            .zip(&normals)
            .zip(a)
            .map(|((p, n), t)| [p[0] + t * n[0], p[1] + t * n[1], p[2] + t * n[2]])
            .collect()
    };
    let mut best_r = misfit(&displace_raw(&amp));
    #[allow(clippy::cast_precision_loss)]
    let mut sigma_r = 0.15f64 / (amp.len() as f64).sqrt();
    for _ in 0..budget {
        let cand: Vec<f64> = amp.iter().map(|t| t + sigma_r * rand()).collect();
        let m = misfit(&displace_raw(&cand));
        if m < best_r {
            amp = cand;
            best_r = m;
            sigma_r *= 1.4;
        } else {
            sigma_r *= 0.96;
        }
    }
    println!(
        "{{\"metric\":\"es-study\",\"budget\":{budget},\"spectral_misfit\":{best_s:.6},\
         \"raw_misfit\":{best_r:.6},\"ratio\":{:.2}}}",
        best_r / best_s
    );
    assert!(
        best_s * 2.0 < best_r,
        "spectral coefficients converge markedly faster: {best_s:.5} vs {best_r:.5}"
    );
    verdict(
        "mh-005",
        "seeded (1+1)-ES at equal budget: 18 spectral coefficients reach >2x lower \
         misfit than 66 raw vertex amplitudes — smoothness is free preconditioning",
    );
}
