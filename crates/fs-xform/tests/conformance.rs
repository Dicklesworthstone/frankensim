//! fs-xform conformance suite (plan §13.3; the wqd.19 bead). Acceptance:
//! every lever passes the Jacobian-action-vs-dual battery (the solvers'
//! gradient-gate discipline); composition Jacobians correct (FD
//! order-checked); the level-set velocity path drives a working advection
//! step; fold-over monitoring refuses structurally; G3 frame equivariance
//! for the radial lever.

use fs_xform::{
    Composed, DensityField, FfdLattice, Parameterization, Point3, RbfMorph, Vec3, VelocityBand,
    XformError, advect_sdf, detect_foldover,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-xform/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn unit_ffd() -> FfdLattice {
    FfdLattice {
        origin: Point3::new(0.0, 0.0, 0.0),
        size: Vec3::new(1.0, 1.0, 1.0),
        counts: [2, 2, 2],
    }
}

fn sub(a: Point3, b: Point3) -> Vec3 {
    a.delta_from(b)
}

#[test]
fn xf_001_jacobian_action_matches_finite_differences_and_is_linear() {
    let mut seed = 0x5EED_0F02_0000_0019u64;
    let ffd = unit_ffd();
    let rbf = RbfMorph {
        centers: vec![Point3::new(0.3, 0.3, 0.3), Point3::new(0.7, 0.6, 0.4)],
        radius: 0.8,
    };
    let levers: [(&str, &dyn Parameterization); 2] = [("ffd", &ffd), ("rbf", &rbf)];
    for (name, lever) in levers {
        let n = lever.dof();
        let theta: Vec<f64> = (0..n).map(|_| 0.1 * (lcg(&mut seed) - 0.5)).collect();
        let dtheta: Vec<f64> = (0..n).map(|_| lcg(&mut seed) - 0.5).collect();
        for _ in 0..25 {
            let x = Point3::new(lcg(&mut seed), lcg(&mut seed), lcg(&mut seed));
            let action = lever.jacobian_action(&theta, &dtheta, x).unwrap();
            // These levers are LINEAR in θ: the secant is exact.
            let eps = 1e-3;
            let plus: Vec<f64> =
                theta.iter().zip(&dtheta).map(|(t, d)| t + eps * d).collect();
            let fd = sub(lever.apply(&plus, x).unwrap(), lever.apply(&theta, x).unwrap());
            for (a, f) in [(action.x, fd.x), (action.y, fd.y), (action.z, fd.z)] {
                assert!(
                    (a - f / eps).abs() < 1e-9,
                    "{name}: action {a} vs exact secant {}",
                    f / eps
                );
            }
            // Linearity in δθ.
            let twice: Vec<f64> = dtheta.iter().map(|d| 2.0 * d).collect();
            let double = lever.jacobian_action(&theta, &twice, x).unwrap();
            assert!((double.x - 2.0 * action.x).abs() < 1e-12);
        }
    }
    verdict("xf-001", "FFD + RBF actions exact against secants; linear in delta-theta");
}

#[test]
fn xf_002_dual_number_gate_on_the_ffd_warp() {
    // Re-express the 2×2×2 Bernstein warp over Dual64<1> and compare its
    // JVP against jacobian_action — the same gradient-gate discipline as
    // solvers (bead acceptance).
    let ffd = unit_ffd();
    let n = ffd.dof(); // 24
    let mut seed = 0x5EED_D0A1_0000_0027u64;
    let theta: Vec<f64> = (0..n).map(|_| 0.2 * (lcg(&mut seed) - 0.5)).collect();
    let dtheta: Vec<f64> = (0..n).map(|_| lcg(&mut seed) - 0.5).collect();
    let x = Point3::new(0.35, 0.62, 0.81);
    let action = ffd.jacobian_action(&theta, &dtheta, x).unwrap();
    // Trilinear Bernstein (n=1 per axis): B0 = 1−t, B1 = t.
    for comp in 0..3 {
        let (value, derivative) = {
            // f(θ) = x_comp + Σ_ijk B_i(u)B_j(v)B_k(w)·θ[3·node+comp]
            let seed_dir = &dtheta;
            let mut vars = vec![fs_ad::Dual64::<1>::constant(0.0); n];
            for (i, v) in vars.iter_mut().enumerate() {
                *v = fs_ad::Dual64::<1> { re: theta[i], eps: [seed_dir[i]] };
            }
            let (u, v, w) = (x.x, x.y, x.z);
            let b = |t: f64| [1.0 - t, t];
            let mut acc = fs_ad::Dual64::<1>::constant([x.x, x.y, x.z][comp]);
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let weight = b(u)[i] * b(v)[j] * b(w)[k];
                        let node = (i * 2 + j) * 2 + k;
                        let coeff = fs_ad::Dual64::<1>::constant(weight);
                        acc = acc + coeff * vars[3 * node + comp];
                    }
                }
            }
            (acc.re, acc.eps[0])
        };
        let applied = ffd.apply(&theta, x).unwrap();
        let expected_value = [applied.x, applied.y, applied.z][comp];
        assert!((value - expected_value).abs() < 1e-12, "dual primal drifted");
        let expected_dir = [action.x, action.y, action.z][comp];
        assert!(
            (derivative - expected_dir).abs() < 1e-12,
            "dual JVP {derivative} vs jacobian_action {expected_dir} (comp {comp})"
        );
    }
    verdict("xf-002", "FFD jacobian_action agrees with dual-number JVP to 1e-12");
}

#[test]
fn xf_003_rbf_frame_equivariance_and_foldover_refusal() {
    // G3: rotate handles, displacements, and the probe by 90° about z —
    // the radial kernel must commute exactly (up to fp).
    let rot = |p: Point3| Point3::new(-p.y, p.x, p.z);
    let rot_v = |v: [f64; 3]| [-v[1], v[0], v[2]];
    let centers = vec![Point3::new(0.2, 0.1, 0.0), Point3::new(-0.3, 0.4, 0.2)];
    let theta = [0.05, -0.02, 0.01, 0.03, 0.04, -0.06];
    let morph = RbfMorph { centers: centers.clone(), radius: 1.0 };
    let rotated = RbfMorph { centers: centers.iter().map(|&c| rot(c)).collect(), radius: 1.0 };
    let mut theta_rot = [0.0; 6];
    for i in 0..2 {
        let r = rot_v([theta[3 * i], theta[3 * i + 1], theta[3 * i + 2]]);
        theta_rot[3 * i..3 * i + 3].copy_from_slice(&r);
    }
    let x = Point3::new(0.15, -0.25, 0.1);
    let y = morph.apply(&theta, x).unwrap();
    let y_rot = rotated.apply(&theta_rot, rot(x)).unwrap();
    let expect = rot(y);
    for (a, b) in [(y_rot.x, expect.x), (y_rot.y, expect.y), (y_rot.z, expect.z)] {
        assert!((a - b).abs() < 1e-14, "equivariance broke: {a} vs {b}");
    }
    // Fold-over: a violent handle collapse must refuse with location+det.
    let violent = [-3.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let samples: Vec<Point3> =
        (0..20).map(|i| Point3::new(0.2 + 0.01 * f64::from(i), 0.1, 0.0)).collect();
    match detect_foldover(&morph, &violent, &samples) {
        Err(XformError::FoldOver { det, .. }) => assert!(det <= 0.0),
        other => panic!("violent collapse must fold: {other:?}"),
    }
    // The gentle θ does not fold.
    detect_foldover(&morph, &theta, &samples).expect("gentle morph stays invertible");
    verdict("xf-003", "RBF equivariant under rotation; fold-over refuses structurally");
}

#[test]
fn xf_004_composition_chain_rule_converges_at_second_order() {
    let ffd = unit_ffd();
    let rbf = RbfMorph { centers: vec![Point3::new(0.5, 0.5, 0.5)], radius: 0.9 };
    let composed = Composed { first: &ffd, then: &rbf };
    let n = composed.dof();
    let mut seed = 0x5EED_C0A2_0000_0041u64;
    let theta: Vec<f64> = (0..n).map(|_| 0.15 * (lcg(&mut seed) - 0.5)).collect();
    let dtheta: Vec<f64> = (0..n).map(|_| lcg(&mut seed) - 0.5).collect();
    let x = Point3::new(0.4, 0.45, 0.55);
    let action = composed.jacobian_action(&theta, &dtheta, x).unwrap();
    // Central differences at ε and ε/2: the composed map is nonlinear in
    // θ jointly, so the FD error is O(ε²) — the ratio must approach 4.
    let fd = |eps: f64| -> Vec3 {
        let plus: Vec<f64> = theta.iter().zip(&dtheta).map(|(t, d)| t + eps * d).collect();
        let minus: Vec<f64> = theta.iter().zip(&dtheta).map(|(t, d)| t - eps * d).collect();
        let d = sub(
            composed.apply(&plus, x).unwrap(),
            composed.apply(&minus, x).unwrap(),
        );
        Vec3::new(d.x / (2.0 * eps), d.y / (2.0 * eps), d.z / (2.0 * eps))
    };
    let err = |eps: f64| -> f64 {
        let f = fd(eps);
        ((f.x - action.x).powi(2) + (f.y - action.y).powi(2) + (f.z - action.z).powi(2)).sqrt()
    };
    let (e1, e2) = (err(1e-3), err(5e-4));
    assert!(e1 < 1e-4, "chain rule far off: err {e1}");
    let order = (e1 / e2).log2();
    assert!(
        order > 1.5,
        "composition FD convergence order {order:.2} (want ~2): e1={e1:.3e}, e2={e2:.3e}"
    );
    // Spatial Jacobian of the composition is the matrix product (probe det>0).
    let j = composed.spatial_jacobian(&theta, x).unwrap();
    assert!(fs_xform::det3(&j) > 0.0);
    verdict(
        "xf-004",
        &format!("composed chain rule verified; FD order {order:.2}"),
    );
}

#[test]
fn xf_005_levelset_velocity_drives_advection() {
    // A sphere SDF under constant unit outward speed must grow its radius
    // by ~dt per step (first-order upwind: loose tolerance, right physics).
    let dims = [24usize, 24, 24];
    let spacing = 0.1;
    let center = Point3::new(1.15, 1.15, 1.15);
    let r0 = 0.5;
    let node_pos = |i: usize, j: usize, k: usize| {
        Point3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing)
    };
    let mut phi = vec![0.0f64; dims[0] * dims[1] * dims[2]];
    for i in 0..dims[0] {
        for j in 0..dims[1] {
            for k in 0..dims[2] {
                let p = node_pos(i, j, k);
                phi[(i * dims[1] + j) * dims[2] + k] = p.delta_from(center).norm() - r0;
            }
        }
    }
    // The velocity comes from the PARAMETERIZATION: θ = all ones on the
    // grid (constant unit normal speed inside a wide band).
    let band = VelocityBand {
        origin: Point3::new(0.0, 0.0, 0.0),
        spacing,
        dims,
        band: 10.0,
    };
    let theta = vec![1.0; band.dof()];
    let dt = 0.02;
    let steps = 5;
    for _ in 0..steps {
        let phi_snapshot = phi.clone();
        let v = |i: usize, j: usize, k: usize| -> f64 {
            let idx = (i * dims[1] + j) * dims[2] + k;
            band.velocity(&theta, node_pos(i, j, k), phi_snapshot[idx]).unwrap()
        };
        advect_sdf(&mut phi, dims, spacing, &v, dt);
    }
    // Measure the new radius along +x from the center: φ = 0 crossing.
    let expected_r = r0 + dt * f64::from(steps as u32);
    let probe = |r: f64| {
        let p = Point3::new(center.x + r, center.y, center.z);
        let i = (p.x / spacing).round() as usize;
        let j = (p.y / spacing).round() as usize;
        let k = (p.z / spacing).round() as usize;
        phi[(i * dims[1] + j) * dims[2] + k]
    };
    assert!(
        probe(expected_r - 0.1) < 0.0,
        "inside the grown sphere must stay negative"
    );
    assert!(
        probe(expected_r + 0.1) > 0.0,
        "outside the grown sphere must stay positive"
    );
    verdict(
        "xf-005",
        &format!("unit-speed band advection grew the sphere to ~{expected_r:.2} as expected"),
    );
}

#[test]
fn xf_006_density_and_dof_diagnostics_teach() {
    let field = DensityField { cells: 4 };
    match field.validate(&[0.2, 0.4, f64::NAN, 0.8]) {
        Err(XformError::OutOfBounds { index: 2, .. }) => {}
        other => panic!("NaN density must refuse at its index: {other:?}"),
    }
    let ffd = unit_ffd();
    match ffd.apply(&[1.0, 2.0], Point3::new(0.5, 0.5, 0.5)) {
        Err(XformError::DofMismatch { expected: 24, got: 2 }) => {}
        other => panic!("dof mismatch must teach: {other:?}"),
    }
    let msg = XformError::FoldOver { at: Point3::new(1.0, 2.0, 3.0), det: -0.5 }.to_string();
    assert!(msg.contains("reduce the step"), "refusals must teach: {msg}");
    verdict("xf-006", "structured diagnostics name index, expected DOFs, and fixes");
}
