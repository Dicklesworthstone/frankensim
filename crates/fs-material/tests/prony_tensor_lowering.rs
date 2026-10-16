//! Integration and Gauntlet (G0 / G1 / G3) tests for `fs_material::visco::tensor::lower`
//! (IRP M0.2-W2b: `frankensim-irp-v5-instrument-realism-pexfp.22`).

use std::time::Instant;

use fs_la::eigen::jacobi_eigh;
use fs_material::OrthotropicElastic;
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::tensor::{rotate_mandel, rotation};
use fs_material::visco::tensor::{
    CERTIFICATE_GRID_POINTS, ComplexEngineering, KelvinFractional, MAX_PRONY_BRANCHES, Mandel6,
    PlaneStressTarget, PronyBranch, TensorComplexModulus, TensorLoweringTolerance, TensorProny,
    certify_lowered_tensor, lower_tensor_to_prony, lower_tensor_to_prony_with_max_iters,
    relative_quadratic_errors,
};
use fs_material::visco::{FractionalZener, ViscoError, lower_to_prony};
use fs_math::det;

fn unit_axis(v: [f64; 3]) -> [f64; 3] {
    let n = det::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Construct a FractionalZener with equilibrium modulus `e0_pa`, fractional order `alpha`,
/// relaxation time `tau_s`, and peak loss factor `eta_peak` at `ω = 1 / tau_s`.
fn zener_with_peak_eta(e0_pa: f64, alpha: f64, tau_s: f64, eta_peak: f64) -> FractionalZener {
    // At ω·τ = 1, E*(1/τ) = E0 + (E∞ - E0)·(1 + c + i s)/(2(1 + c)) where c = cos(απ/2), s = sin(απ/2).
    // Since Re((i)^α / (1 + (i)^α)) = 1/2 and Im((i)^α / (1 + (i)^α)) = s / (2(1 + c)),
    // η(1/τ) = ΔE · (s / (2(1 + c))) / (E0 + ΔE / 2) = eta_peak.
    let half_pi = 0.5 * core::f64::consts::PI;
    let c = det::cos(alpha * half_pi);
    let s = det::sin(alpha * half_pi);
    let k_im = s / (2.0 * (1.0 + c));
    let delta_e = (eta_peak * e0_pa) / (k_im - 0.5 * eta_peak);
    FractionalZener::new(e0_pa, e0_pa + delta_e, alpha, tau_s).expect("valid FractionalZener")
}

/// Fixture 1: 6-mode KelvinFractional in a rotated orthotropic Sitka-spruce eigenbasis.
/// α ∈ [0.25, 0.50], η plateau ∈ [0.005, 0.02].
fn kelvin_wood_fixture() -> KelvinFractional {
    let ortho = OrthotropicElastic::new(
        [11.0e9, 0.90e9, 0.50e9],
        [0.37, 0.47, 0.43],
        [0.70e9, 0.06e9, 0.65e9],
        0.01,
    )
    .expect("admissible orthotropic baseline");
    let q = rotation(unit_axis([1.0, 2.0, -1.5]), 0.48);
    let c_ref: Mandel6 = rotate_mandel(ortho.stiffness(), q);
    let flat: Vec<f64> = c_ref.iter().flatten().copied().collect();
    let (evals, evecs) = jacobi_eigh(&flat, 6);
    let mut basis = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for k in 0..6 {
            basis[i][k] = evecs[i * 6 + k];
        }
    }
    let alphas = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50];
    let etas = [0.005, 0.008, 0.011, 0.014, 0.017, 0.020];
    let taus = [1.2e-3, 9.0e-4, 7.0e-4, 5.0e-4, 3.5e-4, 2.5e-4];
    let modes = core::array::from_fn(|k| zener_with_peak_eta(evals[k], alphas[k], taus[k], etas[k]));
    KelvinFractional::new(basis, modes).expect("orthogonal KelvinFractional")
}

/// Fixture 2: Non-commuting orthotropic ComplexEngineering wood target with
/// E_L* (η ≈ 0.006), E_R* (η ≈ 0.02), G_LR* (η ≈ 0.015), rotated off-axis so all 6
/// Mandel coordinates couple and the eigenbasis drifts with ω.
fn complex_engineering_wood_fixture() -> ComplexEngineering {
    let e_l = zener_with_peak_eta(11.2e9, 0.32, 6.0e-4, 0.006);
    let e_r = zener_with_peak_eta(0.92e9, 0.38, 8.5e-4, 0.020);
    let e_t = zener_with_peak_eta(0.52e9, 0.42, 1.1e-3, 0.024);
    let g_lr = zener_with_peak_eta(0.72e9, 0.35, 7.5e-4, 0.015);
    let g_rt = zener_with_peak_eta(0.065e9, 0.44, 1.2e-3, 0.026);
    let g_tl = zener_with_peak_eta(0.66e9, 0.36, 7.0e-4, 0.016);
    let q = rotation(unit_axis([0.8, -1.1, 0.6]), 0.37);
    ComplexEngineering::new([e_l, e_r, e_t], [0.37, 0.47, 0.43], [g_lr, g_rt, g_tl])
        .expect("admissible ComplexEngineering")
        .with_rotation(q)
        .expect("valid rotation")
}

/// Fixture 3: 3-D anisotropic wood TensorProny law for `PlaneStressTarget` in `MaterialPlane::Axes12`.
fn wood_3d_prony_fixture() -> TensorProny {
    let eng = complex_engineering_wood_fixture();
    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    lower_tensor_to_prony(&eng, 20.0, 20_000.0, 24, tol)
        .expect("3D wood law lowers")
        .law
}

fn assert_all_branches_psd(law: &TensorProny) -> f64 {
    assert!(!law.branches().is_empty(), "law must have at least one branch");
    let mut min_eig = f64::INFINITY;
    for (idx, b) in law.branches().iter().enumerate() {
        let flat: Vec<f64> = b.c_pa.iter().flatten().copied().collect();
        let (evals, _) = jacobi_eigh(&flat, 6);
        min_eig = min_eig.min(evals[0]);
        assert!(
            evals[0] >= 0.0,
            "branch {idx} must have min eigenvalue >= 0.0 after projection, got {:e}",
            evals[0]
        );
    }
    min_eig
}

#[test]
fn g0_scalar_target_reproduces_lower_to_prony_bitwise() {
    let t0 = Instant::now();
    // Same FractionalZener as visco.rs:1279 (`lowering_reproduces_constant_loss_plateau_with_eight_terms`).
    let fz = FractionalZener::new(8.0e9, 16.0e9, 0.35, 1e-3).expect("valid FractionalZener");
    let scalar = lower_to_prony(&fz, 20.0, 2.0e5, 8, 0.02).expect("scalar lower_to_prony");
    let tol = TensorLoweringTolerance {
        storage_rel: 0.02,
        loss_rel: 0.25,
    };
    let tensor = lower_tensor_to_prony(&fz, 20.0, 2.0e5, 8, tol).expect("tensor lower_to_prony");

    assert_eq!(tensor.law.branches().len(), scalar.model.terms.len());
    assert_eq!(
        tensor.law.c_eq()[0][0].to_bits(),
        scalar.model.e_inf.to_bits()
    );
    assert_eq!(
        tensor.certificate.band_rad_s.0.to_bits(),
        scalar.band.0.to_bits()
    );
    assert_eq!(
        tensor.certificate.band_rad_s.1.to_bits(),
        scalar.band.1.to_bits()
    );
    assert_eq!(tensor.certificate.grid_points, scalar.verification_points);
    assert_eq!(
        tensor.certificate.sup_rel_err.to_bits(),
        scalar.sup_rel_err.to_bits()
    );
    assert_eq!(
        tensor.certificate.grid_sup_storage.to_bits(),
        scalar.sup_rel_storage_error(&fz).to_bits()
    );
    assert_eq!(
        tensor.certificate.grid_sup_loss.to_bits(),
        scalar.sup_rel_loss_error(&fz).to_bits()
    );

    for (branch, &(e_scalar, tau_scalar)) in
        tensor.law.branches().iter().zip(&scalar.model.terms)
    {
        assert_eq!(branch.c_pa[0][0].to_bits(), e_scalar.to_bits());
        assert_eq!(branch.tau_s.to_bits(), tau_scalar.to_bits());
    }

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g0_scalar_target_reproduces_lower_to_prony_bitwise\",\"tier\":\"G0\",\"terms\":{},\"sup_rel_err\":{:.6e},\"grid_sup_loss\":{:.6e},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        tensor.law.branches().len(),
        tensor.certificate.sup_rel_err,
        tensor.certificate.grid_sup_loss,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g1_kelvin_fractional_target_lowers_per_mode_within_certificate() {
    let t0 = Instant::now();
    let target = kelvin_wood_fixture();
    let loose_tol = TensorLoweringTolerance {
        storage_rel: 0.5,
        loss_rel: 0.5,
    };
    let rungs = [4usize, 8, 16];
    let mut storage_rungs = [0.0f64; 3];
    let mut loss_rungs = [0.0f64; 3];
    for (idx, &n) in rungs.iter().enumerate() {
        let m = lower_tensor_to_prony(&target, 20.0, 20_000.0, n, loose_tol)
            .expect("rung lowering must succeed");
        storage_rungs[idx] = m.certificate.total_storage_bound();
        loss_rungs[idx] = m.certificate.total_loss_bound();
    }

    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let lowered =
        lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol).expect("8-term Kelvin lowering");
    let min_eig = assert_all_branches_psd(&lowered.law);

    assert!(
        lowered.certificate.total_storage_bound() <= 1e-3,
        "storage bound {:.6e} > 1e-3",
        lowered.certificate.total_storage_bound()
    );
    assert!(
        lowered.certificate.total_loss_bound() <= 1e-2,
        "loss bound {:.6e} > 1e-2",
        lowered.certificate.total_loss_bound()
    );
    assert!(
        loss_rungs[1] <= loss_rungs[0] && loss_rungs[2] <= loss_rungs[1],
        "loss certificate ladder must be monotone non-increasing across [4,8,16]: {loss_rungs:?}"
    );

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g1_kelvin_fractional_target_lowers_per_mode_within_certificate\",\"tier\":\"G1\",\"inputs\":{{\"target\":\"kelvin\",\"band_hz\":[20,20000],\"n_terms\":8}},\"provenance\":\"named synthetic fixture\",\"rungs\":{{\"n_terms\":[4,8,16]}},\"observed\":{{\"storage_sup\":[{:.6e},{:.6e},{:.6e}],\"loss_sup\":[{:.6e},{:.6e},{:.6e}],\"lipschitz\":{:.6e},\"min_branch_eig\":{:.6e}}},\"tolerance\":{{\"storage\":1e-3,\"loss\":1e-2}},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        storage_rungs[0],
        storage_rungs[1],
        storage_rungs[2],
        loss_rungs[0],
        loss_rungs[1],
        loss_rungs[2],
        lowered.certificate.lipschitz_term,
        min_eig,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g1_complex_engineering_constants_target_lowers_with_psd_branches() {
    let t0 = Instant::now();
    let target = complex_engineering_wood_fixture();

    // Verify that the target's storage and loss operators do NOT commute (eigenbasis drifts with ω).
    let (re_w, im_w) = target.complex_modulus(2_000.0).expect("modulus at 2000 rad/s");
    let mut comm_norm_sq = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            let mut ab = 0.0;
            let mut ba = 0.0;
            for k in 0..6 {
                ab += re_w[i][k] * im_w[k][j];
                ba += im_w[i][k] * re_w[k][j];
            }
            let d = ab - ba;
            comm_norm_sq += d * d;
        }
    }
    assert!(
        det::sqrt(comm_norm_sq) > 1e14,
        "ComplexEngineering storage and loss must be non-commuting"
    );

    let loose_tol = TensorLoweringTolerance {
        storage_rel: 0.5,
        loss_rel: 0.5,
    };
    let rungs = [4usize, 8, 16];
    let mut storage_rungs = [0.0f64; 3];
    let mut loss_rungs = [0.0f64; 3];
    for (idx, &n) in rungs.iter().enumerate() {
        let m = lower_tensor_to_prony(&target, 20.0, 20_000.0, n, loose_tol)
            .expect("rung lowering must succeed");
        storage_rungs[idx] = m.certificate.total_storage_bound();
        loss_rungs[idx] = m.certificate.total_loss_bound();
    }

    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let lowered = lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol)
        .expect("8-term ComplexEngineering lowering");
    let min_eig = assert_all_branches_psd(&lowered.law);

    assert!(
        lowered.certificate.total_storage_bound() <= 1e-3,
        "storage bound {:.6e} > 1e-3",
        lowered.certificate.total_storage_bound()
    );
    assert!(
        lowered.certificate.total_loss_bound() <= 1e-2,
        "loss bound {:.6e} > 1e-2",
        lowered.certificate.total_loss_bound()
    );
    assert!(
        loss_rungs[1] <= loss_rungs[0] && loss_rungs[2] <= loss_rungs[1],
        "loss certificate ladder must be monotone non-increasing across [4,8,16]: {loss_rungs:?}"
    );

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g1_complex_engineering_constants_target_lowers_with_psd_branches\",\"tier\":\"G1\",\"inputs\":{{\"target\":\"engineering\",\"band_hz\":[20,20000],\"n_terms\":8}},\"provenance\":\"named synthetic fixture\",\"rungs\":{{\"n_terms\":[4,8,16]}},\"observed\":{{\"storage_sup\":[{:.6e},{:.6e},{:.6e}],\"loss_sup\":[{:.6e},{:.6e},{:.6e}],\"lipschitz\":{:.6e},\"min_branch_eig\":{:.6e}}},\"tolerance\":{{\"storage\":1e-3,\"loss\":1e-2}},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        storage_rungs[0],
        storage_rungs[1],
        storage_rungs[2],
        loss_rungs[0],
        loss_rungs[1],
        loss_rungs[2],
        lowered.certificate.lipschitz_term,
        min_eig,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g1_plane_stress_schur_target_lowers() {
    let t0 = Instant::now();
    let law_3d = wood_3d_prony_fixture();
    let target = PlaneStressTarget::new(&law_3d, MaterialPlane::Axes12);

    let loose_tol = TensorLoweringTolerance {
        storage_rel: 0.5,
        loss_rel: 0.5,
    };
    let rungs = [4usize, 8, 16];
    let mut storage_rungs = [0.0f64; 3];
    let mut loss_rungs = [0.0f64; 3];
    for (idx, &n) in rungs.iter().enumerate() {
        let m = lower_tensor_to_prony(&target, 20.0, 20_000.0, n, loose_tol)
            .expect("rung lowering must succeed");
        storage_rungs[idx] = m.certificate.total_storage_bound();
        loss_rungs[idx] = m.certificate.total_loss_bound();
    }

    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let lowered = lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol)
        .expect("8-term PlaneStressTarget lowering");
    let min_eig = assert_all_branches_psd(&lowered.law);

    assert!(
        lowered.certificate.total_storage_bound() <= 1e-3,
        "storage bound {:.6e} > 1e-3",
        lowered.certificate.total_storage_bound()
    );
    assert!(
        lowered.certificate.total_loss_bound() <= 1e-2,
        "loss bound {:.6e} > 1e-2",
        lowered.certificate.total_loss_bound()
    );
    assert!(
        loss_rungs[1] <= loss_rungs[0] && loss_rungs[2] <= loss_rungs[1],
        "loss certificate ladder must be monotone non-increasing across [4,8,16]: {loss_rungs:?}"
    );

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g1_plane_stress_schur_target_lowers\",\"tier\":\"G1\",\"inputs\":{{\"target\":\"plane\",\"band_hz\":[20,20000],\"n_terms\":8}},\"provenance\":\"named synthetic fixture\",\"rungs\":{{\"n_terms\":[4,8,16]}},\"observed\":{{\"storage_sup\":[{:.6e},{:.6e},{:.6e}],\"loss_sup\":[{:.6e},{:.6e},{:.6e}],\"lipschitz\":{:.6e},\"min_branch_eig\":{:.6e}}},\"tolerance\":{{\"storage\":1e-3,\"loss\":1e-2}},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        storage_rungs[0],
        storage_rungs[1],
        storage_rungs[2],
        loss_rungs[0],
        loss_rungs[1],
        loss_rungs[2],
        lowered.certificate.lipschitz_term,
        min_eig,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g1_certificate_encloses_dense_independent_check() {
    let t0 = Instant::now();
    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let kelvin = kelvin_wood_fixture();
    let eng = complex_engineering_wood_fixture();
    let law_3d = wood_3d_prony_fixture();
    let plane = PlaneStressTarget::new(&law_3d, MaterialPlane::Axes12);

    let targets: [(&str, &dyn TensorComplexModulus); 3] = [
        ("kelvin", &kelvin),
        ("engineering", &eng),
        ("plane", &plane),
    ];

    let dense_points = 65_536usize;
    let mut worst_storage_margin = f64::INFINITY;
    let mut worst_loss_margin = f64::INFINITY;

    for (name, target) in targets {
        let lowered = lower_tensor_to_prony(target, 20.0, 20_000.0, 8, tol)
            .unwrap_or_else(|e| panic!("{name} lowering failed: {e}"));
        let cert = lowered.certificate;
        let (w_lo, w_hi) = cert.band_rad_s;
        let log_lo = det::ln(w_lo);
        let log_hi = det::ln(w_hi);

        let mut dense_sup_storage = 0.0f64;
        let mut dense_sup_loss = 0.0f64;
        for i in 0..dense_points {
            let t = i as f64 / ((dense_points - 1) as f64);
            let w = det::exp(log_lo + t * (log_hi - log_lo)).clamp(w_lo, w_hi);
            let (err_re, err_im) = relative_quadratic_errors(target, &lowered.law, w)
                .expect("dense evaluation must succeed");
            dense_sup_storage = dense_sup_storage.max(err_re);
            dense_sup_loss = dense_sup_loss.max(err_im);
        }

        let storage_margin = cert.total_storage_bound() - dense_sup_storage;
        let loss_margin = cert.total_loss_bound() - dense_sup_loss;
        worst_storage_margin = worst_storage_margin.min(storage_margin);
        worst_loss_margin = worst_loss_margin.min(loss_margin);

        assert!(
            dense_sup_storage <= cert.total_storage_bound(),
            "{name}: 65536-pt storage sup {dense_sup_storage:.6e} exceeded certificate {:.6e}",
            cert.total_storage_bound()
        );
        assert!(
            dense_sup_loss <= cert.total_loss_bound(),
            "{name}: 65536-pt loss sup {dense_sup_loss:.6e} exceeded certificate {:.6e}",
            cert.total_loss_bound()
        );
    }

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g1_certificate_encloses_dense_independent_check\",\"tier\":\"G1\",\"dense_points\":{},\"min_storage_margin\":{:.6e},\"min_loss_margin\":{:.6e},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        dense_points,
        worst_storage_margin,
        worst_loss_margin,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g3_dropping_any_branch_breaks_the_certificate() {
    let t0 = Instant::now();
    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let target = complex_engineering_wood_fixture();
    let lowered =
        lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol).expect("baseline lowering");
    let orig_loss = lowered.certificate.total_loss_bound();
    let orig_branches = lowered.law.branches().to_vec();
    let (w_lo, w_hi) = lowered.law.band_rad_s();

    let mut min_degradation_ratio = f64::INFINITY;
    for drop_idx in 0..orig_branches.len() {
        let mut kept = orig_branches.clone();
        let _ = kept.remove(drop_idx);
        let mutated_law = TensorProny::new(*lowered.law.c_eq(), None, kept, (w_lo, w_hi))
            .expect("mutated law remains passive");
        let mutated_cert =
            certify_lowered_tensor(&target, &mutated_law, w_lo, w_hi, CERTIFICATE_GRID_POINTS)
                .expect("mutated certificate");
        let mutated_loss = mutated_cert.total_loss_bound();
        let ratio = mutated_loss / orig_loss.max(1e-30);
        min_degradation_ratio = min_degradation_ratio.min(ratio);
        assert!(
            mutated_loss > tol.loss_rel
                || mutated_cert.total_storage_bound() > tol.storage_rel
                || ratio >= 2.0,
            "dropping branch {drop_idx} did not break certificate or double error: orig={orig_loss:.6e}, mutated={mutated_loss:.6e}, ratio={ratio:.3}"
        );
    }

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g3_dropping_any_branch_breaks_the_certificate\",\"tier\":\"G3\",\"branches_tested\":{},\"min_degradation_ratio\":{:.6e},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        orig_branches.len(),
        min_degradation_ratio,
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g0_term_count_and_band_refusals() {
    let t0 = Instant::now();
    let target = complex_engineering_wood_fixture();
    let loose_tol = TensorLoweringTolerance {
        storage_rel: 0.05,
        loss_rel: 0.05,
    };

    // n_terms = 64 is accepted:
    let ok64 = lower_tensor_to_prony(&target, 20.0, 20_000.0, MAX_PRONY_BRANCHES, loose_tol)
        .expect("n_terms = 64 must be accepted");
    assert!(ok64.law.branches().len() <= MAX_PRONY_BRANCHES);

    // n_terms = 0 and 65 are refused with Parameters:
    assert!(matches!(
        lower_tensor_to_prony(&target, 20.0, 20_000.0, 0, loose_tol),
        Err(ViscoError::Parameters { .. })
    ));
    assert!(matches!(
        lower_tensor_to_prony(&target, 20.0, 20_000.0, MAX_PRONY_BRANCHES + 1, loose_tol),
        Err(ViscoError::Parameters { .. })
    ));

    // Bad frequency bands are refused with Parameters:
    for (f_lo, f_hi) in [
        (0.0, 1000.0),
        (-20.0, 1000.0),
        (1000.0, 1000.0),
        (2000.0, 100.0),
        (f64::NAN, 1000.0),
        (20.0, f64::INFINITY),
    ] {
        assert!(
            matches!(
                lower_tensor_to_prony(&target, f_lo, f_hi, 8, loose_tol),
                Err(ViscoError::Parameters { .. })
            ),
            "band ({f_lo}, {f_hi}) must be refused"
        );
    }

    // Bad tolerances are refused with Parameters:
    assert!(matches!(
        lower_tensor_to_prony(
            &target,
            20.0,
            20_000.0,
            8,
            TensorLoweringTolerance {
                storage_rel: 0.0,
                loss_rel: 1e-2,
            }
        ),
        Err(ViscoError::Parameters { .. })
    ));

    // Unreachable tolerance (loss_rel = 1e-6 with 2 terms) -> FitTolerance:
    let tight_err = lower_tensor_to_prony(
        &target,
        20.0,
        20_000.0,
        2,
        TensorLoweringTolerance {
            storage_rel: 1e-2,
            loss_rel: 1e-6,
        },
    )
    .expect_err("unreachable loss tolerance with 2 terms must fail");
    assert!(
        matches!(tight_err, ViscoError::FitTolerance { .. }),
        "expected FitTolerance, got {tight_err:?}"
    );

    // Forced iteration cap of 1 on general solver -> FitDiverged:
    let div_err = lower_tensor_to_prony_with_max_iters(&target, 20.0, 20_000.0, 8, loose_tol, 1)
        .expect_err("max_iters = 1 must refuse with FitDiverged");
    assert!(
        matches!(div_err, ViscoError::FitDiverged { iterations: 1, .. }),
        "expected FitDiverged(iterations: 1), got {div_err:?}"
    );

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g0_term_count_and_band_refusals\",\"tier\":\"G0\",\"max_branches\":{},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        ok64.law.branches().len(),
        t0.elapsed().as_secs_f64()
    );
}

#[test]
fn g0_lowered_tensor_is_bitwise_deterministic() {
    let t0 = Instant::now();
    let tol = TensorLoweringTolerance {
        storage_rel: 1e-3,
        loss_rel: 1e-2,
    };
    let target = complex_engineering_wood_fixture();
    let a = lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol).expect("run a");
    let b = lower_tensor_to_prony(&target, 20.0, 20_000.0, 8, tol).expect("run b");

    assert_eq!(a.law.branches().len(), b.law.branches().len());
    for r in 0..6 {
        for c in 0..6 {
            assert_eq!(a.law.c_eq()[r][c].to_bits(), b.law.c_eq()[r][c].to_bits());
        }
    }
    for (ba, bb) in a.law.branches().iter().zip(b.law.branches()) {
        assert_eq!(ba.tau_s.to_bits(), bb.tau_s.to_bits());
        for r in 0..6 {
            for c in 0..6 {
                assert_eq!(ba.c_pa[r][c].to_bits(), bb.c_pa[r][c].to_bits());
            }
        }
    }
    assert_eq!(
        a.certificate.grid_sup_storage.to_bits(),
        b.certificate.grid_sup_storage.to_bits()
    );
    assert_eq!(
        a.certificate.grid_sup_loss.to_bits(),
        b.certificate.grid_sup_loss.to_bits()
    );
    assert_eq!(
        a.certificate.lipschitz_term.to_bits(),
        b.certificate.lipschitz_term.to_bits()
    );

    // Also verify KelvinFractional determinism:
    let k_target = kelvin_wood_fixture();
    let ka = lower_tensor_to_prony(&k_target, 20.0, 20_000.0, 8, tol).expect("kelvin a");
    let kb = lower_tensor_to_prony(&k_target, 20.0, 20_000.0, 8, tol).expect("kelvin b");
    for (ba, bb) in ka.law.branches().iter().zip(kb.law.branches()) {
        assert_eq!(ba.tau_s.to_bits(), bb.tau_s.to_bits());
        for r in 0..6 {
            for c in 0..6 {
                assert_eq!(ba.c_pa[r][c].to_bits(), bb.c_pa[r][c].to_bits());
            }
        }
    }
    let _ = PronyBranch::new(ka.law.branches()[0].c_pa, ka.law.branches()[0].tau_s).expect("psd");

    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2b\",\"test\":\"g0_lowered_tensor_is_bitwise_deterministic\",\"tier\":\"G0\",\"branches\":{},\"verdict\":\"pass\",\"elapsed_s\":{:.6}}}",
        a.law.branches().len(),
        t0.elapsed().as_secs_f64()
    );
}
