//! Integration and Gauntlet verification suite for passive full-tensor Prony
//! viscoelasticity (`fs_material::visco::tensor::TensorProny`, IRP M0.2-W2a).

use fs_la::eigen::jacobi_eigh;
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::tensor::{mandel_rotation, rotate, rotation};
use fs_material::visco::tensor::{MAX_PRONY_BRANCHES, Mandel6, PronyBranch, TensorProny};
use fs_material::visco::{GeneralizedMaxwell, ViscoError};
use fs_material::OrthotropicElastic;
use fs_math::c64::C64;
use fs_math::det;
use fs_rand::StreamKey;
use std::time::Instant;

fn identity6(scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { scale } else { 0.0 }))
}

fn frob6(c: &Mandel6) -> f64 {
    det::sqrt(c.iter().flatten().map(|v| v * v).sum::<f64>())
}

fn diff_frob6(a: &Mandel6, b: &Mandel6) -> f64 {
    let mut sum = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            let d = a[i][j] - b[i][j];
            sum += d * d;
        }
    }
    det::sqrt(sum)
}

fn golub_welsch_laguerre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut j = vec![0.0f64; n * n];
    for k in 0..n {
        j[k * n + k] = (2 * k + 1) as f64;
        if k + 1 < n {
            let off = (k + 1) as f64;
            j[k * n + (k + 1)] = off;
            j[(k + 1) * n + k] = off;
        }
    }
    let (nodes, evecs) = jacobi_eigh(&j, n);
    let weights: Vec<f64> = (0..n).map(|m| evecs[m] * evecs[m]).collect();
    (nodes, weights)
}

#[test]
fn g0_admission_requires_pd_equilibrium_psd_branches_positive_tau() {
    let t0 = Instant::now();
    let eq = identity6(10.0e9);
    let valid_branch = PronyBranch::new(identity6(1.0e8), 0.01).expect("valid branch");

    // 0 branches with KV absent is an elastic law and accepted.
    assert!(TensorProny::new(eq, None, vec![], (1.0, 1.0e5)).is_ok());

    // Branch cap: 64 accepted, 65 refused.
    assert_eq!(MAX_PRONY_BRANCHES, 64);
    let branches_64 = vec![valid_branch; 64];
    assert!(TensorProny::new(eq, None, branches_64, (1.0, 1.0e5)).is_ok());
    let branches_65 = vec![valid_branch; 65];
    assert!(matches!(
        TensorProny::new(eq, None, branches_65, (1.0, 1.0e5)),
        Err(ViscoError::Parameters { .. })
    ));

    // Indefinite C_eq refused.
    let mut indef_eq = eq;
    indef_eq[2][2] = -1.0e6;
    assert!(TensorProny::new(indef_eq, None, vec![], (1.0, 1.0e5)).is_err());
    let mut semidef_eq = eq;
    semidef_eq[0][1] = 10.0e9;
    semidef_eq[1][0] = 10.0e9;
    assert!(TensorProny::new(semidef_eq, None, vec![], (1.0, 1.0e5)).is_err());

    // Branch C_j with eigenvalue -1e-12 * ||C_j|| refused.
    let q = mandel_rotation(rotation(
        [1.0 / det::sqrt(3.0), 1.0 / det::sqrt(3.0), 1.0 / det::sqrt(3.0)],
        0.73,
    ));
    let mut diag_eigs = [1.0e9, 2.0e9, 3.0e9, 1.5e9, 0.5e9, 0.0];
    let norm_diag = det::sqrt(diag_eigs.iter().map(|v| v * v).sum::<f64>());
    diag_eigs[5] = -1.0e-12 * norm_diag;
    let mut neg_branch_c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut s = 0.0;
            for k in 0..6 {
                s += q[i][k] * diag_eigs[k] * q[j][k];
            }
            neg_branch_c[i][j] = s;
            neg_branch_c[j][i] = s;
        }
    }
    assert!(PronyBranch::new(neg_branch_c, 0.01).is_err());
    assert!(
        TensorProny::new(
            eq,
            None,
            vec![PronyBranch::unvalidated(neg_branch_c, 0.01)],
            (1.0, 1.0e5)
        )
        .is_err()
    );

    // tau <= 0 refused.
    for bad_tau in [0.0, -1.0e-3, f64::NAN, f64::INFINITY] {
        assert!(
            TensorProny::new(
                eq,
                None,
                vec![PronyBranch::unvalidated(identity6(1.0e8), bad_tau)],
                (1.0, 1.0e5)
            )
            .is_err()
        );
    }

    // Non-finite entries refused.
    let mut nan_eq = eq;
    nan_eq[0][0] = f64::NAN;
    assert!(TensorProny::new(nan_eq, None, vec![], (1.0, 1.0e5)).is_err());
    let mut inf_branch = identity6(1.0e8);
    inf_branch[3][3] = f64::INFINITY;
    assert!(
        TensorProny::new(
            eq,
            None,
            vec![PronyBranch::unvalidated(inf_branch, 0.01)],
            (1.0, 1.0e5)
        )
        .is_err()
    );

    // Asymmetric entry refused.
    let mut asym_eq = eq;
    asym_eq[0][1] = 1.0e8;
    assert!(TensorProny::new(asym_eq, None, vec![], (1.0, 1.0e5)).is_err());
    let mut asym_branch = identity6(1.0e8);
    asym_branch[1][4] = 1.0e6;
    assert!(
        TensorProny::new(
            eq,
            None,
            vec![PronyBranch::unvalidated(asym_branch, 0.01)],
            (1.0, 1.0e5)
        )
        .is_err()
    );

    // Band with lo > hi or negative lo refused.
    assert!(TensorProny::new(eq, None, vec![], (1.0e4, 1.0e2)).is_err());
    assert!(TensorProny::new(eq, None, vec![], (-1.0, 1.0e2)).is_err());

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g0_admission_requires_pd_equilibrium_psd_branches_positive_tau\",\"tier\":\"G0\",\"inputs\":{{\"law\":\"random-psd\",\"branches\":64,\"seed\":\"0x4d3032573261\"}},\"provenance\":\"synthetic fixture (admission boundary)\",\"rungs\":{{\"branch_cap\":[0,64,65]}},\"observed\":{{\"accepted_64\":true,\"refused_65\":true,\"refused_neg_eig_1e12\":true}},\"tolerance\":1e-12,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g0_random_multiaxial_histories_never_create_energy() {
    let t0 = Instant::now();
    let seed: u64 = 0x4d30_3257_3261;
    let histories = 256u32;
    let steps = 512usize;
    let mut min_step_heat = f64::INFINITY;
    let mut max_ledger_rel = 0.0f64;

    for history in 0..histories {
        let mut rng = StreamKey {
            seed,
            kernel: 1,
            tile: history,
        }
        .stream();

        // Build SPD equilibrium C_eq = 5e9 I + B0 B0ᵀ.
        let mut b0 = [[0.0f64; 6]; 6];
        for r in 0..6 {
            for c in 0..6 {
                b0[r][c] = 1.0e4 * (rng.next_f64() - 0.5);
            }
        }
        let mut c_eq = identity6(5.0e9);
        for i in 0..6 {
            for j in 0..=i {
                let mut s = c_eq[i][j];
                for k in 0..6 {
                    s += b0[i][k] * b0[j][k];
                }
                c_eq[i][j] = s;
                c_eq[j][i] = s;
            }
        }

        let n_branches = 1 + (history as usize % 6);
        let mut branches = Vec::with_capacity(n_branches);
        let mut tau_min = f64::INFINITY;
        for b_idx in 0..n_branches {
            let rank = 1 + ((history as usize + b_idx) % 6);
            let mut g_mat = [[0.0f64; 6]; 6];
            for r in 0..6 {
                for k in 0..rank {
                    g_mat[r][k] = 2.0e4 * (rng.next_f64() - 0.5);
                }
            }
            let mut c_pa = [[0.0f64; 6]; 6];
            for i in 0..6 {
                for j in 0..=i {
                    let mut s = 0.0;
                    for k in 0..rank {
                        s += g_mat[i][k] * g_mat[j][k];
                    }
                    c_pa[i][j] = s;
                    c_pa[j][i] = s;
                }
            }
            let log_tau = -4.0 + 4.0 * rng.next_f64();
            let tau_s = det::exp(log_tau * core::f64::consts::LN_10);
            tau_min = tau_min.min(tau_s);
            branches.push(PronyBranch::unvalidated(c_pa, tau_s));
        }

        let law = TensorProny::new(c_eq, None, branches, (0.0, 1.0e9))
            .expect("random PSD branches must be admitted");
        let mut state = law.state();
        let mut strain = [0.0f64; 6];

        for _ in 0..steps {
            for comp in &mut strain {
                *comp += 1.0e-4 * (rng.next_f64() - 0.5);
            }
            // Log-uniform dt in [1e-4, 1e2] * tau_min.
            let log_mult = -4.0 + 6.0 * rng.next_f64();
            let dt = tau_min * det::exp(log_mult * core::f64::consts::LN_10);
            let prev_heat = state.dissipated;
            let _ = law.step(&mut state, strain, dt).expect("step must succeed");
            let step_heat = state.dissipated - prev_heat;
            assert!(
                step_heat >= 0.0,
                "history {history}: step heat {step_heat:.3e} < 0"
            );
            assert!(
                state.dissipated >= prev_heat,
                "history {history}: cumulative heat decreased"
            );
            min_step_heat = min_step_heat.min(step_heat);
            let denom = state.work.abs().max(state.free_energy()).max(1.0e-30);
            let ledger_rel = (state.work - state.free_energy() - state.dissipated).abs() / denom;
            max_ledger_rel = max_ledger_rel.max(ledger_rel);
            assert!(
                ledger_rel <= 1.0e-12,
                "history {history}: ledger closure relative error {ledger_rel:.3e} > 1e-12"
            );
        }
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g0_random_multiaxial_histories_never_create_energy\",\"tier\":\"G0\",\"inputs\":{{\"law\":\"random-psd\",\"branches\":6,\"seed\":\"0x4d3032573261\"}},\"provenance\":\"synthetic fixture (Philox multiaxial strain histories)\",\"rungs\":{{\"histories\":{histories},\"steps\":{steps}}},\"observed\":{{\"min_step_heat\":{min_step_heat:.6e},\"max_ledger_rel\":{max_ledger_rel:.6e}}},\"tolerance\":1e-12,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g1_isotropic_proportional_law_reduces_to_scalar_maxwell_bitwise() {
    let t0 = Instant::now();
    let e_inf = 10.8e9;
    let nu = 0.28;
    let e_terms = vec![(1.8e9, 0.0025), (0.95e9, 0.018), (0.45e9, 0.14)];
    let e_law = GeneralizedMaxwell::new(e_inf, e_terms.clone()).unwrap();

    let k_scale = 1.0 / (3.0 * (1.0 - 2.0 * nu));
    let g_scale = 1.0 / (2.0 * (1.0 + nu));
    let k_law = GeneralizedMaxwell::new(
        e_inf * k_scale,
        e_terms.iter().map(|&(e, tau)| (e * k_scale, tau)).collect(),
    )
    .unwrap();
    let g_law = GeneralizedMaxwell::new(
        e_inf * g_scale,
        e_terms.iter().map(|&(e, tau)| (e * g_scale, tau)).collect(),
    )
    .unwrap();
    let k3_scalar = GeneralizedMaxwell::new(
        3.0 * k_law.e_inf,
        k_law.terms.iter().map(|&(k, tau)| (3.0 * k, tau)).collect(),
    )
    .unwrap();
    let g2_scalar = GeneralizedMaxwell::new(
        2.0 * g_law.e_inf,
        g_law.terms.iter().map(|&(g, tau)| (2.0 * g, tau)).collect(),
    )
    .unwrap();

    let tensor_law = TensorProny::isotropic(k_law, g_law).unwrap();

    // 1. Pure shear history eps12(t) equals GeneralizedMaxwell{2G...}::step bitwise.
    let mut s_tensor_shear = tensor_law.state();
    let mut s_scalar_shear = g2_scalar.state();
    for step in 0..64 {
        let eps12 = 1.2e-3 * det::sin(0.23 * step as f64) - 0.4e-3 * det::cos(0.07 * step as f64);
        let dt = if step == 0 {
            0.0
        } else {
            1.0e-4 * (1 + (step % 7)) as f64
        };
        let sig_tensor = tensor_law
            .step(&mut s_tensor_shear, [0.0, 0.0, 0.0, eps12, 0.0, 0.0], dt)
            .unwrap();
        let sig_scalar = g2_scalar.step(&mut s_scalar_shear, eps12, dt);
        assert_eq!(
            sig_tensor[3].to_bits(),
            sig_scalar.to_bits(),
            "pure shear step {step} must match GeneralizedMaxwell{{2G}} bitwise"
        );
    }

    // 2. Pure volumetric history [eps_v(t), eps_v(t), eps_v(t), 0, 0, 0] equals GeneralizedMaxwell{3K...}::step bitwise.
    let mut s_tensor_vol = tensor_law.state();
    let mut s_scalar_vol = k3_scalar.state();
    for step in 0..64 {
        let eps_v = 0.8e-3 * det::sin(0.19 * step as f64) + 0.3e-3 * det::cos(0.11 * step as f64);
        let dt = if step == 0 {
            0.0
        } else {
            2.0e-4 * (1 + (step % 5)) as f64
        };
        let sig_tensor = tensor_law
            .step(&mut s_tensor_vol, [eps_v, eps_v, eps_v, 0.0, 0.0, 0.0], dt)
            .unwrap();
        let sig_scalar = k3_scalar.step(&mut s_scalar_vol, eps_v, dt);
        assert_eq!(
            sig_tensor[0].to_bits(),
            sig_scalar.to_bits(),
            "pure volumetric step {step} must match GeneralizedMaxwell{{3K}} bitwise"
        );
        assert_eq!(sig_tensor[1].to_bits(), sig_scalar.to_bits());
        assert_eq!(sig_tensor[2].to_bits(), sig_scalar.to_bits());
    }

    // 3. Uniaxial-stress E*(omega) equals scalar GeneralizedMaxwell::modulus within 1e-14 relative on 64 log-spaced omega.
    let mut max_uniaxial_rel = 0.0f64;
    for idx in 0..64 {
        let frac = idx as f64 / 63.0;
        let omega = det::exp((-2.0 + 7.0 * frac) * core::f64::consts::LN_10);
        let (c_re, c_im) = tensor_law.complex_modulus(omega).unwrap();
        let c00 = C64::new(c_re[0][0], c_im[0][0]);
        let c01 = C64::new(c_re[0][1], c_im[0][1]);
        // Schur complement of normal block for uniaxial stress along x: E* = C00 - 2 C01² / (C00 + C01).
        let e_star = c00 - (c01 * c01).scale(2.0) / (c00 + c01);
        let (ep_ref, epp_ref) = e_law.modulus(omega);
        let e_ref = C64::new(ep_ref, epp_ref);
        let rel = (e_star - e_ref).abs() / e_ref.abs();
        max_uniaxial_rel = max_uniaxial_rel.max(rel);
        assert!(
            rel <= 1.0e-14,
            "omega={omega:.3e}: uniaxial E* relative error {rel:.3e} > 1e-14"
        );
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g1_isotropic_proportional_law_reduces_to_scalar_maxwell_bitwise\",\"tier\":\"G1\",\"inputs\":{{\"law\":\"isotropic\",\"branches\":6,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (isotropic proportional Zener/Maxwell)\",\"rungs\":{{\"steps\":64,\"frequencies\":64}},\"observed\":{{\"shear_ulp_diff\":0,\"volumetric_ulp_diff\":0,\"max_uniaxial_rel\":{max_uniaxial_rel:.6e}}},\"tolerance\":1e-14,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g1_kelvin_constructor_matches_general_form() {
    let t0 = Instant::now();
    let ortho = OrthotropicElastic::new(
        [12.0e9, 1.1e9, 0.7e9],
        [0.37, 0.45, 0.42],
        [0.75e9, 0.09e9, 0.68e9],
        0.02,
    )
    .unwrap();
    let c_eq = ortho.stiffness();
    let flat: Vec<f64> = c_eq.iter().flatten().copied().collect();
    let (evals, evecs) = jacobi_eigh(&flat, 6);

    let taus = [0.003, 0.025];
    let weights = [
        [0.08, 0.04],
        [0.12, 0.06],
        [0.10, 0.05],
        [0.15, 0.07],
        [0.09, 0.03],
        [0.11, 0.05],
    ];
    let per_mode: [GeneralizedMaxwell; 6] = core::array::from_fn(|k| {
        GeneralizedMaxwell::new(
            evals[k],
            vec![
                (evals[k] * weights[k][0], taus[0]),
                (evals[k] * weights[k][1], taus[1]),
            ],
        )
        .unwrap()
    });

    let kelvin_law = TensorProny::kelvin(c_eq, per_mode.clone()).unwrap();

    // Explicit general form constructed from rank-1 projector branches.
    let mut explicit_branches = Vec::new();
    for (k, mode) in per_mode.iter().enumerate() {
        for &(e_jk, tau_s) in &mode.terms {
            let mut c_pa = [[0.0f64; 6]; 6];
            for i in 0..6 {
                for j in 0..=i {
                    let v = e_jk * evecs[i * 6 + k] * evecs[j * 6 + k];
                    c_pa[i][j] = v;
                    c_pa[j][i] = v;
                }
            }
            explicit_branches.push(PronyBranch::unvalidated(c_pa, tau_s));
        }
    }
    let general_law = TensorProny::new(c_eq, None, explicit_branches, (0.0, f64::MAX)).unwrap();

    let mut max_rel = 0.0f64;
    for omega in [0.0, 10.0, 250.0, 4000.0] {
        let (re_k, im_k) = kelvin_law.complex_modulus(omega).unwrap();
        let (re_g, im_g) = general_law.complex_modulus(omega).unwrap();
        let rel_re = diff_frob6(&re_k, &re_g) / frob6(&re_g);
        let rel_im = if frob6(&im_g) > 0.0 {
            diff_frob6(&im_k, &im_g) / frob6(&im_g)
        } else {
            0.0
        };
        max_rel = max_rel.max(rel_re).max(rel_im);
    }
    assert!(max_rel <= 1.0e-15, "kelvin vs general max_rel={max_rel:.3e}");

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g1_kelvin_constructor_matches_general_form\",\"tier\":\"G1\",\"inputs\":{{\"law\":\"kelvin\",\"branches\":12,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (orthotropic Kelvin modes)\",\"rungs\":{{\"modes\":6,\"frequencies\":4}},\"observed\":{{\"max_rel\":{max_rel:.6e}}},\"tolerance\":1e-15,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g3_rotation_preserves_spectrum_and_dissipation() {
    let t0 = Instant::now();
    let ortho = OrthotropicElastic::new(
        [11.5e9, 0.95e9, 0.55e9],
        [0.36, 0.44, 0.41],
        [0.72e9, 0.08e9, 0.64e9],
        0.02,
    )
    .unwrap();
    let c_eq = ortho.stiffness();
    let mut b1 = identity6(0.15e9);
    b1[0][1] = 0.04e9;
    b1[1][0] = 0.04e9;
    b1[0][3] = -0.02e9;
    b1[3][0] = -0.02e9;
    let mut b2 = identity6(0.08e9);
    b2[1][2] = 0.025e9;
    b2[2][1] = 0.025e9;
    b2[2][4] = 0.015e9;
    b2[4][2] = 0.015e9;

    let base_law = TensorProny::new(
        c_eq,
        None,
        vec![
            PronyBranch::new(b1, 0.004).unwrap(),
            PronyBranch::new(b2, 0.025).unwrap(),
        ],
        (0.0, 1.0e7),
    )
    .unwrap();

    let base_eq_eigs = jacobi_eigh(&base_law.c_eq().iter().flatten().copied().collect::<Vec<_>>(), 6).0;
    let base_b_eigs: Vec<Vec<f64>> = base_law
        .branches()
        .iter()
        .map(|b| jacobi_eigh(&b.c_pa.iter().flatten().copied().collect::<Vec<_>>(), 6).0)
        .collect();

    // Reference strain history in Voigt tensor-shear coordinates:
    let sqrt2 = core::f64::consts::SQRT_2;
    let to_mandel = |v: &[f64; 6]| -> [f64; 6] {
        [v[0], v[1], v[2], sqrt2 * v[3], sqrt2 * v[4], sqrt2 * v[5]]
    };

    let mut base_state = base_law.state();
    let history_voigt: Vec<([f64; 6], f64)> = (0..48)
        .map(|k| {
            let s = (k + 1) as f64;
            let v = [
                1.0e-3 * det::sin(0.17 * s),
                -0.7e-3 * det::cos(0.13 * s),
                0.4e-3 * det::sin(0.29 * s),
                0.6e-3 * det::cos(0.19 * s),
                -0.5e-3 * det::sin(0.23 * s),
                0.3e-3 * det::cos(0.31 * s),
            ];
            (v, 0.0015)
        })
        .collect();

    for &(v_eps, dt) in &history_voigt {
        let _ = base_law
            .step(&mut base_state, to_mandel(&v_eps), dt)
            .unwrap();
    }

    let mut rng = StreamKey {
        seed: 0x4d30_3257_3261,
        kernel: 3,
        tile: 0,
    }
    .stream();
    let mut max_eig_rel = 0.0f64;
    let mut max_heat_rel = 0.0f64;

    for _rot_idx in 0..32 {
        let ax = rng.next_f64() - 0.5;
        let ay = rng.next_f64() - 0.5;
        let az = rng.next_f64() - 0.5;
        let norm = det::sqrt(ax * ax + ay * ay + az * az).max(1.0e-12);
        let axis = [ax / norm, ay / norm, az / norm];
        let angle = core::f64::consts::TAU * rng.next_f64();
        let q = rotation(axis, angle);

        let rot_law = base_law.rotated(q).unwrap();
        let rot_eq_eigs =
            jacobi_eigh(&rot_law.c_eq().iter().flatten().copied().collect::<Vec<_>>(), 6).0;
        let eq_scale = base_eq_eigs[5];
        for k in 0..6 {
            let rel = (rot_eq_eigs[k] - base_eq_eigs[k]).abs() / eq_scale;
            max_eig_rel = max_eig_rel.max(rel);
            assert!(rel <= 1.0e-13, "C_eq eigenvalue drift {rel:.3e} > 1e-13");
        }
        for (b_idx, branch) in rot_law.branches().iter().enumerate() {
            let eigs = jacobi_eigh(&branch.c_pa.iter().flatten().copied().collect::<Vec<_>>(), 6).0;
            let b_scale = base_b_eigs[b_idx][5];
            for k in 0..6 {
                let rel = (eigs[k] - base_b_eigs[b_idx][k]).abs() / b_scale;
                max_eig_rel = max_eig_rel.max(rel);
                assert!(rel <= 1.0e-13, "branch {b_idx} eigenvalue drift {rel:.3e} > 1e-13");
            }
        }

        let mut rot_state = rot_law.state();
        for &(v_eps, dt) in &history_voigt {
            let v_rot = rotate(&v_eps, &q);
            let _ = rot_law.step(&mut rot_state, to_mandel(&v_rot), dt).unwrap();
        }
        let heat_rel =
            (rot_state.dissipated - base_state.dissipated).abs() / base_state.dissipated;
        max_heat_rel = max_heat_rel.max(heat_rel);
        assert!(
            heat_rel <= 1.0e-12,
            "rotated heat relative drift {heat_rel:.3e} > 1e-12"
        );
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g3_rotation_preserves_spectrum_and_dissipation\",\"tier\":\"G3\",\"inputs\":{{\"law\":\"random-psd\",\"branches\":2,\"seed\":\"0x4d3032573261\"}},\"provenance\":\"synthetic fixture (32 proper Rodrigues rotations)\",\"rungs\":{{\"rotations\":32,\"steps\":48}},\"observed\":{{\"max_eig_rel\":{max_eig_rel:.6e},\"max_heat_rel\":{max_heat_rel:.6e}}},\"tolerance\":1e-12,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g1_kramers_kronig_consistency() {
    let t0 = Instant::now();
    let band = (1.0e2, 1.0e4);
    // Log grid spanning the band +/- 3 decades: [1e-1, 1e7] rad/s (8 decades, 4096 points).
    let grid_lo = band.0 * 1.0e-3;
    let grid_hi = band.1 * 1.0e3;
    let c_eq = identity6(10.0e9);
    let b1 = PronyBranch::new(identity6(1.2e9), 1.0 / 350.0).unwrap();
    let b2 = PronyBranch::new(identity6(0.8e9), 1.0 / 2800.0).unwrap();
    let law = TensorProny::new(c_eq, None, vec![b1, b2], (grid_lo, grid_hi)).unwrap();

    let n_grid = 4096usize;
    let ln_lo = det::ln(grid_lo);
    let ln_hi = det::ln(grid_hi);
    let d_ln_w = (ln_hi - ln_lo) / n_grid as f64;
    // Midpoint log grid avoidsomega == omega_0 pole directly while maintaining 2nd-order convergence.
    let omegas: Vec<f64> = (0..n_grid)
        .map(|m| det::exp(ln_lo + (m as f64 + 0.5) * d_ln_w))
        .collect();
    let delta_storage_00: Vec<f64> = omegas
        .iter()
        .map(|&w| law.complex_modulus(w).unwrap().0[0][0] - c_eq[0][0])
        .collect();

    let mut max_kk_rel = 0.0f64;
    let mut min_falsifier_rel = f64::INFINITY;
    let non_causal_hysteretic = 0.5e9;

    for probe_idx in 0..16 {
        let frac = (probe_idx as f64 + 0.5) / 16.0;
        let w0 = det::exp(det::ln(band.0) + frac * (det::ln(band.1) - det::ln(band.0)));
        let (re0, im0) = law.complex_modulus(w0).unwrap();
        let dc0 = re0[0][0] - c_eq[0][0];
        let exact_im = im0[0][0];

        let mut integral = 0.0f64;
        for (&w, &dc) in omegas.iter().zip(&delta_storage_00) {
            let denom = w * w - w0 * w0;
            if denom.abs() > 1.0e-12 * w0 * w0 {
                integral += ((dc - dc0) / denom) * w * d_ln_w;
            }
        }
        let kk_im = (2.0 * w0 / core::f64::consts::PI) * integral;
        let rel_err = (kk_im - exact_im).abs() / exact_im;
        max_kk_rel = max_kk_rel.max(rel_err);
        assert!(
            rel_err <= 1.0e-3,
            "w0={w0:.2e}: KK relative error {rel_err:.3e} > 1e-3"
        );

        let perturbed_im = exact_im + non_causal_hysteretic;
        let falsifier_rel = (kk_im - perturbed_im).abs() / perturbed_im;
        min_falsifier_rel = min_falsifier_rel.min(falsifier_rel);
    }
    assert!(
        min_falsifier_rel > 1.0e-1,
        "non-causal constant-loss perturbation must fail KK by > 1e-1 (got {min_falsifier_rel:.3e})"
    );

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g1_kramers_kronig_consistency\",\"tier\":\"G1\",\"inputs\":{{\"law\":\"isotropic\",\"branches\":2,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (4096-point log-grid Kramers-Kronig transform)\",\"rungs\":{{\"grid_points\":{n_grid},\"probes\":16}},\"observed\":{{\"max_kk_rel\":{max_kk_rel:.6e},\"min_falsifier_rel\":{min_falsifier_rel:.6e}}},\"tolerance\":1e-3,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g1_relaxation_function_fourier_matches_complex_modulus() {
    let t0 = Instant::now();
    let ortho = OrthotropicElastic::new(
        [11.0e9, 1.0e9, 0.6e9],
        [0.35, 0.42, 0.40],
        [0.7e9, 0.085e9, 0.62e9],
        0.02,
    )
    .unwrap();
    let c_eq = ortho.stiffness();
    let b1 = PronyBranch::new(identity6(0.8e9), 0.010).unwrap();
    let b2 = PronyBranch::new(identity6(0.4e9), 0.016).unwrap();
    let law = TensorProny::new(c_eq, None, vec![b1, b2], (0.1, 1.0e5)).unwrap();

    // 64-point Gauss-Laguerre quadrature via Golub-Welsch (`jacobi_eigh`) on
    // G(t) = C_eq + Σ C_j exp(-t / τ_j). With t = s * τ_max (τ_max = 0.016 s),
    // exp(s) * (G(s τ_max) - C_eq) = Σ C_j exp(-s (τ_max/τ_j - 1)) is bounded
    // for all s >= 0, and nodes with s > 36 have w_m < 2.3e-16 where C_eq + ΔG(t)
    // would otherwise suffer floating-point absorption.
    let (nodes, weights) = golub_welsch_laguerre(64);
    let tau_max = 0.016f64;
    let mut max_fourier_rel = 0.0f64;

    for omega in [10.0, 25.0, 60.0, 120.0] {
        let mut c_re_quad = c_eq;
        let mut c_im_quad = [[0.0f64; 6]; 6];
        for (&s, &w) in nodes.iter().zip(&weights) {
            if s > 36.0 {
                continue;
            }
            let t = s * tau_max;
            let g_t = law.relaxation_modulus(t).unwrap();
            let phase = -omega * t;
            let factor = w * tau_max * det::exp(s);
            let cos_p = det::cos(phase);
            let sin_p = det::sin(phase);
            for i in 0..6 {
                for j in 0..6 {
                    let dg = g_t[i][j] - c_eq[i][j];
                    c_re_quad[i][j] += factor * (-omega * dg * sin_p);
                    c_im_quad[i][j] += factor * (omega * dg * cos_p);
                }
            }
        }
        let (c_re_exact, c_im_exact) = law.complex_modulus(omega).unwrap();
        let rel_re = diff_frob6(&c_re_quad, &c_re_exact) / frob6(&c_re_exact);
        let rel_im = diff_frob6(&c_im_quad, &c_im_exact) / frob6(&c_im_exact);
        max_fourier_rel = max_fourier_rel.max(rel_re).max(rel_im);
        assert!(
            rel_re <= 1.0e-9 && rel_im <= 1.0e-9,
            "omega={omega}: Fourier-Laguerre error re={rel_re:.3e}, im={rel_im:.3e} > 1e-9"
        );
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g1_relaxation_function_fourier_matches_complex_modulus\",\"tier\":\"G1\",\"inputs\":{{\"law\":\"kelvin\",\"branches\":2,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (64-point Golub-Welsch Gauss-Laguerre quadrature of G(t))\",\"rungs\":{{\"quadrature_order\":64,\"frequencies\":4}},\"observed\":{{\"max_fourier_rel\":{max_fourier_rel:.6e}}},\"tolerance\":1e-9,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g0_plane_stress_schur_complement_is_passive_and_matches_elastic_limit() {
    let t0 = Instant::now();
    let e = [11.88e9, 0.92664e9, 0.51084e9];
    let nu = [0.372, 0.467, 0.435];
    let g = [0.76032e9, 0.03564e9, 0.72468e9];
    let ortho = OrthotropicElastic::new(e, nu, g, 0.01).unwrap();
    let c_eq = ortho.stiffness();

    // 1. Elastic limit (C_j -> 0) matches elastic plane-stress reduction of OrthotropicElastic within 1e-14.
    let elastic_prony = TensorProny::new(c_eq, None, vec![], (0.0, 1.0e6)).unwrap();
    let (q_re_12, q_im_12) = elastic_prony
        .plane_stress_complex(MaterialPlane::Axes12, 100.0)
        .unwrap();
    let nu21 = nu[0] * e[1] / e[0];
    let denom = 1.0 - nu[0] * nu21;
    let expected_12 = [
        e[0] / denom,
        nu[0] * e[1] / denom,
        0.0,
        nu[0] * e[1] / denom,
        e[1] / denom,
        0.0,
        0.0,
        0.0,
        g[0],
    ];
    let mut max_elastic_rel = 0.0f64;
    let norm_exp = det::sqrt(expected_12.iter().map(|v| v * v).sum::<f64>());
    for i in 0..9 {
        let rel = (q_re_12[i] - expected_12[i]).abs() / norm_exp;
        max_elastic_rel = max_elastic_rel.max(rel);
        assert_eq!(q_im_12[i], 0.0);
    }
    assert!(
        max_elastic_rel <= 1.0e-14,
        "elastic limit mismatch {max_elastic_rel:.3e} > 1e-14"
    );

    // 2. Frozen limit (tau -> inf, evaluated at omega -> inf vs C_eq + sum C_j) also matches within 1e-14.
    let b1 = PronyBranch::new(c_eq.map(|row| row.map(|v| 0.05 * v)), 1.0e12).unwrap();
    let frozen_law = TensorProny::new(c_eq, None, vec![b1], (0.0, 1.0e6)).unwrap();
    let (q_re_frozen, _) = frozen_law
        .plane_stress_complex(MaterialPlane::Axes12, 100.0)
        .unwrap();
    for i in 0..9 {
        let rel = (q_re_frozen[i] - 1.05 * expected_12[i]).abs() / (1.05 * norm_exp);
        max_elastic_rel = max_elastic_rel.max(rel);
        assert!(rel <= 1.0e-14);
    }

    // 3. Viscoelastic law with non-proportional branches: Im Q(omega) >= 0 at 64 frequencies.
    let mut b_aniso = identity6(0.12e9);
    b_aniso[0][2] = 0.04e9;
    b_aniso[2][0] = 0.04e9;
    let visco_law = TensorProny::new(
        c_eq,
        None,
        vec![
            PronyBranch::new(b_aniso, 0.002).unwrap(),
            PronyBranch::new(identity6(0.06e9), 0.03).unwrap(),
        ],
        (1.0, 1.0e6),
    )
    .unwrap();

    let mut min_im_eig = f64::INFINITY;
    for idx in 0..64 {
        let frac = idx as f64 / 63.0;
        let omega = det::exp((0.0 + 6.0 * frac) * core::f64::consts::LN_10).clamp(1.0, 1.0e6);
        let (re, im) = visco_law
            .plane_stress_complex(MaterialPlane::Axes12, omega)
            .unwrap();
        let re_eigs = jacobi_eigh(&re, 3).0;
        let im_eigs = jacobi_eigh(&im, 3).0;
        assert!(re_eigs[0] > 0.0);
        assert!(im_eigs[0] >= 0.0, "omega={omega}: Im eig {0:.3e} < 0", im_eigs[0]);
        min_im_eig = min_im_eig.min(im_eigs[0]);
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g0_plane_stress_schur_complement_is_passive_and_matches_elastic_limit\",\"tier\":\"G0\",\"inputs\":{{\"law\":\"random-psd\",\"branches\":2,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (orthotropic plane-stress Schur complement)\",\"rungs\":{{\"frequencies\":64}},\"observed\":{{\"max_elastic_rel\":{max_elastic_rel:.6e},\"min_im_eig\":{min_im_eig:.6e}}},\"tolerance\":1e-14,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}

#[test]
fn g0_out_of_band_refuses() {
    let t0 = Instant::now();
    let eq = identity6(8.0e9);
    let band = (100.0, 20_000.0);
    let law = TensorProny::new(
        eq,
        None,
        vec![PronyBranch::new(identity6(0.5e9), 0.01).unwrap()],
        band,
    )
    .unwrap();

    assert!(law.complex_modulus(100.0).is_ok());
    assert!(law.complex_modulus(20_000.0).is_ok());
    assert_eq!(
        law.complex_modulus(99.9),
        Err(ViscoError::OutOfBand { omega: 99.9, band })
    );
    assert_eq!(
        law.complex_modulus(20_000.1),
        Err(ViscoError::OutOfBand {
            omega: 20_000.1,
            band
        })
    );
    assert_eq!(
        law.plane_stress_complex(MaterialPlane::Axes12, 50.0),
        Err(ViscoError::OutOfBand { omega: 50.0, band })
    );

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W2a\",\"test\":\"g0_out_of_band_refuses\",\"tier\":\"G0\",\"inputs\":{{\"law\":\"isotropic\",\"branches\":1,\"seed\":\"none\"}},\"provenance\":\"synthetic fixture (band boundary refusal)\",\"rungs\":{{\"band_rad_s\":[100.0,20000.0]}},\"observed\":{{\"refused_below\":true,\"refused_above\":true}},\"tolerance\":0.0,\"verdict\":\"pass\",\"elapsed_s\":{elapsed:.6}}}"
    );
}
