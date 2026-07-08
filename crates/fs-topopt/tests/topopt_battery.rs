//! fs-topopt battery (7tv.11 slice 1): G0 filter laws (linearity,
//! transpose adjointness, constant preservation), projection
//! monotonicity + endpoints + slope, FULL-CHAIN sensitivity
//! verification vs FD at multiple continuation stages (the acceptance
//! requirement — SIMP ∘ projection ∘ filter ∘ elasticity-solve), an
//! OC cantilever run with volume control and deterministic bitwise
//! replay (G5: a whole topo run replayable), and the golden hash.

use fs_adjoint::verify_gradient;
use fs_feec::kuhn_cube;
use fs_rand::StreamKey;
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, SimpParams, heaviside, heaviside_derivative,
    optimality_criteria,
};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-topopt\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 51,
        kernel: 0x0770,
        tile,
    }
    .stream();
    (0..n).map(|_| s.next_f64()).collect()
}

/// Cantilever fixture: unit-cube kuhn mesh, x = 0 face fully fixed,
/// downward tip load along the x = 1, z = 0 edge.
fn cantilever(m: usize) -> (fs_rep_mesh::TetComplex, Vec<[f64; 3]>, Vec<f64>, Vec<f64>) {
    let (complex, positions) = kuhn_cube(m);
    // Exact grid coordinates are intentional (kuhn positions are
    // rational multiples of 1/m).
    let fixed = |p: [f64; 3]| p[0].to_bits() == 0.0f64.to_bits();
    let el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &fixed);
    let mut force = vec![0.0f64; el.n()];
    for (v, &p) in positions.iter().enumerate() {
        if p[0].to_bits() == 1.0f64.to_bits() && p[2].to_bits() == 0.0f64.to_bits() {
            force[3 * v + 2] = -1.0;
        }
    }
    let geo = fs_feec::element_geometry(&complex, &positions);
    let vol: Vec<f64> = geo.vol_signed.iter().map(|v| v.abs()).collect();
    (complex, positions, force, vol)
}

#[test]
fn filter_g0_laws() {
    let (complex, positions) = kuhn_cube(2);
    let nc = complex.tets.len();
    let filter = DensityFilter::new(&complex, &positions, 0.15);
    // Linearity: F(a·x + y) = a·F(x) + F(y) to solver tolerance.
    let x = rand_vec(nc, 1);
    let y = rand_vec(nc, 2);
    let a = 1.7f64;
    let combo: Vec<f64> = x.iter().zip(&y).map(|(xi, yi)| a * xi + yi).collect();
    let f_combo = filter.apply(&combo);
    let fx = filter.apply(&x);
    let fy = filter.apply(&y);
    let worst = f_combo
        .iter()
        .zip(fx.iter().zip(&fy))
        .map(|(fc, (fxi, fyi))| (fc - (a * fxi + fyi)).abs())
        .fold(0.0f64, f64::max);
    assert!(worst < 1e-9, "filter not linear: {worst:.3e}");
    // Transpose adjointness: ⟨F·x, w⟩ = ⟨x, Fᵀ·w⟩.
    let w = rand_vec(nc, 3);
    let lhs: f64 = fx.iter().zip(&w).map(|(p, q)| p * q).sum();
    let ft_w = filter.apply_transpose(&w);
    let rhs: f64 = x.iter().zip(&ft_w).map(|(p, q)| p * q).sum();
    let rel = (lhs - rhs).abs() / lhs.abs().max(1e-30);
    assert!(
        rel < 1e-9,
        "filter transpose broken: {lhs:.12e} vs {rhs:.12e}"
    );
    // Constants preserved (natural BCs — no boundary droop).
    let ones = vec![0.7f64; nc];
    let f_ones = filter.apply(&ones);
    let dev = f_ones
        .iter()
        .map(|v| (v - 0.7).abs())
        .fold(0.0f64, f64::max);
    assert!(dev < 1e-9, "filter must preserve constants: dev {dev:.3e}");
    log(
        "filter-g0",
        "pass",
        &format!("linearity {worst:.1e}, adjoint {rel:.1e}, const {dev:.1e}"),
    );
}

#[test]
fn projection_g0_laws() {
    let (beta, eta) = (4.0f64, 0.5f64);
    // Endpoints exact.
    assert!((heaviside(0.0, beta, eta)).abs() < 1e-15);
    assert!((heaviside(1.0, beta, eta) - 1.0).abs() < 1e-15);
    // Monotonicity on a fine sweep.
    let mut prev = heaviside(0.0, beta, eta);
    for k in 1..=100 {
        let r = f64::from(k) / 100.0;
        let cur = heaviside(r, beta, eta);
        assert!(cur >= prev, "projection not monotone at {r}");
        prev = cur;
    }
    // Slope vs FD.
    for &r in &[0.2f64, 0.5, 0.8] {
        let eps = 1e-6;
        let fd = (heaviside(r + eps, beta, eta) - heaviside(r - eps, beta, eta)) / (2.0 * eps);
        let an = heaviside_derivative(r, beta, eta);
        assert!(
            (fd - an).abs() < 1e-8,
            "slope mismatch at {r}: {fd} vs {an}"
        );
    }
    log("projection-g0", "pass", "endpoints, monotone, slope");
}

#[test]
fn full_chain_sensitivity_at_continuation_stages() {
    // The acceptance gate: dc/dρ through SIMP ∘ projection ∘ filter ∘
    // solve, FD-verified at MULTIPLE continuation stages (early,
    // mid, sharp).
    let (complex, positions, force, _vol) = cantilever(2);
    let nc = complex.tets.len();
    let rho0: Vec<f64> = rand_vec(nc, 10).iter().map(|v| 0.3 + 0.5 * v).collect();
    for (stage, (penal, beta)) in [(1.0f64, 1.0f64), (3.0, 2.0), (3.0, 8.0)]
        .iter()
        .enumerate()
    {
        let pipeline = DesignPipeline {
            filter: DensityFilter::new(&complex, &positions, 0.12),
            params: SimpParams {
                e_min: 1e-6,
                penal: *penal,
                beta: *beta,
                eta: 0.5,
            },
        };
        let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
            p[0].to_bits() == 0.0f64.to_bits()
        });
        let (_, _, grad) = pipeline.compliance_and_gradient(&mut el, &rho0, &force);
        let j = |rho: &[f64]| -> f64 {
            let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
                p[0].to_bits() == 0.0f64.to_bits()
            });
            pipeline.compliance_and_gradient(&mut el2, rho, &force).0
        };
        let dirs: Vec<Vec<f64>> = (0..2)
            .map(|k| rand_vec(nc, 20 + u32::try_from(stage).expect("small") * 10 + k))
            .collect();
        let verdict = verify_gradient(&j, &rho0, &grad, &dirs, 1e-6, 2e-4);
        assert!(
            verdict.pass,
            "stage {stage} (p={penal}, beta={beta}): sensitivity failed FD: {:.3e}",
            verdict.max_rel_err
        );
        log(
            "chain-sensitivity",
            "pass",
            &format!("p={penal} beta={beta} rel={:.2e}", verdict.max_rel_err),
        );
    }
}

#[test]
fn oc_cantilever_descends_and_replays() {
    let (complex, positions, force, vol) = cantilever(3);
    let nc = complex.tets.len();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
        },
    };
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let vol_frac = 0.4;
    let rho0 = vec![vol_frac; nc];
    let rep = optimality_criteria(&pipeline, &mut el, &force, &rho0, &vol, vol_frac, 0.2, 12);
    let c0 = rep.compliance[0];
    let c_final = *rep.compliance.last().expect("trace");
    assert!(
        c_final < 0.8 * c0,
        "OC failed to improve compliance: {c0} -> {c_final}"
    );
    let v_final = *rep.volume.last().expect("trace");
    assert!(
        (v_final - vol_frac).abs() < 0.03,
        "volume constraint missed: {v_final} vs {vol_frac}"
    );
    // Design is differentiated (not uniform gray): spread must grow.
    let spread = rep
        .rho
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
    assert!(
        spread.1 - spread.0 > 0.5,
        "design stayed gray: range {spread:?}"
    );
    // G5: a whole run replays bitwise.
    let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let rep2 = optimality_criteria(&pipeline, &mut el2, &force, &rho0, &vol, vol_frac, 0.2, 12);
    assert!(
        rep.rho
            .iter()
            .zip(&rep2.rho)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "topo run not replayable"
    );
    log(
        "oc-cantilever",
        "pass",
        &format!(
            "c: {c0:.4e} -> {c_final:.4e}, vol {v_final:.3}, range [{:.2},{:.2}], change {:.3}",
            spread.0, spread.1, rep.final_change
        ),
    );
}

const GOLDEN_HASH: u64 = 0x772a_2f8c_a720_dd64; // recorded at 7tv.11 slice 1, frozen

#[test]
fn topopt_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions, force, vol) = cantilever(2);
    let nc = complex.tets.len();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
    };
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    // Pipeline forward + gradient fingerprints.
    let rho = rand_vec(nc, 40);
    let (rho_tilde, rho_bar, moduli) = pipeline.forward(&rho);
    for v in rho_tilde.iter().step_by(5).chain(rho_bar.iter().step_by(7)) {
        feed(*v);
    }
    for v in moduli.iter().step_by(3) {
        feed(*v);
    }
    let (c, _, grad) = pipeline.compliance_and_gradient(&mut el, &rho, &force);
    feed(c);
    for v in grad.iter().step_by(3) {
        feed(*v);
    }
    // Short OC fingerprint.
    let rep = optimality_criteria(
        &pipeline,
        &mut el,
        &force,
        &vec![0.4; nc],
        &vol,
        0.4,
        0.2,
        3,
    );
    for v in rep.rho.iter().step_by(5) {
        feed(*v);
    }
    log("topopt-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "topopt bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn robust_three_field_ordering_and_sensitivity() {
    use fs_topopt::RobustPipeline;
    let (complex, positions, force, _vol) = cantilever(2);
    let nc = complex.tets.len();
    let pipeline = RobustPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.12),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 4.0,
            eta: 0.5,
        },
        eta_offset: 0.15,
    };
    // Pointwise ordering: eroded ≤ nominal ≤ dilated for random designs.
    for tile in 60..63u32 {
        let rho = rand_vec(nc, tile);
        let tf = pipeline.three_fields(&rho);
        for i in 0..nc {
            assert!(
                tf.eroded[i] <= tf.nominal[i] + 1e-14 && tf.nominal[i] <= tf.dilated[i] + 1e-14,
                "three-field ordering violated at cell {i}: {} / {} / {}",
                tf.eroded[i],
                tf.nominal[i],
                tf.dilated[i]
            );
        }
    }
    // Eroded-compliance sensitivity FD gate (the robust objective's
    // gradient through the eroded projection chain).
    let rho0: Vec<f64> = rand_vec(nc, 63).iter().map(|v| 0.3 + 0.5 * v).collect();
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let (_, grad) = pipeline.eroded_compliance_and_gradient(&mut el, &rho0, &force);
    let j = |rho: &[f64]| -> f64 {
        let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
            p[0].to_bits() == 0.0f64.to_bits()
        });
        pipeline
            .eroded_compliance_and_gradient(&mut el2, rho, &force)
            .0
    };
    let dirs: Vec<Vec<f64>> = (0..2).map(|k| rand_vec(nc, 70 + k)).collect();
    let verdict = verify_gradient(&j, &rho0, &grad, &dirs, 1e-6, 2e-4);
    assert!(
        verdict.pass,
        "robust sensitivity failed FD: {:.3e}",
        verdict.max_rel_err
    );
    log(
        "robust-fields",
        "pass",
        &format!("ordering ok, eroded-grad rel={:.2e}", verdict.max_rel_err),
    );
}

#[test]
fn robust_oc_improves_erosion_retention() {
    use fs_topopt::{RobustPipeline, robust_optimality_criteria};
    // THE minimum-length-scale claim, measured: the robust design must
    // survive erosion better than the slice-1 (non-robust) design —
    // features thinner than the erode band cannot carry the robust
    // load path, so vol(eroded)/vol(nominal) stays high.
    let (complex, positions, force, vol) = cantilever(3);
    let nc = complex.tets.len();
    let params = SimpParams {
        e_min: 1e-6,
        penal: 3.0,
        beta: 6.0,
        eta: 0.5,
    };
    let vol_frac = 0.4;
    let rho0 = vec![vol_frac; nc];
    // Non-robust baseline (slice-1 OC), audited with the SAME
    // three-field probe.
    let nominal_pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params,
    };
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let base = optimality_criteria(
        &nominal_pipeline,
        &mut el,
        &force,
        &rho0,
        &vol,
        vol_frac,
        0.2,
        25,
    );
    let probe = RobustPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params,
        eta_offset: 0.15,
    };
    let base_tf = probe.three_fields(&base.rho);
    let vf = |field: &[f64]| -> f64 {
        let total: f64 = vol.iter().sum();
        field.iter().zip(&vol).map(|(r, v)| r * v).sum::<f64>() / total
    };
    let base_retention = vf(&base_tf.eroded) / vf(&base_tf.nominal).max(1e-30);
    // Robust run.
    let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let rep = robust_optimality_criteria(&probe, &mut el2, &force, &rho0, &vol, vol_frac, 0.2, 25);
    // Eroded compliance descends.
    let c0 = rep.compliance_eroded[0];
    let c_final = *rep.compliance_eroded.last().expect("trace");
    assert!(
        c_final < 0.8 * c0,
        "robust OC failed to improve eroded compliance: {c0} -> {c_final}"
    );
    // Volumes ordered.
    let (ve, vn, vd) = rep.volumes;
    assert!(
        ve <= vn + 1e-12 && vn <= vd + 1e-12,
        "volumes disordered: {ve} {vn} {vd}"
    );
    // NOMINAL volume at target (the adapted-dilated-constraint
    // contract: the nominal design carries the budget).
    assert!((vn - vol_frac).abs() < 0.05, "nominal volume missed: {vn}");
    // The length-scale signal: robust retention ≥ baseline retention
    // (strictly better in practice; ≥ −1e-9 guards roundoff ties).
    assert!(
        rep.erosion_retention >= base_retention - 1e-9,
        "robust design must survive erosion at least as well: robust {:.3} vs base {:.3}",
        rep.erosion_retention,
        base_retention
    );
    log(
        "robust-oc",
        "pass",
        &format!(
            "eroded c {c0:.3e}->{c_final:.3e}, vols ({ve:.3},{vn:.3},{vd:.3}), retention {:.3} vs base {:.3}",
            rep.erosion_retention, base_retention
        ),
    );
}

const ROBUST_GOLDEN_HASH: u64 = 0x519a_41e3_466e_4b7d; // recorded at 7tv.11 slice 2, frozen

#[test]
fn robust_golden_hash() {
    use fs_topopt::{RobustPipeline, robust_optimality_criteria};
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions, force, vol) = cantilever(2);
    let nc = complex.tets.len();
    let pipeline = RobustPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
        eta_offset: 0.12,
    };
    let rho = rand_vec(nc, 80);
    let tf = pipeline.three_fields(&rho);
    for v in tf
        .eroded
        .iter()
        .step_by(5)
        .chain(tf.dilated.iter().step_by(7))
    {
        feed(*v);
    }
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let (c, grad) = pipeline.eroded_compliance_and_gradient(&mut el, &rho, &force);
    feed(c);
    for v in grad.iter().step_by(4) {
        feed(*v);
    }
    let rep = robust_optimality_criteria(
        &pipeline,
        &mut el,
        &force,
        &vec![0.4; nc],
        &vol,
        0.4,
        0.2,
        2,
    );
    for v in rep.rho.iter().step_by(5) {
        feed(*v);
    }
    log("robust-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, ROBUST_GOLDEN_HASH,
        "robust bits changed: {acc:#018x} vs {ROBUST_GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn eigenfrequency_gradient_fd_gates() {
    use fs_topopt::{eigenfrequency_objective, lowest_eigenpairs, mass_interp};
    let (complex, positions, _force, _vol) = cantilever(2);
    let nt = complex.tets.len();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.12),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
        },
    };
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let rho0: Vec<f64> = rand_vec(nt, 100).iter().map(|v| 0.4 + 0.4 * v).collect();
    // Smooth-min aggregate gradient vs FD through the WHOLE chain
    // (filter, projection, SIMP stiffness, interpolated mass, eigen).
    let (agg0, grad) = eigenfrequency_objective(&pipeline, &mut el, &rho0, 4, 40.0);
    assert!(agg0 > 0.0, "base eigenvalue must be positive: {agg0}");
    let j = |rho: &[f64]| -> f64 {
        let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
            p[0].to_bits() == 0.0f64.to_bits()
        });
        eigenfrequency_objective(&pipeline, &mut el2, rho, 4, 40.0).0
    };
    let dirs: Vec<Vec<f64>> = (0..2).map(|k| rand_vec(nt, 110 + k)).collect();
    let verdict = fs_adjoint::verify_gradient(&j, &rho0, &grad, &dirs, 1e-5, 5e-4);
    assert!(
        verdict.pass,
        "eigenfrequency gradient failed FD: {:.3e}",
        verdict.max_rel_err
    );
    // The mass-interpolation trap gate: a design with VOID cells must
    // not host spurious low modes — λ_min of a mostly-void design with
    // one solid load path stays well above zero.
    let mut rho_void = vec![0.03f64; nt];
    for (c, r) in rho_void.iter_mut().enumerate() {
        if c % 3 == 0 {
            *r = 1.0;
        }
    }
    let (_, rho_bar_v, moduli_v) = pipeline.forward(&rho_void);
    el.moduli = moduli_v;
    let mass_v: Vec<f64> = rho_bar_v.iter().map(|&r| mass_interp(r)).collect();
    let (lv, _, _) = lowest_eigenpairs(&el, &mass_v, 2);
    assert!(
        lv[0] > 1e-4,
        "spurious void mode: lambda_min {:.3e} (the rho^6 mass floor must prevent this)",
        lv[0]
    );
    log(
        "eigenfreq-grad",
        "pass",
        &format!(
            "agg {agg0:.4}, FD rel {:.2e}, void lambda_min {:.2e}",
            verdict.max_rel_err, lv[0]
        ),
    );
}

#[test]
fn eigenfrequency_clustered_and_improvement() {
    use fs_topopt::eigenfrequency_objective;
    let (complex, positions, _force, vol) = cantilever(2);
    let nt = complex.tets.len();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.12),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
        },
    };
    // The cantilever cube has two symmetric bending modes (y/z) that
    // sit CLOSE for symmetric designs — exactly the clustered case the
    // smooth-min handles. Verify the aggregate gradient near that
    // near-crossing via FD (the single-eigenvalue objective would be
    // nonsmooth here; the aggregate is smooth).
    let rho_sym = vec![0.55f64; nt];
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let (_, grad) = eigenfrequency_objective(&pipeline, &mut el, &rho_sym, 4, 60.0);
    let j = |rho: &[f64]| -> f64 {
        let mut el2 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
            p[0].to_bits() == 0.0f64.to_bits()
        });
        eigenfrequency_objective(&pipeline, &mut el2, rho, 4, 60.0).0
    };
    let dirs: Vec<Vec<f64>> = (0..2).map(|k| rand_vec(nt, 120 + k)).collect();
    let verdict = fs_adjoint::verify_gradient(&j, &rho_sym, &grad, &dirs, 1e-5, 5e-4);
    assert!(
        verdict.pass,
        "clustered aggregate gradient failed FD near the crossing: {:.3e}",
        verdict.max_rel_err
    );
    // Improvement demo: a few projected-gradient ascent steps at fixed
    // volume must RAISE the aggregate eigenvalue (measured).
    let vol_frac = 0.5f64;
    let total: f64 = vol.iter().sum();
    let project = |rho: &mut Vec<f64>| {
        // Scale toward the volume budget then clamp (simple exact-ish
        // projection for the demo; OC/AL are the production drivers).
        for _ in 0..40 {
            let v: f64 = rho.iter().zip(&vol).map(|(r, w)| r * w).sum::<f64>() / total;
            let s = vol_frac / v.max(1e-12);
            for r in rho.iter_mut() {
                *r = (*r * s).clamp(1e-3, 1.0);
            }
            if (v - vol_frac).abs() < 1e-6 {
                break;
            }
        }
    };
    let mut rho = vec![vol_frac; nt];
    project(&mut rho);
    let mut el3 = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let (l0, _) = eigenfrequency_objective(&pipeline, &mut el3, &rho, 4, 60.0);
    for _ in 0..8 {
        let (_, g) = eigenfrequency_objective(&pipeline, &mut el3, &rho, 4, 60.0);
        let gnorm = g.iter().map(|x| x.abs()).fold(0.0f64, f64::max).max(1e-30);
        for (r, gi) in rho.iter_mut().zip(&g) {
            *r = (*r + 0.05 * gi / gnorm).clamp(1e-3, 1.0);
        }
        project(&mut rho);
    }
    let (l_final, _) = eigenfrequency_objective(&pipeline, &mut el3, &rho, 4, 60.0);
    assert!(
        l_final > 1.05 * l0,
        "ascent must raise the aggregate eigenvalue: {l0:.5} -> {l_final:.5}"
    );
    log(
        "eigenfreq-opt",
        "pass",
        &format!(
            "clustered FD rel {:.2e}, lambda_agg {l0:.5} -> {l_final:.5} (+{:.0}%)",
            verdict.max_rel_err,
            100.0 * (l_final / l0 - 1.0)
        ),
    );
}

const EIGEN_GOLDEN_HASH: u64 = 0xbb7e_5ad3_851a_2bf1; // recorded at mdx2 slice A, frozen

#[test]
fn eigenfreq_golden_hash() {
    use fs_topopt::{eigenfrequency_objective, lowest_eigenpairs, mass_interp};
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions, _force, _vol) = cantilever(2);
    let nt = complex.tets.len();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
    };
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p: [f64; 3]| {
        p[0].to_bits() == 0.0f64.to_bits()
    });
    let rho = rand_vec(nt, 130)
        .iter()
        .map(|v| 0.4 + 0.4 * v)
        .collect::<Vec<_>>();
    let (_, rho_bar, moduli) = pipeline.forward(&rho);
    el.moduli = moduli;
    let mass: Vec<f64> = rho_bar.iter().map(|&r| mass_interp(r)).collect();
    let (lambdas, _, _) = lowest_eigenpairs(&el, &mass, 4);
    for l in &lambdas {
        feed(*l);
    }
    let (agg, grad) = eigenfrequency_objective(&pipeline, &mut el, &rho, 3, 50.0);
    feed(agg);
    for v in grad.iter().step_by(9) {
        feed(*v);
    }
    log("eigenfreq-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, EIGEN_GOLDEN_HASH,
        "eigenfreq bits changed: {acc:#018x} vs {EIGEN_GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
