//! Regression probe: the UNDAMPED dilated-target adaptation limit-cycled
//! (period 2: uniform 0.4 <-> 0.6, compliance 8.3e5 <-> 2.0e3 forever).
//! This probe mirrors the damped update and asserts the cycle is gone:
//! nominal volume settles near the budget and compliance stays at the
//! converged scale.
use fs_topopt::{DensityElasticity, DensityFilter, RobustPipeline, SimpParams};

#[test]
fn probe_robust_iterations() {
    let (complex, positions) = fs_feec::kuhn_cube(3);
    let fixed = |p: [f64; 3]| p[0].to_bits() == 0.0f64.to_bits();
    let mut el = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &fixed);
    let mut force = vec![0.0f64; el.n()];
    for (v, &p) in positions.iter().enumerate() {
        if p[0].to_bits() == 1.0f64.to_bits() && p[2].to_bits() == 0.0f64.to_bits() {
            force[3 * v + 2] = -1.0;
        }
    }
    let geo = fs_feec::element_geometry(&complex, &positions);
    let vol: Vec<f64> = geo.vol_signed.iter().map(|v| v.abs()).collect();
    let nc = complex.tets.len();
    let pipeline = RobustPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 6.0,
            eta: 0.5,
        },
        eta_offset: 0.15,
    };
    let mut rho = vec![0.4f64; nc];
    let mut target = 0.4f64;
    for it in 0..12 {
        let tf = pipeline.three_fields(&rho);
        let vf = |f: &[f64]| -> f64 {
            let t: f64 = vol.iter().sum();
            f.iter().zip(&vol).map(|(r, v)| r * v).sum::<f64>() / t
        };
        let (c, grad) = pipeline.eroded_compliance_and_gradient(&mut el, &rho, &force);
        let gmax = grad.iter().map(|g| g.abs()).fold(0.0f64, f64::max);
        let gmin = grad.iter().map(|g| g.abs()).fold(f64::INFINITY, f64::min);
        println!(
            "it {it}: c={c:.4e} ve={:.3} vn={:.3} vd={:.3} |g| in [{gmin:.2e},{gmax:.2e}] rho[0]={:.4}",
            vf(&tf.eroded),
            vf(&tf.nominal),
            vf(&tf.dilated),
            rho[0]
        );
        // One hand-rolled OC step mirroring robust_optimality_criteria.
        let vn_now = vf(&tf.nominal);
        let vd_now = vf(&tf.dilated);
        let instant = (0.4 * vd_now / vn_now).clamp(0.4, 1.0);
        target = 0.3f64.mul_add(instant - target, target);
        let sens: Vec<f64> = grad.iter().map(|g| (-g).max(1e-30)).collect();
        let (mut lo, mut hi) = (1e-12f64, 1e12f64);
        let mut cand = rho.clone();
        for _ in 0..80 {
            let lambda = fs_math::det::sqrt(lo * hi);
            for i in 0..nc {
                let scale = fs_math::det::sqrt(sens[i] / (lambda * vol[i]));
                cand[i] = (rho[i] * scale)
                    .clamp(rho[i] - 0.2, rho[i] + 0.2)
                    .clamp(1e-3, 1.0);
            }
            let d = pipeline.three_fields(&cand).dilated;
            if vf(&d) > target {
                lo = lambda;
            } else {
                hi = lambda;
            }
        }
        println!(
            "   target={target:.3} cand[0]={:.4} cand_max={:.4}",
            cand[0],
            cand.iter().copied().fold(0.0f64, f64::max)
        );
        rho = cand;
    }
    // No limit cycle: the final state is converged, not bouncing.
    let tf = pipeline.three_fields(&rho);
    let t: f64 = vol.iter().sum();
    let vn: f64 = tf.nominal.iter().zip(&vol).map(|(r, v)| r * v).sum::<f64>() / t;
    let (c, _) = pipeline.eroded_compliance_and_gradient(&mut el, &rho, &force);
    assert!(
        (vn - 0.4).abs() < 0.08,
        "nominal volume did not settle: {vn}"
    );
    assert!(
        c < 1e4,
        "eroded compliance stuck at the void scale: {c:.3e}"
    );
}
