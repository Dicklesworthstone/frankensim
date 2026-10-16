//! Mixed direction/orientation/frame optimization with exact reverse gradients.
//! Dimensionless analytic example, not an identified physical model.
use fs_ascent::{ReverseManifoldStudy, StopRule};
use fs_opt::{Manifold, Problem, ProblemBuilder, Sense, ReverseProblem, reverse::ReverseLimits};
use fs_qty::Dims;

/// An R2 position, unit direction, SO(3) rotation and two-column frame.
pub fn fixture() -> (Problem, Vec<f64>) {
    let mut b = ProblemBuilder::new();
    let specs = [
        ("position", Manifold::Rn { dim: 2 }, vec![0.0, -0.25], 1.0),
        ("direction", Manifold::Sphere { ambient: 3 }, vec![1.0, 0.0, 0.0], 0.5),
        // Only vector quaternion components enter the loss: q and -q agree.
        ("rotation", Manifold::So3, vec![0.0; 4], 1.0),
        ("frame", Manifold::Stiefel { n: 3, p: 2 }, vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0], 0.5),
    ];
    let mut roots = Vec::new();
    let mut total = b.konst(0.0, Dims::NONE).unwrap();
    for (name, man, target, weight) in specs {
        let variable = b.var(name, man, Dims::NONE).unwrap();
        let root = b.var_ref(variable).unwrap();
        roots.push(root);
        for (i, value) in target.into_iter().enumerate() {
            if matches!(man, Manifold::So3) && i == 0 { continue; }
            let coordinate = b.component(root, i as u32).unwrap();
            let target = b.konst(value, Dims::NONE).unwrap();
            let error = b.sub(coordinate, target).unwrap();
            let squared = b.powi(error, 2).unwrap();
            let weight = b.konst(weight, Dims::NONE).unwrap();
            let term = b.mul(weight, squared).unwrap();
            total = b.add(total, term).unwrap();
        }
    }
    // Couple two different factor types, rather than solving four independent fits.
    let p = b.component(roots[0], 0).unwrap();
    let u = b.component(roots[1], 1).unwrap();
    let difference = b.sub(p, u).unwrap();
    let squared = b.powi(difference, 2).unwrap();
    let weight = b.konst(0.2, Dims::NONE).unwrap();
    let coupled = b.mul(weight, squared).unwrap();
    total = b.add(total, coupled).unwrap();
    b.objective(total, Sense::Minimize, 1.0).unwrap();
    (b.finish(), vec![-1.2, 0.4, 0.6, 0.8, 0.0, 0.8, 0.6, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.4_f64.cos(), 0.4_f64.sin()])
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (problem, point) = fixture();
    let oracle = ReverseProblem::new(&problem, ReverseLimits { max_nodes: 4096, max_scalar_slots: 16384 })?;
    let mut study = ReverseManifoldStudy::new(&oracle, &point, 7, None)?;
    let report = study.run(&StopRule::Any(vec![StopRule::GradNorm(1e-7), StopRule::Budget(200)]), 100, None)?;
    println!("stop={:?} objective={:.12e} gradient={:.12e} iterations={} samples={} point_dim={} parameter_dim={}",
        report.reason, report.f, report.grad_norm, report.iters, report.evals,
        study.point().len(), study.gradient().len());
    println!("accepted={:?}", study.point());
    if report.reason != fs_ascent::StopReason::GradNorm { return Err("study did not meet the gradient threshold".into()); }
    Ok(())
}
