//! Minimum-mass rectangular cantilever sizing through the live IR, not callbacks.
//! cargo run -p fs-ascent --example reverse_sqp_beam -- [allowable_stress_pa]
//!
//! Illustrative linear Euler-Bernoulli algebra, not a validated structural design:
//! I=b*h^3/12, tip deflection=F*L^3/(3*E*I), root stress=6*F*L/(b*h^2).
//! The aspect constraint is h=2*b. All solver coordinates/residuals/objectives
//! are explicitly nondimensionalized by the SI reference values below.
use fs_ascent::{ReverseSqpReport, ReverseSqpStudy};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, EvalLimit, Manifold, OptError, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

const LENGTH_M: f64 = 1.0;
const WIDTH_REFERENCE_M: f64 = 0.05;
const HEIGHT_REFERENCE_M: f64 = 0.10;
const LOAD_N: f64 = 1000.0;
const YOUNG_PA: f64 = 200e9;
const DENSITY_KG_M3: f64 = 7800.0;
const DEFLECTION_LIMIT_M: f64 = 0.001;
const MIN_RATIO: f64 = 0.2;

fn problem(stress_limit_pa: f64) -> Result<Problem, OptError> {
    let mut b = ProblemBuilder::new();
    let v = b.var("width/0.05m,height/0.10m", Manifold::Rn { dim: 2 }, Dims::NONE)?;
    let r = b.var_ref(v)?;
    let x = b.component(r, 0)?;
    let y = b.component(r, 1)?;
    let mass_ratio = b.mul(x, y)?;
    b.objective(mass_ratio, Sense::Minimize, 1.0)?;
    let one = b.konst(1.0, Dims::NONE)?;
    let y2 = b.powi(y, 2)?;
    let y3 = b.powi(y, 3)?;
    let stress_denominator = b.mul(x, y2)?;
    let deflection_denominator = b.mul(x, y3)?;
    let stress = 6.0 * LOAD_N * LENGTH_M
        / (WIDTH_REFERENCE_M * HEIGHT_REFERENCE_M.powi(2) * stress_limit_pa);
    let inertia = WIDTH_REFERENCE_M * HEIGHT_REFERENCE_M.powi(3) / 12.0;
    let deflection = LOAD_N * LENGTH_M.powi(3) / (3.0 * YOUNG_PA * inertia * DEFLECTION_LIMIT_M);
    let ks = b.konst(stress, Dims::NONE)?;
    let kd = b.konst(deflection, Dims::NONE)?;
    let stress_ratio = b.div(ks, stress_denominator)?;
    let deflection_ratio = b.div(kd, deflection_denominator)?;
    let stress_excess = b.sub(stress_ratio, one)?;
    let deflection_excess = b.sub(deflection_ratio, one)?;
    let aspect = b.sub(y, x)?;
    b.constraint(stress_excess, ConstraintKind::LeZero, "stress/allowable-1")?;
    b.constraint(aspect, ConstraintKind::EqZero, "height_ratio-width_ratio")?;
    b.constraint(deflection_excess, ConstraintKind::LeZero, "deflection/limit-1")?;
    let minimum = b.konst(MIN_RATIO, Dims::NONE)?;
    let min_x = b.sub(minimum, x)?;
    let min_y = b.sub(minimum, y)?;
    b.constraint(min_x, ConstraintKind::LeZero, "minimum_width_ratio-x")?;
    b.constraint(min_y, ConstraintKind::LeZero, "minimum_height_ratio-y")?;
    b.set_eval_limit(EvalLimit::Limited(NonZeroU64::new(200).unwrap()));
    Ok(b.finish())
}

fn solve(stress_limit_pa: f64) -> Result<ReverseSqpReport, Box<dyn std::error::Error>> {
    if !stress_limit_pa.is_finite() || stress_limit_pa <= 0.0 {
        return Err("allowable stress must be finite and positive".into());
    }
    let problem = problem(stress_limit_pa)?;
    let oracle = ReverseProblem::new(&problem, ReverseLimits { max_nodes: 100, max_scalar_slots: 1024 })?;
    let mut study = ReverseSqpStudy::new(&oracle, &[1.0, 1.0], 16, None)?;
    Ok(study.run(1e-7, 80, None)?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 1 { return Err("usage: reverse_sqp_beam [allowable_stress_pa]".into()); }
    let stress_limit_pa = args.first().map_or(Ok(150e6), |arg| arg.parse::<f64>())?;
    let report = solve(stress_limit_pa)?;
    let mass_reference = DENSITY_KG_M3 * WIDTH_REFERENCE_M * HEIGHT_REFERENCE_M * LENGTH_M;
    println!("stop={:?}; converged={}; evaluations={}; accepted_steps={}",
        report.stop, report.solution.converged, report.solution.evals, report.solution.iters);
    println!("width_m={:.12}; height_m={:.12}; mass_kg={:.12}",
        WIDTH_REFERENCE_M * report.solution.x[0], HEIGHT_REFERENCE_M * report.solution.x[1],
        mass_reference * report.solution.f);
    println!("normalized_kkt={:?}; declared_constraint_multipliers={:?}",
        report.solution.kkt, report.constraint_multipliers);
    if !report.solution.converged { return Err("study stopped without satisfying its KKT tolerance".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deflection_limited_beam_matches_the_closed_form() {
        let report = solve(150e6).unwrap();
        let ratio = 0.4_f64.sqrt().sqrt();
        assert!(report.solution.converged, "{report:?}");
        assert!((report.solution.x[0] - ratio).abs() < 1e-6);
        assert!((report.solution.x[1] - ratio).abs() < 1e-6);
        assert!(report.constraint_multipliers[2] > 0.0);
        assert!(report.constraint_multipliers[0].abs() < 1e-6);
    }

    #[test]
    fn changing_allowable_stress_changes_the_active_constraint_and_design() {
        let report = solve(10e6).unwrap();
        let ratio = 1.2_f64.cbrt();
        assert!(report.solution.converged, "{report:?}");
        assert!((report.solution.x[0] - ratio).abs() < 1e-6);
        assert!((report.solution.x[1] - ratio).abs() < 1e-6);
        assert!(report.constraint_multipliers[0] > 0.0);
        assert!(report.constraint_multipliers[2].abs() < 1e-6);
    }
}
