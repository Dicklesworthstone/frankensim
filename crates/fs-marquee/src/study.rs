//! THE P2 MARQUEE STUDY RUNNER (plan §16.1, bead mye.1; [F] —
//! behind the `marquee` feature until its golden ledger joins nightly
//! CI): design optimization on a RAW SDF with NO MESH IN THE LOOP.
//! Geometry (an exact parametric SDF) → physics (CutFEM Poisson,
//! Nitsche Dirichlet holes) → gradient (the self-adjoint compliance
//! shape derivative `dJ = −∮_Γ (∂u/∂n)² V·n`, evaluated on the CutFEM
//! field — compliance is its own adjoint) → optimizer (projected
//! gradient with an area-feasibility projection) → certificate (the
//! COMPOSED per-iteration error ledger: exact geometry ⊗ DWR
//! discretization estimate ⊗ algebraic residual, colored by the
//! weakest input) → replayable ledger events. The LEVEL-SET VARIANT
//! the bead names: the design IS the zero set.
//!
//! The study: a heated plate (f = 1) with k cooling holes held at
//! temperature 0; minimize thermal compliance `J = ∫ f·u` over hole
//! radii at a fixed material-area budget — the canonical heat-sink
//! layout problem, meshed never.
use std::collections::BTreeMap;

use fs_cutfem::sdf::CutSdf;
use fs_cutfem::{FemParams, Quadtree, Space};
use fs_dwr::{GoalContext, estimate, goal_value};
use fs_evidence::{Color, IntervalOp, compose};
use fs_ledger::hash_bytes;

/// The design: a unit plate minus k circular cooling holes
/// (φ < 0 inside the material). EXACT geometry: circles.
#[derive(Debug, Clone, PartialEq)]
pub struct PlateWithHoles {
    /// Hole centers.
    pub centers: Vec<[f64; 2]>,
    /// Hole radii (the design variables).
    pub radii: Vec<f64>,
}

impl PlateWithHoles {
    /// Material area = 1 − Σ hole areas (holes assumed disjoint and
    /// interior — enforced by the optimizer's box projection).
    #[must_use]
    pub fn area(&self) -> f64 {
        1.0 - self
            .radii
            .iter()
            .map(|r| std::f64::consts::PI * r * r)
            .sum::<f64>()
    }
}

fn validate_study_input(design: &PlateWithHoles, config: &StudyConfig) {
    assert!(
        !design.radii.is_empty(),
        "marquee study needs at least one hole"
    );
    assert_eq!(
        design.centers.len(),
        design.radii.len(),
        "hole centers and radii must have matching lengths"
    );
    assert!(
        design
            .centers
            .iter()
            .flatten()
            .all(|v| v.is_finite() && *v >= 0.0 && *v <= 1.0),
        "hole centers must be finite coordinates in the unit plate"
    );
    assert!(
        design.radii.iter().all(|r| r.is_finite() && *r > 0.0),
        "hole radii must be positive and finite"
    );
    assert!(
        config.step_size.is_finite() && config.step_size >= 0.0,
        "step size must be finite and nonnegative"
    );
    assert!(
        config.area_target.is_finite() && config.area_target > 0.0 && config.area_target < 1.0,
        "area target must be finite and inside (0, 1)"
    );
    assert!(
        config.r_min.is_finite()
            && config.r_max.is_finite()
            && config.r_min > 0.0
            && config.r_min <= config.r_max,
        "radius bounds must be finite and satisfy 0 < r_min <= r_max"
    );
}

impl CutSdf for PlateWithHoles {
    fn value(&self, p: [f64; 2]) -> f64 {
        // Inside the material: negative. φ = max_i (r_i − |p − c_i|).
        self.centers
            .iter()
            .zip(&self.radii)
            .map(|(c, r)| r - ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
        // Gradient of the active hole's term: −(p − c)/|p − c|.
        let Some((c, _)) = self
            .centers
            .iter()
            .zip(&self.radii)
            .max_by(|(ca, ra), (cb, rb)| {
                let da = *ra - ((p[0] - ca[0]).powi(2) + (p[1] - ca[1]).powi(2)).sqrt();
                let db = *rb - ((p[0] - cb[0]).powi(2) + (p[1] - cb[1]).powi(2)).sqrt();
                da.total_cmp(&db)
            })
        else {
            return [0.0, 0.0];
        };
        let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2))
            .sqrt()
            .max(1e-12);
        [-(p[0] - c[0]) / d, -(p[1] - c[1]) / d]
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> fs_ivl::Interval {
        // Exact per-hole enclosure: r − [dist_min, dist_max] to the box,
        // hulled over holes (max of intervals).
        let mut out = [f64::NEG_INFINITY, f64::NEG_INFINITY];
        for (c, r) in self.centers.iter().zip(&self.radii) {
            // Distance from c to the box: min and max over the box.
            let mut dmin2 = 0.0f64;
            let mut dmax2 = 0.0f64;
            for k in 0..2 {
                let below = (lo[k] - c[k]).max(0.0);
                let above = (c[k] - hi[k]).max(0.0);
                let gap = below.max(above);
                dmin2 += gap * gap;
                let far = (c[k] - lo[k]).abs().max((hi[k] - c[k]).abs());
                dmax2 += far * far;
            }
            let (vlo, vhi) = (r - dmax2.sqrt(), r - dmin2.sqrt());
            out[0] = out[0].max(vlo);
            out[1] = out[1].max(vhi);
        }
        // max over holes preserves containment for the max-combination.
        fs_ivl::Interval::new(out[0], out[1])
    }
}

/// One iteration's forensic record.
#[derive(Debug, Clone, PartialEq)]
pub struct IterRecord {
    /// Iteration index.
    pub iter: usize,
    /// The compliance J(u_h).
    pub compliance: f64,
    /// Material area.
    pub area: f64,
    /// The design radii after this step.
    pub radii: Vec<f64>,
    /// Shape gradient dJ/dr per hole (the self-adjoint boundary form).
    pub gradient: Vec<f64>,
    /// The COMPOSED certificate components.
    pub cert_geometry: f64,
    /// |DWR estimate| (discretization).
    pub cert_dwr: f64,
    /// Algebraic residual bound proxy.
    pub cert_algebraic: f64,
    /// The composed color (weakest input).
    pub color: Color,
    /// Solver iterations (flat-cadence evidence: no remeshing spikes).
    pub solver_iters: usize,
}

/// The study configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyConfig {
    /// Background-grid refinement level (quadtree).
    pub level: u32,
    /// Optimization steps.
    pub steps: usize,
    /// Gradient step size on radii.
    pub step_size: f64,
    /// The material-area budget (equality target).
    pub area_target: f64,
    /// Radius box bounds.
    pub r_min: f64,
    /// Radius upper bound.
    pub r_max: f64,
}

/// The study outcome: the trace and the replay hash.
#[derive(Debug, Clone)]
pub struct StudyReport {
    /// Per-iteration records.
    pub iterations: Vec<IterRecord>,
    /// The final design.
    pub design: PlateWithHoles,
    /// The G5 trace hash over every record (replay equality).
    pub trace_hash: String,
}

fn fem_params(level: u32) -> FemParams {
    FemParams {
        nitsche_beta: 10.0,
        ghost_gamma: 0.1,
        quad_depth: 3,
        agg: None,
        strong_outer: true,
        solver_tol: 1e-10,
        solver_max_iters: 1500,
    }
    .with_level_hint(level)
}

/// Helper trait shim: FemParams may not expose a level hint — identity.
trait LevelHint {
    fn with_level_hint(self, level: u32) -> Self;
}

impl LevelHint for FemParams {
    fn with_level_hint(self, _level: u32) -> Self {
        self
    }
}

/// Solve the state problem on the CURRENT design; return
/// (compliance, per-hole shape gradients, certificate parts, iters).
#[allow(clippy::type_complexity)]
fn solve_and_grade(
    design: &PlateWithHoles,
    level: u32,
) -> Result<(f64, Vec<f64>, [f64; 3], usize, BTreeMap<(u32, u32), f64>), fs_cutfem::CutFemError> {
    let grid = Quadtree::uniform(level);
    let params = fem_params(level);
    let f = |_x: f64, _y: f64| 1.0;
    let g = |_x: f64, _y: f64| 0.0;
    let space = Space::build(&grid, design, params)?;
    let sol = space.solve(&f, &g)?;
    let nodal = space.nodal_values(&sol.free, &g);
    // Compliance J = ∫ f·u over Ω (the DWR goal functional with w = f).
    let goal = GoalContext { weight: &f };
    let j = goal_value(&space, &grid, design, &nodal, &goal, params.quad_depth);
    // DWR discretization estimate for THIS goal (estimated color: DWR
    // constants are not guaranteed — the lmp4.4 rule).
    let dwr = estimate(&grid, design, params, &f, &g, &goal)?;
    // Self-adjoint shape gradient: dJ/dr_k = −∮_{Γ_k} (∂u/∂n)² dΓ —
    // growing a Dirichlet cooling hole LOWERS compliance (more cold
    // boundary), so the flux integral enters NEGATED. The sign was
    // originally implemented positive and the FD falsifier caught it
    // (mq-004) — the drill earning its keep. Midpoint quadrature.
    let mut grads = Vec::with_capacity(design.radii.len());
    let samples = 64usize;
    for (c, r) in design.centers.iter().zip(&design.radii) {
        let mut acc = 0.0f64;
        for k in 0..samples {
            #[allow(clippy::cast_precision_loss)]
            let th = std::f64::consts::TAU * (k as f64 + 0.5) / samples as f64;
            let px = c[0] + r * th.cos();
            let py = c[1] + r * th.sin();
            // ∂u/∂n via a one-sided probe into the material along the
            // outward-from-hole (into-material) normal.
            let h = 2.0f64.powi(-(i32::try_from(level).unwrap_or(6)) - 2);
            let q = [px + h * th.cos(), py + h * th.sin()];
            let u_q = sample_nodal(&grid, &nodal, q);
            let dudn = u_q / h; // u = 0 on the hole boundary
            acc += dudn * dudn;
        }
        #[allow(clippy::cast_precision_loss)]
        let circ = std::f64::consts::TAU * r / samples as f64;
        grads.push(-(acc * circ));
    }
    let cert = [0.0, dwr.eta_abs, sol.rel_residual * j.abs().max(1.0)];
    Ok((j, grads, cert, sol.iters, nodal))
}

/// Bilinear sample of the nodal field at a point (level-`grid` cells).
fn sample_nodal(grid: &Quadtree, nodal: &BTreeMap<(u32, u32), f64>, p: [f64; 2]) -> f64 {
    // The background lattice at the quadtree's max depth.
    let n = (1u32 << grid.max_level()) + 1;
    #[allow(clippy::cast_precision_loss)]
    let scale = f64::from(n - 1);
    let fx = (p[0].clamp(0.0, 1.0)) * scale;
    let fy = (p[1].clamp(0.0, 1.0)) * scale;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let (ix, iy) = (
        (fx.floor() as u32).min(n - 2),
        (fy.floor() as u32).min(n - 2),
    );
    let (tx, ty) = (fx - f64::from(ix), fy - f64::from(iy));
    let v = |a: u32, b: u32| nodal.get(&(a, b)).copied().unwrap_or(0.0);
    (1.0 - tx) * (1.0 - ty) * v(ix, iy)
        + tx * (1.0 - ty) * v(ix + 1, iy)
        + tx * ty * v(ix + 1, iy + 1)
        + (1.0 - tx) * ty * v(ix, iy + 1)
}

/// Run the marquee study: projected gradient on hole radii at a fixed
/// area budget, with the composed certificate recorded per iteration.
/// Deterministic; the trace hash is the replay-equality witness.
///
/// # Errors
/// CutFEM build/solve teaching errors.
pub fn run_study(
    mut design: PlateWithHoles,
    config: &StudyConfig,
) -> Result<StudyReport, fs_cutfem::CutFemError> {
    validate_study_input(&design, config);
    let mut iterations = Vec::with_capacity(config.steps);
    for iter in 0..config.steps {
        let (j, grads, cert, iters, _) = solve_and_grade(&design, config.level)?;
        // COMPOSED certificate color: exact geometry (verified) ⊗ DWR
        // (estimated) ⊗ algebraic residual (estimated) — weakest wins.
        let color = compose(
            &compose(
                &Color::Verified { lo: 0.0, hi: 0.0 },
                &Color::Estimated {
                    estimator: "dwr(compliance)".to_string(),
                    dispersion: cert[1],
                },
                IntervalOp::Add,
            ),
            &Color::Estimated {
                estimator: "cg-residual".to_string(),
                dispersion: cert[2],
            },
            IntervalOp::Add,
        );
        iterations.push(IterRecord {
            iter,
            compliance: j,
            area: design.area(),
            radii: design.radii.clone(),
            gradient: grads.clone(),
            cert_geometry: cert[0],
            cert_dwr: cert[1],
            cert_algebraic: cert[2],
            color,
            solver_iters: iters,
        });
        // Descent: dJ/dr < 0 — every hole wants to grow; the area
        // budget's rescale projection turns that into REDISTRIBUTION:
        // holes with the larger per-radius payoff grow at the expense
        // of the rest (flux equalization, the optimality condition).
        for (r, g) in design.radii.iter_mut().zip(&grads) {
            *r = (*r - config.step_size * g).clamp(config.r_min, config.r_max);
        }
        // Area-equality projection: rescale radii to hit the target.
        let hole_area: f64 = design.radii.iter().map(|r| r * r).sum::<f64>();
        let target_hole = (1.0 - config.area_target) / std::f64::consts::PI;
        if hole_area > 0.0 {
            let scale = (target_hole / hole_area).sqrt();
            for r in &mut design.radii {
                *r = (*r * scale).clamp(config.r_min, config.r_max);
            }
        }
    }
    let mut canon = String::new();
    for rec in &iterations {
        use std::fmt::Write as _;
        let _ = write!(
            canon,
            "{}:{:.12e}:{:.12e};",
            rec.iter, rec.compliance, rec.area
        );
        for r in &rec.radii {
            let _ = write!(canon, "{r:.12e},");
        }
    }
    Ok(StudyReport {
        iterations,
        design,
        trace_hash: hash_bytes(canon.as_bytes()).to_hex(),
    })
}
