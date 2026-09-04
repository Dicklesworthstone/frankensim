//! THE P2 MARQUEE STUDY RUNNER (plan §16.1, bead mye.1; [F] —
//! behind the `marquee` feature until its golden ledger joins nightly
//! CI): design optimization on a RAW SDF with NO MESH IN THE LOOP.
//! Geometry (an exact parametric SDF) → physics (CutFEM Poisson,
//! Nitsche Dirichlet holes) → gradient (the self-adjoint compliance
//! shape derivative `dJ = −∮_Γ (∂u/∂n)² V·n`, evaluated on the CutFEM
//! field — compliance is its own adjoint) → optimizer (projected
//! gradient with an area-feasibility projection) → certificate (the
//! COMPOSED per-iteration error ledger: exact geometry ⊗ DWR
//! discretization estimate ⊗ recomputed Euclidean algebraic residual,
//! colored by the
//! weakest input) → replayable ledger events. The LEVEL-SET VARIANT
//! the bead names: the design IS the zero set.
//!
//! The study: a heated plate (f = 1) with k cooling holes held at
//! temperature 0; minimize thermal compliance `J = ∫ f·u` over hole
//! radii at a fixed material-area budget — the canonical heat-sink
//! layout problem, meshed never.
use fs_cutfem::sdf::CutSdf;
use fs_cutfem::{FemParams, Quadtree, ScalarSample, Space};
use fs_dwr::{GoalContext, estimate, goal_value};
use fs_evidence::{Color, ColorRank, IntervalOp, compose};
use fs_ledger::hash_bytes;

/// The fixed projection iteration count and its admissible scale bracket.
/// With a ratio of at most 2^20, 80 bisections leave less than 2^-60 of
/// unresolved scale width before the radius-to-area conversion.
const PROJECTION_BISECTION_STEPS: usize = 80;
const MAX_RADIUS_RATIO: f64 = 1_048_576.0;
/// Largest permitted material-area residual after radius projection.
pub const AREA_PROJECTION_TOLERANCE: f64 = 1e-12;

/// Declared affine heat source on the normalized unit plate:
/// `f(x,y) = constant + x_slope*x + y_slope*y`.
/// The same source weights compliance and enters the state/adjoint solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalSource {
    /// Spatially uniform part.
    pub constant: f64,
    /// Source change across the unit x extent.
    pub x_slope: f64,
    /// Source change across the unit y extent.
    pub y_slope: f64,
}

impl ThermalSource {
    /// The original normalized unit-source problem.
    pub const UNIT: Self = Self {
        constant: 1.0,
        x_slope: 0.0,
        y_slope: 0.0,
    };

    fn value(self, x: f64, y: f64) -> f64 {
        self.constant + self.x_slope * x + self.y_slope * y
    }

    fn validate(self) -> Result<(), fs_cutfem::CutFemError> {
        let corners = [
            self.value(0.0, 0.0),
            self.value(1.0, 0.0),
            self.value(0.0, 1.0),
            self.value(1.0, 1.0),
        ];
        if [self.constant, self.x_slope, self.y_slope]
            .iter()
            .any(|v| !v.is_finite())
            || corners.iter().any(|v| !v.is_finite() || *v < 0.0)
            || corners.iter().all(|v| *v == 0.0)
        {
            return Err(fs_cutfem::CutFemError::InvalidFemInput { what: "affine thermal source must be finite, nonnegative on the unit plate and not identically zero".into() });
        }
        Ok(())
    }
}

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

    /// Compute the unprojected compliance of this design on a uniform quadtree.
    ///
    /// Evaluates the true PDE compliance functional J(u_h) = ∫ f·u over the
    /// design domain without applying area or radius projections.
    pub fn compliance(&self, level: u32) -> Result<f64, fs_cutfem::CutFemError> {
        let (j, _, _, _, _) = solve_and_grade(self, level)?;
        Ok(j)
    }
}

fn validate_study_input(
    design: &PlateWithHoles,
    config: &StudyConfig,
) -> Result<(), fs_cutfem::CutFemError> {
    macro_rules! require {
        ($condition:expr, $message:literal) => {
            if !$condition {
                return Err(fs_cutfem::CutFemError::InvalidFemInput {
                    what: $message.to_string(),
                });
            }
        };
    }
    require!(
        !design.radii.is_empty(),
        "marquee study needs at least one hole"
    );
    require!(
        design.centers.len() == design.radii.len(),
        "hole centers and radii must have matching lengths"
    );
    require!(
        design
            .centers
            .iter()
            .flatten()
            .all(|v| v.is_finite() && *v >= 0.0 && *v <= 1.0),
        "hole centers must be finite coordinates in the unit plate"
    );
    require!(
        design.radii.iter().all(|r| r.is_finite() && *r > 0.0),
        "hole radii must be positive and finite"
    );
    require!(
        config.step_size.is_finite() && config.step_size >= 0.0,
        "step size must be finite and nonnegative"
    );
    require!(
        config.area_target.is_finite() && config.area_target > 0.0 && config.area_target < 1.0,
        "area target must be finite and inside (0, 1)"
    );
    require!(
        config.r_min.is_finite()
            && config.r_max.is_finite()
            && config.r_min > 0.0
            && config.r_min <= config.r_max,
        "radius bounds must be finite and satisfy 0 < r_min <= r_max"
    );
    let ratio = config.r_max / config.r_min;
    require!(
        ratio.is_finite() && ratio <= MAX_RADIUS_RATIO,
        "radius ratio exceeds the deterministic projection envelope"
    );
    require!(
        hole_geometry_is_valid(design),
        "holes must be strictly interior and pairwise disjoint"
    );
    #[allow(clippy::cast_precision_loss)]
    let count = design.radii.len() as f64;
    let target_hole = (1.0 - config.area_target) / std::f64::consts::PI;
    require!(
        target_hole >= count * config.r_min * config.r_min
            && target_hole <= count * config.r_max * config.r_max,
        "area target is infeasible for the number of holes and radius bounds"
    );
    require!(
        config.level > 0 && config.level <= 12,
        "mesh level must lie in 1..=12"
    );
    Ok(())
}

/// Require the disk layout assumed by the summed-disc area and CutSdf model.
fn hole_geometry_is_valid(design: &PlateWithHoles) -> bool {
    if design.centers.len() != design.radii.len() {
        return false;
    }
    if !design
        .centers
        .iter()
        .zip(&design.radii)
        .all(|(center, radius)| {
            center[0] - radius > 0.0
                && center[0] + radius < 1.0
                && center[1] - radius > 0.0
                && center[1] + radius < 1.0
        })
    {
        return false;
    }
    for left in 0..design.radii.len() {
        for right in (left + 1)..design.radii.len() {
            let dx = design.centers[left][0] - design.centers[right][0];
            let dy = design.centers[left][1] - design.centers[right][1];
            let minimum_separation = design.radii[left] + design.radii[right];
            if dx.mul_add(dx, dy * dy) <= minimum_separation * minimum_separation {
                return false;
            }
        }
    }
    true
}

/// Refuse a projected layout whose fixed-ratio projection is not a valid disk
/// configuration. This study deliberately does not claim to search the
/// center-aware feasible set; callers must supply a target that remains valid
/// under its documented projection.
fn ensure_projected_hole_geometry(design: &PlateWithHoles) -> Result<(), fs_cutfem::CutFemError> {
    if hole_geometry_is_valid(design) {
        Ok(())
    } else {
        Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: "area projection produced a non-interior or overlapping hole layout; \
                   center-aware feasible allocation is unsupported"
                .to_string(),
        })
    }
}

/// Project candidate radii onto their box and the declared hole-area equality.
///
/// A single clamp-then-rescale can leave the sum outside the target when a
/// bound activates. Scaling followed by clamping is monotone, so fixed-count
/// bisection finds the deterministic feasible multiplier without changing the
/// background grid or remeshing the design.
fn project_radii_to_area(radii: &mut [f64], config: &StudyConfig) {
    let target_hole = (1.0 - config.area_target) / std::f64::consts::PI;
    let max_scale = config.r_max / config.r_min;
    for radius in radii.iter_mut() {
        *radius = radius.clamp(config.r_min, config.r_max);
    }
    let mut lo = 0.0;
    let mut hi = max_scale;

    for _ in 0..PROJECTION_BISECTION_STEPS {
        let scale = 0.5 * (lo + hi);
        let hole_area = radii
            .iter()
            .map(|radius| (scale * radius).clamp(config.r_min, config.r_max).powi(2))
            .sum::<f64>();
        if hole_area < target_hole {
            lo = scale;
        } else {
            hi = scale;
        }
    }

    for radius in radii.iter_mut() {
        *radius = (hi * *radius).clamp(config.r_min, config.r_max);
    }
    let projected_material_area =
        1.0 - std::f64::consts::PI * radii.iter().map(|radius| radius * radius).sum::<f64>();
    assert!(
        (projected_material_area - config.area_target).abs() <= AREA_PROJECTION_TOLERANCE,
        "radius projection residual exceeds {AREA_PROJECTION_TOLERANCE}"
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
    /// The compliance J(u_h) at `radii`, before the accepted step.
    pub compliance: f64,
    /// Material area.
    pub area: f64,
    /// Material volume in the unit-thickness 2-D study (equal to `area`).
    pub volume: f64,
    /// The design radii evaluated for `compliance`, before the accepted step.
    pub radii: Vec<f64>,
    /// Candidate radii retained after the Armijo decision.
    pub accepted_radii: Vec<f64>,
    /// Objective evaluated at `accepted_radii` (or `compliance` on exhaustion).
    pub accepted_compliance: f64,
    /// Geometry contribution of the solve at `accepted_radii`.
    pub accepted_cert_geometry: f64,
    /// DWR contribution of the solve at `accepted_radii`.
    pub accepted_cert_dwr: f64,
    /// Algebraic-residual contribution of the solve at `accepted_radii`.
    pub accepted_cert_algebraic: f64,
    /// Composed certificate color of the solve at `accepted_radii`.
    pub accepted_color: Color,
    /// Solver iterations of the solve at `accepted_radii`.
    pub accepted_solver_iters: usize,
    /// Cut-cell rules of the solve at `accepted_radii`.
    pub accepted_cut_cell_count: usize,
    /// Shape gradient dJ/dr per hole (the self-adjoint boundary form).
    pub gradient: Vec<f64>,
    /// Euclidean norm of `gradient`.
    pub gradient_norm: f64,
    /// The COMPOSED certificate components.
    pub cert_geometry: f64,
    /// |DWR estimate| (discretization).
    pub cert_dwr: f64,
    /// The per-iteration DWR estimate retained for JSONL rows.
    pub dwr_estimate: f64,
    /// Estimated goal-error contribution derived from a recomputed Euclidean
    /// algebraic residual (not CG's recursive residual estimate).
    pub cert_algebraic: f64,
    /// The composed color (weakest input).
    pub color: Color,
    /// Solver iterations (flat-cadence evidence: no remeshing spikes).
    pub solver_iters: usize,
    /// Certified cut-cell quadrature rules retained by the state solve.
    pub cut_cell_count: usize,
    /// Accepted Armijo step after the current objective evaluation.
    pub accepted_step: f64,
    /// Number of deterministic half-step retries before acceptance.
    pub backtracks: usize,
}

impl IterRecord {
    /// Serialize one deterministic JSONL object, including its terminating
    /// newline, for the per-iteration study ledger stream.
    #[must_use]
    pub fn jsonl_row(&self) -> String {
        format!(
            concat!(
                "{{\"iter\":{},\"compliance\":{:.17e},\"volume\":{:.17e},",
                "\"gradient_norm\":{:.17e},\"dwr_estimate\":{:.17e},",
                "\"cut_cell_count\":{},\"accepted_step\":{:.17e},",
                "\"backtracks\":{},\"accepted_compliance\":{:.17e},",
                "\"accepted_cert_geometry\":{:.17e},\"accepted_dwr_estimate\":{:.17e},",
                "\"accepted_cert_algebraic\":{:.17e},\"accepted_solver_iters\":{},",
                "\"accepted_cut_cell_count\":{},\"color_rank\":{},\"color_payload\":{},",
                "\"accepted_color_rank\":{},\"accepted_color_payload\":{}}}\n"
            ),
            self.iter,
            self.compliance,
            self.volume,
            self.gradient_norm,
            self.dwr_estimate,
            self.cut_cell_count,
            self.accepted_step,
            self.backtracks,
            self.accepted_compliance,
            self.accepted_cert_geometry,
            self.accepted_cert_dwr,
            self.accepted_cert_algebraic,
            self.accepted_solver_iters,
            self.accepted_cut_cell_count,
            color_rank_tag(self.color.rank()),
            self.color.payload_json(),
            color_rank_tag(self.accepted_color.rank()),
            self.accepted_color.payload_json(),
        )
    }
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
pub fn solve_and_grade(
    design: &PlateWithHoles,
    level: u32,
) -> Result<(f64, Vec<f64>, [f64; 3], usize, usize), fs_cutfem::CutFemError> {
    solve_and_grade_with_source(design, level, ThermalSource::UNIT)
}

/// Solve and grade the declared affine-source scalar problem.
pub fn solve_and_grade_with_source(
    design: &PlateWithHoles,
    level: u32,
    source: ThermalSource,
) -> Result<(f64, Vec<f64>, [f64; 3], usize, usize), fs_cutfem::CutFemError> {
    source.validate()?;
    let grid = Quadtree::uniform(level);
    let params = fem_params(level);
    let f = |x: f64, y: f64| source.value(x, y);
    let g = |_x: f64, _y: f64| 0.0;
    let space = Space::build(&grid, design, params)?;
    let sol = space.solve(&f, &g)?;
    let nodal = space.nodal_values(&sol.free, &g);
    // Compliance J = ∫ f·u over Ω (the DWR goal functional with w = f).
    let goal = GoalContext { weight: &f };
    let j = goal_value(&space, &nodal, &goal)?;
    // DWR discretization estimate for THIS goal (estimated color: DWR
    // constants are not guaranteed — the lmp4.4 rule).
    let dwr = estimate(&grid, design, params, &f, &g, &goal)?;
    // Self-adjoint shape gradient: dJ/dr_k = −∮_{Γ_k} (∂u/∂n)² dΓ.
    // For compliance J(u) = ∫ f·u over Ω with homogeneous Dirichlet data
    // on the hole boundaries Γ_k, enlarging a hole removes material from Ω
    // where u > 0, so the boundary moves opposite to the outward domain
    // normal n_Ω. Thus the radius derivative is negative (mq-004).
    // Midpoint quadrature.
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
            // det-ok: base 2, exact (4xnt)
            let h = 2.0f64.powi(-(i32::try_from(level).unwrap_or(6)) - 2);
            let q = [px + h * th.cos(), py + h * th.sin()];
            // Canonical fail-closed sampling (ay40): missing or
            // non-finite active evidence refuses instead of reading as
            // a plausible zero. A probe that leaves the solid past the
            // cooled hole boundary lands in a certified-Outside leaf,
            // where the homogeneous Dirichlet exterior value u = 0 is
            // the explicit physical meaning, not a fallback. The rim
            // clamp keeps legal probes inside the half-open background
            // box; validated hole layouts keep probes interior anyway.
            let q = [q[0].clamp(1e-9, 1.0 - 1e-9), q[1].clamp(1e-9, 1.0 - 1e-9)];
            let u_q = match space.sample_scalar(&nodal, q)? {
                ScalarSample::Active(v) => v,
                ScalarSample::CertifiedOutside => 0.0,
            };
            let dudn = u_q / h; // u = 0 on the hole boundary
            acc += dudn * dudn;
        }
        #[allow(clippy::cast_precision_loss)]
        let circ = std::f64::consts::TAU * r / samples as f64;
        grads.push(-(acc * circ));
    }
    let euclidean_rel_residual =
        sol.euclidean_rel_residual()
            .ok_or_else(|| fs_cutfem::CutFemError::InvalidFemInput {
                what: "marquee certificate requires a recomputed Euclidean solver residual"
                    .to_string(),
            })?;
    let cert = [0.0, dwr.eta_abs, euclidean_rel_residual * j.abs().max(1.0)];
    if !j.is_finite() || grads.iter().chain(cert.iter()).any(|v| !v.is_finite()) {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: "thermal solve produced non-finite objective, gradient or error estimate".into(),
        });
    }
    Ok((j, grads, cert, sol.iters, space.cut_rules().len()))
}

/// Sufficient decrease parameter (c1) in the Armijo line search condition.
pub const ARMIJO_SUFFICIENT_DECREASE: f64 = 1e-4;
/// Maximum number of deterministic half-step retries in Armijo line search.
pub const MAX_ARMIJO_BACKTRACKS: usize = 8;

/// Return an accepted box-and-area-feasible design, or retain the current one
/// after the bounded Armijo budget is exhausted.
pub fn armijo_next_design(
    design: &PlateWithHoles,
    config: &StudyConfig,
    current_objective: f64,
    gradient: &[f64],
    current_cert: [f64; 3],
    current_solver_iters: usize,
    current_cut_cell_count: usize,
) -> Result<(PlateWithHoles, f64, f64, usize, [f64; 3], usize, usize), fs_cutfem::CutFemError> {
    armijo_with_source(
        design,
        config,
        current_objective,
        gradient,
        current_cert,
        current_solver_iters,
        current_cut_cell_count,
        ThermalSource::UNIT,
    )
}

fn armijo_with_source(
    design: &PlateWithHoles,
    config: &StudyConfig,
    current_objective: f64,
    gradient: &[f64],
    current_cert: [f64; 3],
    current_solver_iters: usize,
    current_cut_cell_count: usize,
    source: ThermalSource,
) -> Result<(PlateWithHoles, f64, f64, usize, [f64; 3], usize, usize), fs_cutfem::CutFemError> {
    if config.step_size == 0.0 {
        return Ok((
            design.clone(),
            current_objective,
            0.0,
            0,
            current_cert,
            current_solver_iters,
            current_cut_cell_count,
        ));
    }

    for backtracks in 0..=MAX_ARMIJO_BACKTRACKS {
        #[allow(clippy::cast_precision_loss)]
        let step = config.step_size * 0.5_f64.powi(backtracks as i32);
        let mut candidate = design.clone();
        for (radius, derivative) in candidate.radii.iter_mut().zip(gradient) {
            *radius = (*radius - step * derivative).clamp(config.r_min, config.r_max);
        }
        project_radii_to_area(&mut candidate.radii, config);
        if !hole_geometry_is_valid(&candidate) {
            continue;
        }

        let directional_derivative = gradient
            .iter()
            .zip(&candidate.radii)
            .zip(&design.radii)
            .map(|((derivative, candidate_radius), current_radius)| {
                derivative * (candidate_radius - current_radius)
            })
            .sum::<f64>();
        if directional_derivative >= 0.0 {
            continue;
        }

        let (
            candidate_objective,
            _,
            candidate_cert,
            candidate_solver_iters,
            candidate_cut_cell_count,
        ) = solve_and_grade_with_source(&candidate, config.level, source)?;
        if candidate_objective
            <= current_objective + ARMIJO_SUFFICIENT_DECREASE * directional_derivative
        {
            return Ok((
                candidate,
                candidate_objective,
                step,
                backtracks,
                candidate_cert,
                candidate_solver_iters,
                candidate_cut_cell_count,
            ));
        }
    }

    Ok((
        design.clone(),
        current_objective,
        0.0,
        MAX_ARMIJO_BACKTRACKS,
        current_cert,
        current_solver_iters,
        current_cut_cell_count,
    ))
}

fn certificate_color(cert: [f64; 3]) -> Color {
    // COMPOSED certificate color: exact geometry (verified) ⊗ DWR
    // (estimated) ⊗ algebraic residual (estimated) — weakest wins.
    compose(
        &compose(
            // declared-color-ok: exact-arithmetic identity leaf inside a compose() whose weakest-input rule keeps the result estimated (6pf9)
            &Color::Verified { lo: 0.0, hi: 0.0 },
            &Color::Estimated {
                estimator: "dwr(compliance)".to_string(),
                dispersion: cert[1],
            },
            IntervalOp::Add,
        ),
        &Color::Estimated {
            estimator: "recomputed-euclidean-cg-residual".to_string(),
            dispersion: cert[2],
        },
        IntervalOp::Add,
    )
}

/// Run the marquee study: projected gradient on hole radii at a fixed
/// area budget, with the composed certificate recorded per iteration.
/// Deterministic; the trace hash is the replay-equality witness.
///
/// # Errors
/// CutFEM build/solve teaching errors.
pub fn run_study(
    design: PlateWithHoles,
    config: &StudyConfig,
) -> Result<StudyReport, fs_cutfem::CutFemError> {
    let mut runner = StudyRunner::new(design, config.clone())?;
    while runner.advance()? {}
    Ok(runner.report())
}

/// Incremental thermal radius optimizer. An unsuccessful step leaves the
/// accepted design and trace unchanged. Callers can meter, persist, or stop
/// between steps; a step includes at most nine bounded line-search solves.
#[derive(Debug, Clone)]
pub struct StudyRunner {
    design: PlateWithHoles,
    config: StudyConfig,
    iterations: Vec<IterRecord>,
    source: ThermalSource,
}

impl StudyRunner {
    /// Admit the geometry and project it onto the requested area before any solve.
    pub fn new(
        design: PlateWithHoles,
        config: StudyConfig,
    ) -> Result<Self, fs_cutfem::CutFemError> {
        Self::new_with_source(design, config, ThermalSource::UNIT)
    }

    /// Admit a study with an explicit spatial heat load. Every line-search
    /// trial uses this same load; it is also bound into the replay identity.
    pub fn new_with_source(
        mut design: PlateWithHoles,
        config: StudyConfig,
        source: ThermalSource,
    ) -> Result<Self, fs_cutfem::CutFemError> {
        source.validate()?;
        validate_study_input(&design, &config)?;
        project_radii_to_area(&mut design.radii, &config);
        ensure_projected_hole_geometry(&design)?;
        Ok(Self {
            design,
            config,
            iterations: Vec::new(),
            source,
        })
    }

    /// The currently accepted geometry, including initial area projection.
    pub fn design(&self) -> &PlateWithHoles {
        &self.design
    }

    /// Accepted transition records, in execution order.
    pub fn iterations(&self) -> &[IterRecord] {
        &self.iterations
    }

    /// Advance once, or return false once the declared iteration count is reached.
    pub fn advance(&mut self) -> Result<bool, fs_cutfem::CutFemError> {
        let iter = self.iterations.len();
        if iter == self.config.steps {
            return Ok(false);
        }
        let design = &self.design;
        let config = &self.config;
        let (j, grads, cert, iters, cut_cell_count) =
            solve_and_grade_with_source(design, config.level, self.source)?;
        let (
            next_design,
            accepted_compliance,
            accepted_step,
            backtracks,
            accepted_cert,
            accepted_solver_iters,
            accepted_cut_cell_count,
        ) = armijo_with_source(
            design,
            config,
            j,
            &grads,
            cert,
            iters,
            cut_cell_count,
            self.source,
        )?;
        let color = certificate_color(cert);
        let accepted_color = certificate_color(accepted_cert);
        let gradient_norm = grads
            .iter()
            .map(|gradient| gradient * gradient)
            .sum::<f64>()
            .sqrt();
        self.iterations.push(IterRecord {
            iter,
            compliance: j,
            area: design.area(),
            volume: design.area(),
            radii: design.radii.clone(),
            accepted_radii: next_design.radii.clone(),
            accepted_compliance,
            accepted_cert_geometry: accepted_cert[0],
            accepted_cert_dwr: accepted_cert[1],
            accepted_cert_algebraic: accepted_cert[2],
            accepted_color,
            accepted_solver_iters,
            accepted_cut_cell_count,
            gradient: grads.clone(),
            gradient_norm,
            cert_geometry: cert[0],
            cert_dwr: cert[1],
            dwr_estimate: cert[1],
            cert_algebraic: cert[2],
            color,
            solver_iters: iters,
            cut_cell_count,
            accepted_step,
            backtracks,
        });
        self.design = next_design;
        Ok(true)
    }

    /// Snapshot the exact accepted trace. No solve or area re-projection occurs.
    pub fn report(&self) -> StudyReport {
        let mut trace_hash = canonical_trace_hash(&self.iterations, &self.design, &self.config);
        if self.source != ThermalSource::UNIT {
            let mut bytes = b"fs-marquee-affine-source-trace-v1\0".to_vec();
            bytes.extend_from_slice(trace_hash.as_bytes());
            append_f64(&mut bytes, self.source.constant);
            append_f64(&mut bytes, self.source.x_slope);
            append_f64(&mut bytes, self.source.y_slope);
            trace_hash = hash_bytes(&bytes).to_hex();
        }
        StudyReport {
            trace_hash,
            iterations: self.iterations.clone(),
            design: self.design.clone(),
        }
    }
}

fn append_usize(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(
        &u64::try_from(value)
            .expect("a Rust allocation length fits u64")
            .to_le_bytes(),
    );
}

fn append_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn append_f64_slice(bytes: &mut Vec<u8>, values: &[f64]) {
    append_usize(bytes, values.len());
    for value in values {
        append_f64(bytes, *value);
    }
}

fn color_rank_tag(rank: ColorRank) -> u8 {
    match rank {
        ColorRank::Estimated => 0,
        ColorRank::Validated => 1,
        ColorRank::Verified => 2,
    }
}

/// Canonical replay identity for every recorded solve/acceptance transition
/// and the final retained geometry. Floats are committed by IEEE bits, and
/// colors use their versioned canonical payload bytes.
fn canonical_trace_hash(
    iterations: &[IterRecord],
    design: &PlateWithHoles,
    config: &StudyConfig,
) -> String {
    let mut bytes = b"fs-marquee-study-trace-v3\0".to_vec();
    bytes.extend_from_slice(&config.level.to_le_bytes());
    append_usize(&mut bytes, config.steps);
    append_f64(&mut bytes, config.step_size);
    append_f64(&mut bytes, config.area_target);
    append_f64(&mut bytes, config.r_min);
    append_f64(&mut bytes, config.r_max);
    append_usize(&mut bytes, iterations.len());
    for record in iterations {
        append_usize(&mut bytes, record.iter);
        append_f64(&mut bytes, record.compliance);
        append_f64(&mut bytes, record.area);
        append_f64(&mut bytes, record.volume);
        append_f64_slice(&mut bytes, &record.radii);
        append_f64_slice(&mut bytes, &record.accepted_radii);
        append_f64(&mut bytes, record.accepted_compliance);
        append_f64(&mut bytes, record.accepted_cert_geometry);
        append_f64(&mut bytes, record.accepted_cert_dwr);
        append_f64(&mut bytes, record.accepted_cert_algebraic);
        bytes.push(color_rank_tag(record.accepted_color.rank()));
        let accepted_color = record.accepted_color.canonical_bytes();
        append_usize(&mut bytes, accepted_color.len());
        bytes.extend_from_slice(&accepted_color);
        append_usize(&mut bytes, record.accepted_solver_iters);
        append_usize(&mut bytes, record.accepted_cut_cell_count);
        append_f64_slice(&mut bytes, &record.gradient);
        append_f64(&mut bytes, record.gradient_norm);
        append_f64(&mut bytes, record.cert_geometry);
        append_f64(&mut bytes, record.cert_dwr);
        append_f64(&mut bytes, record.dwr_estimate);
        append_f64(&mut bytes, record.cert_algebraic);
        bytes.push(color_rank_tag(record.color.rank()));
        let color = record.color.canonical_bytes();
        append_usize(&mut bytes, color.len());
        bytes.extend_from_slice(&color);
        append_usize(&mut bytes, record.solver_iters);
        append_usize(&mut bytes, record.cut_cell_count);
        append_f64(&mut bytes, record.accepted_step);
        append_usize(&mut bytes, record.backtracks);
    }
    append_usize(&mut bytes, design.centers.len());
    for center in &design.centers {
        append_f64(&mut bytes, center[0]);
        append_f64(&mut bytes, center[1]);
    }
    append_f64_slice(&mut bytes, &design.radii);
    hash_bytes(&bytes).to_hex()
}
