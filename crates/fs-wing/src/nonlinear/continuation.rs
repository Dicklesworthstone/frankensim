//! Caller-owned continuation of the unchanged section/induction equations.
//!
//! A previous circulation is a numerical branch anchor, not additional air
//! inertia or physical lag. Incidence continuity does not require identical
//! regime labels: an attached/blended transition can be smooth. The declared
//! incidence, row-cancellation and filament-clearance domain is admission,
//! never clipping, regularization, a uniqueness certificate or calibration.
//!
//! The Newton/FGMRES and Armijo algorithms remain owned by `fs-solver`;
//! see Kelley, *Iterative Methods for Linear and Nonlinear Equations*,
//! chapters 6 and 8: <https://doi.org/10.1137/1.9781611970944.ch6> and
//! <https://doi.org/10.1137/1.9781611970944.ch8>. Here only the original
//! section residual and its geometry/normalization chain rule are assembled.

use super::{
    GAMMA_TOL, InfluenceOperator, NonlinearReport, Panel, ROW_SUM_FLOOR, Refusal, SectionClosure,
    StripRegime, StripSpec, geometry_hash, refuse, solve_nonlinear_impl,
};
use fs_blake3::hash_domain;
use fs_math::det;
use fs_solver::{
    Globalization, LineSearchConfig, NewtonKrylovConfig, NewtonKrylovState, NonlinearProblem,
};
use std::cell::{Cell, RefCell};

pub const MAX_NEWTON_ITERATIONS: usize = 64;
pub const MAX_RESIDUAL_EVALUATIONS: usize = 600;
pub const MAX_JACOBIAN_PRODUCTS: usize = 8192;

/// Explicit local model and continuation domain. All angles are in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContinuationControls {
    /// Absolute effective incidence of the original section equation.
    pub maximum_incidence_rad: f64,
    /// Maximum change from each strip's last accepted incidence, per solve.
    pub maximum_incidence_change_rad: f64,
    /// `sum(abs(linear panel gamma))/abs(sum(linear panel gamma))` per strip.
    /// Signed row ratios are legal; unresolved cancellation is not. Multirow
    /// totals at the legacy 1e-14 m^2/s row-share floor also refuse.
    pub maximum_row_condition: f64,
    /// Minimum distance from an incidence probe to any trailing filament,
    /// divided by the probed strip chord. This does not alter vortex kernels.
    pub minimum_trailing_clearance_chords: f64,
    pub maximum_newton_iterations: usize,
    pub maximum_residual_evaluations: usize,
    pub maximum_jacobian_products: usize,
    /// Centered numerical section-slope stencil, in radians. Geometry and
    /// residual-normalization derivatives are composed analytically.
    pub section_difference_rad: f64,
}

impl ContinuationControls {
    pub fn admit(self) -> Result<(), Refusal> {
        if !self.maximum_incidence_rad.is_finite()
            || self.maximum_incidence_rad <= 0.0
            || self.maximum_incidence_rad >= core::f64::consts::FRAC_PI_2
            || !self.maximum_incidence_change_rad.is_finite()
            || self.maximum_incidence_change_rad <= 0.0
            || self.maximum_incidence_change_rad > self.maximum_incidence_rad
            || !self.maximum_row_condition.is_finite()
            || self.maximum_row_condition < 1.0
            || !self.minimum_trailing_clearance_chords.is_finite()
            || self.minimum_trailing_clearance_chords <= 0.0
            || !(1..=MAX_NEWTON_ITERATIONS).contains(&self.maximum_newton_iterations)
            || !(1..=MAX_RESIDUAL_EVALUATIONS).contains(&self.maximum_residual_evaluations)
            || !(1..=MAX_JACOBIAN_PRODUCTS).contains(&self.maximum_jacobian_products)
            || !self.section_difference_rad.is_finite()
            || self.section_difference_rad <= 0.0
            || self.section_difference_rad > 1e-3
        {
            return Err(refuse(
                "nonlinear-continuation-invalid",
                "invalid model domain, Newton budget or section derivative stencil".into(),
                "declare finite positive model limits, incidence below pi/2, row condition at least one and bounded solver work",
            ));
        }
        Ok(())
    }
}

/// Last accepted physical section solution on an unchanged panel/strip layout.
/// Clone this with the caller's mechanical/wake checkpoint for exact replay.
#[derive(Clone, Debug, PartialEq)]
pub struct WingContinuation {
    layout_id: String,
    gamma: Vec<f64>,
    effective_incidence_rad: Vec<f64>,
    maximum_incidence_change_rad: f64,
    maximum_row_condition: f64,
    minimum_trailing_clearance_chords: f64,
    residual_evaluations: usize,
    jacobian_products: usize,
    domain_backtracks: usize,
}

impl WingContinuation {
    #[must_use]
    pub fn gamma(&self) -> &[f64] {
        &self.gamma
    }

    #[must_use]
    pub fn effective_incidence_rad(&self) -> &[f64] {
        &self.effective_incidence_rad
    }

    #[must_use]
    pub fn maximum_incidence_change_rad(&self) -> f64 {
        self.maximum_incidence_change_rad
    }

    #[must_use]
    pub fn maximum_row_condition(&self) -> f64 {
        self.maximum_row_condition
    }

    #[must_use]
    pub fn minimum_trailing_clearance_chords(&self) -> f64 {
        self.minimum_trailing_clearance_chords
    }

    #[must_use]
    pub fn residual_evaluations(&self) -> usize {
        self.residual_evaluations
    }

    #[must_use]
    pub fn jacobian_products(&self) -> usize {
        self.jacobian_products
    }

    #[must_use]
    pub fn domain_backtracks(&self) -> usize {
        self.domain_backtracks
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ContinuedReport {
    pub solution: NonlinearReport,
    pub continuation: WingContinuation,
}

/// Solve with an explicit accepted branch anchor, or establish an initial
/// anchor with `previous=None`. The existing fs-solver Newton--Krylov/Armijo
/// owner solves the SAME normalized circulation residual with the original
/// `GAMMA_TOL` final gate. Trial points outside the incidence domain are
/// rejected by native backtracking; budgets and source failures refuse.
///
/// The geometry/ownership must be unchanged between anchors; freestream,
/// local incident flow, wash and strip twists may change. Rebuilding the
/// influence operator remains mandatory whenever its own identity changes.
#[allow(clippy::too_many_arguments)]
pub fn solve_with_incident(
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    freestream: [f64; 3],
    rho: f64,
    closure: SectionClosure<'_>,
    incident: &[[f64; 3]],
    axial_wash: Option<&[f64]>,
    previous: Option<&WingContinuation>,
    domain: ContinuationControls,
) -> Result<ContinuedReport, Refusal> {
    domain.admit()?;
    if !rho.is_finite() || rho <= 0.0 {
        return Err(refuse(
            "nonlinear-continuation-invalid",
            "finite positive air density required".into(),
            "supply density in kg/m^3",
        ));
    }
    let mut guard = Guard {
        domain,
        previous,
        layout_id: String::new(),
        incidence: vec![0.0; strips.len()],
        row_condition: 1.0,
        clearance: f64::INFINITY,
        residual_evaluations: 0,
        jacobian_products: 0,
        domain_backtracks: 0,
    };
    let solution = solve_nonlinear_impl(
        op,
        panels,
        strips,
        freestream,
        rho,
        closure,
        Some(incident),
        axial_wash,
        previous.map(WingContinuation::gamma),
        Some(&mut guard),
    )?;
    let maximum_incidence_change_rad = previous.map_or(0.0, |anchor| {
        guard
            .incidence
            .iter()
            .zip(&anchor.effective_incidence_rad)
            .map(|(now, old)| (now - old).abs())
            .fold(0.0, f64::max)
    });
    let continuation = WingContinuation {
        layout_id: guard.layout_id,
        gamma: solution.gamma.clone(),
        effective_incidence_rad: guard.incidence,
        maximum_incidence_change_rad,
        maximum_row_condition: guard.row_condition,
        minimum_trailing_clearance_chords: guard.clearance,
        residual_evaluations: guard.residual_evaluations,
        jacobian_products: guard.jacobian_products,
        domain_backtracks: guard.domain_backtracks,
    };
    Ok(ContinuedReport {
        solution,
        continuation,
    })
}

pub(super) struct Guard<'a> {
    domain: ContinuationControls,
    previous: Option<&'a WingContinuation>,
    layout_id: String,
    incidence: Vec<f64>,
    row_condition: f64,
    clearance: f64,
    residual_evaluations: usize,
    jacobian_products: usize,
    domain_backtracks: usize,
}

impl Guard<'_> {
    pub(super) fn prepare(
        &mut self,
        panels: &[Panel],
        strips: &[StripSpec],
        freestream: [f64; 3],
        linear: &[f64],
    ) -> Result<(), Refusal> {
        let mut identity = geometry_hash(panels, [0.0; 3]).into_bytes();
        for strip in strips {
            identity.extend_from_slice(&strip.chord_m.to_bits().to_le_bytes());
            identity.extend_from_slice(&(strip.panel_indices.len() as u64).to_le_bytes());
            for &j in &strip.panel_indices {
                identity.extend_from_slice(&(j as u64).to_le_bytes());
            }
        }
        self.layout_id =
            hash_domain("org.frankensim.fs-wing.continuation-layout.v1", &identity).to_hex();
        if let Some(previous) = self.previous
            && (previous.layout_id != self.layout_id
                || previous.effective_incidence_rad.len() != strips.len())
        {
            return Err(refuse(
                "nonlinear-continuation-stale",
                "the accepted circulation belongs to another panel/strip layout".into(),
                "establish a new anchor or supply a verified circulation transfer after geometry changes",
            ));
        }
        let speed = det::sqrt(freestream.iter().map(|x| x * x).sum());
        let stream = freestream.map(|x| x / speed);
        for (s, strip) in strips.iter().enumerate() {
            let total: f64 = strip.panel_indices.iter().map(|&j| linear[j]).sum();
            let absolute: f64 = strip.panel_indices.iter().map(|&j| linear[j].abs()).sum();
            let condition = if strip.panel_indices.len() == 1 {
                1.0
            } else {
                absolute / total.abs()
            };
            // The legacy row law substitutes equal shares at |total|<=1e-14.
            // An explicit continuation must not turn that fallback into a
            // supposedly resolved signed distribution. Single-row ownership
            // is exactly one share even at zero circulation.
            if !condition.is_finite()
                || condition > self.domain.maximum_row_condition
                || (strip.panel_indices.len() > 1 && total.abs() <= ROW_SUM_FLOOR)
            {
                return Err(refuse(
                    "nonlinear-row-split-ill-conditioned",
                    format!(
                        "strip {s} signed row-sum condition {condition:e}, total {total:e} m^2/s"
                    ),
                    "the linear row distribution is unresolved in the declared domain; do not clip ratios or reuse a stale operator",
                ));
            }
            self.row_condition = self.row_condition.max(condition);
            let front = &panels[strip.panel_indices[0]];
            let probe: [f64; 3] = core::array::from_fn(|i| 0.5 * (front.a[i] + front.b[i]));
            for panel in panels {
                for endpoint in [panel.a, panel.b] {
                    let offset: [f64; 3] = core::array::from_fn(|i| probe[i] - endpoint[i]);
                    let downstream = -offset.iter().zip(stream).map(|(a, b)| a * b).sum::<f64>();
                    let closest: [f64; 3] =
                        core::array::from_fn(|i| offset[i] + downstream.max(0.0) * stream[i]);
                    let clearance = det::sqrt(closest.iter().map(|x| x * x).sum()) / strip.chord_m;
                    if !clearance.is_finite()
                        || clearance < self.domain.minimum_trailing_clearance_chords
                    {
                        return Err(refuse(
                            "nonlinear-filament-proximity",
                            format!("strip {s} trailing-filament clearance {clearance:e} chords"),
                            "the straight singular-wake model is outside the declared incidence-probe clearance domain",
                        ));
                    }
                    self.clearance = self.clearance.min(clearance);
                }
            }
        }
        Ok(())
    }

    pub(super) fn allows(&self, strip: usize, incidence: f64) -> bool {
        incidence.is_finite()
            && incidence.abs() <= self.domain.maximum_incidence_rad
            && self.previous.is_none_or(|anchor| {
                (incidence - anchor.effective_incidence_rad[strip]).abs()
                    <= self.domain.maximum_incidence_change_rad
            })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve(
        &mut self,
        initial: &[f64],
        strips: &[StripSpec],
        ratios: &[Vec<f64>],
        downwash: &[Vec<f64>],
        flow: &[(f64, f64, f64)],
        closure: SectionClosure<'_>,
        panel_count: usize,
    ) -> Result<Solved, Refusal> {
        let problem = CirculationProblem {
            guard: self,
            strips,
            ratios,
            downwash,
            flow,
            closure,
            scratch: RefCell::new(Scratch {
                panels: vec![0.0; panel_count],
                direction: vec![0.0; panel_count],
                incidence: vec![0.0; strips.len()],
                regimes: vec![StripRegime::Attached; strips.len()],
            }),
            residuals: Cell::new(0),
            products: Cell::new(0),
            domain_backtracks: Cell::new(0),
            failure: Cell::new(None),
        };
        let config = NewtonKrylovConfig {
            absolute_tolerance: 0.5 * GAMMA_TOL,
            relative_tolerance: 0.0,
            linear_restart: strips.len().min(32),
            max_linear_cycles: 4,
            globalization: Globalization::LineSearch(LineSearchConfig::default()),
            ..NewtonKrylovConfig::default()
        };
        let mut solver = NewtonKrylovState::new(&problem, initial.to_vec(), config).map_err(|error| {
            problem.failure_refusal().unwrap_or_else(|| refuse(
                "nonlinear-continuation-outside-domain", error.to_string(),
                "the initial circulation is outside the declared incidence domain; supply a closer accepted anchor",
            ))
        })?;
        let run =
            solver.run_controlled(&problem, self.domain.maximum_newton_iterations, &mut || {
                problem.failure.get().is_some()
            });
        if let Some(error) = problem.failure_refusal() {
            return Err(error);
        }
        if !run.report.converged {
            return Err(refuse(
                "nonlinear-continuation-unresolved",
                format!(
                    "native Newton {:?}, {:?}, {} attempts, residual {:e}",
                    run.progress,
                    run.report.diagnosis,
                    run.report.iterations,
                    run.report.residual_norm
                ),
                "reduce the physical continuation step or refine the declared model; no alternate root or linear fallback was returned",
            ));
        }
        // Independent of the driver's convergence flag, reevaluate the actual
        // returned circulation with the original componentwise physical gate.
        let mut final_residual = vec![0.0; strips.len()];
        problem.residual(&solver.x, &mut final_residual);
        if let Some(error) = problem.failure_refusal() {
            return Err(error);
        }
        let residual = final_residual.iter().map(|x| x.abs()).fold(0.0, f64::max);
        if final_residual.iter().any(|x| !x.is_finite()) || residual >= GAMMA_TOL {
            return Err(refuse(
                "nonlinear-continuation-unresolved",
                format!("actual final componentwise residual {residual:e}"),
                "the unchanged section equations must satisfy the original 1e-10 gate",
            ));
        }
        let scratch = problem.scratch.borrow();
        let incidence = scratch.incidence.clone();
        let regimes = scratch.regimes.clone();
        let counts = (
            problem.residuals.get(),
            problem.products.get(),
            problem.domain_backtracks.get(),
        );
        drop(scratch);
        drop(problem);
        self.incidence = incidence;
        self.residual_evaluations = counts.0;
        self.jacobian_products = counts.1;
        self.domain_backtracks = counts.2;
        Ok(Solved {
            gamma: solver.x,
            residual,
            regimes,
            iterations: run.report.iterations as u32,
        })
    }
}

pub(super) struct Solved {
    pub(super) gamma: Vec<f64>,
    pub(super) residual: f64,
    pub(super) regimes: Vec<StripRegime>,
    pub(super) iterations: u32,
}

struct Scratch {
    panels: Vec<f64>,
    direction: Vec<f64>,
    incidence: Vec<f64>,
    regimes: Vec<StripRegime>,
}

struct CirculationProblem<'a> {
    guard: &'a Guard<'a>,
    strips: &'a [StripSpec],
    ratios: &'a [Vec<f64>],
    downwash: &'a [Vec<f64>],
    flow: &'a [(f64, f64, f64)],
    closure: SectionClosure<'a>,
    scratch: RefCell<Scratch>,
    residuals: Cell<usize>,
    products: Cell<usize>,
    domain_backtracks: Cell<usize>,
    failure: Cell<Option<(&'static str, &'static str)>>,
}

impl CirculationProblem<'_> {
    fn fail(&self, code: &'static str, message: &'static str) {
        if self.failure.get().is_none() {
            self.failure.set(Some((code, message)));
        }
    }

    fn failure_refusal(&self) -> Option<Refusal> {
        self.failure.get().map(|(code, message)| refuse(code,
            format!("{message}; residual evaluations {}, Jacobian products {}",self.residuals.get(),self.products.get()),
            "respect the declared section domain and finite work budget; no accepted state was changed"))
    }

    fn scatter(&self, totals: &[f64], panels: &mut [f64]) {
        for (s, strip) in self.strips.iter().enumerate() {
            for (k, &j) in strip.panel_indices.iter().enumerate() {
                panels[j] = totals[s] * self.ratios[s][k];
            }
        }
    }

    fn induced(&self, s: usize, panels: &[f64]) -> f64 {
        let mut result = 0.0;
        for (j, coefficient) in self.downwash[s].iter().enumerate() {
            let term = panels[j] * coefficient;
            result += term;
        }
        result / self.flow[s].1
    }
}

impl NonlinearProblem for CirculationProblem<'_> {
    fn dimension(&self) -> usize {
        self.strips.len()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.failure.get().is_some() {
            return;
        }
        if self.residuals.get() >= self.guard.domain.maximum_residual_evaluations {
            self.fail(
                "nonlinear-continuation-budget",
                "residual evaluation budget exhausted",
            );
            return;
        }
        self.residuals.set(self.residuals.get() + 1);
        let mut scratch = self.scratch.borrow_mut();
        self.scatter(x, &mut scratch.panels);
        for s in 0..self.dimension() {
            let alpha =
                self.flow[s].0 + self.strips[s].twist_rad + self.induced(s, &scratch.panels);
            if !self.guard.allows(s, alpha) {
                // This is a rejected trial, NOT a latched source failure.
                // Native Armijo must be allowed to backtrack inside the tube.
                self.domain_backtracks.set(self.domain_backtracks.get() + 1);
                output.fill(f64::NAN);
                return;
            }
            let (cl, regime) = (self.closure)(s, alpha);
            if !cl.is_finite() {
                self.fail(
                    "nonlinear-did-not-converge",
                    "section closure returned nonfinite lift",
                );
                return;
            }
            let target = -0.5 * self.flow[s].2 * self.strips[s].chord_m * cl;
            output[s] = (x[s] - target) / (x[s].abs() + 1e-9);
            scratch.incidence[s] = alpha;
            scratch.regimes[s] = regime;
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.failure.get().is_some() {
            return;
        }
        if self.products.get() >= self.guard.domain.maximum_jacobian_products {
            self.fail(
                "nonlinear-continuation-budget",
                "Jacobian product budget exhausted",
            );
            return;
        }
        self.products.set(self.products.get() + 1);
        let mut scratch = self.scratch.borrow_mut();
        self.scatter(x, &mut scratch.panels);
        self.scatter(direction, &mut scratch.direction);
        let h = self.guard.domain.section_difference_rad;
        for s in 0..self.dimension() {
            let alpha =
                self.flow[s].0 + self.strips[s].twist_rad + self.induced(s, &scratch.panels);
            let lower_alpha = alpha - h;
            let upper_alpha = alpha + h;
            if !alpha.is_finite() || !(lower_alpha < alpha && upper_alpha > alpha) {
                self.fail(
                    "nonlinear-continuation-invalid",
                    "unrepresentable section slope stencil",
                );
                return;
            }
            if lower_alpha.abs().max(upper_alpha.abs()) > self.guard.domain.maximum_incidence_rad {
                self.fail(
                    "nonlinear-continuation-outside-domain",
                    "section slope stencil leaves the incidence domain",
                );
                return;
            }
            let cl = (self.closure)(s, alpha).0;
            let lower = (self.closure)(s, lower_alpha).0;
            let upper = (self.closure)(s, upper_alpha).0;
            let slope = (upper - lower) / (2.0 * h);
            if ![cl, lower, upper, slope].into_iter().all(f64::is_finite) {
                self.fail(
                    "nonlinear-did-not-converge",
                    "nonfinite numerical section slope",
                );
                return;
            }
            let multiplier = 0.5 * self.flow[s].2 * self.strips[s].chord_m;
            let mismatch = x[s] + multiplier * cl;
            let derivative =
                direction[s] + multiplier * slope * self.induced(s, &scratch.direction);
            let denominator = x[s].abs() + 1e-9;
            // Differentiate the ACTUAL state-dependent normalization too.
            let sign = if x[s] > 0.0 {
                1.0
            } else if x[s] < 0.0 {
                -1.0
            } else {
                0.0
            };
            output[s] = derivative / denominator
                - (mismatch / denominator) * (sign * direction[s] / denominator);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> ContinuationControls {
        ContinuationControls {
            maximum_incidence_rad: 0.6,
            maximum_incidence_change_rad: 0.5,
            maximum_row_condition: 8.0,
            minimum_trailing_clearance_chords: 0.02,
            maximum_newton_iterations: 32,
            maximum_residual_evaluations: 600,
            maximum_jacobian_products: 4096,
            section_difference_rad: 1e-6,
        }
    }

    fn guard(n: usize) -> Guard<'static> {
        Guard {
            domain: controls(),
            previous: None,
            layout_id: String::new(),
            incidence: vec![0.0; n],
            row_condition: 1.0,
            clearance: 0.1,
            residual_evaluations: 0,
            jacobian_products: 0,
            domain_backtracks: 0,
        }
    }

    #[test]
    fn complete_normalized_jacobian_matches_independent_directional_differences() {
        let guard = guard(2);
        let strips = [
            StripSpec {
                panel_indices: vec![0, 1],
                chord_m: 1.3,
                twist_rad: 0.02,
            },
            StripSpec {
                panel_indices: vec![2, 3],
                chord_m: 0.6,
                twist_rad: -0.03,
            },
        ];
        let ratios = [vec![1.2, -0.2], vec![0.4, 0.6]];
        let downwash = [
            vec![0.11, -0.04, 0.07, -0.05],
            vec![-0.08, 0.02, 0.13, 0.09],
        ];
        let flow = [(0.1, 3.0, 4.0), (-0.1, 2.0, 3.0)];
        let closure = |_: usize, alpha: f64| {
            (
                core::f64::consts::TAU * (alpha + 0.08) + 0.4 * alpha.powi(3),
                StripRegime::Attached,
            )
        };
        let problem = CirculationProblem {
            guard: &guard,
            strips: &strips,
            ratios: &ratios,
            downwash: &downwash,
            flow: &flow,
            closure: &closure,
            scratch: RefCell::new(Scratch {
                panels: vec![0.0; 4],
                direction: vec![0.0; 4],
                incidence: vec![0.0; 2],
                regimes: vec![StripRegime::Attached; 2],
            }),
            residuals: Cell::new(0),
            products: Cell::new(0),
            domain_backtracks: Cell::new(0),
            failure: Cell::new(None),
        };
        for x in [[-0.7, 0.3], [0.4, -0.8]] {
            let direction = [0.6, -0.5];
            let mut action = [0.0; 2];
            problem.jacobian_apply(&x, &direction, &mut action);
            let h = 2e-5;
            let mut plus = [0.0; 2];
            let mut minus = [0.0; 2];
            problem.residual(
                &core::array::from_fn::<_, 2, _>(|i| x[i] + h * direction[i]),
                &mut plus,
            );
            problem.residual(
                &core::array::from_fn::<_, 2, _>(|i| x[i] - h * direction[i]),
                &mut minus,
            );
            for i in 0..2 {
                let numerical = (plus[i] - minus[i]) / (2.0 * h);
                assert!(
                    (action[i] - numerical).abs() < 3e-8 * (1.0 + numerical.abs()),
                    "{} != {}",
                    action[i],
                    numerical
                );
            }
        }
        assert!(problem.failure.get().is_none());
    }

    #[test]
    fn native_armijo_backtracks_domain_rejections_then_solves_the_same_equation() {
        let mut guard = guard(1);
        let strips = [StripSpec {
            panel_indices: vec![0],
            chord_m: 1.0,
            twist_rad: 0.0,
        }];
        let closure = |_: usize, alpha: f64| {
            (
                0.8 * alpha + 10.0 * alpha.powi(3) + 0.1,
                StripRegime::Attached,
            )
        };
        let solved = guard
            .solve(
                &[-4.0],
                &strips,
                &[vec![1.0]],
                &[vec![0.1]],
                &[(0.0, 1.0, 2.0)],
                &closure,
                1,
            )
            .unwrap();
        let x = solved.gamma[0];
        // Independent scalar equation after eliminating the represented flow.
        assert!((1.08 * x + 0.01 * x.powi(3) + 0.1).abs() < 1e-11);
        assert!(solved.residual < GAMMA_TOL);
        assert!(
            guard.domain_backtracks > 0,
            "an out-of-domain trial must actually exercise the native retry path"
        );
        assert!(solved.iterations < 32);
    }
}
