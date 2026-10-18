//! The reference's energy-discrete WoolFelt contact force, prepared per step.
//!
//! `grand_piano/felt.rs::average` integrates fs-material's `WoolFelt` stress
//! exactly over the strain increment of one step; `felt.rs::solve_with_rate`
//! then finds the end overlap. This module evaluates the SAME closed forms
//! (virgin power-law envelope, power-law unload/reload through the residual
//! crush) and solves the SAME scalar equation to the SAME tolerance.
//!
//! Why it exists: the reference path spends about thirty `det::pow` calls per
//! solve, and `det::pow` runs on a software fused multiply-add in wasm
//! (measured 0.28 us per call, 9.2 us per solve). One struck trichord needs
//! several solves per string per sample. Here the start-of-step integrals are
//! computed once, each Newton iterate costs one `f64::powf`, and the iteration
//! is warm-started from the previous force.
//!
//! The material card, its committed crush state and the stored-energy books
//! remain fs-material's `WoolFelt` and `felt.rs::stored`. Equivalence of the
//! resulting motion with the reference engine is measured by
//! `examples/validate.rs`; `tests` below compare the force itself.
use crate::refmods::felt::{State, power_law_moments};
use fs_material::WoolFelt;

/// `felt.rs::stored`: recoverable compression energy density [J/m^3] on a
/// committed branch, evaluated with `powf`.
pub fn stored(law: &WoolFelt, e: f64, state: &State) -> f64 {
    let residual = law.crush_fraction * state.eps_max;
    if e <= residual || state.eps_max <= residual {
        return 0.0;
    }
    let span = state.eps_max - residual;
    state.sig_max * span / (law.q + 1.0) * ((e - residual) / span).powf(law.q + 1.0)
}

/// One contact site for one step: law, committed branch and start overlap.
#[derive(Clone, Copy, Debug)]
pub struct Site {
    p: f64,
    q: f64,
    f_ref: f64,
    eps_ref: f64,
    eps_max: f64,
    sig_max: f64,
    residual: f64,
    span: f64,
    start_strain: f64,
    primitive_start: f64,
    stored_max: f64,
    virgin_max: f64,
    thickness: f64,
    area: f64,
    maximum_m: f64,
    /// RT-0425's R_H [N s/m^p] for this site (0: no rate term), the step [s],
    /// and the start overlap raised to p.
    rate: f64,
    dt: f64,
    start_power: f64,
}

impl Site {
    pub fn new(law: &WoolFelt, state: &State, start_m: f64, thickness_m: f64, area_m2: f64) -> Self {
        let residual = law.crush_fraction * state.eps_max;
        let span = state.eps_max - residual;
        let mut site = Self {
            p: law.p,
            q: law.q,
            f_ref: law.f_ref,
            eps_ref: law.eps_ref,
            eps_max: state.eps_max,
            sig_max: state.sig_max,
            residual,
            span,
            start_strain: start_m / thickness_m,
            primitive_start: 0.0,
            stored_max: 0.0,
            virgin_max: 0.0,
            thickness: thickness_m,
            area: area_m2,
            maximum_m: thickness_m * law.eps_densify,
            rate: 0.0,
            dt: 1.0,
            start_power: if start_m <= 0.0 { 0.0 } else { start_m.powf(law.p) },
        };
        site.stored_max = if span > 0.0 { state.sig_max * span / (law.q + 1.0) } else { 0.0 };
        site.virgin_max = site.virgin(state.eps_max);
        site.primitive_start = site.primitive(site.start_strain);
        site
    }

    /// Add RT-0425's Kelvin-Voigt term R_H d(e^p)/dt over a step of `dt_s`
    /// (`felt.rs::average_with_rate`); R_H = 0 leaves the elastic force alone.
    pub fn with_rate(mut self, rate_n_s_m_p: f64, dt_s: f64) -> Self {
        self.rate = rate_n_s_m_p;
        self.dt = dt_s;
        self
    }

    /// Integral of the virgin envelope from zero strain.
    #[inline]
    fn virgin(&self, e: f64) -> f64 {
        if e <= 0.0 { 0.0 } else { self.f_ref * (e / self.eps_ref).powf(self.p) * e / (self.p + 1.0) }
    }

    /// `felt.rs::stored` on the committed branch.
    #[inline]
    fn stored(&self, e: f64) -> f64 {
        if e <= self.residual || self.eps_max <= self.residual {
            return 0.0;
        }
        self.sig_max * self.span / (self.q + 1.0) * ((e - self.residual) / self.span).powf(self.q + 1.0)
    }

    /// `felt.rs::primitive`: integral of the trial stress, committed branch fixed.
    #[inline]
    fn primitive(&self, e: f64) -> f64 {
        if e <= self.eps_max {
            return self.stored(e);
        }
        self.stored_max + self.virgin(e) - self.virgin_max
    }

    /// (primitive, stress) at strain `e` with ONE power evaluation.
    #[inline]
    fn primitive_and_stress(&self, e: f64) -> (f64, f64) {
        if e >= self.eps_max {
            if e <= 0.0 {
                return (if e <= self.eps_max { self.stored(e) } else { self.stored_max - self.virgin_max }, 0.0);
            }
            let stress = self.f_ref * (e / self.eps_ref).powf(self.p);
            let primitive = if e <= self.eps_max {
                // e == eps_max exactly: the committed branch ends here.
                self.stored(e)
            } else {
                self.stored_max + stress * e / (self.p + 1.0) - self.virgin_max
            };
            return (primitive, stress);
        }
        if e <= self.residual {
            return (0.0, 0.0);
        }
        let x = (e - self.residual) / self.span;
        let xq = x.powf(self.q);
        (self.sig_max * self.span / (self.q + 1.0) * xq * x, self.sig_max * xq)
    }

    /// `WoolFelt::tangent` (only for the degenerate tiny-increment branch).
    fn tangent(&self, e: f64) -> f64 {
        if e >= self.eps_max {
            if e <= 0.0 {
                return 0.0;
            }
            return self.f_ref * self.p * (e / self.eps_ref).powf(self.p - 1.0) / self.eps_ref;
        }
        if e <= self.residual {
            return 0.0;
        }
        let x = (e - self.residual) / self.span;
        self.sig_max * self.q * x.powf(self.q - 1.0) / self.span
    }

    /// Stable same-branch mean and weighted tangent, with one native power.
    /// Coincident, crossing and distant intervals retain their existing path.
    #[inline]
    fn same_branch_average(&self, a: f64, b: f64) -> Option<(f64, f64)> {
        if a == b { return None; }
        let (origin, scale, amplitude, exponent) = if a >= self.eps_max
            && b >= self.eps_max && a > 0.0 && b > 0.0 {
            (0.0, self.eps_ref, self.f_ref, self.p)
        } else if a > self.residual && b > self.residual && a <= self.eps_max && b < self.eps_max {
            // A start at the committed maximum belongs to this unloading
            // interval. Keep END-at-maximum stress conventions on the fallback.
            (self.residual, self.span, self.sig_max, self.q)
        } else { return None; };
        let mid = 0.5 * (a - origin) + 0.5 * (b - origin);
        if mid < f64::MIN_POSITIVE { return None; }
        let t = 0.5 * (b - a) / mid;
        let (mean, tangent) = power_law_moments(exponent, t)?;
        let ratio = mid / scale;
        let power = ratio.powf(exponent - 1.0);
        // Computing the tangent independently of the stress preserves it even
        // when the force underflows on a tiny positive compression interval.
        let sigma = amplitude * power * ratio;
        let derivative = amplitude * exponent * power / scale;
        Some((sigma * mean, 0.5 * derivative * tangent))
    }

    /// Held force [N] and its derivative with respect to END overlap [N/m]:
    /// `felt.rs::average_with_rate` with the start integral prepared. The rate
    /// difference is passive because e^p is monotone; on the unilateral branch a
    /// net tensile force is dropped (its unrecovered potential books as loss).
    #[inline]
    pub fn force(&self, end_m: f64) -> (f64, f64) {
        let (elastic, slope) = self.elastic_force(end_m);
        if self.rate == 0.0 {
            return (elastic, slope);
        }
        let end_power = if end_m <= 0.0 { 0.0 } else { end_m.powf(self.p) };
        let total = elastic + self.rate * (end_power - self.start_power) / self.dt;
        if total <= 0.0 {
            return (0.0, 0.0);
        }
        // R_H p e^(p-1) / dt from the same power.
        let tangent = if end_m > 0.0 { self.rate * self.p * end_power / end_m / self.dt } else { 0.0 };
        (total, slope + tangent)
    }

    /// `felt.rs::average`: the conservative part of `force`.
    #[inline]
    fn elastic_force(&self, end_m: f64) -> (f64, f64) {
        let a = self.start_strain;
        let b = end_m / self.thickness;
        let de = b - a;
        let (stress, slope) = if let Some(pair) = self.same_branch_average(a, b) {
            pair
        } else if de.abs() < 1.0e-7 * (a.abs().max(b.abs()) + 1.0e-12) {
            let mid = 0.5 * (a + b);
            (self.primitive_and_stress(mid).1, 0.5 * self.tangent(mid))
        } else {
            let (primitive, end_stress) = self.primitive_and_stress(b);
            let stress = (primitive - self.primitive_start) / de;
            (stress, (end_stress - stress) / de)
        };
        (self.area * stress, self.area * slope / self.thickness)
    }

    /// Solve `F = force(free - compliance*F)` for the end overlap, keeping a
    /// bracket and stopping at the overlap residual's own rounding floor (16 ulps
    /// of its terms) or a bracket collapsed to adjacent doubles, exactly as
    /// `felt.rs::solve_with_rate`; `warm_force` only seeds the iteration.
    pub fn solve(&self, free_m: f64, compliance_m_n: f64, warm_force: f64) -> Result<f64, &'static str> {
        let residual = |end: f64| {
            let (f, df) = self.force(end);
            (end + compliance_m_n * f - free_m, 1.0 + compliance_m_n * df, f)
        };
        let floor = |end: f64, force: f64|
            16.0 * f64::EPSILON * (end.abs() + (compliance_m_n * force).abs() + free_m.abs());
        let mut hi = free_m.min(self.maximum_m);
        let (rhi, _, fhi) = residual(hi);
        if rhi < -1.0e-12 {
            return Err("felt exceeded its admitted densification strain");
        }
        if rhi.abs() <= floor(hi, fhi) {
            return Ok(fhi);
        }
        let start_m = self.start_strain * self.thickness;
        let mut lo = start_m.min(free_m).min(0.0) - self.thickness;
        let mut bracketed = false;
        for _ in 0..16 {
            if residual(lo).0 <= 0.0 {
                bracketed = true;
                break;
            }
            lo = 2.0 * lo - self.thickness;
        }
        if !bracketed && residual(lo).0 > 0.0 {
            return Err("felt contact could not bracket its end overlap");
        }
        let seed = free_m - compliance_m_n * warm_force;
        let mut end = if warm_force > 0.0 && seed > lo && seed < hi { seed } else { 0.5 * (lo + hi) };
        for _ in 0..64 {
            let (r, derivative, force) = residual(end);
            if !r.is_finite() || !derivative.is_finite() || !force.is_finite() || force < 0.0 {
                return Err("nonfinite or tensile discrete felt force");
            }
            // A bracket collapsed to adjacent doubles is the root at this precision.
            if r.abs() <= floor(end, force) || hi - lo <= 4.0 * f64::EPSILON * lo.abs().max(hi.abs()) {
                return Ok(force);
            }
            if r > 0.0 { hi = end; } else { lo = end; }
            let newton = end - r / derivative;
            end = if newton > lo && newton < hi { newton } else { 0.5 * (lo + hi) };
        }
        Err("felt contact exhausted its safeguarded solve budget")
    }
}
