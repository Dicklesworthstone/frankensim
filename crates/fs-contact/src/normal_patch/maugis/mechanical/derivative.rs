//! Original mean-force derivatives in the regular contact-radius coordinate.
use super::{MaugisAttachedPotential, MaugisContinuationError, MaugisError, PI};
use crate::normal_patch::maugis::kernel;

impl MaugisAttachedPotential {
    /// Derivative of the original path mean force with respect to `to_m`, N/m,
    /// holding `from_m` fixed. No source perturbation or force subtraction is
    /// used, including for tiny motion and paths ending at a contact fold.
    ///
    /// Contact derivatives integrate the analytic radius-coordinate slopes.
    /// Successive five-point Gaussian refinements use the original work
    /// budget: its absolute tolerance is scaled by the squared source approach
    /// scale, and its relative tolerance and evaluation limit are unchanged.
    /// Each quadrature sample is one explicit kernel evaluation, not an inverse
    /// contact solve. This is bounded quadrature, not an interval certificate.
    /// The noncontact disk is closed form; a crossed zero-area junction retains
    /// both original branches. A coincident terminal fold has no finite slope.
    pub fn path_end_derivative(
        &self,
        from_m: f64,
        to_m: f64,
    ) -> Result<f64, MaugisContinuationError> {
        self.evaluate(from_m)?;
        self.evaluate(to_m)?;
        let lo = from_m.min(to_m);
        let hi = from_m.max(to_m);
        let slopes = if let Some(limits) = self.noncontact
            && lo < limits.first_contact_indentation_m
        {
            let touch = limits.first_contact_indentation_m;
            let disk = self.disk_slopes(lo, hi.min(touch), limits.minimum_indentation_m)?;
            if hi <= touch {
                disk
            } else {
                let contact = self.contact_slopes(touch, hi)?;
                let span = hi - lo;
                let left = (touch - lo) / span;
                let right = (hi - touch) / span;
                // Integrate the affine-path weight over the two physical
                // phases. The cross moment uses their resolved force secants,
                // never a subtraction of almost equal endpoint forces.
                Slopes {
                    forward: left * left * disk.forward
                        + right * right * contact.forward
                        + left * right * contact.secant(),
                    reverse: left * left * disk.reverse
                        + right * right * contact.reverse
                        + left * right * disk.secant(),
                }
            }
        } else {
            self.contact_slopes(lo, hi)?
        };
        let value = if to_m >= from_m {
            slopes.forward
        } else {
            slopes.reverse
        } * (self.source.force_scale_n / self.source.indentation_scale_m);
        if !value.is_finite() {
            return Err(MaugisError::Unrepresentable {
                field: "Maugis path endpoint derivative",
            }
            .into());
        }
        Ok(value)
    }

    fn disk_slopes(&self, lo: f64, hi: f64, minimum: f64) -> Result<Slopes, MaugisError> {
        let scale = self.source.indentation_scale_m;
        let x = (2.0 * ((lo - minimum) / scale)).sqrt();
        let y = (2.0 * ((hi - minimum) / scale)).sqrt();
        let sum = x + y;
        if sum == 0.0 {
            return Err(MaugisError::EndpointTangent);
        }
        let fold = (8.0 / (3.0 * PI)) * self.source.lambda;
        let factor = -0.5 * PI * self.source.lambda;
        // delta-delta_fold = x^2/2, F/F0 = -pi*lambda*(c_fold+x)^2/2.
        // Cancel the final-radius factor before differentiating the path
        // average; unequal paths ending at x=0 therefore stay finite.
        Ok(Slopes {
            forward: factor * (1.0 + (4.0 / 3.0) * (fold / sum) * (1.0 + x / sum)),
            reverse: factor * (1.0 + (4.0 / 3.0) * (fold / sum) * (1.0 + y / sum)),
        })
    }

    fn contact_slopes(&self, lo: f64, hi: f64) -> Result<Slopes, MaugisError> {
        let a = self.source.state_at_indentation(lo)?;
        let b = self.source.state_at_indentation(hi)?;
        if a.t == b.t {
            if a.t == self.source.minimum_t && self.source.minimum_t > 0.0 || a.ds <= 0.0 {
                return Err(MaugisError::EndpointTangent);
            }
            // At zero solid area both d(delta)/dt and dF/dt vanish, but their
            // original analytically cancelled ratio ps/ds is finite unless
            // this point is also the displacement fold.
            let half = 0.5 * (a.ps / a.ds);
            return Ok(Slopes {
                forward: half,
                reverse: half,
            });
        }
        let mid = kernel::state(self.source.lambda, a.t.midpoint(b.t))?;
        let dscale = a.dd.max(b.dd).max(mid.dd);
        let pscale = dp(a).abs().max(dp(b).abs()).max(dp(mid).abs());
        if !dscale.is_finite() || dscale <= 0.0 || !pscale.is_finite() || pscale <= 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "Maugis path slope scales",
            });
        }
        let energy_scale = self.source.force_scale_n * self.source.indentation_scale_m;
        let absolute = self.budget.work.absolute_tolerance_j / energy_scale;
        let relative = self.budget.work.relative_tolerance;
        let mut previous: Option<Slopes> = None;
        let mut panels = 1_usize;
        let mut evaluations = 0_usize;
        loop {
            let remaining = self
                .budget
                .work
                .maximum_evaluations
                .saturating_sub(evaluations);
            if panels > remaining / 30 {
                return Err(MaugisError::WorkBudgetExceeded);
            }
            let mut total = Moments::default();
            for panel in 0..panels {
                let left = panel as f64 / panels as f64;
                let width = 1.0 / panels as f64;
                let mut local = Moments::default();
                for (x, weight) in GAUSS {
                    let state =
                        kernel::state(self.source.lambda, a.t + (b.t - a.t) * (left + width * x))?;
                    let force_slope = dp(state) / pscale;
                    let mut inner = 0.0;
                    for (u, inner_weight) in GAUSS {
                        let state = kernel::state(
                            self.source.lambda,
                            a.t + (b.t - a.t) * (left + width * x * u),
                        )?;
                        inner += inner_weight * (state.dd / dscale);
                    }
                    local.approach += weight * (state.dd / dscale);
                    local.force += weight * force_slope;
                    local.forward += weight * force_slope * x * inner;
                }
                local.approach *= width;
                local.force *= width;
                local.forward *= width * width;
                total.forward += local.forward + total.approach * local.force;
                total.approach += local.approach;
                total.force += local.force;
            }
            evaluations += 30 * panels;
            // Fbar_to = integral P'(t) integral_from^t d'(s) ds dt / Delta_d^2.
            // The source-coordinate interval width cancels from numerator and
            // denominator before evaluation. Hence no tiny-force differences,
            // no nearly equal displacement subtraction, and no fold tangent.
            let forward = (pscale / dscale) * (total.forward / total.approach / total.approach);
            let secant = (pscale / dscale) * (total.force / total.approach);
            let current = Slopes {
                forward,
                reverse: secant - forward,
            };
            if total.approach <= 0.0 || !forward.is_finite() || !current.reverse.is_finite() {
                return Err(MaugisError::Unrepresentable {
                    field: "Maugis path derivative moments",
                });
            }
            if let Some(old) = previous {
                let error = (current.forward - old.forward)
                    .abs()
                    .max((current.reverse - old.reverse).abs());
                let tolerance =
                    absolute + relative * current.forward.abs().max(current.reverse.abs());
                if tolerance.is_finite() && tolerance > 0.0 && error <= tolerance {
                    return Ok(current);
                }
            }
            previous = Some(current);
            panels = panels
                .checked_mul(2)
                .ok_or(MaugisError::WorkBudgetExceeded)?;
        }
    }
}

#[derive(Clone, Copy)]
struct Slopes {
    forward: f64,
    reverse: f64,
}
impl Slopes {
    fn secant(self) -> f64 {
        self.forward + self.reverse
    }
}
#[derive(Default)]
struct Moments {
    approach: f64,
    force: f64,
    forward: f64,
}
fn dp(state: kernel::State) -> f64 {
    (state.t / state.t.hypot(1.0)) * state.ps
}
// Unit-interval nodes and weights; each pair of nested rules costs 30 original
// explicit kernel evaluations and allocates no storage.
const GAUSS: [(f64, f64); 5] = [
    (0.046_910_077_030_668, 0.118_463_442_528_094_54),
    (0.230_765_344_947_158_45, 0.239_314_335_249_683_24),
    (0.5, 0.284_444_444_444_444_44),
    (0.769_234_655_052_841_5, 0.239_314_335_249_683_24),
    (0.953_089_922_969_332, 0.118_463_442_528_094_54),
];
