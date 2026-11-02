//! Pressure-weighted sliding during simultaneous normal/shear evolution.
//!
//! With constant local Coulomb coefficient mu>0, physical friction work is
//! `D=mu*integral p(t,r)*abs(ds(t,r)) dA`, including microslip induced by normal
//! unloading and the entry/exit of contacting material. The original analytic
//! path already integrates this work; dividing by its immutable mu yields the
//! load-distance measure without a new solve or spatial/time discretization.
//!
//! At mu=0 all tangential motion slides freely. Indentation and shear are
//! affine, so the measure is `mean_Hertz_force*abs(actual_shear_increment)`.
//! The Hertz force is NOT affine: its mean is evaluated analytically below.
//! See the original model's MDR work equivalence and arbitrary-load history in
//! Popov, Hess and Willert (2019), section 4.4, Eqs. 4.22--4.25 and 4.41--4.42:
//! <https://doi.org/10.1007/978-3-662-58709-6_4>.
//!
//! This global measure is neither macroscopic sliding distance nor normal
//! work. It does not reconstruct local varying-load slip from endpoint fields,
//! evolve wear geometry, or select any wear coefficient, hardness or heat share.

use super::{HertzMindlinContact, HertzMindlinState, MindlinPathError, MindlinPathStep};
use crate::{TangentialStepInput, TriboError, checked_mul, nonnegative_finite};
use fs_math::det;

/// Original accepted path and its integrated load-distance measure. Private
/// construction prevents arbitrary endpoints from standing in for a history.
#[derive(Debug)]
pub struct MindlinPathSlidingAdvance {
    path: MindlinPathStep,
    pressure_weighted_sliding_j: f64,
}

impl HertzMindlinContact {
    /// Resolve the original simultaneous path once and retain its load-distance
    /// measure. The original force/iteration/history budgets, body compliance,
    /// cancellation and owner checks apply unchanged. No additional history is
    /// allocated or cloned, and errors leave the accepted state untouched.
    pub fn advance_coupled_path_with_sliding(
        &self,
        before: &HertzMindlinState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinPathSlidingAdvance, MindlinPathError> {
        let path = self.advance_coupled_path(before, input, cancelled)?;
        let mu = self.parameters().friction_coefficient;
        let pressure_weighted_sliding_j = if mu > 0.0 {
            path.tangential.dissipated_j / mu
        } else {
            checked_mul(
                mean_hertz_force(before, &path.tangential.state)?,
                path.tangential.slip_increment_m[0].abs(),
                "mindlin_path_pressure_weighted_sliding_j",
            )?
        };
        nonnegative_finite(
            pressure_weighted_sliding_j,
            "mindlin_path_pressure_weighted_sliding_j",
        )?;
        Ok(MindlinPathSlidingAdvance {
            path,
            pressure_weighted_sliding_j,
        })
    }
}

impl MindlinPathSlidingAdvance {
    /// Original accepted history, actual shear travel, mechanical effort and
    /// separate normal/shear work; these are unchanged by this observer.
    #[must_use]
    pub const fn path_step(&self) -> &MindlinPathStep {
        &self.path
    }

    /// Consume the observer and return the unchanged original path candidate.
    #[must_use]
    pub fn into_path_step(self) -> MindlinPathStep {
        self.path
    }

    /// `integral p(t,r)*abs(ds(t,r)) dA` over the complete accepted path [N·m].
    /// It can be positive at zero macroscopic travel when normal unloading
    /// releases retained shear through microslip. Open-to-open motion gives
    /// zero; at mu=0 this is still load-distance, despite zero friction work.
    #[must_use]
    pub const fn pressure_weighted_sliding_j(&self) -> f64 {
        self.pressure_weighted_sliding_j
    }
}

fn mean_hertz_force(
    before: &HertzMindlinState,
    after: &HertzMindlinState,
) -> Result<f64, TriboError> {
    let n0 = before.normal_force_n();
    let n1 = after.normal_force_n();
    if n0 == n1 {
        return Ok(n0); // Exact fixed-load and open-to-open reductions.
    }
    let d0 = before.indentation_m();
    let d1 = after.indentation_m();
    let ratio = det::sqrt(d0.min(d1)) / det::sqrt(d0.max(d1));
    // Divide the fifth-power difference by the second-power difference, then
    // normalize by the larger indentation. All terms are positive and bounded;
    // no endpoint energy subtraction, unscaled high powers or 0/0 occur.
    let polynomial = 1.0 + ratio * (1.0 + ratio * (1.0 + ratio * (1.0 + ratio)));
    let shape = 0.4 * polynomial / (1.0 + ratio);
    checked_mul(n0.max(n1), shape, "mindlin_path_mean_normal_force_n")
}
