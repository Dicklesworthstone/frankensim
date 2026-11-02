//! Local sliding during an original fixed-normal, monotone scalar advance.
//!
//! Within the unchanged patch, `ds = du_macro - (u_elastic_after-u_elastic_before)`.
//! The positive Abel displacement kernel preserves the sign of each spring's
//! plastic increment, so local sliding cannot reverse within this segment.
//! Consequently the irreversible friction work per area is `mu*p*abs(ds)`.
//! This is Eq. 4.41 of Popov, Hess and Willert (2019), equivalent to its MDR
//! work integral, Eq. 4.42: <https://doi.org/10.1007/978-3-662-58709-6_4>.
//!
//! Evaluation integrates the positive spring plastic increments directly. It
//! avoids subtracting nearly equal macroscopic and elastic displacements, which
//! would create false sliding in a sticking region. The original scalar,
//! decoupled Cattaneo--Mindlin approximation still applies. Changing normal
//! load, contact entry/exit, intervening reversals, wear evolution and thermal
//! partition need their respective owners; endpoint differences alone do not
//! resolve these paths. Friction work is not an additional wear energy sink.

use super::super::{MemoryInterval, MindlinDeresiewiczPoint, MindlinState};
use super::displacement;
use crate::{
    TangentialLaw, TangentialStep, TangentialStepInput, TriboError, checked_mul, nonnegative_finite,
};
use core::f64::consts::PI;
use fs_math::det;

/// Local irreversible motion and work at a material radius in the fixed patch.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MindlinSurfaceSliding {
    /// Unchanged compressive Hertz pressure [Pa].
    pub normal_pressure_pa: f64,
    /// Signed local plastic sliding along the first tangent [m].
    pub signed_sliding_increment_m: f64,
    /// Local sliding distance [m], equal to the absolute signed increment only
    /// because this receipt represents one monotone scalar segment.
    pub sliding_distance_m: f64,
    /// Nonnegative irreversible friction work per area [J/m²]. Its area
    /// integral equals the original step's `dissipated_j`.
    pub dissipated_energy_j_per_m2: f64,
}

/// An original accepted step and its spatial sliding observer. Construction is
/// private: arbitrary endpoint pairs cannot masquerade as a monotone path.
/// The old state is borrowed, and the original step owns its candidate history.
#[derive(Debug)]
pub struct MindlinSurfaceAdvance<'a> {
    before: &'a MindlinState,
    step: TangentialStep<MindlinState>,
    pressure_weighted_sliding_j: f64,
}

impl MindlinDeresiewiczPoint {
    /// Advance once through the original constitutive/body-compliance solver,
    /// retaining a readonly observer of local sliding over its accepted path.
    /// All original input, force, convergence, history and owner gates apply.
    /// No additional history allocation or clone is needed for the observer.
    pub fn advance_with_surface_sliding<'a>(
        &self,
        before: &'a MindlinState,
        input: &TangentialStepInput,
    ) -> Result<MindlinSurfaceAdvance<'a>, TriboError> {
        let step = self.advance(before, input)?;
        let pressure_weighted_sliding_j = if self.parameters.friction_coefficient == 0.0 {
            checked_mul(
                self.parameters.normal_force_n,
                step.slip_increment_m[0].abs(),
                "mindlin_pressure_weighted_sliding_j",
            )?
        } else {
            step.dissipated_j / self.parameters.friction_coefficient
        };
        nonnegative_finite(
            pressure_weighted_sliding_j,
            "mindlin_pressure_weighted_sliding_j",
        )?;
        Ok(MindlinSurfaceAdvance {
            before,
            step,
            pressure_weighted_sliding_j,
        })
    }
}

impl MindlinSurfaceAdvance<'_> {
    /// Original accepted candidate, force, actual slip and physical work.
    #[must_use]
    pub const fn tangential_step(&self) -> &TangentialStep<MindlinState> {
        &self.step
    }

    /// Consume the observer and return the unchanged original candidate.
    #[must_use]
    pub fn into_tangential_step(self) -> TangentialStep<MindlinState> {
        self.step
    }

    /// Exact original integral `integral p(r)*sliding_distance(r) dA` [N·m].
    /// At positive mu this equals original dissipation divided by mu; at mu=0
    /// all motion slides freely and it equals `N*abs(actual_du)` instead.
    /// This is a load-distance measure, not heat at a frictionless interface.
    #[must_use]
    pub const fn pressure_weighted_sliding_j(&self) -> f64 {
        self.pressure_weighted_sliding_j
    }

    /// Evaluate local sliding in O(old history intervals), without allocation.
    /// The radius must lie in `[0,a]`; the edge is the zero-pressure limit.
    /// Exterior elastic displacement is not sliding of contacting material.
    pub fn surface_sliding(&self, radius_m: f64) -> Result<MindlinSurfaceSliding, TriboError> {
        let law = &self.before.owner;
        let a = law.parameters.contact_radius_m;
        nonnegative_finite(radius_m, "mindlin_sliding_radius_m")?;
        if radius_m > a {
            return Err(TriboError::OutsideApplicability {
                field: "mindlin_sliding_radius_m",
                value: radius_m,
                minimum: 0.0,
                maximum: a,
            });
        }
        let increment = self.step.slip_increment_m[0];
        let distance = if law.capacity_n == 0.0 || increment == 0.0 {
            increment.abs()
        } else {
            plastic_displacement(self.before, increment, radius_m)?
        };
        let radius = radius_m / a;
        let pressure = if radius_m == a {
            0.0
        } else {
            ((law.parameters.normal_force_n / a) * (1.5 / PI) / a)
                * det::sqrt((1.0 - radius) * (1.0 + radius))
        };
        nonnegative_finite(pressure, "mindlin_sliding_pressure_pa")?;
        let work = checked_mul(
            law.parameters.friction_coefficient,
            checked_mul(pressure, distance, "mindlin_local_load_distance_j_per_m2")?,
            "mindlin_local_dissipation_j_per_m2",
        )?;
        Ok(MindlinSurfaceSliding {
            normal_pressure_pa: pressure,
            signed_sliding_increment_m: increment.signum() * distance,
            sliding_distance_m: distance,
            dissipated_energy_j_per_m2: work,
        })
    }
}

fn plastic_displacement(
    before: &MindlinState,
    increment_m: f64,
    radius_m: f64,
) -> Result<f64, TriboError> {
    let law = &before.owner;
    let increment = increment_m / law.critical_slip_m;
    let direction = increment.signum();
    let mut lower = 0.0;
    let mut plastic = 0.0;
    for old in &before.memory {
        // Original clipping event: |dplastic|/d=max(0,A-B*r), B>=0.
        let rate = 1.0 - direction * old.slope;
        let shifted = direction * (old.offset + increment);
        let crossing = if rate == 0.0 {
            if shifted > 0.0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            }
        } else {
            shifted / rate
        };
        let upper = crossing.clamp(lower, old.upper);
        if upper > lower {
            let at_upper = if crossing <= old.upper {
                0.0
            } else {
                shifted - rate * upper
            };
            nonnegative_finite(at_upper, "mindlin_local_plastic_increment")?;
            // Zero outside this old interval; the middle polynomial carries
            // only positive plastic motion. Stack storage uses the same analytic
            // inverse Abel evaluator as the elastic surface-field owner.
            let pieces = [
                MemoryInterval {
                    upper: lower,
                    slope: 0.0,
                    offset: 0.0,
                },
                MemoryInterval {
                    upper,
                    slope: -rate,
                    offset: rate * upper + at_upper,
                },
                MemoryInterval {
                    upper: 1.0,
                    slope: 0.0,
                    offset: 0.0,
                },
            ];
            let start = usize::from(lower == 0.0);
            let end = if upper == 1.0 { 2 } else { 3 };
            plastic += displacement(
                &pieces[start..end],
                radius_m,
                law.parameters.contact_radius_m,
            );
        }
        lower = old.upper;
    }
    let distance = law.critical_slip_m * plastic;
    nonnegative_finite(distance, "mindlin_local_sliding_distance_m")?;
    Ok(distance)
}
