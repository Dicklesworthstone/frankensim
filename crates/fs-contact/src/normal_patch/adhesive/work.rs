//! Exact polynomial work of the original displacement-controlled JKR branch.
use super::{JkrError, SourceJkrSphereContact};

impl SourceJkrSphereContact {
    /// Mean signed force along an admitted attached approach interval, N.
    /// Both endpoints retain the original source and numerical admission. The
    /// fold is admitted without requiring its singular tangent. Equal endpoints
    /// return the original force; reversing a path leaves the mean unchanged.
    ///
    /// In t = sqrt(a/a_s)-1, F is degree six and d(delta)/dt is degree three.
    /// Five Gauss--Legendre nodes therefore integrate their product exactly in
    /// real arithmetic. The ratio of the force and displacement moments avoids
    /// subtracting nearly equal radii or stored energies. This is a bounded,
    /// allocation-free floating-point evaluation, not an interval certificate.
    pub fn attached_path_average(&self, from_m: f64, to_m: f64) -> Result<f64, JkrError> {
        let from_t = self.attached_offset(from_m)?;
        let from = self.response_at_offset(from_t, from_m)?;
        if from_m == to_m {
            return Ok(from.normal_force_n);
        }
        let to_t = self.attached_offset(to_m)?;
        self.response_at_offset(to_t, to_m)?;
        if from_t == to_t {
            return Ok(from.normal_force_n);
        }
        let lo = from_t.min(to_t);
        let hi = from_t.max(to_t);
        let midpoint = lo.midpoint(hi);
        let half = 0.5 * (hi - lo);
        // d(delta)/dt = 4 h t(t^2+3t+3). Divide by its upper-end value,
        // retaining a well-scaled denominator even one ulp from the fold.
        let upper_factor = hi * hi + 3.0 * hi + 3.0;
        let mut force_moment = 0.0;
        let mut displacement_moment = 0.0;
        // Half weights sum to one; a finite force cannot overflow solely from
        // the quadrature's conventional weight sum of two.
        for (node, weight) in [
            (-0.906_179_845_938_664, 0.118_463_442_528_094_54),
            (-0.538_469_310_105_683_1, 0.239_314_335_249_683_24),
            (0.0, 0.284_444_444_444_444_44),
            (0.538_469_310_105_683_1, 0.239_314_335_249_683_24),
            (0.906_179_845_938_664, 0.118_463_442_528_094_54),
        ] {
            let t = midpoint + half * node;
            let jacobian = (t / hi) * ((t * t + 3.0 * t + 3.0) / upper_factor);
            let x = 1.0 + t;
            let force = self.force_from_cubed_radius_parameter(x * x * x);
            let measure = weight * jacobian;
            force_moment += measure * force;
            displacement_moment += measure;
        }
        let mean = force_moment / displacement_moment;
        if !mean.is_finite() || displacement_moment <= 0.0 {
            return Err(JkrError::Unrepresentable {
                field: "attached path mean force",
            });
        }
        Ok(mean)
    }

    /// Derivative of the attached path's mean force with respect to `to_m`,
    /// holding `from_m` fixed, N/m. Both endpoints retain the original source
    /// and numerical admission. At coincident endpoints the limit is half the
    /// force tangent; the coincident fold has no finite derivative.
    ///
    /// For unequal endpoints this is `(F(to) - mean) / (to - from)`. Factored
    /// polynomial differences cancel the interval width before the same exact
    /// five-node quadrature, retaining finite derivatives for paths ending at
    /// the fold without subtracting nearly equal forces or displacements.
    pub fn attached_path_end_derivative(&self, from_m: f64, to_m: f64) -> Result<f64, JkrError> {
        let from_t = self.attached_offset(from_m)?;
        let from = self.response_at_offset(from_t, from_m)?;
        if from_m == to_m {
            return from
                .tangent_n_per_m
                .map(|tangent| 0.5 * tangent)
                .ok_or(JkrError::SingularTangent);
        }
        let to_t = self.attached_offset(to_m)?;
        self.response_at_offset(to_t, to_m)?;
        if from_t == to_t {
            return from
                .tangent_n_per_m
                .map(|tangent| 0.5 * tangent)
                .ok_or(JkrError::SingularTangent);
        }

        let lo = from_t.min(to_t);
        let hi = from_t.max(to_t);
        let midpoint = lo.midpoint(hi);
        let half = 0.5 * (hi - lo);
        let upper_factor = hi * hi + 3.0 * hi + 3.0;
        // (delta(to_t) - delta(from_t)) / (h * (to_t - from_t)).
        let from_t2 = from_t * from_t;
        let to_t2 = to_t * to_t;
        let displacement_factor = 6.0 * (from_t + to_t)
            + 4.0 * (from_t2 + from_t * to_t + to_t2)
            + (from_t + to_t) * (from_t2 + to_t2);
        let to_x = 1.0 + to_t;
        let to_x2 = to_x * to_x;
        let to_x3 = to_x2 * to_x;
        let mut derivative_moment = 0.0;
        let mut displacement_moment = 0.0;
        for (node, weight) in [
            (-0.906_179_845_938_664, 0.118_463_442_528_094_54),
            (-0.538_469_310_105_683_1, 0.239_314_335_249_683_24),
            (0.0, 0.284_444_444_444_444_44),
            (0.538_469_310_105_683_1, 0.239_314_335_249_683_24),
            (0.906_179_845_938_664, 0.118_463_442_528_094_54),
        ] {
            let t = midpoint + half * node;
            let measure = weight * (t / hi) * ((t * t + 3.0 * t + 3.0) / upper_factor);
            let x = 1.0 + t;
            let x2 = x * x;
            let x3 = x2 * x;
            // (to_t - t) / (to_t - from_t), including a reversed path.
            let end_fraction = 0.5
                * if to_t > from_t {
                    1.0 - node
                } else {
                    1.0 + node
                };
            // F(to_t)-F(t) = s*(to_t-t)*(to_x^2+to_x*x+x^2)
            //                           *(to_x^3+x^3-6).
            derivative_moment +=
                measure * end_fraction * (to_x2 + to_x * x + x2) * (to_x3 + x3 - 6.0);
            displacement_moment += measure;
        }
        let derivative = (self.force_scale_n / self.displacement_scale_m)
            * (derivative_moment / displacement_moment / displacement_factor);
        if !derivative.is_finite() || displacement_moment <= 0.0 || displacement_factor <= 0.0 {
            return Err(JkrError::Unrepresentable {
                field: "attached path end derivative",
            });
        }
        Ok(derivative)
    }

    /// Signed reversible work on the original attached branch, J. This does
    /// not cross a snap or supply its released energy; use `advance_monotone`
    /// for an explicitly prescribed path with attachment/detachment history.
    pub fn reversible_work_between(&self, from_m: f64, to_m: f64) -> Result<f64, JkrError> {
        let mean = self.attached_path_average(from_m, to_m)?;
        let work = mean * (to_m - from_m);
        if !work.is_finite() || (work == 0.0 && mean != 0.0 && from_m != to_m) {
            return Err(JkrError::Unrepresentable {
                field: "attached path work",
            });
        }
        Ok(work)
    }
}
