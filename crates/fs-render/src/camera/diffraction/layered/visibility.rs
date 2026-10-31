//! Conservative perspective footprints for full-pupil geometric visibility.
//! A far source theta_b crosses the near plane at
//! theta_b + pupil_xy*(1/z_near - 1/z_far). Intervals include all pupil
//! points and lateral camera positions, rather than testing just the centre.

#[derive(Clone, Copy, Debug)]
pub(crate) struct ProjectedPupilFootprint {
    near_m: f64,
    far_m: f64,
    lower: [f64; 2],
    upper: [f64; 2],
}

impl ProjectedPupilFootprint {
    /// Include coordinate-roundoff bounds before perspective division. For
    /// positive depth, every point of a triangle projects into its vertices'
    /// angular bounding rectangle, even with small depth roundoff.
    pub(crate) fn point(position: [f64; 3], roundoff_m: f64) -> Option<Self> {
        if position.iter().any(|value| !value.is_finite())
            || !roundoff_m.is_finite()
            || roundoff_m < 0.0
        {
            return None;
        }
        let near_m = (position[2] - roundoff_m).next_down();
        let far_m = (position[2] + roundoff_m).next_up();
        if near_m <= 0.0 || !far_m.is_finite() {
            return None;
        }
        let mut lower = [0.0; 2];
        let mut upper = [0.0; 2];
        for axis in 0..2 {
            let lo = (position[axis] - roundoff_m).next_down();
            let hi = (position[axis] + roundoff_m).next_up();
            let quotients = [lo / near_m, lo / far_m, hi / near_m, hi / far_m];
            lower[axis] = quotients
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down();
            upper[axis] = quotients
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up();
        }
        if lower.iter().chain(&upper).any(|value| !value.is_finite()) {
            return None;
        }
        Some(Self {
            near_m,
            far_m,
            lower,
            upper,
        })
    }

    pub(crate) fn include(&mut self, other: Self) {
        self.near_m = self.near_m.min(other.near_m);
        self.far_m = self.far_m.max(other.far_m);
        for axis in 0..2 {
            self.lower[axis] = self.lower[axis].min(other.lower[axis]);
            self.upper[axis] = self.upper[axis].max(other.upper[axis]);
        }
    }

    /// Certify separation of this nearer layer from a farther layer for all
    /// positions in the supplied pupil/trajectory rectangle. Touching or
    /// numerically unresolved bounds refuse. This is geometric clearance,
    /// not a full-wave obstacle-diffraction error enclosure.
    pub(crate) fn clears(
        self,
        farther: Self,
        pupil_lower_m: [f64; 2],
        pupil_upper_m: [f64; 2],
    ) -> bool {
        if self.far_m >= farther.near_m
            || pupil_lower_m
                .into_iter()
                .zip(pupil_upper_m)
                .any(|(lo, hi)| !lo.is_finite() || !hi.is_finite() || lo > hi)
        {
            return false;
        }
        let coefficient = [
            ((1.0 / self.far_m).next_down() - (1.0 / farther.near_m).next_up()).next_down(),
            ((1.0 / self.near_m).next_up() - (1.0 / farther.far_m).next_down()).next_up(),
        ];
        if coefficient[0] <= 0.0 || coefficient.iter().any(|value| !value.is_finite()) {
            return false;
        }
        for axis in 0..2 {
            let lo = pupil_lower_m[axis];
            let hi = pupil_upper_m[axis];
            let products = [
                lo * coefficient[0],
                lo * coefficient[1],
                hi * coefficient[0],
                hi * coefficient[1],
            ];
            let shift_lo = products
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down();
            let shift_hi = products
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up();
            let projected_lo = (farther.lower[axis] + shift_lo).next_down();
            let projected_hi = (farther.upper[axis] + shift_hi).next_up();
            if projected_lo.is_finite()
                && projected_hi.is_finite()
                && (projected_hi < self.lower[axis] || projected_lo > self.upper[axis])
            {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rectangle(depth: f64, lower: [f64; 2], upper: [f64; 2]) -> ProjectedPupilFootprint {
        let mut value =
            ProjectedPupilFootprint::point([depth * lower[0], depth * lower[1], depth], 1.0e-14)
                .unwrap();
        value.include(
            ProjectedPupilFootprint::point([depth * upper[0], depth * upper[1], depth], 1.0e-14)
                .unwrap(),
        );
        value
    }

    #[test]
    fn full_pupil_guard_detects_hidden_edge_clipping_and_signed_lateral_motion() {
        let near = rectangle(1.0, [-0.002; 2], [0.002; 2]);
        let partial = rectangle(2.0, [0.0021, -0.001], [0.0025, 0.001]);
        assert!(
            near.clears(partial, [0.0; 2], [0.0; 2]),
            "centre ray alone misses the foreground"
        );
        assert!(
            !near.clears(partial, [-0.001; 2], [0.001; 2]),
            "the outer pupil sees the occluding edge"
        );
        let separated = rectangle(2.0, [0.003, -0.001], [0.004, 0.001]);
        assert!(near.clears(separated, [-0.001; 2], [0.001; 2]));
        let left = rectangle(2.0, [-0.004, -0.001], [-0.0035, 0.001]);
        assert!(near.clears(left, [-0.001; 2], [0.001; 2]));
        assert!(
            !near.clears(left, [0.003, -0.001], [0.005, 0.001]),
            "positive camera translation sweeps the far footprint into the foreground"
        );
        assert!(near.clears(left, [-0.005, -0.001], [-0.003, 0.001]));
        assert!(!near.clears(near, [-0.001; 2], [0.001; 2]));
    }
}
