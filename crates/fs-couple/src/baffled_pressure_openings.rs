//! Coplanar circular flow ports through the original Rayleigh pressure kernel.
//! Linear surface superposition retains signed flow, distance and phase, not
//! a sum of magnitudes or separately normalized audio. All openings face +z.
//! The geometry is an infinite rigid baffle, not an unbaffled instrument body.
use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium, MAX_RAYLEIGH_INPUTS};
use std::collections::BTreeMap;

/// One physical circular opening in a shared stationary z=0 baffle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularOpening {
    /// Center [m] in the common baffle plane.
    pub center_m: [f64; 2],
    /// Actual aperture radius [m]. Input flow is uniform velocity times pi r^2.
    pub radius_m: f64,
}

impl BaffledPressure {
    /// Prepare multiple flow inputs using the EXISTING circular surface rule.
    /// Each disk is translated into its local frame for cold construction, then
    /// its exact delay coefficients become one column of a single shared FIR.
    /// Every step supplies all signed outward volume flows [m^3/s] together.
    /// No new time integrator, propagation approximation or radiation load.
    ///
    /// Positions, surface resolution and homogeneous exterior gas are physical
    /// inputs. The frequency bound is not a filter. Inter-opening radiation
    /// reaction and free-field body scattering are not modeled by this observer.
    ///
    /// # Errors
    /// Empty/over-budget inputs, nonfinite/overlapping disks, more than 120000
    /// total quadrature points, or any original circular receiver refusal.
    pub fn circular_openings(openings: &[CircularOpening], rate: u32,
        receiver: CircularOutletReceiver, medium: RayleighMedium) -> Result<Self, String>
    {
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&openings.len())
            || !(1..=128).contains(&receiver.radial_rings)
            || !(8..=512).contains(&receiver.angular_points)
            || openings.len() * receiver.radial_rings * receiver.angular_points > 120_000 {
            return Err("circular opening array exceeds the shared input or surface budget".into());
        }
        for (i, opening) in openings.iter().enumerate() {
            if !opening.radius_m.is_finite() || opening.radius_m <= 0.0
                || opening.center_m.iter().any(|x| !x.is_finite()) {
                return Err("circular opening needs finite coordinates and positive radius".into());
            }
            for other in &openings[..i] {
                let dx = opening.center_m[0] - other.center_m[0];
                let dy = opening.center_m[1] - other.center_m[1];
                let distance = fs_math::det::sqrt(dx * dx + dy * dy);
                let minimum = opening.radius_m + other.radius_m;
                if !distance.is_finite() || !minimum.is_finite() || distance < minimum {
                    return Err("coplanar circular openings must not overlap".into());
                }
            }
        }
        let mut bins = BTreeMap::<usize, Vec<f64>>::new();
        let mut first = None;
        for (i, opening) in openings.iter().enumerate() {
            let local = CircularOutletReceiver { position_m: [
                receiver.position_m[0] - opening.center_m[0],
                receiver.position_m[1] - opening.center_m[1], receiver.position_m[2]], ..receiver };
            let mut disk = Self::circular_outlet(opening.radius_m, rate, local, medium)?;
            if openings.len() == 1 {
                disk.position_m = receiver.position_m;
                return Ok(disk); // Preserve the original one-port arithmetic/history exactly.
            }
            for (delay, weights) in disk.kernels.drain(..) {
                bins.entry(delay).or_insert_with(|| vec![0.0; openings.len()])[i] = weights[0];
            }
            if first.is_none() { first = Some(disk); }
        }
        // Reuse the first quiescent receiver's validated clock/derivative mode.
        // Resize only its input basis and bounded zero history before publication.
        let mut result = first.expect("nonempty openings admitted");
        let lo = *bins.first_key_value().expect("circular kernels are nonempty").0;
        let hi = *bins.last_key_value().expect("circular kernels are nonempty").0;
        result.kernels = bins.into_iter().collect();
        result.modes = openings.len();
        result.history_frames = hi + 1;
        let samples = result.history_frames.checked_mul(result.modes).ok_or("opening history size overflow")?;
        result.history.try_reserve_exact(samples.saturating_sub(result.history.len()))
            .map_err(|_| "opening history allocation failed")?;
        result.history.resize(samples, 0.0);
        result.previous_velocity.resize(result.modes, 0.0);
        result.set_input_ratio(1)?;
        result.position_m = receiver.position_m;
        result.delay_samples = (lo, hi);
        Ok(result)
    }
}
