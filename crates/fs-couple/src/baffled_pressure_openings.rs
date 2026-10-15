//! Coplanar circular flow ports through the original Rayleigh pressure kernel.
//! Linear surface superposition retains signed flow, distance and phase, not
//! a sum of magnitudes or separately normalized audio. All openings face +z.
//! The geometry is an infinite rigid baffle, not an unbaffled instrument body.
use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium, MAX_RAYLEIGH_INPUTS};
use std::collections::BTreeMap;

/// Physical duct-port binding and complete substep traces, shared by wind sources.
#[path = "baffled_pressure_duct.rs"]
pub mod duct;


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
        Self::assemble_openings(openings, receiver, 1, |radius, local| {
            Self::circular_outlet(radius, rate, local, medium)
        })
    }

    /// Full-band pressure from EVERY accepted mechanical flow sample of a disk.
    /// Reuses the original surface quadrature, bandlimited differentiator,
    /// anti-alias stages and fractional travel-time kernel. Mechanical flow is
    /// never endpoint-downsampled or converted to a fitted pressure gain.
    ///
    /// `step_trace` accepts one interleaved flow frame per mechanical substep.
    /// Clocks must have an integer ratio in 1..=16. The declared physical band
    /// may extend to 0.45 of the output rate when spatial quadrature admits it.
    /// Filter latency is reported in OUTPUT samples, separately from travel.
    /// Oversampling cannot repair an underresolved jet, tube or structural image.
    /// The original finite-difference constructor remains an independent image.
    ///
    /// # Errors
    /// Invalid clocks, geometry/medium, underresolved band or causal-history budget.
    pub fn circular_outlet_multirate_differentiating(radius_m: f64,
        mechanics_rate: u32, output_rate: u32, receiver: CircularOutletReceiver,
        medium: RayleighMedium) -> Result<Self, String>
    {
        Self::opening_clock_ratio(mechanics_rate, output_rate)?;
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err("differentiating outlet band exceeds 0.45 of its output clock".into());
        }
        let surface = Self::circular_surface(radius_m, receiver, medium)?;
        Self::new_multirate_differentiating(&surface, 1, mechanics_rate,
            output_rate, receiver.position_m, medium)
    }

    /// Multi-opening counterpart of `circular_outlet_multirate_differentiating`.
    /// Inputs are time-major: all signed outlet flows at substep 0, then all at
    /// substep 1, etc. One shared output sample includes every mechanical frame.
    /// No kernel, derivative or numerical latency is approximated a second time.
    /// Single-opening arithmetic is identical to the corresponding disk image.
    ///
    /// # Errors
    /// The original opening count/overlap/surface limits and differentiating
    /// clock/band/history admission. Failure creates no running receiver.
    pub fn circular_openings_multirate_differentiating(openings: &[CircularOpening],
        mechanics_rate: u32, output_rate: u32, receiver: CircularOutletReceiver,
        medium: RayleighMedium) -> Result<Self, String>
    {
        let ratio = Self::opening_clock_ratio(mechanics_rate, output_rate)?;
        Self::assemble_openings(openings, receiver, ratio, |radius, local| {
            Self::circular_outlet_multirate_differentiating(radius, mechanics_rate,
                output_rate, local, medium)
        })
    }

    fn opening_clock_ratio(mechanics: u32, output: u32) -> Result<usize, String> {
        if !(8_000..=192_000).contains(&output) || mechanics == 0
            || mechanics % output != 0 || !(1..=16).contains(&(mechanics / output)) {
            return Err("outlet reconstruction requires an integer mechanics/output ratio in 1..=16".into());
        }
        Ok((mechanics / output) as usize)
    }

    fn assemble_openings(openings: &[CircularOpening], receiver: CircularOutletReceiver,
        ratio: usize, mut disk_receiver: impl FnMut(f64, CircularOutletReceiver) -> Result<Self, String>)
        -> Result<Self, String>
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
            let mut disk = disk_receiver(opening.radius_m, local)?;
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
        result.set_input_ratio(ratio)?;
        result.position_m = receiver.position_m;
        result.delay_samples = (lo, hi);
        // Both the basis and delay rows changed. Rebuild the shared optimized
        // evaluator even when growing the history retained its allocation.
        result.prepare_taps();
        Ok(result)
    }
}
