//! Coplanar circular and rectangular flow ports through the original Rayleigh kernel.
//! Linear surface superposition retains signed flow, distance and phase, not
//! a sum of magnitudes or separately normalized audio. All openings face +z.
//! The geometry is an infinite rigid baffle, not an unbaffled instrument body.
use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium, MAX_RAYLEIGH_INPUTS};
use std::collections::BTreeMap;
use super::gas::{GasPropagation, GasPropagationReport};
use fs_material::gas::GasState;

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

#[cfg(test)]
mod rectangular_tests {
    use super::*;

    fn gain(receiver: &BaffledPressure, frequency_hz: f64) -> f64 {
        let omega = core::f64::consts::TAU * frequency_hz / receiver.rate;
        let mut real = 0.0;
        let mut imag = 0.0;
        for (delay, weights) in &receiver.kernels {
            real += weights[0] * (omega * *delay as f64).cos();
            imag -= weights[0] * (omega * *delay as f64).sin();
        }
        real.hypot(imag)
    }

    #[test]
    fn rectangular_flow_directivity_converges_to_the_independent_far_field_sinc_product() {
        // Rayleigh rectangular piston, MIT 2.067 lectures 4–6, p.6:
        // https://ocw.mit.edu/courses/2-067-advanced-structural-dynamics-and-acoustics-13-811-spring-2004/9b4456251e62368eb020ee03a8b017e2_lect_4_6.pdf
        // The radiated pressure per dQ/dt is rho/(2*pi*r) times this sinc product.
        let medium = RayleighMedium { density: 1.2, sound_speed: 343.0 };
        let size = [0.08, 0.04];
        let frequency = 3500.0;
        let radius = 10.0;
        let sinc = |x: f64| if x == 0.0 { 1.0 } else { x.sin()/x };
        for (theta, phi) in [(0.3_f64, 0.4_f64), (1.0_f64, 0.2_f64)] {
            let direction = [theta.sin()*phi.cos(), theta.sin()*phi.sin(), theta.cos()];
            let wave_number = core::f64::consts::TAU * frequency / medium.sound_speed;
            let exact = medium.density/(core::f64::consts::TAU*radius)
                * sinc(0.5*wave_number*size[0]*direction[0])
                * sinc(0.5*wave_number*size[1]*direction[1]);
            let mut errors = Vec::new();
            // Uniform midpoint cells sum exp(i*beta*x) with the exact relative
            // error (beta*h/2)/sin(beta*h/2) - 1 per axis. At theta=1, phi=0.2
            // that is 7.36e-4 at 32 cells, above 5e-4 by construction, and
            // 1.84e-4 at 64 cells (plus a ~1.3e-5 finite-distance floor).
            for cells in [8, 16, 32, 64] {
                let spec = RectangularOutletReceiver { center_m: [0.0; 2],
                    position_m: direction.map(|d| radius*d), cells: [cells; 2],
                    maximum_frequency_hz: frequency };
                let receiver = BaffledPressure::rectangular_outlet_multirate_differentiating(
                    size, 192_000, 192_000, spec, medium).unwrap();
                errors.push((gain(&receiver, frequency)/exact - 1.0).abs());
            }
            assert!(errors.windows(2).all(|pair| pair[1] < pair[0]), "{errors:?}");
            assert!(errors[3] < 5e-4, "{errors:?}");
        }
    }

    #[test]
    fn rectangular_geometry_band_and_physical_density_are_not_hidden_gains() {
        let medium = RayleighMedium { density: 1.2, sound_speed: 343.0 };
        let spec = RectangularOutletReceiver { center_m: [-0.3, 0.0],
            position_m: [0.0, 0.0, 1.0], cells: [4; 2], maximum_frequency_hz: 1000.0 };
        let a = BaffledPressure::rectangular_outlet([0.01, 0.004], 48_000, spec, medium).unwrap();
        let b = BaffledPressure::rectangular_outlet([0.01, 0.004], 48_000, spec,
            RayleighMedium { density: 2.4, ..medium }).unwrap();
        assert_eq!(gain(&b, 1000.0), 2.0*gain(&a, 1000.0));
        assert_eq!(a.delay_samples, b.delay_samples);
        for refused in [RectangularOutletReceiver { cells: [0, 4], ..spec },
            RectangularOutletReceiver { cells: [512, 512], ..spec },
            RectangularOutletReceiver { maximum_frequency_hz: 4801.0, ..spec },
            RectangularOutletReceiver { center_m: [f64::NAN, 0.0], ..spec }] {
            assert!(BaffledPressure::rectangular_outlet([0.01, 0.004], 48_000, refused, medium).is_err());
        }
        assert!(BaffledPressure::rectangular_outlet([0.0, 0.004], 48_000, spec, medium).is_err());
        assert!(BaffledPressure::rectangular_outlet([1.0, 0.004], 48_000, spec, medium).is_err());
    }
}

/// Geometry and numerical resolution for a uniform rectangular flow aperture.
/// Its edges follow x/y in the same stationary z=0 rigid baffle as disk outlets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularOutletReceiver {
    /// Aperture center [m] in the shared baffle plane.
    pub center_m: [f64; 2],
    /// Receiver position [m] in that same frame.
    pub position_m: [f64; 3],
    /// Positive midpoint quadrature counts along width and height, each <= 512.
    pub cells: [usize; 2],
    /// Resolved physical band [Hz], with the same temporal limits as disks.
    pub maximum_frequency_hz: f64,
}

impl BaffledPressure {
    /// Observe total outward volume flow through an actual rectangular aperture.
    /// Area-normalized velocity shapes retain the rectangle's directivity and
    /// travel-time differences; no equivalent circular radius or compact gain.
    /// The source is unchanged and no reaction impedance is added.
    pub fn rectangular_outlet(size_m: [f64; 2], rate: u32,
        receiver: RectangularOutletReceiver, medium: RayleighMedium) -> Result<Self, String>
    {
        if receiver.maximum_frequency_hz > 0.1 * f64::from(rate) {
            return Err("rectangular receiver band exceeds its finite-difference clock".into());
        }
        let surface = Self::rectangular_surface(size_m, receiver, medium)?;
        Self::new(&surface, 1, rate, receiver.position_m, medium)
    }

    /// The same physical rectangle with the shared differentiator, antialias
    /// stages and causal travel-time kernel. Every mechanics-frame flow enters
    /// `step_trace`; numerical latency is identical to the disk image's.
    pub fn rectangular_outlet_multirate_differentiating(size_m: [f64; 2],
        mechanics_rate: u32, output_rate: u32, receiver: RectangularOutletReceiver,
        medium: RayleighMedium) -> Result<Self, String>
    {
        Self::opening_clock_ratio(mechanics_rate, output_rate)?;
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err("rectangular receiver band exceeds its differentiating output clock".into());
        }
        let surface = Self::rectangular_surface(size_m, receiver, medium)?;
        Self::new_multirate_differentiating(&surface, 1, mechanics_rate,
            output_rate, receiver.position_m, medium)
    }

    /// The same rectangle with material transport on every surface path.
    /// The complete exterior state, signed flow, quadrature and physical delay
    /// enter the shared gas Rayleigh owner; no filter follows the pressure sum.
    #[allow(clippy::too_many_arguments)]
    pub fn rectangular_outlet_multirate_gas_differentiating(
        size_m: [f64; 2], mechanics_rate: u32, output_rate: u32,
        receiver: RectangularOutletReceiver, gas: &GasState,
        propagation: &GasPropagation, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasPropagationReport), String> {
        Self::opening_clock_ratio(mechanics_rate, output_rate)?;
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err("rectangular receiver band exceeds its differentiating output clock".into());
        }
        let surface = Self::rectangular_surface(size_m, receiver, gas.into())?;
        Self::new_multirate_gas_differentiating(&surface, 1, mechanics_rate,
            output_rate, receiver.position_m, gas, propagation, cancelled)
    }

    pub(crate) fn rectangular_surface(size_m: [f64; 2], receiver: RectangularOutletReceiver,
        medium: RayleighMedium) -> Result<Vec<super::SurfaceSample>, String>
    {
        let points = receiver.cells[0].checked_mul(receiver.cells[1]);
        if size_m.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || receiver.center_m.iter().any(|x| !x.is_finite())
            || receiver.cells.iter().any(|n| !(1..=512).contains(n))
            || points.is_none_or(|n| n > 120_000)
            || !receiver.maximum_frequency_hz.is_finite() || receiver.maximum_frequency_hz <= 0.0
            || !medium.sound_speed.is_finite() || medium.sound_speed <= 0.0 {
            return Err("rectangular receiver requires physical geometry and bounded resolution".into());
        }
        let area = size_m[0] * size_m[1];
        let dx = size_m[0] / receiver.cells[0] as f64;
        let dy = size_m[1] / receiver.cells[1] as f64;
        let phase = core::f64::consts::TAU * receiver.maximum_frequency_hz
            * fs_math::det::sqrt(dx * dx + dy * dy) / medium.sound_speed;
        if !area.is_finite() || area <= 0.0 || !(1.0 / area).is_finite()
            || !phase.is_finite() || phase > core::f64::consts::FRAC_PI_4 {
            return Err("rectangular area is unrepresentable or quadrature underresolves the acoustic band".into());
        }
        let points = points.expect("bounded quadrature");
        let mut surface = Vec::with_capacity(points);
        for y in 0..receiver.cells[1] {
            for x in 0..receiver.cells[0] {
                surface.push(super::SurfaceSample {
                    position_m: [receiver.center_m[0] + (x as f64 + 0.5) * dx - 0.5 * size_m[0],
                        receiver.center_m[1] + (y as f64 + 0.5) * dy - 0.5 * size_m[1], 0.0],
                    area_m2: area / points as f64,
                    mode_shape: vec![1.0 / area],
                });
            }
        }
        Ok(surface)
    }

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

    /// Original circular opening surfaces with absorption and dispersion from
    /// the complete exterior gas on each individual quadrature path. All disks
    /// share one physical kernel and causal padding, preserving their relative
    /// phases and signed volume-flow basis.
    #[allow(clippy::too_many_arguments)]
    pub fn circular_openings_multirate_gas_differentiating(
        openings: &[CircularOpening], mechanics_rate: u32, output_rate: u32,
        receiver: CircularOutletReceiver, gas: &GasState,
        propagation: &GasPropagation, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasPropagationReport), String> {
        Self::opening_clock_ratio(mechanics_rate, output_rate)?;
        Self::admit_openings(openings, receiver)?;
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err("differentiating outlet band exceeds 0.45 of its output clock".into());
        }
        let mut surface = Vec::new();
        for (index, opening) in openings.iter().enumerate() {
            if cancelled() { return Err("material opening propagation cancelled".into()); }
            let local = CircularOutletReceiver { position_m: [
                receiver.position_m[0] - opening.center_m[0],
                receiver.position_m[1] - opening.center_m[1], receiver.position_m[2]], ..receiver };
            for mut sample in Self::circular_surface(opening.radius_m, local, gas.into())? {
                sample.position_m[0] += opening.center_m[0];
                sample.position_m[1] += opening.center_m[1];
                let weight = sample.mode_shape[0];
                sample.mode_shape = vec![0.0; openings.len()];
                sample.mode_shape[index] = weight;
                surface.push(sample);
            }
        }
        Self::new_multirate_gas_differentiating(&surface, openings.len(),
            mechanics_rate, output_rate, receiver.position_m, gas, propagation, cancelled)
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
        Self::admit_openings(openings, receiver)?;
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

    fn admit_openings(openings: &[CircularOpening], receiver: CircularOutletReceiver)
        -> Result<(), String>
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
        Ok(())
    }
}
