//! Finite-area ideal capsules observing the original coherent Rayleigh field.
//!
//! A capsule reads the area average of p and u, before applying its common
//! `FirstOrder` pattern. This retains spatial cancellation, proximity response,
//! and all source-to-capsule travel times. Averaging point powers or magnitudes
//! would describe a different measurement. Positive quadrature areas describe
//! the sensing aperture, not a second radiating surface or an acoustic load.
//!
//! Source and capsule quadrature are separate refinement axes. The circular
//! builder is equal-area midpoint quadrature, not an exact disk integral or a
//! bandwidth certificate. There is no capsule-body scattering, diaphragm
//! mechanics, measured frequency response, moving sensor, or backreaction.

use super::{
    BaffledPressure, MAX_RAYLEIGH_INPUTS, RayleighMedium, SurfaceSample, directional, path,
};
use fs_bem::near_field::FirstOrder;
use fs_math::{dd::Dd, det};
use std::f64::consts::{PI, TAU};

/// Maximum explicitly retained capsule quadrature. It does not multiply the
/// number of runtime FIR banks: one pressure bank, at most one velocity bank.
pub const MAX_APERTURE_POINTS: usize = 4096;
/// Hard ceiling on cold source/capsule/filter coefficient contributions.
pub const MAX_KERNEL_TERMS: usize = 1 << 28;

/// A positive physical sensing-area quadrature point in the source frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    pub position_m: [f64; 3],
    pub area_m2: f64,
}

/// Immutable sensing geometry with DD-normalized physical area weights.
#[derive(Clone, Debug)]
pub struct Aperture {
    samples: Vec<Sample>,
    weights: Vec<Dd>,
    area_m2: f64,
    center_m: [f64; 3],
}

impl Aperture {
    /// Prepare a positive quadrature of the actual sensing surface. Coincident
    /// points are allowed, but negative mixing weights are not physical areas.
    /// Refinement, coverage of the surface, and spatial resolution are caller
    /// obligations. Every consumer separately admits its baffle clearance.
    pub fn new(samples: Vec<Sample>) -> Result<Self, String> {
        if !(1..=MAX_APERTURE_POINTS).contains(&samples.len())
            || samples.iter().any(|s| {
                s.position_m.iter().any(|x| !x.is_finite())
                    || !s.area_m2.is_finite()
                    || s.area_m2 <= 0.0
            })
        {
            return Err("capsule aperture needs bounded finite points with positive areas".into());
        }
        let mut area = Dd::ZERO;
        for sample in &samples {
            area = area + Dd::from_f64(sample.area_m2);
        }
        if !finite(area) || !area.to_f64().is_finite() || area.hi <= 0.0 {
            return Err("capsule aperture total area is unrepresentable".into());
        }
        let mut weights = Vec::with_capacity(samples.len());
        let mut center = [Dd::ZERO; 3];
        for sample in &samples {
            let weight = Dd::from_f64(sample.area_m2) / area;
            if !finite(weight) || weight.hi <= 0.0 {
                return Err("capsule aperture area ratio is unrepresentable".into());
            }
            for c in 0..3 {
                center[c] = center[c] + weight * Dd::from_f64(sample.position_m[c]);
            }
            weights.push(weight);
        }
        if center
            .iter()
            .any(|x| !finite(*x) || !x.to_f64().is_finite())
        {
            return Err("capsule aperture center is unrepresentable".into());
        }
        Ok(Self {
            samples,
            weights,
            area_m2: area.to_f64(),
            center_m: center.map(Dd::to_f64),
        })
    }

    /// Uniform circular sensing disk, with independent radial and azimuthal
    /// refinement. `normal` defines its plane; the microphone front axis is
    /// still supplied by `FirstOrder`. No geometry is inferred from a pattern.
    pub fn circular(
        center_m: [f64; 3],
        normal: [f64; 3],
        radius_m: f64,
        radial_rings: usize,
        angular_points: usize,
    ) -> Result<Self, String> {
        let count = radial_rings.checked_mul(angular_points);
        if center_m.iter().chain(normal.iter()).any(|x| !x.is_finite())
            || !radius_m.is_finite()
            || radius_m <= 0.0
            || !(1..=512).contains(&radial_rings)
            || !(4..=1024).contains(&angular_points)
            || count.is_none_or(|n| n > MAX_APERTURE_POINTS)
        {
            return Err(
                "circular capsule requires finite geometry and bounded explicit quadrature".into(),
            );
        }
        let scale = normal.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        if scale == 0.0 {
            return Err("circular capsule needs a nonzero plane normal".into());
        }
        let scaled = normal.map(|x| x / scale);
        let norm = det::sqrt(scaled.iter().map(|x| x * x).sum());
        let n = scaled.map(|x| x / norm);
        // Choose a well-separated axis before forming a deterministic frame.
        let mut axis = 0;
        for c in 1..3 {
            if n[c].abs() < n[axis].abs() {
                axis = c;
            }
        }
        let mut unit = [0.0; 3];
        unit[axis] = 1.0;
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let u = cross(n, unit);
        let length = det::sqrt(u.iter().map(|x| x * x).sum());
        let u = u.map(|x| x / length);
        let v = cross(n, u);
        let count = count.expect("admitted quadrature count");
        let area = PI * radius_m * radius_m;
        let cell_area = area / count as f64;
        if !area.is_finite() || !cell_area.is_finite() || cell_area <= 0.0 {
            return Err("circular capsule area is unrepresentable".into());
        }
        let mut samples = Vec::with_capacity(count);
        for ring in 0..radial_rings {
            let radius = radius_m * det::sqrt((ring as f64 + 0.5) / radial_rings as f64);
            for angle in 0..angular_points {
                let phase = TAU * (angle as f64 + 0.5) / angular_points as f64;
                let (x, y) = (radius * det::cos(phase), radius * det::sin(phase));
                samples.push(Sample {
                    position_m: std::array::from_fn(|c| center_m[c] + x * u[c] + y * v[c]),
                    area_m2: cell_area,
                });
            }
        }
        Self::new(samples)
    }

    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }
    pub fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn center_m(&self) -> [f64; 3] {
        self.center_m
    }
    pub(super) fn weights(&self) -> &[Dd] {
        &self.weights
    }
}

fn finite(x: Dd) -> bool {
    x.hi.is_finite() && x.lo.is_finite()
}

/// Cold work limit over the complete source surface and sensing aperture.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Bounds (capsule points + one template) * source points * input modes *
    /// (65 sinc taps + 2 admission taps) * field banks. The extra template is
    /// the existing point constructor, including its full source admission.
    /// This counts projection contributions, not wall time or all allocator work.
    pub maximum_kernel_terms: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_kernel_terms: 1 << 24,
        }
    }
}

impl BaffledPressure {
    /// Area-averaged ideal capsule at one mechanical/output clock.
    #[allow(clippy::too_many_arguments)]
    pub fn new_aperture_first_order(
        surface: &[SurfaceSample],
        inputs: usize,
        rate: u32,
        aperture: &Aperture,
        medium: RayleighMedium,
        microphone: FirstOrder,
        options: Options,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, String> {
        Self::new_multirate_aperture_first_order(
            surface, inputs, rate, rate, aperture, medium, microphone, options, cancelled,
        )
    }

    /// Coherently compile the complete sensing area into the original one/two
    /// causal FIR banks. Sample stepping, multirate alignment, modal outputs,
    /// and numerical latency remain owned by `BaffledPressure`; runtime work
    /// depends on the merged FIR support, not the number of capsule points.
    ///
    /// Every sensing point must satisfy the original 5 cm baffle clearance.
    /// A singleton delegates bit for bit to the existing point constructor.
    /// DD receiver weights and source products meet before rounding each FIR
    /// coefficient. The centroid is reported in `position_m`; travel delays
    /// always use the actual individual positions, never that centroid.
    /// Cancellation is checked before/after the original point construction,
    /// between sensing points and every 64 sources in the merged projection.
    /// No partial receiver escapes on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn new_multirate_aperture_first_order(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        aperture: &Aperture,
        medium: RayleighMedium,
        microphone: FirstOrder,
        options: Options,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, String> {
        let banks = if microphone.pressure_fraction() == 1.0 {
            1
        } else {
            2
        };
        let work = surface
            .len()
            .checked_mul(inputs)
            .and_then(|n| n.checked_mul(aperture.samples.len() + 1))
            .and_then(|n| n.checked_mul(path::TAPS + 2))
            .and_then(|n| n.checked_mul(banks));
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&inputs)
            || surface.is_empty()
            || !(1..=MAX_KERNEL_TERMS).contains(&options.maximum_kernel_terms)
            || work.is_none_or(|n| n > options.maximum_kernel_terms)
            || aperture.samples.iter().any(|s| s.position_m[2] < 0.05)
        {
            return Err(
                "capsule source/receiver work or 5 cm baffle clearance budget exceeded".into(),
            );
        }
        if cancelled() {
            return Err("Rayleigh aperture preparation cancelled".into());
        }
        let mut receiver = Self::new_multirate_first_order(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            aperture.samples[0].position_m,
            medium,
            microphone,
        )?;
        if cancelled() {
            return Err("Rayleigh aperture preparation cancelled".into());
        }
        if aperture.samples.len() == 1 {
            return Ok(receiver);
        }
        let alignment = receiver
            .directional_velocity
            .as_ref()
            .map_or(0.0, |reactive| {
                receiver.decimator.delay_output_frames() - reactive.decimator.delay_output_frames()
            });
        let mut acceleration = path::KernelBins::default();
        let mut velocity = path::KernelBins::default();
        let window = path::window();
        for (sample, &weight) in aperture.samples.iter().zip(&aperture.weights) {
            directional::add_paths(
                surface,
                output_rate,
                sample.position_m,
                medium,
                microphone,
                alignment,
                weight,
                &window,
                &mut acceleration,
                &mut velocity,
                &mut cancelled,
            )?;
        }
        directional::install(&mut receiver, acceleration)?;
        receiver.position_m = aperture.center_m;
        if let Some(reactive) = receiver.directional_velocity.as_mut() {
            directional::install(reactive, velocity)?;
            reactive.position_m = aperture.center_m;
            receiver.delay_samples = (
                receiver.delay_samples.0.min(reactive.delay_samples.0),
                receiver.delay_samples.1.max(reactive.delay_samples.1),
            );
        }
        if cancelled() {
            return Err("Rayleigh aperture preparation cancelled".into());
        }
        Ok(receiver)
    }
}

#[cfg(test)]
#[path = "baffled_pressure_aperture_tests.rs"]
mod tests;
