//! Dynamic longitudinal strain on the existing moderate-slope string law.
//!
//! For one or two transverse displacements, the material strain is
//! `e = u_x + (y_x² + z_x²)/2`. Longitudinal displacement has fixed endpoints
//! and mass-normalized sine coordinates. A caller supplies the transverse
//! slope coefficients, including the moving-end chord, as
//! `y_x = a_0 + sum_n a_n cos(n pi x/L)` (and likewise for `z`).
//!
//! Write geometric strain as `s_0 + sum_j s_j cos(j pi x/L)`. Retaining a
//! declared set of longitudinal coordinates gives the positive potential
//! `EA L s_0²/2 + EA L/4 sum_j (d_j u_j + s_j)²`. Other zero-mean strain
//! channels are statically equilibrated. They have NO retained inertia or
//! resolved propagation; a small selection is not a full axial-wave claim.
//! With no retained coordinates this is the Kirchhoff--Carrier potential.
//!
//! This module owns energy and its work-conjugate discrete gradient, not a
//! propagator, damping law, force lag, or instrument-specific material.
//! Derivation uses cosine product identities, with no spatial quadrature in
//! the runtime. Boundary motion and modal-coordinate pullback are caller-owned.

use fs_math::det;

use crate::{KcStringParams, NlModalError};

/// Maximum number of transverse sine modes per direction.
pub const MAX_TRANSVERSE_MODES: usize = 512;
/// Maximum number of explicitly retained longitudinal sine coordinates.
pub const MAX_LONGITUDINAL_MODES: usize = 32;
/// Largest axial index, covering the quadratic strain band of 512 sine modes.
pub const MAX_LONGITUDINAL_INDEX: usize = 2 * MAX_TRANSVERSE_MODES;

fn invalid(what: &'static str) -> NlModalError {
    NlModalError::Parameter { what }
}

/// Dot product with four independent partial sums (equal-length slices).
fn dot4(a: &[f64], b: &[f64]) -> f64 {
    let mut sums = [0.0_f64; 4];
    let (ca, cb) = (a.chunks_exact(4), b.chunks_exact(4));
    let tail: f64 = ca.remainder().iter().zip(cb.remainder()).map(|(x, y)| x * y).sum();
    for (x, y) in ca.zip(cb) {
        for i in 0..4 {
            sums[i] += x[i] * y[i];
        }
    }
    (sums[0] + sums[1]) + (sums[2] + sums[3]) + tail
}

/// Conservative axial storage for a prestressed string in modal coordinates.
///
/// Slopes are plane-major arrays `[a_0, a_1, ..., a_N]`, repeated for a
/// second direction if selected. `a_0` is the physical chord slope `b/L`;
/// `a_n = sqrt(2/(mu L)) (n pi/L) q_relative_n`. Axial displacements use
/// mass-normalized sine coordinates, in the supplied mode-index order.
/// Kinetic energy is `sum_j velocity_j²/2`; it is not included by `energy`.
#[derive(Clone, Debug)]
pub struct LongitudinalString {
    transverse_modes: usize,
    planes: usize,
    axial_modes: Vec<usize>,
    axial_omegas: Vec<f64>,
    axial_slope_scales: Vec<f64>,
    axial_displacement_scale: f64,
    mean_energy_coefficient: f64,
}

impl LongitudinalString {
    /// Prepare exact Fourier strain coefficients for explicit axial modes.
    ///
    /// Require finite positive `length`, `lin_density`, `tension` and `ea`,
    /// zero to 512 transverse modes, one or two directions, and at most 32
    /// strictly increasing axial indices in 1..=1024. Empty axial selection is
    /// the statically condensed Kirchhoff--Carrier image. A zero transverse
    /// count still retains the moving chord and any selected axial modes.
    /// Frequency/cutoff admission belongs to the caller's declared clock;
    /// construction never removes a requested high-frequency coordinate.
    ///
    /// # Errors
    /// Refuses invalid dimensions, material, mode order, or unrepresentable
    /// mass normalization, stiffness, or modal angular frequency.
    pub fn new(
        params: KcStringParams,
        transverse_modes: usize,
        planes: usize,
        axial_modes: &[usize],
    ) -> Result<Self, NlModalError> {
        if [params.length, params.lin_density, params.tension, params.ea]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || transverse_modes > MAX_TRANSVERSE_MODES
            || !(1..=2).contains(&planes)
            || axial_modes.len() > MAX_LONGITUDINAL_MODES
            || axial_modes.first() == Some(&0)
            || axial_modes
                .last()
                .is_some_and(|index| *index > MAX_LONGITUDINAL_INDEX)
            || axial_modes.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "longitudinal string material, dimensions or mode order",
            ));
        }
        let normalization = det::sqrt(2.0 / (params.lin_density * params.length));
        let wave_speed = det::sqrt(params.ea / params.lin_density);
        let mean_energy_coefficient = 0.5 * params.ea * params.length;
        if [normalization, wave_speed, mean_energy_coefficient]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid("longitudinal string mass or stiffness scale"));
        }
        let mut axial_omegas = Vec::with_capacity(axial_modes.len());
        let mut axial_slope_scales = Vec::with_capacity(axial_modes.len());
        for &index in axial_modes {
            let wave_number = index as f64 * core::f64::consts::PI / params.length;
            let scale = wave_number * normalization;
            let omega = wave_number * wave_speed;
            if !scale.is_finite()
                || scale <= 0.0
                || !omega.is_finite()
                || omega <= 0.0
                || !(omega * omega).is_finite()
            {
                return Err(invalid(
                    "longitudinal string mode scale or stiffness overflow",
                ));
            }
            axial_omegas.push(omega);
            axial_slope_scales.push(scale);
        }
        Ok(Self {
            transverse_modes,
            planes,
            axial_modes: axial_modes.to_vec(),
            axial_omegas,
            axial_slope_scales,
            axial_displacement_scale: normalization,
            mean_energy_coefficient,
        })
    }

    /// Number of slope coefficients, including each direction's chord.
    #[must_use]
    pub fn slope_count(&self) -> usize {
        self.planes * (self.transverse_modes + 1)
    }

    /// Explicit one-based indices of retained longitudinal sine modes.
    #[must_use]
    pub fn axial_modes(&self) -> &[usize] {
        &self.axial_modes
    }

    /// Small-amplitude angular frequencies [rad/s], with full axial inertia.
    #[must_use]
    pub fn axial_omegas(&self) -> &[f64] {
        &self.axial_omegas
    }

    /// Derivatives of mass-normalized axial sine modes [1/(m sqrt(kg))].
    #[must_use]
    pub fn axial_slope_scales(&self) -> &[f64] {
        &self.axial_slope_scales
    }

    /// Work-conjugate axial displacement/force row [1/sqrt(kg)] at a physical
    /// fraction of the complete fixed-end span. Uses the SAME solid-mass
    /// normalization and retained mode order as strain and kinetic energy.
    /// Endpoints return exact zeros; this does not supply moving-end inertia.
    ///
    /// # Errors
    /// Refuses nonfinite fractions or positions outside the closed span.
    pub fn axial_displacement_shapes(&self, fraction: f64) -> Result<Vec<f64>, NlModalError> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("axial displacement needs a physical station on the complete span"));
        }
        if fraction == 0.0 || fraction == 1.0 {
            return Ok(vec![0.0; self.axial_modes.len()]);
        }
        Ok(self.axial_modes.iter().map(|&mode| {
            self.axial_displacement_scale
                * det::sin(mode as f64 * core::f64::consts::PI * fraction)
        }).collect())
    }

    fn validate_state(&self, slopes: &[f64], axial: &[f64]) -> Result<(), NlModalError> {
        if slopes.len() != self.slope_count()
            || axial.len() != self.axial_modes.len()
            || slopes.iter().chain(axial).any(|x| !x.is_finite())
        {
            return Err(invalid("complete finite longitudinal string coordinates"));
        }
        Ok(())
    }

    fn geometric_mean(&self, slopes: &[f64]) -> f64 {
        slopes
            .chunks_exact(self.transverse_modes + 1)
            .map(|plane| {
                0.5 * plane[0] * plane[0] + 0.25 * plane[1..].iter().map(|x| x * x).sum::<f64>()
            })
            .sum()
    }

    fn geometric_harmonic(&self, slopes: &[f64], mode: usize) -> f64 {
        let n = self.transverse_modes;
        let mut coefficient = 0.0;
        for plane in slopes.chunks_exact(n + 1) {
            if mode <= n {
                coefficient += plane[0] * plane[mode];
                for j in 1..=n - mode {
                    coefficient += 0.5 * plane[j] * plane[j + mode];
                }
            }
            // The complete ordered sum includes the diagonal once. Bounds
            // depend on transverse count, never on an unbounded axial index.
            let low = mode.saturating_sub(n).max(1);
            let high = (mode - 1).min(n);
            for j in low..=high {
                coefficient += 0.25 * plane[j] * plane[mode - j];
            }
        }
        coefficient
    }

    /// Total axial elastic potential [J], including the KC mean and retained
    /// longitudinal quadratic stiffness; transverse linear storage is separate.
    /// This method allocates no memory.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite coordinates or nonfinite derived energy.
    pub fn energy(&self, slopes: &[f64], axial: &[f64]) -> Result<f64, NlModalError> {
        self.validate_state(slopes, axial)?;
        let mean = self.geometric_mean(slopes);
        let mut energy = self.mean_energy_coefficient * mean * mean;
        for (j, &mode) in self.axial_modes.iter().enumerate() {
            let strain =
                self.axial_slope_scales[j] * axial[j] + self.geometric_harmonic(slopes, mode);
            energy += 0.5 * self.mean_energy_coefficient * strain * strain;
        }
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(invalid("longitudinal string energy overflow"))
        }
    }

    /// Mean additional material strain and the sum of absolute nonuniform
    /// strain coefficients. The complete reduced strain lies in
    /// `[mean - nonuniform, mean + nonuniform]` everywhere on the span.
    ///
    /// For installed tension `T`, a conservative local tension interval is
    /// `T + EA * [mean - nonuniform, mean + nonuniform]`. Omitted channels are
    /// statically equilibrated, so they add neither strain nor fictitious
    /// inertial tension. The bound need not be attained at one position.
    /// This method allocates no memory.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite coordinates or nonfinite derived strain.
    pub fn strain_bounds(&self, slopes: &[f64], axial: &[f64]) -> Result<(f64, f64), NlModalError> {
        self.validate_state(slopes, axial)?;
        let mean = self.geometric_mean(slopes);
        let mut nonuniform = 0.0;
        for (j, &mode) in self.axial_modes.iter().enumerate() {
            nonuniform += (self.axial_slope_scales[j] * axial[j]
                + self.geometric_harmonic(slopes, mode))
            .abs();
        }
        if mean.is_finite() && nonuniform.is_finite() {
            Ok((mean, nonuniform))
        } else {
            Err(invalid("longitudinal string strain overflow"))
        }
    }

    /// Additional reduced material strain at a span fraction.
    /// `fraction` is `x/L` in `[0,1]`; fixed axial displacement endpoints can
    /// still carry nonzero material strain and work-conjugate tension.
    ///
    /// # Errors
    /// Refuses invalid coordinates, a fraction outside the span, or overflow.
    pub fn strain_at(
        &self,
        slopes: &[f64],
        axial: &[f64],
        fraction: f64,
    ) -> Result<f64, NlModalError> {
        self.validate_state(slopes, axial)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("longitudinal string span fraction"));
        }
        let mut strain = self.geometric_mean(slopes);
        for (j, &mode) in self.axial_modes.iter().enumerate() {
            strain += (self.axial_slope_scales[j] * axial[j]
                + self.geometric_harmonic(slopes, mode))
                * det::cos(mode as f64 * core::f64::consts::PI * fraction);
        }
        if strain.is_finite() {
            Ok(strain)
        } else {
            Err(invalid("longitudinal string strain overflow"))
        }
    }

    /// `EA L / 2` [N m]: storage is this times `s_0^2` plus half of it times
    /// each retained `(d_j u_j + s_j)^2`.
    #[must_use]
    pub fn strain_energy_coefficient(&self) -> f64 {
        self.mean_energy_coefficient
    }

    /// Geometric strain of one state: returns the mean `s_0` and writes each
    /// retained axial mode's cosine coefficient `s_j`, one pass per harmonic.
    /// With [`Self::storage_from_strain`] and [`Self::geometric_strain_gradient`]
    /// it gives `energy`, `strain_bounds` and the equal-endpoint `discrete_gradient`
    /// of one state, to rounding (the products are summed in a different order, for
    /// a real-time caller), without recomputing a harmonic. Allocates nothing.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite slopes or a wrong output length.
    pub fn geometric_strain(&self, slopes: &[f64], harmonics: &mut [f64]) -> Result<f64, NlModalError> {
        if slopes.len() != self.slope_count()
            || harmonics.len() != self.axial_modes.len()
            || slopes.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("complete finite longitudinal string slopes"));
        }
        let n = self.transverse_modes;
        for (slot, &mode) in harmonics.iter_mut().zip(&self.axial_modes) {
            let mut total = 0.0;
            for plane in slopes.chunks_exact(n + 1) {
                if mode <= n {
                    total += plane[0] * plane[mode] + 0.5 * dot4(&plane[1..=n - mode], &plane[1 + mode..]);
                }
                let low = mode.saturating_sub(n).max(1);
                let high = (mode - 1).min(n);
                for j in low..=high {
                    total += 0.25 * plane[j] * plane[mode - j];
                }
            }
            *slot = total;
        }
        Ok(self.geometric_mean(slopes))
    }

    /// Material strain rate of the SAME retained moderate-slope strain.
    ///
    /// Returns the mean rate and writes the retained cosine coefficients of
    /// `u_xt + y_x*y_xt + z_x*z_xt`. `axial_velocity` is in the original solid
    /// mass basis; slope velocities include any explicitly supplied chord.
    /// Because each geometric coefficient is quadratic, evaluating this at
    /// midpoint slopes with coordinate difference / dt gives the EXACT finite
    /// strain increment / dt, including both transverse polarizations.
    /// This is kinematics only: no viscosity, loss or time step is inferred.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite inputs, wrong output extent or overflow.
    /// Output contents are unspecified after an error. Allocates nothing.
    pub fn strain_rate(
        &self,
        slopes: &[f64],
        slope_velocity: &[f64],
        axial_velocity: &[f64],
        harmonics: &mut [f64],
    ) -> Result<f64, NlModalError> {
        self.validate_state(slopes, axial_velocity)?;
        if slope_velocity.len() != self.slope_count()
            || harmonics.len() != self.axial_modes.len()
            || slope_velocity.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("complete finite longitudinal string strain rates"));
        }
        let n = self.transverse_modes;
        let mut mean = 0.0;
        for (a, v) in slopes.chunks_exact(n + 1).zip(slope_velocity.chunks_exact(n + 1)) {
            mean += a[0] * v[0] + 0.5 * dot4(&a[1..], &v[1..]);
        }
        for (j, (&mode, out)) in self.axial_modes.iter().zip(harmonics.iter_mut()).enumerate() {
            let mut rate = self.axial_slope_scales[j] * axial_velocity[j];
            for (a, v) in slopes.chunks_exact(n + 1).zip(slope_velocity.chunks_exact(n + 1)) {
                if mode <= n {
                    rate += a[0] * v[mode] + v[0] * a[mode];
                    for k in 1..=n - mode {
                        rate += 0.5 * (a[k] * v[k + mode] + v[k] * a[k + mode]);
                    }
                }
                let low = mode.saturating_sub(n).max(1);
                let high = (mode - 1).min(n);
                for k in low..=high {
                    rate += 0.25 * (a[k] * v[mode - k] + v[k] * a[mode - k]);
                }
            }
            *out = rate;
        }
        if !mean.is_finite() || harmonics.iter().any(|v| !v.is_finite()) {
            return Err(invalid("longitudinal string strain-rate overflow"));
        }
        Ok(mean)
    }

    /// Storage and nonuniform strain bound from [`Self::geometric_strain`]'s
    /// output and the axial coordinates; writes each retained mode's total
    /// strain coefficient `d_j u_j + s_j` to `strain`. Same arithmetic as
    /// `energy` and `strain_bounds`.
    ///
    /// # Errors
    /// Refuses wrong lengths, nonfinite input, or overflow.
    pub fn storage_from_strain(
        &self,
        mean: f64,
        harmonics: &[f64],
        axial: &[f64],
        strain: &mut [f64],
    ) -> Result<(f64, f64), NlModalError> {
        let n = self.axial_modes.len();
        if harmonics.len() != n || axial.len() != n || strain.len() != n
            || !mean.is_finite() || harmonics.iter().chain(axial).any(|x| !x.is_finite())
        {
            return Err(invalid("complete finite longitudinal string strain"));
        }
        let mut energy = self.mean_energy_coefficient * mean * mean;
        let mut nonuniform = 0.0;
        for (((slot, &s), &u), &scale) in strain.iter_mut().zip(harmonics).zip(axial).zip(&self.axial_slope_scales) {
            *slot = scale * u + s;
            energy += 0.5 * self.mean_energy_coefficient * *slot * *slot;
            nonuniform += slot.abs();
        }
        if energy.is_finite() && nonuniform.is_finite() {
            Ok((energy, nonuniform))
        } else {
            Err(invalid("longitudinal string energy overflow"))
        }
    }

    /// Gradient with respect to the slope coefficients of
    /// `mean_weight * s_0 + sum_j weights[j] * s_j` at one state. With
    /// `mean_weight = 2 c s_0` and `weights[j] = c (d_j u_j + s_j)`, `c` the
    /// [`Self::strain_energy_coefficient`], this is the slope part of the
    /// equal-endpoint `discrete_gradient`, to rounding. Allocates nothing.
    ///
    /// # Errors
    /// Refuses wrong lengths, nonfinite input, or overflow.
    pub fn geometric_strain_gradient(
        &self,
        slopes: &[f64],
        mean_weight: f64,
        weights: &[f64],
        out: &mut [f64],
    ) -> Result<(), NlModalError> {
        if slopes.len() != self.slope_count()
            || out.len() != self.slope_count()
            || weights.len() != self.axial_modes.len()
            || !mean_weight.is_finite()
            || slopes.iter().chain(weights).any(|x| !x.is_finite())
        {
            return Err(invalid("complete finite longitudinal string gradient input"));
        }
        let n = self.transverse_modes;
        for (a, out) in slopes.chunks_exact(n + 1).zip(out.chunks_exact_mut(n + 1)) {
            out[0] = mean_weight * a[0];
            for j in 1..=n {
                out[j] = 0.5 * mean_weight * a[j];
            }
        }
        for (&factor, &mode) in weights.iter().zip(&self.axial_modes) {
            for (a, out) in slopes.chunks_exact(n + 1).zip(out.chunks_exact_mut(n + 1)) {
                if mode <= n {
                    out[0] += factor * a[mode];
                    out[mode] += factor * a[0];
                    // Two contiguous axpys (vectorizable) instead of one interleaved loop.
                    let half = 0.5 * factor;
                    for (o, &x) in out[1..=n - mode].iter_mut().zip(&a[1 + mode..]) {
                        *o += half * x;
                    }
                    for (o, &x) in out[1 + mode..].iter_mut().zip(&a[1..=n - mode]) {
                        *o += half * x;
                    }
                }
                let low = mode.saturating_sub(n).max(1);
                let high = (mode - 1).min(n);
                for k in low..=high {
                    let other = mode - k;
                    out[k] += 0.25 * factor * a[other];
                    out[other] += 0.25 * factor * a[k];
                }
            }
        }
        if out.iter().any(|x| !x.is_finite()) {
            return Err(invalid("longitudinal string gradient overflow"));
        }
        Ok(())
    }

    /// Exact discrete gradient of total axial storage between two states.
    ///
    /// Dotting the returned slope/axial gradients with their coordinate changes
    /// equals `energy(b) - energy(a)`, to arithmetic roundoff. Equal endpoints
    /// give the ordinary force gradient. Outputs are caller-owned scratch,
    /// overwritten on success; after any error their contents are unspecified.
    /// No allocation, iteration, force lag, or division by state differences.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite states, wrong output dimensions or overflow.
    pub fn discrete_gradient(
        &self,
        slopes_a: &[f64],
        slopes_b: &[f64],
        axial_a: &[f64],
        axial_b: &[f64],
        slope_gradient: &mut [f64],
        axial_gradient: &mut [f64],
    ) -> Result<(), NlModalError> {
        self.gradient_impl(
            slopes_a,
            slopes_b,
            axial_a,
            axial_b,
            slope_gradient,
            axial_gradient,
            false,
        )
    }

    /// Discrete gradient of the KC mean plus nonlinear axial coupling, with
    /// ONLY the retained axial linear stiffness removed from the force.
    ///
    /// Use this when the caller's existing linear propagator already advances
    /// `omega_j² u_j`. The total storage remains `energy`: removing its positive
    /// quadratic part from energy as well would omit physical axial stiffness.
    /// This residual can have either sign; total potential is nonnegative.
    /// Other input/output and allocation rules match `discrete_gradient`.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite states, wrong output dimensions or overflow.
    pub fn coupling_discrete_gradient(
        &self,
        slopes_a: &[f64],
        slopes_b: &[f64],
        axial_a: &[f64],
        axial_b: &[f64],
        slope_gradient: &mut [f64],
        axial_gradient: &mut [f64],
    ) -> Result<(), NlModalError> {
        self.gradient_impl(
            slopes_a,
            slopes_b,
            axial_a,
            axial_b,
            slope_gradient,
            axial_gradient,
            true,
        )
    }

    fn gradient_impl(
        &self,
        slopes_a: &[f64],
        slopes_b: &[f64],
        axial_a: &[f64],
        axial_b: &[f64],
        slope_gradient: &mut [f64],
        axial_gradient: &mut [f64],
        remove_axial_linear: bool,
    ) -> Result<(), NlModalError> {
        self.validate_state(slopes_a, axial_a)?;
        self.validate_state(slopes_b, axial_b)?;
        if slope_gradient.len() != self.slope_count()
            || axial_gradient.len() != self.axial_modes.len()
        {
            return Err(invalid("longitudinal string gradient output dimensions"));
        }
        slope_gradient.fill(0.0);
        axial_gradient.fill(0.0);
        let n = self.transverse_modes;
        let mean_factor = 2.0
            * self.mean_energy_coefficient
            * f64::midpoint(self.geometric_mean(slopes_a), self.geometric_mean(slopes_b));
        for ((a, b), out) in slopes_a
            .chunks_exact(n + 1)
            .zip(slopes_b.chunks_exact(n + 1))
            .zip(slope_gradient.chunks_exact_mut(n + 1))
        {
            out[0] = mean_factor * f64::midpoint(a[0], b[0]);
            for j in 1..=n {
                out[j] = 0.5 * mean_factor * f64::midpoint(a[j], b[j]);
            }
        }
        for (j, &mode) in self.axial_modes.iter().enumerate() {
            let geometric_a = self.geometric_harmonic(slopes_a, mode);
            let geometric_b = self.geometric_harmonic(slopes_b, mode);
            let geometric_mid = f64::midpoint(geometric_a, geometric_b);
            let strain_mid =
                self.axial_slope_scales[j] * f64::midpoint(axial_a[j], axial_b[j]) + geometric_mid;
            let factor = self.mean_energy_coefficient * strain_mid;
            axial_gradient[j] = self.mean_energy_coefficient
                * self.axial_slope_scales[j]
                * if remove_axial_linear {
                    geometric_mid
                } else {
                    strain_mid
                };
            // Each geometric strain coefficient is a quadratic polynomial.
            // Its exact two-point gradient is its derivative at the midpoint.
            for ((a, b), out) in slopes_a
                .chunks_exact(n + 1)
                .zip(slopes_b.chunks_exact(n + 1))
                .zip(slope_gradient.chunks_exact_mut(n + 1))
            {
                if mode <= n {
                    out[0] += factor * f64::midpoint(a[mode], b[mode]);
                    out[mode] += factor * f64::midpoint(a[0], b[0]);
                    for k in 1..=n - mode {
                        out[k] += 0.5 * factor * f64::midpoint(a[k + mode], b[k + mode]);
                        out[k + mode] += 0.5 * factor * f64::midpoint(a[k], b[k]);
                    }
                }
                let low = mode.saturating_sub(n).max(1);
                let high = (mode - 1).min(n);
                for k in low..=high {
                    let other = mode - k;
                    out[k] += 0.25 * factor * f64::midpoint(a[other], b[other]);
                    out[other] += 0.25 * factor * f64::midpoint(a[k], b[k]);
                }
            }
        }
        if slope_gradient
            .iter()
            .chain(axial_gradient.iter())
            .any(|x| !x.is_finite())
        {
            return Err(invalid("longitudinal string discrete gradient overflow"));
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "longitudinal_tests.rs"]
mod tests;
