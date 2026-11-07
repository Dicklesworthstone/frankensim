//! Coherent fields of the original bound velocity-port solutions.
//! Source coefficients combine before observation and outgoing-flux integration.
//! No new boundary solve, pressure-density substitution or modal truncation.

use super::{C64, HelmholtzError, VelocityPortSolution, bad, finite, poll, zeroed};
use super::super::{Prepared, directions};
use fs_math::dd::Dd;

/// Extra resources for one coherent field and its streaming power assessment.
/// Original source columns, prepared receiver rows and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct CoherentPortFieldLimits {
    /// 32*boundary_unknowns + 32*panels + 32*A + 256 numeric bytes,
    /// including the original angular-direction table and streaming wave scratch.
    /// A=64*4^angular_refinement uses the prepared boundary's unchanged rule.
    pub maximum_payload_bytes: usize,
    /// 64*ports*(boundary_unknowns+panels) + 128*panels projection units.
    pub maximum_work: u64,
}

/// Same original P1/P0 boundary field, carrying complete complex interference.
/// All geometry, frequency, medium and densities are immutable and source-bound.
/// Unit-source backward errors are retained explicitly; their maximum is NOT
/// a backward/forward error bound for a nearly cancelling coherent field.
pub struct CoherentVelocityPortField<'p, 's> {
    pub(in crate::reciprocal_radiation) owner: &'p Prepared<'s>,
    pub(in crate::reciprocal_radiation) coefficients: Vec<C64>,
    pub(in crate::reciprocal_radiation) k: f64,
    pub(in crate::reciprocal_radiation) medium: super::Medium,
    pressure: Vec<C64>,
    velocity: Vec<C64>,
    surface_power: C64,
    apparent_power: f64,
    maximum_source_backward_residual: f64,
}

/// Original boundary work and independently integrated coherent outgoing flux.
/// Their discrepancy is numerical evidence for this discrete field, not a
/// continuum, spatial-quadrature or cancellation-conditioned source error bound.
#[derive(Clone, Copy, Debug)]
pub struct CoherentRadiationPower {
    /// Half the integral of pressure times conjugate outward velocity [W].
    pub surface_power_w: C64,
    /// Full angular integral, including all source interference [W].
    pub outgoing_power_w: f64,
    /// Absolute real-work/flux difference divided by apparent boundary work
    /// or outgoing flux, whichever is larger. No unit-bearing absolute floor.
    pub relative_flux_defect: f64,
}

fn sum_products(values: impl Iterator<Item = (C64, C64)>) -> C64 {
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for (a, b) in values {
        let (ar, ai, br, bi) = (
            Dd::from_f64(a.re), Dd::from_f64(a.im),
            Dd::from_f64(b.re), Dd::from_f64(b.im),
        );
        re = re + ar * br - ai * bi;
        im = im + ar * bi + ai * br;
    }
    C64::new(re.to_f64(), im.to_f64())
}

impl<'p, 's> VelocityPortSolution<'p, 's> {
    /// Combine the original complete unit-velocity sources without solving again.
    /// Amplitudes use the same physical port units and exp(-i omega t) convention
    /// as the source batch. Signed, dependent, rotational and zero-flow columns
    /// remain in their original order; sources combine as fields, never powers.
    ///
    /// # Errors
    /// Incomplete/nonfinite amplitudes, extra work/payload exhaustion, a changed
    /// source identity, unrepresentable arithmetic or cancellation.
    pub fn coherent(
        &self,
        amplitudes: &[C64],
        limits: CoherentPortFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentVelocityPortField<'p, 's>, HelmholtzError> {
        poll(&mut cancelled)?;
        let first = self.solutions.first().ok_or_else(|| bad("empty coherent velocity-port batch"))?;
        let (n, p, fields) = (first.owner.unknowns(), first.velocity.len(), self.ports());
        let angular = 64_u128 * (1_u128 << (2 * first.owner.budget.angular_refinement));
        let bytes = 32 * n as u128 + 32 * p as u128 + 32 * angular + 256;
        let work = 64 * fields as u128 * (n as u128 + p as u128) + 128 * p as u128;
        if amplitudes.len() != fields || amplitudes.iter().any(|&v| !finite(v))
            || bytes > limits.maximum_payload_bytes as u128
            || work > limits.maximum_work as u128
        {
            return Err(bad("complete coherent port field exceeds amplitude or resource admission"));
        }
        if self.solutions.iter().any(|source|
            !core::ptr::eq(source.owner, first.owner)
                || source.k.to_bits() != first.k.to_bits() || source.medium != first.medium
                || source.coefficients.len() != n || source.velocity.len() != p)
        {
            return Err(bad("coherent velocity-port sources disagree"));
        }
        let mut coefficients = zeroed(n, C64::ZERO)?;
        let mut velocity = zeroed(p, C64::ZERO)?;
        for (i, value) in coefficients.iter_mut().enumerate() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            *value = sum_products(self.solutions.iter().zip(amplitudes)
                .map(|(source, &amplitude)| (source.coefficients[i], amplitude)));
        }
        for (i, value) in velocity.iter_mut().enumerate() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            *value = sum_products(self.solutions.iter().zip(amplitudes)
                .map(|(source, &amplitude)| (source.velocity[i], amplitude)));
        }
        let mut pressure = zeroed(p, C64::ZERO)?;
        let mut apparent = Dd::ZERO;
        for (i, (value, facet)) in pressure.iter_mut().zip(&first.owner.facets).enumerate() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            *value = (coefficients[facet.ids[0]] + coefficients[facet.ids[1]]
                + coefficients[facet.ids[2]]).scale(1.0 / 3.0);
            apparent = apparent + Dd::from_f64(0.5 * facet.area * value.abs() * velocity[i].abs());
        }
        let surface_power = sum_products(pressure.iter().zip(&velocity).zip(&first.owner.facets)
            .map(|((&p, &v), facet)| (p, v.conj().scale(0.5 * facet.area))));
        let apparent_power = apparent.to_f64();
        if coefficients.iter().chain(&pressure).chain(&velocity).any(|&v| !finite(v))
            || !finite(surface_power) || !apparent_power.is_finite()
        {
            return Err(bad("unrepresentable coherent boundary field or work"));
        }
        poll(&mut cancelled)?;
        Ok(CoherentVelocityPortField {
            owner: first.owner, coefficients, k: first.k, medium: first.medium,
            pressure, velocity, surface_power, apparent_power,
            maximum_source_backward_residual: self.solutions.iter()
                .map(|source| source.backward_residual).fold(0.0, f64::max),
        })
    }
}

impl CoherentVelocityPortField<'_, '_> {
    /// Original physical angular frequency [rad/s].
    pub fn omega_rad_s(&self) -> f64 { self.k * self.medium.sound_speed }
    /// Original homogeneous fluid; no inferred attenuation or thermal loss.
    pub fn medium(&self) -> super::Medium { self.medium }
    /// Area means of the actual P1 pressure trace, for boundary work only [Pa].
    pub fn surface_pressure_pa(&self) -> &[C64] { &self.pressure }
    /// Actual prescribed outward P0 surface velocity [m/s].
    pub fn surface_velocity_m_s(&self) -> &[C64] { &self.velocity }
    /// Complex coherent boundary work, with all self and mutual terms [W].
    pub fn surface_power_w(&self) -> C64 { self.surface_power }
    /// Maximum error of the independent original unit-source solves.
    /// It is not relabelled as a residual of this coherent combination.
    pub fn maximum_source_backward_residual(&self) -> f64 {
        self.maximum_source_backward_residual
    }

    /// Far-field amplitude F in p ~ F exp(ikr)/r, from the same density kernels.
    /// # Errors
    /// Invalid direction, original kernel refusal or nonfinite contraction.
    pub fn far_field(&self, direction: [f64; 3]) -> Result<C64, HelmholtzError> {
        let radius = super::super::norm(direction);
        if direction.iter().any(|v| !v.is_finite()) || !radius.is_finite() || radius <= 0.0 {
            return Err(bad("invalid coherent far-field direction"));
        }
        let wave = self.owner.wave_coefficients(self.k, direction.map(|v| v / radius))?;
        let value = sum_products(wave.into_iter().zip(self.coefficients.iter().copied()))
            .scale(1.0 / (4.0 * core::f64::consts::PI));
        if !finite(value) { return Err(bad("coherent far-field overflow")); }
        Ok(value)
    }

    /// Stream the original antipodal angular rule over this coherent field.
    /// Work includes original plane-wave kernel construction and contractions;
    /// scratch was included in CoherentPortFieldLimits. No source solve occurs.
    /// No individual radiated powers are summed and no discrepancy is clipped.
    ///
    /// # Errors
    /// Exceeded work, native kernel refusal, overflow or cancellation.
    pub fn power(
        &self,
        maximum_work: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentRadiationPower, HelmholtzError> {
        poll(&mut cancelled)?;
        let count = 64_u128 * (1_u128 << (2 * self.owner.budget.angular_refinement));
        let work = 128 * count * (self.owner.unknowns() as u128
            + 16 * self.owner.facets.len() as u128);
        if work > maximum_work as u128 {
            return Err(bad("coherent outgoing flux exceeds its declared work allowance"));
        }
        let mut sum = Dd::ZERO;
        for (direction, weight) in directions(self.owner.budget.angular_refinement) {
            poll(&mut cancelled)?;
            let wave = self.owner.wave_coefficients(self.k, direction)?;
            let forward = sum_products(wave.iter().copied().zip(self.coefficients.iter().copied()));
            let opposite = sum_products(wave.iter().map(|a| a.conj())
                .zip(self.coefficients.iter().copied()));
            sum = sum + Dd::from_f64(weight)
                * (Dd::from_f64(forward.norm_sq()) + Dd::from_f64(opposite.norm_sq()));
        }
        let pi = core::f64::consts::PI;
        let outgoing_power_w = sum.to_f64() / (32.0 * pi * pi * self.medium.density * self.medium.sound_speed);
        if !outgoing_power_w.is_finite() || outgoing_power_w < 0.0 {
            return Err(bad("coherent outgoing flux is unrepresentable"));
        }
        let scale = self.apparent_power.max(outgoing_power_w);
        let relative_flux_defect = if scale == 0.0 {
            if self.surface_power.re != 0.0 {
                return Err(bad("nonzero coherent work has no representable normalization"));
            }
            0.0
        } else {
            (self.surface_power.re / scale - outgoing_power_w / scale).abs()
        };
        if !relative_flux_defect.is_finite() {
            return Err(bad("coherent work defect is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(CoherentRadiationPower {
            surface_power_w: self.surface_power, outgoing_power_w, relative_flux_defect,
        })
    }
}

#[cfg(test)]
mod tests;
