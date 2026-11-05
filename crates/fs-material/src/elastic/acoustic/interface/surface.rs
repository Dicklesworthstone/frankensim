//! Coherent reflection by an infinite free or clamped isotropic solid surface.
//!
//! The original solid occupies z < 0; +z points out of it. The wave, traction,
//! phase and modal velocity conventions are those of the bonded interface.
//! A free surface sets all three tractions to zero; a stationary clamp sets
//! all three velocities to zero. P and SV generally convert into each other,
//! including an evanescent P field for supercritical SV incidence. SH remains
//! independent. Neither boundary absorbs energy or introduces a second medium.
//!
//! The same original tensor supplies every wave trace. One 2x2 boundary solve
//! is prepared per tangential slowness, then coherent drives allocate nothing.
//! This is a bulk plane-wave boundary, not a finite plate, flexural edge,
//! Rayleigh-wave source, compliant support or surface-roughness model.
//! Free-surface analytical coefficients: Sánchez-Sesma et al. (2023), Appendix 2,
//! <https://doi.org/10.1186/s40623-023-01780-0>.

use super::*;

/// Explicit stationary, lossless mechanical boundary condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElasticSurfaceCondition {
    /// Vacuum/free boundary: `sigma.e_z = 0`, with unconstrained velocity.
    TractionFree,
    /// Ideal immovable support: `v = 0`, with unconstrained reaction traction.
    Clamped,
}

/// Original isotropic material and its cached P/shear speeds and impedances.
#[derive(Debug)]
pub struct IsotropicElasticSurface<'a> {
    state: &'a ElasticAcousticState,
    condition: ElasticSurfaceCondition,
    speeds: [f64; 2],
    impedances: [f64; 2],
}

/// Fixed-size reflection operator at a conserved tangential slowness.
#[derive(Debug)]
pub struct PreparedElasticSurface<'a> {
    state: &'a ElasticAcousticState,
    condition: ElasticSurfaceCondition,
    tangential_slowness_s_m: f64,
    incoming: [ElasticWaveTrace; 3],
    reflected: [ElasticWaveTrace; 3],
    scattering: [[C64; 3]; 3],
}

/// Coherent physical fields and independently reconstructed surface work.
/// Modal arrays are ordered P, SV, SH; powers are per physical surface area.
#[derive(Clone, Debug)]
pub struct ElasticSurfaceResponse {
    pub material_state_identity: ContentHash,
    pub incoming_velocity_m_s: [C64; 3],
    pub reflected_velocity_m_s: [C64; 3],
    pub incoming_normal_power_w_m2: [f64; 3],
    pub reflected_normal_power_w_m2: [f64; 3],
    pub total_incoming_power_w_m2: f64,
    pub total_reflected_power_w_m2: f64,
    pub boundary_velocity_m_s: [C64; 3],
    /// Traction `sigma.e_z` on the outward normal, in Pa.
    pub boundary_traction_pa: [C64; 3],
    /// Independently evaluated `-Re(v*.traction)/2`; zero for either condition.
    pub boundary_normal_flux_w_m2: f64,
    /// Conservative Mandel-strain norm bound throughout the solid half-space.
    pub strain_norm_bound: f64,
}

impl<'a> IsotropicElasticSurface<'a> {
    /// Retain the original stiffness, density, identity and strain domain.
    /// An anisotropic source is never projected to an isotropic approximation.
    pub fn new(
        state: &'a ElasticAcousticState,
        condition: ElasticSurfaceCondition,
    ) -> Result<Self, MaterialError> {
        if !state.is_explicitly_isotropic() {
            return Err(invalid(
                "elastic surface requires an explicitly isotropic state",
            ));
        }
        let modes = state.directional_waves([0.0, 0.0, 1.0])?.modes;
        Ok(Self {
            state,
            condition,
            speeds: [modes[2].phase_speed_m_s, modes[0].phase_speed_m_s],
            impedances: [modes[2].impedance_pa_s_m, modes[0].impedance_pa_s_m],
        })
    }

    pub const fn state(&self) -> &'a ElasticAcousticState {
        self.state
    }

    /// P and shear phase speeds [m/s], in that order.
    pub const fn phase_speeds_m_s(&self) -> [f64; 2] {
        self.speeds
    }

    /// Fix tangential slowness using a propagating incident mode. Every other
    /// propagating channel at that slowness can subsequently be driven too.
    pub fn prepare_incidence(
        &self,
        mode: ElasticMode,
        cosine: f64,
    ) -> Result<PreparedElasticSurface<'a>, MaterialError> {
        self.prepare_incidence_cancellable(mode, cosine, &mut || false)
    }

    /// Bounded preparation with cancellation before traces and each solve.
    /// Exact grazing incidence and singular boundary systems refuse.
    pub fn prepare_incidence_cancellable(
        &self,
        mode: ElasticMode,
        cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticSurface<'a>, MaterialError> {
        poll(cancelled)?;
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(invalid("elastic surface incident cosine must be in (0,1]"));
        }
        let incident_speed = self.speeds[mode.index().min(1)];
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let tangential = sine / incident_speed;
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) {
            return Err(invalid(
                "elastic surface tangential slowness is unrepresentable",
            ));
        }
        let mut result = PreparedElasticSurface {
            state: self.state,
            condition: self.condition,
            tangential_slowness_s_m: tangential,
            incoming: [empty_trace(); 3],
            reflected: [empty_trace(); 3],
            scattering: [[C64::ZERO; 3]; 3],
        };
        for kind in 0..3 {
            poll(cancelled)?;
            let speed = self.speeds[kind.min(1)];
            let normal = normal_ratio(incident_speed, speed, cosine)?;
            let ratio = speed / incident_speed;
            let wave_cosine = normal.scale(ratio);
            let wave_sine = sine * ratio;
            let normal_slowness = divided(normal, incident_speed);
            if !finite(wave_cosine)
                || !wave_sine.is_finite()
                || (normal != C64::ZERO && normal_slowness == C64::ZERO)
            {
                return Err(invalid("elastic surface wave direction is unrepresentable"));
            }
            for (direction, traces) in [(1.0, &mut result.incoming), (-1.0, &mut result.reflected)]
            {
                traces[kind] = trace(
                    self.state,
                    self.impedances,
                    kind,
                    direction,
                    wave_sine,
                    wave_cosine,
                    tangential,
                    normal_slowness,
                )?;
            }
        }
        result.prepare_psv(cancelled)?;
        if result.incoming[2].normal_power_factor_pa_s_m > 0.0 {
            result.scattering[2][2] = match self.condition {
                ElasticSurfaceCondition::TractionFree => C64::ONE,
                ElasticSurfaceCondition::Clamped => -C64::ONE,
            };
        }
        poll(cancelled)?;
        Ok(result)
    }
}

impl PreparedElasticSurface<'_> {
    pub const fn condition(&self) -> ElasticSurfaceCondition {
        self.condition
    }

    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential_slowness_s_m
    }

    /// Evanescent reflected modes remain physical fields with zero normal power.
    pub fn reflected_mode(&self, mode: ElasticMode) -> ElasticWaveTrace {
        self.reflected[mode.index()]
    }

    /// Only propagating incoming modes are available for external excitation.
    pub fn incoming_mode(&self, mode: ElasticMode) -> Option<ElasticWaveTrace> {
        let wave = self.incoming[mode.index()];
        (wave.normal_power_factor_pa_s_m > 0.0).then_some(wave)
    }

    /// Reflected P/SV/SH modal velocities per unit incident modal velocity.
    pub fn velocity_coefficients(&self, mode: ElasticMode) -> Option<[C64; 3]> {
        self.incoming_mode(mode)?;
        Some(core::array::from_fn(|i| self.scattering[i][mode.index()]))
    }

    /// Apply coherent peak modal velocities [m/s] at a common frequency/phase
    /// reference. No allocation, factorization, clipping or phase averaging.
    pub fn drive(&self, incoming: [C64; 3]) -> Result<ElasticSurfaceResponse, MaterialError> {
        for (i, &amplitude) in incoming.iter().enumerate() {
            if !finite(amplitude)
                || (amplitude != C64::ZERO && self.incoming[i].normal_power_factor_pa_s_m == 0.0)
            {
                return Err(invalid(
                    "elastic surface drive needs finite propagating amplitudes",
                ));
            }
        }
        let reflected: [C64; 3] = core::array::from_fn(|i| {
            complex_sum((0..3).map(|j| self.scattering[i][j] * incoming[j]))
        });
        if reflected.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "elastic surface reflected amplitude is unrepresentable",
            ));
        }
        let mut input_power = [0.0; 3];
        let mut reflected_power = [0.0; 3];
        let mut strain = Dd::ZERO;
        for mode in 0..3 {
            input_power[mode] = power(self.incoming[mode], incoming[mode])?;
            reflected_power[mode] = power(self.reflected[mode], reflected[mode])?;
            for (wave, amplitude) in [
                (self.incoming[mode], incoming[mode]),
                (self.reflected[mode], reflected[mode]),
            ] {
                let norm = wave
                    .strain_mandel_s_m
                    .iter()
                    .fold(0.0, |n, e| det::hypot(n, e.abs()));
                strain = strain + Dd::from_f64(norm * amplitude.abs());
            }
        }
        let strain_bound = strain.to_f64();
        if !strain_bound.is_finite() || strain_bound > self.state.maximum_absolute_mandel_strain() {
            return Err(invalid(
                "elastic surface drive exceeds original strain domain",
            ));
        }
        let velocity: [C64; 3] = core::array::from_fn(|component| {
            complex_sum((0..3).flat_map(|m| {
                [
                    self.incoming[m].velocity[component] * incoming[m],
                    self.reflected[m].velocity[component] * reflected[m],
                ]
            }))
        });
        let traction: [C64; 3] = core::array::from_fn(|component| {
            complex_sum((0..3).flat_map(|m| {
                [
                    self.incoming[m].traction_pa_s_m[component] * incoming[m],
                    self.reflected[m].traction_pa_s_m[component] * reflected[m],
                ]
            }))
        });
        let total_incoming = real_sum(input_power);
        let total_reflected = real_sum(reflected_power);
        let flux = -0.5 * real_sum((0..3).map(|i| (velocity[i].conj() * traction[i]).re));
        if velocity.iter().chain(&traction).any(|&x| !finite(x))
            || [total_incoming, total_reflected, flux]
                .into_iter()
                .any(|x| !x.is_finite())
            || !balanced(total_incoming, total_reflected, total_incoming)
            || !balanced(flux, 0.0, total_incoming)
        {
            return Err(invalid(
                "elastic surface power is unrepresentable or does not close",
            ));
        }
        Ok(ElasticSurfaceResponse {
            material_state_identity: self.state.material_state_identity(),
            incoming_velocity_m_s: incoming,
            reflected_velocity_m_s: reflected,
            incoming_normal_power_w_m2: input_power,
            reflected_normal_power_w_m2: reflected_power,
            total_incoming_power_w_m2: total_incoming,
            total_reflected_power_w_m2: total_reflected,
            boundary_velocity_m_s: velocity,
            boundary_traction_pa: traction,
            boundary_normal_flux_w_m2: flux,
            strain_norm_bound: strain_bound,
        })
    }

    fn boundary_entry(&self, wave: ElasticWaveTrace, row: usize) -> C64 {
        let component = 2 * row;
        match self.condition {
            ElasticSurfaceCondition::TractionFree => wave.traction_pa_s_m[component],
            ElasticSurfaceCondition::Clamped => wave.velocity[component],
        }
    }

    fn prepare_psv(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialError> {
        let mut matrix = [C64::ZERO; 4];
        let mut scales = [0.0_f64; 2];
        for row in 0..2 {
            for column in 0..2 {
                let entry = self.boundary_entry(self.reflected[column], row);
                matrix[2 * row + column] = entry;
                scales[row] = scales[row].max(entry.abs());
            }
            if scales[row] == 0.0 || !scales[row].is_finite() {
                return Err(invalid(
                    "elastic surface boundary row is singular or unrepresentable",
                ));
            }
            for column in 0..2 {
                matrix[2 * row + column] = divided(matrix[2 * row + column], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 2, &mut *cancelled).map_err(solve_error)?;
        for incident in 0..2 {
            if self.incoming[incident].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let rhs: [C64; 2] = core::array::from_fn(|row| {
                -divided(
                    self.boundary_entry(self.incoming[incident], row),
                    scales[row],
                )
            });
            let mut amplitudes = rhs;
            factor
                .solve_with_cancel(&mut amplitudes, &mut *cancelled)
                .map_err(solve_error)?;
            if amplitudes.iter().any(|&a| !finite(a)) {
                return Err(invalid(
                    "elastic surface reflection solve is unrepresentable",
                ));
            }
            for row in 0..2 {
                let terms: [C64; 2] = core::array::from_fn(|j| matrix[2 * row + j] * amplitudes[j]);
                let residual = (complex_sum(terms) - rhs[row]).abs();
                let magnitude = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                if !residual.is_finite() || (magnitude > 0.0 && residual / magnitude > 1e-11) {
                    return Err(invalid("elastic surface boundary solve did not converge"));
                }
                self.scattering[row][incident] = amplitudes[row];
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
