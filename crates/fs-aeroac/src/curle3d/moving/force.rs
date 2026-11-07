//! Full fields of moving compact force sources with explicit past forcing.
//!
//! For a force F on the fluid, define
//! `B(x,t) = integral_{t0}^{tau(x,t)} F(s)/|x-y(s)| ds`.
//! Away from the complete source worldline, linear momentum gives
//! `u_F = grad(div B)/(4 pi rho)`. The retarded upper limit is differentiated
//! analytically; the only nonlocal input is
//! `J = integral (3 n (n.F) - F)/R^3 ds`, at the SAME fixed receiver position.
//! The supplied J must include the chosen initial force history. Zero J at t0
//! means zero prior force contribution, not an arbitrary unknown initial field.
//!
//! The pressure is the original moving-source Green-function derivative in
//! [`super`], including acceleration and induction. Writing M=v/c and
//! kappa=1-n.M, the local velocity boundary term is
//! `n p_F/(rho c) + (F_n ((1+kappa)n-M) - kappa F)/(4 pi rho c R^2 kappa^2)`.
//! Adding J/(4 pi rho) recovers the full force field. At rest this reduces to
//! the original stationary impulse, near-field and radiation terms.
//!
//! This module does not approximate a moving force by a delayed stationary
//! source or by a plane wave. History integration, continuous worldline
//! clearance and source-clock coverage remain explicit caller obligations.

use super::{
    CompactSourceStrength, MovingCompactObserver, MovingCompactSample, PreparedMovingCompactSample,
    SourceKinematics, invalid, nonfinite, norm,
};
use crate::AeroacError;
use fs_math::dd::Dd;

/// Original pressure plus full instantaneous velocity of volume/force sources.
///
/// The velocity is the acoustic perturbation in the medium's rest frame.
/// Intensity is the signed instantaneous pressure-work flux, without a harmonic
/// factor of one half. Coherent sources and paths must sum p and u before p*u.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingCompactField {
    /// Original retarded metadata and coherent monopole/dipole pressure [Pa].
    pub sample: MovingCompactSample,
    /// Cartesian acoustic particle velocity [m/s].
    pub particle_velocity_m_s: [f64; 3],
    /// Signed instantaneous pressure-work flux [W/m²].
    pub instantaneous_intensity_w_m2: [f64; 3],
}

/// Cached retarded boundary terms; the past force history is supplied at apply.
///
/// Preparation adds no cost or payload to the existing pressure-only or moving
/// monopole owners. Applications perform bounded constant work and allocate
/// nothing. The history integral must belong to this exact receiver and
/// emission event; an integral at a previous microphone position is invalid.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingForceFieldSample {
    pressure: PreparedMovingCompactSample,
    force_near_velocity: [[f64; 3]; 3],
    force_pressure_velocity: [f64; 3],
    history_velocity_scale: f64,
}

impl MovingCompactObserver {
    /// Prepare full-field boundary coefficients at one original emission event.
    ///
    /// Uses the same pressure preparation and instantaneous subsonic/clearance
    /// admission as [`Self::prepare_emission`]. The caller must separately admit
    /// the complete trajectory, retarded root and retained force history.
    pub fn prepare_force_emission(
        &self,
        source_time_s: f64,
        kinematics: SourceKinematics,
    ) -> Result<PreparedMovingForceFieldSample, AeroacError> {
        let pressure = self.prepare_emission(source_time_s, kinematics)?;
        self.prepare_force_sample(pressure, kinematics)
    }

    // The trajectory owners already solved and admitted this exact retarded
    // event. Reuse their original pressure preparation without a second root
    // or another pressure-coefficient evaluation.
    pub(super) fn prepare_force_sample(
        &self,
        pressure: PreparedMovingCompactSample,
        kinematics: SourceKinematics,
    ) -> Result<PreparedMovingForceFieldSample, AeroacError> {
        let radius = pressure.distance_m();
        let c = self.medium().sound_speed_m_s();
        let rho = self.medium().density_kg_m3();
        let direction = core::array::from_fn::<_, 3, _>(|i| {
            (self.observer_position_m()[i] - kinematics.position_m[i]) / radius
        });
        let mach = kinematics.velocity_m_s.map(|v| v / c);
        let kappa = 1.0 - pressure.radial_mach();
        let history_velocity_scale = (0.25 / core::f64::consts::PI) / rho;
        // Preserve split powers and divisions; never form rho*c or R^2.
        let force_scale =
            (0.25 / core::f64::consts::PI) / radius / kappa / kappa / radius / c / rho;
        let pressure_scale = 1.0 / c / rho;
        if [history_velocity_scale, force_scale, pressure_scale]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid(
                "moving force-field velocity scales are unrepresentable",
            ));
        }
        let force_pressure_velocity = direction.map(|n| n * pressure_scale);
        let force_near_velocity = core::array::from_fn(|i| {
            let radial = (1.0 + kappa).mul_add(direction[i], -mach[i]);
            core::array::from_fn(|j| {
                let diagonal = if i == j { -kappa } else { 0.0 };
                radial.mul_add(direction[j], diagonal) * force_scale
            })
        });
        if force_near_velocity.iter().flatten().any(|v| !v.is_finite())
            || force_pressure_velocity.iter().any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving force-field boundary coefficients"));
        }
        Ok(PreparedMovingForceFieldSample {
            pressure,
            force_near_velocity,
            force_pressure_velocity,
            history_velocity_scale,
        })
    }

    /// Integrand of the past-force contribution at this fixed receiver.
    ///
    /// Returns `(3 n (n.F) - F)/R^3` [N/m³]. Integrate on the ORIGINAL source
    /// clock, from the declared force-history origin to the retarded emission.
    /// Directions and R use each past source position and this SAME receiver.
    /// There is no Doppler factor in this history integrand: the retarded upper
    /// limit and its derivatives are already in the prepared boundary terms.
    ///
    /// This point query checks finite geometry and clearance, but cannot prove
    /// clearance between quadrature nodes. The history owner must establish
    /// that continuous condition before integrating. Work is fixed-size and
    /// allocation-free, with double-double force contractions.
    pub fn force_history_integrand(
        &self,
        source_position_m: [f64; 3],
        force_on_fluid_n: [f64; 3],
    ) -> Result<[f64; 3], AeroacError> {
        if source_position_m
            .iter()
            .chain(&force_on_fluid_n)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving force history geometry or strength"));
        }
        let displacement =
            core::array::from_fn(|i| self.observer_position_m()[i] - source_position_m[i]);
        let radius = norm(displacement);
        if !radius.is_finite() || radius < self.limits().minimum_distance_m {
            return Err(invalid("moving force history source/receiver separation"));
        }
        let inverse_cube = 1.0 / radius / radius / radius;
        if !inverse_cube.is_finite() || inverse_cube <= 0.0 {
            return Err(invalid("moving force history geometry is unrepresentable"));
        }
        let direction = displacement.map(|x| x / radius);
        let result = core::array::from_fn(|i| {
            let mut value = Dd::ZERO;
            for j in 0..3 {
                let diagonal = if i == j { -1.0 } else { 0.0 };
                let coefficient =
                    (3.0 * direction[i]).mul_add(direction[j], diagonal) * inverse_cube;
                value = value + Dd::from_f64(coefficient) * Dd::from_f64(force_on_fluid_n[j]);
            }
            value.to_f64()
        });
        if result.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving force history integrand"));
        }
        Ok(result)
    }
}

impl PreparedMovingForceFieldSample {
    /// Original pressure preparation, including emission/arrival metadata.
    pub const fn pressure_sample(&self) -> &PreparedMovingCompactSample {
        &self.pressure
    }

    /// Evaluate mixed volume/force fields with the original past-force integral.
    ///
    /// Strengths and their derivatives belong to the prepared source time.
    /// `force_history_integral_n_s_m3` is J [N s/m³] from this module's
    /// definition, before division by 4*pi*rho. The caller supplies force
    /// prehistory explicitly; additional homogeneous acoustic fields must be
    /// composed separately. Neither is inferred from the instantaneous force.
    ///
    /// Nonfinite strengths, history, velocity or flux refuse the complete field.
    /// The original pressure-only owner's representable domain is unchanged.
    pub fn evaluate(
        &self,
        strength: CompactSourceStrength,
        force_history_integral_n_s_m3: [f64; 3],
    ) -> Result<MovingCompactField, AeroacError> {
        let (sample, pressure, velocity) =
            self.evaluate_with_fields_dd(strength, force_history_integral_n_s_m3)?;
        let particle_velocity_m_s = velocity.map(Dd::to_f64);
        let instantaneous_intensity_w_m2 = velocity.map(|u| (pressure * u).to_f64());
        if particle_velocity_m_s
            .iter()
            .chain(&instantaneous_intensity_w_m2)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving compact particle velocity or intensity"));
        }
        Ok(MovingCompactField {
            sample,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
        })
    }

    // Retain the first-order fields through coherent image-path contraction.
    // Squaring individual paths here would omit cross-path energy transport.
    pub(super) fn evaluate_with_fields_dd(
        &self,
        strength: CompactSourceStrength,
        force_history_integral_n_s_m3: [f64; 3],
    ) -> Result<(MovingCompactSample, Dd, [Dd; 3]), AeroacError> {
        self.evaluate_dd(strength, force_history_integral_n_s_m3.map(Dd::from_f64))
    }

    /// Unrounded pressure and velocity for coherent source/image composition.
    ///
    /// The accepted-history owner may retain the past-force integral in DD,
    /// then accumulate every path before rounding, forming intensity or applying
    /// a directional microphone. This evaluates the SAME field as `evaluate`;
    /// it neither supplies missing prehistory nor squares individual sources.
    pub fn evaluate_dd(
        &self,
        strength: CompactSourceStrength,
        force_history_integral_n_s_m3: [Dd; 3],
    ) -> Result<(MovingCompactSample, Dd, [Dd; 3]), AeroacError> {
        if force_history_integral_n_s_m3
            .iter()
            .any(|v| !v.to_f64().is_finite())
        {
            return Err(nonfinite("moving force history integral"));
        }
        let (sample, pressure) = self.pressure.evaluate_with_pressure_dd(strength)?;
        let monopole = self.pressure.monopole_velocity.ok_or_else(|| {
            invalid("moving force-field monopole coefficients are unrepresentable")
        })?;
        let mut force_pressure = Dd::ZERO;
        for (i, value) in strength
            .force_on_fluid_n
            .iter()
            .chain(&strength.force_rate_on_fluid_n_s)
            .enumerate()
        {
            force_pressure = force_pressure
                + Dd::from_f64(self.pressure.coefficients[2 + i]) * Dd::from_f64(*value);
        }
        let velocity: [Dd; 3] = core::array::from_fn(|i| {
            let mut value = Dd::from_f64(monopole[i][0]) * Dd::from_f64(strength.volume_flow_m3_s)
                + Dd::from_f64(monopole[i][1]) * Dd::from_f64(strength.volume_flow_rate_m3_s2);
            value = value + Dd::from_f64(self.force_pressure_velocity[i]) * force_pressure;
            for j in 0..3 {
                value = value
                    + Dd::from_f64(self.force_near_velocity[i][j])
                        * Dd::from_f64(strength.force_on_fluid_n[j]);
            }
            value + Dd::from_f64(self.history_velocity_scale) * force_history_integral_n_s_m3[i]
        });
        if velocity.iter().any(|v| !v.to_f64().is_finite()) {
            return Err(nonfinite("moving compact particle velocity"));
        }
        Ok((sample, pressure, velocity))
    }
}
