//! Compact 3D source operators in a homogeneous lossless fluid.
//!
//! For actual force **on the fluid** `F` [N], the retarded free-space solution is
//! `p(x,t) = n·F(t-r/c)/(4πr²) + n·Fdot(t-r/c)/(4πcr)`.
//! This is the stationary-source specialization of Farassat, NASA TP-3428,
//! Eq. 4.55: <https://ntrs.nasa.gov/citations/19940029887>.
//! Density already enters the physical force; multiplying this solution by
//! density again is dimensionally incorrect. The material's sound speed still
//! sets propagation time and the radiative coefficient.
//!
//! Inputs are prescribed 3D forces, never 2D force-per-span data. The stationary
//! [`CompactForceObserver`] below
//! does not validate a force model, resolve finite source size, moving sources,
//! boundaries, absorption or feedback, or confer measured-SPL authority.
//! The separate [`moving`] owner observes prescribed subsonic source motion.

/// Prepared harmonic arrays of physical volume sources and force dipoles.
pub mod array;
/// Retarded pressure from prescribed moving compact volume and force sources.
pub mod moving;
/// Analytic whole-sphere power of coherent compact-source arrays.
pub mod power;
mod statistics;
pub use array::{
    CoherentCompactField, CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactSourceObserver, PreparedHalfSpaceObserver,
};
pub use array::convected::{
    ConvectedMonopoleField, ConvectedMonopoleLimits, ConvectedMonopoleStatistics,
    PreparedConvectedMonopoleObserver,
};
pub use power::{
    CompactRadiatedPower, PreparedCompactRadiatedPower, RadiatedPowerBudget,
    RadiatedPowerModeBudget,
};
pub use statistics::{CoherenceBudget, CompactFieldStatistics};

use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{dd::Dd, det};

/// The instantaneous field of one prescribed stationary compact force source.
///
/// These are real time-domain values, not peak phasors or RMS amplitudes. At a
/// pressure null the transverse particle velocity may still be nonzero. The
/// intensity is the signed product `p u`; it is neither `p²/(rho c)` nor a
/// harmonic half-product. For multiple coherent sources, sum their pressure
/// and velocity before forming the total intensity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedForceField {
    /// Acoustic pressure perturbation [Pa].
    pub pressure_pa: f64,
    /// Particle velocity [m/s] in the original Cartesian frame.
    pub particle_velocity_m_s: [f64; 3],
    /// Instantaneous acoustic energy flux `p u` [W/m²], including near fields.
    pub instantaneous_intensity_w_m2: [f64; 3],
}

/// Geometry and material coefficients prepared once for a stationary force.
/// Evaluation is allocation-free with fixed-order three-component reductions.
/// The caller owns source history and any sampled propagation-delay filter.
#[derive(Clone, Copy, Debug)]
pub struct CompactForceObserver {
    acoustic: FluidAcousticState,
    observer_m: [f64; 3],
    source_m: [f64; 3],
    direction: [f64; 3],
    near_scale: f64,
    far_scale: f64,
    // Optional so the wider pressure-only domain is not narrowed by a
    // nonrepresentable 1/r³ velocity coefficient.
    velocity_scales: Option<[f64; 3]>,
    delay_s: f64,
}

impl CompactForceObserver {
    /// Prepare one free-space source/observer pair in metres in the same frame.
    /// Refuses coincident positions, non-finite geometry, and coefficients
    /// outside the representable positive finite range.
    pub fn new(
        acoustic: FluidAcousticState,
        observer_m: [f64; 3],
        source_m: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if observer_m.iter().chain(&source_m).any(|v| !v.is_finite()) {
            return Err(AeroacError::NonFinite {
                what: "compact-force positions",
            });
        }
        let displacement = core::array::from_fn::<_, 3, _>(|i| observer_m[i] - source_m[i]);
        let range = det::hypot(
            det::hypot(displacement[0], displacement[1]),
            displacement[2],
        );
        let inv_four_pi_r = (0.25 / core::f64::consts::PI) / range;
        let near_scale = inv_four_pi_r / range;
        let far_scale = inv_four_pi_r / acoustic.sound_speed_m_s();
        let delay_s = range / acoustic.sound_speed_m_s();
        if [range, near_scale, far_scale, delay_s]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(AeroacError::InvalidParameter {
                what: "compact-force separation and propagation must be positive finite",
            });
        }
        let impedance = acoustic.characteristic_impedance_pa_s_m();
        let velocity_scales = [
            (near_scale / range) / acoustic.density_kg_m3(),
            near_scale / impedance,
            far_scale / impedance,
        ];
        let velocity_scales = velocity_scales
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            .then_some(velocity_scales);
        Ok(Self {
            acoustic,
            observer_m,
            source_m,
            direction: displacement.map(|v| v / range),
            near_scale,
            far_scale,
            velocity_scales,
            delay_s,
        })
    }

    /// Original material response used to prepare the observer.
    #[must_use]
    pub const fn acoustic_state(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Physical propagation time `r/c` [s]. Apply it exactly once: evaluate the
    /// source at `observer_time - delay_s`, or use the existing delay filter.
    #[must_use]
    pub const fn delay_s(&self) -> f64 {
        self.delay_s
    }

    /// Prepare this same material and geometry for a stationary harmonic source.
    ///
    /// The resulting operator retains pressure, all three particle-velocity
    /// components and active/reactive intensity. Supply a
    /// [`CompactSourcePhasor::force_dipole`] with physical force on the fluid [N].
    /// The source phasor uses exp(-i omega t) at the original source clock;
    /// the outgoing transfer already includes the travel phase. Do not delay
    /// the phasor again. Preparation allocates once; subsequent evaluation does
    /// not allocate. This does not validate compactness or a force source model.
    pub fn prepare_harmonic(
        &self,
        frequency_hz: f64,
    ) -> Result<PreparedCompactSourceObserver, AeroacError> {
        PreparedCompactSourceObserver::new(
            self.acoustic,
            frequency_hz,
            self.observer_m,
            &[self.source_m],
            1,
        )
    }

    /// Full pressure [Pa], including the `1/r²` induction term. Both arguments
    /// must describe the SAME retarded source time `t-r/c`, in N and N/s.
    /// No time shift or force differentiation is performed here.
    pub fn pressure_at_retarded(
        &self,
        force_n: [f64; 3],
        force_rate_n_s: [f64; 3],
    ) -> Result<f64, AeroacError> {
        self.pressure_from_radial(self.radial(force_n)?, self.radial(force_rate_n_s)?)
    }

    /// Full transient pressure, vector particle velocity and instantaneous flux.
    ///
    /// All three inputs belong to the SAME retarded source time `t-r/c`:
    /// physical force on the fluid [N], its time derivative [N/s], and its
    /// time integral [N s]. The impulse must retain the original source history:
    /// for a causal source initially at rest it is `integral(-infinity,t-r/c) F`.
    /// A different integration constant declares a different initial near field;
    /// resetting it at block boundaries changes the physics. A steady harmonic
    /// source instead supplies the matching zero-mean harmonic antiderivative.
    /// This observer does not integrate, differentiate, delay or resample inputs.
    ///
    /// Away from the source, integrating `rho du/dt = -grad(p)` gives
    /// `u = [(3nnᵀ-I) impulse/r³ + (3nnᵀ-I) F/(c r²)
    ///       + nnᵀ Fdot/(c² r)]/(4 pi rho)`.
    /// The `1/r³` and `1/r²` terms retain transverse near-field motion;
    /// `u=n p/(rho c)` is valid only for the radiative terms. See de Hoop,
    /// *Handbook of Radiation and Scattering of Waves*, §5.8, Eqs. 5.8-21–25:
    /// <https://www.atdehoop.com/Books/DeHoopAT_book_2008/DeHoopAT1995Handbook_chap05.pdf>.
    ///
    /// The same stationary, compact, linear, homogeneous lossless model limits
    /// as the pressure owner apply. Source validity, nonlinear propagation,
    /// scattering and radiation reaction remain outside this one-way transfer.
    /// Evaluation is allocation-free, uses cached coefficients and fixed-order
    /// double-double reductions, and refuses nonfinite inputs/outputs or
    /// unrepresentable velocity coefficients without publishing a partial field.
    pub fn field_at_retarded(
        &self,
        force_impulse_n_s: [f64; 3],
        force_n: [f64; 3],
        force_rate_n_s: [f64; 3],
    ) -> Result<RetardedForceField, AeroacError> {
        let [near, intermediate, radiative] =
            self.velocity_scales.ok_or(AeroacError::InvalidParameter {
                what: "compact-force velocity coefficients are not representable",
            })?;
        let radial_impulse = self.radial(force_impulse_n_s)?;
        let radial_force = self.radial(force_n)?;
        let radial_rate = self.radial(force_rate_n_s)?;
        let pressure_pa = self.pressure_from_radial(radial_force, radial_rate)?;
        let particle_velocity_m_s = core::array::from_fn(|axis| {
            let n = Dd::from_f64(self.direction[axis]);
            let impulse =
                n * radial_impulse * Dd::from_f64(3.0) - Dd::from_f64(force_impulse_n_s[axis]);
            let force = n * radial_force * Dd::from_f64(3.0) - Dd::from_f64(force_n[axis]);
            (impulse * Dd::from_f64(near)
                + force * Dd::from_f64(intermediate)
                + n * radial_rate * Dd::from_f64(radiative))
            .to_f64()
        });
        let instantaneous_intensity_w_m2 = particle_velocity_m_s.map(|v| pressure_pa * v);
        if particle_velocity_m_s
            .iter()
            .chain(&instantaneous_intensity_w_m2)
            .any(|v| !v.is_finite())
        {
            return Err(AeroacError::NonFinite {
                what: "compact-force particle velocity or instantaneous intensity",
            });
        }
        Ok(RetardedForceField {
            pressure_pa,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
        })
    }

    /// Radiative `1/r` pressure [Pa] with signed dipole directivity. Supply the
    /// force derivative [N/s] at `t-r/c`. Omitting the induction term is a
    /// far-field approximation; this is not an omnidirectional monopole.
    pub fn far_field_at_retarded(&self, force_rate_n_s: [f64; 3]) -> Result<f64, AeroacError> {
        finite_pressure((self.radial(force_rate_n_s)? * Dd::from_f64(self.far_scale)).to_f64())
    }

    fn pressure_from_radial(&self, force: Dd, rate: Dd) -> Result<f64, AeroacError> {
        finite_pressure(
            (force * Dd::from_f64(self.near_scale) + rate * Dd::from_f64(self.far_scale)).to_f64(),
        )
    }

    fn radial(&self, vector: [f64; 3]) -> Result<Dd, AeroacError> {
        if vector.iter().any(|v| !v.is_finite()) {
            return Err(AeroacError::NonFinite {
                what: "compact-force source",
            });
        }
        Ok((0..3).fold(Dd::ZERO, |sum, axis| {
            sum + Dd::from_f64(self.direction[axis]) * Dd::from_f64(vector[axis])
        }))
    }
}

fn finite_pressure(pressure: f64) -> Result<f64, AeroacError> {
    if pressure.is_finite() {
        Ok(pressure)
    } else {
        Err(AeroacError::NonFinite {
            what: "compact-force pressure",
        })
    }
}
