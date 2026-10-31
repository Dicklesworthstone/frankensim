//! Compact, stationary 3D force dipoles in a homogeneous lossless fluid.
//!
//! For actual force **on the fluid** `F` [N], the retarded free-space solution is
//! `p(x,t) = n·F(t-r/c)/(4πr²) + n·Fdot(t-r/c)/(4πcr)`.
//! This is the stationary-source specialization of Farassat, NASA TP-3428,
//! Eq. 4.55: <https://ntrs.nasa.gov/citations/19940029887>.
//! Density already enters the physical force; multiplying this solution by
//! density again is dimensionally incorrect. The material's sound speed still
//! sets propagation time and the radiative coefficient.
//!
//! Inputs are prescribed 3D forces, never 2D force-per-span data. This operator
//! does not validate a force model, resolve finite source size, moving sources,
//! boundaries, absorption or feedback, or confer measured-SPL authority.

use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::det;

/// Geometry and material coefficients prepared once for a stationary force.
/// Evaluation is allocation-free with fixed-order three-component reductions.
/// The caller owns source history and any sampled propagation-delay filter.
#[derive(Clone, Copy, Debug)]
pub struct CompactForceObserver {
    acoustic: FluidAcousticState,
    direction: [f64; 3],
    near_scale: f64,
    far_scale: f64,
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
        Ok(Self {
            acoustic,
            direction: displacement.map(|v| v / range),
            near_scale,
            far_scale,
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

    /// Full pressure [Pa], including the `1/r²` induction term. Both arguments
    /// must describe the SAME retarded source time `t-r/c`, in N and N/s.
    /// No time shift or force differentiation is performed here.
    pub fn pressure_at_retarded(
        &self,
        force_n: [f64; 3],
        force_rate_n_s: [f64; 3],
    ) -> Result<f64, AeroacError> {
        let near = self.radial(force_n)? * self.near_scale;
        finite_pressure(near + self.far_field_at_retarded(force_rate_n_s)?)
    }

    /// Radiative `1/r` pressure [Pa] with signed dipole directivity. Supply the
    /// force derivative [N/s] at `t-r/c`. Omitting the induction term is a
    /// far-field approximation; this is not an omnidirectional monopole.
    pub fn far_field_at_retarded(&self, force_rate_n_s: [f64; 3]) -> Result<f64, AeroacError> {
        finite_pressure(self.radial(force_rate_n_s)? * self.far_scale)
    }

    fn radial(&self, vector: [f64; 3]) -> Result<f64, AeroacError> {
        if vector.iter().any(|v| !v.is_finite()) {
            return Err(AeroacError::NonFinite {
                what: "compact-force source",
            });
        }
        Ok(
            (self.direction[0] * vector[0] + self.direction[1] * vector[1])
                + self.direction[2] * vector[2],
        )
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
