//! Original fixed-channel momentum under exp(-i omega t).
//!
//! This borrowed view exposes the already assembled finite-volume masses,
//! pressure/flow areas and no-slip face conductances. It adds no spatial
//! discretization, time integrator, fitted inertance or prescribed profile.
//! The full descriptor is (s M + G) v = A dp, Q = A^T v. Every original
//! transverse/radial velocity remains a physical state.

use super::{AcousticRealizeError, invalid};
use fs_math::c64::C64;

/// Borrowed continuous physical coefficients of the unchanged channel owner.
/// Valid for fixed geometry, the retained coherent Newtonian gas, and the
/// original incompressible, fully developed, laminar reduction about rest.
/// The nonlinear uniform-exit-head option cannot produce this view.
#[derive(Clone, Copy, Debug)]
pub struct ChannelContinuousOperator<'a> {
    mass: &'a [f64],
    area: &'a [f64],
    conductance: &'a [f64],
}

/// Complete peak velocity phasors and independently observed physical work.
/// There is no time-discretization loss in this continuous frequency query.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelHarmonicResponse {
    pub velocities_m_s: Vec<C64>,
    pub pressure_drop_pa: C64,
    pub flow_m3_s: C64,
    /// Cycle-mean kinetic storage, sum(M_i*|v_i|^2)/4 [J].
    pub mean_stored_energy_j: f64,
    /// Cycle-mean physical no-slip friction heat [W].
    pub viscous_dissipation_w: f64,
    /// Cycle-mean pressure work entering this channel [W].
    pub input_power_w: f64,
    /// Wall heat minus pressure work [W]; mean storage change is zero.
    pub power_residual_w: f64,
    /// Largest original cell-equation residual relative to uncancelled terms.
    pub solve_relative_residual: f64,
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

impl<'a> ChannelContinuousOperator<'a> {
    pub(super) fn from_parts(mass: &'a [f64], area: &'a [f64], conductance: &'a [f64]) -> Self {
        Self {
            mass,
            area,
            conductance,
        }
    }
    pub fn dimension(self) -> usize {
        self.mass.len()
    }
    /// Original physical fluid mass in each velocity cell [kg].
    pub const fn cell_masses_kg(self) -> &'a [f64] {
        self.mass
    }
    /// Original area mapping pressure to force and velocity to flow [m2].
    pub const fn cell_areas_m2(self) -> &'a [f64] {
        self.area
    }
    /// Original face shear conductance [kg/s]. Face 0 and face n meet the
    /// stationary boundaries; the circular axis has exactly zero conductance.
    /// The tridiagonal G has diagonal g_i+g_(i+1) and neighbours -g_i.
    pub const fn face_conductances_kg_s(self) -> &'a [f64] {
        self.conductance
    }

    /// Observe a candidate full field against the ORIGINAL continuous cell
    /// equations and physical energy identity. No state changes. The supplied
    /// pressure is the independently observed pressure-plane difference.
    /// A strictly positive frequency is required for these cycle means;
    /// a constant DC field does not have the peak-phasor averaging factor.
    /// Work admission precedes the returned field allocation; a conservative
    /// 256*n arithmetic-work envelope is charged and cancellation is checked per cell.
    pub fn observe_harmonic(
        self,
        omega_rad_s: f64,
        velocities_m_s: &[C64],
        pressure_drop_pa: C64,
        maximum_work: usize,
        maximum_relative_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelHarmonicResponse, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let n = self.dimension();
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || velocities_m_s.len() != n
            || velocities_m_s.iter().any(|&v| !finite(v))
            || !finite(pressure_drop_pa)
            || !maximum_relative_residual.is_finite()
            || maximum_relative_residual <= 0.0
            || maximum_relative_residual >= 1.0
            || n.checked_mul(256).is_none_or(|work| work > maximum_work)
        {
            return Err(invalid(
                "continuous channel observation needs positive finite frequency, full field, tolerance and declared work",
            ));
        }
        let tolerance = maximum_relative_residual + 64.0 * n as f64 * f64::EPSILON;
        let mut flow = C64::ZERO;
        let mut energy = 0.0;
        let mut heat = 0.0;
        let mut largest_residual = 0.0_f64;
        for i in 0..n {
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
            let velocity = velocities_m_s[i];
            let left = if i == 0 {
                C64::ZERO
            } else {
                velocities_m_s[i - 1]
            };
            let right = if i + 1 == n {
                C64::ZERO
            } else {
                velocities_m_s[i + 1]
            };
            let inertia = C64::new(0.0, -omega_rad_s * self.mass[i]) * velocity;
            let down = (velocity - left).scale(self.conductance[i]);
            let up = (velocity - right).scale(self.conductance[i + 1]);
            let forcing = pressure_drop_pa.scale(self.area[i]);
            let residual = inertia + down + up - forcing;
            let scale = inertia.abs()
                + forcing.abs()
                + self.conductance[i] * (velocity.abs() + left.abs())
                + self.conductance[i + 1] * (velocity.abs() + right.abs());
            if !finite(residual)
                || !scale.is_finite()
                || residual.abs() > tolerance * scale.max(1e-300)
            {
                return Err(invalid(
                    "continuous channel field violates its original cell momentum equations",
                ));
            }
            largest_residual = largest_residual.max(residual.abs() / scale.max(1e-300));
            flow = flow + velocity.scale(self.area[i]);
            energy += 0.25 * self.mass[i] * velocity.norm_sq();
            heat += 0.5 * self.conductance[i] * (velocity - left).norm_sq();
        }
        heat += 0.5 * self.conductance[n] * velocities_m_s[n - 1].norm_sq();
        let complex_power = (pressure_drop_pa * flow.conj()).scale(0.5);
        let power_residual = heat - complex_power.re;
        // Imaginary pressure work is the same original kinetic storage:
        // Im(dp*conj(Q)/2) = -2*omega*mean(E) under exp(-i omega t).
        let reactive_residual = complex_power.im + 2.0 * omega_rad_s * energy;
        let power_scale = heat + complex_power.abs() + 2.0 * omega_rad_s * energy;
        if !finite(flow)
            || !finite(complex_power)
            || ![energy, heat, power_residual, reactive_residual, power_scale]
                .iter()
                .all(|x| x.is_finite())
            || power_residual.abs() > tolerance * power_scale.max(1e-300)
            || reactive_residual.abs() > tolerance * power_scale.max(1e-300)
        {
            return Err(invalid(
                "continuous channel pressure work disagrees with original fluid storage and wall heat",
            ));
        }
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        Ok(ChannelHarmonicResponse {
            velocities_m_s: velocities_m_s.to_vec(),
            pressure_drop_pa,
            flow_m3_s: flow,
            mean_stored_energy_j: energy,
            viscous_dissipation_w: heat,
            input_power_w: complex_power.re,
            power_residual_w: power_residual,
            solve_relative_residual: largest_residual,
        })
    }
}
