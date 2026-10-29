//! Physical terminal flow from the same molecular chain as pressure feedback.

use super::MolecularDrivingPoint;
use crate::driving_point::DrivingPointError;
use crate::driving_point::regional::outlet::{OutletError, OutletPressure};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_duct::Termination;
use fs_math::c64::C64;

impl MolecularDrivingPoint {
    /// Observe the flanged outlet in the original explicitly bound exterior.
    /// Both the bore's molecular storage and its diffusive walls participate
    /// in the actual outlet flow; no equilibrium outlet transfer is substituted.
    /// The caller must use the same flanged load for its pressure feedback.
    ///
    /// `fft_samples` is the already-admitted characteristic line's grid. The
    /// shared outlet owner retains its finite-history and Rayleigh quadrature
    /// approximations, explicit axial/exterior delays and allocation-free step.
    pub fn baffled_outlet(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
    ) -> Result<OutletPressure, OutletError> {
        self.baffled_outlet_delayed_at(
            sample_rate_hz,
            fft_samples,
            receiver,
            self.exterior.sound_speed,
        )
    }

    /// The same outlet with a declared exterior wavefront speed for its delays.
    /// A subsequent dispersive observer law owns any additional excess delay;
    /// the bore material and its frozen wavefront are unchanged.
    pub fn baffled_outlet_delayed_at(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        delay_speed_m_s: f64,
    ) -> Result<OutletPressure, OutletError> {
        let radius = self
            .geometry
            .segments
            .last()
            .expect("admitted nonempty molecular duct")
            .outlet_radius();
        OutletPressure::from_flow_response(
            self.exterior,
            radius,
            sample_rate_hz,
            fft_samples,
            self.one_way_wavefront_s,
            receiver,
            delay_speed_m_s,
            |omega| self.flanged_outlet_flow(omega),
        )
    }

    // Qout/p+ uses the very same real inlet reference as ReedBoreVoice.
    // A flanged radiation load is dissipative at every sampled positive
    // frequency; original composition/conditioning refusals still propagate.
    pub(super) fn flanged_outlet_flow(&self, omega: f64) -> Result<C64, DrivingPointError> {
        let radius = self
            .geometry
            .segments
            .last()
            .expect("admitted nonempty molecular duct")
            .outlet_radius();
        let load = fs_duct::compact_radiation_impedance(
            Termination::FlangedOpen,
            &self.exterior,
            radius,
            omega,
        )?;
        let response = self
            .material
            .response(omega, Some(load), self.composition)?;
        let radius = self.inlet_radius_m();
        let area = core::f64::consts::PI * radius * radius;
        let gas = self.inlet_gas();
        let zc = gas.density * gas.sound_speed / area;
        let denominator = C64::new(response.z_in.re + zc, response.z_in.im);
        Ok(response.u_out_over_u_in.scale(2.0) * denominator.recip())
    }
}

#[cfg(test)]
mod tests;
