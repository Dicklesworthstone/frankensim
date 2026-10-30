//! Physical terminal flow from the same molecular chain as pressure feedback.

use super::MolecularDrivingPoint;
use crate::driving_point::DrivingPointError;
use crate::driving_point::regional::outlet::{MAX_OUTLET_FFT_SAMPLES, OutletError, OutletPressure};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_duct::Termination;
use fs_math::c64::C64;

/// Cold all-mouth observation budget, counting each retained positive-frequency
/// bin once per mouth (including the terminal). It bounds the shared frequency
/// sweep before allocation and never removes a physical opening to fit.
pub const MAX_MOLECULAR_OUTLET_BINS: usize = 1_048_576;

impl MolecularDrivingPoint {
    /// Observe the flanged outlet in the original explicitly bound exterior.
    /// Both the bore's molecular storage and its diffusive walls participate
    /// in the actual outlet flow; no equilibrium outlet transfer is substituted.
    /// The caller must use the same flanged load for its pressure feedback.
    /// This method deliberately observes only the terminal; use
    /// [`Self::baffled_outlet_and_open_holes_delayed_at`] for the full opening set.
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

    /// Terminal and every open/vented chimney, solved in one molecular sweep.
    /// Each mouth uses the same flanged load as feedback and its full authored
    /// radius. Its solved volume flow already contains the vent fraction.
    ///
    /// This retains the existing parallel-baffle geometry convention: a side
    /// receiver is translated by `L - x_hole` along its local x axis. The disk
    /// normals stay parallel and exterior cross-loading is not modeled. These
    /// independent Rayleigh pressures may be summed at the common receiver.
    /// Closed chimneys affect feedback but have no radiating mouth.
    pub fn baffled_outlet_and_open_holes_delayed_at(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        delay_speed_m_s: f64,
    ) -> Result<(OutletPressure, Vec<OutletPressure>), OutletError> {
        self.baffled_outlet_and_open_holes_with_delay_speeds(
            sample_rate_hz,
            fft_samples,
            receiver,
            delay_speed_m_s,
            &vec![delay_speed_m_s; self.open_hole_sites.len()],
        )
    }

    /// The same all-opening response with each exterior propagation clock
    /// explicitly bound. A dispersive observer may have a different FIR delay
    /// at each physical receiver range. Speeds must match the exact authored
    /// open-mouth order; none changes the bore or radiation-load gas state.
    pub fn baffled_outlet_and_open_holes_with_delay_speeds(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        terminal_delay_speed_m_s: f64,
        hole_delay_speeds_m_s: &[f64],
    ) -> Result<(OutletPressure, Vec<OutletPressure>), OutletError> {
        if hole_delay_speeds_m_s.len() != self.open_hole_sites.len()
            || !terminal_delay_speed_m_s.is_finite()
            || terminal_delay_speed_m_s <= 0.0
            || hole_delay_speeds_m_s
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "each molecular opening requires one finite positive exterior delay speed",
            )
            .into());
        }
        if self.open_hole_sites.is_empty() {
            return self
                .baffled_outlet_delayed_at(
                    sample_rate_hz,
                    fft_samples,
                    receiver,
                    terminal_delay_speed_m_s,
                )
                .map(|terminal| (terminal, Vec::new()));
        }
        if sample_rate_hz == 0
            || !(256..=MAX_OUTLET_FFT_SAMPLES).contains(&fft_samples)
            || !fft_samples.is_power_of_two()
        {
            return Err(invalid("molecular opening observation requires a physical clock and power-of-two grid in 256..=65536").into());
        }
        let bins = fft_samples / 2;
        let count = self.open_hole_sites.len() + 1;
        if bins
            .checked_mul(count)
            .is_none_or(|n| n > MAX_MOLECULAR_OUTLET_BINS)
        {
            return Err(invalid(
                "joint molecular opening observation exceeds the retained frequency-bin budget",
            )
            .into());
        }

        // One original chain solve produces every mouth at each frequency.
        // The existing scalar observer then owns differentiation, delay
        // factoring, FFT, Rayleigh quadrature and runtime convolution.
        let mut flows: Vec<Vec<C64>> = (0..count).map(|_| Vec::with_capacity(bins)).collect();
        let rate = f64::from(sample_rate_hz);
        for k in 1..=bins {
            let omega = core::f64::consts::TAU * k as f64 / fft_samples as f64 * rate;
            let (terminal, holes) = self.flanged_opening_flows(omega)?;
            flows[0].push(terminal);
            for (channel, flow) in flows[1..].iter_mut().zip(holes) {
                channel.push(flow);
            }
        }
        let mut channels = flows.into_iter();
        let mut build = |radius, time, receiver, speed| {
            let values = channels.next().expect("one retained channel per opening");
            let mut sample = 0;
            OutletPressure::from_flow_response(
                self.exterior,
                radius,
                sample_rate_hz,
                fft_samples,
                time,
                receiver,
                speed,
                |omega| {
                    sample += 1;
                    let expected =
                        core::f64::consts::TAU * sample as f64 / fft_samples as f64 * rate;
                    if omega != expected {
                        return Err(invalid(
                            "molecular outlet frequency grid changed during shared realization",
                        ));
                    }
                    values.get(sample - 1).copied().ok_or_else(|| {
                        invalid("molecular outlet frequency grid exceeded its retained response")
                    })
                },
            )
        };
        let terminal = build(
            self.geometry
                .segments
                .last()
                .expect("admitted axial terminal")
                .outlet_radius(),
            self.one_way_wavefront_s,
            receiver,
            terminal_delay_speed_m_s,
        )?;
        let length = super::duct_axial_length(&self.geometry);
        let speed = self.material.state().frozen_sound_speed_m_s();
        let mut holes = Vec::with_capacity(self.open_hole_sites.len());
        for (site, &delay_speed) in self.open_hole_sites.iter().zip(hole_delay_speeds_m_s) {
            let mut position = receiver;
            position.position_m[0] += length - site.axial_position_m;
            holes.push(build(
                site.hole_radius_m,
                (site.axial_position_m + site.open_chimney_length_m) / speed,
                position,
                delay_speed,
            )?);
        }
        Ok((terminal, holes))
    }

    // Qout/p+ uses the very same real inlet reference as ReedBoreVoice.
    // A flanged radiation load is dissipative at every sampled positive
    // frequency; original composition/conditioning refusals still propagate.
    pub(super) fn flanged_outlet_flow(&self, omega: f64) -> Result<C64, DrivingPointError> {
        self.flanged_opening_flows(omega)
            .map(|(terminal, _)| terminal)
    }

    // All transfers use the same original scattering solution. In particular,
    // chimney input flow is not substituted for the flow at its exterior mouth.
    fn flanged_opening_flows(&self, omega: f64) -> Result<(C64, Vec<C64>), DrivingPointError> {
        let load = self.terminal_load(omega, Termination::FlangedOpen)?;
        let mouth_loads = self.hole_mouth_loads(omega)?;
        let response =
            self.material
                .response_with_holes(omega, load, &mouth_loads, self.composition)?;
        let radius = self.inlet_radius_m();
        let area = core::f64::consts::PI * radius * radius;
        let gas = self.inlet_gas();
        let zc = gas.density * gas.sound_speed / area;
        let denominator = C64::new(response.duct.z_in.re + zc, response.duct.z_in.im);
        let inlet_flow = denominator.recip().scale(2.0);
        let holes: Vec<_> = response
            .holes
            .iter()
            .filter(|hole| hole.site.open_fraction > 0.0)
            .map(|hole| hole.u_mouth_over_u_in * inlet_flow)
            .collect();
        if holes.len() != self.open_hole_sites.len() {
            return Err(invalid(
                "molecular mouth count changed during frequency evaluation",
            ));
        }
        Ok((
            response.duct.u_out_over_u_in.scale(2.0) * denominator.recip(),
            holes,
        ))
    }
}

fn invalid(what: &'static str) -> DrivingPointError {
    DrivingPointError::Invalid { what }
}

#[cfg(test)]
mod tests;
