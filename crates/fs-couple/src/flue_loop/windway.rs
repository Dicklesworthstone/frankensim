//! Fixed-channel flow and a refinable convection field for jet-card feedback.
//!
//! The channel resolves unsteady transverse no-slip momentum. The feedback
//! field solves dv/dt + c(t) dv/dx = 0 across the supplied cut-up, using the
//! accepted bulk speed times the ORIGINAL card convection ratio. Linear
//! semi-Lagrangian interpolation is refinable in space and time. Zero speed
//! preserves the field; startup has zero initial history, and a pressure release
//! does not erase the moving gas or clamp an ever-growing time delay.
//!
//! This does not turn the surrounding jet-drive card into resolved mouth flow.
//! Its prescribed profile, amplification and receptivity keep their old scope.

use super::{FlueError, FlueGeometry};
use crate::bernoulli_aperture::slit_channel::{
    SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
};
use fs_material::gas::GasState;

/// Maximum spatial cells in the convection field.
pub const MAX_CONVECTION_CELLS: usize = 4096;

/// Actual unsteady source channel and its retained spatial convection state.
#[derive(Debug, Clone)]
pub struct UnsteadyWindway {
    channel: SlitChannel,
    convection: ConvectedVelocity,
    convection_ratio: f64,
}

impl UnsteadyWindway {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(gas: &GasState, geometry: &FlueGeometry, length_m: f64,
        transverse_cells: usize, convection_cells: usize, sample_period_s: f64,
        convection_ratio: f64) -> Result<Self, FlueError>
    {
        let channel = SlitChannel::from_gas(gas, SlitChannelSpec {
            length_m, width_m: geometry.flue_width_m, height_m: geometry.flue_height_m,
            transverse_cells, sample_period_s, exit: SlitExit::UniformKineticHead,
        }).map_err(FlueError::Physics)?;
        let convection = ConvectedVelocity::new(geometry.cut_up_m, convection_cells)?;
        if !convection_ratio.is_finite() || convection_ratio <= 0.0 {
            return Err(FlueError::Invalid { what: "unsteady windway requires positive card convection speed ratio" });
        }
        Ok(Self { channel, convection, convection_ratio })
    }

    /// Source gas, physical geometry, transverse velocity field and energy balance.
    #[must_use]
    pub const fn channel(&self) -> &SlitChannel { &self.channel }

    /// Explicit spatial resolution across the physical mouth cut-up.
    #[must_use]
    pub fn convection_cells(&self) -> usize { self.convection.field.len() }

    /// Cold equilibrium response used for control admission, not live stepping.
    pub(super) fn steady_speed(&self, pressure_pa: f64) -> Result<f64, FlueError> {
        self.channel.steady_bulk_speed_m_s(pressure_pa).map_err(FlueError::Physics)
    }

    pub(super) fn advance(&mut self, pressure_pa: f64, root_velocity_m_s: f64)
        -> Result<(SlitChannelStep, f64), FlueError>
    {
        let frame = self.channel.step(pressure_pa).map_err(FlueError::Physics)?;
        if frame.bulk_speed_m_s < 0.0 {
            return Err(FlueError::Invalid { what: "flue jet requires nonnegative accepted windway speed" });
        }
        let convected = self.convection.advance(frame.bulk_speed_m_s * self.convection_ratio,
            self.channel.spec().sample_period_s, root_velocity_m_s)?;
        Ok((frame, convected))
    }
}

#[derive(Debug, Clone)]
struct ConvectedVelocity {
    dx_m: f64,
    // Nodes dx, 2dx, ... length. The inlet at zero is the supplied boundary.
    field: Vec<f64>,
}

impl ConvectedVelocity {
    fn new(length_m: f64, cells: usize) -> Result<Self, FlueError> {
        let dx_m = length_m / cells as f64;
        if !(2..=MAX_CONVECTION_CELLS).contains(&cells) || !dx_m.is_finite() || dx_m <= 0.0 {
            return Err(FlueError::Invalid { what: "jet convection requires a positive physical length and 2..=4096 cells" });
        }
        Ok(Self { dx_m, field: vec![0.0; cells] })
    }

    fn advance(&mut self, speed_m_s: f64, dt_s: f64, boundary: f64) -> Result<f64, FlueError> {
        let shift = speed_m_s * dt_s / self.dx_m;
        if ![speed_m_s, dt_s, boundary, shift].iter().all(|v| v.is_finite())
            || speed_m_s < 0.0 || dt_s <= 0.0
        {
            return Err(FlueError::Invalid { what: "jet convection requires finite nonnegative travel and finite boundary velocity" });
        }
        // At rest no old parcel is overwritten merely because a callback ran.
        if shift == 0.0 { return Ok(self.field[self.field.len() - 1]); }
        // Upstream source indices are still old state in this descending pass.
        for i in (0..self.field.len()).rev() {
            let from = (i + 1) as f64 - shift;
            self.field[i] = if from <= 0.0 {
                boundary
            } else {
                let lower = from.floor() as usize;
                let fraction = from - lower as f64;
                let left = if lower == 0 { boundary } else { self.field[lower - 1] };
                // from <= i+1; a zero fractional part at the final node needs
                // no out-of-bounds upper sample (very small representable shift).
                if fraction == 0.0 { left } else {
                    let right = self.field[lower];
                    (1.0 - fraction) * left + fraction * right
                }
            };
        }
        Ok(self.field[self.field.len() - 1])
    }
}

#[cfg(test)]
mod tests {
    use super::ConvectedVelocity;

    #[test]
    fn integer_characteristics_cross_exactly_and_stop_preserves_parcels() {
        let mut field = ConvectedVelocity::new(0.004, 8).unwrap();
        // One cell per step: pure characteristic translation, no numerical blur.
        for i in 0..8 {
            let value = field.advance(1.0, 0.0005, 2.0).unwrap();
            assert_eq!(value, if i == 7 { 2.0 } else { 0.0 });
        }
        let retained = field.field.clone();
        for _ in 0..100 { assert_eq!(field.advance(0.0, 0.001, -5.0).unwrap(), 2.0); }
        assert_eq!(field.field, retained);
        // A fresh characteristic must travel the whole actual cut-up on restart.
        for i in 0..8 {
            assert_eq!(field.advance(1.0, 0.0005, -3.0).unwrap(), if i == 7 { -3.0 } else { 2.0 });
        }
    }

    #[test]
    fn fractional_characteristics_converge_to_smooth_transport() {
        // Exact v(x,t)=sin(omega*(t-x/c)) behind a sinusoidal boundary.
        // Refine dt together with dx at fixed fractional Courant number 0.4.
        let mut errors = Vec::new();
        for cells in [16, 32, 64] {
            let mut field = ConvectedVelocity::new(1.0, cells).unwrap();
            let dt = 0.4 / cells as f64;
            let omega = 1.1;
            for n in 0..5 * cells {
                let time = n as f64 * dt;
                field.advance(1.0, dt, (omega * time).sin()).unwrap();
            }
            let exact = (omega * (2.0 - 1.0)).sin();
            errors.push((field.field[cells - 1] - exact).abs());
        }
        assert!(errors[0] > 1e-4, "coarse-grid error must be live");
        assert!(errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1], "{errors:?}");
    }
}
