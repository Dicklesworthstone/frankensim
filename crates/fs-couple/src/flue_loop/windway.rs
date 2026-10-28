//! Fixed-channel flow and a refinable convection field for jet-card feedback.
//!
//! The channel resolves unsteady transverse no-slip momentum. The feedback
//! field solves dv/dt + c(t) dv/dx = 0 across the supplied cut-up, using the
//! accepted bulk speed times the ORIGINAL card convection ratio. Linear
//! semi-Lagrangian interpolation is refinable in space and time. Characteristics
//! entering during a step interpolate the previous and current mouth velocity;
//! their current coefficient participates in the original joint mouth solve.
//! Source momentum and transport publish only after that solve is accepted. Zero speed
//! preserves the field; startup has zero initial history, and a pressure release
//! does not erase the moving gas or clamp an ever-growing time delay.
//! The labium feedback is read at a supplied fraction of the cut-up (the jet
//! card's band-lag compensation), linearly between the bracketing nodes.
//!
//! This does not turn the surrounding jet-drive card into resolved mouth flow.
//! Its prescribed profile, amplification and receptivity keep their old scope.

use super::{FlueError, FlueGeometry, RetardedVelocity};
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
    // Fraction of the cut-up at which the last accepted sample read the feedback.
    tap_fraction: f64,
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
        Ok(Self { channel, convection, convection_ratio, tap_fraction: 1.0 })
    }

    /// Distance from the flue exit at which the last accepted sample read the
    /// labium feedback [m]; the cut-up before any sample.
    #[must_use]
    pub fn feedback_tap_m(&self) -> f64 {
        self.tap_fraction * self.convection.dx_m * self.convection.field.len() as f64
    }

    /// Source gas, physical geometry, transverse velocity field and energy balance.
    #[must_use]
    pub const fn channel(&self) -> &SlitChannel { &self.channel }

    /// Explicit spatial resolution across the physical mouth cut-up.
    #[must_use]
    pub fn convection_cells(&self) -> usize { self.convection.field.len() }

    /// Accepted transverse feedback velocities at x=dx, 2dx, ... cut-up [m/s].
    /// These are the field actually transported to the labium, in physical order.
    #[must_use]
    pub fn convected_velocity_m_s(&self) -> &[f64] { &self.convection.field }

    /// Last accepted mouth velocity at the flue boundary of the transport field.
    #[must_use]
    pub const fn boundary_velocity_m_s(&self) -> f64 { self.convection.boundary_m_s }

    /// Cold equilibrium response used for control admission, not live stepping.
    pub(super) fn steady_speed(&self, pressure_pa: f64) -> Result<f64, FlueError> {
        self.channel.steady_bulk_speed_m_s(pressure_pa).map_err(FlueError::Physics)
    }

    pub(super) fn prepare(&self, pressure_pa: f64) -> Result<WindwayTrial, FlueError>
    {
        // The existing channel is the only momentum owner. Its complete trial
        // remains private until mouth, wave and transport admission succeed.
        let mut channel = self.channel.clone();
        let frame = channel.step(pressure_pa).map_err(FlueError::Physics)?;
        if frame.bulk_speed_m_s < 0.0 {
            return Err(FlueError::Invalid { what: "flue jet requires nonnegative accepted windway speed" });
        }
        let convection = self.convection.prepare(frame.bulk_speed_m_s * self.convection_ratio,
            self.channel.spec().sample_period_s)?;
        Ok(WindwayTrial { channel, convection, convection_ratio: self.convection_ratio,
            tap_fraction: self.tap_fraction })
    }
}

/// A source-fluid candidate and the affine current-mouth trace it needs.
/// Finishing it validates every transported value before any participant commits.
pub(super) struct WindwayTrial {
    channel: SlitChannel,
    convection: ConvectionTrial,
    convection_ratio: f64,
    tap_fraction: f64,
}
impl WindwayTrial {
    pub(super) fn frame(&self) -> SlitChannelStep { self.channel.last_step() }
    /// The feedback at `tap_fraction` of the cut-up (0 is the flue exit, 1 the
    /// labium). The transport step's own linear interpolant between the
    /// bracketing nodes; below the first node the other end is the current mouth
    /// velocity. Every node is affine in that velocity, so the read is too.
    pub(super) fn mouth_response(&mut self, tap_fraction: f64)
        -> Result<RetardedVelocity, FlueError>
    {
        if !(0.0..=1.0).contains(&tap_fraction) {
            return Err(FlueError::Invalid { what: "jet feedback tap must lie on the cut-up" });
        }
        let values = &self.convection.values;
        let position = tap_fraction * values.len() as f64;
        let lower = position.floor() as usize;
        let fraction = position - lower as f64;
        let node = |n: usize| if n == 0 {
            RetardedVelocity { past_m_s: 0.0, current_weight: 1.0 }
        } else { values[n - 1] };
        let left = node(lower);
        let value = if fraction == 0.0 { left } else {
            let right = node(lower + 1);
            RetardedVelocity {
                past_m_s: (1.0 - fraction) * left.past_m_s + fraction * right.past_m_s,
                current_weight: (1.0 - fraction) * left.current_weight
                    + fraction * right.current_weight,
            }
        };
        self.tap_fraction = tap_fraction;
        Ok(value)
    }
    pub(super) fn finish(self, boundary_m_s: f64) -> Result<UnsteadyWindway, FlueError> {
        Ok(UnsteadyWindway { channel: self.channel,
            convection: self.convection.finish(boundary_m_s)?,
            convection_ratio: self.convection_ratio, tap_fraction: self.tap_fraction })
    }
}

#[derive(Debug, Clone)]
struct ConvectedVelocity {
    dx_m: f64,
    // Nodes dx, 2dx, ... length. The inlet at zero is the supplied boundary.
    field: Vec<f64>,
    boundary_m_s: f64,
}

struct ConvectionTrial {
    dx_m: f64,
    values: Vec<RetardedVelocity>,
}
impl ConvectionTrial {
    fn finish(self, boundary_m_s: f64) -> Result<ConvectedVelocity, FlueError> {
        if !boundary_m_s.is_finite() {
            return Err(FlueError::Invalid { what: "jet convection needs a finite accepted mouth velocity" });
        }
        let field: Vec<_> = self.values.iter()
            .map(|value| value.past_m_s + value.current_weight * boundary_m_s).collect();
        if field.iter().any(|value| !value.is_finite()) {
            return Err(FlueError::Invalid { what: "candidate jet convection field is not finite" });
        }
        Ok(ConvectedVelocity { dx_m: self.dx_m, field, boundary_m_s })
    }
}

impl ConvectedVelocity {
    fn new(length_m: f64, cells: usize) -> Result<Self, FlueError> {
        let dx_m = length_m / cells as f64;
        if !(2..=MAX_CONVECTION_CELLS).contains(&cells) || !dx_m.is_finite() || dx_m <= 0.0 {
            return Err(FlueError::Invalid { what: "jet convection requires a positive physical length and 2..=4096 cells" });
        }
        Ok(Self { dx_m, field: vec![0.0; cells], boundary_m_s: 0.0 })
    }

    fn prepare(&self, speed_m_s: f64, dt_s: f64) -> Result<ConvectionTrial, FlueError> {
        let shift = speed_m_s * dt_s / self.dx_m;
        if ![speed_m_s, dt_s, self.boundary_m_s, shift].iter().all(|v| v.is_finite())
            || speed_m_s < 0.0 || dt_s <= 0.0
        {
            return Err(FlueError::Invalid { what: "jet convection requires finite nonnegative travel and finite boundary velocity" });
        }
        let mut values = Vec::with_capacity(self.field.len());
        for i in 0..self.field.len() {
            let position = (i + 1) as f64;
            let from = position - shift;
            let value = if shift == 0.0 {
                // At rest the field stays intact; the boundary still has its
                // own accepted acoustic clock for a later restart.
                RetardedVelocity { past_m_s: self.field[i], current_weight: 0.0 }
            } else if from <= 0.0 {
                // This characteristic entered at t_new-x/c within the step.
                // Linear time interpolation gives an affine CURRENT mouth
                // term, not an extra sample of delay or a boundary held early.
                let current_weight = (shift-position)/shift;
                RetardedVelocity { past_m_s: (1.0-current_weight)*self.boundary_m_s,
                    current_weight }
            } else {
                let lower = from.floor() as usize;
                let fraction = from - lower as f64;
                let left = if lower == 0 { self.boundary_m_s } else { self.field[lower - 1] };
                // from <= i+1; a zero fractional part at the final node needs
                // no out-of-bounds upper sample (very small representable shift).
                let past_m_s = if fraction == 0.0 { left } else {
                    let right = self.field[lower];
                    (1.0 - fraction) * left + fraction * right
                };
                RetardedVelocity { past_m_s, current_weight: 0.0 }
            };
            if !value.past_m_s.is_finite() || !value.current_weight.is_finite()
                || !(0.0..=1.0).contains(&value.current_weight) {
                return Err(FlueError::Invalid { what: "affine jet convection trace is not finite and causal" });
            }
            values.push(value);
        }
        Ok(ConvectionTrial { dx_m: self.dx_m, values })
    }

    #[cfg(test)]
    fn advance(&mut self, speed_m_s: f64, dt_s: f64, boundary_m_s: f64)
        -> Result<f64, FlueError>
    {
        *self = self.prepare(speed_m_s, dt_s)?.finish(boundary_m_s)?;
        Ok(self.field[self.field.len()-1])
    }
}

#[cfg(test)]
mod tests {
    use super::ConvectedVelocity;

    #[test]
    fn integer_characteristics_cross_exactly_and_stop_preserves_parcels() {
        let mut field = ConvectedVelocity::new(0.004, 8).unwrap();
        // Establish the constant inflow at rest, before its first transit.
        field.advance(0.0, 0.0005, 2.0).unwrap();
        // One cell per step: pure characteristic translation, no numerical blur.
        for i in 0..8 {
            let value = field.advance(1.0, 0.0005, 2.0).unwrap();
            assert_eq!(value, if i == 7 { 2.0 } else { 0.0 });
        }
        let retained = field.field.clone();
        for _ in 0..100 { assert_eq!(field.advance(0.0, 0.001, -5.0).unwrap(), 2.0); }
        assert_eq!(field.field, retained);
        // The boundary evolved to -5 while transport stood still. That physical
        // restart parcel arrives first, followed by the new -3 endpoint value.
        for i in 0..8 {
            assert_eq!(field.advance(1.0, 0.0005, -3.0).unwrap(), if i == 7 { -5.0 } else { 2.0 });
        }
        assert_eq!(field.advance(1.0, 0.0005, -3.0).unwrap(), -3.0);
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
                let time = (n+1) as f64 * dt;
                field.advance(1.0, dt, (omega * time).sin()).unwrap();
            }
            let exact = (omega * (2.0 - 1.0)).sin();
            errors.push((field.field[cells - 1] - exact).abs());
        }
        assert!(errors[0] > 1e-4, "coarse-grid error must be live");
        assert!(errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1], "{errors:?}");
    }

    #[test]
    fn linear_characteristics_use_the_actual_boundary_time_and_joint_current_coefficient() {
        for (cells, speed, dt) in [(4, 4.0, 0.3), (16, 1.0, 0.1)] {
            let mut field = ConvectedVelocity::new(1.0, cells).unwrap();
            field.boundary_m_s = 3.0;
            // Exact v(x,t) = 3 + 2*(t-x/c), including the initial interior.
            for (i, value) in field.field.iter_mut().enumerate() {
                *value = 3.0 - 2.0*(i+1) as f64*field.dx_m/speed;
            }
            let old = field.field.clone();
            let candidate = field.prepare(speed, dt).unwrap();
            let outlet = candidate.values[cells-1];
            let current = 3.0+2.0*dt;
            let exact = 3.0+2.0*(dt-1.0/speed);
            assert!((outlet.past_m_s+outlet.current_weight*current-exact).abs() < 1e-14);
            assert!((outlet.current_weight-(1.0-1.0/(speed*dt)).max(0.0)).abs() < 1e-14);
            assert_eq!(field.field, old, "preparing a joint solve cannot advance accepted history");
            let accepted = candidate.finish(current).unwrap();
            for (i, value) in accepted.field.iter().enumerate() {
                let exact = 3.0+2.0*(dt-(i+1) as f64*field.dx_m/speed);
                assert!((value-exact).abs() < 1e-14);
            }
            assert_eq!(accepted.boundary_m_s, current);
        }
    }

    #[test]
    fn rejected_source_and_transport_candidates_leave_both_physical_histories_untouched() {
        use super::UnsteadyWindway;
        use crate::flue_loop::FlueGeometry;
        use fs_material::gas::{GasSpec, GasState};
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let geometry = FlueGeometry { flue_width_m: 0.01, flue_height_m: 0.001,
            cut_up_m: 0.004, labium_offset_m: 0.0, jet_angle_rad: 0.0,
            bore_area_m2: 1e-4, mouth_loss_rel: 0.02 };
        let original = UnsteadyWindway::new(&gas, &geometry, 0.015, 16, 32,
            1.0/48_000.0, 0.4).unwrap();
        let candidate = original.prepare(300.0).unwrap();
        assert!(candidate.frame().bulk_speed_m_s > 0.0);
        assert!(candidate.frame().kinetic_energy_j > 0.0);
        assert!(candidate.finish(f64::INFINITY).is_err());
        assert!(original.channel().velocities_m_s().iter().all(|v| *v == 0.0));
        assert!(original.convected_velocity_m_s().iter().all(|v| *v == 0.0));
        assert_eq!(original.boundary_velocity_m_s(), 0.0);
        assert_eq!(original.channel().last_step().kinetic_energy_j, 0.0);
        let accepted = original.prepare(300.0).unwrap().finish(0.3).unwrap();
        assert!(accepted.channel().last_step().kinetic_energy_j > 0.0);
        assert_eq!(accepted.boundary_velocity_m_s(), 0.3);
    }

    #[test]
    fn feedback_tap_reads_the_affine_field_anywhere_on_the_cut_up() {
        // The frozen path's band-lag read is exact on an affine signal (flue_loop's
        // g1_retarded_affine_signal_has_no_extra_sample_at_any_ring_position, 64 epsilon).
        // Exact v(x,t) = 3 + 2*(t-x/c): a tap at x must return 3 + 2*(dt-x/c) at the exit,
        // below the first node, at and between nodes and at the labium, with the current
        // mouth coefficient the characteristic geometry gives.
        use super::UnsteadyWindway;
        use crate::flue_loop::FlueGeometry;
        use fs_material::gas::{GasSpec, GasState};
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let geometry = FlueGeometry { flue_width_m: 0.01, flue_height_m: 0.001,
            cut_up_m: 0.004, labium_offset_m: 0.0, jet_angle_rad: 0.0,
            bore_area_m2: 1e-4, mouth_loss_rel: 0.02 };
        let mut windway = UnsteadyWindway::new(&gas, &geometry, 0.015, 16, 16,
            1.0/48_000.0, 0.4).unwrap();
        assert_eq!(windway.feedback_tap_m(), windway.convection.dx_m * 16.0);
        // dx = 0.25 mm. At 40 m/s and 48 kHz a characteristic crosses 3.3 cells per step,
        // so the exit and the first three nodes carry a current-mouth coefficient.
        let (speed, dt) = (40.0, 1.0/48_000.0);
        windway.convection.boundary_m_s = 3.0;
        let dx = windway.convection.dx_m;
        for (i, value) in windway.convection.field.iter_mut().enumerate() {
            *value = 3.0 - 2.0*(i+1) as f64*dx/speed;
        }
        let current = 3.0 + 2.0*dt;
        for fraction in [0.0, 0.03125, 0.0625, 0.125, 0.15, 0.775, 0.765, 1.0] {
            let mut trial = super::WindwayTrial { channel: windway.channel.clone(),
                convection: windway.convection.prepare(speed, dt).unwrap(),
                convection_ratio: 0.4, tap_fraction: 1.0 };
            let read = trial.mouth_response(fraction).unwrap();
            let x = fraction * geometry.cut_up_m;
            let exact = 3.0 + 2.0*(dt - x/speed);
            assert!((read.past_m_s + read.current_weight*current - exact).abs()
                < 64.0*f64::EPSILON*exact.abs().max(1.0), "fraction {fraction}");
            assert!((read.current_weight - (1.0 - x/(speed*dt)).max(0.0)).abs() < 1e-14,
                "fraction {fraction}: {}", read.current_weight);
            let accepted = trial.finish(current).unwrap();
            assert!((accepted.feedback_tap_m() - x).abs() < 1e-18);
        }
        let mut trial = super::WindwayTrial { channel: windway.channel.clone(),
            convection: windway.convection.prepare(speed, dt).unwrap(),
            convection_ratio: 0.4, tap_fraction: 1.0 };
        for fraction in [-1e-12, 1.0 + 1e-12, f64::NAN] {
            assert!(trial.mouth_response(fraction).is_err(), "fraction {fraction}");
        }
    }
}
