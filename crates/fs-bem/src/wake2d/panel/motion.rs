//! Prescribed planar rigid motion through the original coupled wake owner.

use super::{
    BemError, IncidentField, PanelWakeAdaptiveOptions, PanelWakeAdaptiveReport,
    PanelWakeTransportOptions, PanelWakeTransportReport, WakeSim, WakeVortex, checkpoint, reserved,
};
use crate::panel2d::Airfoil2d;
use fs_ga::{So3, So3Tangent, Vec3};

/// Pose and velocity of the normalized body chord frame at one instant.
///
/// Lengths use a fixed reference chord `c`, time uses `c/U_ref`, and linear
/// velocity uses `U_ref`. The origin is the leading edge, +x points to the
/// trailing edge, and positive angle/rate is counterclockwise. The supplied
/// velocity and rate must be derivatives of this same prescribed pose; the
/// wake solver neither integrates nor certifies the motion law.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PanelFrame2d {
    /// Leading-edge position in fixed world coordinates, divided by `c`.
    pub origin_world: [f64; 2],
    /// Unwrapped body-to-world angle [rad].
    pub angle_rad: f64,
    /// World derivative of the leading-edge position, divided by `U_ref`.
    pub origin_velocity_world: [f64; 2],
    /// Angle derivative [rad per nondimensional time], physical omega*c/U_ref.
    pub angular_velocity: f64,
}

fn invalid(operation: &'static str) -> BemError {
    BemError::NonFiniteResult { operation }
}

#[derive(Clone, Copy)]
struct PreparedFrame {
    state: PanelFrame2d,
    rotation: So3,
    inverse: So3,
}

impl PreparedFrame {
    fn new(state: PanelFrame2d) -> Result<Self, BemError> {
        if !state
            .origin_world
            .into_iter()
            .chain(state.origin_velocity_world)
            .chain([state.angle_rad, state.angular_velocity])
            .all(f64::is_finite)
        {
            return Err(invalid("panel wake prescribed motion"));
        }
        let rotation = So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, state.angle_rad)))
            .map_err(|_| invalid("panel wake prescribed rotation"))?;
        Ok(Self {
            state,
            rotation,
            inverse: rotation.inverse(),
        })
    }

    fn world_point(self, point: [f64; 2]) -> Result<[f64; 2], BemError> {
        let rotated = if self.state.angle_rad == 0.0 {
            point
        } else {
            let value = self
                .rotation
                .rotate(Vec3::new(point[0], point[1], 0.0))
                .map_err(|_| invalid("panel wake world point"))?;
            [value.x, value.y]
        };
        let world = core::array::from_fn(|axis| {
            if self.state.origin_world[axis] == 0.0 {
                rotated[axis]
            } else {
                self.state.origin_world[axis] + rotated[axis]
            }
        });
        if !world.into_iter().all(f64::is_finite) {
            return Err(invalid("panel wake world point"));
        }
        Ok(world)
    }

    fn relative_velocity(self, point: [f64; 2], air: [f64; 2]) -> Result<[f64; 2], BemError> {
        // Subtract common translation before rotating, so a large common
        // world boost cannot erase a small, physically resolved slip.
        let slip = core::array::from_fn::<_, 2, _>(|axis| {
            if self.state.origin_velocity_world[axis] == 0.0 {
                air[axis]
            } else {
                air[axis] - self.state.origin_velocity_world[axis]
            }
        });
        let mut local = if self.state.angle_rad == 0.0 {
            slip
        } else {
            let value = self
                .inverse
                .rotate(Vec3::new(slip[0], slip[1], 0.0))
                .map_err(|_| invalid("panel wake frame-relative velocity"))?;
            [value.x, value.y]
        };
        if self.state.angular_velocity != 0.0 {
            local[0] += self.state.angular_velocity * point[1];
            local[1] -= self.state.angular_velocity * point[0];
        }
        if !local.into_iter().all(f64::is_finite) {
            return Err(invalid("panel wake frame-relative velocity"));
        }
        Ok(local)
    }
}

type Motion<'a> = dyn FnMut(f64) -> Result<PanelFrame2d, BemError> + 'a;
type WorldAir<'a> = dyn FnMut([f64; 2], f64) -> Result<[f64; 2], BemError> + 'a;

pub(super) struct MovingField<'a> {
    motion: &'a mut Motion<'a>,
    air: &'a mut WorldAir<'a>,
    start: Option<(f64, PanelFrame2d)>,
    cached: Option<(f64, PreparedFrame)>,
}

impl<'a> MovingField<'a> {
    pub(super) fn new(motion: &'a mut Motion<'a>, air: &'a mut WorldAir<'a>) -> Self {
        Self {
            motion,
            air,
            start: None,
            cached: None,
        }
    }

    pub(super) fn begin(
        &mut self,
        time: f64,
        accepted: PanelFrame2d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), BemError> {
        self.start = Some((time, accepted));
        self.frame(time, cancelled).map(|_| ())
    }

    fn frame(
        &mut self,
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedFrame, BemError> {
        if let Some((cached_time, frame)) = self.cached {
            if time == cached_time {
                return Ok(frame);
            }
        }
        if !time.is_finite() {
            return Err(invalid("panel wake motion time"));
        }
        checkpoint(cancelled)?;
        let frame = PreparedFrame::new((self.motion)(time)?)?;
        checkpoint(cancelled)?;
        if let Some((start_time, accepted)) = self.start {
            if time == start_time
                && (frame.state.origin_world != accepted.origin_world
                    || frame.state.angle_rad != accepted.angle_rad)
            {
                return Err(BemError::InvalidScalar {
                    name: "panel wake motion continuity",
                    value: time,
                    requirement: "start origin and unwrapped angle must match the accepted frame",
                });
            }
        }
        self.cached = Some((time, frame));
        Ok(frame)
    }

    pub(super) fn velocity(
        &mut self,
        point: [f64; 2],
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], BemError> {
        let frame = self.frame(time, cancelled)?;
        let world = frame.world_point(point)?;
        checkpoint(cancelled)?;
        let air = (self.air)(world, time)?;
        checkpoint(cancelled)?;
        if !air.into_iter().all(f64::is_finite) {
            return Err(invalid("panel wake world air velocity"));
        }
        frame.relative_velocity(point, air)
    }

    pub(super) fn accepted_frame(&self) -> Option<PanelFrame2d> {
        self.cached.map(|(_, frame)| frame.state)
    }
}

impl WakeSim {
    /// Prepare the existing panel/Kelvin wake in a prescribed moving chord frame.
    ///
    /// The foil is normalized exactly as in [`Self::new_panel_coupled`]; its
    /// geometry is then placed by `initial_frame`. Reference chord/speed remain
    /// fixed. The nominal steady diagnostic uses unit chordwise flow; actual
    /// flow comes entirely from the total world-air callback during motion steps.
    /// Existing fixed-frame stepping APIs refuse this mode.
    pub fn new_panel_coupled_with_motion(
        foil: &Airfoil2d,
        dt: f64,
        core: f64,
        initial_frame: PanelFrame2d,
    ) -> Result<Self, BemError> {
        PreparedFrame::new(initial_frame)?;
        let (_, chord, _) = foil.chord_frame();
        let angle = fs_math::det::atan2(chord[1], chord[0]);
        let mut wake = Self::new_panel_coupled(foil, angle, dt, core)?;
        let coupling = wake.panel_coupling.as_mut().expect("native constructor");
        // The supplied world pose replaces the input geometry's orientation;
        // retain an exact unit-chordwise reference even after normalization.
        wake.alpha = 0.0;
        wake.steady_gamma = coupling.steady_basis[0];
        coupling.frame = Some(initial_frame);
        Ok(wake)
    }

    /// The last accepted body pose/twist, or `None` in a fixed-frame wake.
    #[must_use]
    pub fn panel_frame(&self) -> Option<PanelFrame2d> {
        self.panel_coupling.as_ref()?.frame
    }

    /// Export the accepted vortices in the prescribed fixed world frame.
    ///
    /// Positions remain in reference-chord units; circulation retains its
    /// clockwise-positive `U_ref*c` normalization. This does not advance state.
    pub fn wake_in_world(&self) -> Result<Vec<WakeVortex>, BemError> {
        let frame = self.panel_frame().ok_or(BemError::InvalidScalar {
            name: "panel wake frame mode",
            value: 0.0,
            requirement: "construct with WakeSim::new_panel_coupled_with_motion",
        })?;
        let prepared = PreparedFrame::new(frame)?;
        let mut result = reserved(self.wake.len(), "world panel wake export")?;
        for vortex in &self.wake {
            result.push(WakeVortex {
                pos: prepared.world_point(vortex.pos)?,
                gamma: vortex.gamma,
            });
        }
        Ok(result)
    }

    /// Advance real prescribed translation/heave/pitch through the native owner.
    ///
    /// `motion(time)` supplies a consistent pose and its derivatives. One frame
    /// is prepared per distinct consecutive stage time, shared by every panel
    /// and wake point. `world_air(point,time)` supplies TOTAL external material
    /// air velocity in the fixed world frame; exclude body/wake induction.
    /// All quantities use [`PanelFrame2d`]'s fixed reference scales.
    ///
    /// At body coordinate r, sample x=X+R*r and use
    /// `R^-1*(U_world-X_dot) - omega*perp(r)` in BOTH the surface boundary solve
    /// and wake transport. The rotation term is coordinate motion, not external
    /// physical vorticity. Physical panel/wake circulation and the original
    /// Kutta/Kelvin closure are unchanged. Wake coordinates, panel velocities,
    /// peak speed and displacement guards are relative to the moving chord
    /// frame; use [`Self::wake_in_world`] for accepted world trajectories.
    ///
    /// The start origin and unwrapped angle must exactly continue the accepted
    /// pose. Callback samples must replay consistently during retries. Smooth
    /// motion/flow supports midpoint convergence; split discontinuities at step
    /// boundaries. The prescribed actuator's dynamics/work and unsteady fluid
    /// pressures/loads are not supplied. This is planar rigid motion, with the
    /// original attached-flow, core and finite packet-release approximations.
    /// Any error/cancellation withholds both the trial wake and endpoint frame.
    pub fn step_panel_coupled_with_motion(
        &mut self,
        options: PanelWakeTransportOptions,
        mut motion: impl FnMut(f64) -> Result<PanelFrame2d, BemError>,
        mut world_air: impl FnMut([f64; 2], f64) -> Result<[f64; 2], BemError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PanelWakeTransportReport, BemError> {
        self.advance_panel_transport(
            self.alpha,
            options,
            &mut IncidentField::moving(&mut motion, &mut world_air),
            &mut cancelled,
        )
    }

    /// Adapt the same prescribed moving-body step without changing shedding.
    ///
    /// Frame, source, scale and callback semantics are identical to
    /// [`Self::step_panel_coupled_with_motion`]. Rejected trials return to the
    /// same accepted body-coordinate wake and absolute forcing times. The
    /// endpoint-position error is measured in the common endpoint body frame,
    /// so rigid world placement cannot change its Euclidean norm. The relative
    /// tolerance and displacement guard retain their body-frame meaning.
    pub fn step_panel_coupled_adaptive_with_motion(
        &mut self,
        options: PanelWakeAdaptiveOptions,
        mut motion: impl FnMut(f64) -> Result<PanelFrame2d, BemError>,
        mut world_air: impl FnMut([f64; 2], f64) -> Result<[f64; 2], BemError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PanelWakeAdaptiveReport, BemError> {
        self.advance_panel_adaptive(
            self.alpha,
            options,
            &mut IncidentField::moving(&mut motion, &mut world_air),
            &mut cancelled,
        )
    }
}
