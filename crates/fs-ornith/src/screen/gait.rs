//! Physical sectional gait screening through the original moving-body wake.
//!
//! The candidate supplies section thickness and mean incidence. Its legacy
//! circulation-gain genes have no role here. Gait kinematics drive the native
//! boundary/Kutta/Kelvin solve and material wake transport; no shed circulation
//! is rescaled. This is a two-dimensional inviscid attached-flow screen, without
//! a thrust, propulsive-efficiency or unsteady pressure/load prediction.

use crate::param::OrnithCandidate;
use fs_bem::BemError;
use fs_bem::panel2d::naca4_symmetric;
use fs_bem::wake2d::{
    MAX_DIRECT_WAKE_VORTICES, PanelFrame2d, PanelWakeAdaptiveOptions, PanelWakeTransportOptions,
    WakeSim,
};
use fs_math::det;

/// Harmonic heave and pitch in the fixed world frame of a unit-chord section.
///
/// `t` uses `c/U_ref`, and the external air is `[1,0] * U_ref`. The material
/// pivot `[pivot_chord,0]` moves vertically as `h=A_h*sin(2*k*t+phase_h)`.
/// The body chord angle is `-candidate.alpha + A_p*sin(2*k*t+phase_p)`;
/// positive pitch rotates the chord counterclockwise. This supplies exact
/// prescribed pose derivatives, not an actuator dynamics or power model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionGait {
    /// Nonnegative heave amplitude divided by the fixed reference chord.
    pub heave_amplitude_chord: f64,
    /// Nonnegative amplitude of body-chord angular displacement [rad].
    pub pitch_amplitude_rad: f64,
    /// Nonnegative reduced frequency `k = omega*c/(2*U_ref)`.
    pub reduced_frequency: f64,
    /// Heave phase [rad].
    pub heave_phase_rad: f64,
    /// Pitch phase [rad], independent of the heave phase.
    pub pitch_phase_rad: f64,
    /// Material pivot along the chord, from leading edge 0 to trailing edge 1.
    pub pivot_chord: f64,
}

impl SectionGait {
    /// Pose/twist for the original BEM moving-frame callback, in its fixed
    /// reference units. This can also drive that owner directly.
    ///
    /// Zero reduced frequency holds the initial pose; zero amplitudes hold
    /// the candidate's mean incidence. Invalid or unrepresentable data refuse.
    pub fn frame_at(
        self,
        mean_incidence_rad: f64,
        time_nondimensional: f64,
    ) -> Result<PanelFrame2d, BemError> {
        if ![
            self.heave_amplitude_chord,
            self.pitch_amplitude_rad,
            self.reduced_frequency,
            self.heave_phase_rad,
            self.pitch_phase_rad,
            self.pivot_chord,
            mean_incidence_rad,
            time_nondimensional,
        ]
        .into_iter()
        .all(f64::is_finite)
            || self.heave_amplitude_chord < 0.0
            || self.pitch_amplitude_rad < 0.0
            || self.reduced_frequency < 0.0
            || !(0.0..=1.0).contains(&self.pivot_chord)
            || mean_incidence_rad.abs() > std::f64::consts::FRAC_PI_2
        {
            return Err(BemError::InvalidScalar {
                name: "physical sectional gait",
                value: time_nondimensional,
                requirement: "finite gait/phase/time, nonnegative amplitudes/frequency, pivot in [0,1] and mean incidence in [-pi/2,pi/2]",
            });
        }
        let omega = 2.0 * self.reduced_frequency;
        let heave_phase = omega * time_nondimensional + self.heave_phase_rad;
        let pitch_phase = omega * time_nondimensional + self.pitch_phase_rad;
        if ![omega, heave_phase, pitch_phase]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(nonfinite());
        }
        let heave = self.heave_amplitude_chord * det::sin(heave_phase);
        let heave_rate = self.heave_amplitude_chord * omega * det::cos(heave_phase);
        let angle = -mean_incidence_rad + self.pitch_amplitude_rad * det::sin(pitch_phase);
        let rate = self.pitch_amplitude_rad * omega * det::cos(pitch_phase);
        let pivot = [
            self.pivot_chord * det::cos(angle),
            self.pivot_chord * det::sin(angle),
        ];
        let frame = PanelFrame2d {
            origin_world: [self.pivot_chord - pivot[0], heave - pivot[1]],
            angle_rad: angle,
            origin_velocity_world: [rate * pivot[1], heave_rate - rate * pivot[0]],
            angular_velocity: rate,
        };
        if !frame
            .origin_world
            .into_iter()
            .chain(frame.origin_velocity_world)
            .chain([frame.angle_rad, frame.angular_velocity])
            .all(f64::is_finite)
        {
            return Err(nonfinite());
        }
        Ok(frame)
    }
}

/// Original native transport selection; no additional timestep owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GaitTransport {
    Fixed(PanelWakeTransportOptions),
    Adaptive(PanelWakeAdaptiveOptions),
}

/// Finite screening interval and native accuracy/work controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaitScreenControls {
    /// Section panel count, admitted by the original NACA/BEM owner.
    pub panels: usize,
    /// Number of original shedding intervals, in 1..=MAX_DIRECT_WAKE_VORTICES.
    /// This conservative cap also bounds the report's complete retained wake.
    pub intervals: usize,
    /// Fixed shedding interval in `c/U_ref`; refinement does not emit extra rows.
    pub interval_nondimensional: f64,
    /// Vortex core divided by `c`, with the original native admission.
    pub core_radius_chord: f64,
    /// Fixed midpoint refinement or original adaptive position-error control.
    pub transport: GaitTransport,
    /// Whole-screen coupled-solve cap, INCLUDING adaptive rejected trials.
    /// Each interval also retains its selected native per-interval cap.
    pub maximum_coupled_solves: usize,
}

/// Actual accepted native wake and an explicitly kinematic screening metric.
#[derive(Clone, Debug)]
pub struct GaitScreenReport {
    /// Original BEM owner: accepted circulation/history, body-frame wake,
    /// accepted pose and `wake_in_world()` export. No strengths are modified.
    pub wake: WakeSim,
    /// `sum |Gamma_i|*|Delta x_world_i| / sum |Gamma_i|*Delta t`, in `U_ref`.
    /// Each term follows an existing physical packet across one accepted
    /// shedding interval. Newly emitted packets join on the NEXT interval,
    /// because this consumer does not reconstruct the owner's release path.
    /// This measures streamwise travel, including ambient convection; it is
    /// neither a force, signed momentum flux nor a propulsive efficiency.
    /// Zero represented circulation-time gives zero.
    pub circulation_weighted_streamwise_speed: f64,
    /// Denominator of that metric, in `U_ref*c` times `c/U_ref`.
    pub circulation_time_weight: f64,
    /// Largest accepted `|bound + sum physical wake Gamma|`, in `U_ref*c`.
    /// This is a numerical residual, not a validation of the model's accuracy.
    pub maximum_kelvin_residual: f64,
    /// All actual native coupled solves, including rejected adaptive work.
    pub coupled_solves: usize,
    /// Accepted native transport substeps over all shedding intervals.
    pub transport_substeps: usize,
}

fn nonfinite() -> BemError {
    BemError::NonFiniteResult {
        operation: "physical sectional gait",
    }
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), BemError> {
    if cancelled() {
        Err(BemError::Cancelled {
            operation: "sectional gait screen",
        })
    } else {
        Ok(())
    }
}

/// Screen a candidate section under an explicit physical heave/pitch gait.
///
/// `candidate.thickness` and `candidate.alpha` supply the shape and mean
/// incidence. Legacy `flap_amp`/`flap_freq` circulation-gain genes and the inlet
/// proxy do not alter this prescribed gait. The world flow is uniform `[1,0]`
/// in fixed reference-speed units. One retained panel factorization serves
/// every stage, including adaptive retries, and all old circulation is retained.
///
/// Controls bound panels, wake population, intervals and total coupled solves.
/// Cancellation reaches the original native stage checks as well as this
/// consumer's reductions. Any refusal returns no partial report. The caller
/// retains the original candidate/gait; cancellation callback effects are external.
///
/// This exposes a real 2-D moving-section screen. It does not resolve 3-D
/// flapping, dynamic stall, viscosity, actuator work or unsteady pressure forces.
/// Refine shedding intervals, panels and core as well as transport when judging
/// physical predictions; adaptive transport alone does not bound those errors.
pub fn screen_section_gait(
    candidate: &OrnithCandidate,
    gait: SectionGait,
    controls: GaitScreenControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GaitScreenReport, BemError> {
    checkpoint(&mut cancelled)?;
    if controls.intervals == 0 || controls.intervals > MAX_DIRECT_WAKE_VORTICES {
        return Err(BemError::InvalidScalar {
            name: "sectional gait intervals",
            value: controls.intervals as f64,
            requirement: "positive and no greater than the native direct-wake ceiling",
        });
    }
    let frame = gait.frame_at(candidate.alpha, 0.0)?;
    let foil = naca4_symmetric(candidate.thickness, controls.panels)?;
    let mut wake = WakeSim::new_panel_coupled_with_motion(
        &foil,
        controls.interval_nondimensional,
        controls.core_radius_chord,
        frame,
    )?;
    let mut previous = wake.wake_in_world()?;
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut kelvin: f64 = 0.0;
    let mut coupled_solves = 0;
    let mut transport_substeps = 0;
    for _ in 0..controls.intervals {
        checkpoint(&mut cancelled)?;
        let remaining = controls.maximum_coupled_solves - coupled_solves;
        let motion = |time| gait.frame_at(candidate.alpha, time);
        let (solves, substeps) = match controls.transport {
            GaitTransport::Fixed(mut options) => {
                options.max_coupled_solves = options.max_coupled_solves.min(remaining);
                let report = wake.step_panel_coupled_with_motion(
                    options,
                    motion,
                    |_, _| Ok([1.0, 0.0]),
                    &mut cancelled,
                )?;
                (report.coupled_solves, report.substeps)
            }
            GaitTransport::Adaptive(mut options) => {
                options.max_coupled_solves = options.max_coupled_solves.min(remaining);
                let report = wake.step_panel_coupled_adaptive_with_motion(
                    options,
                    motion,
                    |_, _| Ok([1.0, 0.0]),
                    &mut cancelled,
                )?;
                (report.coupled_solves, report.substeps)
            }
        };
        coupled_solves += solves;
        transport_substeps += substeps;
        let world = wake.wake_in_world()?;
        let mut circulation = wake
            .panel_state()
            .expect("accepted native interval")
            .bound_circulation;
        for vortex in &world {
            checkpoint(&mut cancelled)?;
            circulation += vortex.gamma;
        }
        kelvin = kelvin.max(circulation.abs());
        for (old, new) in previous.iter().zip(&world) {
            checkpoint(&mut cancelled)?;
            numerator += old.gamma.abs() * (new.pos[0] - old.pos[0]).abs();
            denominator += old.gamma.abs() * controls.interval_nondimensional;
        }
        if ![numerator, denominator, circulation]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(nonfinite());
        }
        previous = world;
    }
    let speed = if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    };
    if !speed.is_finite() {
        return Err(nonfinite());
    }
    checkpoint(&mut cancelled)?;
    Ok(GaitScreenReport {
        wake,
        circulation_weighted_streamwise_speed: speed,
        circulation_time_weight: denominator,
        maximum_kelvin_residual: kelvin,
        coupled_solves,
        transport_substeps,
    })
}
