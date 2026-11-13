//! Inertial chamber motion with the original finite IF97 water inventory.
//!
//! Total liquid-plus-vapor volume is the ONLY piston-position coordinate.
//! Each trial asks `EquilibriumWaterCell::propose_volume` for its original
//! conserved M,U, phase inventory and reciprocal pressure work. The existing
//! gas-piston mechanical kernel uses that SAME mean pressure in its impulse:
//! dV=A dt (v0+v1)/2, M_wall (v1-v0)=dt [A (p_mean-p_external)+F].
//! Latent calorics occur once in the water source; vapor volume is never
//! substituted for total chamber volume or interpreted as a spherical bubble.
//!
//! This is a sealed, constant-area, uniform-pressure chamber with a rigid
//! translating wall. Uniform temperature and rapid liquid/vapor equilibration
//! are explicit caller conditions. It does not resolve water-hammer waves,
//! nucleation delay, moving phase fronts, noncondensable gas or wall impacts.
//! Heat is an integrated external energy port; pressure and additional force
//! are held over each interval. Split known load changes at their actual times.
//!
//! Water work uses its original trapezoidal endpoint pressure. Wall momentum
//! uses native midpoint kinematics, so their discrete energy exchange is
//! reciprocal. Time, volume and temperature limits are admission bounds, not
//! an error estimator; refine the actual interval to resolve nonlinear motion
//! and phase crossings. Both conserved and original-source total-energy
//! residuals are checked independently before either history is published.

use super::volume::{WaterVolumeControls, WaterVolumeError, WaterVolumeStep};
use super::{EquilibriumWaterCell, WaterCellError};
pub use crate::cavity_phs::piston::mechanics::PistonSupport;
use crate::cavity_phs::piston::mechanics::{self, MechanicalTrial};
use crate::cavity_phs::piston::thermal::feedback;

/// External ports over one actual physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPistonInput {
    /// Integrated heat INTO the finite water inventory [J].
    pub heat_j: f64,
    /// Held absolute pressure on the external wall face [Pa], zero for vacuum.
    pub exterior_pressure_pa: f64,
    /// Held additional force [N], positive towards increasing chamber volume.
    /// Water and exterior-pressure loads are already accounted separately.
    pub external_force_n: f64,
}

/// Native material-work controls and independent mechanical admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPistonControls {
    pub water: WaterVolumeControls,
    /// Positive endpoint-volume bracket containing the accepted TOTAL volume.
    /// Every queried state must fit the original source/volume/temperature
    /// bounds. Fixed support may use [V,V] and queries only the original V.
    pub volume_bracket_m3: [f64; 2],
    pub momentum_tolerance_n_s: f64,
    /// Applies separately to mechanical, conserved-total and source-total work.
    pub energy_tolerance_j: f64,
    /// Bounded native scalar-solve iterations, in 1..=128. At most this count
    /// plus two complete, individually bounded water proposals are evaluated.
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum WaterPistonError {
    Invalid(&'static str),
    Water(WaterVolumeError),
    OutsideBracket {
        residual_n_s: [f64; 2],
    },
    NotConverged {
        iterations: usize,
        residual_n_s: f64,
    },
    Budget {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl From<WaterVolumeError> for WaterPistonError {
    fn from(error: WaterVolumeError) -> Self {
        if error == WaterVolumeError::Water(WaterCellError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Water(error)
        }
    }
}
impl core::fmt::Display for WaterPistonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite water/piston coupling: {self:?}")
    }
}
impl std::error::Error for WaterPistonError {}

/// One original water cell and one native rigid-wall momentum coordinate.
/// The water's accepted clock and total volume remain the sole time/position
/// owners; kinetic energy is derived from wall mass and velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPiston {
    water: EquilibriumWaterCell,
    wall_mass_kg: f64,
    area_m2: f64,
    velocity_m_s: f64,
    support: PistonSupport,
}

/// Actual material, momentum, phase and reciprocal-work accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPistonStep {
    pub water: WaterVolumeStep,
    pub displacement_m: f64,
    pub velocity_m_s: [f64; 2],
    pub kinetic_energy_change_j: f64,
    pub exterior_pressure_work_j: f64,
    pub external_force_work_j: f64,
    pub support_impulse_n_s: f64,
    pub support_work_j: f64,
    /// Delta V minus A dt times actual midpoint velocity [m3].
    pub kinematic_residual_m3: f64,
    pub momentum_residual_n_s: f64,
    /// Delta K minus original water work plus external pressure work minus
    /// additional-force and support work [J].
    pub mechanical_energy_residual_j: f64,
    /// Conserved water Delta U plus Delta K and external pressure work,
    /// minus heat, additional-force and support work [J].
    pub total_energy_residual_j: f64,
    /// Same total account using the original recovered IF97 energy change.
    pub source_total_energy_residual_j: f64,
    /// Positive evaporation, negative condensation, from original phase masses.
    pub vapor_mass_change_kg: f64,
    /// Actual original vapor-volume change, distinct from chamber displacement.
    pub vapor_volume_change_m3: f64,
    pub mechanical_trials: usize,
}

/// Inspect a complete mechanical/material proposal before atomic publication.
/// Dropping it leaves the original water, wall velocity and clock unchanged.
#[must_use = "inspect and commit to publish the complete water/piston interval"]
pub struct WaterPistonCandidate<'a> {
    owner: &'a mut WaterPiston,
    endpoint: WaterPiston,
    step: WaterPistonStep,
}
impl WaterPistonCandidate<'_> {
    pub const fn endpoint(&self) -> &WaterPiston {
        &self.endpoint
    }
    pub const fn step(&self) -> &WaterPistonStep {
        &self.step
    }
    pub fn commit(self) -> WaterPistonStep {
        *self.owner = self.endpoint;
        self.step
    }
}

struct Trial {
    water: EquilibriumWaterCell,
    step: WaterVolumeStep,
    motion: MechanicalTrial,
    trials: usize,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterPistonError> {
    if cancelled() {
        Err(WaterPistonError::Cancelled)
    } else {
        Ok(())
    }
}
fn budget(quantity: &'static str, residual: f64, tolerance: f64) -> Result<(), WaterPistonError> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(WaterPistonError::Budget {
            quantity,
            residual,
            tolerance,
        })
    } else {
        Ok(())
    }
}

impl WaterPiston {
    pub fn new(
        water: EquilibriumWaterCell,
        wall_mass_kg: f64,
        area_m2: f64,
        velocity_m_s: f64,
    ) -> Result<Self, WaterPistonError> {
        if !positive(wall_mass_kg)
            || !positive(area_m2)
            || !velocity_m_s.is_finite()
            || !positive(water.state().volume_m3() / area_m2)
            || !(0.5 * wall_mass_kg * velocity_m_s * velocity_m_s).is_finite()
        {
            return Err(WaterPistonError::Invalid(
                "positive represented wall mass/area/position and finite velocity/kinetic energy required",
            ));
        }
        Ok(Self {
            water,
            wall_mass_kg,
            area_m2,
            velocity_m_s,
            support: PistonSupport::Free,
        })
    }

    pub const fn water(&self) -> &EquilibriumWaterCell {
        &self.water
    }
    pub const fn wall_mass_kg(&self) -> f64 {
        self.wall_mass_kg
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn velocity_m_s(&self) -> f64 {
        self.velocity_m_s
    }
    pub fn position_m(&self) -> f64 {
        self.water.state().volume_m3() / self.area_m2
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        0.5 * self.wall_mass_kg * self.velocity_m_s * self.velocity_m_s
    }
    pub const fn support(&self) -> PistonSupport {
        self.support
    }

    /// Attach an ideal support to an already stationary wall. Stopping motion
    /// requires a separate accounted impulse; it cannot silently delete energy.
    pub fn with_fixed_support(mut self) -> Result<Self, WaterPistonError> {
        if self.velocity_m_s != 0.0 {
            return Err(WaterPistonError::Invalid(
                "fixed support requires a wall already at rest",
            ));
        }
        self.support = PistonSupport::Fixed;
        Ok(self)
    }

    /// Solve and publish water pressure/phase and wall motion together.
    pub fn step(
        &mut self,
        end_time_s: f64,
        input: WaterPistonInput,
        controls: WaterPistonControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<WaterPistonStep, WaterPistonError> {
        Ok(self
            .propose_step(end_time_s, input, controls, cancelled)?
            .commit())
    }

    /// Prepare one complete physical interval. All source evaluations start
    /// from the SAME accepted inventory. Rejected trials never advance time,
    /// alter phase mass, or consume heat. Source, work, impulse or cancellation
    /// failures preserve the full original owner. Callback side effects are
    /// outside this transaction; cancellation must remain bounded and pure.
    pub fn propose_step(
        &mut self,
        end_time_s: f64,
        input: WaterPistonInput,
        controls: WaterPistonControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterPistonCandidate<'_>, WaterPistonError> {
        poll(&mut cancelled)?;
        let volume = self.water.state().volume_m3();
        let duration = end_time_s - self.water.state().time_s();
        let [lo, hi] = controls.volume_bracket_m3;
        if !end_time_s.is_finite()
            || !positive(duration)
            || !input.heat_j.is_finite()
            || !input.external_force_n.is_finite()
            || !input.exterior_pressure_pa.is_finite()
            || input.exterior_pressure_pa < 0.0
            || !positive(lo)
            || !positive(hi)
            || lo > volume
            || hi < volume
            || (lo == hi && self.support == PistonSupport::Free)
            || !positive(controls.momentum_tolerance_n_s)
            || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
        {
            return Err(WaterPistonError::Invalid(
                "forward physical clock, finite loads, admitted total-volume bracket and positive bounded work controls required",
            ));
        }
        let trial = if self.support == PistonSupport::Fixed {
            self.trial(volume, end_time_s, input, controls.water, 1, &mut cancelled)?
        } else {
            feedback::solve_scalar(
                [lo - volume, hi - volume],
                controls.momentum_tolerance_n_s,
                controls.maximum_iterations,
                |increment, trials, cancelled| {
                    let next = volume + increment;
                    if increment != 0.0 && next == volume {
                        return Err(WaterPistonError::Invalid(
                            "piston volume increment is unrepresentable",
                        ));
                    }
                    let trial =
                        self.trial(next, end_time_s, input, controls.water, trials, cancelled)?;
                    let residual = trial.motion.residual;
                    Ok((trial, residual))
                },
                &mut cancelled,
            )
            .map_err(|error| match error {
                feedback::ScalarError::Evaluation(error) => error,
                feedback::ScalarError::Invalid(message) => WaterPistonError::Invalid(message),
                feedback::ScalarError::Cancelled => WaterPistonError::Cancelled,
                feedback::ScalarError::OutsideBracket { lower, upper } => {
                    WaterPistonError::OutsideBracket {
                        residual_n_s: [lower, upper],
                    }
                }
                feedback::ScalarError::NotConverged {
                    iterations,
                    residual,
                } => WaterPistonError::NotConverged {
                    iterations,
                    residual_n_s: residual,
                },
            })?
        };
        let step = self.finish(&trial, input, duration, controls)?;
        let mut endpoint = *self;
        endpoint.water = trial.water;
        endpoint.velocity_m_s = trial.motion.velocity;
        poll(&mut cancelled)?;
        Ok(WaterPistonCandidate {
            owner: self,
            endpoint,
            step,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        volume_m3: f64,
        end_time_s: f64,
        input: WaterPistonInput,
        controls: WaterVolumeControls,
        trials: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, WaterPistonError> {
        poll(cancelled)?;
        // Copying this small source-bound owner retains its original controls
        // and M,U; only the tentative copy can be changed by a rejected solve.
        let mut water = self.water;
        let candidate = water.propose_volume(
            volume_m3,
            input.heat_j,
            end_time_s,
            controls,
            &mut *cancelled,
        )?;
        let motion = mechanics::trial(
            self.wall_mass_kg,
            self.area_m2,
            self.water.state().volume_m3(),
            self.velocity_m_s,
            self.support,
            volume_m3,
            candidate.step().mean_boundary_pressure_pa,
            end_time_s - self.water.state().time_s(),
            input.exterior_pressure_pa,
            input.external_force_n,
        )
        .map_err(WaterPistonError::Invalid)?;
        let step = candidate.commit();
        poll(cancelled)?;
        Ok(Trial {
            water,
            step,
            motion,
            trials,
        })
    }

    fn finish(
        &self,
        trial: &Trial,
        input: WaterPistonInput,
        duration: f64,
        controls: WaterPistonControls,
    ) -> Result<WaterPistonStep, WaterPistonError> {
        let before = self.velocity_m_s;
        let after = trial.motion.velocity;
        let kinetic = (0.5 * self.wall_mass_kg * (after - before)) * (after + before);
        let exterior = input.exterior_pressure_pa * trial.step.volume_change_m3;
        let external = input.external_force_n * trial.motion.displacement;
        let support = trial.motion.support_impulse * f64::midpoint(before, after);
        let kinematic =
            trial.step.volume_change_m3 - self.area_m2 * duration * f64::midpoint(before, after);
        let mechanical = kinetic - trial.step.water.boundary_work_j + exterior - external - support;
        let total = trial.step.water.internal_energy_change_j + kinetic + exterior
            - external
            - support
            - input.heat_j;
        let source_total = trial.step.water.source_first_law_residual_j + mechanical;
        budget(
            "wall kinematics [m3]",
            kinematic,
            self.water.controls().volume_tolerance_m3,
        )?;
        budget(
            "wall impulse [N s]",
            trial.motion.residual,
            controls.momentum_tolerance_n_s,
        )?;
        budget(
            "mechanical energy [J]",
            mechanical,
            controls.energy_tolerance_j,
        )?;
        budget(
            "conserved total energy [J]",
            total,
            controls.energy_tolerance_j,
        )?;
        budget(
            "original-source total energy [J]",
            source_total,
            controls.energy_tolerance_j,
        )?;
        if ![kinetic, exterior, external, support]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(WaterPistonError::Invalid("unrepresentable piston work"));
        }
        let vapor_mass = trial.step.water.state.phase_masses_kg()[1]
            - trial.step.water.initial.phase_masses_kg()[1];
        let vapor_volume = trial.step.water.state.phase_volumes_m3()[1]
            - trial.step.water.initial.phase_volumes_m3()[1];
        Ok(WaterPistonStep {
            water: trial.step,
            displacement_m: trial.motion.displacement,
            velocity_m_s: [before, after],
            kinetic_energy_change_j: kinetic,
            exterior_pressure_work_j: exterior,
            external_force_work_j: external,
            support_impulse_n_s: trial.motion.support_impulse,
            support_work_j: support,
            kinematic_residual_m3: kinematic,
            momentum_residual_n_s: trial.motion.residual,
            mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total,
            source_total_energy_residual_j: source_total,
            vapor_mass_change_kg: vapor_mass,
            vapor_volume_change_m3: vapor_volume,
            mechanical_trials: trial.trials,
        })
    }
}

#[cfg(test)]
#[path = "piston/tests.rs"]
mod tests;
