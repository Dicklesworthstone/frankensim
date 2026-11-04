//! Exact linear time observations of the original OU inlet on native CFD steps.

use super::*;
use crate::nsf3::{NsfBoundaryDrive, NsfCharacteristicTrial};
use fs_atmo::mann_field::fixed_step::MannStepAverage;

/// Whole-interval work ceilings, including catch-up and the final profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NsfEvolvingMannAverageBudget {
    pub maximum_cell_steps: usize,
    /// Four normals per mode for catch-up steps; eight per averaged step
    /// (the unchanged endpoint innovation and an independent integral bridge).
    pub maximum_normal_draws: usize,
    /// Tiles times modes times `(steps + 1)`, including the final endpoint.
    pub maximum_mode_samples: usize,
}

/// Provenance and work for linear velocity averages over native steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfEvolvingMannAverageReport {
    pub steps: usize,
    pub profile_updates: usize,
    pub mode_samples: usize,
    pub normal_draws: usize,
    /// Source endpoints of the averaged interval, excluding prior catch-up.
    pub atmosphere_step_range: [u64; 2],
    pub atmosphere_time_range_s: [f64; 2],
    pub cfd_time_range_s: [f64; 2],
    /// Instantaneous profile restored at the final accepted fluid/source clock.
    pub endpoint: NsfEvolvingMannApplied,
}

/// One unpublished native fluid interval and its original OU endpoint owner.
/// Dropping the trial leaves both accepted owners untouched. The supplied
/// `MannStepAverage` is observation scratch, not a third accepted participant.
pub struct NsfEvolvingMannAverageTrial<'a, L: Lattice3<27>> {
    inlet: &'a mut NsfEvolvingMannInflow,
    fluid: NsfCharacteristicTrial<'a, L>,
    report: NsfEvolvingMannAverageReport,
}

impl<'a, L: Lattice3<27>> NsfEvolvingMannAverageTrial<'a, L> {
    pub const fn fluid(&self) -> &NsfCharacteristicTrial<'a, L> {
        &self.fluid
    }

    pub const fn field(&self) -> &FixedStepMannField {
        &self.inlet.trial
    }

    pub const fn report(&self) -> NsfEvolvingMannAverageReport {
        self.report
    }

    /// Publish the already admitted native fluid and atmosphere endpoints.
    pub fn accept(self) {
        self.fluid.accept();
        core::mem::swap(&mut self.inlet.field, &mut self.inlet.trial);
    }
}

fn atmosphere_error(error: fs_atmo::Refusal) -> NsfEvolvingMannError {
    if error.code == "mann-fixed-step-cancelled" {
        NsfEvolvingMannError::Cancelled
    } else {
        NsfEvolvingMannError::Atmosphere(error)
    }
}

impl NsfEvolvingMannInflow {
    /// Use the exact Gaussian time-average of the advected OU velocity during
    /// every native CFD step, jointly sampled with its unchanged endpoint.
    /// The authored point/face/cell observation and optional slip plane are
    /// applied to that same linear average; sampling uses no additional RNG.
    ///
    /// This reduces temporal undersampling of the represented velocity without
    /// increasing the fluid step count. It is a linear velocity observation:
    /// the native reservoir still holds that velocity during its step, so its
    /// nonlinear momentum and energy fluxes are not exact time-averaged fluxes.
    /// Existing spatial and advective Nyquist admissions remain in force.
    /// No subgrid kinetic energy or acoustic source is invented.
    ///
    /// Prepare `average` once from `self.field()` with its explicit byte budget.
    /// This method reuses it, the inlet's trial source/profile, and the original
    /// single fluid candidate. The last average can remain candidate scratch
    /// after refusal or a dropped trial. Accepted fluid/source states remain
    /// unchanged, and retry reproduces the same endpoint and bridge draws.
    /// A successful interval restores the instantaneous final profile, making
    /// ordinary same-time `apply` idempotent and preserving the checkpoint law.
    pub fn preview_time_averaged_boundary_interval<'a, L: Lattice3<27>>(
        &'a mut self,
        grid: &'a mut NsfGrid3<L>,
        average: &mut MannStepAverage,
        drives: [NsfBoundaryDrive; 2],
        steps: usize,
        budget: NsfEvolvingMannAverageBudget,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfEvolvingMannAverageTrial<'a, L>, NsfEvolvingMannError> {
        if cancelled() {
            return Err(NsfEvolvingMannError::Cancelled);
        }
        if steps == 0 {
            return Err(NsfEvolvingMannError::Input(
                "time-averaged inlet interval requires positive native steps",
            ));
        }
        average
            .validate_source(&self.field)
            .map_err(atmosphere_error)?;
        let first_cfd_step = grid.steps();
        let first_target = self.target_step(grid)?;
        let catch_up = self.transition_work(self.field.step(), first_target)?;
        let final_cfd_step =
            first_cfd_step
                .checked_add(steps)
                .ok_or(NsfEvolvingMannError::Input(
                    "time-averaged inlet CFD clock overflows",
                ))?;
        if final_cfd_step as u128 >= (1u128 << 53) {
            return Err(NsfEvolvingMannError::Input(
                "CFD step clock exceeds exact f64 integer range",
            ));
        }
        let final_target = self.atmosphere_step_at(final_cfd_step)?;
        let draws_per_step =
            self.field
                .model()
                .mode_count()
                .checked_mul(8)
                .ok_or(NsfEvolvingMannError::Input(
                    "time-averaged inlet draw work overflows",
                ))?;
        let normal_draws = steps
            .checked_mul(draws_per_step)
            .and_then(|draws| draws.checked_add(catch_up))
            .ok_or(NsfEvolvingMannError::Input(
                "time-averaged inlet draw work overflows",
            ))?;
        let profiles = steps.checked_add(1).ok_or(NsfEvolvingMannError::Input(
            "time-averaged inlet sample work overflows",
        ))?;
        let mode_samples = profiles
            .checked_mul(self.profile.len())
            .and_then(|work| work.checked_mul(self.field.model().mode_count()))
            .ok_or(NsfEvolvingMannError::Input(
                "time-averaged inlet sample work overflows",
            ))?;
        for (requested, maximum) in [
            (normal_draws, budget.maximum_normal_draws),
            (mode_samples, budget.maximum_mode_samples),
        ] {
            if requested > maximum {
                return Err(NsfEvolvingMannError::Work { requested, maximum });
            }
        }
        self.trial.copy_temporal_state_from(&self.field)?;
        let mut first_atmosphere_time = None;
        let mut updates = 0;
        let mut endpoint = None;
        let fluid = grid
            .preview_boundary_interval_with_updates(
                drives,
                [None; 2],
                &[],
                steps,
                budget.maximum_cell_steps,
                &cancelled,
                |candidate| {
                    // Catch up only after the native interval has admitted its
                    // cell work and eligibility. Refused native budgets must
                    // not consume a potentially long atmosphere prehistory.
                    if updates == 0 {
                        self.trial
                            .advance_to(first_target, catch_up, &cancelled)
                            .map_err(atmosphere_error)?;
                        first_atmosphere_time = Some(self.trial.time_s());
                    }
                    if candidate.steps() == final_cfd_step {
                        endpoint = Some(self.stage_profile(candidate, 0, &cancelled)?);
                    } else {
                        // Each mean covers exactly the step about to be taken,
                        // while the phase/source endpoint advances to its end.
                        if self.target_step(candidate)? != self.trial.step() {
                            return Err(NsfEvolvingMannError::Input(
                                "time-averaged inlet source and fluid clocks differ",
                            ));
                        }
                        self.trial
                            .advance_one_with_average(average, draws_per_step, &cancelled)
                            .map_err(atmosphere_error)?;
                        sample_profile(
                            candidate,
                            self.binding.aperture,
                            SamplingPlan {
                                config: self.config.frame(),
                                maximum_wavenumber_rad_m: self.maximum_wavenumber_rad_m,
                                convection_speed_m_s: self.advection_speed_m_s,
                                mode_count: self.trial.model().mode_count(),
                                averaging: self.averaging.as_ref(),
                            },
                            &mut self.profile,
                            |point, filter| match filter {
                                Some(SpatialFilter::Wall(filter)) => {
                                    average.sample_wall_blocked_filtered(point, filter)
                                }
                                Some(SpatialFilter::Homogeneous(filter)) => {
                                    average.sample_filtered(point, filter)
                                }
                                None => match self.wall {
                                    Some(wall) => average.sample_wall_blocked(point, wall),
                                    None => average.sample(point),
                                },
                            },
                            &mut &cancelled,
                        )?;
                        candidate
                            .set_reservoir_aperture_profile(self.binding.aperture, &self.profile)?;
                    }
                    updates += 1;
                    Ok(())
                },
            )
            .map_err(|error| match error {
                NsfEvolvingMannError::Boundary(NsfError::Input(
                    "characteristic interval cancelled",
                )) => NsfEvolvingMannError::Cancelled,
                other => other,
            })?;
        debug_assert_eq!(updates, profiles);
        debug_assert_eq!(self.trial.step(), final_target);
        let report = NsfEvolvingMannAverageReport {
            steps,
            profile_updates: updates,
            mode_samples,
            normal_draws,
            atmosphere_step_range: [first_target, final_target],
            atmosphere_time_range_s: [
                first_atmosphere_time.expect("positive interval has a first native step"),
                self.trial.time_s(),
            ],
            cfd_time_range_s: [
                first_cfd_step as f64 * self.binding.dt_s,
                final_cfd_step as f64 * self.binding.dt_s,
            ],
            endpoint: endpoint.expect("native interval refreshes its final endpoint"),
        };
        Ok(NsfEvolvingMannAverageTrial {
            inlet: self,
            fluid,
            report,
        })
    }
}

#[cfg(test)]
mod tests;
