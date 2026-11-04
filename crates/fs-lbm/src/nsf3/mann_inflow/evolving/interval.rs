//! Native CFD intervals with a distinct evolving inlet at every step boundary.

use super::*;
use crate::nsf3::{NsfBoundaryDrive, NsfCharacteristicTrial};

/// Work ceilings for the complete interval, not allowances per substep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NsfEvolvingMannIntervalBudget {
    pub maximum_cell_steps: usize,
    /// Four innovations per represented mode per elapsed atmosphere step,
    /// including any catch-up before the first fluid step and the final endpoint.
    pub maximum_ou_mode_transitions: usize,
    /// Tiles times modes times `(steps + 1)`, including the final profile.
    pub maximum_mode_samples: usize,
}

/// Both endpoint clocks and the work actually consumed by a complete interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfEvolvingMannIntervalReport {
    pub steps: usize,
    pub profile_updates: usize,
    pub mode_samples: usize,
    pub ou_mode_transitions: usize,
    /// Source/profile used for the first original native step.
    pub first: NsfEvolvingMannApplied,
    /// Source/profile synchronized to the final CFD clock, after the last step.
    pub endpoint: NsfEvolvingMannApplied,
}

/// Unpublished fluid and atmosphere interval. Dropping it rolls back both
/// accepted owners; the inlet retains only reusable, unpublished scratch state.
/// The native fluid trial carries all original mass, energy and work accounts.
pub struct NsfEvolvingMannTrial<'a, L: Lattice3<27>> {
    inlet: &'a mut NsfEvolvingMannInflow,
    fluid: NsfCharacteristicTrial<'a, L>,
    report: NsfEvolvingMannIntervalReport,
}

impl<'a, L: Lattice3<27>> NsfEvolvingMannTrial<'a, L> {
    /// Original candidate grid and original chronological interval accounts.
    pub const fn fluid(&self) -> &NsfCharacteristicTrial<'a, L> {
        &self.fluid
    }

    /// Candidate atmosphere at the final CFD endpoint, for inspection/restart.
    pub const fn field(&self) -> &FixedStepMannField {
        &self.inlet.trial
    }

    pub const fn report(&self) -> NsfEvolvingMannIntervalReport {
        self.report
    }

    /// Publish both already admitted owners with no further fallible work.
    pub fn accept(self) {
        self.fluid.accept();
        core::mem::swap(&mut self.inlet.field, &mut self.inlet.trial);
    }
}

impl NsfEvolvingMannInflow {
    /// Evolve turbulence at every native CFD step boundary while previewing
    /// the original characteristic interval. The mean, frame, gas source,
    /// spatial filter and optional plane are the inlet's existing selections.
    /// Other reservoir and characteristic-patch controls remain held.
    ///
    /// At CFD time `n*dt`, the corresponding atmosphere supplies the profile
    /// held through step `n -> n+1`. The final endpoint is also synchronized,
    /// so accepting the trial publishes matching clocks and a ready next-step
    /// profile. There are `steps + 1` profiles, but only `steps` fluid transfers.
    /// This preserves the native held-control time discretization within each
    /// step; it does not interpolate OU noise inside predictor/corrector passes.
    ///
    /// The existing native interval owns the single fluid candidate and all
    /// original mass/energy aggregation. Source evolution uses the inlet's
    /// retained trial buffer. Refusal, cancellation or dropping the result
    /// changes neither accepted participant; retry draws the same innovations.
    pub fn preview_boundary_interval<'a, L: Lattice3<27>>(
        &'a mut self,
        grid: &'a mut NsfGrid3<L>,
        drives: [NsfBoundaryDrive; 2],
        steps: usize,
        budget: NsfEvolvingMannIntervalBudget,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfEvolvingMannTrial<'a, L>, NsfEvolvingMannError> {
        if cancelled() {
            return Err(NsfEvolvingMannError::Cancelled);
        }
        if steps == 0 {
            return Err(NsfEvolvingMannError::Input(
                "evolving inlet interval requires positive native steps",
            ));
        }
        // Check the bound geometry and both ends of the relative clock before
        // allocating the native fluid candidate or consuming any innovations.
        let first_target = self.target_step(grid)?;
        self.transition_work(self.field.step(), first_target)?;
        let final_cfd_step = grid
            .steps()
            .checked_add(steps)
            .ok_or(NsfEvolvingMannError::Input(
                "evolving inlet interval CFD clock overflows",
            ))?;
        if final_cfd_step as u128 >= (1u128 << 53) {
            return Err(NsfEvolvingMannError::Input(
                "CFD step clock exceeds exact f64 integer range",
            ));
        }
        let final_target = self.atmosphere_step_at(final_cfd_step)?;
        let transitions = self.transition_work(self.field.step(), final_target)?;
        let profiles = steps.checked_add(1).ok_or(NsfEvolvingMannError::Input(
            "evolving inlet interval sample work overflows",
        ))?;
        let mode_samples = profiles
            .checked_mul(self.profile.len())
            .and_then(|work| work.checked_mul(self.field.model().mode_count()))
            .ok_or(NsfEvolvingMannError::Input(
                "evolving inlet interval sample work overflows",
            ))?;
        for (requested, maximum) in [
            (transitions, budget.maximum_ou_mode_transitions),
            (mode_samples, budget.maximum_mode_samples),
        ] {
            if requested > maximum {
                return Err(NsfEvolvingMannError::Work { requested, maximum });
            }
        }
        self.trial.copy_temporal_state_from(&self.field)?;
        let mut first = None;
        let mut endpoint = None;
        let mut used_transitions = 0;
        let mut used_samples = 0;
        let mut profile_updates = 0;
        let fluid = grid
            .preview_boundary_interval_with_updates(
                drives,
                [None; 2],
                &[],
                steps,
                budget.maximum_cell_steps,
                &cancelled,
                |candidate| {
                    let applied = self.stage_profile(
                        candidate,
                        budget.maximum_ou_mode_transitions - used_transitions,
                        &cancelled,
                    )?;
                    used_transitions += applied.ou_mode_transitions;
                    used_samples += applied.mode_samples;
                    profile_updates += 1;
                    first.get_or_insert(applied);
                    endpoint = Some(applied);
                    Ok(())
                },
            )
            .map_err(|error| match error {
                NsfEvolvingMannError::Boundary(NsfError::Input(
                    "characteristic interval cancelled",
                )) => NsfEvolvingMannError::Cancelled,
                other => other,
            })?;
        debug_assert_eq!(fluid.grid().steps(), final_cfd_step);
        debug_assert_eq!(used_transitions, transitions);
        debug_assert_eq!(used_samples, mode_samples);
        debug_assert_eq!(profile_updates, profiles);
        Ok(NsfEvolvingMannTrial {
            inlet: self,
            fluid,
            report: NsfEvolvingMannIntervalReport {
                steps,
                profile_updates,
                mode_samples: used_samples,
                ou_mode_transitions: used_transitions,
                first: first.expect("positive interval has a first native step"),
                endpoint: endpoint.expect("native interval refreshes its endpoint"),
            },
        })
    }
}
