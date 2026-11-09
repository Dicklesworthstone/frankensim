//! Bounded step doubling of one Kelvin-coupled shedding interval.

use super::{
    BemError, CoupledSolveBudget, IncidentField, PreparedTransport, TransportCandidate, WakeSim,
    admit_displacement_fraction, admit_substeps, checkpoint,
};

/// Temporal position-error control with an unchanged shedding discretization.
/// Length and time retain the native wake's unit-chord/unit-freestream units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PanelWakeAdaptiveOptions {
    /// Positive midpoint count for the first coarse integration.
    pub initial_substeps: usize,
    /// Largest permitted count in ANY integration, including the fine result.
    /// At least twice `initial_substeps`; every retry doubles the count.
    pub max_substeps: usize,
    /// Total coupled boundary solves across all trials, including work spent
    /// before a displacement-guard refusal. Charged immediately before solving.
    pub max_coupled_solves: usize,
    /// Finite positive absolute endpoint position tolerance [chord lengths].
    pub absolute_position_tolerance: f64,
    /// Finite nonnegative relative tolerance against the larger of the core
    /// and coarse/fine displacement from the SAME original particle position.
    /// This scale is independent of the coordinate origin.
    pub relative_position_tolerance: f64,
    /// Optional stage-displacement fraction in (0, 1], with the same meaning
    /// as `PanelWakeTransportOptions::max_displacement_fraction`. Exceeding it
    /// retries from the original wake at twice the substep count.
    pub max_displacement_fraction: Option<f64>,
}

/// Accuracy and actual work for one accepted adaptive shedding interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PanelWakeAdaptiveReport {
    /// Midpoint substeps of the accepted fine candidate.
    pub substeps: usize,
    /// All coupled solves performed, including coarse and refused trials.
    pub coupled_solves: usize,
    /// Completed integrations, each including its endpoint boundary solve.
    pub completed_trials: usize,
    /// Incomplete integrations discarded after exceeding the stage guard.
    pub guard_refusals: usize,
    /// Completed coarse/fine comparisons that exceeded the error tolerance.
    pub rejected_comparisons: usize,
    /// Maximum `|fine - coarse| / 3` over endpoint positions [chord lengths].
    pub estimated_max_position_error: f64,
    /// Maximum estimated error divided by its absolute/relative allowance.
    /// Every accepted result has a ratio no larger than one.
    pub max_error_ratio: f64,
    /// Largest sampled displacement fraction in the accepted fine integration.
    pub max_displacement_fraction: f64,
}

impl WakeSim {
    /// Adapt midpoint transport within ONE original shedding interval. Both
    /// coarse and fine trials start from the same retained wake and newborn
    /// packet. Every stage retains older strengths and uses the original
    /// previous-bound Kelvin target; retries never shed additional packets.
    ///
    /// Compare `N` and `2*N` midpoint integrations using the second-order
    /// estimate `|x_fine - x_coarse| / 3` for every endpoint position. Its
    /// allowance is `abs_tol + rel_tol * max(core, coarse_motion, fine_motion)`.
    /// Accept the fine geometry and its endpoint panel solution without
    /// extrapolation. Otherwise reuse that fine candidate as the next coarse
    /// result. A guard refusal discards the incomplete trial and starts a new
    /// adjacent-resolution comparison after refinement.
    ///
    /// This controls an asymptotic temporal position-error estimate, not the
    /// error of discrete shedding, panel/core models, or continuous collision
    /// detection. Incidence stays fixed throughout this interval. Completed
    /// trials share the existing panel factorization and transport scratch.
    /// Solve and substep caps bound all retries; cancellation and any final
    /// refusal leave the accepted physical state and history unchanged.
    ///
    /// # Errors
    /// Native panel-wake refusals, invalid tolerances/settings, exhausted solve
    /// or refinement budgets, nonrepresentable error estimates or cancellation.
    /// Body contact/interior stages and failed boundary solves refuse directly;
    /// only the declared displacement guard and error estimate trigger retries.
    pub fn step_panel_coupled_adaptive(
        &mut self,
        alpha_chord: f64,
        options: PanelWakeAdaptiveOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PanelWakeAdaptiveReport, BemError> {
        self.advance_panel_adaptive(
            alpha_chord,
            options,
            &mut IncidentField::uniform(),
            &mut cancelled,
        )
    }

    /// Adapt the same native wake transport in a prescribed incident field.
    ///
    /// The callback and fixed reference scales have exactly the semantics of
    /// [`Self::step_panel_coupled_with_incident`]. Every coarse/fine trial samples
    /// the same field at its actual panel/wake positions and absolute stage
    /// times. Retrying never advances a forcing clock or sheds another packet.
    /// The callback must reproduce the field for repeated position/time queries;
    /// its external side effects are not rolled back.
    ///
    /// The existing error estimate controls the smooth prescribed-field ODE's
    /// temporal position error, not callback uncertainty, panel/core/shedding
    /// error, external-vorticity evolution or unsteady pressure loads. The
    /// original coupled-solve cap includes all trials, and their panel/vortex
    /// ceilings bound the number of incident samples.
    ///
    /// # Errors
    /// Original adaptive refusals plus callback/nonfinite-flow refusals. The
    /// accepted wake, body solution and history survive every failed trial.
    pub fn step_panel_coupled_adaptive_with_incident(
        &mut self,
        alpha_chord: f64,
        options: PanelWakeAdaptiveOptions,
        mut incident_variation: impl FnMut([f64; 2], f64) -> Result<[f64; 2], BemError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PanelWakeAdaptiveReport, BemError> {
        self.advance_panel_adaptive(
            alpha_chord,
            options,
            &mut IncidentField::new(&mut incident_variation),
            &mut cancelled,
        )
    }

    pub(super) fn advance_panel_adaptive(
        &mut self,
        alpha_chord: f64,
        options: PanelWakeAdaptiveOptions,
        incident: &mut IncidentField<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PanelWakeAdaptiveReport, BemError> {
        checkpoint(cancelled)?;
        self.admit_panel_incidence(alpha_chord)?;
        admit_options(self.dt, options)?;
        let mut transport = PreparedTransport::new(self, alpha_chord, incident, cancelled)?;
        let mut budget = CoupledSolveBudget::new(options.max_coupled_solves);
        let mut substeps = options.initial_substeps;
        let mut coarse = None;
        let mut completed_trials = 0;
        let mut guard_refusals = 0;
        let mut rejected_comparisons = 0;
        loop {
            checkpoint(cancelled)?;
            let candidate = match transport.integrate(
                substeps,
                options.max_displacement_fraction,
                &mut budget,
                incident,
                cancelled,
            ) {
                Ok(candidate) => candidate,
                Err(BemError::InvalidScalar {
                    name: "panel wake stage displacement fraction",
                    ..
                }) => {
                    guard_refusals += 1;
                    // Never compare nonadjacent resolutions with the factor 3.
                    coarse = None;
                    substeps = doubled_substeps(transport.dt, substeps, options.max_substeps)?;
                    continue;
                }
                Err(error) => return Err(error),
            };
            completed_trials += 1;
            if let Some(previous) = coarse.take() {
                let (error, ratio) =
                    compare(&transport, &previous, &candidate, options, cancelled)?;
                if ratio <= 1.0 {
                    let report = PanelWakeAdaptiveReport {
                        substeps,
                        coupled_solves: budget.used,
                        completed_trials,
                        guard_refusals,
                        rejected_comparisons,
                        estimated_max_position_error: error,
                        max_error_ratio: ratio,
                        max_displacement_fraction: candidate.report.max_displacement_fraction,
                    };
                    drop(transport);
                    self.publish_panel_transport(candidate, cancelled)?;
                    return Ok(report);
                }
                rejected_comparisons += 1;
            }
            coarse = Some(candidate);
            substeps = doubled_substeps(transport.dt, substeps, options.max_substeps)?;
        }
    }
}

fn admit_options(dt: f64, options: PanelWakeAdaptiveOptions) -> Result<(), BemError> {
    admit_substeps(dt, options.initial_substeps)?;
    admit_displacement_fraction(options.max_displacement_fraction)?;
    if !options.absolute_position_tolerance.is_finite()
        || options.absolute_position_tolerance <= 0.0
    {
        return Err(BemError::InvalidScalar {
            name: "panel wake absolute position tolerance",
            value: options.absolute_position_tolerance,
            requirement: "finite and positive in chord lengths",
        });
    }
    if !options.relative_position_tolerance.is_finite() || options.relative_position_tolerance < 0.0
    {
        return Err(BemError::InvalidScalar {
            name: "panel wake relative position tolerance",
            value: options.relative_position_tolerance,
            requirement: "finite and nonnegative",
        });
    }
    doubled_substeps(dt, options.initial_substeps, options.max_substeps)?;
    Ok(())
}

fn doubled_substeps(dt: f64, count: usize, maximum: usize) -> Result<usize, BemError> {
    let next = count.checked_mul(2);
    if next.is_none_or(|next| next > maximum) {
        return Err(BemError::WorkEnvelopeExceeded {
            operation: "panel wake adaptive substeps",
            requested: next.unwrap_or(usize::MAX),
            max: maximum,
        });
    }
    let next = next.expect("admitted doubled count");
    admit_substeps(dt, next)?;
    Ok(next)
}

fn compare(
    transport: &PreparedTransport<'_>,
    coarse: &TransportCandidate,
    fine: &TransportCandidate,
    options: PanelWakeAdaptiveOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), BemError> {
    let initial = transport
        .retained
        .iter()
        .chain(std::iter::once(&transport.newborn));
    let mut error = 0.0_f64;
    let mut ratio = 0.0_f64;
    for ((start, coarse), fine) in initial.zip(&coarse.wake).zip(&fine.wake) {
        checkpoint(cancelled)?;
        let (local_error, local_ratio) = position_error(
            start.pos,
            coarse.pos,
            fine.pos,
            transport.core,
            options.absolute_position_tolerance,
            options.relative_position_tolerance,
        )?;
        error = error.max(local_error);
        ratio = ratio.max(local_ratio);
    }
    Ok((error, ratio))
}

fn position_error(
    initial: [f64; 2],
    coarse: [f64; 2],
    fine: [f64; 2],
    core: f64,
    absolute: f64,
    relative: f64,
) -> Result<(f64, f64), BemError> {
    let distance = |a: [f64; 2], b: [f64; 2]| fs_math::det::hypot(a[0] - b[0], a[1] - b[1]);
    let error = distance(coarse, fine) / 3.0;
    let motion = distance(initial, coarse).max(distance(initial, fine));
    let allowance = absolute + relative * core.max(motion);
    if !error.is_finite() || !allowance.is_finite() || allowance <= 0.0 {
        return Err(BemError::NonFiniteResult {
            operation: "panel wake adaptive error estimate",
        });
    }
    // A very small positive tolerance may yield +infinity; that rejects the
    // comparison and refines under the same bounded work/substep envelope.
    Ok((error, error / allowance))
}

#[cfg(test)]
mod tests {
    use super::position_error;

    #[test]
    fn relative_error_scale_depends_on_motion_and_core_not_coordinate_origin() {
        let initial = [1.0, -2.0];
        let coarse = [1.25, -1.5];
        let fine = [1.125, -1.625];
        let reference = position_error(initial, coarse, fine, 0.5, 1e-6, 1e-3).unwrap();
        for offset in [[128.0, -256.0], [-64.0, 512.0]] {
            let shift = |p: [f64; 2]| [p[0] + offset[0], p[1] + offset[1]];
            assert_eq!(
                position_error(shift(initial), shift(coarse), shift(fine), 0.5, 1e-6, 1e-3)
                    .unwrap(),
                reference,
            );
        }
        let near = position_error([0.0; 2], [0.0; 2], [1e-6, 0.0], 0.5, 1e-6, 1e-3).unwrap();
        assert_eq!(near.1, near.0 / (1e-6 + 1e-3 * 0.5));
    }
}
