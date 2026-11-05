//! Accuracy-controlled RKL2 on the original planar field and face operators.

use super::super::stabilized::{Action, reserve_evaluation};
use super::super::{
    DissipationError, MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler2d, MolecularTransport2d,
    PlanarTransport, StabilizedMolecularReport, TransportLaw, copy_values, poll,
};
use super::{
    AdaptiveDissipationControl, PhysicalState, Workspace, advance, face_error, state_error,
};
use crate::compressible::{AxisBoundaries, BoundarySide, SIDES};
use fs_qty::Time;

/// Original balances, actual full/two-half work and accepted local accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMolecularReport {
    /// Original physical transport. Accepted steps count pairs of half steps;
    /// all coarse, fine and refused RHS evaluations enter the work account.
    pub stabilized: StabilizedMolecularReport,
    /// Largest accepted component or individual-face error ratio, at most one.
    pub maximum_estimated_local_error: f64,
    /// Probe refusals caused by estimated temporal error rather than physics.
    pub temporal_rejections: usize,
}

impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Complete flow/operator RKL2 storage plus two field-sized comparison
    /// buffers and two boundary-quadrature buffers. Pass this allowance to
    /// `new`; accuracy scratch is allocated only for a nonzero first interval.
    pub fn required_adaptive_workspace_bytes(
        flow: &MixtureEuler2d,
        source: &PlanarTransport<L>,
    ) -> Result<usize, DissipationError> {
        let [nx, ny] = flow.grid.shape;
        let boundary = nx
            .checked_add(ny)
            .and_then(|n| n.checked_mul(2 * (flow.gas.species_count() + 3)))
            .ok_or(DissipationError::Invalid(
                "adaptive planar workspace size overflow",
            ))?;
        Self::required_stabilized_workspace_bytes(flow, source)?
            .checked_add(Workspace::required_bytes(flow.states.len(), boundary)?)
            .ok_or(DissipationError::Invalid(
                "adaptive planar workspace size overflow",
            ))
    }

    /// Advance a complete molecular interval to explicit local physical accuracy.
    ///
    /// One original RKL2 scheme serves the full step and both half steps. The
    /// admitted two-half field and its own signed face transfers are retained,
    /// without extrapolation, clipping or balance repair. Every original species,
    /// temperature, velocity component and individual boundary transfer is
    /// checked, so large formation energies and cancelling faces cannot hide
    /// local error. The initial and midpoint RHS are reused: an unrefused
    /// s-stage probe costs 3*s+1 source evaluations.
    ///
    /// These are estimated local molecular errors, not global, spatial or
    /// splitting error bounds. Original caloric/constitutive admission remains
    /// mandatory at every stage. Failed evolution preserves the accepted field,
    /// Euler clock, reports and published face histories.
    pub fn advance_adaptive(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        controls: AdaptiveDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveMolecularReport, DissipationError> {
        let report = self.advance_adaptive_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            &mut 0,
            &mut cancelled,
        )?;
        self.publish(report.stabilized.transport);
        Ok(report)
    }

    /// The composing owner may retain this counter across several intervals.
    pub(in crate::compressible::mixture) fn advance_adaptive_unpublished(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        controls: AdaptiveDissipationControl,
        maximum_workspace_bytes: usize,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AdaptiveMolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        if duration.value() > 0.0 {
            self.reserve_adaptive_workspace(flow, maximum_workspace_bytes, cancelled)?;
        }
        let stabilized = controls.stabilized;
        let result = advance(
            duration.value(),
            controls,
            evaluations,
            cancelled,
            |action, elapsed, evaluations, cancelled| match action {
                Action::First => {
                    copy_values(
                        &self.work.current,
                        &mut self.work.adaptive.anchor,
                        cancelled,
                    )?;
                    Ok((self.stabilized_anchor(flow, stabilized, cancelled)?, 0.0))
                }
                Action::Trial {
                    dt,
                    scheme,
                    first,
                    cached,
                } => {
                    let attempt = (|| {
                        let anchor = if cached {
                            first
                        } else {
                            reserve_evaluation(
                                evaluations,
                                stabilized.maximum_operator_evaluations,
                                elapsed,
                            )?;
                            self.stabilized_anchor(flow, stabilized, cancelled)?
                        };
                        let coarse = self.stabilized_trial(
                            flow,
                            dt,
                            scheme,
                            anchor,
                            stabilized,
                            elapsed,
                            evaluations,
                            cancelled,
                        )?;
                        self.capture_coarse_coordinates(cancelled)?;
                        copy_values(
                            &self.work.boundary_step,
                            &mut self.work.adaptive.coarse_boundary,
                            cancelled,
                        )?;
                        // The original anchor flux was retained by the coarse
                        // recurrence; the first fine half can use it directly.
                        let first = self.stabilized_trial(
                            flow,
                            0.5 * dt,
                            scheme,
                            anchor,
                            stabilized,
                            elapsed,
                            evaluations,
                            cancelled,
                        )?;
                        copy_values(
                            &self.work.boundary_step,
                            &mut self.work.adaptive.first_boundary,
                            cancelled,
                        )?;
                        core::mem::swap(&mut self.work.current, &mut self.work.stage);
                        self.cache_stabilized_anchor(cancelled)?;
                        let second = self.stabilized_trial(
                            flow,
                            0.5 * dt,
                            scheme,
                            first,
                            stabilized,
                            elapsed,
                            evaluations,
                            cancelled,
                        )?;
                        let error = self.adaptive_error(controls, cancelled)?;
                        let mut peak = coarse.merge(first).merge(second);
                        peak.rate = coarse.rate.max(0.5 * first.rate).max(0.5 * second.rate);
                        Ok((peak, error))
                    })();
                    match attempt {
                        Ok((peak, error)) if error <= 1.0 => {
                            // Preserve both oriented quadratures in the original
                            // compensated total; do not subtract cancelling faces.
                            for i in 0..self.work.boundary_trial.len() {
                                if i % 4096 == 0 {
                                    poll(cancelled)?;
                                }
                                for change in [
                                    self.work.adaptive.first_boundary[i],
                                    self.work.boundary_step[i],
                                ] {
                                    let adjusted = change - self.work.boundary_correction[i];
                                    let next = self.work.boundary_trial[i] + adjusted;
                                    self.work.boundary_correction[i] =
                                        (next - self.work.boundary_trial[i]) - adjusted;
                                    self.work.boundary_trial[i] = next;
                                }
                            }
                            core::mem::swap(&mut self.work.current, &mut self.work.stage);
                            Ok((peak, error))
                        }
                        refused => {
                            copy_values(
                                &self.work.adaptive.anchor,
                                &mut self.work.current,
                                cancelled,
                            )?;
                            refused
                        }
                    }
                }
            },
        )?;
        let accepted = result.stabilized;
        let transport =
            self.complete_interval(flow, duration, before, accepted.interval, cancelled)?;
        Ok(AdaptiveMolecularReport {
            stabilized: StabilizedMolecularReport {
                transport,
                operator_evaluations: accepted.operator_evaluations,
                maximum_stages_used: accepted.maximum_stages_used,
                maximum_scaled_parabolic_number: accepted.maximum_scaled_parabolic_number,
            },
            maximum_estimated_local_error: result.maximum_estimated_local_error,
            temporal_rejections: result.temporal_rejections,
        })
    }

    fn reserve_adaptive_workspace(
        &mut self,
        flow: &MixtureEuler2d,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let states = self.work.current.len();
        let boundary = self.work.boundary_step.len();
        let adaptive = self.work.adaptive.additional_bytes(states, boundary)?;
        let extra = adaptive
            .checked_add(self.additional_stabilized_workspace_bytes()?)
            .ok_or(DissipationError::Invalid(
                "adaptive planar workspace size overflow",
            ))?;
        let limit = maximum_workspace_bytes.min(MAX_MIXTURE_WORKSPACE_BYTES);
        let required = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .and_then(|n| n.checked_add(extra))
            .ok_or(DissipationError::Invalid(
                "adaptive planar workspace size overflow",
            ))?;
        if required > limit {
            return Err(DissipationError::Invalid(
                "adaptive planar molecular total workspace budget",
            ));
        }
        // Reserve both pending groups before allocation, then retain room for
        // adaptive buffers while checking the actual RKL allocation capacities.
        self.reserve_stabilized_workspace(flow, limit - adaptive, cancelled)?;
        self.work.adaptive.reserve(states, boundary, cancelled)?;
        if self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .is_none_or(|n| n > limit)
        {
            return Err(DissipationError::Invalid(
                "adaptive planar molecular actual workspace budget",
            ));
        }
        poll(cancelled)
    }

    fn capture_coarse_coordinates(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let n = self.gas().species_count();
        let width = n + 3;
        for (index, (source, target)) in self
            .work
            .stage
            .chunks_exact(width)
            .zip(self.work.adaptive.coarse.chunks_exact_mut(width))
            .enumerate()
        {
            poll(cancelled)?;
            let primitive = self.work.primitive[index];
            target[..n].copy_from_slice(&source[..n]);
            target[n..n + 2].copy_from_slice(&primitive.velocity);
            target[n + 2] = primitive.temperature;
        }
        Ok(())
    }

    fn adaptive_error(
        &self,
        controls: AdaptiveDissipationControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DissipationError> {
        let gas = self.gas();
        let n = gas.species_count();
        let width = n + 3;
        let mut error = 0.0_f64;
        for (index, (coarse, fine)) in self
            .work
            .adaptive
            .coarse
            .chunks_exact(width)
            .zip(self.work.stage.chunks_exact(width))
            .enumerate()
        {
            poll(cancelled)?;
            let p = self.work.primitive[index];
            error = error.max(state_error(
                gas,
                coarse,
                fine,
                PhysicalState {
                    density: coarse[..n].iter().sum(),
                    velocity: [coarse[n], coarse[n + 1]],
                    temperature: coarse[n + 2],
                },
                PhysicalState {
                    density: p.density,
                    velocity: p.velocity,
                    temperature: p.temperature,
                },
                controls,
            )?);
        }
        let [nx, ny] = self.grid.shape;
        for side in SIDES {
            if self.boundaries.axis(side.axis()) == AxisBoundaries::Periodic {
                continue;
            }
            for (line, face) in self.grid.boundary_range(side).enumerate() {
                poll(cancelled)?;
                let cell = match side {
                    BoundarySide::XMin => nx * line,
                    BoundarySide::XMax => nx * line + nx - 1,
                    BoundarySide::YMin => line,
                    BoundarySide::YMax => nx * (ny - 1) + line,
                };
                let at = face * width;
                error = error.max(face_error::<2>(
                    gas,
                    &self.work.adaptive.coarse_boundary[at..at + width],
                    &self.work.adaptive.first_boundary[at..at + width],
                    &self.work.boundary_step[at..at + width],
                    &self.work.stage[cell * width..cell * width + n],
                    self.work.primitive[cell].temperature,
                    self.grid.cell_volume(),
                    controls,
                    cancelled,
                )?);
            }
        }
        Ok(error)
    }
}
