//! RKL2 for the original three-dimensional frozen-composition operator.
//! Full stress/work, Fourier heat, authored walls and selected SGS closures
//! retain their original actual-state evaluations and source admission.

use super::*;
use crate::compressible::mixture::molecular::{
    interval::Stage,
    stabilized::{self as controller, Action, Rkl2Scheme, reserve_evaluation},
};
use crate::compressible::molecular::stabilized::{admit_workspace, copy_values};
use fs_gas::multicomponent::dissipation::StabilizedDissipationControl;
use std::cell::Cell;
use std::mem::size_of;

/// Complete physical balances and actual three-dimensional RKL2 work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularReport3d {
    /// Original conservative report; accepted steps are RKL2 macrosteps.
    pub transport: MolecularReport3d,
    /// All attempted spatial evaluations, including retries and endpoints.
    pub operator_evaluations: usize,
    /// Largest stage count among accepted macrosteps.
    pub maximum_stages_used: usize,
    /// Largest accepted `dt * rate / stability_extent`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

pub(super) struct Workspace {
    anchor_flux: [Vec<SpatialFaceFlux>; 3],
    boundary_stage: Vec<BoundaryImpulse>,
}
impl Workspace {
    fn bytes(&self) -> usize {
        self.anchor_flux.iter().map(Vec::capacity).sum::<usize>() * size_of::<SpatialFaceFlux>()
            + self.boundary_stage.capacity() * size_of::<BoundaryImpulse>()
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Required EXTRA anchor-face/boundary bytes, excluding existing storage.
    pub fn required_stabilized_workspace_bytes(&self) -> Result<usize, MolecularError> {
        self.work.flux[0]
            .len()
            .checked_add(self.work.flux[1].len())
            .and_then(|n| n.checked_add(self.work.flux[2].len()))
            .and_then(|n| n.checked_mul(size_of::<SpatialFaceFlux>()))
            .and_then(|n| {
                self.work
                    .boundary_step
                    .len()
                    .checked_mul(size_of::<BoundaryImpulse>())
                    .and_then(|b| n.checked_add(b))
            })
            .ok_or(MolecularError::Invalid(
                "volumetric stabilized workspace size overflow",
            ))
    }

    /// Actual retained capacity of the additional RKL2 buffers.
    pub fn stabilized_workspace_bytes(&self) -> usize {
        self.stabilized.as_ref().map_or(0, Workspace::bytes)
    }

    /// Advance the entire duration atomically with bounded second-order RKL2.
    ///
    /// The native source, caloric, full stress/work, Fourier and wall operators
    /// qualify every stage and its final endpoint. Selected SGS is reevaluated
    /// at those same states. Boundary impulses use the field's exact Legendre
    /// recurrence. Three existing state arrays and reusable extra face/boundary
    /// buffers suffice; no per-stage field allocation is introduced.
    ///
    /// Linear diffusion stability is not a positivity or entropy certificate.
    /// Refused physical trials halve without clipping. The caller controls
    /// time accuracy with `controls.maximum_step_s`; no error estimate is
    /// implied. Accepted Euler time and advection metadata remain unchanged.
    pub fn advance_stabilized(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        stabilization: MolecularStabilization,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport3d, MolecularError> {
        self.advance_stabilized_with_budget(
            flow,
            duration,
            controls,
            stabilization,
            &mut 0,
            timed::Clock::interval(0.0, duration.value()),
            None,
            &mut cancelled,
        )
    }

    pub(super) fn prepare_stabilized(
        &mut self,
        maximum_bytes: usize,
    ) -> Result<(), MolecularError> {
        if self.stabilized.is_none() {
            admit_workspace(self.required_stabilized_workspace_bytes()?, maximum_bytes)?;
            let work = Workspace {
                anchor_flux: [
                    allocated(self.work.flux[0].len(), SpatialFaceFlux::default())?,
                    allocated(self.work.flux[1].len(), SpatialFaceFlux::default())?,
                    allocated(self.work.flux[2].len(), SpatialFaceFlux::default())?,
                ],
                boundary_stage: allocated(self.work.boundary_step.len(), Totals::default())?,
            };
            admit_workspace(work.bytes(), maximum_bytes)?;
            self.stabilized = Some(work);
        }
        admit_workspace(self.stabilized_workspace_bytes(), maximum_bytes)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_stabilized_with_budget(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        stabilization: MolecularStabilization,
        evaluations: &mut usize,
        clock: timed::Clock,
        mut profiles: Option<&mut timed::TimedSource<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport3d, MolecularError> {
        poll(cancelled)?;
        let selected = stabilization.controller(controls)?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid(
                "molecular operator and flow differ",
            ));
        }
        if !positive(duration.value()) {
            return Err(MolecularError::Invalid("molecular duration"));
        }
        if profiles.is_some() {
            clock.admit(duration.value())?;
        }
        self.prepare_stabilized(stabilization.maximum_extra_workspace_bytes)?;
        let before = checked_totals(&flow.states, self.grid.cell_volume(), cancelled)?;
        copy_values(&flow.states, &mut self.work.current, cancelled)?;
        for tile in self.work.boundary_trial.chunks_mut(256) {
            poll(cancelled)?;
            tile.fill(Totals::default());
        }
        let accepted = controller::advance(
            duration.value(),
            selected,
            evaluations,
            cancelled,
            |action, elapsed, evaluations, cancelled| {
                let result = (|| {
                    let rate = match action {
                        Action::First => {
                            let reservoirs = timed::stage(
                                &mut flow.reservoirs,
                                profiles.as_deref_mut(),
                                &self.gas,
                                &self.grid,
                                clock.at(elapsed, duration.value()),
                                self.reservoir_exchange,
                                cancelled,
                            )?;
                            return spatial(
                                &self.gas,
                                &self.grid,
                                self.boundaries,
                                &self.walls,
                                reservoirs,
                                &self.work.current,
                                &mut self.law,
                                &mut self.work.primitive,
                                &mut self.work.gradient,
                                &mut self.work.viscosity,
                                &mut self.work.flux,
                                &mut self.work.binary,
                                #[cfg(feature = "experimental-sgs")]
                                self.sgs.as_mut(),
                                cancelled,
                            )
                            .map(|rate| Stage {
                                rate,
                                ..Stage::default()
                            });
                        }
                        Action::Trial {
                            dt,
                            scheme,
                            first,
                            cached,
                        } => self.rkl_trial(
                            flow,
                            profiles.as_deref_mut(),
                            clock,
                            duration.value(),
                            dt,
                            scheme,
                            first.rate,
                            cached,
                            selected,
                            elapsed,
                            evaluations,
                            cancelled,
                        )?,
                    };
                    for (index, (total, change)) in self
                        .work
                        .boundary_trial
                        .iter_mut()
                        .zip(&self.work.boundary_step)
                        .enumerate()
                    {
                        if index % 256 == 0 {
                            poll(cancelled)?;
                        }
                        *total = add_totals(*total, *change);
                    }
                    std::mem::swap(&mut self.work.current, &mut self.work.stage);
                    Ok(Stage {
                        rate,
                        ..Stage::default()
                    })
                })();
                result.map_err(|error| timed::AttemptFailure {
                    error,
                    source_failed: profiles.as_deref().is_some_and(timed::TimedSource::failed),
                })
            },
        )
        .map_err(|failure| failure.error)?;
        let after = checked_totals(&self.work.current, self.grid.cell_volume(), cancelled)?;
        let outward = sum_arrays(
            self.work.boundary_trial.iter().map(|q| q.array()),
            cancelled,
        )?;
        let boundary_change = Totals::from_array(outward.map(|v| -v));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(cancelled)?;
        flow.states.copy_from_slice(&self.work.current);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        if profiles.is_some() {
            flow.reservoirs.profiles.publish();
        }
        Ok(StabilizedMolecularReport3d {
            transport: MolecularReport3d {
                elapsed: duration,
                accepted_steps: accepted.interval.accepted_steps,
                rejected_steps: accepted.interval.rejected_steps,
                maximum_parabolic_number: accepted.interval.maximum_parabolic_number,
                before,
                after,
                boundary_change,
                conservation_residual,
            },
            operator_evaluations: accepted.operator_evaluations,
            maximum_stages_used: accepted.maximum_stages_used,
            maximum_scaled_parabolic_number: accepted.maximum_scaled_parabolic_number,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn rkl_trial(
        &mut self,
        flow: &mut Euler3d,
        mut profiles: Option<&mut timed::TimedSource<'_>>,
        clock: timed::Clock,
        duration: f64,
        dt: f64,
        scheme: Rkl2Scheme,
        first: f64,
        cached: bool,
        controls: StabilizedDissipationControl,
        elapsed: f64,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MolecularError> {
        let first = if cached {
            first
        } else {
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            let reservoirs = timed::stage(
                &mut flow.reservoirs,
                profiles.as_deref_mut(),
                &self.gas,
                &self.grid,
                clock.at(elapsed, duration),
                self.reservoir_exchange,
                cancelled,
            )?;
            spatial(
                &self.gas,
                &self.grid,
                self.boundaries,
                &self.walls,
                reservoirs,
                &self.work.current,
                &mut self.law,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.viscosity,
                &mut self.work.flux,
                &mut self.work.binary,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                cancelled,
            )?
        };
        let maximum = controls.transport.parabolic_cfl * scheme.stability_extent();
        parabolic(dt * first, maximum)?;
        let stable = self.stabilized.as_mut().ok_or(MolecularError::Invalid(
            "missing prepared volumetric stabilized workspace",
        ))?;
        for (source, target) in self.work.flux.iter().zip(&mut stable.anchor_flux) {
            copy_values(source, target, cancelled)?;
        }
        evolve::update(
            &self.gas,
            &self.grid,
            &self.work.current,
            &self.work.current,
            &self.work.flux,
            &mut self.work.stage,
            dt * scheme.first_weight(),
            false,
            None,
            None, // Do not apply advective gravity again in the diffusion split.
            cancelled,
        )?;
        evolve::capture_boundaries(
            &self.grid,
            &self.work.flux,
            dt * scheme.first_weight(),
            &mut self.work.boundary_step,
            false,
            cancelled,
        )?;
        for tile in stable.boundary_stage.chunks_mut(256) {
            poll(cancelled)?;
            tile.fill(Totals::default());
        }
        let mut rate = first;
        for j in 2..=scheme.stage_count() {
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            let reservoirs = timed::stage(
                &mut flow.reservoirs,
                profiles.as_deref_mut(),
                &self.gas,
                &self.grid,
                clock.at(elapsed + dt * scheme.abscissa(j - 1)?, duration),
                self.reservoir_exchange,
                cancelled,
            )?;
            let stage_rate = spatial(
                &self.gas,
                &self.grid,
                self.boundaries,
                &self.walls,
                reservoirs,
                &self.work.stage,
                &mut self.law,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.viscosity,
                &mut self.work.flux,
                &mut self.work.binary,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                cancelled,
            )?;
            parabolic(dt * stage_rate, maximum)?;
            rate = rate.max(stage_rate);
            let [mu, nu, mu_tilde, gamma_tilde] = scheme.coefficients(j)?;
            let older = Cell::from_mut(self.work.candidate.as_mut_slice()).as_slice_of_cells();
            cartesian::update::<5, Error, _>(
                &self.grid,
                5,
                dt,
                false,
                |_, _| 0.0,
                |k, d| {
                    let anchor = self.work.current[k].array()[d];
                    let old = if j == 2 {
                        anchor
                    } else {
                        older[k].get().array()[d]
                    };
                    anchor + mu * (self.work.stage[k].array()[d] - anchor) + nu * (old - anchor)
                },
                |axis, face, d| {
                    mu_tilde * evolve::flux_array(self.work.flux[axis][face], axis)[d]
                        + gamma_tilde * evolve::flux_array(stable.anchor_flux[axis][face], axis)[d]
                },
                cancelled,
                |k, row, cancelled| {
                    let q = State::from_array([row[0], row[1], row[2], row[3], row[4]]);
                    q.primitive(&self.gas, &mut *cancelled)?;
                    older[k].set(q);
                    Ok(())
                },
            )?;
            cartesian::capture::<Error>(
                &self.grid,
                5,
                dt,
                false,
                |axis, face, d| {
                    mu_tilde * evolve::flux_array(self.work.flux[axis][face], axis)[d]
                        + gamma_tilde * evolve::flux_array(stable.anchor_flux[axis][face], axis)[d]
                },
                |face, d, value, _| {
                    let previous = self.work.boundary_step[face].array()[d];
                    let mut target = stable.boundary_stage[face].array();
                    target[d] = mu * previous + nu * target[d] + value;
                    stable.boundary_stage[face] = Totals::from_array(target);
                },
                cancelled,
            )?;
            std::mem::swap(&mut self.work.stage, &mut self.work.candidate);
            std::mem::swap(&mut self.work.boundary_step, &mut stable.boundary_stage);
        }
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let endpoint = if dt == duration - elapsed {
            duration
        } else {
            elapsed + dt
        };
        let reservoirs = timed::stage(
            &mut flow.reservoirs,
            profiles.as_deref_mut(),
            &self.gas,
            &self.grid,
            clock.at(endpoint, duration),
            self.reservoir_exchange,
            cancelled,
        )?;
        let endpoint_rate = spatial(
            &self.gas,
            &self.grid,
            self.boundaries,
            &self.walls,
            reservoirs,
            &self.work.stage,
            &mut self.law,
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.viscosity,
            &mut self.work.flux,
            &mut self.work.binary,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            cancelled,
        )?;
        parabolic(dt * endpoint_rate, maximum)?;
        Ok(rate.max(endpoint_rate))
    }
}
