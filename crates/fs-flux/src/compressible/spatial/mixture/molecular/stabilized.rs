//! RKL2 on the original three-dimensional molecular/SGS spatial operator.
//!
//! The shared gas owner supplies the Legendre coefficients and the shared
//! molecular interval controller owns bounded time/work selection. Only the
//! recurrence of existing conservative fields and boundary transfers belongs
//! here. Every internal state and endpoint retains the original source,
//! Maxwell--Stefan, partial-enthalpy, heat and complete stress evaluations.

use super::*;
use crate::compressible::mixture::molecular::stabilized::{
    self as controller, Action, Rkl2Scheme, reserve_evaluation,
};
use std::cell::Cell;

/// Complete-duration physical transport and actual RKL2 work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularReport {
    /// Original species, momentum and referenced/absolute energy balances.
    /// Accepted steps are RKL2 macrosteps; the parabolic number is unscaled.
    pub transport: MolecularReport,
    /// All actual spatial evaluations, including rejected trials and endpoints.
    pub operator_evaluations: usize,
    /// Largest stage count among accepted macrosteps; zero for zero duration.
    pub maximum_stages_used: usize,
    /// Accepted `dt * rate / stability_extent`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

/// Only the anchor faces and a second boundary recurrence need extra storage.
pub(super) struct Workspace {
    anchor_flux: [Vec<f64>; 3],
    boundary_stage: Vec<f64>,
}
impl Workspace {
    fn required_bytes(grid: &CartesianGrid, width: usize) -> Result<usize, DissipationError> {
        grid.face_len(0)
            .checked_add(grid.face_len(1))
            .and_then(|n| n.checked_add(grid.face_len(2)))
            .and_then(|n| n.checked_add(grid.boundary_len()))
            .and_then(|n| n.checked_mul(width))
            .and_then(|n| n.checked_mul(size_of::<f64>()))
            .ok_or(DissipationError::Invalid(
                "volumetric stabilized workspace size overflow",
            ))
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.anchor_flux.iter().map(Vec::capacity).sum::<usize>() + self.boundary_stage.capacity())
            * size_of::<f64>()
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Base flow/operator bytes plus the lazily prepared RKL2 face/boundary
    /// buffers. Separately selected SGS reserves keep their original charge.
    pub fn required_stabilized_workspace_bytes(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_workspace_bytes(flow, source)?
            .checked_add(Workspace::required_bytes(
                &flow.grid,
                flow.gas.species_count() + 4,
            )?)
            .ok_or(DissipationError::Invalid(
                "volumetric stabilized workspace size overflow",
            ))
    }

    /// Advance the complete interval atomically using second-order RKL2.
    ///
    /// Periodic faces and static prescribed slip/no-slip thermal walls retain
    /// their original transport laws. Optional SGS is recomputed by the same
    /// spatial callback at every state. Finite wall inertia is refused because
    /// its reciprocal carrier solve is not this fixed-wall recurrence.
    ///
    /// `transport.maximum_step_s` controls physical macrostep accuracy; this
    /// path does not estimate temporal error. The enlarged linear parabolic
    /// stability interval proves neither nonlinear positivity nor entropy
    /// decay. Original physical checks admit every stage, with halving on
    /// refusal and no clipping. Euler time and advection reports are unchanged.
    pub fn advance_stabilized(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: StabilizedDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, DissipationError> {
        let report = self.advance_stabilized_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            &mut 0,
            &mut cancelled,
        )?;
        self.publish(report.transport, None);
        Ok(report)
    }

    /// The composing owner retains this work counter through discarded probes.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_stabilized_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: StabilizedDissipationControl,
        maximum_workspace_bytes: usize,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        if self.inertia.is_some() {
            return Err(DissipationError::Invalid(
                "stabilized molecular transport requires prescribed walls without finite inertia",
            ));
        }
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        self.prepare_stabilized(flow, duration.value() > 0.0, maximum_workspace_bytes)?;
        let accepted = controller::advance(
            duration.value(),
            controls,
            evaluations,
            cancelled,
            |action, elapsed, evaluations, cancelled| {
                let stage = match action {
                    Action::First => {
                        return spatial(
                            &mut self.source,
                            &self.grid,
                            self.boundaries,
                            &self.walls,
                            &self.work.current,
                            &mut self.work.primitive,
                            &mut self.work.gradient,
                            &mut self.work.flux,
                            #[cfg(feature = "experimental-sgs")]
                            self.sgs.as_mut(),
                            controls.transport,
                            cancelled,
                        );
                    }
                    Action::Trial {
                        dt,
                        scheme,
                        first,
                        cached,
                    } => self.rkl_trial(
                        dt,
                        scheme,
                        first,
                        cached,
                        controls,
                        elapsed,
                        evaluations,
                        cancelled,
                    )?,
                };
                for (index, ((total, correction), change)) in self
                    .work
                    .boundary_trial
                    .iter_mut()
                    .zip(&mut self.work.boundary_correction)
                    .zip(&self.work.boundary_step)
                    .enumerate()
                {
                    if index % 4096 == 0 {
                        poll(cancelled)?;
                    }
                    let adjusted = change - *correction;
                    let next = *total + adjusted;
                    *correction = (next - *total) - adjusted;
                    *total = next;
                }
                std::mem::swap(&mut self.work.current, &mut self.work.stage);
                Ok(stage)
            },
        )?;
        let transport =
            self.complete_interval(flow, duration, before, accepted.interval, None, cancelled)?;
        Ok(StabilizedMolecularReport {
            transport,
            operator_evaluations: accepted.operator_evaluations,
            maximum_stages_used: accepted.maximum_stages_used,
            maximum_scaled_parabolic_number: accepted.maximum_scaled_parabolic_number,
        })
    }

    fn prepare_stabilized(
        &mut self,
        flow: &MixtureEuler3d,
        needed: bool,
        maximum_workspace_bytes: usize,
    ) -> Result<(), DissipationError> {
        let owned = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "volumetric stabilized workspace size overflow",
            ))?;
        let width = self.gas().species_count() + 4;
        let extra = if needed && self.stabilized.is_none() {
            Workspace::required_bytes(&self.grid, width)?
        } else {
            0
        };
        admit_workspace(
            owned.checked_add(extra).ok_or(DissipationError::Invalid(
                "volumetric stabilized workspace size overflow",
            ))?,
            maximum_workspace_bytes,
        )?;
        if needed && self.stabilized.is_none() {
            let work = Workspace {
                anchor_flux: [
                    allocated(self.grid.face_len(0) * width, 0.0)?,
                    allocated(self.grid.face_len(1) * width, 0.0)?,
                    allocated(self.grid.face_len(2) * width, 0.0)?,
                ],
                boundary_stage: allocated(self.grid.boundary_len() * width, 0.0)?,
            };
            admit_workspace(
                owned
                    .checked_add(work.workspace_bytes())
                    .ok_or(DissipationError::Invalid(
                        "volumetric stabilized workspace size overflow",
                    ))?,
                maximum_workspace_bytes,
            )?;
            self.stabilized = Some(work);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn rkl_trial(
        &mut self,
        dt: f64,
        scheme: Rkl2Scheme,
        first: Stage,
        cached: bool,
        controls: StabilizedDissipationControl,
        elapsed: f64,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        let first = if cached {
            first
        } else {
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            spatial(
                &mut self.source,
                &self.grid,
                self.boundaries,
                &self.walls,
                &self.work.current,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.flux,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                controls.transport,
                cancelled,
            )?
        };
        let maximum = controls.transport.parabolic_cfl * scheme.stability_extent();
        parabolic(dt * first.rate, maximum)?;
        let stable = self.stabilized.as_mut().ok_or(DissipationError::Invalid(
            "missing prepared volumetric stabilized workspace",
        ))?;
        for (source, target) in self.work.flux.iter().zip(&mut stable.anchor_flux) {
            copy_values(source, target, cancelled)?;
        }
        let width = self.source.gas().species_count() + 4;
        update(
            self.source.gas(),
            &self.grid,
            &self.work.current,
            &self.work.current,
            &self.work.flux,
            &mut self.work.stage,
            dt * scheme.first_weight(),
            false,
            None,
            cancelled,
        )?;
        capture(
            &self.grid,
            width,
            &self.work.flux,
            dt * scheme.first_weight(),
            false,
            &mut self.work.boundary_step,
            cancelled,
        )?;
        for tile in stable.boundary_stage.chunks_mut(4096) {
            poll(cancelled)?;
            tile.fill(0.0);
        }
        let mut stats = first;
        for j in 2..=scheme.stage_count() {
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            let stage = spatial(
                &mut self.source,
                &self.grid,
                self.boundaries,
                &self.walls,
                &self.work.stage,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.flux,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                controls.transport,
                cancelled,
            )?;
            parabolic(dt * stage.rate, maximum)?;
            stats = stats.merge(stage);
            let [mu, nu, mu_tilde, gamma_tilde] = scheme.coefficients(j)?;
            let older = Cell::from_mut(self.work.candidate.as_mut_slice()).as_slice_of_cells();
            // Native Cartesian arithmetic owns the divergence. The Cell view
            // safely lets its callbacks consume Y_(j-2) before replacing each
            // admitted row with Y_j, without a fourth whole-field array.
            cartesian::update::<WIDTH_LIMIT, DissipationError, _>(
                &self.grid,
                width,
                dt,
                false,
                |_, _| 0.0,
                |k, d| {
                    let index = k * width + d;
                    let anchor = self.work.current[index];
                    let old = if j == 2 { anchor } else { older[index].get() };
                    anchor + mu * (self.work.stage[index] - anchor) + nu * (old - anchor)
                },
                |axis, face, d| {
                    let index = face * width + d;
                    mu_tilde * self.work.flux[axis][index]
                        + gamma_tilde * stable.anchor_flux[axis][index]
                },
                cancelled,
                |k, row, cancelled| {
                    cell(row, width - 4).primitive(self.source.gas(), &mut *cancelled)?;
                    for (d, value) in row.iter().enumerate() {
                        older[k * width + d].set(*value);
                    }
                    Ok(())
                },
            )?;
            cartesian::capture::<DissipationError>(
                &self.grid,
                width,
                dt,
                false,
                |axis, face, d| {
                    let index = face * width + d;
                    mu_tilde * self.work.flux[axis][index]
                        + gamma_tilde * stable.anchor_flux[axis][index]
                },
                |face, d, value, _| {
                    let index = face * width + d;
                    stable.boundary_stage[index] = mu * self.work.boundary_step[index]
                        + nu * stable.boundary_stage[index]
                        + value;
                },
                cancelled,
            )?;
            std::mem::swap(&mut self.work.stage, &mut self.work.candidate);
            std::mem::swap(&mut self.work.boundary_step, &mut stable.boundary_stage);
        }
        // The final recurrence state must itself qualify actual transport data
        // and its enlarged parabolic bound before private acceptance.
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let endpoint = spatial(
            &mut self.source,
            &self.grid,
            self.boundaries,
            &self.walls,
            &self.work.stage,
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.flux,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            controls.transport,
            cancelled,
        )?;
        parabolic(dt * endpoint.rate, maximum)?;
        Ok(stats.merge(endpoint))
    }
}
