//! Bounded RKL2 integration of the original frozen-composition face operator.
//!
//! The gas owner supplies Legendre coefficients, the shared interval controller
//! selects work/time, and the original Cartesian divergence owns conservation.
//! No transport law, heat capacity, source gate or stress closure is replaced.

use super::*;
use crate::compressible::mixture::molecular::{
    interval::Stage,
    stabilized::{self as controller, Action, IntervalFailure, Rkl2Scheme, reserve_evaluation},
};
use fs_gas::multicomponent::dissipation::{DissipationControl, StabilizedDissipationControl};
use std::cell::Cell;
use std::mem::size_of;

/// Additional RKL2 controls. Physical time, CFL and retry limits remain in
/// `MolecularControl`, also when composing molecular/advection/molecular steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStabilization {
    /// Maximum Legendre stages in one macrostep, in [2, 256].
    pub maximum_stages: usize,
    /// Actual spatial calls, including failed prefixes and endpoint checks.
    /// In [1, 1,000,000]; a successful s-stage macrostep costs s+1 calls.
    pub maximum_operator_evaluations: usize,
    /// Byte cap for the EXTRA retained anchor faces and boundary recurrence.
    /// This is not a cap on the existing flow, law, molecular or SGS storage.
    pub maximum_extra_workspace_bytes: usize,
}
impl Default for MolecularStabilization {
    fn default() -> Self {
        Self {
            maximum_stages: 32,
            maximum_operator_evaluations: 100_000,
            maximum_extra_workspace_bytes: 64 * 1024 * 1024,
        }
    }
}
impl MolecularStabilization {
    pub(in crate::compressible) fn controller(
        self,
        transport: MolecularControl,
    ) -> Result<StabilizedDissipationControl, MolecularError> {
        transport.validate()?;
        let controls = StabilizedDissipationControl {
            transport: DissipationControl {
                parabolic_cfl: transport.parabolic_cfl,
                maximum_step_s: transport.maximum_step_s,
                minimum_step_s: transport.minimum_step_s,
                maximum_steps: transport.maximum_steps,
                maximum_halvings: transport.maximum_halvings,
                // The common interval controller reads time/work controls only.
                // Frozen composition never invokes a Maxwell--Stefan solve.
                ..DissipationControl::default()
            },
            maximum_stages: self.maximum_stages,
            maximum_operator_evaluations: self.maximum_operator_evaluations,
        };
        controls.validate()?;
        Ok(controls)
    }
}

/// Private composing-step account. Only extrema reset when a whole split is
/// retried; spent operator work is never refunded.
pub(in crate::compressible) struct StepWork {
    pub controls: MolecularStabilization,
    pub evaluations: usize,
    pub maximum_stages_used: usize,
    pub maximum_scaled_parabolic_number: f64,
}
impl StepWork {
    pub fn new(controls: MolecularStabilization) -> Self {
        Self {
            controls,
            evaluations: 0,
            maximum_stages_used: 0,
            maximum_scaled_parabolic_number: 0.0,
        }
    }

    pub fn reset_extrema(&mut self) {
        self.maximum_stages_used = 0;
        self.maximum_scaled_parabolic_number = 0.0;
    }

    pub fn record(&mut self, stages: usize, number: f64) {
        self.maximum_stages_used = self.maximum_stages_used.max(stages);
        self.maximum_scaled_parabolic_number = self.maximum_scaled_parabolic_number.max(number);
    }
}

impl IntervalFailure for MolecularError {
    /// A cancellation, an operator work limit, or a refused subgrid transport ratio ends the
    /// interval with its own error. The ratio refusals (`SubgridViscosityRatio`,
    /// `SubgridThermalDiffusivityRatio`) compare the fitted eddy transport with the molecular
    /// transport; the held coefficient and the face strain set that ratio, not the step length,
    /// so halving the step cannot admit it and every coupling retry would be futile work.
    fn is_terminal(&self) -> bool {
        #[cfg(feature = "experimental-sgs")]
        if matches!(
            self,
            Self::SubgridViscosityRatio { .. } | Self::SubgridThermalDiffusivityRatio { .. }
        ) {
            return true;
        }
        matches!(
            self,
            Self::Cancelled
                | Self::Law(
                    DissipationError::OperatorWorkLimit { .. } | DissipationError::WorkLimit { .. }
                )
        )
    }
    fn step_failure(self, elapsed: f64, dt: f64) -> Self {
        Self::StepFailure {
            elapsed_s: elapsed,
            step_s: dt,
            last: Box::new(self),
        }
    }
}

/// Complete physical balances and actual work in a stabilized interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularReport {
    /// Original conservative transport report; accepted steps are macrosteps.
    pub transport: MolecularReport,
    /// All attempted operator evaluations, including rejected stages.
    pub operator_evaluations: usize,
    /// Largest accepted macrostep stage count.
    pub maximum_stages_used: usize,
    /// Largest accepted `dt * rate / stability_extent`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

pub(super) struct Workspace {
    anchor_flux: [Vec<NormalFaceFlux>; 2],
    boundary_stage: Vec<BoundaryImpulse>,
}
impl Workspace {
    fn bytes(&self) -> usize {
        self.anchor_flux.iter().map(Vec::capacity).sum::<usize>() * size_of::<NormalFaceFlux>()
            + self.boundary_stage.capacity() * size_of::<BoundaryImpulse>()
    }
}

pub(in crate::compressible) fn admit_workspace(
    required: usize,
    maximum: usize,
) -> Result<(), MolecularError> {
    if required > maximum {
        Err(MolecularError::StabilizedWorkspace {
            required_bytes: required,
            maximum_bytes: maximum,
        })
    } else {
        Ok(())
    }
}

pub(in crate::compressible) fn copy_values<T: Copy>(
    source: &[T],
    target: &mut [T],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MolecularError> {
    for (source, target) in source.chunks(256).zip(target.chunks_mut(256)) {
        poll(cancelled)?;
        target.copy_from_slice(source);
    }
    Ok(())
}

impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Additional RKL2 storage required for this prepared operator, in bytes.
    /// Existing flow/law/operator storage is outside this explicit extra cap.
    pub fn required_stabilized_workspace_bytes(&self) -> Result<usize, MolecularError> {
        self.work.x_flux.len().checked_add(self.work.y_flux.len())
            .and_then(|n| n.checked_mul(size_of::<NormalFaceFlux>()))
            .and_then(|n| self.work.boundary_step.len()
                .checked_mul(size_of::<BoundaryImpulse>())
                .and_then(|b| n.checked_add(b)))
            .ok_or(MolecularError::Invalid("planar stabilized workspace size overflow"))
    }

    /// Actual retained capacity of the additional RKL2 buffers.
    pub fn stabilized_workspace_bytes(&self) -> usize {
        self.stabilized.as_ref().map_or(0, Workspace::bytes)
    }

    /// Advance the complete duration atomically with second-order RKL2.
    ///
    /// Original source laws, full viscous stress/work, Fourier heat and wall
    /// conditions are evaluated at every stage and at the final endpoint.
    /// Boundary transfers use the SAME Legendre weights as the field.
    /// Three existing state buffers are reused; extra face/boundary buffers
    /// are prepared once, with no per-stage field allocation.
    ///
    /// The enlarged linear stability interval does not prove nonlinear
    /// positivity or entropy decay. Physical failures halve the private trial;
    /// no state is clipped. Choose `controls.maximum_step_s` for temporal
    /// accuracy: this entry does not estimate time error. Euler time and its
    /// advection report are unchanged, also on a late refusal.
    pub fn advance_stabilized(
        &mut self,
        flow: &mut Euler2d,
        duration: Time,
        controls: MolecularControl,
        stabilization: MolecularStabilization,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, MolecularError> {
        self.advance_stabilized_with_budget(
            flow, duration, controls, stabilization, &mut 0, &mut cancelled,
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
                    allocated(self.work.x_flux.len(), NormalFaceFlux::default())?,
                    allocated(self.work.y_flux.len(), NormalFaceFlux::default())?,
                ],
                boundary_stage: allocated(self.work.boundary_step.len(), Totals::default())?,
            };
            admit_workspace(work.bytes(), maximum_bytes)?;
            self.stabilized = Some(work);
        }
        admit_workspace(self.stabilized_workspace_bytes(), maximum_bytes)
    }

    /// One cumulative counter is shared by composing owners, including retries.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_stabilized_with_budget(
        &mut self,
        flow: &mut Euler2d,
        duration: Time,
        controls: MolecularControl,
        stabilization: MolecularStabilization,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, MolecularError> {
        poll(cancelled)?;
        let selected = stabilization.controller(controls)?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid("molecular operator and flow differ"));
        }
        if !positive(duration.value()) {
            return Err(MolecularError::Invalid("molecular duration"));
        }
        self.prepare_stabilized(stabilization.maximum_extra_workspace_bytes)?;
        let before = checked_totals(&flow.states, self.grid.cell_volume(), cancelled)?;
        copy_values(&flow.states, &mut self.work.current, cancelled)?;
        for tile in self.work.boundary_trial.chunks_mut(256) {
            poll(cancelled)?;
            tile.fill(Totals::default());
        }
        let accepted = controller::advance(
            duration.value(), selected, evaluations, cancelled,
            |action, elapsed, evaluations, cancelled| {
                let rate = match action {
                    Action::First => {
                        return spatial(
                            &self.gas, &self.grid, self.boundaries, &self.walls,
                            &self.work.current, &mut self.law, &mut self.work.primitive,
                            &mut self.work.gradient, &mut self.work.viscosity,
                            &mut self.work.x_flux, &mut self.work.y_flux,
                            &mut self.work.binary, cancelled,
                        ).map(|rate| Stage { rate, ..Stage::default() });
                    }
                    Action::Trial { dt, scheme, first, cached } => self.rkl_trial(
                        dt, scheme, first.rate, cached, selected, elapsed, evaluations, cancelled,
                    )?,
                };
                for (index, (total, change)) in self.work.boundary_trial.iter_mut()
                    .zip(&self.work.boundary_step).enumerate()
                {
                    if index % 256 == 0 { poll(cancelled)?; }
                    *total = add_totals(*total, *change);
                }
                std::mem::swap(&mut self.work.current, &mut self.work.stage);
                Ok(Stage { rate, ..Stage::default() })
            },
        )?;
        let after = checked_totals(&self.work.current, self.grid.cell_volume(), cancelled)?;
        let outward = sum_arrays(self.work.boundary_trial.iter().map(|q| q.array()), cancelled)?;
        let boundary_change = Totals::from_array(outward.map(|v| -v));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(cancelled)?;
        flow.states.copy_from_slice(&self.work.current);
        self.boundary_impulses.copy_from_slice(&self.work.boundary_trial);
        Ok(StabilizedMolecularReport {
            transport: MolecularReport {
                elapsed: duration,
                accepted_steps: accepted.interval.accepted_steps,
                rejected_steps: accepted.interval.rejected_steps,
                maximum_parabolic_number: accepted.interval.maximum_parabolic_number,
                before, after, boundary_change, conservation_residual,
            },
            operator_evaluations: accepted.operator_evaluations,
            maximum_stages_used: accepted.maximum_stages_used,
            maximum_scaled_parabolic_number: accepted.maximum_scaled_parabolic_number,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn rkl_trial(
        &mut self,
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
            spatial(
                &self.gas, &self.grid, self.boundaries, &self.walls,
                &self.work.current, &mut self.law, &mut self.work.primitive,
                &mut self.work.gradient, &mut self.work.viscosity,
                &mut self.work.x_flux, &mut self.work.y_flux,
                &mut self.work.binary, cancelled,
            )?
        };
        let maximum = controls.transport.parabolic_cfl * scheme.stability_extent();
        parabolic(dt * first, maximum)?;
        let stable = self.stabilized.as_mut().ok_or(MolecularError::Invalid(
            "missing prepared planar stabilized workspace",
        ))?;
        copy_values(&self.work.x_flux, &mut stable.anchor_flux[0], cancelled)?;
        copy_values(&self.work.y_flux, &mut stable.anchor_flux[1], cancelled)?;
        update(
            &self.gas, &self.grid, &self.work.current, &self.work.current,
            &self.work.x_flux, &self.work.y_flux, &mut self.work.stage,
            dt * scheme.first_weight(), false, None, None, cancelled,
        )?;
        capture_boundaries(
            &self.grid, &self.work.x_flux, &self.work.y_flux,
            dt * scheme.first_weight(), &mut self.work.boundary_step, false, cancelled,
        )?;
        for tile in stable.boundary_stage.chunks_mut(256) {
            poll(cancelled)?;
            tile.fill(Totals::default());
        }
        let mut rate = first;
        for j in 2..=scheme.stage_count() {
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            let stage_rate = spatial(
                &self.gas, &self.grid, self.boundaries, &self.walls,
                &self.work.stage, &mut self.law, &mut self.work.primitive,
                &mut self.work.gradient, &mut self.work.viscosity,
                &mut self.work.x_flux, &mut self.work.y_flux,
                &mut self.work.binary, cancelled,
            )?;
            parabolic(dt * stage_rate, maximum)?;
            rate = rate.max(stage_rate);
            let [mu, nu, mu_tilde, gamma_tilde] = scheme.coefficients(j)?;
            let older = Cell::from_mut(self.work.candidate.as_mut_slice()).as_slice_of_cells();
            cartesian::update::<4, Error, _>(
                &self.grid, 4, dt, false,
                |_, _| 0.0,
                |k, d| {
                    let anchor = self.work.current[k].array()[d];
                    let old = if j == 2 { anchor } else { older[k].get().array()[d] };
                    anchor + mu * (self.work.stage[k].array()[d] - anchor) + nu * (old - anchor)
                },
                |axis, face, d| {
                    let current = if axis == 0 { self.work.x_flux[face] } else { self.work.y_flux[face] };
                    mu_tilde * flux_array(current, axis)[d]
                        + gamma_tilde * flux_array(stable.anchor_flux[axis][face], axis)[d]
                },
                cancelled,
                |k, row, cancelled| {
                    let q = State::from_array([row[0], row[1], row[2], row[3]]);
                    q.primitive(&self.gas, &mut *cancelled)?;
                    older[k].set(q);
                    Ok(())
                },
            )?;
            cartesian::capture::<Error>(
                &self.grid, 4, dt, false,
                |axis, face, d| {
                    let current = if axis == 0 { self.work.x_flux[face] } else { self.work.y_flux[face] };
                    mu_tilde * flux_array(current, axis)[d]
                        + gamma_tilde * flux_array(stable.anchor_flux[axis][face], axis)[d]
                },
                |face, d, value, _| {
                    let previous = self.work.boundary_step[face].array()[d];
                    let target = stable.boundary_stage[face].component_mut(d);
                    *target = mu * previous + nu * *target + value;
                },
                cancelled,
            )?;
            std::mem::swap(&mut self.work.stage, &mut self.work.candidate);
            std::mem::swap(&mut self.work.boundary_step, &mut stable.boundary_stage);
        }
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let endpoint_rate = spatial(
            &self.gas, &self.grid, self.boundaries, &self.walls,
            &self.work.stage, &mut self.law, &mut self.work.primitive,
            &mut self.work.gradient, &mut self.work.viscosity,
            &mut self.work.x_flux, &mut self.work.y_flux,
            &mut self.work.binary, cancelled,
        )?;
        parabolic(dt * endpoint_rate, maximum)?;
        Ok(rate.max(endpoint_rate))
    }
}
