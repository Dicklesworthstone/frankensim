//! RKL2 stages on the original planar source, divergence and boundary owners.

use super::super::{
    DissipationError, MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler2d, MolecularReport,
    MolecularTransport2d, PlanarTransport, Stage, TransportLaw, WIDTH_LIMIT, allocated, capture,
    cartesian, cell, copy_values, parabolic, poll, spatial,
};
use super::{Action, Rkl2Scheme, StabilizedDissipationControl, advance, reserve_evaluation};
use fs_qty::Time;

/// Complete original planar balances and actual stabilized integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularReport {
    /// Original physical transfers and balances. Accepted steps are macrosteps;
    /// its maximum parabolic number is unscaled by the stability extension.
    pub transport: MolecularReport,
    /// Actual RHS evaluations, including rejected trials and endpoint admission.
    pub operator_evaluations: usize,
    /// Largest accepted RKL2 stage count; zero for a zero-duration interval.
    pub maximum_stages_used: usize,
    /// Largest accepted `dt * rate / ((s*s+s-2)/4)`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Combined field/operator numerical bytes required by the RKL2 path.
    ///
    /// This adds two anchor-face arrays and one boundary recurrence array to
    /// the original allocation plan. Source/provider exclusions and excess
    /// input state capacity follow [`Self::required_workspace_bytes`]. Pass the
    /// returned allowance to [`Self::new`] when preparing stabilized transport.
    pub fn required_stabilized_workspace_bytes(
        flow: &MixtureEuler2d,
        source: &PlanarTransport<L>,
    ) -> Result<usize, DissipationError> {
        let base = Self::required_workspace_bytes(flow, source)?;
        let [nx, ny] = flow.grid.shape;
        let checked = || -> Option<usize> {
            let faces = (nx + 1)
                .checked_mul(ny)?
                .checked_add(nx.checked_mul(ny + 1)?)?;
            let boundary = nx.checked_add(ny)?.checked_mul(2)?;
            let extra = faces
                .checked_add(boundary)?
                .checked_mul(flow.gas.species_count() + 3)?
                .checked_mul(size_of::<f64>())?;
            base.checked_add(extra)
        };
        checked().ok_or(DissipationError::Invalid(
            "stabilized planar molecular workspace size overflow",
        ))
    }

    /// Advance a complete molecular interval with the original bounded RKL2.
    ///
    /// Every intermediate state and the endpoint evaluate the same source law,
    /// Maxwell--Stefan, enthalpy, heat and full planar stress face owner. Held
    /// baths retain their actual half-cell distance and boundary energy effort.
    /// Signed face transfers use the same recurrence as the conservative field.
    /// The original three field arrays are reused; extra face/boundary scratch
    /// is reserved lazily under the combined numerical byte budget.
    ///
    /// `controls.transport.maximum_step_s` caps the physical macrostep for
    /// accuracy. This path has no temporal error estimator or unconditional
    /// nonlinear positivity guarantee. Source or stage refusal halves the whole
    /// macrostep without clipping. Cancellation or exhausted work preserves the
    /// complete accepted field, reports, face loads and Euler clock.
    pub fn advance_stabilized(
        &mut self,
        flow: &mut MixtureEuler2d,
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
        self.publish(report.transport);
        Ok(report)
    }

    // The composing owner carries the same work counter through both molecular
    // halves and every discarded complete source/chemistry/advection attempt.
    pub(in crate::compressible::mixture) fn advance_stabilized_unpublished(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        controls: StabilizedDissipationControl,
        maximum_workspace_bytes: usize,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        if duration.value() > 0.0 {
            self.reserve_stabilized_workspace(flow, maximum_workspace_bytes, cancelled)?;
        }
        let result = advance(
            duration.value(),
            controls,
            evaluations,
            cancelled,
            |action, elapsed, evaluations, cancelled| match action {
                Action::First => self.stabilized_anchor(flow, controls, cancelled),
                Action::Trial {
                    dt,
                    scheme,
                    first,
                    cached,
                } => {
                    let first = if cached {
                        first
                    } else {
                        reserve_evaluation(
                            evaluations,
                            controls.maximum_operator_evaluations,
                            elapsed,
                        )?;
                        self.stabilized_anchor(flow, controls, cancelled)?
                    };
                    let peak = self.stabilized_trial(
                        flow,
                        dt,
                        scheme,
                        first,
                        controls,
                        elapsed,
                        evaluations,
                        cancelled,
                    )?;
                    self.accept_stabilized_trial(cancelled)?;
                    Ok(peak)
                }
            },
        )?;
        let transport =
            self.complete_interval(flow, duration, before, result.interval, cancelled)?;
        Ok(StabilizedMolecularReport {
            transport,
            operator_evaluations: result.operator_evaluations,
            maximum_stages_used: result.maximum_stages_used,
            maximum_scaled_parabolic_number: result.maximum_scaled_parabolic_number,
        })
    }

    fn reserve_stabilized_workspace(
        &mut self,
        flow: &MixtureEuler2d,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let mut required = self.workspace_bytes().checked_add(flow.workspace_bytes());
        for (saved, count) in [
            (&self.work.stabilized_x_flux, self.work.x_flux.len()),
            (&self.work.stabilized_y_flux, self.work.y_flux.len()),
            (
                &self.work.stabilized_boundary,
                self.work.boundary_step.len(),
            ),
        ] {
            if saved.len() != count {
                required = required
                    .and_then(|bytes| count.checked_mul(size_of::<f64>())?.checked_add(bytes));
            }
        }
        let required = required.ok_or(DissipationError::Invalid(
            "stabilized planar molecular workspace size overflow",
        ))?;
        if required > maximum_workspace_bytes || required > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(DissipationError::Invalid(
                "stabilized planar molecular total workspace budget",
            ));
        }
        poll(cancelled)?;
        if self.work.stabilized_x_flux.len() != self.work.x_flux.len() {
            self.work.stabilized_x_flux = allocated(self.work.x_flux.len(), 0.0)?;
        }
        if self.work.stabilized_y_flux.len() != self.work.y_flux.len() {
            self.work.stabilized_y_flux = allocated(self.work.y_flux.len(), 0.0)?;
        }
        if self.work.stabilized_boundary.len() != self.work.boundary_step.len() {
            self.work.stabilized_boundary = allocated(self.work.boundary_step.len(), 0.0)?;
        }
        // Count real capacities as well as the admission plan; allocator excess
        // cannot disappear from subsequent source/flow composition accounting.
        let actual = self.workspace_bytes().checked_add(flow.workspace_bytes());
        if actual.is_none_or(|bytes| {
            bytes > maximum_workspace_bytes || bytes > MAX_MIXTURE_WORKSPACE_BYTES
        }) {
            return Err(DissipationError::Invalid(
                "stabilized planar molecular actual workspace budget",
            ));
        }
        poll(cancelled)
    }

    pub(in crate::compressible::mixture::molecular) fn stabilized_anchor(
        &mut self,
        flow: &MixtureEuler2d,
        controls: StabilizedDissipationControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        let first = spatial(
            &mut self.source,
            &self.grid,
            self.boundaries,
            &flow.reservoirs,
            &self.work.current,
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.x_flux,
            &mut self.work.y_flux,
            controls.transport,
            cancelled,
        )?;
        copy_values(
            &self.work.x_flux,
            &mut self.work.stabilized_x_flux,
            cancelled,
        )?;
        copy_values(
            &self.work.y_flux,
            &mut self.work.stabilized_y_flux,
            cancelled,
        )?;
        Ok(first)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::mixture::molecular) fn stabilized_trial(
        &mut self,
        flow: &MixtureEuler2d,
        dt: f64,
        scheme: Rkl2Scheme,
        first: Stage,
        controls: StabilizedDissipationControl,
        elapsed: f64,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        let limit = controls.transport.parabolic_cfl * scheme.stability_extent();
        parabolic(dt * first.rate, limit)?;
        super::super::update(
            self.source.gas(),
            &self.grid,
            &self.work.current,
            &self.work.current,
            &self.work.stabilized_x_flux,
            &self.work.stabilized_y_flux,
            &mut self.work.stage,
            dt * scheme.first_weight(),
            false,
            None,
            cancelled,
        )?;
        copy_values(&self.work.current, &mut self.work.candidate, cancelled)?;
        capture(
            &self.grid,
            self.gas().species_count(),
            &self.work.stabilized_x_flux,
            &self.work.stabilized_y_flux,
            dt * scheme.first_weight(),
            false,
            &mut self.work.boundary_step,
            cancelled,
        )?;
        for tile in self.work.stabilized_boundary.chunks_mut(4096) {
            poll(cancelled)?;
            tile.fill(0.0);
        }
        let mut peak = first;
        for stage in 2..=scheme.stage_count() + 1 {
            poll(cancelled)?;
            reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
            let evaluated = spatial(
                &mut self.source,
                &self.grid,
                self.boundaries,
                &flow.reservoirs,
                &self.work.stage,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.x_flux,
                &mut self.work.y_flux,
                controls.transport,
                cancelled,
            )?;
            parabolic(dt * evaluated.rate, limit)?;
            peak = peak.merge(evaluated);
            if stage > scheme.stage_count() {
                break; // final state passed full original source/face admission
            }
            self.stabilized_recurrence(dt, scheme.coefficients(stage)?, cancelled)?;
            core::mem::swap(&mut self.work.stage, &mut self.work.candidate);
        }
        for tile in self.work.boundary_step.chunks(4096) {
            poll(cancelled)?;
            if tile.iter().any(|value| !value.is_finite()) {
                return Err(DissipationError::Invalid(
                    "stabilized planar molecular boundary exchange overflow",
                ));
            }
        }
        Ok(peak)
    }

    fn accept_stabilized_trial(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
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
        core::mem::swap(&mut self.work.current, &mut self.work.stage);
        Ok(())
    }

    fn stabilized_recurrence(
        &mut self,
        dt: f64,
        coefficients: [f64; 4],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let [mu, nu, mu_tilde, gamma_tilde] = coefficients;
        let gas = self.source.gas();
        let species = gas.species_count();
        let width = species + 3;
        let work = &mut self.work;
        // Cartesian update reads all components of a single cell before its
        // write callback. Interior mutability therefore reuses the older-stage
        // array safely without adding a field or duplicating the divergence.
        let destination =
            core::cell::Cell::from_mut(work.candidate.as_mut_slice()).as_slice_of_cells();
        let flux = |axis, face, component| {
            let offset = face * width + component;
            if axis == 0 {
                mu_tilde * work.x_flux[offset] + gamma_tilde * work.stabilized_x_flux[offset]
            } else {
                mu_tilde * work.y_flux[offset] + gamma_tilde * work.stabilized_y_flux[offset]
            }
        };
        cartesian::update::<WIDTH_LIMIT, DissipationError, _>(
            &self.grid,
            width,
            dt,
            false,
            |index, component| work.current[index * width + component],
            |index, component| {
                let offset = index * width + component;
                let anchor = work.current[offset];
                anchor
                    + mu * (work.stage[offset] - anchor)
                    + nu * (destination[offset].get() - anchor)
            },
            &flux,
            cancelled,
            |index, row, cancelled| {
                cell(row, species).primitive(gas, &mut *cancelled)?;
                for (target, value) in destination[index * width..(index + 1) * width]
                    .iter()
                    .zip(row)
                {
                    target.set(*value);
                }
                Ok(())
            },
        )?;
        cartesian::capture::<DissipationError>(
            &self.grid,
            width,
            dt,
            false,
            &flux,
            |face, component, value, _| {
                let offset = face * width + component;
                work.stabilized_boundary[offset] =
                    mu * work.boundary_step[offset] + nu * work.stabilized_boundary[offset] + value;
            },
            cancelled,
        )?;
        core::mem::swap(&mut work.boundary_step, &mut work.stabilized_boundary);
        Ok(())
    }
}
