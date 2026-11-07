//! Physical local-error control on the original volumetric RKL2 operator.
//!
//! The shared molecular controller selects bounded full/two-half comparisons.
//! Every probe retains the native spatial source, held reservoirs, thermal
//! walls, conservative divergence and signed boundary quadrature. Only the
//! admitted two-half solution is accepted; no extrapolated state is created.

use super::*;
pub use crate::compressible::mixture::molecular::adaptive::AdaptiveDissipationControl;
pub(super) use crate::compressible::mixture::molecular::adaptive::Workspace;
use crate::compressible::mixture::molecular::adaptive::{
    self as controller, PhysicalState, face_error, state_error,
};
use crate::compressible::mixture::molecular::stabilized::{Action, Rkl2Scheme};

/// Complete original physical balances, actual probe work and local accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMolecularReport {
    /// Native RKL2 balances and all RHS calls, including rejected comparisons.
    /// Accepted steps count full/two-half comparisons; each publishes two halves.
    pub stabilized: StabilizedMolecularReport,
    /// Maximum accepted component/face error ratio, at most one.
    pub maximum_estimated_local_error: f64,
    /// Complete probes rejected by the local temporal-error comparison.
    pub temporal_rejections: usize,
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Complete flow/operator storage including all lazy RKL2 and error buffers.
    /// Separately configured SGS and finite-wall reserves retain their charges.
    pub fn required_adaptive_workspace_bytes(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_stabilized_workspace_bytes(flow, source)?
            .checked_add(Workspace::required_bytes(
                flow.states.len(),
                flow.grid.boundary_len() * (flow.gas().species_count() + 4),
            )?)
            .ok_or(DissipationError::Invalid(
                "volumetric adaptive workspace size overflow",
            ))
    }

    /// Evolve a complete molecular interval to declared physical local accuracy.
    ///
    /// One full step and two half steps use the SAME Legendre stage count.
    /// Species densities, temperature, all three velocity components and every
    /// physical face's signed transfer are compared separately. These estimated
    /// local errors do not certify spatial, splitting or whole-trajectory error.
    ///
    /// Reservoir traces, wall temperatures and prescribed wall velocities stay
    /// fixed throughout the interval. Finite wall inertia retains its native
    /// stabilized-path refusal. Euler time and advection reports are unchanged.
    /// All four additional scratch arrays are retained and reused across calls;
    /// no per-probe whole-field allocation, clipping or flux correction occurs.
    pub fn advance_adaptive(
        &mut self,
        flow: &mut MixtureEuler3d,
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
        self.publish(report.stabilized.transport, None);
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::mixture) fn advance_adaptive_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: AdaptiveDissipationControl,
        maximum_workspace_bytes: usize,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AdaptiveMolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        self.admit_stabilized()?;
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        self.prepare_adaptive(
            flow,
            duration.value() > 0.0,
            maximum_workspace_bytes,
            cancelled,
        )?;
        let accepted = controller::advance(
            duration.value(),
            controls,
            evaluations,
            cancelled,
            |action, elapsed, evaluations, cancelled| match action {
                Action::First => spatial(
                    &mut self.source,
                    &self.grid,
                    self.boundaries,
                    self.pressure_outlets,
                    &flow.reservoirs,
                    &self.walls,
                    &self.work.current,
                    &mut self.work.primitive,
                    &mut self.work.gradient,
                    &mut self.work.viscosity,
                    &mut self.work.flux,
                    #[cfg(feature = "experimental-sgs")]
                    self.sgs.as_mut(),
                    controls.stabilized.transport,
                    cancelled,
                )
                .map(|stage| (stage, 0.0)),
                Action::Trial {
                    dt,
                    scheme,
                    first,
                    cached,
                } => self.adaptive_trial(
                    &flow.reservoirs,
                    dt,
                    scheme,
                    first,
                    cached,
                    controls,
                    elapsed,
                    evaluations,
                    cancelled,
                ),
            },
        )?;
        let transport = self.complete_interval(
            flow,
            duration,
            before,
            accepted.stabilized.interval,
            None,
            cancelled,
        )?;
        Ok(AdaptiveMolecularReport {
            stabilized: StabilizedMolecularReport {
                transport,
                operator_evaluations: accepted.stabilized.operator_evaluations,
                maximum_stages_used: accepted.stabilized.maximum_stages_used,
                maximum_scaled_parabolic_number: accepted
                    .stabilized
                    .maximum_scaled_parabolic_number,
            },
            maximum_estimated_local_error: accepted.maximum_estimated_local_error,
            temporal_rejections: accepted.temporal_rejections,
        })
    }

    fn prepare_adaptive(
        &mut self,
        flow: &MixtureEuler3d,
        needed: bool,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        if !needed {
            return Ok(());
        }
        let boundary = self.grid.boundary_len() * (self.gas().species_count() + 4);
        let extra = if let Some(work) = &self.adaptive {
            work.additional_bytes(flow.states.len(), boundary)?
        } else {
            Workspace::required_bytes(flow.states.len(), boundary)?
        };
        let stable_extra = self.additional_stabilized_workspace_bytes()?;
        let required = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .and_then(|bytes| bytes.checked_add(extra))
            .and_then(|bytes| bytes.checked_add(stable_extra))
            .ok_or(DissipationError::Invalid(
                "volumetric adaptive workspace size overflow",
            ))?;
        // Admit the complete pending assembly before preparing either buffer set.
        admit_workspace(required, maximum_workspace_bytes)?;
        self.prepare_stabilized(flow, true, maximum_workspace_bytes)?;
        if self.adaptive.is_none() {
            let mut work = Workspace::default();
            work.reserve(flow.states.len(), boundary, cancelled)?;
            let actual = self
                .workspace_bytes()
                .checked_add(flow.workspace_bytes())
                .and_then(|bytes| bytes.checked_add(work.workspace_bytes()))
                .ok_or(DissipationError::Invalid(
                    "volumetric adaptive workspace size overflow",
                ))?;
            admit_workspace(actual, maximum_workspace_bytes)?;
            self.adaptive = Some(work);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn adaptive_trial(
        &mut self,
        reservoirs: &Reservoirs,
        dt: f64,
        scheme: Rkl2Scheme,
        first: Stage,
        cached: bool,
        controls: AdaptiveDissipationControl,
        elapsed: f64,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Stage, f64), DissipationError> {
        let work = self.adaptive.as_mut().ok_or(DissipationError::Invalid(
            "missing prepared volumetric adaptive workspace",
        ))?;
        copy_values(&self.work.current, &mut work.anchor, cancelled)?;
        let result = self.adaptive_probe(
            reservoirs,
            dt,
            scheme,
            first,
            cached,
            controls,
            elapsed,
            evaluations,
            cancelled,
        );
        if let Ok((stage, error)) = result {
            if (0.0..=1.0).contains(&error) {
                let work = self.adaptive.as_ref().ok_or(DissipationError::Invalid(
                    "missing prepared volumetric adaptive workspace",
                ))?;
                for (index, (((total, correction), first), second)) in self
                    .work
                    .boundary_trial
                    .iter_mut()
                    .zip(&mut self.work.boundary_correction)
                    .zip(&work.first_boundary)
                    .zip(&self.work.boundary_step)
                    .enumerate()
                {
                    if index % 4096 == 0 {
                        poll(cancelled)?;
                    }
                    let adjusted = (first + second) - *correction;
                    let next = *total + adjusted;
                    *correction = (next - *total) - adjusted;
                    *total = next;
                }
                std::mem::swap(&mut self.work.current, &mut self.work.stage);
                return Ok((stage, error));
            }
        }
        // A refused fine probe may have replaced current with its first half.
        // Every retry must restart at the same complete accepted anchor.
        let work = self.adaptive.as_ref().ok_or(DissipationError::Invalid(
            "missing prepared volumetric adaptive workspace",
        ))?;
        copy_values(&work.anchor, &mut self.work.current, cancelled)?;
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn adaptive_probe(
        &mut self,
        reservoirs: &Reservoirs,
        dt: f64,
        scheme: Rkl2Scheme,
        first: Stage,
        cached: bool,
        controls: AdaptiveDissipationControl,
        elapsed: f64,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Stage, f64), DissipationError> {
        let coarse = self.rkl_trial(
            reservoirs,
            dt,
            scheme,
            first,
            cached,
            controls.stabilized,
            elapsed,
            evaluations,
            cancelled,
        )?;
        let work = self.adaptive.as_mut().ok_or(DissipationError::Invalid(
            "missing prepared volumetric adaptive workspace",
        ))?;
        let species = self.source.gas().species_count();
        let width = species + 4;
        for (index, (row, saved)) in self
            .work
            .stage
            .chunks_exact(width)
            .zip(work.coarse.chunks_exact_mut(width))
            .enumerate()
        {
            poll(cancelled)?;
            saved[..species].copy_from_slice(&row[..species]);
            saved[species..species + 3].copy_from_slice(&self.work.primitive[index].velocity);
            saved[species + 3] = self.work.primitive[index].temperature;
        }
        copy_values(
            &self.work.boundary_step,
            &mut work.coarse_boundary,
            cancelled,
        )?;
        // The full probe retained its original anchor faces. Reuse those exact
        // source evaluations for the first half instead of solving MS again.
        let stable = self.stabilized.as_ref().ok_or(DissipationError::Invalid(
            "missing prepared volumetric stabilized workspace",
        ))?;
        for (source, target) in stable.anchor_flux.iter().zip(&mut self.work.flux) {
            copy_values(source, target, cancelled)?;
        }
        let first_half = self.rkl_trial(
            reservoirs,
            0.5 * dt,
            scheme,
            first,
            true,
            controls.stabilized,
            elapsed,
            evaluations,
            cancelled,
        )?;
        let work = self.adaptive.as_mut().ok_or(DissipationError::Invalid(
            "missing prepared volumetric adaptive workspace",
        ))?;
        copy_values(
            &self.work.boundary_step,
            &mut work.first_boundary,
            cancelled,
        )?;
        std::mem::swap(&mut self.work.current, &mut self.work.stage);
        // The first half's qualified endpoint is this half's exact anchor;
        // its original face fluxes and source/SGS scratch are still prepared.
        let second_half = self.rkl_trial(
            reservoirs,
            0.5 * dt,
            scheme,
            first_half,
            true,
            controls.stabilized,
            elapsed,
            evaluations,
            cancelled,
        )?;
        let mut stats = coarse.merge(first_half).merge(second_half);
        stats.rate = coarse
            .rate
            .max(0.5 * first_half.rate)
            .max(0.5 * second_half.rate);
        Ok((stats, self.adaptive_error(controls, cancelled)?))
    }

    fn adaptive_error(
        &self,
        controls: AdaptiveDissipationControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DissipationError> {
        let work = self.adaptive.as_ref().ok_or(DissipationError::Invalid(
            "missing prepared volumetric adaptive workspace",
        ))?;
        let gas = self.gas();
        let species = gas.species_count();
        let width = species + 4;
        let mut error = 0.0_f64;
        for (index, (coarse, fine)) in work
            .coarse
            .chunks_exact(width)
            .zip(self.work.stage.chunks_exact(width))
            .enumerate()
        {
            poll(cancelled)?;
            let fine_state = self.work.primitive[index];
            error = error.max(state_error(
                gas,
                coarse,
                fine,
                PhysicalState {
                    density: coarse[..species].iter().sum(),
                    temperature: coarse[species + 3],
                    velocity: [coarse[species], coarse[species + 1], coarse[species + 2]],
                },
                PhysicalState {
                    density: fine_state.density,
                    temperature: fine_state.temperature,
                    velocity: fine_state.velocity,
                },
                controls,
            )?);
        }
        for side in BoundarySide::ALL {
            if physical_condition(self.boundaries, side).is_none() {
                continue;
            }
            let axis = side.axis();
            let [a, b] = tangents(axis);
            for (line, face) in self.grid.boundary_range(side).enumerate() {
                poll(cancelled)?;
                let mut coordinates = [0; 3];
                coordinates[axis] = if side.sign() < 0.0 {
                    0
                } else {
                    self.grid.shape[axis] - 1
                };
                coordinates[a] = line % self.grid.shape[a];
                coordinates[b] = line / self.grid.shape[a];
                let index = self.grid.index(coordinates);
                let row = face * width..(face + 1) * width;
                error = error.max(face_error::<3>(
                    gas,
                    &work.coarse_boundary[row.clone()],
                    &work.first_boundary[row.clone()],
                    &self.work.boundary_step[row],
                    &self.work.stage[index * width..index * width + species],
                    self.work.primitive[index].temperature,
                    self.grid.cell_volume(),
                    controls,
                    cancelled,
                )?);
            }
        }
        Ok(error)
    }
}
