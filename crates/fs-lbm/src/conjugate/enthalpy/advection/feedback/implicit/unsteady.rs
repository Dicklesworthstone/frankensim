//! Endpoint thermal/material closure with the ORIGINAL retained face momentum.
//! Every native trial starts from the same accepted velocity and pressure;
//! neither constitutive iteration nor a Newton probe advances physical time.

use super::super::rheology::{RheologyControls, RheologyReport};
use super::*;
use crate::conjugate::UnsteadyHeatingAudit;

/// Independent thermal residual and optional physical-rate closure controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitUnsteadyFlowControls {
    pub thermal: ImplicitFlowControls,
    /// Required exactly when `with_native_rheology` bound a source. This closes
    /// its viscosity/conductivity/Cp at the actual endpoint strain invariant.
    pub rheology: Option<RheologyControls>,
}

/// One accepted coupled interval; these accounts are observations, never new
/// sources. The returned owner retains the exact source epoch and native faces.
#[derive(Debug)]
pub struct ImplicitUnsteadyFlowReport {
    pub thermal: EnergyAudit,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    /// All native solves across Newton, rate closure and final acceptance.
    pub momentum_solves: usize,
    pub rheology: Option<RheologyReport>,
    /// Native face-mass kinetic storage, physical heat, signed actuator work
    /// and nonnegative backward-Euler numerical loss remain separate.
    pub mechanical: UnsteadyHeatingAudit,
    pub deposited_heat_j: f64,
    pub net_stream_enthalpy_j: f64,
    pub other_thermal_input_j: f64,
    /// Delta thermal + delta native kinetic + numerical loss - physical supply.
    /// The independently bounded pressure-constraint defect is not a heat source.
    pub combined_residual_j: f64,
}

struct Endpoint<P> {
    owner: P,
    rheology: Option<RheologyReport>,
    deposited_w: f64,
    stream_w: f64,
    other_w: f64,
}

impl<F: Clone, R, T: FlowTransportSnapshot> MaterialFlowFeedback<'_, F, R, (), T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Close one backward-Euler interval of retained momentum and total h.
    ///
    /// Select the existing unsteady Stokes physical-heating owner first. The
    /// original face mass, viscosity/stress, pressure, porous drag, Boussinesq
    /// and actuator owners solve native inertia from the SAME accepted history
    /// for every trial. Source temperature, latent fraction and optional native
    /// shear all belong to that trial's endpoint. The converged ENDPOINT fluxes
    /// advect h through the original donor graph; conduction uses the identical
    /// source samples. Caloric storage and local/global residual acceptance stay
    /// in BackwardEulerStep. Numerical BE loss never becomes material heat.
    ///
    /// Returns the accepted thermal state AND its retained feedback owner.
    /// Resume with those two together, including after changing a prescribed
    /// actuator. Source callbacks must be repeatable and cannot mutate shared
    /// history. Nonconvergence, source/budget refusal and cancellation publish
    /// neither state. Interior events refuse; target events remain pending.
    /// `applied()` now describes the endpoint used by this implicit interval;
    /// `audit()` remains the separate explicit-advection observer account.
    ///
    /// Closed wall/symmetry boundaries, fixed reference mass/grid and unchanged
    /// phase occupancy are required. Mobile mushy fractions and sourced porous
    /// coefficients can evolve; opening/locking faces needs a separate admitted
    /// event/kinetic handoff and refuses here. No mutable wrapped observer,
    /// moving skeleton, free surface, compressibility, finite actuator reservoir,
    /// viscoelastic memory or gravitational-potential store is introduced.
    /// The existing Reynolds/source bounds still apply. No explicit thermal CFL
    /// is imposed; nonlinear convergence is not a temporal/mesh accuracy claim.
    #[allow(clippy::type_complexity)]
    pub fn step_unsteady_implicit<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: ImplicitUnsteadyFlowControls,
        initial_energy_increments_j: &[f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'t>, Self, ImplicitUnsteadyFlowReport), ImplicitFlowError>
    {
        let checkpoint = RefCell::new(checkpoint);
        let check = || check_gate(self.gate, &mut *checkpoint.borrow_mut());
        check()?;
        controls
            .thermal
            .admit(thermal.cells().len(), initial_energy_increments_j)?;
        let heat = self
            .unsteady
            .as_ref()
            .and_then(MaterialMomentum::heating_controls)
            .ok_or(ImplicitFlowError::Invalid(
                "implicit inertia needs the original unsteady Stokes heating owner",
            ))?;
        if self.viscous_heating.is_some()
            || self.advected.is_none()
            || self.support_capture_endpoint_s.is_some()
            || self
                .rheology
                .as_ref()
                .is_some_and(|binding| binding.quasisteady.is_some())
            || controls.rheology.is_some() != self.rheology.is_some()
        {
            return Err(ImplicitFlowError::Invalid(
                "incompatible momentum history or missing/extraneous rheology controls",
            ));
        }
        if let Some(rheology) = controls.rheology {
            rheology.admit()?;
        }
        self.check_map(thermal, &mut || check())?;
        self.check_unsteady_start(thermal, &mut || check())?;
        let step = BackwardEulerStep::new(thermal, target_s, || check())?;
        let dt = step.duration_s();
        let work = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate_with_transport(
                x,
                |state| {
                    let (mut proposal, rheology) = self.resolve_unsteady_endpoint(
                        state,
                        thermal.time_s(),
                        dt,
                        controls.rheology,
                        &work,
                        &mut || check(),
                    )?;
                    let epoch = proposal.applied.as_ref().expect("resolved native endpoint");
                    let active = proposal.momentum_domain(&epoch.samples, &mut || check())?;
                    let domain = active.as_ref().unwrap_or(proposal.domain);
                    let (flows, boundaries) = mass_graph(
                        domain,
                        &epoch.flow.field,
                        proposal.density_kg_m3,
                        proposal.input,
                        &mut || check(),
                    )?;
                    let graph = Advected::new((), flows, proposal.input.max_connections)?
                        .with_boundaries(boundaries, proposal.input.max_connections)?;
                    let transport = graph.endpoint_transport(state, || check())?;
                    let mut sources_w = transport.heat.bulk.sources_w;
                    let stream_w = sources_w.values().try_fold(0.0, |s, &p| finite(s + p))?;
                    let power = &proposal
                        .unsteady_mechanical_step()
                        .ok_or(invalid("implicit inertia lacks endpoint physical power"))?
                        .power;
                    let mut deposited_w = 0.0;
                    for sample in &epoch.samples {
                        check()?;
                        let watts = power.cell_heat_w[sample.voxel];
                        deposited_w = finite(deposited_w + watts)?;
                        let source = sources_w.entry(sample.cell).or_insert(0.0);
                        *source = finite(*source + watts)?;
                    }
                    proposal.advected = Some(graph);
                    let conductances_w_k = proposal.conductances(state, &mut || check())?;
                    let other_w =
                        state.external_step_power_w(&BTreeMap::new(), &conductances_w_k, || {
                            check()
                        })?;
                    Ok::<_, MaterialFlowError<DepositionError>>(TransportStepCoefficients {
                        heat: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: Endpoint {
                                owner: proposal,
                                rheology,
                                deposited_w,
                                stream_w,
                                other_w,
                            },
                        }
                        .into(),
                        transfers: transport.transfers,
                    })
                },
                || check(),
            )
        };
        let AcceptedEndpoint {
            thermal: next,
            payload: endpoint,
            audit: thermal_audit,
            nonlinear,
            maximum_cell_residual_j,
        } = solve_endpoint(
            step.cell_ids().len(),
            controls.thermal,
            initial_energy_increments_j,
            evaluate,
            self.gate,
            &check,
        )?;
        let mut accepted = endpoint.owner;
        accepted
            .unsteady
            .as_mut()
            .expect("retained endpoint")
            .accept_heating_step(thermal.time_s(), dt, target_s)
            .map_err(MaterialFlowError::from)?;
        let mechanical = accepted
            .unsteady_heating_audit()
            .expect("accepted physical heat");
        let deposited_heat_j = finite(dt * endpoint.deposited_w)?;
        let net_stream_enthalpy_j = finite(dt * endpoint.stream_w)?;
        let other_thermal_input_j = finite(dt * endpoint.other_w)?;
        let combined_residual_j = finite(
            thermal_audit.stored_energy_change_j
                + finite(mechanical.final_kinetic_energy_j - mechanical.initial_kinetic_energy_j)?
                + mechanical.backward_euler_loss_j
                - net_stream_enthalpy_j
                - other_thermal_input_j
                - mechanical.supplied_wall_work_j
                - mechanical.supplied_body_work_j,
        )?;
        let power_allowance = finite(dt * heat.power_tolerance_w)?;
        // Native mechanical residual and pressure-constraint work each have
        // their own independent power bound. Neither defect is deposited.
        let allowance =
            finite(controls.thermal.total_energy_tolerance_j + finite(2.0 * power_allowance)?)?;
        if (deposited_heat_j - mechanical.generated_heat_j).abs() > power_allowance
            || combined_residual_j.abs() > allowance
        {
            return Err(invalid(
                "implicit thermal, native kinetic, physical heat and work do not close",
            )
            .into());
        }
        accepted.bound_time_s = target_s;
        accepted.audit = AdvectionAudit::default();
        accepted.solves = work.get();
        accepted.check_unsteady_start(&next, &mut || check())?;
        check()?;
        Ok((
            next,
            accepted,
            ImplicitUnsteadyFlowReport {
                thermal: thermal_audit,
                nonlinear,
                maximum_cell_residual_j,
                momentum_solves: work.get(),
                rheology: endpoint.rheology,
                mechanical,
                deposited_heat_j,
                net_stream_enthalpy_j,
                other_thermal_input_j,
                combined_residual_j,
            },
        ))
    }

    fn resolve_unsteady_endpoint(
        &self,
        state: &DepositionSimulation<'_>,
        start_s: f64,
        dt: f64,
        controls: Option<RheologyControls>,
        work: &Cell<usize>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Option<RheologyReport>), MaterialFlowError<DepositionError>> {
        let mut rates = self.endpoint_shear_rate_per_s().map(<[f64]>::to_vec);
        for iteration in 1..=controls.map_or(1, |c| c.max_iterations) {
            check_gate(self.gate, check)?;
            // Every trial restores ORIGINAL accepted native inertia, never the
            // preceding fixed-point result. Only source-rate coordinates relax.
            let mut proposal = self.clone();
            proposal.solves = work.get();
            proposal
                .unsteady
                .as_mut()
                .expect("admitted native inertia")
                .begin_advance()?;
            let resolved = proposal.solve_material_flow_epoch_at_rates(state, rates.take(), check);
            work.set(proposal.solves);
            let (epoch, reference) = resolved?;
            for (old, new) in self.endpoint.iter().zip(&epoch.samples) {
                check_gate(self.gate, check)?;
                if old.voxel != new.voxel
                    || old.cell != new.cell
                    || old.transport.mobile_transport().is_some()
                        != new.transport.mobile_transport().is_some()
                {
                    return Err(invalid("implicit inertia requires unchanged phase occupancy; split at the phase handoff").into());
                }
            }
            proposal.reference = reference;
            proposal.applied = Some(epoch);
            proposal.advance_resolved_momentum(start_s, dt, state.time_s(), check)?;
            let mut epoch = proposal.applied.take().expect("resolved epoch");
            epoch.flow = proposal
                .unsteady_flow()
                .expect("advanced native inertia")
                .clone();
            let Some(controls) = controls else {
                proposal.endpoint = epoch.samples.clone();
                proposal.applied = Some(epoch);
                proposal.check_porous_endpoint(check)?;
                return Ok((proposal, None));
            };
            let active = proposal.momentum_domain(&epoch.samples, check)?;
            let domain = active.as_ref().unwrap_or(proposal.domain);
            let fluid = proposal
                .reference
                .map(|_| proposal.reference_fluid())
                .transpose()?;
            let actual = proposal
                .unsteady
                .as_ref()
                .expect("advanced native inertia")
                .shear_rates(domain, fluid.as_ref(), proposal.config, proposal.gate)?;
            let (report, closed, observed) =
                proposal.observe_rheology(state, &epoch, actual, iteration, controls, check)?;
            if closed {
                proposal.endpoint = observed;
                proposal.applied = Some(epoch);
                proposal.check_porous_endpoint(check)?;
                return Ok((proposal, Some(report)));
            }
            let mut next = proposal
                .endpoint_shear_rate_per_s()
                .expect("bound native rheology")
                .to_vec();
            let queried = epoch
                .shear_rate_per_s
                .as_ref()
                .expect("queried native rheology");
            if controls.relaxation != 1.0 {
                for (actual, &query) in next.iter_mut().zip(queried) {
                    *actual = controls
                        .relaxation
                        .mul_add(*actual, (1.0 - controls.relaxation) * query);
                }
            }
            rates = Some(next);
        }
        Err(DepositionError::Budget("implicit unsteady constitutive iterations").into())
    }
}
