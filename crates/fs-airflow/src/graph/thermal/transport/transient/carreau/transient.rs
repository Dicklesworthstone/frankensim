//! Retained original radial momentum with regional native enthalpy transport.
//!
//! Every thermal endpoint trial restarts the SAME accepted radial histories.
//! Original local viscosity is held at each trial endpoint temperature during
//! the native midpoint momentum interval. Its midpoint quadrature of mean flow
//! drives native upwind enthalpy transport, and ONLY its integrated viscous heat
//! is deposited locally. This endpoint-thermal/midpoint-momentum composition is
//! first order in coupled time; it is not a midpoint-temperature source model.
//!
//! Fixed regional masses cannot absorb hydraulic compliance. This owner admits
//! only rigid incompressible pressure junctions; the separate one-inventory
//! thermal owner remains available for explicitly compliant networks.

use std::borrow::Borrow;
use std::cell::{Cell, RefCell};

/// Native time refinement that resolves radial history and reversing transport.
pub mod adaptive;
/// Original pump/drive losses deposited at explicitly authored native cells.
pub mod equipment;

use adaptive::RegionalCarreauDrive;
use equipment::EquipmentThermalReport;
use fs_exec::Cx;
use fs_qty::Pressure;
use fs_thermal::enthalpy::DepositionSimulation;

use super::super::nonlinear::{
    ImplicitTransportControls, ImplicitTransportError, ImplicitTransportReport,
    advance_endpoint_transport,
};
use super::{CarreauTransportExchange, RegionalCarreauTransport, thermal_checkpoint};
use crate::graph::rheology::FlowInjection;
use crate::graph::rheology::transient::head::{FlowHeadControls, FlowHeadError, FlowHeadSource};
use crate::graph::rheology::transient::{
    TransientCarreauControls, TransientCarreauError, TransientCarreauGraph, TransientCarreauState,
    TransientCarreauStep, TransientCarreauWorkspace,
};
use crate::graph::{FixedPressure, Sum};

/// Original rigid radial owner and the unchanged finite thermal inventory.
#[derive(Clone)]
pub struct TransientRegionalCarreauTransport {
    momentum: TransientCarreauGraph,
    regional: RegionalCarreauTransport,
}

/// Native accuracy controls and aggregate nonlinear work/energy bounds.
#[derive(Clone, Copy, Debug)]
pub struct TransientRegionalCarreauControls {
    /// The head source-query limit also bounds the aggregate of all thermal
    /// residuals, derivative probes and final acceptance for this interval.
    pub hydraulics: TransientCarreauControls,
    pub thermal: ImplicitTransportControls,
    /// Counts every native branch solve across all thermal residuals, numerical
    /// derivative probes and final acceptance, including rejected pressure trials.
    pub maximum_total_branch_solves: u64,
    /// Independent bound on pressure/head work minus kinetic change and native
    /// viscous heat. Numerical mechanical residuals never become thermal heat.
    pub mechanical_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum TransientCarreauTransportError<E> {
    Momentum(TransientCarreauError),
    /// The authored drive source refused a prepared thermal endpoint.
    Drive(E),
    Exchange(E),
    MechanicalEnergy {
        residual_j: f64,
        tolerance_j: f64,
    },
    CombinedEnergy {
        residual_j: f64,
        tolerance_j: f64,
    },
    EquipmentEnergy {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for TransientCarreauTransportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient regional Carreau transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for TransientCarreauTransportError<E> {}

/// Original work and heat accounts from one jointly accepted endpoint.
#[derive(Debug)]
pub struct TransientRegionalCarreauReport<T> {
    /// Original annular velocities, mean-flow quadrature, kinetic change,
    /// viscous heat, pressure/head work and numerical residuals stay distinct.
    pub hydraulics: TransientCarreauStep,
    /// Native stored energy, external heat and independent local/global gates.
    /// Its heat payload is the caller's exchange result at this same endpoint.
    pub thermal: ImplicitTransportReport<T>,
    pub external_stream_enthalpy_j: f64,
    /// Native ambient, reservoir and authored/callback source heat, evaluated
    /// directly before adding equipment losses, advection or viscous heat;
    /// internal links cancel. Automatically routed losses are never relabeled
    /// as caller-authored heat.
    pub other_heat_j: f64,
    pub equipment: EquipmentThermalReport,
    /// Delta(thermal + kinetic) - reservoir - injection - head work
    /// - equipment heat - other heat - external stream enthalpy.
    /// No residual is deposited.
    pub combined_energy_residual_j: f64,
    /// Native thermal and mechanical tolerances plus arithmetic roundoff.
    pub combined_energy_tolerance_j: f64,
    pub total_branch_solves: u64,
    /// Original equipment queries across all completed thermal trials. A
    /// terminal failed trial remains bounded by its clipped remaining budget.
    pub total_head_evaluations: u64,
}

/// All physical owners advance together, or no replacement is returned.
pub struct TransientRegionalCarreauStep<'m, T> {
    pub next: TransientRegionalCarreauTransport,
    pub momentum: TransientCarreauState,
    pub thermal: DepositionSimulation<'m>,
    pub report: TransientRegionalCarreauReport<T>,
}

struct Prepared<T, D> {
    drive: D,
    source: TransientCarreauGraph,
    state: TransientCarreauState,
    hydraulics: TransientCarreauStep,
    exchange: T,
    external_stream_enthalpy_j: f64,
    other_heat_j: f64,
    equipment: EquipmentThermalReport,
}

impl TransientRegionalCarreauTransport {
    /// Compose the already admitted regional inventory with its original
    /// transient pipe discretization. No hydraulic node may carry compliance.
    /// Original topology, pipe names/orientations, source laws, geometry, rate
    /// ceilings and reference density must agree. Existing radial histories
    /// and workspaces keep their original compatibility identities.
    pub fn new(
        cx: &Cx<'_>,
        regional: RegionalCarreauTransport,
        momentum: TransientCarreauGraph,
    ) -> Result<Self, ImplicitTransportError<TransientCarreauError>> {
        thermal_checkpoint(cx)?;
        if momentum.storage().iter().any(Option::is_some) {
            return Err(ImplicitTransportError::Invalid(
                "fixed regional thermal masses require zero hydraulic compliance",
            ));
        }
        if momentum.graph().node_count() != regional.graph().node_count()
            || momentum.pipes().len() != regional.pipe_cells().len()
            || momentum
                .pipes()
                .iter()
                .any(|pipe| pipe.density_kg_m3() != regional.reference_density().value())
        {
            return Err(ImplicitTransportError::Invalid(
                "regional and radial topology or reference density disagree",
            ));
        }
        let temperatures: Vec<_> = regional
            .graph()
            .branches()
            .iter()
            .map(|branch| branch.pipe.temperature_k())
            .collect();
        let momentum = momentum
            .at_temperatures(&temperatures, cx)
            .map_err(ImplicitTransportError::Resolve)?;
        if momentum.graph().branches() != regional.graph().branches() {
            return Err(ImplicitTransportError::Invalid(
                "regional and original radial pipe sources or geometry disagree",
            ));
        }
        thermal_checkpoint(cx)?;
        Ok(Self { momentum, regional })
    }

    pub fn momentum(&self) -> &TransientCarreauGraph {
        &self.momentum
    }
    pub fn regional(&self) -> &RegionalCarreauTransport {
        &self.regional
    }

    /// Advance the original radial histories and native regional heat together.
    /// Initial momentum and thermal clocks must match exactly. The original
    /// pressure/injection data are held for this interval. Zero-compliance
    /// nodes retain native initial, mean and final flow-continuity admission.
    ///
    /// Each thermal trial freezes its OWN endpoint pipe temperatures, clones
    /// the SAME accepted momentum history, then runs native midpoint momentum.
    /// The exchange callback is pure and receives actual native mean/old/new
    /// flows and source heat for that candidate. It supplies explicit balanced
    /// streams and nonadvective heat only. No pressure datum supplies enthalpy.
    ///
    /// The mean rate is the original midpoint quadrature, not an exact integral
    /// of the continuous trajectory. Native upwinding uses that signed mean and
    /// endpoint donor enthalpy. A reversal within the interval is not located;
    /// refine time to resolve counterflow and residence. This is first order in
    /// coupled time and retains the developed radial model's omission of axial
    /// momentum/kinetic-energy convection and junction kinetic losses.
    ///
    /// Native caloric, rate, work, energy and cancellation gates remain active.
    /// Source refusal or failed acceptance preserves BOTH physical inputs and
    /// permits retry; the caller's reusable workspace is scratch and may change.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        TransientRegionalCarreauStep<'m, T>,
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    > {
        self.advance_impl(
            cx,
            thermal,
            state,
            workspace,
            target_s,
            boundaries,
            injections,
            None,
            controls,
            energy_guess_j,
            exchange,
        )
    }

    /// Same atomic step with original held branch pressure rises. Recovering
    /// head work remains signed mechanical exchange, never negative viscous heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heads<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        branch_heads: &[Pressure],
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        TransientRegionalCarreauStep<'m, T>,
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    > {
        self.advance_impl(
            cx,
            thermal,
            state,
            workspace,
            target_s,
            boundaries,
            injections,
            Some(branch_heads),
            controls,
            energy_guess_j,
            exchange,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_impl<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        branch_heads: Option<&[Pressure]>,
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        TransientRegionalCarreauStep<'m, T>,
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    > {
        let drive = RegionalCarreauDrive {
            boundaries: boundaries.to_vec(),
            injections: injections.to_vec(),
            branch_heads: branch_heads.map(<[Pressure]>::to_vec),
            branch_head_laws: None,
            branch_heat_destinations: None,
        };
        self.advance_with_resolved_drive(
            cx,
            thermal,
            state,
            workspace,
            target_s,
            controls,
            energy_guess_j,
            |_, _, _| Ok::<_, E>(&drive),
            |cx, candidate, flow, _| exchange(cx, candidate, flow),
        )
        .map(|(step, _)| step)
    }

    /// Advance with original flow-dependent equipment sources held over this
    /// interval. Every thermal trial solves the operating point against its own
    /// prepared pipe viscosity. The actual head work remains mechanical; only
    /// original radial viscous heat enters the caloric inventory. Use the
    /// adaptive thermal-drive callback when the equipment source itself must
    /// be prepared from the same candidate thermal endpoint. Sources declaring
    /// positive routed equipment losses require advance_with_drive and explicit
    /// heat recipients; this entry point never guesses those destinations.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_head_laws<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        held_heads: Option<&[Pressure]>,
        head_laws: &[Option<FlowHeadSource>],
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        TransientRegionalCarreauStep<'m, T>,
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    > {
        let drive = RegionalCarreauDrive {
            boundaries: boundaries.to_vec(),
            injections: injections.to_vec(),
            branch_heads: held_heads.map(<[Pressure]>::to_vec),
            branch_head_laws: Some(head_laws.to_vec()),
            branch_heat_destinations: None,
        };
        self.advance_with_resolved_drive(
            cx,
            thermal,
            state,
            workspace,
            target_s,
            controls,
            energy_guess_j,
            |_, _, _| Ok::<_, E>(&drive),
            |cx, candidate, flow, _| exchange(cx, candidate, flow),
        )
        .map(|(step, _)| step)
    }

    /// Advance with the complete original held drive, including explicit
    /// recipients for any original pump/drive routed losses. The same nonlinear
    /// endpoint solve includes those heat sources and local viscosity feedback.
    /// The callback supplies additional independent heat/streams only: it must
    /// not add the automatically routed equipment losses a second time.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_drive<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        drive: &RegionalCarreauDrive,
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        TransientRegionalCarreauStep<'m, T>,
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    > {
        self.advance_with_resolved_drive(
            cx,
            thermal,
            state,
            workspace,
            target_s,
            controls,
            energy_guess_j,
            |_, _, _| Ok::<_, E>(drive),
            |cx, candidate, flow, _| exchange(cx, candidate, flow),
        )
        .map(|(step, _)| step)
    }

    /// One original thermal residual owner for prescribed and endpoint-resolved
    /// drive laws. Returning the accepted drive together with the physical state
    /// keeps discarded nonlinear source samples out of the accepted trajectory.
    /// A borrowed drive preserves prescribed forcing without cloning each trial.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_with_resolved_drive<'m, T, E, D>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        target_s: f64,
        controls: TransientRegionalCarreauControls,
        energy_guess_j: &[f64],
        resolve_drive: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauGraph,
        ) -> Result<D, E>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &TransientCarreauStep,
            &D,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        (TransientRegionalCarreauStep<'m, T>, D),
        ImplicitTransportError<TransientCarreauTransportError<E>>,
    >
    where
        D: Borrow<RegionalCarreauDrive>,
    {
        thermal_checkpoint(cx)?;
        let duration_s = target_s - thermal.time_s();
        if state.accepted_time_s() != thermal.time_s()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || controls.maximum_total_branch_solves == 0
            || !controls.mechanical_energy_tolerance_j.is_finite()
            || controls.mechanical_energy_tolerance_j < 0.0
        {
            return Err(ImplicitTransportError::Invalid(
                "matching physical clocks, positive interval and finite transient regional budgets required",
            ));
        }
        self.regional.inventory.validate_storage(cx, thermal)?;
        let total_branch_solves = Cell::new(0_u64);
        let total_head_evaluations = Cell::new(0_u64);
        let workspace = RefCell::new(workspace);
        let (next_thermal, thermal_report) = advance_endpoint_transport(
            cx,
            thermal,
            target_s,
            controls.thermal,
            energy_guess_j,
            |cx, candidate| {
                self.regional.inventory.validate_storage(cx, candidate)?;
                let momentum_error = |error| {
                    ImplicitTransportError::Resolve(TransientCarreauTransportError::Momentum(error))
                };
                let remaining = controls.maximum_total_branch_solves - total_branch_solves.get();
                if remaining == 0 {
                    return Err(momentum_error(TransientCarreauError::BranchSolveBudget));
                }
                let mut temperatures = Vec::with_capacity(self.regional.pipe_cells.len());
                for &id in self.regional.pipe_cells() {
                    thermal_checkpoint(cx)?;
                    temperatures.push(candidate.temperature_k(id)?);
                }
                // Endpoint-held rheology is part of this first-order coupled
                // residual. It never resets or replaces the accepted velocities.
                let source = self
                    .momentum
                    .at_temperatures(&temperatures, cx)
                    .map_err(momentum_error)?;
                let drive = resolve_drive(cx, candidate, &source).map_err(|error| {
                    ImplicitTransportError::Resolve(TransientCarreauTransportError::Drive(error))
                })?;
                thermal_checkpoint(cx)?;
                let forcing = drive.borrow();
                if forcing
                    .branch_head_laws
                    .as_ref()
                    .is_some_and(|laws| laws.len() != source.pipes().len())
                {
                    return Err(ImplicitTransportError::Invalid(
                        "one optional original equipment source per regional branch required",
                    ));
                }
                let head_remaining = controls.hydraulics.head.maximum_source_evaluations
                    - total_head_evaluations.get();
                let first_head = forcing
                    .branch_head_laws
                    .as_ref()
                    .and_then(|laws| laws.iter().position(Option::is_some));
                if let Some(branch) = first_head {
                    if head_remaining == 0 {
                        return Err(momentum_error(TransientCarreauError::Head {
                            branch,
                            source: FlowHeadError::SourceBudget,
                        }));
                    }
                }
                let mut proposed = state.clone();
                let hydraulic_controls = TransientCarreauControls {
                    maximum_branch_solves: controls.hydraulics.maximum_branch_solves.min(remaining),
                    head: FlowHeadControls {
                        maximum_source_evaluations: head_remaining,
                        ..controls.hydraulics.head
                    },
                    ..controls.hydraulics
                };
                let hydraulics = {
                    let mut scratch = workspace.borrow_mut();
                    match (first_head, forcing.branch_heads.as_deref()) {
                        (Some(_), _) => source.step_with_head_laws(
                            &mut proposed,
                            &mut scratch,
                            &forcing.boundaries,
                            &forcing.injections,
                            forcing.branch_heads.as_deref(),
                            forcing
                                .branch_head_laws
                                .as_deref()
                                .expect("source roster admitted"),
                            duration_s,
                            hydraulic_controls,
                            cx,
                        ),
                        (None, Some(heads)) => source.step_with_heads(
                            &mut proposed,
                            &mut scratch,
                            &forcing.boundaries,
                            &forcing.injections,
                            heads,
                            duration_s,
                            hydraulic_controls,
                            cx,
                        ),
                        (None, None) => source.step(
                            &mut proposed,
                            &mut scratch,
                            &forcing.boundaries,
                            &forcing.injections,
                            duration_s,
                            hydraulic_controls,
                            cx,
                        ),
                    }
                    .map_err(momentum_error)?
                };
                if hydraulics.branch_solves > remaining {
                    return Err(momentum_error(TransientCarreauError::BranchSolveBudget));
                }
                total_branch_solves.set(total_branch_solves.get() + hydraulics.branch_solves);
                if hydraulics.head_evaluations > head_remaining {
                    return Err(ImplicitTransportError::Invalid(
                        "original equipment work exceeded remaining thermal budget",
                    ));
                }
                total_head_evaluations
                    .set(total_head_evaluations.get() + hydraulics.head_evaluations);
                if proposed.accepted_time_s() != target_s {
                    return Err(ImplicitTransportError::Invalid(
                        "native momentum and thermal endpoint clocks disagree",
                    ));
                }
                if hydraulics.physical_work_residual_j.abs()
                    > controls.mechanical_energy_tolerance_j
                {
                    return Err(ImplicitTransportError::Resolve(
                        TransientCarreauTransportError::MechanicalEnergy {
                            residual_j: hydraulics.physical_work_residual_j,
                            tolerance_j: controls.mechanical_energy_tolerance_j,
                        },
                    ));
                }
                let equipment = equipment::prepare(
                    cx,
                    candidate,
                    &hydraulics,
                    forcing.branch_heads.as_deref(),
                    forcing.branch_heat_destinations.as_deref(),
                )?;
                let mut response =
                    exchange(cx, candidate, &hydraulics, &drive).map_err(|error| {
                        ImplicitTransportError::Resolve(TransientCarreauTransportError::Exchange(
                            error,
                        ))
                    })?;
                thermal_checkpoint(cx)?;
                // Measure the original nonadvective external source before
                // local dissipation/stream enthalpy are installed. Internal
                // reciprocal heat links remain internal to the native timeline.
                let other_power_w = candidate.external_step_power_with_boundary_conductances(
                    &response.heat.bulk.sources_w,
                    &response.heat.bulk.conductances_w_k,
                    &response.heat.boundaries,
                    || thermal_checkpoint(cx),
                )?;
                let other_heat_j = finite(other_power_w * duration_s)?;
                // Preserve independent callback/native sources at the same
                // recipients. These powers come only from the retained original
                // evaluation; no second equipment query or loss inference.
                for (&cell, &power_w) in &equipment.sources_w {
                    thermal_checkpoint(cx)?;
                    let source = response.heat.bulk.sources_w.entry(cell).or_default();
                    *source = sum(&[*source, power_w])?;
                }
                let rates: Vec<_> = hydraulics
                    .branches
                    .iter()
                    .map(|branch| branch.midpoint_volume_flow_m3_s)
                    .collect();
                let heat_w: Vec<_> = hydraulics
                    .branches
                    .iter()
                    .map(|branch| branch.viscous_heat_j / duration_s)
                    .collect();
                self.regional.endpoint_coefficients(
                    cx,
                    candidate,
                    &rates,
                    heat_w,
                    response,
                    |payload, stream_w, _| Prepared {
                        drive,
                        source,
                        state: proposed,
                        hydraulics,
                        exchange: payload,
                        external_stream_enthalpy_j: stream_w * duration_s,
                        other_heat_j,
                        equipment: equipment.report,
                    },
                )
            },
        )?;
        let ImplicitTransportReport {
            heat: mut prepared,
            thermal: audit,
            maximum_cell_residual_j,
            endpoint_evaluations,
            nonlinear,
        } = thermal_report;
        let terms = [
            audit.stored_energy_change_j,
            prepared.hydraulics.kinetic_energy_change_j,
            -prepared.hydraulics.reservoir_work_j,
            -prepared.hydraulics.injection_work_j,
            -prepared.hydraulics.branch_head_work_j,
            -prepared.external_stream_enthalpy_j,
            -prepared.other_heat_j,
            -prepared.equipment.deposited_heat_j,
        ];
        let combined = sum(&terms)?;
        let scale = sum(&terms.map(f64::abs))?;
        let allowance = finite(
            controls.thermal.total_tolerance_j
                + controls.mechanical_energy_tolerance_j
                + 64.0 * f64::EPSILON * scale,
        )?;
        if combined.abs() > allowance {
            return Err(ImplicitTransportError::Resolve(
                TransientCarreauTransportError::CombinedEnergy {
                    residual_j: combined,
                    tolerance_j: allowance,
                },
            ));
        }
        // Independently replace ONLY sourced hydraulic work plus routed losses
        // with original electrical input minus original surroundings heat.
        // Unknown-power/held heads remain explicit mechanical inputs.
        let equipment_terms = [
            audit.stored_energy_change_j,
            prepared.hydraulics.kinetic_energy_change_j,
            -prepared.hydraulics.reservoir_work_j,
            -prepared.hydraulics.injection_work_j,
            -prepared.equipment.remaining_head_work_j,
            -prepared.equipment.electrical_input_j,
            prepared.equipment.surroundings_heat_j,
            -prepared.other_heat_j,
            -prepared.external_stream_enthalpy_j,
        ];
        let equipment_residual = sum(&equipment_terms)?;
        let equipment_tolerance = finite(
            controls.thermal.total_tolerance_j
                + controls.mechanical_energy_tolerance_j
                + 128.0 * f64::EPSILON * sum(&equipment_terms.map(f64::abs))?,
        )?;
        if equipment_residual.abs() > equipment_tolerance {
            return Err(ImplicitTransportError::Resolve(
                TransientCarreauTransportError::EquipmentEnergy {
                    residual_j: equipment_residual,
                    tolerance_j: equipment_tolerance,
                },
            ));
        }
        prepared.equipment.electrical_boundary_residual_j = equipment_residual;
        prepared.equipment.electrical_boundary_tolerance_j = equipment_tolerance;
        let mut next_regional = self.regional.clone();
        next_regional.graph = prepared.source.graph().clone();
        thermal_checkpoint(cx)?;
        Ok((
            TransientRegionalCarreauStep {
                next: Self {
                    momentum: prepared.source,
                    regional: next_regional,
                },
                momentum: prepared.state,
                thermal: next_thermal,
                report: TransientRegionalCarreauReport {
                    hydraulics: prepared.hydraulics,
                    thermal: ImplicitTransportReport {
                        heat: prepared.exchange,
                        thermal: audit,
                        maximum_cell_residual_j,
                        endpoint_evaluations,
                        nonlinear,
                    },
                    external_stream_enthalpy_j: prepared.external_stream_enthalpy_j,
                    other_heat_j: prepared.other_heat_j,
                    equipment: prepared.equipment,
                    combined_energy_residual_j: combined,
                    combined_energy_tolerance_j: allowance,
                    total_branch_solves: total_branch_solves.get(),
                    total_head_evaluations: total_head_evaluations.get(),
                },
            },
            prepared.drive,
        ))
    }
}

fn finite<E>(value: f64) -> Result<f64, ImplicitTransportError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ImplicitTransportError::Invalid(
            "unrepresentable transient regional Carreau energy",
        ))
    }
}
fn sum<E>(terms: &[f64]) -> Result<f64, ImplicitTransportError<E>> {
    let mut sum = Sum::default();
    let error = |_| {
        ImplicitTransportError::Invalid("unrepresentable transient regional Carreau energy sum")
    };
    for &term in terms {
        sum.add(term).map_err(error)?;
    }
    sum.value().map_err(error)
}

#[cfg(test)]
mod tests;
