//! Nonlinear conductor-frame response and complete finite-enthalpy feedback.
//! The native field supplies irreversible heat, terminal work and mechanical
//! work separately; only its original local dissipation is deposited.

use super::*;
use crate::circuit_enthalpy::{
    conduction::heat::ConductionHeatError,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
};
use fs_em::conduction::{
    motion::ConductionMotion,
    response::{motion::MovingConductionResponseReport, terminals::ConductionTerminalDrive},
};
use fs_exec::Cx;

/// Held voltage/current-density traces and prescribed material motion on one
/// fixed-geometry interval. Empty current traces insulate the total free current.
#[derive(Clone, Copy, Debug)]
pub struct MovingConductionDrive<'a> {
    pub boundary_voltage_v: &'a [(usize, f64)],
    pub current_traces: &'a [ConductionCurrentTrace],
    pub motion: ConductionMotion<'a>,
}

/// Original held controls for moving nonlinear conduction and finite heat.
/// Direct terminals enforce total current and equipotentiality on their
/// original exterior traces; zero-current electrodes may transfer current
/// locally. Source laws are evaluated at every actual trial field.
#[derive(Clone, Copy, Debug)]
pub enum MovingConductionResponseDrive<'a> {
    Field(MovingConductionDrive<'a>),
    Terminals {
        drive: ConductionTerminalDrive<'a>,
        motion: ConductionMotion<'a>,
    },
}
impl<'a> From<MovingConductionDrive<'a>> for MovingConductionResponseDrive<'a> {
    fn from(drive: MovingConductionDrive<'a>) -> Self {
        Self::Field(drive)
    }
}

/// The accepted native source/field state and its independent work accounts.
/// Positive mechanical work leaves the electrical/thermal system. Electrical
/// supply may be negative during generation; it is never deposited as heat.
#[derive(Clone, Debug)]
pub struct MovingConductionResponseThermalStep {
    pub start_s: f64,
    pub end_s: f64,
    /// Native fields and source closure, including ordered direct-electrode
    /// observations in `field.terminals` when terminal controls are used.
    pub field: MovingConductionResponseReport,
    pub contacts: Vec<contacts::ConductionContactHeat>,
    pub heat_j: f64,
    pub terminal_work_j: f64,
    pub mechanical_work_j: f64,
    /// Terminal work minus mechanical output, paired with Joule heat.
    pub net_supplied_work_j: f64,
}

struct ResponseCandidate {
    problem: PlanarConductionResponse,
    step: MovingConductionResponseThermalStep,
}
impl ThermalDriver for ResponseCandidate {
    type Model = PlanarConductionResponse;
    type Step = MovingConductionResponseThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.step.net_supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.problem, self.step)
    }
}

/// Accept the original enthalpy history and native source state together.
/// Failure or cancellation leaves both supplied owners reusable.
pub struct MovingConductionResponseEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    problem: PlanarConductionResponse,
    report: CircuitEnthalpyReport<S, MovingConductionResponseThermalStep>,
}
impl<'t, S> MovingConductionResponseEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn problem(&self) -> &PlanarConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, MovingConductionResponseThermalStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarConductionResponse,
        CircuitEnthalpyReport<S, MovingConductionResponseThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Close original nonlinear sources at T, phase, |E + u x B| and |J| on
    /// the same backward-Euler endpoint. Each native field solve consumes the
    /// shared circuit budget, including rejected trials and final admission.
    ///
    /// Motion and imposed B are prescribed on the original triangles. Material
    /// occupancy and geometry must remain fixed over the interval. This does
    /// not evolve induction, momentum, advection or a moving surface-work law.
    /// Contact-adjacent cells must have zero velocity, as enforced by the
    /// native owner. Each stationary contact retains its original heat split.
    /// Direct grouped voltage/current/floating/mixed electrodes use one native
    /// field solve per constitutive trial. Their independently authored source
    /// work and motional voltage offset remain in the accepted terminal report.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_response_candidate<'t, 'd>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: impl Into<MovingConductionResponseDrive<'d>>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_moving_response_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive.into(),
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_moving_response_resolved_candidate<'t, S, E: From<ConductionThermalError>>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: MovingConductionResponseDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MovingConductionResponseEnthalpyCandidate<'t, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(map_conduction_error)?;
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&response.maximum_iterations) {
            return Err(map_conduction_error(invalid(
                "invalid bounded moving response iteration limit",
            )));
        }
        let dt = self
            .admit_mesh_step(
                thermal,
                target_s,
                problem.mesh(),
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        let routes = self
            .admit_mesh_routes(thermal, problem.mesh(), cancel)
            .map_err(map_conduction_error)?;
        self.admit_response_caloric(problem, cancel)
            .map_err(map_conduction_error)?;
        let resolve = |state: &DepositionSimulation<'_>| {
            let coefficients = resolve_sources(state)?;
            self.admit_response_heat(sample(&coefficients.source), cancel)
                .map_err(map_conduction_error)?;
            Ok::<_, CircuitEnthalpyError<E>>(coefficients)
        };
        let accepted = resolve(thermal)?;
        let contact_axes = contacts::response_axes(problem.contact_response_sources());
        let contacts = contacts::PreparedContacts::new_for_response_contacts(
            self,
            problem.mesh(),
            problem.contacts(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            &contact_axes,
            cancel,
        )
        .map_err(map_conduction_error)?;
        let initial = self
            .response_materials_at(problem, sample(&accepted.source), &routes, cancel)
            .map_err(map_conduction_error)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(map_conduction_error(invalid(
                "accepted moving conductor sources must match the actual thermal state",
            )));
        }
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                check(cancel).map_err(map_conduction_error)?;
                if remaining == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = resolve(state)?;
                let sources = self
                    .response_materials_at(problem, sample(&coefficients.source), &routes, cancel)
                    .map_err(map_conduction_error)?;
                let contact_sources = contacts
                    .refresh_contacts(problem.contacts(), sample(&coefficients.source), cancel)
                    .map_err(map_conduction_error)?;
                let bounded = ConductionResponseControls {
                    maximum_iterations: response.maximum_iterations.min(remaining),
                    ..response
                };
                let candidate = match drive {
                    MovingConductionResponseDrive::Field(drive) => problem
                        .solve_with_motion_and_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            drive.boundary_voltage_v,
                            drive.current_traces,
                            drive.motion,
                            assembly,
                            field_controls,
                            bounded,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                    MovingConductionResponseDrive::Terminals { drive, motion } => problem
                        .solve_terminals_with_motion_and_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            drive,
                            motion,
                            assembly,
                            field_controls,
                            bounded,
                            remaining,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                }
                .map_err(|error| match error {
                    ConductionResponseError::FieldBudget { .. } => {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    ConductionResponseError::NotConverged { .. }
                        if bounded.maximum_iterations < response.maximum_iterations =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => map_conduction_error(response_error(other)),
                })?;
                let (next, field) = candidate.into_parts();
                let heat_j = finite(dt * field.field.power_w)?;
                let terminal_power_w = field
                    .terminals
                    .as_ref()
                    .map_or(field.motion.terminal_power_w, |terminals| {
                        terminals.supplied_power_w
                    });
                let terminal_work_j = finite(dt * terminal_power_w)?;
                let mechanical_work_j = finite(dt * field.motion.mechanical_power_w)?;
                let net_supplied_work_j = finite(terminal_work_j - mechanical_work_j)?;
                if (field.field.power_w > 0.0 && heat_j == 0.0)
                    || (heat_j - net_supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(map_conduction_error(invalid(
                        "moving Joule heat, terminal supply and mechanical output do not close",
                    )));
                }
                let circuit_trials = field.field.field_solves;
                let contact_heat = contacts
                    .heat(&field.field.contact_response)
                    .map_err(map_conduction_error)?;
                Ok(DrivenCoefficients {
                    circuit: ResponseCandidate {
                        problem: next,
                        step: MovingConductionResponseThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            contacts: contact_heat,
                            heat_j,
                            terminal_work_j,
                            mechanical_work_j,
                            net_supplied_work_j,
                        },
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers = self
                    .deposit_power(
                        &candidate.step.field.field.cell_power_w,
                        &routes,
                        checkpoint,
                    )
                    .map_err(map_conduction_error)?;
                contacts::deposit(&candidate.step.contacts, &mut powers, checkpoint)
                    .map_err(map_conduction_error)?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(map_conduction_error)?;
        Ok(MovingConductionResponseEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}

fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

/// Close nonlinear planar moving conduction through any complete heat assembly.
///
/// Every original thermal contact, conductive layer, reservoir and radiation
/// source sees the same proposed endpoint as the original electrical response.
/// Conductor-frame bulk loss and stationary electrical-contact loss enter their
/// finite receiving cells once; terminal work and Lorentz output stay separate.
/// The complete source sample and typed errors are retained. Direct grouped
/// voltage/current/floating/mixed electrodes and original boundary traces share
/// this endpoint; terminal current controls determine their actual nonlinear
/// voltage including motional EMF, without rescaling an accepted field.
///
/// The original fixed-occupancy, source-pin, confined-phase and field-independent
/// caloric/thermal-source requirements apply. Contact-adjacent triangles must
/// remain stationary. Graph electrical work cannot be silently omitted. Closed
/// enclosure radiation is internal; external non-electrical heat equals
/// report.non_electrical_heat_j - dt * report.source.as_ref().radiation_net_power_w().
/// Cancellation, source refusal or a failed energy gate accepts neither owner.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_moving_conduction_response_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &PlanarConductionResponse,
    drive: impl Into<MovingConductionResponseDrive<'d>>,
    contact_bindings: &[contacts::ConductionContactThermal<'_>],
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    MovingConductionResponseEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    if regions.len() != problem.sources().len() || regions.len() > assembly.max_cells {
        return Err(map_conduction_error(invalid(
            "one bounded heat route per original moving triangle required",
        )));
    }
    let mut routes = BTreeMap::new();
    let mut field_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(map_conduction_error)?;
        let source = problem.sources().get(region.triangle).ok_or_else(|| {
            map_conduction_error(invalid("heat route names an absent moving conductor"))
        })?;
        if routes.insert(region.triangle, region.cell).is_some() {
            return Err(map_conduction_error(invalid(
                "original moving conductor heat routes must be unique",
            )));
        }
        for (axis, _) in source.axes() {
            field_axes.push((region.cell, axis.as_str()));
        }
    }
    if let Some(sources) = problem.contact_response_sources() {
        let trace = problem.mesh().boundary_trace().map_err(|_| {
            map_conduction_error(invalid("moving contact lost its original exterior trace"))
        })?;
        let edge_key = |mut vertices: [usize; 2]| {
            vertices.sort_unstable();
            vertices
        };
        let mut adjacent = BTreeMap::new();
        for edge in &trace.edges {
            check(&cancel).map_err(map_conduction_error)?;
            adjacent.insert(
                edge_key(edge.oriented_vertices.map(|vertex| vertex as usize)),
                edge.source_face,
            );
        }
        for (law, contact) in sources.laws().iter().zip(problem.contacts()) {
            check(&cancel).map_err(map_conduction_error)?;
            for vertices in [contact.first_vertices, contact.second_vertices] {
                let cell = adjacent
                    .get(&edge_key(vertices))
                    .and_then(|face| routes.get(face))
                    .ok_or_else(|| {
                        map_conduction_error(invalid(
                            "moving nonlinear contact lacks its original heat route",
                        ))
                    })?;
                for (axis, _) in law.axes() {
                    field_axes.push((*cell, axis.as_str()));
                }
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &field_axes)
        .map_err(heat_error)?;
    let layout = ConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let next = layout.step_moving_response_resolved_candidate(
        thermal,
        target_s,
        problem,
        drive.into(),
        contact_bindings,
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "moving nonlinear current cannot omit graph electrical work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    let report = next.report();
    closed_balance(
        target_s - thermal.time_s(),
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}
