//! The native nonlinear spatial conductor inside the original heat residual.
//!
//! Temperature comes from the same retained enthalpy chart, E/J from the same
//! FEEC field, and heat from its actual local E.J integrals. There is no second
//! thermal timeline, effective uniform resistance, or iterative-error heating.

use super::*;
use fs_em::conduction::response::{
    ConductionResponseControls, ConductionResponseError, ConductionResponseReport,
    MAX_RESPONSE_ITERATIONS, PlanarConductionResponse,
};
use fs_material::conductor::response::ConductorMagnitude;
use fs_qty::{Dims, QuantitySpec};

fn response_error(error: ConductionResponseError) -> Error {
    match error {
        ConductionResponseError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Response(other)),
    }
}

impl From<ConductionResponseError> for Error {
    fn from(error: ConductionResponseError) -> Self {
        response_error(error)
    }
}

/// Actual native fields and endpoint material evaluations, with independent
/// terminal work. Applied and evaluated source defects stay in `field` and are
/// never deposited as physical heat.
#[derive(Clone, Debug)]
pub struct ConductionResponseThermalStep {
    pub start_s: f64,
    pub end_s: f64,
    pub field: ConductionResponseReport,
    /// Actual native surface losses and their original two-body destinations.
    pub contacts: Vec<contacts::ConductionContactHeat>,
    pub heat_j: f64,
    pub supplied_work_j: f64,
}

struct ResponseCandidate {
    problem: PlanarConductionResponse,
    step: ConductionResponseThermalStep,
}
impl ThermalDriver for ResponseCandidate {
    type Model = PlanarConductionResponse;
    type Step = ConductionResponseThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.problem, self.step)
    }
}

/// Consume the accepted enthalpy history and its source/field state together.
/// Failed stencils, source-domain queries, budgets and cancellation publish none
/// of either candidate; the original histories remain available for retry.
pub struct ConductionResponseEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    problem: PlanarConductionResponse,
    report: CircuitEnthalpyReport<S, ConductionResponseThermalStep>,
}
impl<'t, S> ConductionResponseEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn problem(&self) -> &PlanarConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, ConductionResponseThermalStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarConductionResponse,
        CircuitEnthalpyReport<S, ConductionResponseThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Close source sigma/rho(T, phase, |E|, |J|), spatial current and finite
    /// heat at one common backward-Euler endpoint. Every thermal stencil starts
    /// from the original accepted state and refreshes the SAME complete pinned
    /// conductor bundle; the native response owner closes its E/J coordinates.
    /// Local native E.J powers enter the existing per-cell enthalpy residual.
    /// `maximum_circuit_trials` counts EVERY native field solve, including
    /// constitutive iterations and the final fresh thermal admission.
    ///
    /// Voltage and distributed inward current density retain their native trace
    /// meanings. Total-current electrodes use the original equipotential solve
    /// inside EVERY constitutive iteration. Source work is authored I times
    /// solved electrode voltage; accepted nonlinear fields are never rescaled.
    /// Original geometry, chart identity, source roles/pins and fixed external
    /// conditions remain bound. Source-derived caloric charts retain their
    /// original Cp/density/transition receipts; E/J-dependent response or
    /// applicability refuses. Explicitly authored h(T, phase) charts declare
    /// field independence. This path supplies no reversible field/caloric work.
    /// Thermal conductivity is refreshed at each
    /// proposed T/phase, with field-dependent response/applicability refused.
    /// A phase-aware scalar response retains its declared mass-fraction axis,
    /// supplied by the SAME enthalpy chart at every trial. Its source support
    /// and interpolation remain authoritative. This bounded memoryless scalar
    /// closure is not a convergence, uniqueness or material-validation claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseEnthalpyCandidate<'t>, Error> {
        self.step_response_contact_candidate(
            thermal,
            target_s,
            problem,
            drive,
            &[],
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Close nonlinear bulk E/J response, native surface conduction and both
    /// receiving enthalpies in the same residual. Each thermal trial refreshes
    /// original contact sources and each constitutive trial retains that SAME
    /// surface operator. Explicit interface coordinates and heat fractions keep
    /// their original owners; field-dependent interface laws are not inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_contact_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseEnthalpyCandidate<'t>, Error> {
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
        self.step_response_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
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

    // Both bulk-only and complete material heat assemblies enter the SAME
    // spatial constitutive and enthalpy residual, including the final fresh
    // solve. Retain the caller's source samples rather than discard its surface
    // physics or introduce a second thermal timeline.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_response_resolved_candidate<'t, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, Error>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ConductionResponseEnthalpyCandidate<'t, S>, Error> {
        check(cancel)?;
        // Validate the authored limit BEFORE the outer budget narrows it. The
        // native response owner validates every remaining constitutive control.
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&response.maximum_iterations) {
            return Err(invalid("invalid bounded spatial response iteration limit"));
        }
        let dt = self.admit_mesh_step(
            thermal,
            target_s,
            problem.mesh(),
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let routes = self.admit_mesh_routes(thermal, problem.mesh(), cancel)?;
        self.admit_response_caloric(problem, cancel)?;
        let resolve = |state: &DepositionSimulation<'_>| {
            let coefficients = resolve_sources(state)?;
            self.admit_response_heat(sample(&coefficients.source), cancel)?;
            Ok::<_, Error>(coefficients)
        };
        let accepted = resolve(thermal)?;
        let contacts = contacts::PreparedContacts::new_for_contacts(
            self,
            problem.mesh(),
            problem.contacts(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            cancel,
        )?;
        let initial =
            self.response_materials_at(problem, sample(&accepted.source), &routes, cancel)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted conductor sources must match the actual thermal state",
            ));
        }
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                check(cancel)?;
                let coefficients = resolve(state)?;
                let sources = self.response_materials_at(
                    problem,
                    sample(&coefficients.source),
                    &routes,
                    cancel,
                )?;
                let contact_sources = contacts.refresh_contacts(
                    problem.contacts(),
                    sample(&coefficients.source),
                    cancel,
                )?;
                let bounded = ConductionResponseControls {
                    maximum_iterations: response.maximum_iterations.min(remaining),
                    ..response
                };
                let candidate = match drive {
                    ConductionDrive::Voltage(voltage) => problem
                        .solve_with_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            voltage,
                            &[],
                            assembly,
                            field_controls,
                            bounded,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                    ConductionDrive::CurrentDensity { voltage, currents } => problem
                        .solve_with_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            voltage,
                            currents,
                            assembly,
                            field_controls,
                            bounded,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                    ConductionDrive::Current {
                        positive,
                        negative,
                        current_a,
                    } => problem.solve_current_with_materials_and_contacts_candidate(
                        sources,
                        contact_sources,
                        positive,
                        negative,
                        current_a,
                        assembly,
                        field_controls,
                        bounded,
                        &mut || (*cancel.borrow_mut())(),
                    ),
                }
                .map_err(|error| match error {
                    ConductionResponseError::NotConverged { .. }
                        if bounded.maximum_iterations < response.maximum_iterations =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => response_error(other),
                })?;
                let (next, field) = candidate.into_parts();
                let heat_j = finite(dt * field.power_w)?;
                let supplied_work_j = finite(
                    dt * drive.supplied_power_w(
                        &field.voltage_v,
                        field.power_w,
                        field.power_balance_residual_w,
                    )?,
                )?;
                if (field.power_w > 0.0 && heat_j == 0.0)
                    || (heat_j - supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(invalid(
                        "spatial Joule heat and terminal source work do not close",
                    ));
                }
                let circuit_trials = field.field_solves;
                let contact_heat = contacts.heat(&field.contact_response)?;
                Ok(DrivenCoefficients {
                    circuit: ResponseCandidate {
                        problem: next,
                        step: ConductionResponseThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            contacts: contact_heat,
                            heat_j,
                            supplied_work_j,
                        },
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers =
                    self.deposit_power(&candidate.step.field.cell_power_w, &routes, checkpoint)?;
                contacts::deposit(&candidate.step.contacts, &mut powers, checkpoint)?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel)?;
        Ok(ConductionResponseEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }

    pub(super) fn response_materials_at(
        &self,
        problem: &PlanarConductionResponse,
        sample: &MaterialNetworkSample,
        routes: &[CellId],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ResolvedConductorState>, Error> {
        self.materials_at_with(
            problem.conductors(),
            sample,
            routes,
            cancel,
            |triangle, binding, original, mut point| {
                let source = &problem.sources()[triangle];
                if source.card().content_hash() != binding.point.material.content_hash() {
                    return Err(invalid(
                        "electrical and caloric bindings must use the same card",
                    ));
                }
                let accepted = original.material().query();
                let mut owned = Vec::with_capacity(source.axes().len());
                for (axis, magnitude) in source.axes() {
                    // The thermal reference context may contain its ORIGINAL
                    // zero-field value. Continuation must retain the accepted
                    // spatial operating point before solving the new one.
                    point = point.with(axis, accepted.axes()[axis]).map_err(|_| {
                        invalid("accepted electrical coordinate cannot be retained")
                    })?;
                    owned.push((axis.as_str(), magnitude.quantity()));
                }
                let conductor = source.at_query_point(&point).map_err(|e| {
                    CircuitEnthalpyError::Source(ConductionThermalError::Conductor(e))
                })?;
                for property in conductor.source_properties() {
                    binding
                        .point
                        .admit_static_with_axes(property, &owned)
                        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                }
                Ok(conductor)
            },
        )
    }

    pub(super) fn admit_response_heat(
        &self,
        sample: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        for cell in &sample.cells {
            check(cancel)?;
            let binding = self
                .cells
                .iter()
                .find(|binding| binding.cell == cell.cell)
                .ok_or_else(|| invalid("thermal transport lacks its original binding"))?;
            let owned: Vec<_> = binding
                .point
                .liquid_fraction_axis
                .into_iter()
                .map(|axis| (axis, QuantitySpec::dimensional(Dims([0; 6]))))
                .collect();
            if cell.conductivity.as_ref().is_some_and(|material| {
                material.properties().iter().any(|property| {
                    !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
                })
            }) {
                return Err(invalid(
                    "thermal conductivity response and applicability must be independent of E/J",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn admit_response_caloric(
        &self,
        problem: &PlanarConductionResponse,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        let electrical = |dims| {
            [
                ConductorMagnitude::ElectricField.quantity().dims(),
                ConductorMagnitude::CurrentDensity.quantity().dims(),
            ]
            .contains(&dims)
        };
        for binding in self.cells {
            for receipt in binding.curve.source_receipts() {
                check(cancel)?;
                let claim = binding
                    .point
                    .material
                    .claims()
                    .claim(receipt.selected)
                    .ok_or_else(|| invalid("original caloric source claim is absent"))?;
                let field_axis = |axis: &str| {
                    self.regions
                        .iter()
                        .filter(|region| region.cell == binding.cell)
                        .filter_map(|region| problem.sources().get(region.triangle))
                        .any(|source| source.axes().iter().any(|(name, _)| name == axis))
                        || [
                            claim.validity.axis_quantities().get(axis),
                            receipt.axis_quantities.get(axis),
                            binding.point.context.axis_quantities().get(axis),
                        ]
                        .into_iter()
                        .flatten()
                        .any(|quantity| electrical(quantity.dims()))
                };
                // Inspect the exact original Cp/rho/transition claims. Branch
                // selectors remain legal fixed source-table coordinates; none
                // of these claims is re-selected or frozen at a new field.
                if claim
                    .value
                    .axes()
                    .any(|(axis, dims)| electrical(dims) || field_axis(axis))
                    || claim.validity.bounds().keys().any(|axis| field_axis(axis))
                {
                    return Err(invalid(
                        "caloric source response and applicability must be independent of E/J",
                    ));
                }
            }
        }
        Ok(())
    }
}
