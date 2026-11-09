//! Actual tetrahedral source response and finite heat at one common endpoint.
//!
//! Every thermal proposal closes the original bulk E/J and explicitly bound
//! contact laws on the native field. Only its admitted local volume/surface
//! Joule work enters the original enthalpy residual; discarded constitutive
//! trials neither advance thermal time nor deposit numerical defects.

mod heat;
/// Nonlinear current, rigid momentum and finite heat at one endpoint.
pub mod mechanical;
/// Nonlinear conductor-frame response under prescribed motion and magnetic field.
pub mod motion;
pub use heat::{
    step_spatial_conduction_response_candidate, step_spatial_moving_conduction_response_candidate,
};

use super::*;
use fs_em::spatial::conduction::response::{
    SpatialConductionResponse, SpatialConductionResponseControls, SpatialConductionResponseDrive,
    SpatialConductionResponseError, SpatialConductionResponseReport,
};

fn response_error(error: SpatialConductionResponseError) -> Error {
    match error {
        SpatialConductionResponseError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::SpatialResponse(other)),
    }
}

fn response_solve_error<E: From<ConductionThermalError>>(
    error: SpatialConductionResponseError,
    bounded: SpatialConductionResponseControls,
    requested: SpatialConductionResponseControls,
) -> CircuitEnthalpyError<E> {
    match error {
        SpatialConductionResponseError::FieldSolveBudget { .. }
            if bounded.maximum_field_solves < requested.maximum_field_solves =>
        {
            CircuitEnthalpyError::CircuitBudget
        }
        other => map_conduction_error(response_error(other)),
    }
}

/// Actual endpoint field/source evaluations and independently supplied work.
/// Constitutive defects remain diagnostics, never additional thermal sources.
#[derive(Clone, Debug)]
pub struct SpatialConductionResponseThermalStep<F = SpatialConductionField> {
    pub start_s: f64,
    pub end_s: f64,
    pub field: SpatialConductionResponseReport<F>,
    pub contacts: Vec<SpatialConductionContactHeat>,
    pub heat_j: f64,
    /// Signed electrical work supplied through the native boundary/electrodes.
    pub supplied_work_j: f64,
    /// Signed Lorentz work delivered to prescribed motion; zero at rest.
    pub mechanical_work_j: f64,
    /// Electrical input minus mechanical output, available to finite heat.
    pub net_supplied_work_j: f64,
}

// Both frames use one thermal residual and one original-source refresh path.
// The native field owner supplies local Joule heat and signed mechanical work.
trait ResponseHeatField {
    fn joule_power_w(&self) -> f64;
    fn cell_joule_power_w(&self) -> &[f64];
    fn contact_response(&self) -> &[SpatialConductionContactResponse];
    fn mechanical_power_w(&self) -> f64;
}

impl ResponseHeatField for SpatialConductionField {
    fn joule_power_w(&self) -> f64 {
        self.power_w
    }
    fn cell_joule_power_w(&self) -> &[f64] {
        &self.cell_power_w
    }
    fn contact_response(&self) -> &[SpatialConductionContactResponse] {
        &self.contacts
    }
    fn mechanical_power_w(&self) -> f64 {
        0.0
    }
}

struct ResponseCandidate<F> {
    problem: SpatialConductionResponse,
    step: SpatialConductionResponseThermalStep<F>,
}
impl<F> ThermalDriver for ResponseCandidate<F> {
    type Model = SpatialConductionResponse;
    type Step = SpatialConductionResponseThermalStep<F>;
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

/// Publish the source/field state and its original finite heat history together.
/// A failed source query, trial, work limit or cancellation publishes neither.
pub struct SpatialConductionResponseEnthalpyCandidate<
    't,
    S = MaterialNetworkSample,
    F = SpatialConductionField,
> {
    thermal: DepositionSimulation<'t>,
    problem: SpatialConductionResponse,
    report: CircuitEnthalpyReport<S, SpatialConductionResponseThermalStep<F>>,
}
impl<'t, S, F> SpatialConductionResponseEnthalpyCandidate<'t, S, F> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn problem(&self) -> &SpatialConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, SpatialConductionResponseThermalStep<F>> {
        &self.report
    }
    /// Local conductor-frame Joule and contact heat entered the returned state.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialConductionResponse,
        CircuitEnthalpyReport<S, SpatialConductionResponseThermalStep<F>>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Close nonlinear scalar sigma/rho(T, phase, |E|, |J|), native tetrahedral
    /// current and finite enthalpy at the SAME backward-Euler endpoint. The
    /// accepted source laws, geometry, caloric charts, pins and non-electrical
    /// context stay bound. Every thermal stencil starts from the original
    /// history; native source-coordinate relaxation is retained during feedback
    /// and continuation without weakening actual-field/current/power gates.
    ///
    /// Boundary voltage/current density and mixed equipotential electrodes keep
    /// the native meanings. Direct `Terminals` drives charge one actual-state
    /// field solve per constitutive iteration; matrix-based mixed drives charge
    /// their N-1 basis solves. Every solve consumes the original
    /// `maximum_circuit_trials` allowance.
    /// Caloric and thermal transport sources must be independent of electrical
    /// response/applicability; this path supplies no reversible field storage.
    /// Fixed geometry and confined phase-chart limits still apply. This bounded
    /// first-order closure makes no global convergence or accuracy guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_candidate<'t, 'd>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseEnthalpyCandidate<'t>, Error> {
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

    /// Include source-bound matching-face resistance and its explicit two-body
    /// heat partition in this same endpoint. Both contact temperatures/phases
    /// refresh at each thermal proposal. Explicitly bound nonlinear trace axes
    /// are closed by the native response owner, with its original uniform-trace
    /// restrictions; an unowned electric interface law still refuses.
    /// Bulk and surface heat enter their original receivers exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_contact_candidate<'t, 'd>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit_mesh_step(
            thermal,
            target_s,
            problem.mesh(),
            problem.contacts().len(),
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
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
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    // Complete heat assemblies supply their original source sample to the same
    // residual. Geometry-independent material/volume/heat admission is shared
    // with the existing planar owner; no surrogate planar field is constructed.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_response_resolved_candidate<
        't,
        'd,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'t>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionResponseEnthalpyCandidate<'t, S>, CircuitEnthalpyError<E>> {
        let drive = drive.into();
        let work_drive = match drive {
            SpatialConductionResponseDrive::Terminals(drive) => {
                SpatialConductionDrive::Terminals(drive)
            }
            SpatialConductionResponseDrive::Power {
                positive_faces,
                negative_faces,
                source,
                reference_voltage_v,
            } => SpatialConductionDrive::Power {
                positive_faces,
                negative_faces,
                source,
                reference_voltage_v,
            },
            SpatialConductionResponseDrive::Boundary(boundary) => {
                SpatialConductionDrive::Boundary(boundary)
            }
            SpatialConductionResponseDrive::MixedTerminals {
                terminals,
                current_terminals,
                voltages_v,
                currents_a,
                response,
                mixed,
            } => SpatialConductionDrive::MixedTerminals {
                terminals,
                current_terminals,
                voltages_v,
                currents_a,
                response,
                mixed,
            },
        };
        self.step_response_driver_candidate(
            thermal,
            target_s,
            problem,
            contact_bindings,
            assembly,
            response,
            controls,
            initial_energy_increments_j,
            resolve_sources,
            sample,
            |sources, contacts, bounded| {
                let candidate = problem
                    .solve_with_materials_and_contacts_candidate(
                        sources,
                        contacts,
                        drive,
                        assembly,
                        field_controls,
                        bounded,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| response_solve_error(error, bounded, response))?;
                let (next, field) = candidate.into_parts();
                let supplied_power_w = match &field.terminals {
                    Some(terminals) => terminals.supplied_power_w,
                    None => work_drive
                        .supplied_power_w(next.mesh(), &field.field, cancel)
                        .map_err(map_conduction_error)?,
                };
                Ok((next, field, supplied_power_w))
            },
            cancel,
        )
    }

    // The shared endpoint owns every source/phase refresh, thermal trial,
    // physical work gate and heat partition. Frame adapters only solve the
    // native field and expose its independently recovered terminal work.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_response_driver_candidate<
        't,
        S,
        E: From<ConductionThermalError>,
        F: ResponseHeatField,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'t>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        solve_field: impl Fn(
            Vec<ResolvedConductorState>,
            Vec<SpatialConductionContact>,
            SpatialConductionResponseControls,
        ) -> Result<
            (
                SpatialConductionResponse,
                SpatialConductionResponseReport<F>,
                f64,
            ),
            CircuitEnthalpyError<E>,
        >,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionResponseEnthalpyCandidate<'t, S, F>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(map_conduction_error)?;
        if !(1..=fs_em::conduction::response::MAX_RESPONSE_ITERATIONS)
            .contains(&response.maximum_iterations)
            || response.maximum_field_solves == 0
        {
            return Err(map_conduction_error(invalid(
                "invalid bounded 3D source response limits",
            )));
        }
        let dt = self
            .admit_mesh_step(
                thermal,
                target_s,
                problem.mesh(),
                problem.contacts().len(),
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        let regions: Vec<_> = self
            .regions
            .iter()
            .map(|r| ConductionRegion {
                triangle: r.tetrahedron,
                cell: r.cell,
            })
            .collect();
        let route_owner = ConductionThermalNetwork {
            cells: self.cells,
            links: self.links,
            regions: &regions,
            volume_tolerance_m3: self.volume_tolerance_m3,
        };
        let volumes = (0..problem.mesh().tets.len())
            .map(|tet| {
                check(cancel)?;
                problem
                    .cell_volume_m3(tet)
                    .ok_or_else(|| invalid("original tetrahedral volume is absent"))
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(map_conduction_error)?;
        let routes = route_owner
            .admit_volume_routes(thermal, &volumes, cancel)
            .map_err(map_conduction_error)?;
        route_owner
            .admit_response_caloric_sources(problem.sources(), cancel)
            .map_err(map_conduction_error)?;
        let resolve = |state: &DepositionSimulation<'t>| {
            let coefficients = resolve_sources(state)?;
            if !coefficients.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "3D response cannot omit graph electrical work",
                )));
            }
            route_owner
                .admit_response_heat(sample(&coefficients.source), cancel)
                .map_err(map_conduction_error)?;
            Ok::<_, CircuitEnthalpyError<E>>(coefficients)
        };
        let refresh = |state: &MaterialNetworkSample| {
            route_owner
                .response_materials_from(
                    problem.conductors(),
                    problem.sources(),
                    state,
                    &routes,
                    cancel,
                )
                .map_err(map_conduction_error)
        };
        let accepted = resolve(thermal)?;
        let initial = refresh(sample(&accepted.source))?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(map_conduction_error(invalid(
                "accepted tetrahedral response sources must match the actual thermal state",
            )));
        }
        let contact_sources = problem
            .contacts()
            .iter()
            .map(|c| &c.source)
            .collect::<Vec<_>>();
        let contact_axes = contacts::response_axes(problem.contact_response_sources());
        let contact_owner = contacts::PreparedContacts::new_for_response_sources(
            &route_owner,
            &contact_sources,
            problem.contact_cells(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            &[],
            &contact_axes,
            cancel,
        )
        .map_err(map_conduction_error)?;
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(cancel).map_err(map_conduction_error)?;
                if remaining == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = resolve(trial)?;
                let sources = refresh(sample(&coefficients.source))?;
                let contacts = contact_owner
                    .refresh_sources(&contact_sources, sample(&coefficients.source), cancel)
                    .map_err(map_conduction_error)?;
                let contacts = problem
                    .contacts()
                    .iter()
                    .zip(contacts)
                    .map(|(original, source)| SpatialConductionContact {
                        first_vertices: original.first_vertices,
                        second_vertices: original.second_vertices,
                        source,
                    })
                    .collect();
                let bounded = SpatialConductionResponseControls {
                    maximum_field_solves: response.maximum_field_solves.min(remaining),
                    ..response
                };
                let (next, field, supplied_power_w) = solve_field(sources, contacts, bounded)?;
                let joule_power_w = field.field.joule_power_w();
                let heat_j = finite(dt * joule_power_w)?;
                let supplied_work_j = finite(dt * supplied_power_w)?;
                let mechanical_work_j = finite(dt * field.field.mechanical_power_w())?;
                let net_supplied_work_j = finite(supplied_work_j - mechanical_work_j)?;
                if (joule_power_w > 0.0 && heat_j == 0.0)
                    || (heat_j - net_supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(map_conduction_error(invalid(
                        "nonlinear 3D Joule heat and terminal work minus Lorentz work do not close",
                    )));
                }
                let contact_heat = contact_owner
                    .heat(field.field.contact_response())
                    .map_err(map_conduction_error)?;
                let circuit_trials = field.field_solves;
                Ok(DrivenCoefficients {
                    circuit: ResponseCandidate {
                        problem: next,
                        step: SpatialConductionResponseThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            contacts: contact_heat,
                            heat_j,
                            supplied_work_j,
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
                let mut powers = route_owner
                    .deposit_power(
                        candidate.step.field.field.cell_joule_power_w(),
                        &routes,
                        checkpoint,
                    )
                    .map_err(map_conduction_error)?;
                for contact in &candidate.step.contacts {
                    checkpoint()?;
                    contacts::deposit_partition(
                        contact.cells,
                        contact.heat_fraction_to_first,
                        contact.response.power_w,
                        &mut powers,
                    )
                    .map_err(map_conduction_error)?;
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(map_conduction_error)?;
        Ok(SpatialConductionResponseEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}
