//! Tetrahedral Ohmic transport inside the original finite-enthalpy residual.
//!
//! Each proposed endpoint refreshes the complete original source bundles at
//! the actual T/liquid mass fraction, rebuilds only the conductivity Hodge, and
//! deposits native bulk and explicitly partitioned interface Joule work once.
//! Native thermal storage,
//! source/conductance resolution, Newton work limits, events and acceptance
//! remain with the existing circuit/enthalpy and material-network owners.

use super::contacts::{self, ConductionContactThermal};
use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, ConductionLink,
    ConductionRegion, ConductionThermalError, ConductionThermalNetwork, DrivenCoefficients,
    DriverAccount, Error, MaterialCell, MaterialNetwork, MaterialNetworkError,
    MaterialNetworkSample, ThermalDriver, cell_driver_candidate, check, finite, invalid,
    map_conduction_error, material, same_material,
};
use crate::circuit_enthalpy::CircuitThermalCoefficients;
use fs_em::{
    conduction::ConductionControls,
    spatial::conduction::{
        SpatialConduction, SpatialConductionBoundary, SpatialConductionContact,
        SpatialConductionContactResponse, SpatialConductionError, SpatialConductionField,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial},
    thermal::source::ThermalMaterialError,
};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialConductionRegion {
    pub tetrahedron: usize,
    pub cell: CellId,
}

/// Every original tetrahedron has exactly one finite thermal receiver. Several
/// tets may cover one receiver; their volumes must cover that cell completely.
pub struct SpatialConductionThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub regions: &'a [SpatialConductionRegion],
    pub volume_tolerance_m3: f64,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionThermalStep {
    pub start_s: f64,
    pub end_s: f64,
    pub field: SpatialConductionField,
    /// Exact triangular contact response and the explicit A/B heat partition.
    /// Its heat is already included in the finite cells and `heat_j` once.
    pub contacts: Vec<SpatialConductionContactHeat>,
    pub heat_j: f64,
    pub supplied_work_j: f64,
}

pub type SpatialConductionContactHeat =
    contacts::ConductionContactHeat<SpatialConductionContactResponse>;

struct Candidate<C> {
    problem: SpatialConduction<C>,
    step: SpatialConductionThermalStep,
}
impl<C: SpatialConductorMaterial> ThermalDriver for Candidate<C> {
    type Model = SpatialConduction<C>;
    type Step = SpatialConductionThermalStep;
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

pub struct SpatialConductionEnthalpyCandidate<
    't,
    C = ResolvedConductorState,
    S = MaterialNetworkSample,
> {
    thermal: DepositionSimulation<'t>,
    problem: SpatialConduction<C>,
    report: CircuitEnthalpyReport<S, SpatialConductionThermalStep>,
}
impl<'t, C: SpatialConductorMaterial, S> SpatialConductionEnthalpyCandidate<'t, C, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn problem(&self) -> &SpatialConduction<C> {
        &self.problem
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, SpatialConductionThermalStep> {
        &self.report
    }
    /// Publish the admitted source operator and finite enthalpy together. Joule
    /// heat already entered this thermal state and must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialConduction<C>,
        CircuitEnthalpyReport<S, SpatialConductionThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

fn field_error(error: SpatialConductionError) -> Error {
    match error {
        SpatialConductionError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::SpatialField(other)),
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Close the actual 3D current field, source sigma/rho(T,f), finite-cell
    /// enthalpy and all original material heat links at one common endpoint.
    /// Scalar and fixed-orientation principal transport share the same solve.
    /// A confined phase chart needs equal, constant density and an explicit
    /// liquid-MASS-fraction coordinate; its original source owns sigma(f), not
    /// a manufactured mixing rule. Geometry/mass and the source roles stay fixed.
    ///
    /// Every residual/finite-difference/final evaluation uses the same accepted
    /// thermal history. The original maximum_circuit_trials budget counts each
    /// native spatial solve. Failures, cancellation and unsupported event windows
    /// leave both inputs untouched. Backward Euler and the declared tolerances
    /// do not establish unconditional convergence or temporal/spatial accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionEnthalpyCandidate<'t, C>, Error> {
        self.step_contact_candidate(
            thermal,
            target_s,
            problem,
            boundary,
            &[],
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Advance the original tetrahedral field and every finite caloric cell
    /// with source-bound resistance on matching triangular material contacts.
    /// Every interface is refreshed at BOTH actual endpoint temperatures and
    /// owned liquid mass fractions while its pressure/history and selected
    /// source claims stay pinned. Its exact surface work is divided between the
    /// original adjacent volume receivers using an explicit physical fraction;
    /// no equal split, effective contact temperature or surrogate resistor is
    /// inferred. Bulk E.J and surface heat enter the same native residual once.
    ///
    /// All contacts require a binding. `step_candidate` is the empty-binding
    /// case and refuses contacts without their thermal source/heat ownership.
    /// Fixed matching traces, zero interface storage and the confined-phase
    /// density rules still apply. Rejected trials publish neither state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_contact_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
        contacts: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit_step(
            thermal,
            target_s,
            problem,
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
        self.step_resolved_candidate(
            thermal,
            target_s,
            problem,
            boundary,
            contacts,
            assembly,
            field_controls,
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

    // Admit before material preparation or any numerical trial. The complete
    // heat assembly and the bare Fourier entry point share these limits.
    fn admit_step<C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        assembly: WeightedAssemblyLimits,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<f64, Error> {
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if dt <= 0.0
            || self.cells.len() != thermal.cells().len()
            || self.regions.len() != problem.mesh().tets.len()
            || problem.mesh().tets.len() > assembly.max_cells
            || problem.mesh().vertex_count > assembly.max_dofs
            || problem.mesh().edges.len() > assembly.max_dofs
            || problem
                .mesh()
                .tets
                .len()
                .checked_mul(36)
                .and_then(|bulk| {
                    problem
                        .contacts()
                        .len()
                        .checked_mul(36)
                        .and_then(|contact| bulk.checked_add(contact))
                })
                .is_none_or(|n| n > assembly.max_triplets)
        {
            return Err(invalid(
                "positive bounded 3D interval and complete caloric/volume correspondence required",
            ));
        }
        Ok(dt)
    }

    /// Compose a complete native heat assembly with the original tetrahedral
    /// field/enthalpy residual. Every source evaluation observes the SAME trial
    /// endpoint; no additional field or thermal time integrator is introduced.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_resolved_candidate<
        't,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionEnthalpyCandidate<'t, C, S>, CircuitEnthalpyError<E>> {
        self.step_resolved_candidate_with_geometry(
            thermal, target_s, problem, boundary, contact_bindings, assembly,
            field_controls, controls, initial_energy_increments_j, resolve, sample,
            None, |_| Ok(None), cancel,
        )
    }

    // A native mechanical owner may supply current positions while caloric
    // inventory remains on its original reference volumes. Source resolution,
    // actual field work, heat routes and the one enthalpy residual stay here.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_resolved_candidate_with_geometry<
        't,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        boundary: &SpatialConductionBoundary,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        reference_volumes_m3: Option<&[f64]>,
        trial_positions: impl Fn(&S) -> Result<Option<Vec<[f64; 3]>>, CircuitEnthalpyError<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionEnthalpyCandidate<'t, C, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(map_conduction_error)?;
        let dt = self
            .admit_step(
                thermal,
                target_s,
                problem,
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        if reference_volumes_m3.is_some()
            && (!problem.contacts().is_empty() || !contact_bindings.is_empty())
        {
            return Err(map_conduction_error(invalid(
                "deforming current geometry requires its own moving electrical-interface owner",
            )));
        }
        if contact_bindings.len() != problem.contacts().len() {
            return Err(map_conduction_error(invalid(
                "every tetrahedral electrical contact needs explicit source and heat ownership",
            )));
        }
        // The existing geometry-independent inventory and heat routing owner
        // consumes the actual tetrahedral volumes. No surrogate planar mesh is
        // constructed and no second mass/density admission rule is introduced.
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
        let volumes = if let Some(volumes) = reference_volumes_m3 {
            if volumes.len() != problem.mesh().tets.len()
                || volumes.iter().any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(map_conduction_error(invalid(
                    "deforming current requires every original positive reference volume",
                )));
            }
            volumes.to_vec()
        } else {
            (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("original tetrahedral volume is missing"))
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(map_conduction_error)?
        };
        let routes = route_owner
            .admit_volume_routes(thermal, &volumes, cancel)
            .map_err(map_conduction_error)?;
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|c| (c.cell, c)).collect();
        let refresh = |sample: &MaterialNetworkSample| -> Result<Vec<C>, Error> {
            let phases: BTreeMap<_, _> = sample.cells.iter().map(|c| (c.cell, c.phase)).collect();
            routes.iter().zip(problem.conductors()).map(|(&cell, source)| {
                check(cancel)?;
                let binding = bindings.get(&cell)
                    .ok_or_else(|| invalid("3D conductor lacks its original thermal material"))?;
                if source.liquid_mass_fraction_axis().is_some()
                    && source.liquid_mass_fraction_axis() != binding.point.liquid_fraction_axis
                {
                    return Err(invalid("3D electrical phase coordinate must be driven by the original caloric fraction"));
                }
                let phase = *phases.get(&cell)
                    .ok_or_else(|| invalid("3D conductor lacks the actual endpoint phase"))?;
                let point = binding.point.at(phase)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                let snapshot = source.material().at_query_point(binding.point.material, &point)
                    .map_err(|e| material(MaterialNetworkError::Material(ThermalMaterialError::Material(e))))?;
                let refreshed = source.with_material(snapshot)
                    .map_err(|e| CircuitEnthalpyError::Source(ConductionThermalError::Conductor(e)))?;
                for key in refreshed.source_keys() {
                    let property = refreshed.material().property_by_key(key)
                        .ok_or_else(|| invalid("3D source lost an original principal electrical role"))?;
                    binding.point.admit_static(property)
                        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                }
                Ok(refreshed)
            }).collect()
        };
        let accepted = resolve(thermal)?;
        if !accepted.resistors.is_empty() {
            return Err(map_conduction_error(invalid(
                "3D current cannot silently omit graph electrical work",
            )));
        }
        let initial = refresh(sample(&accepted.source)).map_err(map_conduction_error)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(map_conduction_error(invalid(
                "accepted tetrahedral conductor differs from its actual thermal state",
            )));
        }
        let contact_sources = problem
            .contacts()
            .iter()
            .map(|contact| &contact.source)
            .collect::<Vec<_>>();
        let contact_owner = contacts::PreparedContacts::new_for_sources(
            &route_owner,
            &contact_sources,
            problem.contact_cells(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            cancel,
        )
        .map_err(map_conduction_error)?;
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _| {
                check(cancel).map_err(map_conduction_error)?;
                let coefficients = resolve(trial)?;
                if !coefficients.resistors.is_empty() {
                    return Err(map_conduction_error(invalid(
                        "3D current cannot silently omit graph electrical work",
                    )));
                }
                let refreshed_contacts = contact_owner
                    .refresh_sources(&contact_sources, sample(&coefficients.source), cancel)
                    .map_err(map_conduction_error)?;
                let refreshed_contacts = problem
                    .contacts()
                    .iter()
                    .zip(refreshed_contacts)
                    .map(|(original, source)| SpatialConductionContact {
                        first_vertices: original.first_vertices,
                        second_vertices: original.second_vertices,
                        source,
                    })
                    .collect();
                let conductors = refresh(sample(&coefficients.source))
                    .map_err(map_conduction_error)?;
                let positions = trial_positions(&coefficients.source)?;
                if positions.is_some() != reference_volumes_m3.is_some() {
                    return Err(map_conduction_error(invalid(
                        "current geometry and original caloric reference must have one owner",
                    )));
                }
                let next = if let Some(positions) = positions {
                    SpatialConduction::new(
                        positions, problem.mesh().tets.clone(), conductors, assembly,
                        &mut || (*cancel.borrow_mut())(),
                    )
                } else {
                    problem.with_conductors_and_contacts(
                        conductors, refreshed_contacts, assembly,
                        &mut || (*cancel.borrow_mut())(),
                    )
                }
                .map_err(field_error)
                .map_err(map_conduction_error)?;
                let field = next
                    .solve(boundary, field_controls, &mut || (*cancel.borrow_mut())())
                    .map_err(field_error)
                    .map_err(map_conduction_error)?;
                let heat_j = finite(dt * field.power_w)?;
                let supplied_work_j = finite(dt * field.supplied_power_w)?;
                if (field.power_w > 0.0 && heat_j == 0.0)
                    || (heat_j - supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(map_conduction_error(invalid(
                        "3D Joule heat and independently supplied terminal work do not close",
                    )));
                }
                let contacts = contact_owner
                    .heat(&field.contacts)
                    .map_err(map_conduction_error)?;
                Ok(DrivenCoefficients {
                    circuit: Candidate {
                        problem: next,
                        step: SpatialConductionThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            contacts,
                            heat_j,
                            supplied_work_j,
                        },
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers = route_owner
                    .deposit_power(&candidate.step.field.cell_power_w, &routes, checkpoint)
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
        Ok(SpatialConductionEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}
