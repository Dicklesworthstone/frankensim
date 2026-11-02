//! Tetrahedral Ohmic transport inside the original finite-enthalpy residual.
//!
//! Each proposed endpoint refreshes the complete original source bundles at
//! the actual T/liquid mass fraction, rebuilds only the conductivity Hodge, and
//! deposits native bulk and explicitly partitioned interface Joule work once.
//! Equipotential total-current drive is resolved anew on that actual operator,
//! including the trial bulk and contact resistance, before its heat is deposited.
//! Native thermal storage,
//! source/conductance resolution, Newton work limits, events and acceptance
//! remain with the existing circuit/enthalpy and material-network owners.

use super::contacts::{self, ConductionContactThermal, SolvedContactPressure};
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
        SpatialConductanceControls, SpatialConduction, SpatialConductionBoundary,
        SpatialConductionContact, SpatialConductionContactResponse, SpatialConductionError,
        SpatialConductionField,
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

/// Electrical drive held fixed over one coupled enthalpy interval.
///
/// Existing voltage/current-density traces remain available through `From`.
/// Total current instead fixes equipotential electrodes and lets the native
/// heterogeneous field determine their current-density distribution. Every
/// enthalpy trial rebuilds its terminal response after refreshing the original
/// bulk and contact sources; no accepted-state resistance is frozen in time.
#[derive(Clone, Copy, Debug)]
pub enum SpatialConductionDrive<'a> {
    Boundary(&'a SpatialConductionBoundary),
    /// Signed amperes entering `positive_faces`; `negative_faces` is grounded.
    /// Both sets contain complete exterior triangles and must be disjoint and
    /// electrically connected under the native terminal admission rules.
    /// `reference_voltage_v` is a finite positive numerical probe amplitude,
    /// not the imposed terminal voltage or a second source of Joule work.
    TotalCurrent {
        positive_faces: &'a [usize],
        negative_faces: &'a [usize],
        current_a: f64,
        reference_voltage_v: f64,
    },
    /// Balanced signed amperes INTO every original equipotential electrode.
    /// Terminal zero fixes the voltage gauge; a zero-net-current electrode
    /// floats to its actual potential and may carry current between its faces.
    /// N-1 native basis solves are charged to the shared work budget per trial.
    TerminalCurrents {
        terminals: &'a [&'a [usize]],
        currents_a: &'a [f64],
        response: SpatialConductanceControls,
    },
}

impl<'a> SpatialConductionDrive<'a> {
    /// Use the native one-volt reference probe. Current may be negative or zero;
    /// zero drive still requires an admitted, conducting electrode topology.
    pub fn total_current(
        positive_faces: &'a [usize],
        negative_faces: &'a [usize],
        current_a: f64,
    ) -> Self {
        Self::TotalCurrent {
            positive_faces,
            negative_faces,
            current_a,
            reference_voltage_v: 1.0,
        }
    }

    // The ideal source fixes current, not the recovered reaction. Account its
    // work independently from the Joule integral and native reaction work so
    // the shared energy gate also detects a terminal-current mismatch.
    fn supplied_power_w<C: SpatialConductorMaterial>(
        self,
        problem: &SpatialConduction<C>,
        field: &SpatialConductionField,
    ) -> Result<f64, Error> {
        // Native admission and the completed solve establish nonempty,
        // equipotential exterior face sets for either current-driven mode.
        let potential = |faces: &[usize]| {
            faces
                .first()
                .and_then(|&face| problem.mesh().faces.get(face))
                .and_then(|vertices| field.voltage_v.get(vertices[0] as usize))
                .copied()
                .ok_or_else(|| invalid("current electrode lost its solved potential"))
        };
        let power = match self {
            Self::Boundary(_) => field.supplied_power_w,
            Self::TotalCurrent {
                positive_faces,
                negative_faces,
                current_a,
                ..
            } => current_a * finite(potential(positive_faces)? - potential(negative_faces)?)?,
            Self::TerminalCurrents {
                terminals,
                currents_a,
                ..
            } => {
                if terminals.is_empty() || terminals.len() != currents_a.len() {
                    return Err(invalid(
                        "current electrodes lost their source correspondence",
                    ));
                }
                let reference = potential(terminals[0])?;
                let (mut sum, mut correction) = (0.0, 0.0);
                for (&faces, &current) in terminals.iter().zip(currents_a) {
                    // Independently authored I*V, with the same voltage gauge
                    // as the grounded terminal owner. A floating electrode's
                    // zero net current contributes no source work.
                    let work = finite(current * finite(potential(faces)? - reference)?)?;
                    let adjusted = finite(work - correction)?;
                    let next = finite(sum + adjusted)?;
                    correction = finite((next - sum) - adjusted)?;
                    sum = next;
                }
                sum
            }
        };
        Ok(finite(power)?)
    }
}

impl<'a> From<&'a SpatialConductionBoundary> for SpatialConductionDrive<'a> {
    fn from(boundary: &'a SpatialConductionBoundary) -> Self {
        Self::Boundary(boundary)
    }
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
    /// Pass an existing boundary by reference or [`SpatialConductionDrive`] for
    /// equipotential total-current electrodes. Current is held fixed over the
    /// interval while the actual endpoint source laws determine terminal voltage.
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
    pub fn step_candidate<'t, 'd, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
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
            drive,
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
    pub fn step_contact_candidate<'t, 'd, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
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
            drive,
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
        'd,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
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
            thermal,
            target_s,
            problem,
            drive,
            contact_bindings,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            None,
            |_| Ok(None),
            &[],
            |_| Ok(Vec::new()),
            cancel,
        )
    }

    // A native mechanical owner may supply current positions while caloric
    // inventory remains on its original reference volumes. Matching electrical
    // interfaces retain their original paired vertices and source/heat owners;
    // the native field re-admits their coincidence and actual surface metric at
    // every proposal. An admitted mechanical pressure owner can also refresh
    // original interface coordinates. Accepted pressures come from the stored
    // mechanical state, separately from a target-load solve at old temperatures.
    // Source resolution, field work, heat routes and the one enthalpy residual
    // stay here.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_resolved_candidate_with_geometry<
        't,
        'd,
        'p,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConduction<C>,
        drive: impl Into<SpatialConductionDrive<'d>>,
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
        accepted_pressures: &[SolvedContactPressure<'p>],
        trial_pressures: impl Fn(
            &S,
        ) -> Result<Vec<SolvedContactPressure<'p>>, CircuitEnthalpyError<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionEnthalpyCandidate<'t, C, S>, CircuitEnthalpyError<E>> {
        let drive = drive.into();
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
            accepted_pressures,
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
                let field_solves = match drive {
                    SpatialConductionDrive::TerminalCurrents { terminals, .. } => terminals
                        .len()
                        .checked_sub(1)
                        .filter(|&n| n > 0)
                        .ok_or_else(|| {
                            map_conduction_error(invalid(
                                "multi-terminal current needs at least two electrodes",
                            ))
                        })?,
                    _ => 1,
                };
                if field_solves > remaining {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = resolve(trial)?;
                if !coefficients.resistors.is_empty() {
                    return Err(map_conduction_error(invalid(
                        "3D current cannot silently omit graph electrical work",
                    )));
                }
                let pressures = trial_pressures(&coefficients.source)?;
                let refreshed_contacts = contact_owner
                    .refresh_sources_with_pressures(
                        &contact_sources,
                        sample(&coefficients.source),
                        &pressures,
                        cancel,
                    )
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
                let conductors =
                    refresh(sample(&coefficients.source)).map_err(map_conduction_error)?;
                let positions = trial_positions(&coefficients.source)?;
                if positions.is_some() != reference_volumes_m3.is_some() {
                    return Err(map_conduction_error(invalid(
                        "current geometry and original caloric reference must have one owner",
                    )));
                }
                let next = if let Some(positions) = positions {
                    // Preserve material trace identities, ordered sources and
                    // phase roles. The existing native contact assembly checks
                    // exact coincidence and opposite adjacent volumes on X+u,
                    // then integrates the voltage jump on the CURRENT triangle
                    // area. Opening/sliding therefore refuses atomically; it is
                    // never projected back onto the accepted interface.
                    SpatialConduction::new_with_phase_contacts(
                        positions,
                        problem.mesh().tets.clone(),
                        conductors,
                        refreshed_contacts,
                        problem.contact_liquid_mass_fraction_axes().to_vec(),
                        assembly,
                        &mut || (*cancel.borrow_mut())(),
                    )
                } else {
                    problem.with_conductors_and_contacts(
                        conductors,
                        refreshed_contacts,
                        assembly,
                        &mut || (*cancel.borrow_mut())(),
                    )
                }
                .map_err(field_error)
                .map_err(map_conduction_error)?;
                // A total-current source changes its terminal voltage at every
                // nonlinear trial. Solve on the refreshed native operator so
                // bulk/contact I²R(T,f) feeds back into this same endpoint.
                let field = match drive {
                    SpatialConductionDrive::Boundary(boundary) => {
                        next.solve(boundary, field_controls, &mut || (*cancel.borrow_mut())())
                    }
                    SpatialConductionDrive::TotalCurrent {
                        positive_faces,
                        negative_faces,
                        current_a,
                        reference_voltage_v,
                    } => next.solve_current_with_reference_voltage(
                        positive_faces,
                        negative_faces,
                        current_a,
                        reference_voltage_v,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    ),
                    SpatialConductionDrive::TerminalCurrents {
                        terminals,
                        currents_a,
                        response,
                    } => {
                        let matrix = next
                            .conductance_matrix(terminals, response, field_controls, &mut || {
                                (*cancel.borrow_mut())()
                            })
                            .map_err(field_error)
                            .map_err(map_conduction_error)?;
                        matrix.field_at_currents(
                            currents_a,
                            field_controls.current_tolerance_a,
                            field_controls.power_tolerance_w,
                            &mut || (*cancel.borrow_mut())(),
                        )
                    }
                }
                .map_err(field_error)
                .map_err(map_conduction_error)?;
                let heat_j = finite(dt * field.power_w)?;
                let supplied_work_j = finite(
                    dt * drive
                        .supplied_power_w(&next, &field)
                        .map_err(map_conduction_error)?,
                )?;
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
                    circuit_trials: field_solves,
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
