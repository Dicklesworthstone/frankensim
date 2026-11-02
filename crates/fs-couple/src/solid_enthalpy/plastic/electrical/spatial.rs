//! Native tetrahedral current on the SAME proposed thermoplastic configuration.
//!
//! Current is solved on X+u, including shear and nonuniform P1 deformation.
//! The original reference-mass caloric cells, native J2 history and complete heat
//! assembly stay authoritative. This adds no resistor reduction or integrator.

use super::*;
use crate::circuit_enthalpy::conduction::spatial::{
    SpatialConductionDrive, SpatialConductionRegion, SpatialConductionThermalNetwork,
    SpatialConductionThermalStep,
};
use fs_em::{conduction::ConductionControls, spatial::conduction::SpatialConduction};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::ResolvedConductorState;
use std::cell::RefCell;

/// Existing independent mechanical, field, thermal and assembly bounds.
#[derive(Clone, Copy, Debug)]
pub struct SpatialElectricalSolidSolve {
    pub coupled: ElectricalSolidSolve,
    pub field: ConductionControls,
    pub assembly: WeightedAssemblyLimits,
    /// Absolute tolerance for original electrical/solid/caloric reference volume.
    pub reference_volume_tolerance_m3: f64,
}

/// All accepted histories and the actual deformed electrical operator together.
pub struct SpatialElectricalSolidCandidate<'m, 't, S = HeatNetworkSample> {
    thermal: DepositionSimulation<'t>,
    solid: ThermoplasticSolid<'m>,
    problem: SpatialConduction,
    report: ElectricalSolidReport<S, SpatialConductionThermalStep>,
}
impl<'m, 't, S> SpatialElectricalSolidCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &ThermoplasticSolid<'m> {
        &self.solid
    }
    pub const fn problem(&self) -> &SpatialConduction {
        &self.problem
    }
    pub const fn report(&self) -> &ElectricalSolidReport<S, SpatialConductionThermalStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticSolid<'m>,
        SpatialConduction,
        ElectricalSolidReport<S, SpatialConductionThermalStep>,
    ) {
        (self.thermal, self.solid, self.problem, self.report)
    }
}

fn deformed_positions<E>(
    reference: &[[f64; 3]],
    displacement: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, ElectricalSolidSourceError<E>> {
    if reference.len() != displacement.len() {
        return Err(ElectricalSolidSourceError::Invalid(
            "native current must retain every original solid vertex",
        ));
    }
    reference
        .iter()
        .zip(displacement)
        .map(|(x, u)| {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
            let position = [
                number(x[0] + u[0])?,
                number(x[1] + u[1])?,
                number(x[2] + u[2])?,
            ];
            Ok(position)
        })
        .collect()
}

fn source_error<E>(error: ElectricalSolidSourceError<E>) -> ElectricalSolidError<E> {
    match error {
        ElectricalSolidSourceError::Mechanical(ThermoplasticError::Cancelled)
        | ElectricalSolidSourceError::Thermal(DepositionError::Cancelled) => {
            ElectricalSolidError::Cancelled
        }
        other => ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(other)),
    }
}
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> ElectricalSolidError<E> {
    match error {
        HeatNetworkError::Cancelled => ElectricalSolidError::Cancelled,
        other => ElectricalSolidError::Heat(other.into()),
    }
}

impl<'m> ThermoplasticSolid<'m> {
    /// Close actual 3-D current, Joule heat and native loaded-solid response at
    /// one endpoint of the original finite-enthalpy residual.
    ///
    /// The electrical operator must use this solid's exact vertex/tetrahedron
    /// ordering and accepted X+u positions. Its original scalar conductor in
    /// each tet retains the SAME material card, source selection and complete
    /// thermal query. Every trial starts from accepted plastic history, solves
    /// native mechanical equilibrium, refreshes rho/sigma at the same proposed
    /// temperature, then assembles and solves current on that trial's X+u mesh.
    /// No principal-path or affine wire approximation is made.
    ///
    /// Face identities and supplied electrical traces remain material-bound.
    /// Nodal voltage values follow those vertices; prescribed current density
    /// acts on the actual deformed face area and normal. The original native
    /// boundary owner computes source work and charge balance. Per-tet E.J uses
    /// the actual deformed volume and enters its original reference-mass cell
    /// once. Native plastic/reversible heat shares that finite caloric endpoint;
    /// mechanical storage and numerical loss keep their original accounts.
    ///
    /// Accepts the original boundary reference or shared total-current control.
    /// Ideal equipotential terminals follow their original faces on X+u while
    /// holding TOTAL current independently of the changing area. The original
    /// spatial conduction owner solves terminal voltage and the nonuniform
    /// reaction-current distribution at each material/geometry trial. Terminal
    /// work I times the solved voltage drop and local E.J enter the same finite
    /// heat residual without an equivalent resistor or uniform-flux reduction.
    ///
    /// This bounded quasistatic composition retains fully solid isotropic DC
    /// transport and the native small-strain/rotation limits. Thermal geometry,
    /// reference mass, heat-transfer areas and topology remain fixed. Electrical
    /// contacts, source strain-dependence, material-frame tensor evolution, melt,
    /// motional EMF and Lorentz force need their own owners and are not inferred.
    /// In particular, this is not finite-deformation constitutive mechanics.
    ///
    /// Retain all returned owners and the same prepared heat assembly when
    /// continuing. Discarded trials, source exit, cancellation and failed work
    /// gates publish no partial thermal, electrical or mechanical endpoint.
    /// Endpoint geometry/source evaluation gives first-order coupled accuracy;
    /// timestep and spatial refinement remain independent requirements.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial_electrical_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
        &self,
        thermal: &DepositionSimulation<'t>,
        problem: &SpatialConduction<ResolvedConductorState>,
        assembly: &A,
        drive: impl Into<SpatialConductionDrive<'d>>,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: SpatialElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialElectricalSolidCandidate<'m, 't, A::Sample>, ElectricalSolidError<A::Error>>
    {
        self.advance_spatial_electrical_with_motion_candidate(
            thermal,
            problem,
            assembly,
            drive,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            None,
            cx,
        )
        .map(|(candidate, _)| candidate)
    }

    // Quasistatic and inertial consumers use this same field/material/enthalpy
    // residual. Optional velocity selects the original mechanical mass owner;
    // it does not replace the source material or introduce another heat solve.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(in crate::solid_enthalpy::plastic) fn advance_spatial_electrical_with_motion_candidate<
        't,
        'd,
        A: HeatAssembly + ?Sized,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        problem: &SpatialConduction<ResolvedConductorState>,
        assembly: &A,
        drive: impl Into<SpatialConductionDrive<'d>>,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: SpatialElectricalSolidSolve,
        guess_j: &[f64],
        velocity_mps: Option<&[[f64; 3]]>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            SpatialElectricalSolidCandidate<'m, 't, A::Sample>,
            Option<super::super::dynamics::Motion>,
        ),
        ElectricalSolidError<A::Error>,
    > {
        let drive = drive.into();
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        let dt = target_s - self.time_s;
        let tolerance = controls.coupled.coupled.total_energy_tolerance_j;
        if self.time_s != thermal.time_s()
            || !positive(dt)
            || !target_s.is_finite()
            || !positive(tolerance)
            || !positive(controls.coupled.yield_tolerance_pa)
            || controls.coupled.maximum_mechanical_trials == 0
            || self.nodes_m.len() != problem.positions_m().len()
            || self.tetrahedra.len() != problem.mesh().tets.len()
            || self.conductor_geometry.is_some()
            || !problem.contacts().is_empty()
        {
            return Err(ElectricalSolidError::Invalid(
                "matching bounded solid/current histories without lumped paths or moving electrical contacts required",
            ));
        }
        for (solid, electrical) in self.tetrahedra.iter().zip(&problem.mesh().tets) {
            checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
            if *solid != electrical.map(|v| v as usize) {
                return Err(ElectricalSolidError::Invalid(
                    "current and solid must retain the exact original tetrahedron ordering",
                ));
            }
        }
        let accepted_positions =
            deformed_positions::<A::Error>(self.nodes_m, &self.solution.displacement_m, cx)
                .map_err(source_error)?;
        if accepted_positions != problem.positions_m() {
            return Err(ElectricalSolidError::Invalid(
                "accepted current mesh must be the original solid's actual X+u configuration",
            ));
        }
        self.enthalpies::<Infallible>(thermal, true, cx)
            .map_err(|e| match e {
                ThermoplasticError::Cancelled => ElectricalSolidError::Cancelled,
                other => ElectricalSolidError::Mechanical(other),
            })?;
        let prepared = assembly.prepared_heat_network();
        self.admit_heat_network::<Infallible>(thermal, prepared, cx)
            .map_err(|e| match e {
                SolidHeatError::Cancelled => ElectricalSolidError::Cancelled,
                other => ElectricalSolidError::Binding(other),
            })?;
        let bulk = prepared.thermal_material_network().map_err(heat_error)?;
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        let initial = assembly.resolve_heat(cx, thermal).map_err(|e| match e {
            CircuitEnthalpyError::Cancelled => ElectricalSolidError::Cancelled,
            other => {
                ElectricalSolidError::Coupled(other.map_source(ElectricalSolidSourceError::Heat))
            }
        })?;
        if !initial.contacts.is_empty() || !initial.bulk.resistors.is_empty() {
            return Err(ElectricalSolidError::Invalid(
                "spatial solid current cannot discard graph electrical work",
            ));
        }
        self.admit_electrical_transport::<A::Error>(
            thermal,
            initial.bulk.source.as_ref(),
            &[],
            &BTreeMap::new(),
            cx,
        )
        .map_err(source_error)?;
        for (cell, conductor) in self.cells.iter().zip(problem.conductors()) {
            if conductor.liquid_mass_fraction_axis().is_some() {
                return Err(ElectricalSolidError::Invalid(
                    "deforming native solid current requires fully solid transport",
                ));
            }
            self.admit_transport_state::<A::Error>(
                thermal,
                cell.accepted.id,
                conductor.material(),
                cx,
            )
            .map_err(source_error)?;
        }
        let regions: Vec<_> = self
            .cells
            .iter()
            .enumerate()
            .map(|(tetrahedron, c)| SpatialConductionRegion {
                tetrahedron,
                cell: c.accepted.id,
            })
            .collect();
        let volumes: Vec<_> = self.cells.iter().map(|c| c.volume_m3).collect();
        let layout = SpatialConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions: &regions,
            volume_tolerance_m3: controls.reference_volume_tolerance_m3,
        };
        let mechanics = ThermoplasticConfig {
            maximum_iterations: 1,
            temperature_tolerance_k: f64::EPSILON,
            energy_tolerance_j: tolerance,
            yield_tolerance_pa: controls.coupled.yield_tolerance_pa,
            enthalpy_relaxation: 1.0,
        };
        let trials = Cell::new(0usize);
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let result = layout.step_resolved_candidate_with_geometry(
            thermal,
            target_s,
            problem,
            drive,
            &[],
            controls.assembly,
            controls.field,
            controls.coupled.coupled,
            guess_j,
            |endpoint| {
                checkpoint().map_err(CircuitEnthalpyError::from)?;
                if trials.get() >= controls.coupled.maximum_mechanical_trials {
                    return Err(CircuitEnthalpyError::Source(
                        ElectricalSolidSourceError::MechanicalBudget,
                    ));
                }
                trials.set(trials.get() + 1);
                let ContactThermalCoefficients {
                    bulk: mut coefficients,
                    contacts,
                } = assembly
                    .resolve_heat(cx, endpoint)
                    .map_err(|e| e.map_source(ElectricalSolidSourceError::Heat))?;
                if !contacts.is_empty() || !coefficients.resistors.is_empty() {
                    return Err(CircuitEnthalpyError::Invalid(
                        "spatial solid current cannot omit graph electrical work",
                    ));
                }
                (|| {
                    let h = self
                        .enthalpies::<Infallible>(endpoint, false, cx)
                        .map_err(ElectricalSolidSourceError::Mechanical)?;
                    let mechanical = self
                        .mechanics_with_motion::<Infallible>(
                            &h,
                            prescribed_m,
                            nodal_forces_n,
                            settings,
                            mechanics,
                            dt,
                            velocity_mps,
                            cx,
                        )
                        .map_err(ElectricalSolidSourceError::Mechanical)?;
                    let radiation_net_j =
                        number(dt * coefficients.source.as_ref().radiation_net_power_w())?;
                    let external_heat_j = number(
                        number(
                            dt * endpoint.external_step_power_w(
                                &coefficients.sources_w,
                                &coefficients.conductances_w_k,
                                checkpoint,
                            )?,
                        )? - radiation_net_j,
                    )?;
                    let applied_heat_j = self.deposit_mechanical_heat(
                        &mechanical,
                        &mut coefficients.sources_w,
                        dt,
                        tolerance,
                        cx,
                    )?;
                    Ok(CircuitThermalCoefficients {
                        resistors: coefficients.resistors,
                        sources_w: coefficients.sources_w,
                        conductances_w_k: coefficients.conductances_w_k,
                        source: Trial {
                            mechanical,
                            material: coefficients.source,
                            support_contacts: Vec::new(),
                            applied_heat_j,
                            external_heat_j,
                            radiation_net_j,
                        },
                    })
                })()
                .map_err(CircuitEnthalpyError::Source)
            },
            |trial: &Trial<A::Sample>| &trial.material.as_ref().material,
            Some(&volumes),
            |trial: &Trial<A::Sample>| {
                for gradient in &trial.mechanical.displacement_gradients {
                    checkpoint().map_err(CircuitEnthalpyError::from)?;
                    let mut f = *gradient;
                    for (i, row) in f.iter_mut().enumerate() {
                        row[i] += 1.0;
                    }
                    let j = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
                        - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
                        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
                    if !positive(j) {
                        return Err(CircuitEnthalpyError::Invalid(
                            "deforming current cannot invert or collapse a native solid element",
                        ));
                    }
                }
                deformed_positions::<A::Error>(
                    self.nodes_m,
                    &trial.mechanical.solution.displacement_m,
                    cx,
                )
                .map(Some)
                .map_err(CircuitEnthalpyError::Source)
            },
            &cancel,
        );
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        let candidate = result.map_err(|error| match error {
            CircuitEnthalpyError::Cancelled => ElectricalSolidError::Cancelled,
            other => ElectricalSolidError::Coupled(other),
        })?;
        let (thermal, problem, r) = candidate.into_parts();
        let trial = r.source;
        let mechanical = trial.mechanical;
        let thermal_residual_j = r.thermal.stored_energy_change_j
            - r.thermal.incoming_energy_j
            - trial.external_heat_j
            - trial.applied_heat_j
            - r.deposited_joule_j;
        let combined_residual_j = mechanical.stored_change_j
            + mechanical.thermoelastic_storage_change_j
            + r.thermal.stored_energy_change_j
            + mechanical.numerical_loss_j
            + mechanical.motion.as_ref().map_or(0.0, |motion| {
                motion.energy.kinetic_change_j + motion.energy.numerical_loss_j
            })
            - r.electrical.supplied_work_j
            - mechanical.external_work_j
            - trial.external_heat_j
            - r.thermal.incoming_energy_j;
        let source_residual_j = r.thermal.external_heat_j
            - trial.external_heat_j
            - trial.radiation_net_j
            - trial.applied_heat_j
            - r.deposited_joule_j;
        if [thermal_residual_j, combined_residual_j, source_residual_j]
            .iter()
            .any(|r| !r.is_finite() || r.abs() > tolerance)
        {
            return Err(ElectricalSolidError::Invalid(
                "closed spatial electrical, native solid and caloric balances do not agree",
            ));
        }
        let mut solid = self.clone();
        solid.time_s = target_s;
        solid.materials = mechanical.materials;
        solid.solution = mechanical.solution;
        for cell in &mut solid.cells {
            cell.accepted = thermal.cells()[&cell.accepted.id].clone();
        }
        let report = ElectricalSolidReport {
            support_contacts: trial.support_contacts,
            circuit_thermal: CircuitEnthalpyReport {
                electrical: r.electrical,
                thermal: r.thermal,
                source: ContactThermalSample {
                    material: trial.material,
                    destinations: Vec::new(),
                },
                nonlinear: r.nonlinear,
                maximum_cell_residual_j: r.maximum_cell_residual_j,
                circuit_trials: r.circuit_trials,
                deposited_joule_j: r.deposited_joule_j,
                non_electrical_heat_j: r.non_electrical_heat_j,
                total_energy_residual_j: r.total_energy_residual_j,
            },
            element_energy: mechanical.energy,
            mechanical_trials: trials.get(),
            physical_plastic_heat_j: mechanical.physical_heat_j,
            reversible_heat_j: mechanical.reversible_heat_j,
            deposited_mechanical_heat_j: trial.applied_heat_j,
            stored_mechanical_change_j: mechanical.stored_change_j,
            thermoelastic_storage_change_j: mechanical.thermoelastic_storage_change_j,
            numerical_loss_j: mechanical.numerical_loss_j,
            mechanical_work_j: mechanical.external_work_j,
            mechanical_residual_j: mechanical.mechanical_residual_j,
            closed_radiation_net_j: trial.radiation_net_j,
            external_heat_j: trial.external_heat_j,
            thermal_residual_j,
            combined_residual_j,
        };
        checkpoint().map_err(|_| ElectricalSolidError::Cancelled)?;
        Ok((
            SpatialElectricalSolidCandidate {
                thermal,
                solid,
                problem,
                report,
            },
            mechanical.motion,
        ))
    }
}
