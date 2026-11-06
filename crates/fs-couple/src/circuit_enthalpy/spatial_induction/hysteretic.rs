//! Retained three-dimensional magnetic memory and finite caloric feedback.
//!
//! The original spatial A–phi owner advances every vector play state from the
//! same accepted history at each thermal trial. Its exact local Joule heat and
//! irreversible pinning heat enter their original finite receivers once;
//! numerical endpoint dissipation remains a separate energy export. Conductivity
//! follows the actual endpoint temperature, while each magnetic fit stays fixed
//! over its complete admitted caloric interval. Geometry and material occupancy
//! are stationary; this adapter supplies no transport of magnetic memory.

/// Native winding circuits on the same retained vector-memory thermal solve.
pub mod winding;

use super::{
    CircuitEnthalpyError, CircuitEnthalpySolve, CircuitThermalCoefficients, DrivenCoefficients,
    DriverAccount, Error, MagneticThermalError, MagneticThermalNetwork, MaterialNetwork,
    MaterialNetworkSample, RefCell, SpatialInductionThermalReport, SpatialInductionThermalSource,
    SpatialMagneticMaterialSample, SpatialMagneticThermalNetwork, ThermalDriver,
    cell_driver_candidate, check, deposit_cell_heat, finite, invalid, material, native_error,
    same_material,
};
use crate::circuit_enthalpy::magnetic::MagneticRoutes;
use fs_em::spatial::{
    SpatialMagneticError,
    vector_potential::induction::hysteretic::{
        SpatialHystereticInduction, SpatialHystereticInductionCandidate,
        SpatialHystereticInductionControls, SpatialHystereticInductionError,
        SpatialHystereticInductionStep,
    },
};
use fs_material::{
    conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial},
    field_storage::hysteresis::vector::VectorHysteresisError,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// One accepted A–phi/vector-memory field and its original finite thermal state.
pub struct SpatialHystereticInductionThermalCandidate<
    't,
    'm,
    S = MaterialNetworkSample,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    thermal: DepositionSimulation<'t>,
    field: SpatialHystereticInduction<'m, C>,
    pub(super) report: SpatialInductionThermalReport<SpatialHystereticInductionStep, S>,
}

impl<'t, 'm, S, C: SpatialConductorMaterial>
    SpatialHystereticInductionThermalCandidate<'t, 'm, S, C>
{
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &SpatialHystereticInduction<'m, C> {
        &self.field
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<SpatialHystereticInductionStep, S> {
        &self.report
    }
    /// Physical Joule plus pinning heat already deposited into the finite cells.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.coupled.deposited_joule_j
    }
    /// Publish both histories together. Reported physical heat is already paid;
    /// the numerical loss is not heat and must not be deposited by the caller.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialHystereticInduction<'m, C>,
        SpatialInductionThermalReport<SpatialHystereticInductionStep, S>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver for SpatialHystereticInductionCandidate<'m, C> {
    type Model = SpatialHystereticInduction<'m, C>;
    type Step = SpatialHystereticInductionStep;
    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialHystereticInductionCandidate::into_parts(self)
    }
}

// Current excitation and native circuit terminals share the original magnetic
// source admission, finite-enthalpy residual and per-tetrahedron heat accounts.
trait HystereticInductionDriver: ThermalDriver {
    fn induction_step(&self) -> &SpatialHystereticInductionStep;
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_>;
    fn field_trials(&self) -> usize;
}

impl<C: SpatialConductorMaterial> HystereticInductionDriver
    for SpatialHystereticInductionCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialHystereticInductionStep {
        self.step()
    }
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_change_j,
            supplied_j: step.impressed_work_j,
        }
    }
    fn field_trials(&self) -> usize {
        1
    }
}

fn induction(error: SpatialHystereticInductionError) -> Error {
    match error {
        SpatialHystereticInductionError::Spatial(SpatialMagneticError::Cancelled)
        | SpatialHystereticInductionError::Material(VectorHysteresisError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => {
            CircuitEnthalpyError::Source(MagneticThermalError::SpatialHystereticInduction(other))
        }
    }
}

fn deposit_heat(
    step: &SpatialHystereticInductionStep,
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    let mut powers = deposit_cell_heat(&step.cells, routes, dt, checkpoint)?;
    if step.magnetic_cells.len() != routes.magnetic.len() {
        return Err(invalid(
            "pinning heat changed original tetrahedron correspondence",
        ));
    }
    // Electrically insulating magnetic tetrahedra can dissipate pinning heat.
    // The conductor list alone therefore cannot own this deposition.
    for loss in &step.magnetic_cells {
        checkpoint()?;
        let receiver = routes
            .magnetic
            .get(&loss.cell)
            .ok_or_else(|| invalid("pinning heat lacks its original spatial magnetic region"))?;
        let heat = finite(loss.irreversible_heat_j)?;
        if heat < 0.0 {
            return Err(invalid("negative native spatial pinning heat"));
        }
        if heat == 0.0 {
            continue;
        }
        let receiver = receiver.ok_or_else(|| {
            invalid(
                "every dissipating magnetic tetrahedron needs an actual thermal heat destination",
            )
        })?;
        let power = finite(heat / dt)?;
        let previous = *powers.get(&receiver).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if power == 0.0 || combined == previous {
            return Err(invalid(
                "spatial pinning heat disappeared in receiving-cell sum",
            ));
        }
        powers.insert(receiver, combined);
    }
    Ok(powers)
}

impl SpatialMagneticThermalNetwork<'_> {
    /// Continue spatial vector hysteresis, charge-conserving eddy currents and
    /// finite enthalpy at one jointly accepted endpoint.
    ///
    /// Original scalar or oriented tensor conductivity responds to each actual
    /// caloric temperature. Magnetic source roles, pinning thresholds, branch
    /// slopes and saturation bounds remain constant over the complete source
    /// interval. A changing hysteresis fit needs an additional thermodynamic
    /// law and refuses instead of retaining incompatible memory.
    ///
    /// Conducting tetrahedra and all dissipating magnetic tetrahedra require
    /// real whole-volume, fully solid caloric receivers. Other magnetic regions
    /// may be fixed only when they dissipate no pinning heat. Local Joule and
    /// pinning heat is deposited once; numerical dissipation is never heated.
    /// Every trial freezes the accepted A, vector polarizations and thermal
    /// history. Cancellation or any failed source/balance check publishes none.
    /// Geometry, material occupancy and the original first-order scheme remain
    /// fixed. Source switches and physical reversals require resolved endpoints.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialHystereticInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialHystereticInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialHystereticInductionThermalCandidate<'t, 'm, MaterialNetworkSample, C>, Error>
    {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let material_network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = material_network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_hysteretic_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            controls,
            field_controls,
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

    // Complete contact/radiation/coolant assemblies enter this same residual;
    // there is no second temperature history or separately stepped heat path.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_hysteretic_resolved_candidate<
        't,
        'm,
        S,
        E: From<MagneticThermalError>,
        C: SpatialConductorMaterial,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialHystereticInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialHystereticInductionControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialHystereticInductionThermalCandidate<'t, 'm, S, C>, CircuitEnthalpyError<E>>
    {
        let (thermal, field, report) = self.step_hysteretic_driver_candidate(
            thermal,
            target_s,
            field,
            controls,
            initial_energy_increments_j,
            resolve,
            snapshot,
            |_, _, conductors, _| {
                field
                    .step_with_conductors_candidate(
                        endpoint_current_a_m2,
                        target_s - thermal.time_s(),
                        conductors,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)
                    .map_err(native_error)
            },
            |_, _, _, _| Ok(()),
            cancel,
        )?;
        Ok(SpatialHystereticInductionThermalCandidate {
            thermal,
            field,
            report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_hysteretic_driver_candidate<
        't,
        'm,
        S,
        E: From<MagneticThermalError>,
        C: SpatialConductorMaterial,
        D: HystereticInductionDriver,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialHystereticInduction<'m, C>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        drive: impl Fn(
            &DepositionSimulation<'_>,
            &CircuitThermalCoefficients<S>,
            Vec<Option<C>>,
            usize,
        ) -> Result<D, CircuitEnthalpyError<E>>,
        additional_heat: impl Fn(
            &D,
            &mut BTreeMap<CellId, f64>,
            f64,
            &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<(), CircuitEnthalpyError<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            DepositionSimulation<'t>,
            D::Model,
            SpatialInductionThermalReport<D::Step, S>,
        ),
        CircuitEnthalpyError<E>,
    > {
        check(cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(native_error(invalid(
                "matching hysteretic spatial field/thermal clocks and complete cells required",
            )));
        }
        let regions = self.regions();
        let network = MagneticThermalNetwork {
            cells: self.cells,
            links: self.links,
            magnetic_regions: &regions,
            volume_tolerance_m3: self.volume_tolerance_m3,
        };
        let problem = field.problem();
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem.cell_volume_m3(cell).ok_or_else(|| {
                    invalid("original hysteretic tetrahedron has no physical volume")
                })
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(native_error)?;
        let conductors = field
            .conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, source)| source.as_ref().map(|source| (cell, source)))
            .collect::<BTreeMap<_, _>>();
        let routes = network
            .routes_for_volumes(
                thermal,
                &volumes,
                &conductors.keys().copied().collect::<Vec<_>>(),
                cancel,
            )
            .map_err(native_error)?;
        network
            .admit_hysteretic_temperature_support(problem.materials(), &routes, cancel)
            .map_err(native_error)?;
        let resolve_materials = |state: &DepositionSimulation<'_>,
                                 sample: &MaterialNetworkSample| {
            network.hysteretic_materials_at(
                state,
                problem.materials(),
                &conductors,
                sample,
                &routes,
                cancel,
            )
        };
        let initial = resolve(thermal)?;
        let (accepted, _) =
            resolve_materials(thermal, snapshot(&initial.source)).map_err(native_error)?;
        if accepted.len() != conductors.len()
            || accepted.iter().any(|(cell, source)| {
                conductors.get(cell).is_none_or(|old| {
                    !old.has_same_electrical_source(source)
                        || !same_material(old.material(), source.material())
                })
            })
        {
            return Err(native_error(invalid(
                "accepted hysteretic spatial conductor sources differ from their thermal state",
            )));
        }
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(cancel).map_err(native_error)?;
                let coefficients = resolve(trial)?;
                let (sources, magnetic) = resolve_materials(trial, snapshot(&coefficients.source))
                    .map_err(native_error)?;
                let mut refreshed = vec![None; volumes.len()];
                for (cell, source) in sources {
                    refreshed[cell] = Some(source);
                }
                let proposal = drive(trial, &coefficients, refreshed, remaining)?;
                let field_internal_energy_j = [
                    proposal.induction_step().stored_before_j,
                    proposal.induction_step().stored_after_j,
                ];
                let circuit_trials = proposal.field_trials();
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials,
                    source: SpatialInductionThermalSource {
                        network: coefficients.source,
                        magnetic: magnetic
                            .into_iter()
                            .map(|sample| SpatialMagneticMaterialSample {
                                tetrahedron: sample.triangle,
                                cell: sample.cell,
                                material: sample.material,
                            })
                            .collect(),
                        fields: Vec::new(),
                        field_internal_energy_j,
                        reversible_heat_j: 0.0,
                    },
                })
            },
            |proposal: &D, dt, checkpoint| {
                let mut powers = deposit_heat(proposal.induction_step(), &routes, dt, checkpoint)
                    .map_err(native_error)?;
                additional_heat(proposal, &mut powers, dt, checkpoint)?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        // Fixed magnetic coefficients give U=psi. Preserve the constitutive
        // stable increment instead of subtracting two large biased totals.
        let external_heat_j = report.non_electrical_heat_j;
        let account = D::receipt_account(&report.electrical);
        let total_energy_residual_j = finite(
            account.stored_change_j
                + report.thermal.stored_energy_change_j
                + account.numerical_dissipation_j
                - account.supplied_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(native_error(invalid(
                "spatial hysteresis, finite enthalpy and external work do not close",
            )));
        }
        check(cancel).map_err(native_error)?;
        Ok((
            thermal,
            field,
            SpatialInductionThermalReport {
                coupled: report,
                external_heat_j,
                total_energy_residual_j,
            },
        ))
    }
}
