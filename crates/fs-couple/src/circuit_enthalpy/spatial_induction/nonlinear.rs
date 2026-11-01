//! Retained three-dimensional saturation, Ohmic heat and finite enthalpy.
//!
//! The original nonlinear A–phi owner retains its field and current history.
//! The same thermal residual refreshes conductivity and, with explicit caloric
//! laws, nonlinear magnetic storage from the original material sources. Native
//! Joule heat and reversible entropy exchange enter their finite receivers;
//! backward-Euler magnetic loss stays a separate numerical energy export.

/// Source-separable nonlinear magnetism in the original spatial heat solve.
pub mod thermal;
/// Native distributed winding circuits on the same finite thermal residual.
pub mod winding;

use super::{
    CircuitEnthalpyError, CircuitEnthalpySolve, CircuitThermalCoefficients, DrivenCoefficients,
    DriverAccount, Error, MagneticThermalNetwork, MaterialNetwork, MaterialNetworkSample, RefCell,
    SpatialInductionThermalReport, SpatialInductionThermalSource, SpatialMagneticCaloricSample,
    SpatialMagneticMaterialSample, SpatialMagneticThermalLaw, SpatialMagneticThermalNetwork,
    ThermalDriver, cell_driver_candidate, check, deposit_cell_heat, finite, invalid, material,
    native_error, same_material, MagneticThermalError,
};
use crate::circuit_enthalpy::magnetic::caloric::{
    MagneticCaloricExchange, MagneticConstitutiveThermalLaw, MagneticEnergyCell,
    MagneticThermalLaw, PreparedMagneticCaloric, ResolvedMagneticConstitutiveThermalLaw,
};
use fs_em::spatial::{
    SpatialMagneticError,
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        induction::nonlinear::{
            SpatialNonlinearInduction, SpatialNonlinearInductionCandidate,
            SpatialNonlinearInductionControls, SpatialNonlinearInductionError,
            SpatialNonlinearInductionStep,
        },
        nonlinear::NonlinearVectorError,
    },
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_material::field_storage::{
    admit_magnetic_temperature_refresh, permanent::FLUX_DENSITY_DIMS,
};
use fs_qty::QuantitySpec;
use fs_thermal::enthalpy::DepositionSimulation;
use fs_thermal::enthalpy::{CellId, DepositionError};
use std::collections::BTreeMap;

/// One jointly accepted magnetic/conducting state and its finite heat history.
pub struct SpatialNonlinearInductionThermalCandidate<
    't, 'm, S = MaterialNetworkSample, C: SpatialConductorMaterial = ResolvedConductorState,
> {
    thermal: DepositionSimulation<'t>,
    field: SpatialNonlinearInduction<'m, C>,
    pub(super) report: SpatialInductionThermalReport<SpatialNonlinearInductionStep, S>,
}
impl<'t, 'm, S, C: SpatialConductorMaterial>
    SpatialNonlinearInductionThermalCandidate<'t, 'm, S, C>
{
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &SpatialNonlinearInduction<'m, C> {
        &self.field
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<SpatialNonlinearInductionStep, S> {
        &self.report
    }
    /// Publish both histories together. Reported Joule heat is already paid;
    /// numerical magnetic dissipation must never be deposited a second time.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialNonlinearInduction<'m, C>,
        SpatialInductionThermalReport<SpatialNonlinearInductionStep, S>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver for SpatialNonlinearInductionCandidate<'m, C> {
    type Model = SpatialNonlinearInduction<'m, C>;
    type Step = SpatialNonlinearInductionStep;
    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialNonlinearInductionCandidate::into_parts(self)
    }
}

// A native terminal circuit changes the external driver, not the material
// query, entropy law, finite thermal residual or local field heat integral.
trait InductionThermalDriver: ThermalDriver {
    fn induction_step(&self) -> &SpatialNonlinearInductionStep;
    fn induction_receipt(step: &Self::Step) -> &SpatialNonlinearInductionStep;
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_>;
    fn field_trials(&self) -> usize;
}

impl<C: SpatialConductorMaterial> InductionThermalDriver
    for SpatialNonlinearInductionCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialNonlinearInductionStep {
        self.step()
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialNonlinearInductionStep {
        step
    }
    fn receipt_account(s: &Self::Step) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: s.heat_j,
            numerical_dissipation_j: s.numerical_dissipation_j,
            stored_change_j: s.stored_change_j - s.parameter_work_j,
            supplied_j: s.impressed_work_j - s.mechanical_work_j,
        }
    }
    fn field_trials(&self) -> usize {
        1
    }
}

fn induction(error: SpatialNonlinearInductionError) -> Error {
    match error {
        SpatialNonlinearInductionError::Spatial(SpatialMagneticError::Cancelled)
        | SpatialNonlinearInductionError::Magnetic(NonlinearVectorError::Spatial(
            SpatialMagneticError::Cancelled,
        )) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(
            crate::circuit_enthalpy::magnetic::MagneticThermalError::SpatialNonlinearInduction(
                other,
            ),
        ),
    }
}

impl SpatialMagneticThermalNetwork<'_> {
    /// Close 3D saturation, local source conductivity and finite thermal energy.
    ///
    /// Every conducting tetrahedron has an explicit whole-volume caloric
    /// receiver. Its original temperature-dependent conductivity/resistivity
    /// is re-evaluated at the joint endpoint; the native field supplies the
    /// exact integrated local current and E_conductor dot sigma E_conductor heat. Heat links,
    /// reservoirs, independent sources, budgets and Newton iteration remain
    /// with the same shared thermal owner as linear induction and circuits.
    ///
    /// Original magnetic H(B), linear and recoil laws must be temperature
    /// independent along this trajectory. Their existing source admission
    /// checks the complete field range and operating envelope; coincident
    /// endpoint coefficients cannot hide a temperature-dependent law. For
    /// these fixed laws U_field=psi, so there is no additional magnetic entropy
    /// source. Temperature-dependent magnetic storage belongs to the separate
    /// caloric-derivative owner, not an invented loss term.
    ///
    /// Geometry and conducting occupancy stay fixed and fully solid. Held
    /// Eulerian velocity transfers signed mechanical work to an external
    /// receiver; it does not evolve mechanics. The first-order implicit scheme
    /// reports its positive numerical magnetic loss but never heats the cells
    /// with it. Every trial restarts from the same accepted A and thermal
    /// history; cancellation, source refusal or any failed balance publishes
    /// neither history. Source/event changes require interval boundaries.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, MaterialNetworkSample, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let material_network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = material_network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_nonlinear_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
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

    // Full sourced heat assemblies keep these original field and enthalpy
    // equations and carry their complete interface/radiation source payload.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_nonlinear_resolved_candidate<'t, 'm, S, E: From<MagneticThermalError>, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, S, C>, CircuitEnthalpyError<E>> {
        self.step_constitutive_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            None,
            controls,
            field_controls,
            initial_energy_increments_j,
            resolve_sources,
            sample,
            cancel,
        )
    }

    // Fixed and temperature-dependent laws share every routing, field, heat,
    // cancellation and acceptance operation. Only constitutive admission and
    // its explicitly owned entropy exchange differ.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_constitutive_resolved_candidate<'t, 'm, S, E: From<MagneticThermalError>, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        bindings: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, S, C>, CircuitEnthalpyError<E>> {
        let (thermal, field, report) = self.step_constitutive_driver_candidate(
            thermal,
            target_s,
            field,
            bindings,
            controls,
            initial_energy_increments_j,
            resolve_sources,
            sample,
            |_, _, materials, conductors, _| {
                let dt = target_s - thermal.time_s();
                if let Some(materials) = materials {
                    field.step_with_materials_candidate(
                        endpoint_current_a_m2,
                        velocity_m_s,
                        dt,
                        materials,
                        conductors,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                } else {
                    field.step_with_conductors_candidate(
                        endpoint_current_a_m2,
                        velocity_m_s,
                        dt,
                        conductors,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                }
                .map_err(induction)
                .map_err(native_error)
            },
            |_, _, _, _| Ok(()),
            cancel,
        )?;
        Ok(SpatialNonlinearInductionThermalCandidate {
            thermal,
            field,
            report,
        })
    }

    // The public current drive and native winding terminals use identical
    // original material/caloric preparation and one finite-enthalpy Newton
    // solve. Only the candidate driver and its additional native heat differ.
    #[allow(clippy::too_many_arguments)]
    fn step_constitutive_driver_candidate<
        't, 'm, S, E: From<MagneticThermalError>,
        C: SpatialConductorMaterial, D: InductionThermalDriver,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m, C>,
        bindings: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        drive: impl Fn(
            &DepositionSimulation<'_>,
            &MaterialNetworkSample,
            Option<Vec<SpatialNonlinearMaterial>>,
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
        (DepositionSimulation<'t>, D::Model, SpatialInductionThermalReport<D::Step, S>),
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
                "matching nonlinear spatial field/thermal clocks and complete cells required",
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
        for original in problem.materials() {
            check(cancel).map_err(native_error)?;
            if bindings.is_none()
                && matches!(original, SpatialNonlinearMaterial::Reversible(law)
                if law.field_scale_slope_per_k().is_some_and(|s| s != 0.0))
            {
                return Err(native_error(invalid(
                    "temperature-scaled nonlinear magnetism requires an explicit spatial field-entropy receiver",
                )));
            }
        }
        let originals = problem
            .materials()
            .iter()
            .map(SpatialNonlinearMaterial::material)
            .collect::<Vec<_>>();
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(&cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("original nonlinear tetrahedron has no physical volume"))
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(native_error)?;
        let conductors = field
            .conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, source)| source.as_ref().map(|s| (cell, s)))
            .collect::<BTreeMap<_, _>>();
        let routes = network.routes_for_volumes(
            thermal,
            &volumes,
            &conductors.keys().copied().collect::<Vec<_>>(),
            &cancel,
        )
        .map_err(native_error)?;
        let laws = bindings.map(|bindings| {
            bindings
                .iter()
                .map(|b| MagneticThermalLaw {
                    triangle: b.tetrahedron,
                    law: b.law.clone(),
                    maximum_flux_density_t: b.maximum_flux_density_t,
                })
                .collect::<Vec<_>>()
        });
        let caloric = if let Some(laws) = &laws {
            for (&cell, &receiver) in &routes.magnetic {
                check(cancel).map_err(native_error)?;
                if matches!(
                    problem.materials()[cell],
                    SpatialNonlinearMaterial::Recoil(_)
                ) {
                    return Err(native_error(invalid(
                        "spatial recoil temperature response requires its own remanence entropy law",
                    )));
                }
                if receiver.is_none()
                    && matches!(&problem.materials()[cell],
                    SpatialNonlinearMaterial::Reversible(law)
                        if law.field_scale_slope_per_k().is_some_and(|s| s != 0.0))
                {
                    return Err(native_error(invalid(
                        "temperature-scaled spatial magnetic region requires a finite entropy receiver",
                    )));
                }
            }
            let native = field
                .field(&mut || (*cancel.borrow_mut())())
                .map_err(induction)
                .map_err(native_error)?;
            Some(PreparedMagneticCaloric::new_for_volumes(
                &network,
                thermal,
                &volumes,
                &originals,
                native.flux_density_t(),
                laws,
                &routes,
                cancel,
            )
            .map_err(native_error)?)
        } else {
            None
        };
        let resolve = |state: &DepositionSimulation<'_>,
                       sample: &MaterialNetworkSample|
         -> Result<_, Error> {
            if let Some(laws) = &laws {
                let error = |e| {
                    CircuitEnthalpyError::Source(
                        crate::circuit_enthalpy::magnetic::MagneticThermalError::Magnetic(e),
                    )
                };
                let (sources, samples, resolved) = network.constitutive_materials_at(
                    state,
                    &originals,
                    &conductors,
                    sample,
                    laws,
                    &routes,
                    cancel,
                    |cell, at| match (&problem.materials()[cell], at) {
                        (
                            SpatialNonlinearMaterial::Linear(original),
                            ResolvedMagneticConstitutiveThermalLaw::Linear(at),
                        ) => admit_magnetic_temperature_refresh(original, at.material())
                            .map_err(error),
                        (
                            SpatialNonlinearMaterial::Reversible(original),
                            ResolvedMagneticConstitutiveThermalLaw::Separable(at),
                        ) => original
                            .admit_temperature_refresh(at.magnetic())
                            .map_err(error),
                        _ => Err(invalid(
                            "spatial thermal law must match its original magnetic material variant",
                        )),
                    },
                )?;
                let materials = problem
                    .materials()
                    .iter()
                    .enumerate()
                    .map(|(cell, original)| match resolved.get(&cell) {
                        Some(ResolvedMagneticConstitutiveThermalLaw::Linear(at)) => {
                            SpatialNonlinearMaterial::Linear(at.material().clone())
                        }
                        Some(ResolvedMagneticConstitutiveThermalLaw::Separable(at)) => {
                            SpatialNonlinearMaterial::Reversible(at.magnetic().clone())
                        }
                        None => original.clone(),
                    })
                    .collect::<Vec<_>>();
                return Ok((sources, samples, Some(materials)));
            }
            let (sources, samples) = network.materials_at_using(
                state,
                &originals,
                &conductors,
                sample,
                &routes,
                &cancel,
                |cell, original, snapshot, binding| {
                    problem.materials()[cell]
                        .admit_fixed_material_refresh(snapshot)
                        .map_err(|e| {
                            CircuitEnthalpyError::Source(
                                crate::circuit_enthalpy::magnetic::MagneticThermalError::Magnetic(
                                    e,
                                ),
                            )
                        })?;
                    let axes = match &problem.materials()[cell] {
                        SpatialNonlinearMaterial::Reversible(law) => vec![(
                            law.source_path().axis.as_str(),
                            QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                        )],
                        _ => Vec::new(),
                    };
                    for property in original.properties() {
                        check(&cancel)?;
                        let refreshed = snapshot
                            .properties()
                            .iter()
                            .find(|p| p.requirement() == property.requirement())
                            .ok_or_else(|| {
                                invalid("refreshed spatial magnetic source is missing")
                            })?;
                        binding
                            .point
                            .admit_static_with_axes(refreshed, &axes)
                            .map_err(|e| {
                                material(super::super::material::MaterialNetworkError::Material(e))
                            })?;
                    }
                    Ok(())
                },
            )?;
            Ok((sources, samples, None))
        };
        let initial = resolve_sources(thermal)?;
        let (accepted, samples, _) = resolve(thermal, sample(&initial.source))
            .map_err(native_error)?;
        if let Some(caloric) = &caloric {
            caloric.validate_sources(thermal, &samples, cancel)
                .map_err(native_error)?;
        }
        if accepted.len() != conductors.len()
            || accepted.iter().any(|(cell, source)| {
                conductors.get(cell).is_none_or(|old| {
                    !old.has_same_electrical_source(source)
                        || !same_material(old.material(), source.material())
                })
            })
        {
            return Err(native_error(invalid(
                "accepted nonlinear spatial conductor sources differ from their thermal state",
            )));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel).map_err(native_error)?;
                let mut coefficients = resolve_sources(trial)?;
                let (sources, magnetic, materials) = resolve(trial, sample(&coefficients.source))
                    .map_err(native_error)?;
                if let Some(caloric) = &caloric {
                    caloric.validate_sources(trial, &magnetic, cancel)
                        .map_err(native_error)?;
                }
                let mut refreshed = vec![None; originals.len()];
                for (cell, source) in sources {
                    refreshed[cell] = Some(source);
                }
                let proposal = drive(
                    trial,
                    sample(&coefficients.source),
                    materials,
                    refreshed,
                    remaining,
                )?;
                let exchange = if let Some(caloric) = &caloric {
                    caloric.deposit(
                        trial,
                        proposal
                            .induction_step()
                            .magnetic_cells
                            .iter()
                            .map(|cell| MagneticEnergyCell {
                                triangle: cell.cell,
                                flux_density_t: cell.flux_density_t,
                                free_energy_j: cell.free_energy_j,
                                parameter_work_j: cell.parameter_work_j,
                            }),
                        &mut coefficients.sources_w,
                        dt,
                        controls.total_energy_tolerance_j,
                        cancel,
                    )
                    .map_err(native_error)?
                } else {
                    MagneticCaloricExchange {
                        fields: Vec::new(),
                        field_internal_energy_j: [
                            proposal.induction_step().stored_before_j,
                            proposal.induction_step().stored_after_j,
                        ],
                        reversible_heat_j: 0.0,
                    }
                };
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
                            .map(|s| SpatialMagneticMaterialSample {
                                tetrahedron: s.triangle,
                                cell: s.cell,
                                material: s.material,
                            })
                            .collect(),
                        fields: exchange
                            .fields
                            .into_iter()
                            .map(|s| SpatialMagneticCaloricSample {
                                tetrahedron: s.triangle,
                                cell: s.cell,
                                temperature_k: s.temperature_k,
                                flux_density_t: s.flux_density_t,
                                state: s.state,
                                parameter_work_j: s.parameter_work_j,
                                reversible_heat_j: s.reversible_heat_j,
                            })
                            .collect(),
                        field_internal_energy_j: exchange.field_internal_energy_j,
                        reversible_heat_j: exchange.reversible_heat_j,
                    },
                })
            },
            |proposal: &D, dt, checkpoint| {
                let mut powers = deposit_cell_heat(
                    &proposal.induction_step().cells, &routes, dt, checkpoint,
                ).map_err(native_error)?;
                additional_heat(proposal, &mut powers, dt, checkpoint)?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j =
            finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let induction = D::induction_receipt(&report.electrical);
        let account = D::receipt_account(&report.electrical);
        // Substitute U_field for the field's free-energy contribution while
        // retaining every native circuit/leakage store and external work term.
        let total_energy_residual_j = finite(
            report.source.field_internal_energy_j[1] - report.source.field_internal_energy_j[0]
                + account.stored_change_j - induction.stored_change_j
                + induction.parameter_work_j
                + report.thermal.stored_energy_change_j
                + account.numerical_dissipation_j
                - account.supplied_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(native_error(invalid(
                "spatial nonlinear field, finite enthalpy and external work do not close",
            )));
        }
        check(&cancel).map_err(native_error)?;
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
