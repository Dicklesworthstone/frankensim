//! Source-driven nonlinear induction and temperature-dependent Ohmic heating.
//!
//! Every thermal trial restarts from the original accepted A_z. The same
//! material card supplies conductivity at that trial's temperature, while the
//! complete source H(B) range is re-admitted as an unchanged constitutive law.
//! The original finite enthalpy owner receives ONLY physical sigma E^2 heat.
//! Backward-Euler magnetic convexity loss stays in the driver work account.

/// Explicitly separable nonlinear magnetic entropy and finite thermal feedback.
pub mod thermal;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients,
    DriverAccount, EddyThermalSample, ThermalDriver, cell_driver_candidate, deposit_cell_heat,
    native_error,
};
use crate::circuit_enthalpy::{
    magnetic::{
        Error, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork, MagneticThermalSample,
        check, invalid, material, same_material,
    },
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use fs_em::planar::{
    eddy::{
        EddyConductor,
        nonlinear::{
            NonlinearEddyControls, NonlinearEddyError, NonlinearEddyStep, PlanarNonlinearEddy,
            PlanarNonlinearEddyCandidate,
        },
    },
    nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
};
use fs_material::{
    conductor::ResolvedConductorState,
    field_storage::{
        MAGNETIC_PERMEABILITY_PROPERTY, admit_fixed_magnetic_refresh,
        nonlinear::REVERSIBLE_FIELD_STRENGTH_PROPERTY, permanent::FLUX_DENSITY_DIMS,
    },
};
use fs_qty::QuantitySpec;
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};

/// Both accepted histories and their exact physical/numerical energy account.
pub struct NonlinearEddyEnthalpyCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: PlanarNonlinearEddy<'m>,
    pub(in crate::circuit_enthalpy) report: CircuitEnthalpyReport<EddyThermalSample<S>, NonlinearEddyStep>,
}
impl<'t, 'm, S> NonlinearEddyEnthalpyCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarNonlinearEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<EddyThermalSample<S>, NonlinearEddyStep> {
        &self.report
    }
    /// Commit both histories together. Reported Joule heat is already deposited;
    /// reported numerical dissipation is never deposited into any material.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarNonlinearEddy<'m>,
        CircuitEnthalpyReport<EddyThermalSample<S>, NonlinearEddyStep>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl<'m> ThermalDriver for PlanarNonlinearEddyCandidate<'m> {
    type Model = PlanarNonlinearEddy<'m>;
    type Step = NonlinearEddyStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            supplied_j: step.supplied_work_j + step.parameter_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        PlanarNonlinearEddyCandidate::into_parts(self)
    }
}
fn induction(error: NonlinearEddyError) -> Error {
    match error {
        NonlinearEddyError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::NonlinearInduction(other)),
    }
}

impl MagneticThermalNetwork<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn nonlinear_materials_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarNonlinearMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        // Fixed H refresh cannot close entropy from changing B, even if T is
        // held and its scale value happens to stay unchanged. This shared gate
        // also protects the original nonlinear winding/enthalpy consumer and
        // explicitly fixed external regions.
        for original in problem.materials() {
            check(cancelled)?;
            if matches!(original, IsotropicMagneticMaterial::Reversible(law)
                if law.field_scale_slope_per_k().is_some_and(|s| s != 0.0)) {
                return Err(invalid("temperature-scaled nonlinear magnetism requires the explicit field-entropy thermal owner"));
            }
        }
        let originals: Vec<_> = problem
            .materials()
            .iter()
            .map(IsotropicMagneticMaterial::material)
            .collect();
        self.materials_at_using(
            thermal,
            &originals,
            conductors,
            sample,
            routes,
            cancelled,
            |triangle, original, snapshot, binding| {
                let (property_name, additional) = match &problem.materials()[triangle] {
                    IsotropicMagneticMaterial::Linear(_) => {
                        admit_fixed_magnetic_refresh(original, snapshot).map_err(|e| {
                            CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e))
                        })?;
                        (MAGNETIC_PERMEABILITY_PROPERTY, None)
                    }
                    IsotropicMagneticMaterial::Reversible(law) => {
                        law.admit_fixed_material_refresh(snapshot).map_err(|e| {
                            CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e))
                        })?;
                        (
                            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
                            Some(law.source_path().axis.as_str()),
                        )
                    }
                    IsotropicMagneticMaterial::Recoil(_) => {
                        return Err(invalid(
                            "magnetothermal recoil refresh requires remanence work and entropy data",
                        ));
                    }
                };
                let property = snapshot
                    .property(property_name)
                    .ok_or_else(|| invalid("refreshed nonlinear magnetic property is missing"))?;
                // B=0 is the retained constitutive reference coordinate. The
                // original nonlinear field owner evolves B and admits the full
                // source path; electrical transport remains Ohmic at fixed B.
                let axes: Vec<_> = additional
                    .into_iter()
                    .map(|axis| (axis, QuantitySpec::dimensional(FLUX_DENSITY_DIMS)))
                    .collect();
                binding
                    .point
                    .admit_static_with_axes(property, &axes)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                Ok(())
            },
        )
    }

    /// Evolve nonlinear induction and the SAME finite thermal cells together.
    /// Conductivity is re-resolved at the thermal endpoint; every field trial
    /// starts from the same accepted A/current history. Physical Joule heat,
    /// caloric storage, conduction and authored boundary sources use the shared
    /// original thermal residual and candidate publication path.
    ///
    /// The law H(B), geometry and caloric charts remain fixed. A temperature-
    /// dependent H(B) source refuses thermal refresh even when endpoint values
    /// coincide: magnetic entropy/parameter work needs its own coupled model.
    /// First-order backward Euler numerical dissipation is reported separately
    /// in `report.electrical`; it participates in global energy admission but
    /// never heats the cells. No empirical hysteresis/core-loss map is added.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_eddy_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarNonlinearEddy<'m>,
        impressed_jz_a_m2: &[f64],
        controls: CircuitEnthalpySolve,
        field_controls: NonlinearEddyControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearEddyEnthalpyCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal, controls, initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_nonlinear_eddy_resolved_candidate(
            thermal, target_s, field, impressed_jz_a_m2, controls, field_controls,
            initial_energy_increments_j,
            |trial| prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(material),
            |sample| sample, &cancel,
        )
    }

    // The original fixed-law backward-Euler field step owns storage and
    // numerical loss; only its complete thermal coefficient source is generic.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn step_nonlinear_eddy_resolved_candidate<
        't, 'm, S, E: From<MagneticThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarNonlinearEddy<'m>,
        impressed_jz_a_m2: &[f64],
        controls: CircuitEnthalpySolve,
        field_controls: NonlinearEddyControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'_>) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<NonlinearEddyEnthalpyCandidate<'t, 'm, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal, controls, initial_energy_increments_j,
        )?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(native_error(invalid(
                "matching nonlinear magnetic/thermal clocks and complete cell bindings required",
            )));
        }
        let triangles: Vec<_> = field.conductors().iter().map(|s| s.triangle).collect();
        let routes = self.routes_for_mesh(thermal, field.problem().mesh(), &triangles, cancel)
            .map_err(native_error)?;
        let conductors: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, &s.conductor))
            .collect();
        let regions: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, s.region))
            .collect();
        let initial = resolve(thermal)?;
        let (sources, _) = self.nonlinear_materials_at(
            thermal,
            field.problem(),
            &conductors,
            snapshot(&initial.source),
            &routes,
            cancel,
        ).map_err(native_error)?;
        if sources.len() != field.conductors().len()
            || sources
                .iter()
                .zip(field.conductors())
                .any(|((triangle, conductor), original)| {
                    *triangle != original.triangle
                        || conductor.source_basis() != original.conductor.source_basis()
                        || conductor.source_key() != original.conductor.source_key()
                        || !same_material(conductor.material(), original.conductor.material())
                })
        {
            return Err(native_error(invalid(
                "accepted nonlinear conductor sources must match the actual thermal state",
            )));
        }
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _| {
                check(cancel).map_err(native_error)?;
                let coefficients = resolve(state)?;
                let (sources, magnetic) = self.nonlinear_materials_at(
                    state,
                    field.problem(),
                    &conductors,
                    snapshot(&coefficients.source),
                    &routes,
                    cancel,
                ).map_err(native_error)?;
                let sources = sources
                    .into_iter()
                    .map(|(triangle, conductor)| EddyConductor {
                        triangle,
                        region: regions[&triangle],
                        conductor,
                    })
                    .collect();
                let proposal = field
                    .step_to_candidate(
                        target_s,
                        sources,
                        impressed_jz_a_m2,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| native_error(induction(error)))?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: EddyThermalSample {
                        network: coefficients.source,
                        magnetic,
                    },
                    circuit_trials: 1,
                })
            },
            |proposal: &PlanarNonlinearEddyCandidate<'_>, dt, checkpoint| {
                deposit_cell_heat(&proposal.step().cells, &routes, dt, checkpoint)
                    .map_err(native_error)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(native_error)?;
        Ok(NonlinearEddyEnthalpyCandidate {
            thermal,
            field,
            report,
        })
    }
}
