//! Source-separable H(B,T), actual eddy currents and finite caloric receivers.
//!
//! H=s(T)H_ref(B) is an explicit constitutive declaration, not an inferred
//! general surface. The original field solver owns A, nonlinear reluctivity,
//! Joule heat, parameter work and numerical loss. The shared magnetic caloric
//! owner supplies entropy and routes -P-delta(T*S) exactly once.

use super::*;
use crate::circuit_enthalpy::{
    CircuitThermalCoefficients,
    eddy::thermal::{EddyCaloricReport, EddyCaloricSource},
    finite,
    magnetic::caloric::{
        MagneticConstitutiveThermalLaw, MagneticEnergyCell, MagneticThermalLaw,
        PreparedMagneticCaloric,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::field_storage::admit_magnetic_temperature_refresh;

/// One atomic accepted nonlinear field/thermal history and its physical audit.
pub struct NonlinearEddyCaloricCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: PlanarNonlinearEddy<'m>,
    pub(in crate::circuit_enthalpy) report: EddyCaloricReport<S, NonlinearEddyStep>,
}
impl<'t, 'm, S> NonlinearEddyCaloricCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarNonlinearEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &EddyCaloricReport<S, NonlinearEddyStep> {
        &self.report
    }
    /// Publish both histories together. Physical heat has already been
    /// deposited; backward-Euler numerical loss never enters a thermal cell.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarNonlinearEddy<'m>,
        EddyCaloricReport<S, NonlinearEddyStep>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl MagneticThermalNetwork<'_> {
    #[allow(clippy::too_many_arguments)]
    fn separable_materials_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarNonlinearMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
            Vec<IsotropicMagneticMaterial>,
        ),
        Error,
    > {
        let mut resolved = BTreeMap::new();
        for binding in bindings {
            check(cancel)?;
            let cell = routes.magnetic[&binding.triangle].expect("admitted thermal law");
            let temperature = thermal.temperature_k(cell)?;
            let error = |e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e));
            let at = match (&binding.law, &problem.materials()[binding.triangle]) {
                (
                    MagneticConstitutiveThermalLaw::Linear(law),
                    IsotropicMagneticMaterial::Linear(original),
                ) => {
                    let at = law.at_temperature(temperature).map_err(error)?;
                    admit_magnetic_temperature_refresh(original, at.material()).map_err(error)?;
                    IsotropicMagneticMaterial::Linear(at.material().clone())
                }
                (
                    MagneticConstitutiveThermalLaw::Separable(law),
                    IsotropicMagneticMaterial::Reversible(original),
                ) => {
                    let at = law.at_temperature(temperature).map_err(error)?;
                    original
                        .admit_temperature_refresh(at.magnetic())
                        .map_err(error)?;
                    IsotropicMagneticMaterial::Reversible(at.magnetic().clone())
                }
                _ => {
                    return Err(invalid(
                        "thermal magnetic law must match its original linear or explicitly source-scaled nonlinear material",
                    ));
                }
            };
            resolved.insert(binding.triangle, at);
        }
        let originals: Vec<_> = problem
            .materials()
            .iter()
            .map(IsotropicMagneticMaterial::material)
            .collect();
        let (sources, samples) = self.materials_at_using(
            thermal, &originals, conductors, sample, routes, cancel,
            |triangle, _, snapshot, binding| {
                let at = &resolved[&triangle];
                if !same_material(snapshot, at.material()) {
                    return Err(invalid("nonlinear field and conductor must retain the complete original separable source bundle"));
                }
                let (property_name, axes) = match at {
                    IsotropicMagneticMaterial::Linear(_) => (MAGNETIC_PERMEABILITY_PROPERTY, Vec::new()),
                    IsotropicMagneticMaterial::Reversible(law) => (REVERSIBLE_FIELD_STRENGTH_PROPERTY,
                        vec![(law.source_path().axis.as_str(), QuantitySpec::dimensional(FLUX_DENSITY_DIMS))]),
                    IsotropicMagneticMaterial::Recoil(_) => unreachable!("admitted thermal variants"),
                };
                let property = snapshot.property(property_name)
                    .ok_or_else(|| invalid("missing original nonlinear reference source"))?;
                binding.point.admit_static_with_axes(property, &axes)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                Ok(())
            },
        )?;
        let materials = problem
            .materials()
            .iter()
            .enumerate()
            .map(|(triangle, original)| resolved.get(&triangle).unwrap_or(original).clone())
            .collect();
        Ok((sources, samples, materials))
    }

    /// Evolve mixed linear and explicitly source-separable nonlinear magnetic
    /// materials and the same original finite heat cells. Every thermal region needs its complete
    /// original law/bundle and field envelope; initial field materials must be
    /// constructed with that law's `at_temperature(...).magnetic()` result.
    /// H_ref is fixed and monotone; its separate positive scale is affine in
    /// absolute temperature. Density, solid geometry and all source pins stay
    /// fixed. General H(B,T), hysteresis, phase change and magnetostriction are
    /// outside this model. First-order field/thermal time refinement is needed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_eddy_thermal_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarNonlinearEddy<'m>,
        impressed_jz_a_m2: &[f64],
        bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: NonlinearEddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearEddyCaloricCandidate<'t, 'm>, Error> {
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
        self.step_nonlinear_eddy_thermal_resolved_candidate(
            thermal,
            target_s,
            field,
            impressed_jz_a_m2,
            bindings,
            controls,
            field_controls,
            assembly,
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

    /// Shared residual for ordinary conduction or the complete prepared heat
    /// assembly. The callback supplies original heat-owner coefficients only.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn step_nonlinear_eddy_thermal_resolved_candidate<'t, 'm, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarNonlinearEddy<'m>,
        impressed_jz_a_m2: &[f64],
        bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: NonlinearEddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'_>) -> Result<CircuitThermalCoefficients<S>, Error>,
        snapshot: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<NonlinearEddyCaloricCandidate<'t, 'm, S>, Error> {
        check(cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching nonlinear field/thermal clocks and complete original cells required",
            ));
        }
        let triangles: Vec<_> = field.conductors().iter().map(|s| s.triangle).collect();
        let routes = self.routes_for_mesh(thermal, field.problem().mesh(), &triangles, cancel)?;
        // A fixed-temperature region can still exchange entropy as B changes.
        // It cannot silently claim that transfer without a physical receiver.
        for (&triangle, &receiver) in &routes.magnetic {
            check(cancel)?;
            if matches!(
                field.problem().materials()[triangle],
                IsotropicMagneticMaterial::Recoil(_)
            ) {
                return Err(invalid(
                    "recoil magnetothermal coupling requires its own source remanence entropy law",
                ));
            }
            if receiver.is_none()
                && matches!(&field.problem().materials()[triangle],
                IsotropicMagneticMaterial::Reversible(law) if law.field_scale_slope_per_k().is_some_and(|s| s != 0.0))
            {
                return Err(invalid(
                    "temperature-scaled nonlinear magnetic region needs a finite entropy receiver even at fixed temperature",
                ));
            }
        }
        let native = field
            .field(&mut || (*cancel.borrow_mut())())
            .map_err(induction)?;
        let originals: Vec<_> = field
            .problem()
            .materials()
            .iter()
            .map(IsotropicMagneticMaterial::material)
            .collect();
        let caloric = PreparedMagneticCaloric::new(
            self,
            thermal,
            field.problem().mesh(),
            &originals,
            native.magnetic_flux_density_t(),
            bindings,
            &routes,
            cancel,
        )?;
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
        let (accepted, samples, _) = self.separable_materials_at(
            thermal,
            field.problem(),
            &conductors,
            snapshot(&initial.source),
            bindings,
            &routes,
            cancel,
        )?;
        caloric.validate_sources(thermal, &samples, cancel)?;
        if accepted.len() != field.conductors().len()
            || accepted
                .iter()
                .zip(field.conductors())
                .any(|((triangle, source), original)| {
                    *triangle != original.triangle
                        || source.source_basis() != original.conductor.source_basis()
                        || source.source_key() != original.conductor.source_key()
                        || !same_material(source.material(), original.conductor.material())
                })
        {
            return Err(invalid(
                "accepted nonlinear conductor sources must match the actual thermal state",
            ));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _| {
                check(cancel)?;
                let mut coefficients = resolve(trial)?;
                let (sources, magnetic, materials) = self.separable_materials_at(
                    trial,
                    field.problem(),
                    &conductors,
                    snapshot(&coefficients.source),
                    bindings,
                    &routes,
                    cancel,
                )?;
                caloric.validate_sources(trial, &magnetic, cancel)?;
                let proposal = field
                    .step_with_materials_candidate(
                        target_s,
                        materials,
                        sources
                            .into_iter()
                            .map(|(triangle, conductor)| EddyConductor {
                                triangle,
                                region: regions[&triangle],
                                conductor,
                            })
                            .collect(),
                        impressed_jz_a_m2,
                        assembly,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)?;
                let exchange = caloric.deposit(
                    trial,
                    proposal
                        .step()
                        .magnetic_cells
                        .iter()
                        .map(|cell| MagneticEnergyCell {
                            triangle: cell.triangle,
                            flux_density_t: cell.flux_density_t,
                            free_energy_j: cell.free_energy_j,
                            parameter_work_j: cell.parameter_work_j,
                        }),
                    &mut coefficients.sources_w,
                    dt,
                    controls.total_energy_tolerance_j,
                    cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials: 1,
                    source: EddyCaloricSource {
                        material: EddyThermalSample {
                            network: coefficients.source,
                            magnetic,
                        },
                        fields: exchange.fields,
                        field_internal_energy_j: exchange.field_internal_energy_j,
                        reversible_heat_j: exchange.reversible_heat_j,
                    },
                })
            },
            |proposal: &PlanarNonlinearEddyCandidate<'_>, dt, checkpoint| {
                deposit_cell_heat(&proposal.step().cells, &routes, dt, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j =
            finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(
            report.source.field_internal_energy_j[1] - report.source.field_internal_energy_j[0]
                + report.thermal.stored_energy_change_j
                + report.electrical.numerical_dissipation_j
                - report.electrical.supplied_work_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(invalid(
                "nonlinear magnetic internal energy, caloric storage and separate numerical loss do not close external supply",
            ));
        }
        check(cancel)?;
        Ok(NonlinearEddyCaloricCandidate {
            thermal,
            field,
            report: EddyCaloricReport {
                coupled: report,
                external_heat_j,
                total_energy_residual_j,
            },
        })
    }
}
