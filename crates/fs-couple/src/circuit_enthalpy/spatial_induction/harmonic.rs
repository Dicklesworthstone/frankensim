//! Cycle-averaged induction heating with self-consistent source temperature.
//!
//! Each proposed thermal endpoint resolves its original conductivity sources,
//! then solves the native three-dimensional harmonic A–phi equations. Local
//! cycle-mean Joule powers enter the original finite-cell enthalpy residual,
//! together with the complete heat assembly at that SAME temperature.
//!
//! This is a settled-carrier, slow-thermal model. Electromagnetic transients and
//! envelope magnetic-energy changes are excluded; the thermal balance uses the
//! independently audited steady AC supply. Geometry and material occupancy are
//! stationary and each mapped magnetic permeability source is temperature
//! independent. Use retained transient induction for startup, rapid modulation,
//! magnetic entropy, moving material, saturation or hysteretic losses.

use super::{
    SpatialMagneticMaterialSample, SpatialMagneticRegion, SpatialMagneticThermalNetwork,
    heat::{HeatError, heat_error},
    induction, native_error,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
    contacts::heat::{HeatAssembly, closed_balance},
    finite,
    magnetic::{MagneticThermalError, MagneticThermalNetwork, check, invalid, material},
    material::MaterialNetworkError,
    spatial_induction::SpatialInductionHeatError,
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::induction::{
        SpatialInduction,
        harmonic::{SpatialHarmonicInductionBalance, SpatialHarmonicInductionControls},
    },
};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::spatial::SpatialConductorMaterial, field_storage::admit_source_temperature_refresh,
};
use fs_math::{c64::C64, dd::Dd};
use fs_solver::SolveReport;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Original tetrahedron's actual cycle-mean Joule heat and caloric destination.
#[derive(Clone, Debug)]
pub struct SpatialHarmonicHeatingCell {
    pub tetrahedron: usize,
    pub cell: CellId,
    pub volume_m3: f64,
    pub mean_ohmic_power_w: f64,
    pub heat_j: f64,
}

/// Accepted source-temperature AC response, evaluated freshly after convergence.
///
/// Magnetic storage is the native cycle-mean diagnostic in `balance`; it is
/// not a transient storage increment. Work and heat are steady mean powers
/// times this thermal interval, without a second peak-phasor factor of 1/2.
#[derive(Clone, Debug)]
pub struct SpatialHarmonicHeatingStep {
    pub omega_rad_s: f64,
    pub duration_s: f64,
    pub cells: Vec<SpatialHarmonicHeatingCell>,
    pub balance: SpatialHarmonicInductionBalance,
    pub field: SolveReport,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub maximum_cell_skin_ratio: f64,
    pub heat_j: f64,
    pub supplied_work_j: f64,
}

/// Complete original thermal assembly and source-pinned field material states.
pub struct SpatialHarmonicThermalSample<S> {
    pub network: S,
    pub magnetic: Vec<SpatialMagneticMaterialSample>,
}

/// Unpublished caloric endpoint and its accepted settled AC heating account.
pub struct SpatialHarmonicThermalCandidate<'t, S> {
    thermal: DepositionSimulation<'t>,
    report: CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicHeatingStep>,
}
impl<'t, S> SpatialHarmonicThermalCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicHeatingStep> {
        &self.report
    }
    /// All reported Joule heat is already deposited exactly once.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicHeatingStep>,
    ) {
        (self.thermal, self.report)
    }
}

struct CycleMeanDriver(SpatialHarmonicHeatingStep);
impl ThermalDriver for CycleMeanDriver {
    type Model = ();
    type Step = SpatialHarmonicHeatingStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.0.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.0.supplied_work_j,
        }
    }
    fn into_parts(self) -> ((), Self::Step) {
        ((), self.0)
    }
}

/// Advance finite caloric cells under a self-consistent settled AC excitation.
///
/// `template` supplies the original mesh, gauge, trace and pinned material
/// identities, not a transient magnetic initial condition. It is unchanged and
/// may be reused on consecutive accepted thermal intervals; its clock need not
/// equal the thermal clock. Peak impressed current and angular frequency stay
/// fixed during this interval. All conductors require whole-volume solid
/// caloric receivers. Scalar and oriented conductivity refresh through their
/// original source laws; mapped scalar or oriented permeability must be
/// temperature independent, with its source domain rechecked at every trial.
///
/// Every Newton proposal and the final acceptance re-resolve the complete heat
/// assembly and solve the native field. `circuit_trials` counts those actual
/// field solves and obeys the existing work budget. No AC response matrix is
/// reused after changing conductivity. Contacts, layers, radiation and finite
/// coolant retain their original source payload; numerical closed-radiation
/// net is excluded from external supply. Unowned graph electrical work refuses.
///
/// The caller must justify separation of carrier settling and thermal time
/// scales. This API supplies no automatic time-scale or discretization-error
/// certificate. Thermal rates are backward Euler, so time refinement remains
/// necessary. It resolves neither finite-window phase/beat heating nor the
/// transient energy exchanged while one settled AC state changes to another.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    template: &SpatialInduction<'_, C>,
    impressed_current_a_m2: &[[C64; 3]],
    omega_rad_s: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHarmonicInductionControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialHarmonicThermalCandidate<'t, A::Sample>, HeatError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let dt = finite(target_s - thermal.time_s())?;
    if dt <= 0.0 || !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
        return Err(native_error(invalid(
            "positive thermal interval and AC angular frequency required",
        )));
    }
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    if bulk.cells.len() != thermal.cells().len() {
        return Err(native_error(invalid("complete AC thermal cells required")));
    }
    let layout = SpatialMagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: regions,
        volume_tolerance_m3,
    };
    let regions = layout.regions();
    let network = MagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: &regions,
        volume_tolerance_m3,
    };
    let problem = template.problem();
    let originals = problem
        .materials()
        .iter()
        .map(|m| match m {
            SpatialMagneticMaterial::Linear(m) => Ok(m),
            SpatialMagneticMaterial::Recoil(_) => Err(native_error(invalid(
                "settled AC heating requires original induced-linear magnetic materials",
            ))),
        })
        .collect::<Result<Vec<_>, HeatError<A::Error>>>()?;
    let volumes = (0..problem.mesh().tets.len())
        .map(|cell| {
            check(&cancel).map_err(native_error)?;
            problem.cell_volume_m3(cell).ok_or_else(|| {
                native_error(invalid("original AC tetrahedron has no physical volume"))
            })
        })
        .collect::<Result<Vec<_>, HeatError<A::Error>>>()?;
    let conductors = template
        .conductors()
        .iter()
        .enumerate()
        .filter_map(|(cell, source)| source.as_ref().map(|source| (cell, source)))
        .collect::<BTreeMap<_, _>>();
    let conducting_cells = conductors.keys().copied().collect::<Vec<_>>();
    let routes = network
        .routes_for_volumes(thermal, &volumes, &conducting_cells, &cancel)
        .map_err(native_error)?;
    let resolve = |state: &DepositionSimulation<'_>| {
        check(&cancel).map_err(native_error)?;
        let coefficients = provider
            .resolve_heat(cx, state)
            .map_err(|error| error.map_source(SpatialInductionHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(native_error(invalid(
                "settled AC heating cannot discard graph electrical work",
            )));
        }
        let (sources, magnetic) = network
            .materials_at_using(
                state,
                &originals,
                &conductors,
                &coefficients.bulk.source.as_ref().material,
                &routes,
                &cancel,
                |cell, before, after, binding| {
                    let field_error =
                        |e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e));
                    admit_source_temperature_refresh(before, after).map_err(field_error)?;
                    let law = problem.magnetic_laws()[cell];
                    if law.resolve(before).map_err(field_error)?
                        != law.resolve(after).map_err(field_error)?
                    {
                        return Err(invalid(
                            "settled AC heating cannot change magnetic storage laws",
                        ));
                    }
                    for name in law.source_property_names() {
                        let property = after.property(name).ok_or_else(|| {
                            invalid("settled AC heating lost its original magnetic source")
                        })?;
                        if property
                            .source_axes()
                            .iter()
                            .any(|(axis, _)| axis == binding.point.temperature_axis)
                        {
                            return Err(invalid(
                                "settled AC permeability must be independent of temperature",
                            ));
                        }
                        binding
                            .point
                            .admit_static(property)
                            .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                    }
                    Ok(())
                },
            )
            .map_err(native_error)?;
        Ok((coefficients.bulk, sources, magnetic))
    };
    // Admit the accepted caloric/source correspondence even if a later Newton
    // proposal would happen to return to the template's authored temperature.
    resolve(thermal)?;
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let (thermal, (), mut report) = cell_driver_candidate(
        thermal,
        target_s,
        controls,
        initial_energy_increments_j,
        |trial, _remaining| {
            let (coefficients, sources, magnetic) = resolve(trial)?;
            let mut refreshed = vec![None; originals.len()];
            for (cell, source) in sources {
                refreshed[cell] = Some(source);
            }
            let field = template
                .prepare_materials_trial(
                    magnetic
                        .iter()
                        .map(|s| SpatialMagneticMaterial::Linear(s.material.clone()))
                        .collect(),
                    refreshed,
                    assembly,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(induction)
                .map_err(native_error)?
                .into_harmonic(&mut || (*cancel.borrow_mut())())
                .map_err(induction)
                .map_err(native_error)?;
            let response = field
                .solve_harmonic(
                    impressed_current_a_m2,
                    &velocity,
                    omega_rad_s,
                    field_controls,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(induction)
                .map_err(native_error)?;
            if response.cells().len() != routes.conductors.len() {
                return Err(native_error(invalid(
                    "AC heat omitted original conducting cells",
                )));
            }
            let mut cells = Vec::with_capacity(response.cells().len());
            for loss in response.cells() {
                check(&cancel).map_err(native_error)?;
                let cell = *routes.conductors.get(&loss.cell).ok_or_else(|| {
                    native_error(invalid("AC heat lost its original caloric receiver"))
                })?;
                let heat_j = finite(dt * loss.mean_ohmic_power_w)?;
                if heat_j < 0.0 || loss.mean_ohmic_power_w != 0.0 && heat_j == 0.0 {
                    return Err(native_error(invalid(
                        "AC cell heat is negative or unrepresentable",
                    )));
                }
                cells.push(SpatialHarmonicHeatingCell {
                    tetrahedron: loss.cell,
                    cell,
                    volume_m3: loss.volume_m3,
                    mean_ohmic_power_w: loss.mean_ohmic_power_w,
                    heat_j,
                });
            }
            let balance = *response.balance();
            let step = SpatialHarmonicHeatingStep {
                omega_rad_s,
                duration_s: dt,
                cells,
                balance,
                field: response.report().clone(),
                maximum_ampere_residual_a: response.maximum_ampere_residual_a(),
                maximum_charge_residual_a: response.maximum_charge_residual_a(),
                maximum_cell_skin_ratio: response.maximum_cell_skin_ratio(),
                heat_j: finite(dt * balance.mean_ohmic_power_w)?,
                supplied_work_j: finite(dt * balance.mean_supplied_power_w)?,
            };
            Ok(DrivenCoefficients {
                circuit: CycleMeanDriver(step),
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                circuit_trials: 1,
                source: SpatialHarmonicThermalSample {
                    network: coefficients.source,
                    magnetic: magnetic
                        .into_iter()
                        .map(|s| SpatialMagneticMaterialSample {
                            tetrahedron: s.triangle,
                            cell: s.cell,
                            material: s.material,
                        })
                        .collect(),
                },
            })
        },
        |driver: &CycleMeanDriver, _, checkpoint| {
            let mut powers = BTreeMap::<CellId, Dd>::new();
            for cell in &driver.0.cells {
                checkpoint()?;
                let power = powers.entry(cell.cell).or_insert(Dd::ZERO);
                *power = *power + Dd::from_f64(cell.mean_ohmic_power_w);
            }
            powers
                .into_iter()
                .map(|(id, power)| Ok((id, finite(power.to_f64())?)))
                .collect()
        },
        || (*cancel.borrow_mut())(),
    )?;
    let sample = report.source.network.as_ref();
    closed_balance(
        dt,
        report.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.non_electrical_heat_j = finite(report.non_electrical_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    check(&cancel).map_err(native_error)?;
    Ok(SpatialHarmonicThermalCandidate { thermal, report })
}
