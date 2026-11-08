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

/// Source-bound loaded circuits and conductor/resistor thermal feedback.
pub mod circuit;
/// Physical voltage/current electrodes on the same settled thermal residual.
pub mod terminals;

mod field;

use super::{
    SpatialMagneticMaterialSample, SpatialMagneticRegion, SpatialMagneticThermalNetwork,
    heat::{HeatError, heat_error},
    induction, native_error,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    admit_driver_unknowns, cell_driver_candidate,
    contacts::heat::{HeatAssembly, closed_balance},
    finite,
    magnetic::{MagneticThermalNetwork, check, invalid},
    spatial_induction::SpatialInductionHeatError,
};
use field::{HarmonicFieldMaterials, recover_field_heat};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::induction::{
        SpatialInduction,
        harmonic::{
            SpatialHarmonicInductionBalance, SpatialHarmonicInductionControls,
            SpatialHarmonicInductionTerminalTransfer,
        },
    },
};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_math::{c64::C64, dd::Dd};
use fs_solver::SolveReport;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};
use terminals::HarmonicHeatingDrive;

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
    /// Actual reactions in original electrode order; empty for isolated drive.
    /// Currents are positive into the conductor and voltages use its reference.
    pub terminals: Vec<SpatialHarmonicInductionTerminalTransfer>,
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
    heat_assembly_with_drive(
        provider,
        cx,
        thermal,
        target_s,
        HarmonicHeatingDrive::Isolated(template),
        impressed_current_a_m2,
        omega_rad_s,
        regions,
        volume_tolerance_m3,
        controls,
        field_controls,
        assembly,
        initial_energy_increments_j,
    )
}

// Electrode causality changes only native field preparation and excitation.
// Material admission, physical cell heat, the complete heat assembly, nonlinear
// thermal residual and acceptance all remain in this one original owner.
#[allow(clippy::too_many_arguments)]
fn heat_assembly_with_drive<'t, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    template: HarmonicHeatingDrive<'_, '_, C>,
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
    let field_materials = HarmonicFieldMaterials::new(problem, template.conductors(), &cancel)
        .map_err(native_error)?;
    let routes = field_materials
        .routes(&network, thermal, &cancel)
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
        let (sources, magnetic) = field_materials
            .resolve(
                &network,
                state,
                &coefficients.bulk.source.as_ref().material,
                &routes,
                &cancel,
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
            let mut refreshed = vec![None; field_materials.cell_count()];
            for (cell, source) in sources {
                refreshed[cell] = Some(source);
            }
            let step = template
                .with_materials(
                    magnetic
                        .iter()
                        .map(|s| SpatialMagneticMaterial::Linear(s.material.clone()))
                        .collect(),
                    refreshed,
                    impressed_current_a_m2,
                    &velocity,
                    omega_rad_s,
                    field_controls,
                    assembly,
                    &mut || (*cancel.borrow_mut())(),
                    |response| recover_field_heat(response, dt, &routes, &cancel),
                )
                .map_err(induction)
                .map_err(native_error)?
                .map_err(native_error)?;
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
