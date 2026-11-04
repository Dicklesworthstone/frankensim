//! Coherent finite-time winding heat on the original enthalpy residual.
//!
//! The native AC owner integrates every carrier and beat term once on the
//! actual thermal interval. Its spatial Joule heat enters the existing cells;
//! magnetic energy change and signed terminal work enter the same total audit.
//! Repeated thermal/Newton trials borrow that observation instead of solving
//! fields or integrating the spectrum again. Thermal conduction and caloric
//! storage still use their actual proposed endpoint.
//!
//! This is a frozen-electromagnetic-law receiver, not electrothermal feedback.
//! Mapped electrical laws must be independent of evolving temperature, and
//! magnetic source refresh must retain the original constitutive law. Fixed
//! geometry, finite source support and whole-volume correspondence remain
//! mandatory. The original infinitesimal settled AC/DC-bias interpretation is
//! unchanged: operating DC heat, switch-on history, hysteresis, finite-amplitude
//! saturation and temperature-dependent field re-solves are not supplied here.

use super::WindingThermalSample;
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate, finite,
    magnetic::{self, MagneticThermalError, MagneticThermalNetwork},
    material::{MaterialNetwork, MaterialNetworkSample},
};
use fs_em::planar::eddy::harmonic::{
    HarmonicWindingError,
    exposure::{HarmonicWindingExposure, HarmonicWindingExposureResult, HarmonicWindingHeatSource},
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::field_storage::admit_source_temperature_refresh;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Native source/geometry refusal or original harmonic-field refusal.
#[derive(Debug)]
pub enum HarmonicWindingHeatError {
    Material(MagneticThermalError),
    Field(HarmonicWindingError),
}
impl core::fmt::Display for HarmonicWindingHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent winding/enthalpy: {self:?}")
    }
}
impl std::error::Error for HarmonicWindingHeatError {}
type Error = CircuitEnthalpyError<HarmonicWindingHeatError>;
fn source(error: magnetic::Error) -> Error {
    error.map_source(HarmonicWindingHeatError::Material)
}
fn field(error: HarmonicWindingError) -> Error {
    match error {
        HarmonicWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HarmonicWindingHeatError::Field(other)),
    }
}
fn invalid(message: &'static str) -> Error {
    CircuitEnthalpyError::Invalid(message)
}

/// Accepted thermal endpoint and the exact native AC observation that heated it.
/// No field history is advanced or replaced; the spectral owner stays borrowed.
pub struct HarmonicWindingEnthalpyCandidate<'t, 's> {
    thermal: DepositionSimulation<'t>,
    report: CircuitEnthalpyReport<WindingThermalSample, HarmonicWindingExposureResult<'s>>,
}
impl<'t, 's> HarmonicWindingEnthalpyCandidate<'t, 's> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<WindingThermalSample, HarmonicWindingExposureResult<'s>> {
        &self.report
    }
    /// Publish the complete thermal endpoint. The returned AC heat has already
    /// been deposited exactly once; magnetic storage is never deposited as heat.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CircuitEnthalpyReport<WindingThermalSample, HarmonicWindingExposureResult<'s>>,
    ) {
        (self.thermal, self.report)
    }
}

// The shared private driver carries the original work account. It borrows the
// large spatial observation throughout Jacobian trials; only the final report
// takes ownership, with no per-trial copy of fields or cell exposures.
struct WindowDriver<'o, 's>(&'o HarmonicWindingExposureResult<'s>);
impl ThermalDriver for WindowDriver<'_, '_> {
    type Model = ();
    type Step = ();
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.0.dissipated_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.0.magnetic_energy_change_j,
            supplied_j: self.0.supplied_work_j,
        }
    }
    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}

impl MagneticThermalNetwork<'_> {
    /// Advance the original thermal cells using one settled, coherent AC bank.
    ///
    /// The physical integration window is exactly `[thermal.time_s(), target_s]`
    /// on the AC source's phase clock. Every resolved/stranded conducting cell
    /// requires a destination; mapped magnetic volumes cover whole receiving
    /// cells. Caloric curves, cards, selected electrical/magnetic sources and
    /// nonthermal coordinates remain pinned through every material evaluation.
    ///
    /// Electrical laws depending on the owned temperature refuse, including
    /// coincident endpoints. Changing magnetic laws retain the native refusal.
    /// Temperature-dependent thermal transport remains supported by the existing
    /// material network. Reusing this same bank on consecutive intervals retains
    /// phase and beat heating without another AC field solve.
    /// `report.circuit_trials` counts bounded thermal/material evaluations here;
    /// the original AC observation is computed only once before those trials.
    ///
    /// Thermal rates are backward Euler with interval-integrated AC heat;
    /// general thermal accuracy is first order and requires time refinement.
    /// Process start events first and split at interior topology/source events.
    /// Every refusal leaves the accepted thermal and spectral owners unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_harmonic_winding_candidate<'t, 's>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        exposure: &'s HarmonicWindingExposure<'_, '_, '_>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWindingEnthalpyCandidate<'t, 's>, Error> {
        let cancel = RefCell::new(cancelled);
        magnetic::check(&cancel).map_err(source)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if self.cells.len() != thermal.cells().len() {
            return Err(invalid(
                "complete original harmonic thermal bindings required",
            ));
        }
        let window =
            HarmonicWindow::new(thermal.time_s(), target_s - thermal.time_s()).map_err(|_| {
                invalid("harmonic heating needs a representable positive physical interval")
            })?;
        let model = exposure.model();
        let problem = model.basis().operating_field().problem();
        let conductors: BTreeMap<_, _> = model
            .eddy_conductors()
            .iter()
            .map(|s| (s.triangle, &s.conductor))
            .chain(
                model
                    .winding_conductors()
                    .iter()
                    .map(|s| (s.triangle, &s.conductor)),
            )
            .collect();
        let triangles: Vec<_> = conductors.keys().copied().collect();
        let routes = self
            .routes_for_mesh(thermal, problem.mesh(), &triangles, &cancel)
            .map_err(source)?;
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|s| (s.cell, s)).collect();
        // Same values at two endpoints cannot justify freezing a law that
        // depends on the evolving temperature between them.
        for (&triangle, original) in &conductors {
            magnetic::check(&cancel).map_err(source)?;
            let cell = routes.conductors[&triangle];
            let binding = bindings[&cell];
            if original
                .source_property()
                .source_axes()
                .iter()
                .any(|(axis, _)| axis == binding.point.temperature_axis)
            {
                return Err(invalid(
                    "temperature-dependent conductivity needs a new electromagnetic response and thermal feedback",
                ));
            }
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(magnetic::material)
            .map_err(source)?;
        let material_sample = |state: &DepositionSimulation<'_>, sample: &MaterialNetworkSample| {
            let (refreshed, magnetic) = self
                .nonlinear_materials_at(state, problem, &conductors, sample, &routes, &cancel)
                .map_err(source)?;
            for (triangle, conductor) in refreshed {
                magnetic::check(&cancel).map_err(source)?;
                let original = conductors[&triangle];
                admit_source_temperature_refresh(original.material(), conductor.material())
                    .map_err(|e| {
                        CircuitEnthalpyError::Source(HarmonicWindingHeatError::Material(
                            MagneticThermalError::Magnetic(e),
                        ))
                    })?;
                if original.source_key() != conductor.source_key()
                    || original.source_basis() != conductor.source_basis()
                    || original.resistivity_ohm_m().to_bits()
                        != conductor.resistivity_ohm_m().to_bits()
                {
                    return Err(invalid(
                        "thermal state changed the prepared AC conductivity",
                    ));
                }
            }
            Ok(magnetic)
        };
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(magnetic::material)
            .map_err(source)?;
        material_sample(thermal, &accepted.source)?;
        let observation = exposure
            .observe(window, &mut || (*cancel.borrow_mut())())
            .map_err(field)?;
        // Route physical heat only. Compensated positive sums preserve weak
        // deposits when many original elements share one finite thermal body.
        let mut sums = BTreeMap::<CellId, Dd>::new();
        for cell in &observation.cells {
            magnetic::check(&cancel).map_err(source)?;
            let triangle = match cell.source {
                HarmonicWindingHeatSource::Floating(source)
                | HarmonicWindingHeatSource::Solid { source, .. } => source.triangle,
                HarmonicWindingHeatSource::Stranded { source, .. } => source.triangle,
            };
            let id = routes
                .conductors
                .get(&triangle)
                .ok_or_else(|| invalid("AC heat lost its original material receiver"))?;
            let sum = sums.entry(*id).or_insert(Dd::ZERO);
            *sum = *sum + Dd::from_f64(cell.mean_power_w);
        }
        let powers: BTreeMap<_, _> = sums
            .into_iter()
            .map(|(id, power)| Ok((id, finite(power.to_f64())?)))
            .collect::<Result<_, Error>>()?;
        let (thermal, (), report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                magnetic::check(&cancel).map_err(source)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(magnetic::material)
                    .map_err(source)?;
                let magnetic = material_sample(state, &coefficients.source)?;
                Ok(DrivenCoefficients {
                    circuit: WindowDriver(&observation),
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: WindingThermalSample {
                        network: coefficients.source,
                        magnetic,
                        magnetic_caloric: None,
                        magnetic_balance: None,
                    },
                    circuit_trials: 1,
                })
            },
            |_driver, dt, checkpoint| {
                checkpoint()?;
                if dt.to_bits() != window.duration_s().to_bits() {
                    return Err(invalid("AC exposure and native thermal intervals disagree"));
                }
                Ok(powers.clone())
            },
            || (*cancel.borrow_mut())(),
        )?;
        magnetic::check(&cancel).map_err(source)?;
        Ok(HarmonicWindingEnthalpyCandidate {
            thermal,
            report: CircuitEnthalpyReport {
                electrical: observation,
                thermal: report.thermal,
                source: report.source,
                nonlinear: report.nonlinear,
                maximum_cell_residual_j: report.maximum_cell_residual_j,
                circuit_trials: report.circuit_trials,
                deposited_joule_j: report.deposited_joule_j,
                non_electrical_heat_j: report.non_electrical_heat_j,
                total_energy_residual_j: report.total_energy_residual_j,
            },
        })
    }
}
