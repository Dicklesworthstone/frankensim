//! Photon-conserving transport through ordered, individually mixed gas cells.
//!
//! The beam has one photon energy and one area throughout the column. Each
//! cell may have its own inventory, volume, temperatures and atomic source.
//! Only photoionization changes its opacity during the exposure.
//!
//! For any incident flux history F(t), a cell obeys
//! `d tau / d Phi = -sigma * (1-exp(-tau))`, with `Phi = integral F(t) dt`.
//! Its endpoint therefore depends only on incident fluence, and its outgoing
//! fluence is the incident fluence minus the absorbed photon column. Passing
//! this integral to the next cell gives the exact serial endpoint under the
//! same model, including changing upstream opacity. A downstream cell's
//! reported exposure uses the equivalent mean flux; it is not the actual
//! time-varying flux history. There is one analytic update per cell, without
//! a reaction time-step or a frozen-opacity approximation.
//!
//! These statements require fixed geometry, fixed source cross sections,
//! quasi-static propagation and the parent module's immediate thermalization
//! assumptions. Recombination, collisional reactions, scattering, other
//! absorbers, heat exchange, motion and light-travel delays are absent during
//! the interval. Spatial structure within a cell is not resolved.

use super::{
    PhotoionizationBalance, PhotoionizationControls, PhotoionizationError, PhotoionizationExposure,
    PhotoionizationResult, advance_monochromatic_photoionization, nonnegative, positive,
    validate_exposure,
};
use fs_math::dd::Dd;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::AtomicPhotoionizationV1;

/// An immutable cell and its original atomic source, ordered toward the exit.
#[derive(Clone, Copy, Debug)]
pub struct PhotoionizationCell<'a> {
    /// Initial gas inventory and thermodynamic state.
    pub initial: &'a FrozenIonizedGasStateV1,
    /// Bound-free cross section and source chemical gap for this cell.
    pub source: &'a AtomicPhotoionizationV1,
}

/// Explicit column work and numerical conservation allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoionizationColumnControls {
    /// Cell admission controls. The energy allowance also applies once to the
    /// whole column; it is not multiplied by the number of cells.
    pub cell: PhotoionizationControls,
    /// Maximum cells accepted before any update or result allocation.
    pub maximum_cells: usize,
    /// Absolute photon-account residual divided by incident photons, in (0,1).
    /// This is a numerical allowance, not a physical model error estimate.
    pub relative_photon_tolerance: f64,
}

/// Atomic column refusal; no caller-owned cell is ever changed.
#[derive(Clone, Debug, PartialEq)]
pub enum PhotoionizationColumnError {
    /// Invalid column controls or common exposure.
    Invalid(PhotoionizationError),
    /// The caller's cell-count budget was exceeded.
    CellLimit { requested: usize, maximum: usize },
    /// Original source or thermodynamic update refused at this zero-based cell.
    Cell {
        index: usize,
        error: PhotoionizationError,
    },
    /// Column accounts, mean flux, or result storage cannot be represented.
    Arithmetic(&'static str),
    /// Summed absorptions plus outgoing photons do not close against incoming.
    PhotonBalance { residual_photons: f64 },
    /// The whole column exceeds the unchanged absolute joule allowance.
    EnergyBalance { residual_j: f64 },
    /// Cancellation before publishing a complete column.
    Cancelled,
}
impl core::fmt::Display for PhotoionizationColumnError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "monochromatic photoionization column: {self:?}")
    }
}
impl std::error::Error for PhotoionizationColumnError {}

/// Complete source-retained cell endpoints and the column's physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotoionizationColumnResult {
    cells: Vec<PhotoionizationResult>,
    exposure: PhotoionizationExposure,
    controls: PhotoionizationColumnControls,
    balance: PhotoionizationBalance,
    photon_residual: f64,
}
impl PhotoionizationColumnResult {
    /// Endpoints in incident-to-exit order; each exposure is fluence-equivalent.
    pub fn cells(&self) -> &[PhotoionizationResult] {
        &self.cells
    }
    /// Original constant illumination entering the first cell.
    pub const fn exposure(&self) -> PhotoionizationExposure {
        self.exposure
    }
    /// Applied work and numerical conservation controls.
    pub const fn controls(&self) -> PhotoionizationColumnControls {
        self.controls
    }
    /// Incoming/outgoing photons and summed material changes over all cells.
    pub const fn balance(&self) -> PhotoionizationBalance {
        self.balance
    }
    /// Compensated sum of all absorptions plus outgoing minus incoming photons.
    pub const fn photon_balance_residual(&self) -> f64 {
        self.photon_residual
    }
}

/// Advance a fixed-area monochromatic beam through stratified photoionizing gas.
///
/// Empty columns transmit every photon. Source bindings and temperature domains
/// are checked even in dark cells. Cancellation is polled within each analytic
/// cell update and before final publication. Any late source, budget, arithmetic
/// or conservation refusal discards all local endpoints without modifying the
/// inputs. Time-varying downstream fluxes are handled through their exact
/// integrated fluence, not by freezing an upstream optical depth.
pub fn advance_monochromatic_column(
    cells: &[PhotoionizationCell<'_>],
    exposure: PhotoionizationExposure,
    controls: PhotoionizationColumnControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhotoionizationColumnResult, PhotoionizationColumnError> {
    if cancelled() {
        return Err(PhotoionizationColumnError::Cancelled);
    }
    if cells.len() > controls.maximum_cells {
        return Err(PhotoionizationColumnError::CellLimit {
            requested: cells.len(),
            maximum: controls.maximum_cells,
        });
    }
    if !positive(controls.relative_photon_tolerance) || controls.relative_photon_tolerance >= 1.0 {
        return Err(PhotoionizationColumnError::Invalid(
            PhotoionizationError::Invalid("relative column photon tolerance in (0,1)"),
        ));
    }
    let (_, incident) =
        validate_exposure(exposure, controls.cell).map_err(PhotoionizationColumnError::Invalid)?;
    let mut endpoints = Vec::new();
    endpoints
        .try_reserve_exact(cells.len())
        .map_err(|_| PhotoionizationColumnError::Arithmetic("column result allocation"))?;
    let mut local_exposure = exposure;
    let mut transmitted = incident;
    // Fixed incident-to-exit reduction order. Keep small cell contributions
    // during aggregation instead of absorbing them into a large leading sum.
    let mut sums = [Dd::ZERO; 4];
    for (index, cell) in cells.iter().enumerate() {
        let result = advance_monochromatic_photoionization(
            cell.initial,
            cell.source,
            local_exposure,
            controls.cell,
            &mut cancelled,
        )
        .map_err(|error| match error {
            PhotoionizationError::Cancelled => PhotoionizationColumnError::Cancelled,
            error => PhotoionizationColumnError::Cell { index, error },
        })?;
        let b = result.balance();
        for (sum, value) in sums.iter_mut().zip([
            b.absorbed_photons,
            b.chemical_storage_change_j,
            b.electron_sensible_energy_change_j,
            b.absorbed_energy_j,
        ]) {
            *sum = *sum + Dd::from_f64(value);
        }
        transmitted = b.transmitted_photons;
        if index + 1 < cells.len() && transmitted != b.incident_photons {
            let mean_flux = if exposure.duration.value() == 0.0 {
                0.0
            } else {
                (transmitted / exposure.illuminated_area.value()) / exposure.duration.value()
            };
            if !nonnegative(mean_flux) || (transmitted > 0.0 && mean_flux == 0.0) {
                return Err(PhotoionizationColumnError::Arithmetic(
                    "downstream mean photon flux",
                ));
            }
            local_exposure.beam.number_flux_per_m2_s = mean_flux;
        }
        endpoints.push(result);
    }
    let photon_residual = (sums[0] + Dd::from_f64(transmitted) - Dd::from_f64(incident)).to_f64();
    let energy_residual = (sums[1] + sums[2] - sums[3]).to_f64();
    if sums.iter().any(|sum| !sum.to_f64().is_finite())
        || !photon_residual.is_finite()
        || !energy_residual.is_finite()
    {
        return Err(PhotoionizationColumnError::Arithmetic(
            "column photon or energy accounts",
        ));
    }
    if photon_residual.abs() > controls.relative_photon_tolerance * incident {
        return Err(PhotoionizationColumnError::PhotonBalance {
            residual_photons: photon_residual,
        });
    }
    if energy_residual.abs() > controls.cell.energy_tolerance_j {
        return Err(PhotoionizationColumnError::EnergyBalance {
            residual_j: energy_residual,
        });
    }
    if cancelled() {
        return Err(PhotoionizationColumnError::Cancelled);
    }
    Ok(PhotoionizationColumnResult {
        cells: endpoints,
        exposure,
        controls,
        photon_residual,
        balance: PhotoionizationBalance {
            incident_photons: incident,
            absorbed_photons: sums[0].to_f64(),
            transmitted_photons: transmitted,
            chemical_storage_change_j: sums[1].to_f64(),
            electron_sensible_energy_change_j: sums[2].to_f64(),
            absorbed_energy_j: sums[3].to_f64(),
            energy_residual_j: energy_residual,
        },
    })
}
