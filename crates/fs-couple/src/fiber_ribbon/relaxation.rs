//! Hereditary axial and bending response of the original material fibers.
//!
//! The same isotropic Young law acts on extension and curvature. For a layer
//! with a Maxwell increment Delta E, mass-normalized sine projection gives
//! c_axial^2 = Delta E A k^2 / mu and c_bending^2 = Delta E I k^4 / mu.
//! The shared owner integrates H=(c q-z)^2/2 and z'=(c q-z)/tau exactly with
//! the strings, body and contacts. No per-mode damping estimate replaces z.
//! Fiber count per strip scales A, I and mu equally, so width refinement
//! partitions the physical material without changing its constitutive poles.
//!
//! Composed with the original nonlinear fiber strain, axial arms act on
//! u_x+y_x^2/2 through the shared mean/cosine Maxwell owner. Bending arms keep
//! their original curvature coordinates. The maintained preload does not relax.
//! Source state is fixed; no temperature evolution or material slip is inferred.

use super::strain::{FiberStrain, relaxation::FiberStrainHistory};
use super::{FiberRibbonError, ModalAcousticMode, TensionedFiberRibbon};
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    ModalCoordinate,
    history::{MappedModalHistoryBinding, ModalMaterialHistory},
};
use fs_exec::CancelGate;
use fs_material::section::relaxation::ConcentricPronySection;
use fs_math::det;
use fs_phs::RelaxationBranch;

/// One original material arm acting on one mass-normalized cell mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberRelaxationChannel {
    pub layer: usize,
    pub arm: usize,
    /// Axial modes first, followed by the normal modes, as in `cell_model`.
    pub mode: usize,
    /// Square root of the projected relaxing stiffness [1/s].
    pub coupling_per_s: f64,
    /// Original material time [s], never averaged over layers or modes.
    pub relaxation_time_s: f64,
}

/// Original section spectrum and its complete projection into every fiber cell.
#[derive(Clone, Debug)]
pub struct FiberRibbonRelaxation {
    source: ConcentricPronySection,
    channels: Vec<FiberRelaxationChannel>,
    cells: usize,
    modes_per_cell: usize,
    axial_modes: usize,
    max_memory_states: usize,
    omega_band_rad_s: (f64, f64),
    finite_strain: bool,
    strain_history: Option<FiberStrainHistory>,
}

impl TensionedFiberRibbon {
    /// Bind one source Young spectrum to BOTH extension and curvature.
    ///
    /// Every material layer must be the original resolved section, including
    /// geometry, query points and receipts. Its equilibrium modulus is the
    /// existing elastic modulus. Existing explicitly sourced Kelvin--Voigt
    /// viscosity is an additional parallel dashpot; the caller must select
    /// independent loss mechanisms rather than counting one measurement twice.
    /// No material or arm is selected by a desired audible result.
    ///
    /// The original exact coupled owner retains all new histories from cold
    /// construction through normal/tangential contact and support motion.
    /// When nonlinear fiber strain is selected, every axial arm retains its
    /// geometric mean and original axial harmonics in that same transaction.
    /// The two builder calls may be made in either order. Missing data, a
    /// second binding, insufficient
    /// memory, incomplete source bandwidth or unresolved material times refuse.
    pub fn with_material_relaxation(
        mut self,
        source: ConcentricPronySection,
        max_memory_states: usize,
        gate: &CancelGate,
    ) -> Result<Self, FiberRibbonError> {
        poll(gate)?;
        if self.relaxation.is_some()
            || !(1..=256).contains(&max_memory_states)
            || source.elastic() != self.axial.elastic()
        {
            return Err(FiberRibbonError::Invalid(
                "fiber memory needs one cold binding of the identical material section and 1..=256 states",
            ));
        }
        let positive_arms = source
            .spectra()
            .iter()
            .flat_map(|s| &s.terms)
            .filter(|(delta, _)| *delta > 0.0)
            .count();
        let count = positive_arms
            .checked_mul(self.cell.modes().len())
            .and_then(|n| n.checked_mul(self.resolution.cells));
        if count.is_none_or(|n| n > max_memory_states) {
            return Err(FiberRibbonError::Invalid(
                "fiber material history exceeds its explicit state budget",
            ));
        }
        let band = [
            source.omega_band_rad_s(),
            self.axial.omega_band_rad_s(),
            self.bending.omega_band_rad_s(),
        ]
        .into_iter()
        .fold((0.0_f64, f64::MAX), |b, s| (b.0.max(s.0), b.1.min(s.1)));
        let mu = source.elastic().line_density_kg_m();
        let mut channels = Vec::with_capacity(count.unwrap() / self.resolution.cells);
        let mut instantaneous = self
            .cell
            .modes()
            .iter()
            .map(|m| m.angular_frequency_rad_s.powi(2))
            .collect::<Vec<_>>();
        for (layer, (spectrum, geometry)) in source
            .spectra()
            .iter()
            .zip(source.elastic().layer_integrals())
            .enumerate()
        {
            for (arm, &(delta, tau)) in spectrum.terms.iter().enumerate() {
                poll(gate)?;
                if delta == 0.0 {
                    continue;
                }
                for (mode, high) in instantaneous.iter_mut().enumerate() {
                    let axial = mode < self.resolution.axial_modes;
                    let n = if axial {
                        mode + 1
                    } else {
                        mode + 1 - self.resolution.axial_modes
                    };
                    let k = n as f64 * core::f64::consts::PI / self.geometry.length_m;
                    let stiffness = if axial {
                        delta * geometry.area_m2() / mu * k.powi(2)
                    } else {
                        delta * geometry.second_moment_m4() / mu * k.powi(4)
                    };
                    *high += stiffness;
                    if !stiffness.is_finite() || stiffness <= 0.0 || !high.is_finite() {
                        return Err(FiberRibbonError::Invalid(
                            "fiber material stiffness is unrepresentable",
                        ));
                    }
                    channels.push(FiberRelaxationChannel {
                        layer,
                        arm,
                        mode,
                        coupling_per_s: det::sqrt(stiffness),
                        relaxation_time_s: tau,
                    });
                }
            }
        }
        let mut image = FiberRibbonRelaxation {
            source,
            channels,
            cells: self.resolution.cells,
            modes_per_cell: self.cell.modes().len(),
            axial_modes: self.resolution.axial_modes,
            max_memory_states,
            omega_band_rad_s: band,
            finite_strain: false,
            strain_history: None,
        };
        image.validate_band(
            det::sqrt(instantaneous.into_iter().fold(0.0, f64::max)),
            self.cell.sample_period_s(),
        )?;
        if let Some(strain) = &self.strain {
            let share = self.geometry.fiber_count as f64 / self.resolution.cells as f64;
            image.with_finite_strain(
                strain,
                self.geometry.length_m,
                share,
                self.cell.modes(),
                self.cell.sample_period_s(),
            )?;
        }
        poll(gate)?;
        self.relaxation = Some(image);
        Ok(self)
    }

    /// Original source data, independent of evolving coupled material history.
    pub const fn material_relaxation(&self) -> Option<&FiberRibbonRelaxation> {
        self.relaxation.as_ref()
    }
}

fn poll(gate: &CancelGate) -> Result<(), FiberRibbonError> {
    if gate.is_requested() {
        Err(FiberRibbonError::Cancelled)
    } else {
        Ok(())
    }
}

impl FiberRibbonRelaxation {
    pub const fn source(&self) -> &ConcentricPronySection {
        &self.source
    }
    /// One cell's active linear pattern, ordered by layer, arm and mode.
    /// With full strain selected this contains bending channels only.
    pub fn channels(&self) -> &[FiberRelaxationChannel] {
        &self.channels
    }
    /// Complete history count over ALL original quadrature cells.
    pub fn memory_states(&self) -> usize {
        self.cells
            * (self.channels.len()
                + self
                    .strain_history
                    .as_ref()
                    .map_or(0, ModalMaterialHistory::memory_count))
    }
    /// Whether axial memory follows the same nonlinear strain as elastic storage.
    pub const fn has_finite_strain(&self) -> bool {
        self.finite_strain
    }
    pub const fn max_memory_states(&self) -> usize {
        self.max_memory_states
    }
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }

    pub(super) fn with_finite_strain(
        &mut self,
        strain: &FiberStrain,
        length_m: f64,
        share: f64,
        modes: &[ModalAcousticMode],
        dt: f64,
    ) -> Result<(), FiberRibbonError> {
        let history = FiberStrainHistory::from_source(strain, &self.source, length_m, share)?;
        let bending = self
            .channels
            .iter()
            .filter(|c| c.mode >= self.axial_modes)
            .count();
        let per_cell = bending.checked_add(
            history
                .as_ref()
                .map_or(0, ModalMaterialHistory::memory_count),
        );
        let total = per_cell.and_then(|n| n.checked_mul(self.cells));
        if total.is_none_or(|n| n > self.max_memory_states) {
            return Err(FiberRibbonError::Invalid(
                "fiber full-strain history needs mean, axial and bending states within its original budget",
            ));
        }
        let mut frequencies: Vec<_> = modes
            .iter()
            .map(|m| m.angular_frequency_rad_s.powi(2))
            .collect();
        for channel in self.channels.iter().filter(|c| c.mode >= self.axial_modes) {
            frequencies[channel.mode] += channel.coupling_per_s.powi(2);
        }
        let bound = history
            .as_ref()
            .map_or(0.0, |law| law.stiffness_bound(None));
        self.validate_band(
            det::sqrt(frequencies.into_iter().fold(0.0, f64::max) + bound),
            dt,
        )?;
        self.channels.retain(|c| c.mode >= self.axial_modes);
        self.strain_history = history;
        self.finite_strain = true;
        Ok(())
    }

    pub(crate) fn strain_bindings(&self, first_component: usize) -> Vec<MappedModalHistoryBinding> {
        self.strain_history.as_ref().map_or_else(Vec::new, |law| {
            (0..self.cells)
                .map(|cell| MappedModalHistoryBinding {
                    coordinates: (0..self.modes_per_cell)
                        .map(|mode| ModalCoordinate {
                            component: first_component + cell,
                            mode,
                        })
                        .collect(),
                    law: Box::new(law.clone()),
                    initial_memory: None,
                })
                .collect()
        })
    }

    pub(crate) fn validate_band(&self, high: f64, dt: f64) -> Result<(), FiberRibbonError> {
        if !high.is_finite()
            || high < 0.0
            || !dt.is_finite()
            || dt <= 0.0
            || self.omega_band_rad_s.0 > 0.0
            || high > self.omega_band_rad_s.1
            || self.channels.iter().any(|c| dt / c.relaxation_time_s > 2.0)
        {
            return Err(FiberRibbonError::Invalid(
                "complete fiber band or physical clock exceeds original material-memory support",
            ));
        }
        Ok(())
    }

    pub(crate) fn branches(
        &self,
        first_mode: usize,
        total_modes: usize,
        gate: &CancelGate,
    ) -> Result<Vec<RelaxationBranch>, FiberRibbonError> {
        let dimension = total_modes.checked_mul(2).ok_or(FiberRibbonError::Invalid(
            "fiber history dimension overflow",
        ))?;
        let end = self
            .cells
            .checked_mul(self.modes_per_cell)
            .and_then(|n| first_mode.checked_add(n));
        if end.is_none_or(|end| end > total_modes) {
            return Err(FiberRibbonError::Invalid(
                "fiber history lost its original modal coordinates",
            ));
        }
        let mut branches = Vec::with_capacity(self.memory_states());
        for cell in 0..self.cells {
            for channel in &self.channels {
                poll(gate)?;
                let mut projection = vec![0.0; dimension];
                projection[2 * (first_mode + cell * self.modes_per_cell + channel.mode)] =
                    channel.coupling_per_s;
                branches.push(RelaxationBranch {
                    projection,
                    stiffness: 1.0,
                    relaxation_time_s: channel.relaxation_time_s,
                });
            }
        }
        Ok(branches)
    }
}
