//! Resolved spatial plate loads through the original material and observer.

use super::*;
use crate::curle3d::array::PreparedCompactSourceObserver;
use crate::curle3d::{CoherenceBudget, CompactFieldStatistics};
use crate::slot_jet_3d::distributed::{
    SlotJet3dPatchBudget, SlotJet3dPatchLayout, SlotJet3dPatchRun, run_slot_jet_3d_patches,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchError};

/// The flow record supplies wall dipoles only. A finite periodic-cell source
/// distribution is an authored radiation interpretation, not an isolated-jet
/// experiment or an automatic periodic-array replication.
pub const SPATIAL_PLATE_FIELD_SCOPE: &str = "Resolved stationary splitter-plate force on fluid \
from one periodic 3D flow cell; coherent compact-patch force dipoles in the same homogeneous \
stationary lossless material and free space. Authored geometric compactness/range limits \
are not an acoustic error estimate. No unresolved volume quadrupoles, stochastic SGS sound, \
finite-body scattering, source feedback, periodic-array replication or measured-SPL qualification.";

/// Authored geometric limits on each point-force patch approximation.
/// These bound phase variation and patch size relative to range, not a relative
/// field error: cancellation can amplify the error of arbitrarily small patches.
#[derive(Clone, Copy, Debug)]
pub struct CompactPatchLimits {
    /// Maximum `k * enclosing_radius`, positive and at most one radian.
    pub maximum_ka: f64,
    /// Maximum `enclosing_radius / source_observer_distance`, strictly in (0,1).
    pub maximum_radius_over_range: f64,
}

/// Same-material spatial force history, with a common physical clock.
pub struct MaterialSlotJet3dPatchRun<'a, F: FluidTransport + FluidAcoustics + ?Sized = GasState> {
    material: &'a F,
    acoustic: FluidAcousticState,
    timestep_s: f64,
    positions_m: Vec<[f64; 3]>,
    radii_m: Vec<f64>,
    force_on_fluid_n: Vec<[f64; 3]>,
    /// Maximum final-plane speed divided by this fluid's acoustic speed.
    pub final_plane_physical_mach_max: f64,
    pub lattice_diagnostics: SlotJet3dDiagnostics,
    /// Original closure and recorded lattice-unit shear diagnostics, if enabled.
    pub sgs: Option<SlotJet3dSgsDiagnostics>,
}

impl<'a, F: FluidTransport + FluidAcoustics + ?Sized> MaterialSlotJet3d<'a, F> {
    /// Run native spatial splitter loads in the same material and physical units.
    /// Optional SGS is explicit and affects the original collision trajectory;
    /// it adds no separate acoustic source. Positions use the grid's original
    /// Cartesian frame, with coordinate zero at its minimum cell boundary.
    ///
    /// # Errors
    /// Original source/closure/step, record-budget, cancellation, material speed
    /// envelope and unrepresentable force/geometry errors refuse the whole run.
    pub fn run_patches(
        &self,
        layout: SlotJet3dPatchLayout,
        budget: SlotJet3dPatchBudget,
        closure: Option<SgsConfig>,
        gate: &CancelGate,
    ) -> Result<MaterialSlotJet3dPatchRun<'_, F>, MaterialSlotJet3dSgsError> {
        let run = run_slot_jet_3d_patches(&self.config, layout, budget, closure, gate)
            .map_err(MaterialSlotJet3dSgsError::Source)?;
        self.physical_patch_run(run)
            .map_err(MaterialSlotJet3dSgsError::Material)
    }

    fn physical_patch_run(
        &self,
        mut run: SlotJet3dPatchRun,
    ) -> Result<MaterialSlotJet3dPatchRun<'_, F>, MaterialAeroacError> {
        let speed = self.scaling.velocity_si(run.diagnostics.mach_max_lattice)?;
        if speed > self.max_speed_m_s {
            return Err(invalid(
                "final-plane 3-D jet speed exceeds the declared physical envelope",
            ));
        }
        for force in &mut run.force_on_fluid_lu {
            for component in force {
                *component *= self.force_unit_n;
                if !component.is_finite() {
                    return Err(invalid("physical spatial jet force is unrepresentable"));
                }
            }
        }
        let mut positions = Vec::new();
        let mut radii = Vec::new();
        positions
            .try_reserve_exact(run.patches.len())
            .map_err(|_| invalid("physical spatial jet geometry allocation refused"))?;
        radii
            .try_reserve_exact(run.patches.len())
            .map_err(|_| invalid("physical spatial jet geometry allocation refused"))?;
        let dx = self.scaling.cell_size_m();
        for patch in run.patches {
            let position = patch.center_lu.map(|x| x * dx);
            let radius = patch.radius_lu * dx;
            if position.into_iter().any(|x| !x.is_finite()) || !radius.is_finite() || radius <= 0.0
            {
                return Err(invalid("physical spatial jet geometry is unrepresentable"));
            }
            positions.push(position);
            radii.push(radius);
        }
        Ok(MaterialSlotJet3dPatchRun {
            material: self.scaling.material(),
            acoustic: self.acoustic,
            timestep_s: self.scaling.timestep_s(),
            positions_m: positions,
            radii_m: radii,
            force_on_fluid_n: run.force_on_fluid_lu,
            final_plane_physical_mach_max: speed / self.acoustic.sound_speed_m_s(),
            lattice_diagnostics: run.diagnostics,
            sgs: run.sgs,
        })
    }
}

impl<F: FluidTransport + FluidAcoustics + ?Sized> MaterialSlotJet3dPatchRun<'_, F> {
    #[must_use]
    pub const fn material(&self) -> &F {
        self.material
    }

    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.timestep_s
    }

    /// Physical patch centers in the original flow grid's Cartesian frame.
    #[must_use]
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }

    #[must_use]
    pub fn enclosing_radii_m(&self) -> &[f64] {
        &self.radii_m
    }

    /// Sample-major force ON THE FLUID in N: `sample * patch_count + patch`.
    /// There is no division by span or patch area and no force-axis projection.
    #[must_use]
    pub fn force_on_fluid_n(&self) -> &[[f64; 3]] {
        &self.force_on_fluid_n
    }

    /// Joint Welch factors retain all patch/axis cross terms on the original
    /// physical clock. The original FFT owner reads the dynamic layout directly;
    /// no record transpose or allocation of monopole-zero channels is needed.
    ///
    /// The stored channel order is source-major `[0, Fx, Fy, Fz]`: this stationary
    /// wall measurement supplies no volume-injection term. Leakage, incomplete
    /// tails, detrending and overlapping-segment limitations retain their owners.
    pub fn estimate_welch(
        &self,
        config: WelchConfig,
        budget: WelchBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SlotJet3dPatchSpectrum<'_>, WelchError> {
        let patches = self.positions_m.len();
        let channels = patches.checked_mul(4).ok_or(WelchError::Unrepresentable)?;
        let spectrum = WelchCrossSpectrum::estimate_with(
            channels,
            self.force_on_fluid_n.len() / patches,
            |channel, sample| {
                let axis = channel % 4;
                if axis == 0 {
                    0.0
                } else {
                    self.force_on_fluid_n[sample * patches + channel / 4][axis - 1]
                }
            },
            self.timestep_s,
            config,
            budget,
            cancelled,
        )?;
        Ok(SlotJet3dPatchSpectrum {
            acoustic: self.acoustic,
            positions_m: &self.positions_m,
            radii_m: &self.radii_m,
            spectrum,
        })
    }
}

/// Measured patch coherence and the geometry/material that produced its channels.
/// This borrows geometry from its source record and owns the joint FFT factors.
pub struct SlotJet3dPatchSpectrum<'a> {
    acoustic: FluidAcousticState,
    positions_m: &'a [[f64; 3]],
    radii_m: &'a [f64],
    spectrum: WelchCrossSpectrum,
}

impl SlotJet3dPatchSpectrum<'_> {
    /// Original synchronized factors for other explicitly authored propagation
    /// owners. Channels are always `[0,Fx,Fy,Fz]` per original physical patch.
    #[must_use]
    pub const fn cross_spectrum(&self) -> &WelchCrossSpectrum {
        &self.spectrum
    }

    /// Observe one resolved interior bin with the existing three-dimensional
    /// pressure/velocity Green transfer. Every segment's patch forces propagate
    /// coherently before forming pressure variance and vector energy flux.
    /// Separation can span many wavelengths; only each individual patch must
    /// satisfy the authored limits. No coherent span sum or independent-patch
    /// assumption replaces the measured cross-patch factors.
    ///
    /// The observer is in the source grid's physical Cartesian frame. Limits
    /// are geometry admission only, not acoustic-convergence evidence or an
    /// isolated-radiator qualification. Output is bin-integrated, not PSD/Hz.
    pub fn observe_bin(
        &self,
        bin: usize,
        observer_m: [f64; 3],
        limits: CompactPatchLimits,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        let frequency = self
            .spectrum
            .frequency_hz(bin)
            .map_err(|_| refusal("spatial plate spectral bin is invalid"))?;
        if !limits.maximum_ka.is_finite()
            || limits.maximum_ka <= 0.0
            || limits.maximum_ka > 1.0
            || !limits.maximum_radius_over_range.is_finite()
            || !(0.0..1.0).contains(&limits.maximum_radius_over_range)
            || limits.maximum_radius_over_range == 0.0
            || observer_m.into_iter().any(|x| !x.is_finite())
        {
            return Err(refusal(
                "spatial plate compactness limits or observer are invalid",
            ));
        }
        let k = core::f64::consts::TAU * frequency / self.acoustic.sound_speed_m_s();
        for (position, &radius) in self.positions_m.iter().zip(self.radii_m) {
            if cancelled() {
                return Err(AeroacError::Cancelled {
                    what: "spatial plate observation",
                });
            }
            let range = det::sqrt(
                (0..3)
                    .map(|a| {
                        let displacement = observer_m[a] - position[a];
                        displacement * displacement
                    })
                    .sum(),
            );
            if !k.is_finite()
                || k <= 0.0
                || !range.is_finite()
                || range <= 0.0
                || k * radius > limits.maximum_ka
                || radius / range > limits.maximum_radius_over_range
            {
                return Err(refusal(
                    "spatial plate patch exceeds wavelength or range limit",
                ));
            }
        }
        let observer = PreparedCompactSourceObserver::new(
            self.acoustic,
            frequency,
            observer_m,
            self.positions_m,
            self.positions_m.len(),
        )?;
        let mut field = observer.evaluate_welch_bin(&self.spectrum, bin, budget, cancelled)?;
        field.scope = SPATIAL_PLATE_FIELD_SCOPE;
        Ok(field)
    }
}

fn refusal(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

#[cfg(test)]
mod tests;
