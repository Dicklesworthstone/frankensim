//! Existing single, paired and modal mechanical owners on one exact DG load.
//!
//! LoadedAperture owns dispatch, contact, spatial traction and material memory.
//! This adapter only supplies the same-step DG mean-pressure relation and uses
//! the original shared acoustic/work candidate. It never projects independent
//! blades or modes onto the scalar diagnostic opening coordinate.

/// Finite authored pressure, force and spatial-pressure playback.
pub mod performance;

use super::{ApertureDg, ApertureDgLoad, DgApertureFrame, accept_acoustic_candidate, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::paired::BladePropagator;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;

/// Retains every original blade/mode, closure and material history while using
/// the same DG owner as the single-coordinate aperture and flute source.
pub type LoadedApertureDg = ApertureDg<LoadedAperture>;

impl ApertureDg<LoadedAperture> {
    /// Bind cold original mechanics to its actual DG mean impedance and gas.
    /// Every family must use ExactZoh. Material-bound paired/modal gas states
    /// must match the inlet completely, including transport properties.
    pub fn new_loaded(
        aperture: impl Into<LoadedAperture>,
        load: ApertureDgLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let aperture = aperture.into();
        let (density, impedance) = aperture.load();
        let gas = load.duct().inlet_gas();
        if load.has_driven_mouth()
            || aperture.accepted_steps() != 0
            || load.state().step_count() != 0
            || aperture.blade_propagator() != BladePropagator::ExactZoh
            || aperture.time_step_s().to_bits() != load.propagator().dt().to_bits()
            || density.to_bits() != gas.density.to_bits()
            || impedance.to_bits() != load.impedance_pa_s_m3().to_bits()
            || matches!(&aperture, LoadedAperture::Paired(a) if a.gas().is_some_and(|g| g != gas))
            || matches!(&aperture, LoadedAperture::Modal(a) if a.gas().is_some_and(|g| g != gas))
        {
            return Err(invalid(
                "loaded DG coupling requires cold ExactZoh mechanics bound to its inlet gas, mean inlet-pressure impedance and clock; driven exterior mouth effort is unsupported",
            ));
        }
        Ok(Self { aperture, load, upstream: None })
    }

    /// Mechanical, contact, material and DG storage [J]. The original rational
    /// terminal's signed exchange stays separate from its unmodeled storage.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        let energy = self.aperture.energy() + self.load.stored_energy_j();
        self.upstream.as_ref().map_or(energy, |tract| energy + tract.stored_energy_j())
    }

    /// Continue all original mechanical, acoustic and material histories.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        self.aperture.extend_step_budget(total)
    }

    /// Advance the physical upstream pressure and optional supplied body flow.
    pub fn step(&mut self, drive: TubeDrive) -> Result<DgApertureFrame, AcousticRealizeError> {
        self.step_with_forces_and_pressures(drive, [0.0; 2], [0.0; 2])
    }

    /// Independent blade forces or the original modal force-projection columns.
    /// A Single source refuses a nonzero second force instead of dropping it.
    pub fn step_with_forces(
        &mut self,
        drive: TubeDrive,
        forces_n: [f64; 2],
    ) -> Result<DgApertureFrame, AcousticRealizeError> {
        self.step_with_forces_and_pressures(drive, forces_n, [0.0; 2])
    }

    /// Solve all original mechanical coordinates against the same-step DG load.
    /// Prescribed modal pressure tractions retain their source-bound projections
    /// and are credited once in the mechanical work. Non-modal sources refuse
    /// nonzero spatial pressure. Any fallible trial leaves every owner unchanged.
    pub fn step_with_forces_and_pressures(
        &mut self,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
    ) -> Result<DgApertureFrame, AcousticRealizeError> {
        if self.upstream.is_some() {
            return self.step_with_upstream(drive, forces_n, pressures_pa);
        }
        let relation = self.load.port_response()?;
        let trial = self.aperture.preview(ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: relation.history_pressure_pa / 2.0,
            body_flow_m3_s: drive.body_flow_m3_s,
        }, forces_n, pressures_pa)?;
        let frame = accept_acoustic_candidate(&mut self.load, relation, trial.frame(),
            trial.spatial_pressure(), drive)?;
        self.aperture.accept(trial);
        Ok(frame)
    }
}
