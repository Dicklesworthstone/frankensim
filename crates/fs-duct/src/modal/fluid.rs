//! Original gas/liquid material response in the native multimodal duct.
//!
//! One immutable material snapshot supplies every station and the mouth.
//! Density and sound speed set each radial mode's axial wavelength and cutoff;
//! optional transport/caloric data supply the existing wall-diffusion law.
//! Radius-junction conversion, stable reflection recursion, exact-cutoff
//! pressure/flow transfer and physical outlet observations use the same owner
//! as the gas APIs. No gas state or ideal-gas relation is manufactured.
//!
//! This retains the hole-free, axisymmetric cylindrical/conical staircase and
//! its mode/slice refinement axes. Higher-mode wall losses retain the original
//! plane-ZK approximation; per-mode viscothermal eigenfunctions, asymmetric
//! modes, compliant walls and evolving phase boundaries are not supplied here.

use super::{
    ModalResponse, PreparedModalDuct, mm_propagate_prepared, modal_wave_from_plane_wave,
    physical_response, reflection_from_impedance,
};
use crate::fluid::{FluidDuctError, thermoviscous_wave_from_state, wave_from_state};
use crate::modal::outlet::MouthResponse;
use crate::{DuctError, LossModel, SegmentWave};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;

/// Boundary on the complete retained physical mouth basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MouthBoundary {
    /// Zero normal flow in every retained mode.
    Rigid,
    /// Zero pressure in every retained mode, an ideal pressure-release boundary.
    PressureRelease,
    /// Full coupled radiation matrix into the identical homogeneous fluid,
    /// through a circular mouth in an infinite rigid baffle. All retained modal
    /// cross terms participate; no condensed scalar or matched tail is used.
    Flanged,
}

impl PreparedModalDuct<'_> {
    /// Physical multimodal input impedance from one original acoustic material.
    ///
    /// The material is evaluated once. Its density and adiabatic sound speed
    /// drive both propagation and, for [`MouthBoundary::Flanged`], the original
    /// full radiation matrix. Liquid, real-gas and explicitly declared mixture
    /// sources retain their own phase and relaxation assumptions. The exterior
    /// is stationary and nondispersive; no bulk absorption is inferred.
    ///
    /// # Errors
    /// Original material failures, nonpositive/nonfinite frequency, propagation
    /// or input-impedance poles, and original flanged-radiation admission limits.
    /// Geometry and retained modes were admitted by [`Self::new`].
    pub fn input_impedance_from_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
    ) -> Result<ModalResponse, FluidDuctError> {
        self.evaluate_material(material, omega, boundary, false)
            .map(|response| response.input)
    }

    /// Complete physical mouth-flow map and input impedance from the same
    /// material solve as [`Self::input_impedance_from_material`].
    ///
    /// `U_mouth = flow_map U_inlet` retains every modal phase and cross term.
    /// Exact transverse cutoffs use finite pressure/flow transfer; genuine
    /// physical poles still refuse. No singular native-wave reference escapes.
    ///
    /// # Errors
    /// As [`Self::input_impedance_from_material`].
    pub fn mouth_response_from_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
    ) -> Result<MouthResponse, FluidDuctError> {
        self.evaluate_material(material, omega, boundary, true)
    }

    /// The same modal input impedance with material-driven viscous/thermal wall
    /// diffusion. One complete [`ThermoviscousFluid`] snapshot supplies rho, c,
    /// viscosity, conductivity, Cp and actual Cp/Cv at every station radius.
    ///
    /// [`LossModel::Bessel`] retains the native all-shear-number plane response;
    /// [`LossModel::WideTube`] retains its thin-boundary-layer admission. Higher
    /// modes use the existing plane-ZK-corrected axial dispersion and impedance,
    /// exact in the lossless limit. Neither law supplies bulk relaxation loss.
    ///
    /// # Errors
    /// As [`Self::input_impedance_from_material`], plus original transport,
    /// caloric and wall-diffusion admission failures.
    pub fn input_impedance_from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
    ) -> Result<ModalResponse, FluidDuctError> {
        self.evaluate_material_with_loss(material, omega, loss, boundary, false)
            .map(|response| response.input)
    }

    /// Complete physical mouth-flow map from the same material and wall-loss
    /// solve as [`Self::input_impedance_from_material_with_loss`].
    ///
    /// # Errors
    /// As [`Self::input_impedance_from_material_with_loss`].
    pub fn mouth_response_from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
    ) -> Result<MouthResponse, FluidDuctError> {
        self.evaluate_material_with_loss(material, omega, loss, boundary, true)
    }

    fn evaluate_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
        observe_mouth: bool,
    ) -> Result<MouthResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.acoustic_state()?;
        evaluate_snapshot(self, &state, omega, boundary, observe_mouth, |radius| {
            wave_from_state(&state, radius, omega)
        })
    }

    fn evaluate_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
        observe_mouth: bool,
    ) -> Result<MouthResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.thermoviscous_state()?;
        evaluate_snapshot(
            self,
            &state.acoustic(),
            omega,
            boundary,
            observe_mouth,
            |radius| thermoviscous_wave_from_state(&state, radius, omega, loss),
        )
    }
}

fn admit_frequency(omega: f64) -> Result<(), DuctError> {
    if !omega.is_finite() || omega <= 0.0 {
        return Err(DuctError::BadParameter {
            what: "omega must be finite and positive",
        });
    }
    Ok(())
}

fn evaluate_snapshot(
    geometry: &PreparedModalDuct<'_>,
    state: &FluidAcousticState,
    omega: f64,
    boundary: MouthBoundary,
    observe_mouth: bool,
    mut plane_wave: impl FnMut(f64) -> Result<SegmentWave, DuctError>,
) -> Result<MouthResponse, FluidDuctError> {
    let n = geometry.n_modes;
    let radius = geometry.mouth_radius_m();
    let mut wave_at = |radius, admit_cutoff| {
        modal_wave_from_plane_wave(plane_wave(radius)?, radius, &geometry.gammas, admit_cutoff)
    };
    let mouth_wave = wave_at(radius, true)?;
    let r_load = match boundary {
        MouthBoundary::Flanged => {
            // The snapshot implements FluidAcoustics directly: the original
            // material is not re-evaluated for its exterior radiation.
            let matrix =
                crate::flanged::flanged_modal_impedance_from_material(state, radius, omega, n)?;
            reflection_from_impedance(&matrix, &mouth_wave.zc, n)?
        }
        MouthBoundary::Rigid | MouthBoundary::PressureRelease => {
            // Preserve the existing gas boundary arithmetic, including with
            // complex thermoviscous references; apply it to the whole basis.
            let reflection = if boundary == MouthBoundary::Rigid {
                C64::ONE
            } else {
                let zc0 = mouth_wave.zc[0];
                if zc0.abs() < 1e-300 {
                    return Err(DuctError::Singular.into());
                }
                (C64::ZERO - zc0) * (C64::ZERO + zc0).recip()
            };
            let mut matrix = vec![C64::ZERO; n * n];
            for mode in 0..n {
                matrix[mode * n + mode] = reflection;
            }
            matrix
        }
    };
    let response = mm_propagate_prepared(
        geometry,
        omega,
        state.sound_speed_m_s(),
        omega * radius / state.sound_speed_m_s(),
        mouth_wave,
        r_load,
        observe_mouth,
        true,
        wave_at,
    )?;
    Ok(physical_response(response)?)
}

#[cfg(test)]
mod tests;
