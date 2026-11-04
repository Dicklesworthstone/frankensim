//! Original gas/liquid material response in the native multimodal duct.
//!
//! One immutable interior snapshot supplies every station. Same-fluid mouths
//! reuse it; an explicitly supplied receiving fluid supplies its own radiation.
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
    ModalReflectionResponse, ModalResponse, ModalWaveResponse, PreparedModalDuct,
    mm_propagate_prepared, modal_wave_from_plane_wave, physical_response,
    reflection_from_impedance,
};
use crate::fluid::{FluidDuctError, thermoviscous_wave_from_state, wave_from_state};
use crate::modal::outlet::{self, MouthResponse, MouthSourceResponse, MouthWaveResponse};
use crate::network::PortExcitation;
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
        let response = self.evaluate_material(material, omega, boundary, false, true)?;
        Ok(physical_response(response)?.input)
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
        let response = self.evaluate_material(material, omega, boundary, true, true)?;
        Ok(physical_response(response)?)
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
        let response =
            self.evaluate_material_with_loss(material, omega, loss, boundary, false, true)?;
        Ok(physical_response(response)?.input)
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
        let response =
            self.evaluate_material_with_loss(material, omega, loss, boundary, true, true)?;
        Ok(physical_response(response)?)
    }

    /// Native input reflection from one original acoustic material, without
    /// forming an input impedance. Inlet flow nodes remain finite reflected
    /// waves. The returned response supports the original real-reference and
    /// coherent [`ModalReflectionResponse::drive_sources`] boundary solves.
    ///
    /// Native characteristic coordinates still refuse an exact endpoint
    /// cutoff. Use [`Self::drive_sources_from_material`] for physical source
    /// fields at that cutoff; its finite propagation reference stays private.
    ///
    /// # Errors
    /// As [`Self::input_impedance_from_material`], except an input-impedance
    /// pole alone is not a refusal; singular native endpoint coordinates are.
    pub fn input_reflection_from_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
    ) -> Result<ModalReflectionResponse, FluidDuctError> {
        Ok(self
            .evaluate_material(material, omega, boundary, false, false)?
            .input)
    }

    /// Native input reflection with the same single material snapshot and
    /// wall law as [`Self::input_impedance_from_material_with_loss`].
    ///
    /// # Errors
    /// As [`Self::input_reflection_from_material`], plus original transport,
    /// caloric and wall-diffusion admission failures.
    pub fn input_reflection_from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
    ) -> Result<ModalReflectionResponse, FluidDuctError> {
        Ok(self
            .evaluate_material_with_loss(material, omega, loss, boundary, false, false)?
            .input)
    }

    /// Native reflected waves and physical mouth-flow observation from one
    /// material solve. Retain this response to drive multiple coherent source
    /// vectors without repeating propagation. Input-impedance poles do not
    /// require an inverse or prevent observing the mouth motion.
    ///
    /// # Errors
    /// As [`Self::input_reflection_from_material`], including native endpoint
    /// cutoff refusals. Physical cutoff drives use [`Self::drive_sources_from_material`].
    pub fn mouth_waves_from_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
    ) -> Result<MouthWaveResponse, FluidDuctError> {
        let response = self.evaluate_material(material, omega, boundary, true, false)?;
        Ok(outlet::wave_response(response)?)
    }

    /// Native reflected waves and physical mouth motion with the original
    /// material's wall diffusion. Source work and mouth observation reuse the
    /// same wave solution as [`Self::mouth_waves_from_material`].
    ///
    /// # Errors
    /// As [`Self::input_reflection_from_material_with_loss`].
    pub fn mouth_waves_from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
    ) -> Result<MouthWaveResponse, FluidDuctError> {
        let response =
            self.evaluate_material_with_loss(material, omega, loss, boundary, true, false)?;
        Ok(outlet::wave_response(response)?)
    }

    /// Drive the original material duct with physical inlet source laws and
    /// observe every mouth mode, including at an exact transverse cutoff.
    ///
    /// Each [`PortExcitation`] prescribes volume flow or `P + Zs U = p_open`
    /// with a finite passive complex source impedance. The original coherent
    /// source solver supplies inlet pressure/flow, signed generator work,
    /// source resistance dissipation and mouth flow. It never forms the input
    /// impedance. Finite private references at cutoff do not escape as native
    /// characteristic impedances or change the physical power pairing.
    ///
    /// An unforced, exactly disconnected free mode remains at zero amplitude;
    /// an inconsistent or singular source-connected system still refuses.
    /// One original material snapshot supplies propagation and radiation.
    ///
    /// # Errors
    /// Original material/propagation refusals or source-boundary failures from
    /// [`ModalReflectionResponse::drive_sources`]. A pole of the undriven
    /// input impedance alone does not prevent a finite driven solution.
    pub fn drive_sources_from_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
        sources: &[PortExcitation],
    ) -> Result<MouthSourceResponse, FluidDuctError> {
        let response = self.evaluate_material(material, omega, boundary, true, true)?;
        Ok(outlet::wave_response(response)?.drive_sources(sources)?)
    }

    /// Physical source drive with one complete material snapshot and its
    /// admitted wall diffusion. The native source equations and source work
    /// accounting are exactly those of [`Self::drive_sources_from_material`].
    ///
    /// # Errors
    /// As [`Self::drive_sources_from_material`], plus original transport,
    /// caloric and wall-diffusion admission failures.
    pub fn drive_sources_from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
        sources: &[PortExcitation],
    ) -> Result<MouthSourceResponse, FluidDuctError> {
        let response =
            self.evaluate_material_with_loss(material, omega, loss, boundary, true, true)?;
        Ok(outlet::wave_response(response)?.drive_sources(sources)?)
    }

    /// Drive a material duct radiating into a separately declared exterior
    /// fluid through its circular mouth in an infinite rigid baffle.
    ///
    /// Interior density and sound speed alone set the duct's transverse
    /// cutoffs and axial propagation. The exterior's own rho and c set every
    /// entry of the coupled radiation impedance, including evanescent-mode
    /// coupling. Pressure and outward volume flow are continuous at the mouth;
    /// no averaged sound speed or impedance is assigned to the interface.
    /// Both original materials are evaluated once, before propagation.
    ///
    /// The stationary interface has no membrane, mean flow, mass transfer or
    /// capillary dynamics. The caller owns equilibrium of the two background
    /// states and the validity of their small-signal acoustic material laws.
    /// The receiving half-space is homogeneous and nondispersive. Radiation
    /// quadrature bounds apply to its ka, while duct cutoffs use the interior.
    ///
    /// # Errors
    /// Either original material refusal, the original flanged-radiation
    /// admission limits, and source/propagation failures. As with
    /// [`Self::drive_sources_from_material`], an undriven impedance pole does
    /// not prevent a finite physical source solution.
    pub fn drive_sources_from_material_with_exterior<
        M: FluidAcoustics + ?Sized,
        E: FluidAcoustics + ?Sized,
    >(
        &self,
        material: &M,
        exterior: &E,
        omega: f64,
        sources: &[PortExcitation],
    ) -> Result<MouthSourceResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.acoustic_state()?;
        let exterior = exterior.acoustic_state()?;
        let response = evaluate_snapshot(
            self,
            &state,
            &exterior,
            omega,
            MouthBoundary::Flanged,
            true,
            true,
            |radius| wave_from_state(&state, radius, omega),
        )?;
        Ok(outlet::wave_response(response)?.drive_sources(sources)?)
    }

    /// The same independent receiving medium with the interior material's
    /// original viscothermal wall diffusion. Only the interior must supply
    /// viscosity, conductivity and Cp/Cv; those coefficients never leak into
    /// the exterior radiation law. Loss ownership and higher-mode wall-loss
    /// approximation remain those of
    /// [`Self::drive_sources_from_material_with_loss`].
    ///
    /// # Errors
    /// As [`Self::drive_sources_from_material_with_exterior`], plus the original
    /// interior caloric, transport and wall-diffusion refusals.
    pub fn drive_sources_from_material_with_exterior_and_loss<
        M: ThermoviscousFluid + ?Sized,
        E: FluidAcoustics + ?Sized,
    >(
        &self,
        material: &M,
        exterior: &E,
        omega: f64,
        loss: LossModel,
        sources: &[PortExcitation],
    ) -> Result<MouthSourceResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.thermoviscous_state()?;
        let exterior = exterior.acoustic_state()?;
        let response = evaluate_snapshot(
            self,
            &state.acoustic(),
            &exterior,
            omega,
            MouthBoundary::Flanged,
            true,
            true,
            |radius| thermoviscous_wave_from_state(&state, radius, omega, loss),
        )?;
        Ok(outlet::wave_response(response)?.drive_sources(sources)?)
    }

    fn evaluate_material<M: FluidAcoustics + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        boundary: MouthBoundary,
        observe_mouth: bool,
        admit_boundary_cutoff: bool,
    ) -> Result<ModalWaveResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.acoustic_state()?;
        evaluate_snapshot(
            self,
            &state,
            &state,
            omega,
            boundary,
            observe_mouth,
            admit_boundary_cutoff,
            |radius| wave_from_state(&state, radius, omega),
        )
    }

    fn evaluate_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        loss: LossModel,
        boundary: MouthBoundary,
        observe_mouth: bool,
        admit_boundary_cutoff: bool,
    ) -> Result<ModalWaveResponse, FluidDuctError> {
        admit_frequency(omega)?;
        let state = material.thermoviscous_state()?;
        evaluate_snapshot(
            self,
            &state.acoustic(),
            &state.acoustic(),
            omega,
            boundary,
            observe_mouth,
            admit_boundary_cutoff,
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

#[allow(clippy::too_many_arguments)]
fn evaluate_snapshot(
    geometry: &PreparedModalDuct<'_>,
    state: &FluidAcousticState,
    exterior: &FluidAcousticState,
    omega: f64,
    boundary: MouthBoundary,
    observe_mouth: bool,
    admit_boundary_cutoff: bool,
    mut plane_wave: impl FnMut(f64) -> Result<SegmentWave, DuctError>,
) -> Result<ModalWaveResponse, FluidDuctError> {
    let n = geometry.n_modes;
    let radius = geometry.mouth_radius_m();
    let mut wave_at = |radius, admit_cutoff| {
        modal_wave_from_plane_wave(plane_wave(radius)?, radius, &geometry.gammas, admit_cutoff)
    };
    let mouth_wave = wave_at(radius, admit_boundary_cutoff)?;
    let r_load = match boundary {
        MouthBoundary::Flanged => {
            // The snapshot implements FluidAcoustics directly: the original
            // material is not re-evaluated for its exterior radiation.
            let matrix =
                crate::flanged::flanged_modal_impedance_from_material(exterior, radius, omega, n)?;
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
    Ok(mm_propagate_prepared(
        geometry,
        omega,
        state.sound_speed_m_s(),
        omega * radius / state.sound_speed_m_s(),
        mouth_wave,
        r_load,
        observe_mouth,
        admit_boundary_cutoff,
        wave_at,
    )?)
}

#[cfg(test)]
mod tests;
