//! Source-bound direct axial termination work in the existing support image.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::{
    LongitudinalStringResolution,LongitudinalMaterialLoss,
};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::plate::impact::supported::longitudinal::{LongitudinalSupportModes, StringSupportGeometry};
use fs_math::det;

/// A material string's actual spatial attachment and explicit axial resolution.
/// Axial mass and EA are read from its original material, not entered again.
#[derive(Clone, Debug)]
pub struct MaterialLongitudinalSupport {
    pub moving_end: StringEndMotion,
    pub geometry: StringSupportGeometry,
    pub resolution: LongitudinalStringResolution,
    /// Original mass-normalized solid axial sine states in selected-mode order.
    pub initial_axial: Vec<ModalAcousticState>,
}
impl MaterialLongitudinalSupport {
    pub(super) fn lower(&self,source:&NonlinearMaterialString)
        -> Result<(LongitudinalSupportModes,Option<(f64,f64)>),ImpactError> {
        self.resolution.validate()?;
        if self.initial_axial.len()!=self.resolution.axial_modes.len()
            || self.initial_axial.iter().any(|s|!s.displacement_m_sqrt_kg.is_finite()
                || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
            return Err(invalid("supported axial selection needs every finite initial solid mode"));
        }
        let MaterialStringSource::WithPolarizations(p)=source.source() else {
            return Err(invalid("spatial axial support requires both original transverse polarizations"));
        };
        let direction=self.geometry.components();let norm=det::hypot(direction[0],direction[1]);
        if (self.geometry.length_m()-source.card().length_m).abs()>1e-12*source.card().length_m
            || (0..2).any(|i|(direction[i]-norm*p.directions().bridge()[i]).abs()>1e-12) {
            return Err(invalid("spatial attachment must retain the source span and declared transverse bridge direction"));
        }
        let mut band=None;
        let solid=source.source().longitudinal_material(&mut band)?;
        Ok((LongitudinalSupportModes {solid_linear_density_kg_m:solid,
            modes:self.resolution.axial_modes.clone(),maximum_strain:self.resolution.maximum_strain},band))
    }
}
impl SupportedMaterialStrings {
    /// Retain material axial waves and direct endpoint strain on an oblique
    /// translating support. Actual span geometry determines the projection;
    /// changing density or EA changes the inertia, waves and support reaction.
    /// The same owner retains both transverse planes, initial torsion and
    /// stationary-gas resistance. Axial waves are elastic and conservative:
    /// no axial viscosity or axial gas law is inferred from bending loss.
    /// The reference is already tensioned/equilibrated, so its static traction
    /// is not reapplied. This is one support translation, not a rotating bridge.
    pub fn in_stationary_gas_with_longitudinal_waves(sources:Vec<NonlinearMaterialString>,
        support:TranslatingSupport, attachments:Vec<MaterialLongitudinalSupport>,
        remainder_frequency_hz:f64, gate:&CancelGate) -> Result<Self,ImpactError> {
        Self::in_stationary_gas_with_longitudinal_waves_and_loss(sources,support,attachments,
            LongitudinalMaterialLoss::Elastic,remainder_frequency_hz,gate)
    }

    /// The same complete source-bound support image with explicit axial loss.
    /// Every homogeneous/bonded solid layer supplies its own uniaxial viscosity.
    /// The original fixed-end bank and this moving termination share the same
    /// material binding and solid axial inertia.
    ///
    /// With `KelvinVoigt`, relative axial motion and the moving endpoint's uniform
    /// extension dissipate through the original reciprocal resistance matrix. Nonzero
    /// transverse/axial/torsion states, stationary-gas channels and the support
    /// input survive. There is no second damping step or copied audio engine.
    ///
    /// A moving support mixes axial and lower-frequency mechanical motion. Each
    /// source band must therefore cover DC through the existing conservative
    /// bound on the coupled support/string SMALL-SIGNAL spectrum, not just axial
    /// poles. This does not certify every finite-amplitude generated harmonic.
    /// `FiniteStrainKelvinVoigt` instead applies the complete retained strain
    /// rate, including both transverse planes and affine endpoint motion, in the
    /// SAME nonlinear solve. Evolving heat and axial fluid loss remain absent.
    /// Missing, frequency-dependent or out-of-band source data refuses.
    pub fn in_stationary_gas_with_longitudinal_waves_and_loss(sources:Vec<NonlinearMaterialString>,
        support:TranslatingSupport, attachments:Vec<MaterialLongitudinalSupport>,
        loss:LongitudinalMaterialLoss, remainder_frequency_hz:f64, gate:&CancelGate)
        -> Result<Self,ImpactError>
    {
        poll(gate)?;
        let finite_strain=matches!(&loss,LongitudinalMaterialLoss::FiniteStrainKelvinVoigt {..});
        if sources.is_empty() || sources.len()>32 || attachments.len()!=sources.len() {
            return Err(invalid("material axial support requires 1..=32 complete source attachments"));
        }
        let ends:Vec<_>=attachments.iter().map(|a|a.moving_end).collect();
        let mut image=Self::build(sources,support,Some(&ends),SupportGasFrame::StationaryModal {remainder_frequency_hz},
            Some(&attachments),gate)?;
        if matches!(loss,LongitudinalMaterialLoss::Elastic) {return Ok(image);}
        let upper=image.supported.clock_bounds().iter().copied().fold(0.0_f64,f64::max);
        let mut laws=Vec::with_capacity(image.sources.len());
        let mut rigidities=Vec::with_capacity(image.sources.len());
        for (index,source) in image.sources.iter().enumerate() {
            poll(gate)?;
            let law=loss.bind(source.source())?
                .ok_or_else(||invalid("selected axial material loss is missing its constitutive source"))?;
            for &omega in image.supported.axial_angular_frequencies(index)? {
                law.longitudinal_damping_ratio(omega).map_err(|e|ImpactError::Owner(e.to_string()))?;
            }
            let (lo,hi)=law.omega_band_rad_s();
            if lo>0.0 || hi<upper {
                return Err(invalid("moving-support axial viscosity must cover DC through the complete small-signal spectral bound"));
            }
            rigidities.push(law.viscous_rigidity_n_s());laws.push(law);
        }
        image.supported=if finite_strain {
            image.supported.with_finite_strain_axial_viscosity(rigidities)?
        } else {image.supported.with_axial_viscosity(rigidities)?};
        image.bodies[0].potential=BodyPotential::SupportedStrings(image.supported.clone());
        image.axial_material_loss=Some(laws);
        poll(gate)?;
        Ok(image)
    }
    /// Original material loss receipts, in source order. Elastic construction
    /// retains None; an explicit zero-viscosity source remains a retained law.
    pub fn axial_material_loss(&self)->Option<&[fs_material::section::axial::ConcentricAxialViscosity]> {
        self.axial_material_loss.as_deref()
    }
    /// Original geometric attachments and all explicit axial resolutions/states.
    pub fn longitudinal_supports(&self)->Option<&[MaterialLongitudinalSupport]> {
        self.longitudinal.as_deref()
    }
    /// Work-conjugate absolute axial point motion and force in the full image.
    pub fn axial_shape(&self,string:usize,station:f64)->Result<Vec<f64>,ImpactError> {
        let mut row=self.supported.absolute_axial_shape(string,station)?;
        row.resize(self.mode_count,0.0);Ok(row)
    }
}
