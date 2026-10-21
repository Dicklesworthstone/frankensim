//! Source-bound direct axial termination work in the existing support image.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
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
        let ends:Vec<_>=attachments.iter().map(|a|a.moving_end).collect();
        Self::build(sources,support,Some(&ends),SupportGasFrame::StationaryModal {remainder_frequency_hz},
            Some(&attachments),gate)
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
