//! Close-coiled material windings lowered to the existing string mechanics.
use fs_material::section::winding::{CloseWoundSection, WoundCoreSection};
use fs_scenario::PrestressedString;
use crate::acoustic_realize::AcousticRealizeError;
use crate::string_specimen::StringPrestress;
use super::concentric::with_section_coefficients;

/// An immutable wound specimen with independent core and wrap source receipts.
#[derive(Clone, Debug)]
pub struct ResolvedWoundStringSpecimen {
    string: PrestressedString,
    section: CloseWoundSection,
    prestress: StringPrestress,
}
impl ResolvedWoundStringSpecimen {
    /// Bind mass, spring extension, core flexure/loss and finished diameter in
    /// one transaction. No coefficients or damping survive from another string.
    /// This is the explicit sliding, core-dominated reduction documented on
    /// `CloseWoundSection`, not a resolved winding/contact calculation.
    pub fn try_new(template: PrestressedString, section: CloseWoundSection, prestress: StringPrestress)
        -> Result<Self, AcousticRealizeError>
    {
        section.admit_wavenumber(std::f64::consts::PI * template.n_modes as f64 / template.length_m)
            .map_err(|_| AcousticRealizeError::InvalidDescription {
                what: "wound string requires wavelengths long compared with winding pitch",
            })?;
        let mut string = with_section_coefficients(template, section.line_density_kg_m(),
            section.axial_stiffness_n(), section.core().elastic().bending_stiffness_n_m2(),
            section.outer_radius_m(), prestress)?;
        match section.core() {
            WoundCoreSection::KelvinVoigt(core) => {
                string.kelvin_voigt_bending = Some(fs_scenario::acoustic::KelvinVoigtBending {
                    viscous_stiffness_n_m2_s: core.viscous_stiffness_n_m2_s(),
                    omega_band_rad_s: section.omega_band_rad_s(),
                    material_state_identity: None,
                });
            }
            WoundCoreSection::Prony(core) => {
                string.relaxation_bending = Some(fs_scenario::acoustic::PronyBending {
                    branches: core.branches().iter().map(|arm| fs_scenario::acoustic::BendingRelaxationBranch {
                        relaxing_stiffness_n_m2: arm.relaxing_stiffness_n_m2,
                        relaxation_time_s: arm.relaxation_time_s,
                    }).collect(),
                    // The wrap is elastic, but its coefficients still have to
                    // support the instantaneous frequencies of the whole owner.
                    omega_band_rad_s: section.omega_band_rad_s(),
                    // Both original source bundles and per-layer selectors
                    // remain in section; there is no equivalent material ID.
                    material_state_identity: None,
                    source_properties: None,
                });
            }
        }
        Ok(Self { string, section, prestress })
    }
    /// Ordinary generic string coefficients; source data stay on this specimen.
    #[must_use] pub fn string(&self) -> PrestressedString { self.string.clone() }
    /// Original core, winding geometry and independent material identities.
    #[must_use] pub const fn section(&self) -> &CloseWoundSection { &self.section }
    /// Applied force, extension or explicitly requested setup tuning.
    #[must_use] pub const fn prestress(&self) -> StringPrestress { self.prestress }
}
