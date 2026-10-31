//! Source-bound bowed strings through one modal-loss lowering.
//!
//! Uniform specimens and perfectly bonded concentric sections reuse the same
//! section integrals, support prescription, beam dispersion, Kelvin-Voigt
//! dilution and oscillating-cylinder gas resistance as the general assembly.
//! There is no instrument-name switch, equivalent-density material, copied
//! decay envelope, or second bow stepper. Rebuild the card after a data change;
//! this is frozen-state compilation, not a live material-history transfer.
//!
//! Concentric annuli are NOT helical windings: winding angle, interlayer slip,
//! torsion and winding contact loss are not supplied by this reduction. The
//! gas model supplies small-motion drag, not fluid added mass or turbulent flow.
//! Body modes, friction, and observation geometry remain independent inputs.

use super::BowedStringCard;
use crate::acoustic_realize::{AcousticRealizeError, mode_zeta, validate_incremental_string};
use crate::material_assembly::concentric::with_concentric_kelvin_voigt_section;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::string_specimen::{ResolvedStringSpecimen, StringPrestress};
use fs_exec::Cx;
use fs_material::gas::GasState;
use fs_material::section::ConcentricKelvinVoigtSection;
use fs_material::section::modal_loss::ConcentricModalLossSection;
use fs_scenario::PrestressedString;

fn poll(cx: &Cx<'_>) -> Result<(), AcousticRealizeError> {
    cx.checkpoint().map_err(|_| AcousticRealizeError::Cancelled)
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl BowedStringCard {
    /// Lower direct whole-string measurements through the original transverse
    /// dispersion, cylinder resistance and time-domain admission. Intrinsic
    /// material loss is absent, explicitly unmeasured, unless the specimen
    /// selected its source's own loss factors; no Q is converted into a
    /// viscosity, torsional law or invented homogeneous solid.
    pub fn from_measured_string_specimen(
        cx: &Cx<'_>,
        specimen: &crate::material_assembly::measured_string::ResolvedMeasuredStringSpecimen,
        gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = specimen.string();
        validate_incremental_string(&string)?;
        for mode in 1..=string.n_modes {
            poll(cx)?;
            specimen.admit_angular_frequency(fs_nlmodal::prestressed_beam_omega(
                string.length_m, string.tension_n, string.lin_density_kg_m,
                string.bending_stiffness_n_m2, mode,
            ))?;
        }
        lower_backbone(cx, &string, gas, sample_rate_hz, 0.0)
    }

    /// Compile the admitted close-coiled material winding through the same
    /// string dispersion, gas drag, bridge shear and source-band checks.
    pub fn from_wound_material_specimen(
        cx: &Cx<'_>, specimen: &crate::material_assembly::wound::ResolvedWoundStringSpecimen,
        gas: &GasState, sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        lower(cx, &specimen.string(), gas, sample_rate_hz)
    }

    // Only the joint bank may lower a core with memory; the public card path
    // keeps refusing rather than silently dropping the internal strains.
    pub(crate) fn wound_relaxation_backbone_with_gas(
        cx: &Cx<'_>, specimen: &crate::material_assembly::wound::ResolvedWoundStringSpecimen,
        gas: &GasState, sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = specimen.string();
        if !matches!(specimen.section().core(), fs_material::section::winding::WoundCoreSection::Prony(_))
            || string.relaxation_bending.is_none() {
            return Err(invalid("wound bank backbone requires the original core bending spectrum"));
        }
        let mut elastic = string.clone();
        elastic.relaxation_bending = None;
        validate_incremental_string(&elastic)?;
        lower_backbone(cx, &string, gas, sample_rate_hz, 0.0)
    }

    // A multilayer Prony source takes the same elastic/gas backbone as a
    // homogeneous source. Only the bank may use it: that owner immediately
    // installs every original layer arm in its exact material-memory state.
    pub(crate) fn concentric_relaxation_backbone_with_gas(
        cx: &Cx<'_>, template: PrestressedString,
        section: &fs_material::section::relaxation::ConcentricPronySection,
        prestress: StringPrestress, gas: &GasState, sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = crate::material_assembly::concentric::with_concentric_prony_section(
            template, section, prestress)?;
        let mut elastic = string.clone();
        elastic.relaxation_bending = None;
        validate_incremental_string(&elastic)?;
        lower_backbone(cx, &string, gas, sample_rate_hz, 0.0)
    }

    // Private to a bank constructor which immediately binds the original
    // source's material histories. Memoryless cards retain their refusal.
    pub(crate) fn material_relaxation_backbone_with_gas(
        cx: &Cx<'_>, specimen: &ResolvedStringSpecimen, gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = specimen.string();
        let law = string.relaxation_bending.as_ref()
            .ok_or_else(|| invalid("material bank backbone requires original bending memory"))?;
        if law.material_state_identity != Some(specimen.material().identity())
            || law.source_properties.as_ref().is_none_or(|p| p.len() != law.branches.len())
            || law.branches.len() > 8 || string.kelvin_voigt_bending.is_some()
            || string.rayleigh.is_some() || string.damping_ratio != 0.0 {
            return Err(invalid("material bank memory needs its original exclusive bounded source spectrum"));
        }
        // Validate the unchanged mechanical geometry without asking a
        // memoryless consumer to admit memory. The bank keeps the full source,
        // including explicit zero arms, and checks its complete coupled band.
        let mut elastic = string.clone();
        elastic.relaxation_bending = None;
        validate_incremental_string(&elastic)?;
        lower_backbone(cx, &string, gas, sample_rate_hz, 0.0)
    }

    /// Lower a homogeneous material specimen in an explicitly supplied gas.
    ///
    /// The supplied coherent state may describe air, another ideal gas, or an
    /// admitted mixture. Its actual density and viscosity feed the SAME drag
    /// owner used by other string excitations. Gas sound speed does not retune
    /// the dry structural modes. Material temperature, source band and all
    /// solid coefficients remain those of the original immutable specimen.
    ///
    /// Pass this same gas to any independently constructed acoustic body and
    /// observer; the card alone does not construct or certify a complete body.
    /// The caller retains the specimen's original source receipts.
    ///
    /// # Errors
    /// Refuses inconsistent gas, missing/foreign material loss, unsupported
    /// moving ends/polarizations/memory, out-of-band modes and cancellation.
    pub fn from_material_specimen_with_gas(
        cx: &Cx<'_>,
        specimen: &ResolvedStringSpecimen,
        gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = specimen.string();
        let law = string.kelvin_voigt_bending
            .ok_or_else(|| invalid("material bow card requires a resolved Kelvin-Voigt bending law"))?;
        if law.material_state_identity != Some(specimen.material().identity()) {
            return Err(invalid("bow bending loss must belong to the resolved material state"));
        }
        lower(cx, &string, gas, sample_rate_hz)
    }

    /// Lower the original section's intrinsic loss-angle or circular thermal
    /// loss through the same string modes, gas drag and support prescription.
    ///
    /// Each actual mode receives the cycle-equivalent rate from its own
    /// frequency and mass. Original thermal and optional excess sources remain
    /// in the section; missing/out-of-band selected loss refuses. No effective
    /// viscosity is substituted into elastic or viscous support shear.
    ///
    /// # Errors
    /// Inherited internal loss, unsupported geometry/state, source-band or
    /// time-domain admission failures, inconsistent gas, and cancellation.
    pub fn from_concentric_modal_loss_section(
        cx: &Cx<'_>,
        template: PrestressedString,
        section: &ConcentricModalLossSection,
        prestress: StringPrestress,
        gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = crate::material_assembly::concentric::with_concentric_section(
            template, section.elastic(), prestress)?;
        validate_incremental_string(&string)?;
        lower_backbone_with_modal_loss(cx, &string, gas, sample_rate_hz, 0.0, Some(section))
    }

    /// Lower a material-swappable layered member into the same bow runtime.
    ///
    /// Every annulus contributes its original density, elastic modulus and
    /// bending viscosity. The common section/support adapter derives linear
    /// mass, EI, viscous EI, diameter and tension together; the actual outer
    /// radius drives gas drag. The intersected source band admits every mode.
    /// The returned card uses the existing runtime's elastic AND viscous bridge
    /// shear, so changing a coating cannot leave the old support force behind.
    ///
    /// The template supplies span/resolution, while prestress explicitly chooses
    /// fixed tension, fixed extension, or target frequency. It must not carry
    /// an inherited internal loss law. This is a bonded-section bending model,
    /// not an assertion that a wrapped instrument string is a bonded cylinder.
    /// Source receipts remain in the borrowed section; no averaged source is minted.
    ///
    /// # Errors
    /// Shared section/support/loss refusals, unsupported string images,
    /// inconsistent gas, mode/source bandwidth limits or cancellation.
    pub fn from_concentric_material_section(
        cx: &Cx<'_>,
        template: PrestressedString,
        section: &ConcentricKelvinVoigtSection,
        prestress: StringPrestress,
        gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        poll(cx)?;
        let string = with_concentric_kelvin_voigt_section(template, section, prestress)?;
        lower(cx, &string, gas, sample_rate_hz)
    }
}

// All material sources lower through the same bounded numerical path. Kept
// private: a raw string carrying invented coefficients is not a sourced input.
fn lower(
    cx: &Cx<'_>, string: &PrestressedString, gas: &GasState, sample_rate_hz: u32,
) -> Result<BowedStringCard, AcousticRealizeError> {
    poll(cx)?;
    gas.validate_coherence().map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
    validate_incremental_string(string)?;
    let law = string.kelvin_voigt_bending
        .ok_or_else(|| invalid("material bow card requires a resolved Kelvin-Voigt bending law"))?;
    lower_backbone(cx, string, gas, sample_rate_hz, law.viscous_stiffness_n_m2_s)
}

fn lower_backbone(
    cx: &Cx<'_>, string: &PrestressedString, gas: &GasState, sample_rate_hz: u32,
    viscous_bending_n_m2_s: f64,
) -> Result<BowedStringCard, AcousticRealizeError> {
    lower_backbone_with_modal_loss(cx, string, gas, sample_rate_hz, viscous_bending_n_m2_s, None)
}

fn lower_backbone_with_modal_loss(
    cx: &Cx<'_>, string: &PrestressedString, gas: &GasState, sample_rate_hz: u32,
    viscous_bending_n_m2_s: f64, modal_loss: Option<&ConcentricModalLossSection>,
) -> Result<BowedStringCard, AcousticRealizeError> {
    poll(cx)?;
    gas.validate_coherence().map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
    let mut card = BowedStringCard {
        length_m: string.length_m,
        tension_n: string.tension_n,
        linear_density_kg_m: string.lin_density_kg_m,
        bending_stiffness_n_m2: string.bending_stiffness_n_m2,
        viscous_bending_n_m2_s,
        mode_count: string.n_modes,
        zetas: Vec::with_capacity(string.n_modes),
        sample_rate_hz,
    };
    let mut modes = Vec::with_capacity(string.n_modes);
    for k in 1..=string.n_modes {
        poll(cx)?;
        let omega = card.mode_omega_rad_s(k);
        let wave_number = k as f64 * core::f64::consts::PI / string.length_m;
        let rate = modal_loss.map(|section|
            section.damping_rate_per_s(omega, wave_number, card.linear_density_kg_m))
            .transpose().map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?
            .unwrap_or(0.0);
        let zeta = mode_zeta(string, omega, wave_number, gas)? + 0.5 * rate / omega;
        card.zetas.push(zeta);
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        });
    }
    ModalAcousticTimeModel::try_new(sample_rate_hz, modes,
        ModalAcousticTimeBudget::audible_reference())
        .map_err(|_| invalid("material bow modes violate time-domain admission"))?;
    card.validate().map_err(|_| invalid("material bow card is not physically representable"))?;
    poll(cx)?;
    Ok(card)
}
