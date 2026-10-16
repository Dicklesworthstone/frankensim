//! Twist and transverse motion at the SAME material-string surface contact.
//!
//! For rotary line inertia Ip and circular stiffness GJ, fixed-rotation ends
//! give omega_n = (n*pi/L)*sqrt(GJ/Ip). Shear loss gives 2*zeta*omega =
//! (eta_s J/Ip)*(n*pi/L)^2. The original modal owner advances these coordinates.
//! Bow slip sees v + r*theta_dot; its transpose applies BOTH F and r*F. The
//! bridge here transmits transverse displacement only. No acoustic gain is
//! assigned to twist, and no extra friction work or loss is counted separately.
//!
//! This is small-angle twist of bonded circular layers, not helical winding,
//! interlayer slip, torque-bearing bridge contact, rotary gas drag, or large
//! axial strain. All bows use the same declared positive surface side. The
//! material band screens the retained linear modes, not nonlinear harmonics.
//! See Bavu, Smith & Wolfe, 'Torsional waves in a bowed string' (2005),
//! arXiv:physics/0505059, for the importance of surface rather than centre speed.

use super::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::torsion::ConcentricTorsionalSection;
use fs_material::state_point::MaterialStatePointError;

/// Frozen original bending specimen plus torsion from its SAME material states.
/// Private fields prevent attaching an unrelated section or nesting reductions.
#[derive(Clone, Debug)]
pub struct TorsionalStringSource {
    bending: MaterialStringSource,
    section: ConcentricTorsionalSection,
    mode_count: usize,
    total_modes: usize,
}
impl MaterialStringSource {
    /// Add small-angle torsion using named G [Pa] and shear viscosity [Pa s].
    ///
    /// The original homogeneous or concentric material bundles supply all new
    /// coefficients; there is no additional material query or independent mass.
    /// Source values must already be present. The original bending model and
    /// prestress are retained, not retuned to match a torsional resonance.
    /// All torsional endpoints are fixed in rotation over the complete span.
    /// Gas continues to load transverse motion; rotational gas loading is absent.
    ///
    /// # Errors
    /// Zero/oversized/nested mode selection, unrepresentable total dimension, or
    /// original source-law refusal. Frequency/clock admission occurs at assembly.
    pub fn with_torsion(self, shear_property: &str, viscosity_property: &str, modes: usize)
        -> Result<Self, MaterialStatePointError>
    {
        let invalid = |quantity| MaterialStatePointError::InvalidDerived { quantity };
        if modes == 0 || modes > 4096 || matches!(&self, Self::WithTorsion(_)) {
            return Err(invalid("torsional extension requires 1..=4096 modes and a non-nested source"));
        }
        let total_modes = self.modes().checked_add(modes)
            .ok_or_else(|| invalid("combined bending/torsion mode count overflow"))?;
        let elastic = match &self {
            Self::Homogeneous(s) => ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
                outer_radius_m: s.radius_m(), material: s.material().clone(),
            }])?,
            Self::Concentric { section, .. } => section.elastic().clone(),
            Self::WithTorsion(_) => return Err(invalid("nested torsion is not a new physical section")),
        };
        // The correction multiplies transverse GA only; it never changes GJ.
        let section = elastic.with_shear_modulus(shear_property, 1.0)?
            .with_torsional_viscosity(viscosity_property)?;
        Ok(Self::WithTorsion(Box::new(TorsionalStringSource {
            bending: self, section, mode_count: modes, total_modes,
        })))
    }

    pub(super) fn surface_shape(&self, card: &BowedStringCard, station: f64)
        -> Result<Vec<f64>, ReciprocalStringError>
    {
        let mut shapes = shape(card, station);
        if let Self::WithTorsion(twist) = self {
            let inertia = twist.section.rotary_inertia_kg_m();
            let radius = twist.section.shear().elastic().outer_radius_m();
            for k in 0..twist.mode_count {
                let value = radius * unit_shape(k, station, card.length_m, inertia);
                if !value.is_finite() { return Err(invalid("torsional surface shape is not representable")); }
                shapes.push(value);
            }
        }
        Ok(shapes)
    }
}
impl TorsionalStringSource {
    /// Unmodified original specimen/layers, bending loss and support prescription.
    #[must_use]
    pub const fn bending(&self) -> &MaterialStringSource { &self.bending }
    /// Original material section with independent shear loss and rotary inertia.
    #[must_use]
    pub const fn section(&self) -> &ConcentricTorsionalSection { &self.section }
    /// Number of appended twist modes; bending coordinates always come first.
    #[must_use]
    pub const fn mode_count(&self) -> usize { self.mode_count }
    pub(super) const fn total_modes(&self) -> usize { self.total_modes }

    pub(super) fn append_modes(&self, card: &BowedStringCard, modes: &mut Vec<ModalAcousticMode>, gate: &CancelGate)
        -> Result<(), ReciprocalStringError>
    {
        let inertia = self.section.rotary_inertia_kg_m();
        let speed = (self.section.stiffness_n_m2() / inertia).sqrt();
        let norm = (inertia * card.length_m * 0.5).sqrt();
        let radius = self.section.shear().elastic().outer_radius_m();
        if !speed.is_finite() || speed <= 0.0 || !norm.is_finite() || norm <= 0.0
            || !(radius / norm).is_finite() || radius / norm == 0.0 {
            return Err(invalid("torsion frequency and angular normalization must be representable"));
        }
        let (lo, hi) = self.section.omega_band_rad_s();
        for n in 1..=self.mode_count {
            poll(gate)?;
            let k = n as f64 * core::f64::consts::PI / card.length_m;
            let omega = k * speed;
            let zeta = (self.section.viscous_stiffness_n_m2_s() / inertia) * k * k / (2.0 * omega);
            if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi
                || !zeta.is_finite() || zeta < 0.0
                || (self.section.viscous_stiffness_n_m2_s() > 0.0 && zeta == 0.0) {
                return Err(invalid("actual torsional modes exceed material support or coefficient representability"));
            }
            modes.push(ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: zeta,
                pressure_per_modal_velocity: fs_math::c64::C64::ZERO });
        }
        Ok(())
    }
}
