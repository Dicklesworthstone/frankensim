//! Source-bound isotropic Lorentz polarization at fixed material conditions.
//!
//! Bound-charge inertia produces resonance, phase reversal and negative real
//! permittivity above a line. Every original pole retains its own positive
//! mechanical storage and nonnegative frictional heat. Debye relaxation and
//! Ohmic leakage remain separate laws; no fitted loss tangent substitutes for
//! the resonant dynamics. Optical wave consumers use `exp(+i omega t)`.
//!
//! Pole parameters are explicit SI material properties. Natural frequencies
//! must carry the angular-frequency quantity kind; inverse seconds alone do
//! not permit silently interpreting a frequency in Hz as radians per second.
//! This is a linear, isotropic, fixed-condition law, not measured validation,
//! nonlinear saturation, thermal parameter work, or a transient Maxwell solver.

use super::{
    ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError, MAX_DIELECTRIC_POLES, coefficient, finite,
    invalid,
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::c64::C64;
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

pub mod oscillator;
pub use oscillator::{
    LorentzError, LorentzHarmonicTransfer, LorentzPole, LorentzState, LorentzStep,
    PreparedLorentzHarmonic, PreparedLorentzStep,
};

impl From<LorentzError> for FieldStorageError {
    fn from(value: LorentzError) -> Self {
        Self::Invalid(match value {
            LorentzError::Invalid(what) => what,
            LorentzError::SingularResonance => {
                "lossless Lorentz resonance has no finite periodic response"
            }
            LorentzError::Evolution(_) => {
                "exact Lorentz polarization evolution is not representable"
            }
        })
    }
}

/// Source roles for one bound-charge process. No default frequency or loss.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LorentzPoleProperties {
    pub permittivity_increment: String,
    pub angular_frequency: String,
    pub damping_rate: String,
}

/// Complete original source bundle and its admitted finite oscillator bank.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLorentzMaterial {
    material: ResolvedMaterialStatePoint,
    high_frequency_property: String,
    properties: Vec<LorentzPoleProperties>,
    epsilon_infinity: f64,
    poles: Vec<LorentzPole>,
}

fn angular_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ))
}
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn admit_names(high: &str, poles: &[LorentzPoleProperties]) -> Result<(), FieldStorageError> {
    if high.trim().is_empty() || poles.len() > MAX_DIELECTRIC_POLES {
        return Err(invalid(
            "named high-frequency source and bounded Lorentz pole count required",
        ));
    }
    let mut names = vec![high];
    for pole in poles {
        for name in [
            &pole.permittivity_increment,
            &pole.angular_frequency,
            &pole.damping_rate,
        ] {
            if name.trim().is_empty() || names.contains(&name.as_str()) {
                return Err(invalid(
                    "Lorentz coefficients require distinct original source roles",
                ));
            }
            names.push(name);
        }
    }
    Ok(())
}

fn parameter(
    material: &ResolvedMaterialStatePoint,
    name: &str,
    quantity: QuantitySpec,
) -> Result<f64, FieldStorageError> {
    let property = material
        .property(name)
        .ok_or(invalid("missing original Lorentz parameter"))?;
    if property.requirement().quantity() != quantity
        || !property.is_static_at_fixed_external_conditions()
        || !property.value_si().is_finite()
        || property.value_si() < 0.0
    {
        return Err(invalid(
            "Lorentz parameters require original nonnegative SI sources at fixed external conditions",
        ));
    }
    Ok(property.value_si())
}

impl ResolvedLorentzMaterial {
    /// Keep the complete source bundle, including claims used by other physics.
    /// The empty bank is exactly the instantaneous lossless dielectric.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        high_frequency_property: &str,
        properties: &[LorentzPoleProperties],
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, properties)?;
        let epsilon_infinity = coefficient(
            &material,
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?;
        let mut poles = Vec::with_capacity(properties.len());
        let mut static_epsilon = epsilon_infinity;
        for names in properties {
            let increment = coefficient(
                &material,
                &names.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?;
            let omega = parameter(&material, &names.angular_frequency, angular_quantity())?;
            let gamma = parameter(
                &material,
                &names.damping_rate,
                QuantitySpec::dimensional(RATE),
            )?;
            poles.push(LorentzPole::try_new(increment, omega, gamma)?);
            static_epsilon = finite(static_epsilon + increment)?;
        }
        Ok(Self {
            material,
            high_frequency_property: high_frequency_property.into(),
            properties: properties.to_vec(),
            epsilon_infinity,
            poles,
        })
    }

    /// Resolve original material data. Angular/cyclic frequency mismatches,
    /// missing loss data and frequency-dependent coefficients refuse.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        high_frequency_property: &str,
        properties: &[LorentzPoleProperties],
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, properties)?;
        let mut requirements = vec![ScalarPropertyRequirement::try_new(
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?];
        for names in properties {
            requirements.push(ScalarPropertyRequirement::try_new(
                &names.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
            requirements.push(ScalarPropertyRequirement::try_with_quantity(
                &names.angular_frequency,
                angular_quantity(),
                ScalarAdmissibility::StrictlyPositive,
            )?);
            requirements.push(ScalarPropertyRequirement::try_new(
                &names.damping_rate,
                RATE,
                ScalarAdmissibility::NonNegative,
            )?);
        }
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            high_frequency_property,
            properties,
        )
    }

    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub fn high_frequency_property(&self) -> &str {
        &self.high_frequency_property
    }
    pub fn pole_properties(&self) -> &[LorentzPoleProperties] {
        &self.properties
    }
    pub const fn high_frequency_permittivity_f_m(&self) -> f64 {
        self.epsilon_infinity
    }
    pub fn poles(&self) -> &[LorentzPole] {
        &self.poles
    }

    /// Complex absolute permittivity under `exp(+i omega t)`. Negative real
    /// permittivity is physical; negative damping and lossless poles refuse.
    pub fn permittivity_f_m(&self, omega: f64) -> Result<C64, FieldStorageError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(invalid("finite nonnegative angular frequency required"));
        }
        let mut real = fs_math::dd::Dd::from_f64(self.epsilon_infinity);
        let mut imaginary = fs_math::dd::Dd::ZERO;
        for &pole in &self.poles {
            let contribution = pole.susceptibility_f_m(omega)?;
            real = real + fs_math::dd::Dd::from_f64(contribution.re);
            imaginary = imaginary + fs_math::dd::Dd::from_f64(contribution.im);
        }
        Ok(C64::new(
            finite(real.to_f64())?,
            finite(imaginary.to_f64())?,
        ))
    }

    /// Prepare the complete original pole bank for repeated field intervals.
    /// The returned transitions reuse the exact shared scalar propagator.
    pub fn prepare_step(
        &self,
        duration_s: f64,
    ) -> Result<Vec<PreparedLorentzStep>, FieldStorageError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid("positive finite Lorentz interval required"));
        }
        self.poles
            .iter()
            .map(|p| p.prepare_step(duration_s).map_err(Into::into))
            .collect()
    }
}
