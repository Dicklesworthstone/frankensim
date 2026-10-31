//! Shared SI membrane coefficients. Henry solubility here is concentration
//! divided by pressure; neither inverse Henry coefficients nor STP-volume
//! conventions are silently converted. These values contain no source claim.
use fs_qty::Dims;

pub const MOLAR_PERMEABILITY_DIMS: Dims = Dims([0, -1, 1, 0, 0, 1]);
pub const MOLAR_PERMEANCE_DIMS: Dims = Dims([-1, -1, 1, 0, 0, 1]);
pub const MEMBRANE_DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);
pub const MOLAR_HENRY_SOLUBILITY_DIMS: Dims = Dims([-2, -1, 2, 0, 0, 1]);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MembraneCoefficientError(&'static str);
impl MembraneCoefficientError {
    pub const fn message(self) -> &'static str { self.0 }
}
impl core::fmt::Display for MembraneCoefficientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for MembraneCoefficientError {}

fn nonnegative(value: f64, message: &'static str) -> Result<f64, MembraneCoefficientError> {
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(MembraneCoefficientError(message)) }
}

/// Dense-film diffusivity times molar Henry solubility, in mol/(m s Pa).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolarPermeability(pub(super) f64);
impl MolarPermeability {
    /// Zero is an explicitly impermeable species.
    pub fn try_new(value_mol_m_s_pa: f64) -> Result<Self, MembraneCoefficientError> {
        nonnegative(value_mol_m_s_pa, "molar permeability must be finite and nonnegative").map(Self)
    }
    pub const fn value_si(self) -> f64 { self.0 }
    /// P = D S in the explicitly declared dilute, concentration-based Henry law.
    pub fn from_diffusion(diffusivity: MembraneDiffusivity, solubility: MolarHenrySolubility)
        -> Result<Self, MembraneCoefficientError>
    {
        let value = diffusivity.0 * solubility.0;
        if !value.is_finite() || (diffusivity.0 > 0.0 && solubility.0 > 0.0 && value == 0.0) {
            return Err(MembraneCoefficientError("diffusivity times Henry solubility is unrepresentable"));
        }
        Ok(Self(value))
    }
    /// Thickness conversion is explicit; intrinsic porous k cannot enter here.
    pub fn permeance(self, thickness_m: f64) -> Result<MolarPermeance, MembraneCoefficientError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(MembraneCoefficientError("membrane thickness must be positive and finite"));
        }
        let value = self.0 / thickness_m;
        if !value.is_finite() || (self.0 > 0.0 && value == 0.0) {
            return Err(MembraneCoefficientError("membrane permeance is unrepresentable"));
        }
        Ok(MolarPermeance(value))
    }
}

/// Pressure-normalized molar flux, in mol/(m2 s Pa).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolarPermeance(pub(super) f64);
impl MolarPermeance {
    /// Explicit molar coefficient; no bulk, mass or STP-volume conversion.
    pub fn try_new(value_mol_m2_s_pa: f64) -> Result<Self, MembraneCoefficientError> {
        nonnegative(value_mol_m2_s_pa, "molar permeance must be finite and nonnegative").map(Self)
    }
    pub const fn value_si(self) -> f64 { self.0 }
}

/// Fickian diffusivity of one sorbed species through the dense phase [m2/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneDiffusivity(f64);
impl MembraneDiffusivity {
    pub fn try_new(value_m2_s: f64) -> Result<Self, MembraneCoefficientError> {
        nonnegative(value_m2_s, "membrane diffusivity must be finite and nonnegative").map(Self)
    }
    pub const fn value_si(self) -> f64 { self.0 }
}

/// Molar concentration / equilibrium partial pressure [mol/(m3 Pa)].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolarHenrySolubility(f64);
impl MolarHenrySolubility {
    pub fn try_new(value_mol_m3_pa: f64) -> Result<Self, MembraneCoefficientError> {
        nonnegative(value_mol_m3_pa, "molar Henry solubility must be finite and nonnegative").map(Self)
    }
    pub const fn value_si(self) -> f64 { self.0 }
}
