//! One source-bound thin-membrane law for discrete and continuous ducts.
//! Geometry reduction is shared with the existing graph wall. Frequency-domain
//! memory uses GeneralizedMaxwell::modulus; no fitted loss factor is introduced.
use super::{AxialConstraint, WallMaterialLoss, WallPatch, property, invalid,
    YOUNG_VISCOSITY_PROPERTY};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::state_point::{ResolvedMaterialStatePoint, DENSITY_PROPERTY,
    YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY};
use fs_math::c64::C64;
use fs_phs::{WallPin, wall_specific_impedance};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_vfit::relaxation::RelaxationImpedanceSpec;

/// Immutable material, thickness, restraint and complete local loss law.
///
/// This is the existing isotropic, frozen-state, small-displacement membrane
/// limit. It does not include axial shell waves, pressure prestress, thermal
/// evolution, or nonlocal wall radiation. Every queried radius must satisfy
/// the original thickness/radius <= 0.05 admission.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalMaterialWall {
    material: ResolvedMaterialStatePoint,
    thickness_m: f64,
    axial: AxialConstraint,
    loss: WallMaterialLoss,
    band_hz: [f64; 2],
    density: f64,
    young: f64,
    viscosity: f64,
    axial_factor: f64,
}

impl LocalMaterialWall {
    /// Resolve the same original claims used by graph-based material walls.
    ///
    /// # Errors
    /// Missing/wrong-quantity/out-of-band material data, unsupported material
    /// dependence, nonpositive geometry or invalid isotropic coefficients.
    pub fn new(material: ResolvedMaterialStatePoint, thickness_m: f64,
        axial: AxialConstraint, loss: WallMaterialLoss, band_hz: [f64; 2])
        -> Result<Self, AcousticRealizeError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0
            || !band_hz.iter().all(|f| f.is_finite() && *f > 0.0)
            || band_hz[0] >= band_hz[1] || !(core::f64::consts::TAU * band_hz[1]).is_finite() {
            return Err(invalid("local material wall requires positive thickness and an ordered positive frequency band"));
        }
        let density = property(&material, DENSITY_PROPERTY, Density::DIMS, band_hz)?;
        let young = property(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS, band_hz)?;
        let nu = property(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]), band_hz)?;
        if density <= 0.0 || young <= 0.0 || nu <= -1.0 || nu >= 0.5 {
            return Err(invalid("isotropic wall needs rho,E>0 and -1<nu<0.5"));
        }
        let viscosity = match &loss {
            WallMaterialLoss::Elastic | WallMaterialLoss::Maxwell(_) => 0.0,
            WallMaterialLoss::KelvinVoigt =>
                property(&material, YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, band_hz)?,
        };
        if viscosity < 0.0 { return Err(invalid("wall Young viscosity must be nonnegative")); }
        if let WallMaterialLoss::Maxwell(spectrum) = &loss {
            spectrum.validate_source(&material)?;
            spectrum.validate_material(young, band_hz)?;
        }
        let axial_factor = match axial {
            AxialConstraint::StressFree => 1.0,
            AxialConstraint::StrainFixed => 1.0 / (1.0 - nu * nu),
        };
        Ok(Self { material, thickness_m, axial, loss, band_hz, density, young,
            viscosity, axial_factor })
    }

    /// Original material identity, queries and property receipts.
    #[must_use] pub const fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    #[must_use] pub const fn thickness_m(&self) -> f64 { self.thickness_m }
    #[must_use] pub const fn axial_constraint(&self) -> AxialConstraint { self.axial }
    #[must_use] pub const fn loss(&self) -> &WallMaterialLoss { &self.loss }
    #[must_use] pub const fn band_hz(&self) -> [f64; 2] { self.band_hz }

    /// Base mass/equilibrium/viscosity image for an actual cylindrical strip.
    /// Hereditary arms remain in `loss()`; this is not its complete impedance.
    pub fn base_patch(&self, radius: f64, length: f64) -> Result<WallPatch, AcousticRealizeError> {
        let h = self.thickness_m;
        if ![radius, length, h].iter().all(|x| x.is_finite() && *x > 0.0) || h / radius > 0.05 {
            return Err(invalid("local membrane wall requires positive geometry and thickness/radius <= 0.05"));
        }
        let factor = (h / radius) / radius * self.axial_factor;
        let patch = WallPatch { area_m2: core::f64::consts::TAU * radius * length,
            wall: WallPin { surface_density: self.density * h,
                stiffness_per_area: self.young * factor, resistance: self.viscosity * factor } };
        patch.impedance()?;
        if self.viscosity > 0.0 && patch.wall.resistance == 0.0 {
            return Err(invalid("wall material reduction lost its nonzero viscosity"));
        }
        Ok(patch)
    }

    pub(super) fn load(&self, patch: WallPatch) -> Result<RelaxationImpedanceSpec, AcousticRealizeError> {
        match &self.loss {
            WallMaterialLoss::Maxwell(spectrum) => spectrum.load(patch, self.young, self.band_hz),
            _ => RelaxationImpedanceSpec::new(patch.impedance()?, &[])
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string())),
        }
    }

    /// Complete surface impedance [Pa s/m] under exp(-i omega t).
    /// Storage and loss moduli come from the same source law as the graph's
    /// hereditary states, evaluated at the requested physical frequency.
    pub fn specific_impedance(&self, radius: f64, omega: f64)
        -> Result<C64, AcousticRealizeError> {
        if !omega.is_finite() || omega < core::f64::consts::TAU * self.band_hz[0]
            || omega > core::f64::consts::TAU * self.band_hz[1] {
            return Err(invalid("material wall frequency lies outside its declared source band"));
        }
        let patch = self.base_patch(radius, 1.0)?;
        let mut pin = patch.wall;
        if let WallMaterialLoss::Maxwell(spectrum) = &self.loss {
            // Retain the graph owner's admission of every positive physical arm.
            spectrum.load(patch, self.young, self.band_hz)?;
            // Reuse the complete generalized-Maxwell constitutive owner.
            let (storage, loss) = spectrum.material().modulus(omega);
            let factor = pin.stiffness_per_area / self.young;
            pin.stiffness_per_area = factor * storage;
            pin.resistance += factor * loss / omega;
        }
        wall_specific_impedance(&pin, omega)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    }
}
