//! Independent angular ensembles on laterally homogeneous isotropic slabs.
//!
//! A diffuse hemisphere has constant incident radiance, so its flux measure is
//! `2 mu dmu`, with `mu=cos(theta)` in [0,1]. Azimuth is already integrated out.
//! This is the plane-wave random-incidence model, not a finite room, finite
//! absorber, nearby monopole, or correlated cross-angle field. See Sgard,
//! Atalla and Robin (2024), section 2.2.3,
//! <https://doi.org/10.3389/facou.2024.1414356>.

use super::{
    MAX_POROUS_SPECTRUM_BINS, PorousCovarianceBin, PorousHeatingError, PressureCovariance,
};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, dd::Dd, det};
use fs_thermal::enthalpy::DepositionError;

/// One stationary independent spatial/frequency mode on the complete specimen.
/// Its two exterior faces may be partially coherent with each other. Distinct
/// entries have ZERO cross-covariance, even when their frequencies are equal.
/// Isotropic, parallel layers require only |k_t|; independent equal-|k_t|
/// azimuths must be combined in one covariance before native preparation.
#[derive(Clone, Copy, Debug)]
pub struct PorousAngularMode {
    /// Integrated RMS pressure covariance and original fixed exterior reference.
    pub bin: PorousCovarianceBin,
    /// Conserved nonnegative tangential wavenumber [rad/m].
    pub tangential_wavenumber_rad_m: f64,
}

/// One positive quadrature node for the hemisphere's `dmu`, not `dtheta` or
/// `2 mu dmu`. The builder applies the projected-flux cosine exactly once.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseAngularNode {
    /// `mu=cos(theta)`, strictly above zero and at most one.
    pub cosine: f64,
    /// Positive integration weight in `dmu` over [0,1].
    pub weight_mu: f64,
}

/// Fixed lossless exterior and prescribed incoming diffuse irradiance.
/// Both exterior faces use the same medium; their fields are independent.
/// Irradiance is normal incident acoustic power per FULL specimen area [W/m²]
/// at this frequency/band, not RMS pressure, along-ray intensity or radiance.
/// It stays fixed while the INTERNAL pore gas evolves thermally.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseHemisphere {
    exterior: FluidAcousticState,
    omega: f64,
    area_m2: f64,
    irradiance_w_m2: [f64; 2],
}

impl DiffuseHemisphere {
    /// Admit one original exterior acoustic state and already band-integrated
    /// left/right normal irradiances. The exterior is modeled as lossless for
    /// this plane-wave excitation; the slab still uses its complete pore losses.
    pub fn try_new<M: FluidAcoustics + ?Sized>(
        exterior: &M,
        angular_frequency_rad_s: f64,
        area_m2: f64,
        irradiance_w_m2: [f64; 2],
    ) -> Result<Self, PorousHeatingError> {
        if ![angular_frequency_rad_s, area_m2]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
            || !irradiance_w_m2
                .into_iter()
                .all(|x| x.is_finite() && x >= 0.0)
        {
            return Err(PorousHeatingError::Invalid(
                "diffuse illumination needs positive frequency/area and nonnegative finite irradiance",
            ));
        }
        let exterior = exterior.acoustic_state()?;
        let wavenumber = angular_frequency_rad_s / exterior.sound_speed_m_s();
        if !wavenumber.is_finite() || wavenumber <= 0.0 {
            return Err(PorousHeatingError::Invalid(
                "diffuse exterior wavenumber is unrepresentable",
            ));
        }
        Ok(Self {
            exterior,
            omega: angular_frequency_rad_s,
            area_m2,
            irradiance_w_m2,
        })
    }

    /// Original checked exterior acoustic state, independent of internal gas.
    pub const fn exterior(self) -> FluidAcousticState {
        self.exterior
    }
    /// Prescribed normal incoming irradiances [W/m²], left/right respectively.
    pub const fn irradiance_w_m2(self) -> [f64; 2] {
        self.irradiance_w_m2
    }

    /// Build a bounded full-hemisphere angular quadrature. Positive `dmu` weights
    /// must integrate both 1 and 2*mu to unity within 64 binary64 epsilons;
    /// no empirical grazing cutoff or silent renormalization is applied.
    /// For example, midpoint nodes `(i+0.5)/n` with weights `1/n` obey both
    /// moments; higher-order positive rules may be supplied by the caller.
    /// Angular refinement is still required for the actual absorption function.
    ///
    /// With `Z0=rho0*c0`, normal reference `R=Z0/(area*mu)` and irradiance E,
    /// a node carries `area*E*2*mu*w` watts. Its RMS pressure covariance is
    /// therefore `R * watts = 2*Z0*E*w`: the grazing cosine cancels ANALYTICALLY.
    /// Pressure variance stays finite as mu approaches zero, and power tends
    /// to zero without dividing a tiny assigned flux by mu. Exact grazing
    /// (`mu=0`) has a singular real-wave reference and refuses, even at zero E.
    /// These independent hemispheres have no left/right cross covariance.
    /// Unrepresentable references or lost positive pressure variance refuse;
    /// numerical underflow is never interpreted as a zero-power input.
    pub fn modes(
        self,
        nodes: &[DiffuseAngularNode],
        maximum_modes: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<PorousAngularMode>, PorousHeatingError> {
        checkpoint()?;
        if maximum_modes == 0
            || maximum_modes > MAX_POROUS_SPECTRUM_BINS
            || nodes.is_empty()
            || nodes.len() > maximum_modes
        {
            return Err(DepositionError::Budget("diffuse angular quadrature modes").into());
        }
        let mut solid_angle = Dd::ZERO;
        let mut flux = Dd::ZERO;
        for node in nodes {
            checkpoint()?;
            if !node.cosine.is_finite()
                || node.cosine <= 0.0
                || node.cosine > 1.0
                || !node.weight_mu.is_finite()
                || node.weight_mu <= 0.0
            {
                return Err(PorousHeatingError::Invalid(
                    "diffuse nodes need 0 < cosine <= 1 and positive finite dmu weight",
                ));
            }
            solid_angle = solid_angle + Dd::from_f64(node.weight_mu);
            flux = flux + Dd::from_f64(2.0 * node.cosine) * Dd::from_f64(node.weight_mu);
        }
        if [solid_angle.to_f64(), flux.to_f64()]
            .into_iter()
            .any(|value| !value.is_finite() || (value - 1.0).abs() > 64.0 * f64::EPSILON)
        {
            return Err(PorousHeatingError::Invalid(
                "diffuse full-hemisphere dmu and cosine-flux weights must each normalize to one",
            ));
        }
        let impedance = self.exterior.characteristic_impedance_pa_s_m();
        let wavenumber = self.omega / self.exterior.sound_speed_m_s();
        let mut modes = Vec::with_capacity(nodes.len());
        for node in nodes {
            checkpoint()?;
            let reference = (impedance / self.area_m2) / node.cosine;
            let kt = wavenumber * det::sqrt((1.0 - node.cosine) * (1.0 + node.cosine));
            if !reference.is_finite() || reference <= 0.0 || !kt.is_finite() {
                return Err(PorousHeatingError::Invalid(
                    "diffuse normal wave reference or tangential mode is unrepresentable",
                ));
            }
            let variance = self
                .irradiance_w_m2
                .map(|e| impedance * (2.0 * node.weight_mu * e));
            if variance
                .iter()
                .zip(self.irradiance_w_m2)
                .any(|(&value, irradiance)| irradiance > 0.0 && value == 0.0)
            {
                return Err(PorousHeatingError::Invalid(
                    "positive diffuse pressure variance was lost to numerical underflow",
                ));
            }
            let covariance = PressureCovariance::try_new([
                C64::from_re(variance[0]),
                C64::ZERO,
                C64::ZERO,
                C64::from_re(variance[1]),
            ])?;
            modes.push(PorousAngularMode {
                bin: PorousCovarianceBin {
                    angular_frequency_rad_s: self.omega,
                    area_m2: self.area_m2,
                    reference_impedance_pa_s_m3: reference,
                    covariance,
                },
                tangential_wavenumber_rad_m: kt,
            });
        }
        checkpoint()?;
        Ok(modes)
    }
}
