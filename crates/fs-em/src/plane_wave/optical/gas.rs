//! Original molecular refraction in homogeneous monochromatic Maxwell waves.
//!
//! The material owner evaluates its molecular mixture or atmospheric law
//! once. With exp(+i omega t), epsilon=n_phase²/(mu*c0²); group
//! index is a distinct propagation derivative, never the phase constitutive
//! coefficient. The caller explicitly chooses a lossless coherent-medium
//! approximation: molecular scattering is NOT converted into absorption.

use super::{FREQUENCY_WAVELENGTH_NM, ReflectionError, VACUUM_SPEED_M_S, dispersion, positive};
use fs_material::{
    air_refraction::{CiddorAir, WAVELENGTH_RANGE_NM as AIR_WAVELENGTH_RANGE_NM},
    field_storage::linear_magnetic_permeability,
    gas_scattering::{
        RayleighGasState,
        composition::{MolecularRefractionModel, PartialDensityBasis},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, det};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};

/// Original immutable gas sources, with no reconstructed species or fluid EOS.
#[derive(Clone, Copy, Debug)]
pub enum MolecularGasOpticalSource<'a> {
    /// The exact accepted gas state and its fixed composition/optical sources.
    FixedComposition(&'a RayleighGasState),
    /// Actual homogeneous partial densities in the original model's row order.
    /// The borrowed values are densities, not fractions to be normalized.
    PartialDensities {
        model: &'a MolecularRefractionModel,
        partial_densities: &'a [f64],
        basis: PartialDensityBasis,
    },
    /// Original atmospheric composition and homogeneous thermodynamic state.
    CiddorAir {
        composition: CiddorAir,
        temperature_k: f64,
        pressure_pa: f64,
    },
}

/// One wavelength of an original molecular gas plus explicitly selected mu.
///
/// The magnetic material is a SEPARATE source whose applicability the caller
/// declares; it is not claimed to own the gas chemistry or optical references.
/// No vacuum permeability, new EOS, fitted optical index or extinction is
/// inserted. Existing native interfaces, finite films and angular spectra
/// retain this response. Its finite-band phase/group data do not establish a
/// causal transient reservoir, molecular attenuation or measured accuracy.
///
/// Preparation is bounded by the material owner's component cap and allocates
/// nothing. All subsequent queries use cached epsilon and exact frequency
/// admission, so each angular mode reuses the same molecular evaluation.
/// Local group delay/GVD additionally reuse the cached source derivatives;
/// two-sided frequency derivatives require an interior point of the source band.
#[derive(Clone, Copy, Debug)]
pub struct MolecularGasOpticalResponse<'a> {
    source: MolecularGasOpticalSource<'a>,
    magnetic: &'a ResolvedMaterialStatePoint,
    omega_rad_s: f64,
    wavelength_nm: f64,
    phase_index: f64,
    group_index: f64,
    scaled_wavelength_derivatives: [f64; 2],
    permeability_h_m: f64,
    permittivity_f_m: C64,
}

impl<'a> MolecularGasOpticalResponse<'a> {
    /// Choose the lossless coherent approximation for the original gas state.
    /// Magnetic-source conditions must match the gas's actual T, p and rho.
    /// Typed static absolute temperature, pressure and mass concentration are
    /// mapped by meaning; legacy untyped/dimensional `T` remains kelvin.
    /// Other coordinates refuse. The original constant-polarizability visible
    /// model has equal phase/group index and keeps its 380--780 nm domain.
    pub fn try_lossless_fixed_composition(
        gas: &'a RayleighGasState,
        magnetic: &'a ResolvedMaterialStatePoint,
        omega_rad_s: f64,
    ) -> Result<Self, ReflectionError> {
        admit_fixed_conditions(gas, magnetic)?;
        let wavelength_nm = wavelength(omega_rad_s, [380.0, 780.0])?;
        let sample = gas
            .refractive_index(wavelength_nm)
            .map_err(|_| bad("molecular gas is outside its original refraction domain"))?;
        Self::prepare(
            MolecularGasOpticalSource::FixedComposition(gas),
            magnetic,
            omega_rad_s,
            wavelength_nm,
            sample.phase_index,
            sample.phase_index,
            [0.0; 2],
        )
    }

    /// Choose the lossless coherent approximation for these original partial
    /// densities, including the model's full source-supported dispersion.
    ///
    /// This source does not own temperature or pressure. The magnetic source
    /// must therefore be unconditioned: an unknown condition is never inferred
    /// from an ideal-gas guess or from the model's reference mixture. Density
    /// basis, row count, dilute-use bounds and spectral support are checked by
    /// the original molecular owner, including vacuum and infrared sources.
    pub fn try_lossless_partial_densities(
        model: &'a MolecularRefractionModel,
        partial_densities: &'a [f64],
        basis: PartialDensityBasis,
        magnetic: &'a ResolvedMaterialStatePoint,
        omega_rad_s: f64,
    ) -> Result<Self, ReflectionError> {
        if magnetic.query_point().into_iter().next().is_some() {
            return Err(bad(
                "partial-density gas optics requires an unconditioned magnetic source",
            ));
        }
        let wavelength_nm = wavelength(omega_rad_s, model.wavelength_range_nm())?;
        let sample = model
            .at_wavelength(wavelength_nm, basis)
            .and_then(|spectral| spectral.refractive_index_jet(partial_densities))
            .map_err(|_| {
                bad("molecular partial densities or wavelength are outside their source domain")
            })?;
        Self::prepare(
            MolecularGasOpticalSource::PartialDensities {
                model,
                partial_densities,
                basis,
            },
            magnetic,
            omega_rad_s,
            wavelength_nm,
            sample.index.phase_index,
            sample.index.group_index,
            sample.scaled_wavelength_derivatives,
        )
    }

    /// Choose the lossless coherent approximation for the original Ciddor
    /// atmospheric-air law. Temperature [K], absolute pressure [Pa], water
    /// mole fraction and dry-air CO2 remain fixed during spectral derivatives.
    /// The named material law supplies its own optical density correction;
    /// no fluid EOS, mass density, absorption or scattering is inferred.
    /// Magnetic conditions may name only matching static temperature/pressure
    /// (or legacy kelvin `T`); composition remains owned by the air source.
    pub fn try_lossless_ciddor_air(
        composition: CiddorAir,
        temperature_k: f64,
        pressure_pa: f64,
        magnetic: &'a ResolvedMaterialStatePoint,
        omega_rad_s: f64,
    ) -> Result<Self, ReflectionError> {
        let wavelength_nm = wavelength(omega_rad_s, AIR_WAVELENGTH_RANGE_NM)?;
        let sample = composition
            .prepare(wavelength_nm)
            .and_then(|prepared| prepared.sample(temperature_k, pressure_pa))
            .map_err(|_| bad("Ciddor air state is outside its original refraction domain"))?;
        admit_air_conditions(temperature_k, pressure_pa, magnetic)?;
        Self::prepare(
            MolecularGasOpticalSource::CiddorAir {
                composition,
                temperature_k,
                pressure_pa,
            },
            magnetic,
            omega_rad_s,
            wavelength_nm,
            sample.phase_index,
            sample.group_index,
            [
                wavelength_nm * sample.wavelength_derivative_per_nm,
                wavelength_nm.powi(2) * sample.wavelength_second_derivative_per_nm2,
            ],
        )
    }

    fn prepare(
        source: MolecularGasOpticalSource<'a>,
        magnetic: &'a ResolvedMaterialStatePoint,
        omega_rad_s: f64,
        wavelength_nm: f64,
        phase_index: f64,
        group_index: f64,
        scaled_wavelength_derivatives: [f64; 2],
    ) -> Result<Self, ReflectionError> {
        let permeability_h_m = linear_magnetic_permeability(magnetic)?;
        let scale = positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        let permittivity_f_m = dispersion::permittivity(C64::from_re(phase_index), scale)
            .ok_or_else(|| bad("molecular gas optical permittivity is unrepresentable"))?;
        Ok(Self {
            source,
            magnetic,
            omega_rad_s,
            wavelength_nm,
            phase_index,
            group_index,
            scaled_wavelength_derivatives,
            permeability_h_m,
            permittivity_f_m,
        })
    }

    pub const fn source(self) -> MolecularGasOpticalSource<'a> {
        self.source
    }

    /// The explicitly selected magnetic source ONLY. Molecular optics and the
    /// actual gas state/partial densities remain owned by [`Self::source`].
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.magnetic
    }

    pub const fn omega_rad_s(self) -> f64 {
        self.omega_rad_s
    }
    pub const fn wavelength_nm(self) -> f64 {
        self.wavelength_nm
    }
    pub const fn phase_index(self) -> f64 {
        self.phase_index
    }
    pub const fn group_index(self) -> f64 {
        self.group_index
    }
    pub const fn permeability_h_m(self) -> f64 {
        self.permeability_h_m
    }

    /// Cached absolute phase n-i*k, k=0 under the explicit approximation.
    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        self.admit_frequency(omega_rad_s)?;
        Ok(C64::from_re(self.phase_index))
    }

    /// Cached absolute SI phase permittivity, including the selected mu.
    pub fn permittivity_f_m(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        self.admit_frequency(omega_rad_s)?;
        Ok(self.permittivity_f_m)
    }

    pub(in crate::plane_wave) fn frequency_jet(
        self,
        omega: f64,
    ) -> Result<super::frequency::OpticalFrequencyJet, ReflectionError> {
        self.admit_frequency(omega)?;
        let band = match self.source {
            MolecularGasOpticalSource::FixedComposition(_) => [380.0, 780.0],
            MolecularGasOpticalSource::CiddorAir { .. } => AIR_WAVELENGTH_RANGE_NM,
            MolecularGasOpticalSource::PartialDensities { model, .. } => {
                model.wavelength_range_nm()
            }
        };
        super::frequency::interior(omega, FREQUENCY_WAVELENGTH_NM, band)?;
        let scale = positive((VACUUM_SPEED_M_S * det::sqrt(self.permeability_h_m)).recip())?;
        let [first, second] = self.scaled_wavelength_derivatives;
        super::frequency::finish(
            C64::from_re(self.phase_index),
            scale,
            C64::from_re(first),
            C64::from_re(second),
            C64::from_re(self.group_index),
        )
    }

    fn admit_frequency(self, omega: f64) -> Result<(), ReflectionError> {
        if omega != self.omega_rad_s {
            return Err(bad(
                "molecular gas optics was prepared at another angular frequency",
            ));
        }
        Ok(())
    }

    pub(in crate::plane_wave) fn same_source(self, other: Self) -> bool {
        // Frequency-dependent caches differ across a polychromatic exposure;
        // only the original material and homogeneous state must stay fixed.
        core::ptr::eq(self.magnetic, other.magnetic)
            && match (self.source, other.source) {
                (
                    MolecularGasOpticalSource::FixedComposition(a),
                    MolecularGasOpticalSource::FixedComposition(b),
                ) => core::ptr::eq(a, b),
                (
                    MolecularGasOpticalSource::PartialDensities {
                        model: a,
                        partial_densities: ar,
                        basis: ab,
                    },
                    MolecularGasOpticalSource::PartialDensities {
                        model: b,
                        partial_densities: br,
                        basis: bb,
                    },
                ) => core::ptr::eq(a, b) && ab == bb && ar == br,
                (
                    MolecularGasOpticalSource::CiddorAir {
                        composition: a,
                        temperature_k: at,
                        pressure_pa: ap,
                    },
                    MolecularGasOpticalSource::CiddorAir {
                        composition: b,
                        temperature_k: bt,
                        pressure_pa: bp,
                    },
                ) => a == b && at == bt && ap == bp,
                _ => false,
            }
    }
}

fn wavelength(omega: f64, band: [f64; 2]) -> Result<f64, ReflectionError> {
    dispersion::wavelength_from_frequency(omega, FREQUENCY_WAVELENGTH_NM, band).ok_or_else(|| {
        bad("molecular gas frequency is outside its original vacuum-wavelength band")
    })
}

fn admit_fixed_conditions(
    gas: &RayleighGasState,
    magnetic: &ResolvedMaterialStatePoint,
) -> Result<(), ReflectionError> {
    let state = gas.gas();
    for (axis, value) in magnetic.query_point() {
        let quantity = magnetic.axis_quantities().get(axis).copied();
        let expected = match quantity.and_then(QuantitySpec::semantic_type) {
            Some(ty) if ty.form() == ValueForm::Static => match ty.kind() {
                QuantityKind::AbsoluteTemperature => state.temperature,
                QuantityKind::Pressure => state.pressure,
                QuantityKind::MassConcentration => state.density,
                _ => {
                    return Err(bad(
                        "molecular gas magnetic source has an unmapped condition",
                    ));
                }
            },
            None if axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS)) =>
            {
                state.temperature
            }
            _ => {
                return Err(bad(
                    "molecular gas magnetic source needs original static conditions",
                ));
            }
        };
        if *value != expected {
            return Err(bad("molecular gas and magnetic source conditions disagree"));
        }
    }
    Ok(())
}

fn admit_air_conditions(
    temperature_k: f64,
    pressure_pa: f64,
    magnetic: &ResolvedMaterialStatePoint,
) -> Result<(), ReflectionError> {
    for (axis, value) in magnetic.query_point() {
        let quantity = magnetic.axis_quantities().get(axis).copied();
        let expected = match quantity.and_then(QuantitySpec::semantic_type) {
            Some(ty) if ty.form() == ValueForm::Static => match ty.kind() {
                QuantityKind::AbsoluteTemperature => temperature_k,
                QuantityKind::Pressure => pressure_pa,
                _ => return Err(bad("Ciddor air magnetic source has an unmapped condition")),
            },
            None if axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS)) =>
            {
                temperature_k
            }
            _ => {
                return Err(bad(
                    "Ciddor air magnetic source needs original static conditions",
                ));
            }
        };
        if *value != expected {
            return Err(bad("Ciddor air and magnetic source conditions disagree"));
        }
    }
    Ok(())
}

fn bad(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}

#[cfg(test)]
mod tests;
