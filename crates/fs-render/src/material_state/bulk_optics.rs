//! One conversion of material optical sources shared by surfaces and particles.

use super::*;

pub(super) fn conductor(
    source: &VisibleConductorStatePoint,
) -> Result<ConductorOptics, MaterialOpticalError> {
    let mut samples = [ConductorIorSample::try_new(380.0, 1.0, 1.0)?; 9];
    for (slot, sample) in samples.iter_mut().zip(source.samples()) {
        *slot = ConductorIorSample::try_new(sample.wavelength_nm, sample.eta, sample.k)?;
    }
    Ok(ConductorOptics::try_new(
        samples,
        ConductorSource::try_new(
            source.resolved().identity(),
            ConductorDataStatus::MaterialStateClaims,
        )?,
    )?)
}

pub(super) fn cauchy(
    source: &VisibleDielectricStatePoint,
) -> Result<DielectricGlass, MaterialOpticalError> {
    let [a, b_m2, c_m4] = source.cauchy_coefficients_si();
    Ok(DielectricGlass::new(
        CauchyIor::try_new(a, b_m2 * 1.0e12, c_m4 * 1.0e24)?,
        MaterialOpticalBinding::dielectric_absorption(source.absorption())?,
        GlassProvenance::Custom,
    ))
}

pub(super) fn tabulated(
    source: &VisibleTabulatedDielectricStatePoint,
) -> Result<DielectricGlass, MaterialOpticalError> {
    let mut indices = [DielectricIorSample::try_new(380.0, 1.0)?; 9];
    let mut absorption = [(380.0, 0.0); 9];
    for ((index, coefficient), sample) in indices
        .iter_mut()
        .zip(absorption.iter_mut())
        .zip(source.samples())
    {
        *index = DielectricIorSample::try_new(sample.wavelength_nm, sample.phase_index)?;
        *coefficient = (sample.wavelength_nm, sample.absorption_per_m);
    }
    Ok(DielectricGlass::with_ior(
        DielectricIor::Tabulated(TabulatedIor::try_new(indices)?),
        BeerLambertAbsorption::try_from_spectral_extinction(&absorption)?,
        GlassProvenance::Custom,
    ))
}

pub(super) fn sellmeier(
    source: &VisibleSellmeierStatePoint,
) -> Result<DielectricGlass, MaterialOpticalError> {
    let (strengths, resonances_m2) = source.sellmeier_coefficients_si();
    Ok(DielectricGlass::new(
        SellmeierIor::try_new(strengths, resonances_m2.map(|value| value * 1.0e12))?,
        MaterialOpticalBinding::dielectric_absorption(source.absorption())?,
        GlassProvenance::Custom,
    ))
}

#[cfg(feature = "volumes")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum BulkOpticalLaw {
    Conductor(ConductorOptics),
    Dielectric(DielectricGlass),
}

#[cfg(feature = "volumes")]
impl BulkOpticalLaw {
    pub(super) fn from_state(
        source: &VisibleOpticalStatePoint,
    ) -> Result<Self, MaterialOpticalError> {
        match source {
            VisibleOpticalStatePoint::Conductor(value) => conductor(value).map(Self::Conductor),
            VisibleOpticalStatePoint::Dielectric(value) => cauchy(value).map(Self::Dielectric),
            VisibleOpticalStatePoint::TabulatedDielectric(value) => {
                tabulated(value).map(Self::Dielectric)
            }
            VisibleOpticalStatePoint::DielectricSellmeier(value) => {
                sellmeier(value).map(Self::Dielectric)
            }
        }
    }

    /// Bulk intensity attenuation is alpha = 4*pi*k/lambda_vacuum, with
    /// lambda in metres. This does not infer alpha from scattering losses.
    pub(super) fn complex_index(
        self,
        wavelength_nm: f64,
    ) -> Result<fs_math::c64::C64, MaterialOpticalError> {
        match self {
            Self::Conductor(optics) => {
                let value = optics.complex_ior(wavelength_nm)?;
                Ok(fs_math::c64::C64::new(value.eta(), value.k()))
            }
            Self::Dielectric(glass) => {
                let n = glass.ior().eval(wavelength_nm)?;
                let alpha = glass.absorption().extinction_per_m(wavelength_nm)?;
                let k = alpha * (wavelength_nm * 1.0e-9 / (4.0 * core::f64::consts::PI));
                if !k.is_finite() || (alpha > 0.0 && k == 0.0) {
                    return Err(DielectricError::InvalidAbsorption.into());
                }
                Ok(fs_math::c64::C64::new(n, k))
            }
        }
    }
}
