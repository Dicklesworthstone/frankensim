//! Two original spectral index axes of a passive uniaxial crystal.
//!
//! Each eigenaxis retains its explicitly bound wavelength/property roles from
//! SpectralComplexIndexStatePoint. Both MUST refer to the same complete
//! material snapshot, including absolute scalar permeability. Ordinary means
//! E perpendicular to the optic axis; extraordinary is the principal index
//! for E along it. The latter is not the angle-dependent extraordinary phase
//! index, which the Maxwell interface computes from the index ellipsoid.
//!
//! Selecting these two roles explicitly asserts a uniaxial constitutive model;
//! no material name, surface roughness, or single isotropic spectrum supplies
//! it. The geometry consumer owns optic-axis orientation independently.

use super::{ReflectionError, SpectralOpticalResponse, SpectralOpticalSource, positive};
use fs_material::state_point::{ResolvedMaterialStatePoint, SpectralComplexIndexStatePoint};
use fs_math::c64::C64;

/// Borrowed original ordinary/extraordinary spectra with same-state SI mu.
/// Source support, interpolation gaps and optical coordinate admission reuse
/// the existing spectral Maxwell response, without resampling or allocation.
#[derive(Clone, Copy, Debug)]
pub struct UniaxialOpticalResponse<'a> {
    ordinary: &'a SpectralComplexIndexStatePoint,
    extraordinary: &'a SpectralComplexIndexStatePoint,
    ordinary_response: SpectralOpticalResponse<'a>,
    extraordinary_response: SpectralOpticalResponse<'a>,
    wavelength_bounds_nm: [f64; 2],
}

impl<'a> UniaxialOpticalResponse<'a> {
    /// Bind explicit ordered optical roles from one original complete source.
    /// The maximum interpolation-knot span is authored in vacuum metres.
    /// Equal spectra are admitted as an explicit isotropic limit.
    ///
    /// Positive extinction is never discarded: every phase_indices query
    /// requires both original axes to be exactly lossless at that frequency.
    /// No extrapolation or second independently authored permeability is used.
    pub fn try_new(
        ordinary: &'a SpectralComplexIndexStatePoint,
        extraordinary: &'a SpectralComplexIndexStatePoint,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        if ordinary.resolved().identity() != extraordinary.resolved().identity() {
            return Err(ReflectionError::Invalid(
                "uniaxial optical axes require the same complete material snapshot",
            ));
        }
        let ordinary_response = SpectralOpticalResponse::try_new(
            SpectralOpticalSource::ComplexIndex(ordinary),
            maximum_knot_span_m,
        )?;
        let extraordinary_response = SpectralOpticalResponse::try_new(
            SpectralOpticalSource::ComplexIndex(extraordinary),
            maximum_knot_span_m,
        )?;
        let o = ordinary.wavelength_bounds_nm();
        let e = extraordinary.wavelength_bounds_nm();
        let wavelength_bounds_nm = [o[0].max(e[0]), o[1].min(e[1])];
        if wavelength_bounds_nm[0] > wavelength_bounds_nm[1] {
            return Err(ReflectionError::Invalid(
                "uniaxial optical axes have no common original wavelength support",
            ));
        }
        Ok(Self {
            ordinary,
            extraordinary,
            ordinary_response,
            extraordinary_response,
            wavelength_bounds_nm,
        })
    }

    /// Original explicitly authored ordinary wavelength/key roles and receipts.
    pub const fn ordinary_source(self) -> &'a SpectralComplexIndexStatePoint {
        self.ordinary
    }

    /// Original explicitly authored extraordinary wavelength/key roles and receipts.
    pub const fn extraordinary_source(self) -> &'a SpectralComplexIndexStatePoint {
        self.extraordinary
    }

    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.ordinary.resolved()
    }

    pub const fn permeability_h_m(self) -> f64 {
        self.ordinary_response.permeability_h_m()
    }

    pub const fn maximum_knot_span_m(self) -> f64 {
        self.ordinary_response.maximum_knot_span_m()
    }

    /// Closed intersection of both original vacuum-wavelength supports [nm].
    pub const fn wavelength_bounds_nm(self) -> [f64; 2] {
        self.wavelength_bounds_nm
    }

    /// Original positive principal phase indices [ordinary, extraordinary].
    /// Both owners admit frequency before any endpoint roundoff correction;
    /// positive loss in either axis refuses the transparent-crystal model.
    pub fn phase_indices(self, omega_rad_s: f64) -> Result<[f64; 2], ReflectionError> {
        let [ordinary, extraordinary] = self.complex_principal_indices(omega_rad_s)?;
        if ordinary.im != 0.0 || extraordinary.im != 0.0 {
            return Err(ReflectionError::Invalid(
                "transparent uniaxial refraction cannot discard source extinction in either axis",
            ));
        }
        Ok([positive(ordinary.re)?, positive(extraordinary.re)?])
    }

    /// Original ordered [ordinary, extraordinary] principal indices n-i*k.
    /// Both original wavelength supports, knot-gap limits and state bindings
    /// remain authoritative. A consumer admitting loss must retain both
    /// complex values in Maxwell's tensor; taking only their real parts is
    /// not an absorbing-crystal model.
    pub fn complex_principal_indices(self, omega_rad_s: f64) -> Result<[C64; 2], ReflectionError> {
        Ok([
            self.ordinary_response.complex_index(omega_rad_s)?,
            self.extraordinary_response.complex_index(omega_rad_s)?,
        ])
    }

    /// Original principal group indices [ordinary, extraordinary], n+omega*dn/domega.
    /// Both original frequency jets must remain lossless. Their wavelength
    /// coordinates, interpolation gaps, endpoint and cusp refusals are retained.
    /// Zero and negative values are kept: a local slope is not a signal-speed
    /// or causal energy-storage model. This query allocates no storage.
    pub fn group_indices(self, omega_rad_s: f64) -> Result<[f64; 2], ReflectionError> {
        let jets = [
            self.ordinary_response.frequency_jet(omega_rad_s)?,
            self.extraordinary_response.frequency_jet(omega_rad_s)?,
        ];
        if jets
            .iter()
            .any(|jet| jet.epsilon.im != 0.0 || jet.first_s_m.im != 0.0)
        {
            return Err(ReflectionError::Invalid(
                "transparent uniaxial group transport requires lossless original frequency jets",
            ));
        }
        let indices = jets.map(|jet| jet.first_s_m.re * super::VACUUM_SPEED_M_S);
        if indices.iter().any(|value| !value.is_finite()) {
            return Err(ReflectionError::Invalid(
                "uniaxial principal group indices are unrepresentable",
            ));
        }
        Ok(indices)
    }

    // The normal-wavevector owner needs both original derivative orders.
    // Keep this internal: public results retain the complete crystal/interface
    // source, rather than exposing detachable constitutive derivative tables.
    pub(in crate::plane_wave) fn transparent_frequency_jets(
        self,
        omega_rad_s: f64,
    ) -> Result<[super::frequency::OpticalFrequencyJet; 2], ReflectionError> {
        let jets = [
            self.ordinary_response.frequency_jet(omega_rad_s)?,
            self.extraordinary_response.frequency_jet(omega_rad_s)?,
        ];
        if jets.iter().any(|jet| {
            jet.epsilon.im != 0.0
                || jet.epsilon_derivatives.iter().any(|value| value.im != 0.0)
                || jet.first_s_m.im != 0.0
                || jet.second_scaled_s_m.im != 0.0
        }) {
            return Err(ReflectionError::Invalid(
                "transparent crystal dispersion requires lossless original second-order frequency jets",
            ));
        }
        Ok(jets)
    }

    /// Original principal SI permittivities N²/(mu*c0²) [F/m]. The shared
    /// scalar mu belongs to the same complete source snapshot. Negative real
    /// epsilon is admitted; positive extinction may never round to zero loss.
    pub fn principal_permittivities_f_m(
        self,
        omega_rad_s: f64,
    ) -> Result<[C64; 2], ReflectionError> {
        Ok([
            self.ordinary_response.permittivity_f_m(omega_rad_s)?,
            self.extraordinary_response.permittivity_f_m(omega_rad_s)?,
        ])
    }
}

#[cfg(test)]
#[path = "uniaxial_tests.rs"]
mod tests;
