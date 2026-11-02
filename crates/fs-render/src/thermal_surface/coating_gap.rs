//! Mutual far-field exchange between two opaque, parallel coherent coatings.
//!
//! A vacuum gap with unit view factor returns every reflected ray to the other
//! surface at the same wavelength, cosine and s/p polarization. In the
//! phase-averaged gap limit, summing those reflected POWERS gives the matched-
//! mode transfer fraction `A1*A2 / (1-R1*R2)`. Opaque reciprocal surfaces have
//! `R_i=1-A_i`; positive native work therefore gives the equivalent denominator
//! `A1+A2-A1*A2` without subtracting nearly unit reflectances. The classic
//! parallel-specular-surface limit is Howell and Siegel, NASA SP-164, vol. II
//! (1969), section 4.3.1, equation (4-1); sections 5.3 and 5.6 retain spectral
//! and directional dependence: <https://ntrs.nasa.gov/citations/19690020939>.
//!
//! Each body's films remain optically coherent and internally isothermal;
//! interference between separate gap crossings has been averaged away. This
//! requires a physically justified broad-phase/far-field approximation. The
//! owner does not select that approximation from a supplied distance. It does
//! not model evanescent exchange, gap interference, finite-edge escape, rough
//! reflection, other exterior ports or internal exchange of nonisothermal films.

use super::coating_exchange::{admit, validate_isothermal};
use super::{
    CoatingExchangeBudget, CoatingExchangeError, CoatingExchangeQuadrature, CoatingExchangeRates,
    CoatingExchangeResponse, CoatingThermalEmission, ThermalSurfaceError,
};
use crate::conductor::MAX_COHERENT_STACK_LAYERS;
use fs_blake3::{ContentHash, DomainHasher};

/// Prepared mutual exchange for the two facing sides of opaque planar bodies.
///
/// `area_m2` is their common facing area in the infinite/coextensive-plate
/// approximation: every departing mode reaches the other surface, with no
/// distance attenuation or invented geometric view factor. Both sides have
/// vacuum incidence. Each body has one internal temperature, but the bodies
/// may differ. The original source bands constrain their common integration
/// interval; this response supplies only that finite-band exchange.
///
/// Preparation performs two original coherent absorption solves per quadrature
/// node. Subsequent temperature queries reuse the original Planck exchange
/// owner and allocate nothing. Optical constants remain fixed in the response;
/// native T/h-dependent laws must first be resolved at the current state.
#[derive(Clone, Debug)]
pub struct CoatingGapResponse {
    response: CoatingExchangeResponse,
    optical_identities: [ContentHash; 2],
}

impl CoatingGapResponse {
    /// Prepare both original opaque emitters without assigning either a black
    /// environment. Exact internal isothermality is checked separately on each
    /// side. The budget counts BOTH original optical solves per angular node;
    /// spectral nodes are shared. Failure or cancellation publishes no response.
    pub fn from_emitters(
        first: &CoatingThermalEmission,
        second: &CoatingThermalEmission,
        area_m2: f64,
        quadrature: CoatingExchangeQuadrature,
        budget: CoatingExchangeBudget,
        checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<Self, CoatingExchangeError> {
        let optical_evaluations = admit(quadrature, budget, area_m2, 2)?;
        if first.substrate_temperature_k().is_none() || second.substrate_temperature_k().is_none() {
            return Err(CoatingExchangeError::OpaqueRequired);
        }
        validate_isothermal(first)?;
        validate_isothermal(second)?;
        let first_band = first.wavelength_bounds_nm();
        let second_band = second.wavelength_bounds_nm();
        let wavelength_bounds_nm = (
            first_band.0.max(second_band.0),
            first_band.1.min(second_band.1),
        );
        if wavelength_bounds_nm.0 >= wavelength_bounds_nm.1 {
            return Err(CoatingExchangeError::WavelengthDomain);
        }
        let optical_identities = [first.optical_identity(), second.optical_identity()];
        let mut hash = DomainHasher::new("org.frankensim.render.coating-gap-optics.v1");
        for identity in optical_identities {
            hash.update(identity.as_bytes());
        }
        let response = CoatingExchangeResponse::from_absorptance(
            wavelength_bounds_nm,
            hash.finalize(),
            area_m2,
            quadrature,
            optical_evaluations,
            checkpoint,
            |wavelength, cosine, checkpoint| {
                let a = first.reciprocal_absorbed_fraction(wavelength, cosine)?;
                checkpoint()?;
                let b = second.reciprocal_absorbed_fraction(wavelength, cosine)?;
                Ok([
                    multiple_reflection_exchange(a[0], b[0])?,
                    multiple_reflection_exchange(a[1], b[1])?,
                ])
            },
        )?;
        Ok(Self {
            response,
            optical_identities,
        })
    }

    /// Signed watts from FIRST to SECOND and the positive frozen-optics W/K
    /// secant. A conservative thermal owner debits the first and credits the
    /// second by this same power. Equal temperatures give exact zero heat and
    /// the analytic tangent. Swapping temperatures negates heat with the same
    /// secant; a body temperature query performs no optical solves.
    pub fn evaluate(
        &self,
        first_temperature_k: f64,
        second_temperature_k: f64,
        checkpoint: impl FnMut() -> Result<(), CoatingExchangeError>,
    ) -> Result<CoatingExchangeRates, CoatingExchangeError> {
        self.response
            .evaluate(first_temperature_k, Some(second_temperature_k), checkpoint)
    }

    /// First and second original optical identities, excluding temperature.
    #[must_use]
    pub const fn optical_identities(&self) -> [ContentHash; 2] {
        self.optical_identities
    }

    /// Ordered optical pair and phase-averaged parallel-gap model identity.
    #[must_use]
    pub const fn optical_identity(&self) -> ContentHash {
        self.response.optical_identity()
    }

    /// Pair, common band, physical area, quadrature and retained spectral work.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.response.identity()
    }

    /// Common original source support, in vacuum nanometres.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.response.wavelength_bounds_nm()
    }

    /// Shared physical facing area in square metres.
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.response.area_m2()
    }

    /// Explicit spectral/angular midpoint approximation controls.
    #[must_use]
    pub const fn quadrature(&self) -> CoatingExchangeQuadrature {
        self.response.quadrature()
    }

    /// Original coherent absorption solves, twice the quadrature node count.
    #[must_use]
    pub const fn optical_evaluations(&self) -> usize {
        self.response.optical_evaluations()
    }
}

fn multiple_reflection_exchange(first: f64, second: f64) -> Result<f64, ThermalSurfaceError> {
    // Original positive work may exceed unity by its admitted final rounding
    // error. Refuse a material passivity violation; do not clip weak sources.
    let allowance = 512.0 * f64::EPSILON * (MAX_COHERENT_STACK_LAYERS + 1) as f64;
    if [first, second]
        .into_iter()
        .any(|a| !a.is_finite() || a < 0.0 || a > 1.0 + allowance)
    {
        return Err(ThermalSurfaceError::Unrepresentable);
    }
    let small = first.min(second);
    let large = first.max(second);
    if small == 0.0 {
        return Ok(0.0);
    }
    // Dividing out the larger absorptance keeps even two tiny absorptances
    // from disappearing through their product or a rounded 1-R1*R2 residual.
    let fraction = small / (1.0 + small / large - small);
    if !fraction.is_finite() || fraction <= 0.0 {
        Err(ThermalSurfaceError::Unrepresentable)
    } else {
        Ok(fraction)
    }
}

#[cfg(test)]
mod tests;
