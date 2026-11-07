//! Endpoint tensile activation retaining the exact evaluated source bundle.

use super::{
    ResolvedBrittleCrackBand, ResolvedTensileInitiationError, ResolvedTensileInitiator,
    resolve_plane_band,
};
use crate::Voigt;
use crate::brittle::initiation::TensileActivation;

/// An endpoint observation under one complete original resolved source.
/// It carries no loading path, event time or permission to seed current damage.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTensileActivation<'a> {
    source: &'a ResolvedBrittleCrackBand,
    native: TensileActivation<'a>,
}

impl<'a> ResolvedTensileInitiator<'a> {
    /// Observe this source's elastic stress and strength at the supplied strain.
    /// Thermal consumers supply MECHANICAL strain. The source's dormant plane
    /// does not substitute for a required physical principal-direction hint.
    pub fn activation_at(
        &self,
        strain: &Voigt,
        direction_hint: Option<[f64; 3]>,
    ) -> Result<Option<ResolvedTensileActivation<'a>>, ResolvedTensileInitiationError> {
        self.activation_at_controlled(strain, direction_hint, &mut || false)
    }

    /// Preserve original typed admission and cancellation through the native
    /// bounded endpoint observation, without interpolating source coefficients.
    pub fn activation_at_controlled<C: FnMut() -> bool>(
        &self,
        strain: &Voigt,
        direction_hint: Option<[f64; 3]>,
        cancelled: &mut C,
    ) -> Result<Option<ResolvedTensileActivation<'a>>, ResolvedTensileInitiationError> {
        Ok(self
            .native
            .activation_at_controlled(strain, direction_hint, cancelled)?
            .map(|native| ResolvedTensileActivation {
                source: self.source,
                native,
            }))
    }
}

impl<'a> ResolvedTensileActivation<'a> {
    #[must_use]
    pub const fn source(&self) -> &'a ResolvedBrittleCrackBand {
        self.source
    }
    #[must_use]
    pub const fn native(&self) -> &TensileActivation<'a> {
        &self.native
    }

    /// Prepare the observed plane at a width supplied by receiving geometry.
    /// All source claims, coefficients and controls remain identical. This
    /// creates no accepted state: old-source history transfer and new-source
    /// replay belong to the consumer, not this already-loaded observation.
    pub fn resolve_band(
        &self,
        band_width_m: f64,
    ) -> Result<ResolvedBrittleCrackBand, ResolvedTensileInitiationError> {
        resolve_plane_band(self.source, self.native.normal(), band_width_m)
    }
}
