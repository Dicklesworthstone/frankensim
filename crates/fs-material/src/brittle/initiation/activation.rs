//! Principal tensile activation at an explicitly supplied elastic strain.
//!
//! Changing source coefficients can make an accepted strain critical without
//! supplying a new strain path. This observation evaluates that endpoint only;
//! it invents neither a zero-to-current loading path nor an onset time. The
//! original affine detector continues to refuse a nonelastic path start.

use super::{
    TensileInitiationError, TensileInitiator, normalize, poll, principal_plane, spectrum,
    threshold_side,
};
use crate::Voigt;
use crate::elastic::IsotropicElastic;

/// Original elastic stress at or above the declared tensile criterion.
///
/// This is a numerical endpoint observation, not a located crossing, accepted
/// damage state or interval certificate. It cannot seed an undamaged history
/// at this already-critical strain. A consumer must preserve its OLD accepted
/// source/state and separately transfer and replay through the native band.
#[derive(Clone, Debug, PartialEq)]
pub struct TensileActivation<'a> {
    elastic: &'a IsotropicElastic,
    strength_pa: f64,
    strain: Voigt,
    stress_pa: Voigt,
    normal: [f64; 3],
    principal_stresses_pa: [f64; 3],
    multiplicity: usize,
    eigen_residual_pa: f64,
}

impl<'a> TensileInitiator<'a> {
    /// Observe the strict Rankine criterion at this supplied strain: `None`
    /// below strength, otherwise the actual principal tensile plane/stress.
    /// Strain is MECHANICAL for thermal consumers. No temperature, strain-path
    /// fraction, first-crossing time or prior elastic-history claim is inferred.
    pub fn activation_at(
        &self,
        strain: &Voigt,
        direction_hint: Option<[f64; 3]>,
    ) -> Result<Option<TensileActivation<'a>>, TensileInitiationError> {
        self.activation_at_controlled(strain, direction_hint, &mut || false)
    }

    /// Cancellation is polled before/after bounded native spectral work.
    /// The guarded below-threshold screen avoids eigensolver allocation for
    /// ordinary elastic points; a loaded point uses one original 3x3 spectrum.
    /// Repeated principal tension retains the native physical-hint refusal.
    pub fn activation_at_controlled<C: FnMut() -> bool>(
        &self,
        strain: &Voigt,
        direction_hint: Option<[f64; 3]>,
        cancelled: &mut C,
    ) -> Result<Option<TensileActivation<'a>>, TensileInitiationError> {
        poll(cancelled)?;
        let hint = direction_hint.map(normalize).transpose()?;
        self.check_strain(strain)?;
        let stress_pa = self.stress(strain)?;
        let side = threshold_side(&stress_pa, self.strength_pa);
        if side == Some(true) {
            poll(cancelled)?;
            return Ok(None);
        }
        let tolerance = self.strength_pa * self.settings.relative_stress_tolerance;
        let spectral = spectrum(&stress_pa, tolerance, cancelled)?;
        if spectral.values[2] < self.strength_pa {
            if side == Some(false) {
                return Err(TensileInitiationError::UnresolvedSpectrum);
            }
            poll(cancelled)?;
            return Ok(None);
        }
        let (normal, multiplicity) = principal_plane(&spectral, tolerance, hint)?;
        poll(cancelled)?;
        Ok(Some(TensileActivation {
            elastic: self.elastic,
            strength_pa: self.strength_pa,
            strain: *strain,
            stress_pa,
            normal,
            principal_stresses_pa: spectral.values,
            multiplicity,
            eigen_residual_pa: spectral.residual_pa,
        }))
    }
}

impl<'a> TensileActivation<'a> {
    #[must_use]
    pub const fn elastic(&self) -> &'a IsotropicElastic {
        self.elastic
    }
    #[must_use]
    pub const fn tensile_strength_pa(&self) -> f64 {
        self.strength_pa
    }
    /// Exact supplied elastic strain, not a synthesized loading endpoint.
    #[must_use]
    pub const fn strain(&self) -> &Voigt {
        &self.strain
    }
    #[must_use]
    pub const fn stress_pa(&self) -> Voigt {
        self.stress_pa
    }
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }
    /// Ascending principal elastic stresses at the observed endpoint, Pa.
    #[must_use]
    pub const fn principal_stresses_pa(&self) -> [f64; 3] {
        self.principal_stresses_pa
    }
    #[must_use]
    pub const fn principal_multiplicity(&self) -> usize {
        self.multiplicity
    }
    #[must_use]
    pub const fn eigen_residual_pa(&self) -> f64 {
        self.eigen_residual_pa
    }
}
