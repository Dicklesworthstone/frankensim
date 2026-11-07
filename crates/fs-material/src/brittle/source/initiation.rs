//! Source-preserving tensile initiation and elastic-to-crack-band handoff.
//!
//! The original resolved material supplies the elastic backbone and tensile
//! strength. The native detector locates the plane and onset on a previously
//! elastic affine strain path. Resolving that plane for receiving geometry
//! retains the complete original source bundle, every coefficient and the
//! original strain/work bounds; no new card query or parameter fit occurs.
//! Endpoint activation retains the same source for an explicitly evaluated
//! loaded strain; that observation supplies neither an event time nor a seed.
//!
//! This initializer never accepts an existing damaged history. After handoff,
//! callers replay the actual increment through the receiving native law. Its
//! virtual preload represents existing elastic storage, not additional work.
//! The native tensile-only criterion and fixed-plane fracture limitations apply.

mod activation;
pub use activation::ResolvedTensileActivation;

use super::{BrittleSourceError, ResolvedBrittleCrackBand};
use crate::Voigt;
use crate::brittle::initiation::{
    ElasticBrittleHandoff, TensileInitiationError, TensileInitiationSettings, TensileInitiator,
    TensileOnset,
};

/// A native tensile detector retaining the original resolved material owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTensileInitiator<'a> {
    source: &'a ResolvedBrittleCrackBand,
    native: TensileInitiator<'a>,
}

/// A located event retaining its complete original source and native result.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTensileOnset<'a> {
    source: &'a ResolvedBrittleCrackBand,
    native: TensileOnset<'a>,
}

/// Source admission, native initiation, or source-preserving handoff refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ResolvedTensileInitiationError {
    Source(BrittleSourceError),
    Native(TensileInitiationError),
    /// Card, state coordinates, selected claims, receipts or coefficients differ.
    SourceMismatch,
    /// The receiver changed the original strain domain or material-work gate.
    ControlsMismatch,
}

impl core::fmt::Display for ResolvedTensileInitiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound tensile initiation refused: {self:?}")
    }
}
impl std::error::Error for ResolvedTensileInitiationError {}
impl From<BrittleSourceError> for ResolvedTensileInitiationError {
    fn from(error: BrittleSourceError) -> Self {
        Self::Source(error)
    }
}
impl From<TensileInitiationError> for ResolvedTensileInitiationError {
    fn from(error: TensileInitiationError) -> Self {
        Self::Native(error)
    }
}

impl ResolvedBrittleCrackBand {
    /// Discover tensile onset from this exact admitted elastic/strength source.
    ///
    /// The source's existing plane and width are not an onset-direction hint:
    /// the native detector selects a physical principal plane from the stress.
    /// Only numerical event resolution is supplied separately.
    pub fn tensile_initiator(
        &self,
        settings: TensileInitiationSettings,
    ) -> Result<ResolvedTensileInitiator<'_>, ResolvedTensileInitiationError> {
        let native = TensileInitiator::new(
            self.law().elastic(),
            self.law().parameters().tensile_strength_pa,
            settings,
        )?;
        Ok(ResolvedTensileInitiator {
            source: self,
            native,
        })
    }
}

impl<'a> ResolvedTensileInitiator<'a> {
    /// Original material bundle and original prepared constitutive owner.
    #[must_use]
    pub const fn source(&self) -> &'a ResolvedBrittleCrackBand {
        self.source
    }

    /// Original detector, including its admitted numerical resolution.
    #[must_use]
    pub const fn native(&self) -> &TensileInitiator<'a> {
        &self.native
    }

    /// Locate the first tensile crossing on an affine, previously elastic path.
    /// A repeated principal stress requires an explicit physical direction hint.
    pub fn between(
        &self,
        start: &Voigt,
        end: &Voigt,
        direction_hint: Option<[f64; 3]>,
    ) -> Result<Option<ResolvedTensileOnset<'a>>, ResolvedTensileInitiationError> {
        self.between_controlled(start, end, direction_hint, &mut || false)
    }

    /// Poll cancellation through the original bounded native onset detector.
    /// Batch consumers remain responsible for their count and allocation bounds.
    pub fn between_controlled<C: FnMut() -> bool>(
        &self,
        start: &Voigt,
        end: &Voigt,
        direction_hint: Option<[f64; 3]>,
        cancelled: &mut C,
    ) -> Result<Option<ResolvedTensileOnset<'a>>, ResolvedTensileInitiationError> {
        Ok(self
            .native
            .between_controlled(start, end, direction_hint, cancelled)?
            .map(|native| ResolvedTensileOnset {
                source: self.source,
                native,
            }))
    }
}

impl<'a> ResolvedTensileOnset<'a> {
    /// Original source, not a reconstructed material with matching numbers.
    #[must_use]
    pub const fn source(&self) -> &'a ResolvedBrittleCrackBand {
        self.source
    }

    /// Original numerical event: bracket, plane, strains and spectral residual.
    #[must_use]
    pub const fn native(&self) -> &TensileOnset<'a> {
        &self.native
    }

    /// Prepare the discovered plane at a receiving geometry's declared width.
    ///
    /// The width must come from the receiving geometric owner. The material
    /// checks both fracture envelopes, but cannot establish geometric V/A here.
    /// All other original settings, coefficients, source claims and additional
    /// multiphysics receipts are preserved. No accepted history is transferred.
    pub fn resolve_band(
        &self,
        band_width_m: f64,
    ) -> Result<ResolvedBrittleCrackBand, ResolvedTensileInitiationError> {
        resolve_plane_band(self.source, self.native.normal(), band_width_m)
    }

    /// Seed the exact receiving source/law at the original elastic path start.
    ///
    /// Identical E, nu and tensile strength alone are insufficient: original
    /// card, state, selected claims, extra receipts and all fracture parameters
    /// must agree. The receiver may change only the selected plane and its
    /// geometric width. Native seeding checks plane, domains, undamaged state
    /// and stress/storage continuity. It never charges virtual preload work.
    pub fn seed<'b>(
        &self,
        receiving: &'b ResolvedBrittleCrackBand,
    ) -> Result<ElasticBrittleHandoff<'b>, ResolvedTensileInitiationError> {
        if receiving.material() != self.source.material()
            || receiving.law().parameters() != self.source.law().parameters()
        {
            return Err(ResolvedTensileInitiationError::SourceMismatch);
        }
        let original = self.source.law().settings();
        let target = receiving.law().settings();
        if target.maximum_strain != original.maximum_strain
            || target.energy_tolerance_j_m3 != original.energy_tolerance_j_m3
        {
            return Err(ResolvedTensileInitiationError::ControlsMismatch);
        }
        Ok(self.native.seed(receiving.law())?)
    }
}

fn resolve_plane_band(
    source: &ResolvedBrittleCrackBand,
    normal: [f64; 3],
    band_width_m: f64,
) -> Result<ResolvedBrittleCrackBand, ResolvedTensileInitiationError> {
    let mut settings = source.law().settings();
    settings.normal = normal;
    settings.band_width_m = band_width_m;
    Ok(ResolvedBrittleCrackBand::try_from_material(
        source.material().clone(),
        settings,
    )?)
}
