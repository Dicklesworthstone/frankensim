//! Source-prepared natural lengths for the existing finite-strain membrane.
//!
//! Every endpoint evaluates the ORIGINAL engineering-length path and datum.
//! Mass, thickness and their query receipts remain in the original sheet;
//! only its natural reference and held model context change. Preparation owns
//! all source work; the existing mechanical response performs no source query.
//! This is a prescribed thermal reference, not entropy, caloric feedback,
//! relaxation kinetics or permission to overwrite an accepted material history.

use std::sync::Arc;

use super::{HyperelasticMembrane, NaturalReference, ResolvedHyperelastic};
use crate::hyper::reference::ReferenceDeformation;
use crate::sheet::{ResolvedSheet, SheetError};
use crate::state_point::is_fixed_external_coordinate;
use crate::thermal::transport::expansion::{
    FreeExpansionBasis, FreeExpansionGeometrySample, PreparedFreeExpansion,
};
use fs_blake3::ContentHash;
use fs_matdb::{MaterialCard, QueryPoint};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Declared alignment with the fixed warp/weft/reference-normal sheet axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SheetExpansionFrame {
    /// Requires a path explicitly prepared as isotropic.
    Isotropic,
    FixedPrincipal {
        /// Exact frame identity retained by the directional source path.
        principal_frame: ContentHash,
        /// Nonzero identity of the sheet's fixed material axes.
        sheet_frame: ContentHash,
        /// Columns are source axes in sheet coordinates. Only proper signed
        /// axis permutations are supported, preserving diagonal G at ALL T.
        /// Arbitrary rotation or shear is refused, never projected away.
        principal_to_sheet: [[f64; 3]; 3],
    },
}

impl SheetExpansionFrame {
    fn reference(
        self,
        source: &PreparedFreeExpansion,
        geometry: FreeExpansionGeometrySample,
    ) -> Result<ReferenceDeformation, SheetError> {
        let q = match (source.basis(), self) {
            (FreeExpansionBasis::Isotropic, Self::Isotropic) => IDENTITY,
            (
                FreeExpansionBasis::FixedPrincipal { frame },
                Self::FixedPrincipal {
                    principal_frame,
                    sheet_frame,
                    principal_to_sheet,
                },
            ) if frame == principal_frame && sheet_frame != ContentHash([0; 32]) => {
                if principal_to_sheet
                    .iter()
                    .flatten()
                    .any(|&x| x != 0.0 && x.abs() != 1.0)
                    || (principal_frame == sheet_frame && principal_to_sheet != IDENTITY)
                {
                    return Err(SheetError::Invalid(
                        "sheet expansion requires aligned principal axes; thermal shear and inconsistent frame maps are unsupported",
                    ));
                }
                principal_to_sheet
            }
            _ => {
                return Err(SheetError::Invalid(
                    "sheet expansion requires its original source basis and named aligned material frames",
                ));
            }
        };
        // The native reference owner verifies positive stretches, orthonormal
        // proper axes and representability. Signed permutations keep its
        // diagonal image exact, including when directional lengths coincide.
        Ok(ReferenceDeformation::from_principal_stretches(
            geometry.principal_stretches(),
            q,
        )?)
    }
}

/// Immutable thermal endpoint retaining the complete prepared source path.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetExpansionReference {
    source: Arc<PreparedFreeExpansion>,
    geometry: FreeExpansionGeometrySample,
    frame: SheetExpansionFrame,
    reference_context: QueryPoint,
    current_context: QueryPoint,
    reference: ReferenceDeformation,
}

fn context_at(source: &PreparedFreeExpansion, temperature: f64) -> Result<QueryPoint, SheetError> {
    let axis = source.temperature_axis();
    let original = source.reference_query_point();
    if axis != "T" && original.axes().contains_key("T") {
        return Err(SheetError::Invalid(
            "ambiguous additional sheet temperature coordinate",
        ));
    }
    let mut context = QueryPoint::new();
    for (name, value) in original.axes() {
        if name == axis {
            continue;
        }
        context = match original.axis_quantities().get(name) {
            Some(quantity) => context.with_quantity(name, *quantity, *value),
            None => context.with(name, *value),
        }
        .map_err(|_| SheetError::Invalid("unrepresentable fixed expansion context"))?;
    }
    context
        .with_quantity("T", original.axis_quantities()[axis], temperature)
        .map_err(|_| SheetError::Invalid("unrepresentable expansion endpoint temperature"))
}

impl SheetExpansionReference {
    fn try_new(
        card: &MaterialCard,
        source: PreparedFreeExpansion,
        temperature_k: f64,
        frame: SheetExpansionFrame,
    ) -> Result<Self, SheetError> {
        if source.material_card_identity() != card.content_hash() {
            return Err(SheetError::Invalid(
                "sheet expansion requires its original material card",
            ));
        }
        let axis = source.temperature_axis();
        let original = source.reference_query_point();
        let admit_axis = |name: &str| -> Result<(), SheetError> {
            if name != axis
                && !is_fixed_external_coordinate(
                    name,
                    original.axis_quantities().get(name).copied(),
                    original.axes().get(name).copied(),
                )
            {
                return Err(SheetError::Invalid(
                    "sheet expansion requires fixed external source conditions, not evolving mechanical or history coordinates",
                ));
            }
            Ok(())
        };
        for name in original.axes().keys() {
            admit_axis(name)?;
        }
        for line in source.source_lines() {
            let claim = card
                .claims()
                .claim(line.knots[0].1.receipt.selected)
                .ok_or(SheetError::Invalid(
                    "sheet expansion lost its original selected source",
                ))?;
            for name in claim.validity.bounds().keys() {
                admit_axis(name)?;
            }
            for (name, _) in claim.value.axes() {
                admit_axis(name)?;
            }
        }
        let geometry = source.at_temperature(temperature_k)?;
        let reference = frame.reference(&source, geometry)?;
        let reference_context = context_at(&source, original.axes()[axis])?;
        let current_context = context_at(&source, temperature_k)?;
        Ok(Self {
            source: Arc::new(source),
            geometry,
            frame,
            reference_context,
            current_context,
            reference,
        })
    }

    fn at_temperature(&self, card: &MaterialCard, temperature_k: f64) -> Result<Self, SheetError> {
        if self.source.material_card_identity() != card.content_hash() {
            return Err(SheetError::Invalid(
                "sheet expansion requires its original material card",
            ));
        }
        let geometry = self.source.at_temperature(temperature_k)?;
        Ok(Self {
            source: Arc::clone(&self.source),
            geometry,
            frame: self.frame,
            reference_context: self.reference_context.clone(),
            current_context: context_at(&self.source, temperature_k)?,
            reference: self.frame.reference(&self.source, geometry)?,
        })
    }

    pub fn source(&self) -> &PreparedFreeExpansion {
        &self.source
    }
    /// Sample in ORIGINAL source principal-axis order, before the declared permutation.
    pub const fn geometry(&self) -> &FreeExpansionGeometrySample {
        &self.geometry
    }
    pub const fn frame(&self) -> SheetExpansionFrame {
        self.frame
    }
    pub const fn reference_context(&self) -> &QueryPoint {
        &self.reference_context
    }
    pub const fn current_context(&self) -> &QueryPoint {
        &self.current_context
    }
    /// Diagonal natural map in the fixed sheet material axes.
    pub const fn reference_deformation(&self) -> &ReferenceDeformation {
        &self.reference
    }
    pub const fn temperature_k(&self) -> f64 {
        self.geometry.temperature_k()
    }

    fn into_natural_reference(self) -> NaturalReference {
        let g = self.reference.deformation_gradient();
        let stretches = [g[0], g[4], g[8]];
        NaturalReference {
            stretches,
            inverse: stretches.map(f64::recip),
            volume_ratio: self.reference.volume_ratio(),
            thermal: None,
            expansion: Some(Box::new(self)),
        }
    }
}

impl HyperelasticMembrane {
    /// Original engineering-length path and admitted current thermal endpoint.
    /// Explicit stretches and isolated mean-expansion sources return `None`.
    pub fn free_expansion_reference(&self) -> Option<&SheetExpansionReference> {
        self.natural_reference.as_ref()?.expansion.as_deref()
    }
}

impl ResolvedSheet {
    /// Select a reusable, source-supported thermal reference exactly once.
    /// The source datum must equal the ORIGINAL mass/thickness/model context.
    /// Directional maps must explicitly align with the sheet axes. Original
    /// inventory, independent properties and the selected model card survive.
    /// Native energy retains `Psi = det(G) W(F G^-1)` per original volume.
    pub fn with_hyperelastic_free_expansion(
        mut self,
        card: &MaterialCard,
        source: PreparedFreeExpansion,
        temperature_k: f64,
        frame: SheetExpansionFrame,
    ) -> Result<Self, SheetError> {
        let previous = self.hyperelastic.as_deref().ok_or(SheetError::Invalid(
            "select a hyperelastic membrane before its expansion path",
        ))?;
        if previous.natural_reference.is_some()
            || self.material.card_identity() != card.content_hash()
        {
            return Err(SheetError::Invalid(
                "sheet expansion requires an unmodified same-card natural reference",
            ));
        }
        let endpoint = SheetExpansionReference::try_new(card, source, temperature_k, frame)?;
        if self.material.query() != endpoint.reference_context
            || *previous.model.context() != endpoint.reference_context
        {
            return Err(SheetError::Invalid(
                "sheet inventory and original model must belong to the expansion datum context",
            ));
        }
        let model = ResolvedHyperelastic::try_from_model_card_at(
            previous.model.model_card(),
            &endpoint.current_context,
        )?;
        let mut law = previous.clone();
        law.model = model;
        law.natural_reference = Some(endpoint.into_natural_reference());
        self.hyperelastic = Some(Arc::new(law));
        Ok(self)
    }

    /// Prepare a new held-temperature endpoint from the SAME original path,
    /// datum, fixed frame and selected model. No expansion ratios compound and
    /// no mass/thickness source is refreshed. Both source and model domains
    /// must admit the temperature. The old sheet remains unchanged on refusal.
    ///
    /// This does not update an accepted viscoelastic history or supply entropy
    /// exchange. Use the native reference transaction, with its signed work.
    /// That transaction requires unchanged equilibrium and branch energies;
    /// a source-prepared endpoint is not permission for modulus replacement.
    pub fn at_hyperelastic_temperature(
        &self,
        card: &MaterialCard,
        temperature_k: f64,
    ) -> Result<Self, SheetError> {
        let previous = self.hyperelastic.as_deref().ok_or(SheetError::Invalid(
            "a hyperelastic expansion path is required",
        ))?;
        let endpoint = previous.free_expansion_reference()
            .ok_or(SheetError::Invalid("a reusable engineering-length path is required; endpoint-only evidence cannot evolve temperature"))?
            .at_temperature(card, temperature_k)?;
        let model = ResolvedHyperelastic::try_from_model_card_at(
            previous.model.model_card(),
            &endpoint.current_context,
        )?;
        let mut law = previous.clone();
        law.model = model;
        law.natural_reference = Some(endpoint.into_natural_reference());
        let mut result = self.clone();
        result.hyperelastic = Some(Arc::new(law));
        Ok(result)
    }
}
