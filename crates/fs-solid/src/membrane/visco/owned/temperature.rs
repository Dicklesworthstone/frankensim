//! Source-prepared temperature endpoints on the retained spatial membrane.
//!
//! Each original triangle evaluates its own retained engineering-length path,
//! datum and material frame. Only natural strain and held source context change;
//! geometry, mass and all accepted memory continue through the native owned
//! reference transaction. Preparation never globally scales a heated mesh.

use std::sync::Arc;

use fs_matdb::MaterialCard;
use fs_material::sheet::{MembraneLaw, SheetError};

use super::{OwnedFiniteViscoMembraneMotion, OwnedFiniteViscoMembraneMotionReferenceChange};
use crate::membrane::visco::{
    FiniteViscoMembraneError, FiniteViscoStepControls, ViscoelasticLawField, invalid,
    material_error, poll,
};

impl OwnedFiniteViscoMembraneMotion {
    /// Propose one source-supported absolute temperature per original triangle.
    ///
    /// Every current sheet must already carry a reusable `PreparedFreeExpansion`
    /// reference selected through `with_hyperelastic_free_expansion`. Each card
    /// must be that triangle's original card; both its retained expansion path
    /// and selected model must admit the temperature. Regional values preserve
    /// triangle order and are never averaged or interpreted as nodal samples.
    ///
    /// Expansion always uses the original datum and fixed source/frame mapping,
    /// including on cooling or temperature revisits. The current equilibrium,
    /// branch energies, kinetics and local controls are retained; only their G
    /// is replaced by the admitted sheet reference. A target requiring different
    /// energy coefficients is refused by the original native binding checks.
    ///
    /// Reference coordinates, prepared gradients, connectivity, grain, corner
    /// inventory, cached mass and strain envelope are copied without rebuilding
    /// them. The native reference transaction then restores transverse plane
    /// stress with retained viscous metrics and unchanged displacement/velocity.
    /// Subsequent native advancement resolves the changed global forces.
    ///
    /// Returned reference work and pressure-volume actuation are prescribed
    /// inputs, separate from numerical defects and physical relaxation heat.
    /// This method advances neither time nor a thermal/entropy balance. Any
    /// source, work or cancellation refusal leaves the accepted owner intact.
    pub fn with_temperatures(
        &self,
        cards: &[&MaterialCard],
        temperatures_k: &[f64],
        controls: FiniteViscoStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OwnedFiniteViscoMembraneMotionReferenceChange, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        let count = self.patch.sheets.len();
        if cards.len() != count || temperatures_k.len() != count {
            return Err(invalid(
                "one original source card and absolute temperature per membrane triangle required",
            ));
        }
        let mut target_patch = (*self.patch).clone();
        let mut target_laws = Vec::with_capacity(count);
        let current = self.native_laws();
        for (index, (card, temperature)) in cards.iter().zip(temperatures_k).enumerate() {
            poll(&mut cancelled)?;
            let sheet = self.patch.sheets[index]
                .at_hyperelastic_temperature(card, *temperature)
                .map_err(|source| material_error(index, source))?;
            let MembraneLaw::Hyperelastic(elastic) = sheet
                .membrane_law()
                .map_err(|source| material_error(index, source))?
            else {
                return Err(material_error(
                    index,
                    SheetError::Invalid("source-prepared hyperelastic membrane required"),
                ));
            };
            let reference = elastic.free_expansion_reference().ok_or_else(|| {
                material_error(
                    index,
                    SheetError::Invalid("retained engineering-length source path required"),
                )
            })?;
            target_laws.push(
                current
                    .at(index)
                    .clone()
                    .with_reference_deformation(*reference.reference_deformation()),
            );
            target_patch.sheets[index] = sheet;
            target_patch.elements[index].law = MembraneLaw::Hyperelastic(elastic);
        }
        self.change_reference(
            Arc::new(target_patch),
            ViscoelasticLawField::PerElement(&target_laws),
            controls,
            &mut cancelled,
        )
    }
}
