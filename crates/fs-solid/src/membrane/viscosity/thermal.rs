//! Bind thermal endpoint viscosity through the original pinned source claims.
//!
//! The mechanical patch already selected its endpoint elastic model and fixed
//! natural reference. This module supplies no temperature evolution, thermal
//! coefficient fit, hot rebind, or mechanical parameter-work transition.

use super::{
    MaterialPlane, MembraneLaw, MembranePatch, MembraneViscosity, MembraneViscosityError,
    PlaneStressViscosity, Result, invalid, poll,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};

impl<'a> MembraneViscosity<'a> {
    /// Bind the SAME viscosity claims at each selected thermal elastic endpoint.
    ///
    /// Each input law must still belong to the sheet's EXACT complete reference
    /// bundle, in Axes12 warp/weft coordinates. Each card must be the original
    /// immutable card, in triangle order; references may share a card without
    /// cloning it. Nonthermal facets retain their original loss laws unchanged.
    ///
    /// Thermal facets refresh only the six selected viscosity properties using
    /// the original requirements and ClaimIds. The material state owner admits
    /// the endpoint query at the hyperelastic model's exact context, preserving
    /// every axis name and quantity convention. Missing source support refuses;
    /// another claim cannot substitute. The resulting tensor must still be
    /// passive and cover the complete declared frequency band.
    ///
    /// Reference mass, thickness, caloric receipts and natural G are NEVER
    /// re-evaluated. The original natural-volume pullback remains in force and
    /// [`Self::reference_laws`] retains all reference receipts. [`Self::laws`] holds
    /// the actual endpoint loss coefficients used by force, heat and projection.
    ///
    /// This is one fixed endpoint on the SAME immutable mechanical patch. It
    /// does not qualify an intervening temperature path, update an existing
    /// history, account parameter work, or evolve coefficients as heat is added.
    /// Changing the elastic endpoint/reference requires its mechanical owner.
    pub fn new_with_thermal_endpoints(
        patch: &'a MembranePatch,
        reference_laws: Vec<PlaneStressViscosity>,
        cards: &[&MaterialCard],
        resolved_band_rad_s: (f64, f64),
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let (lo, hi) = resolved_band_rad_s;
        if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi < lo {
            return Err(invalid("finite ordered membrane viscosity band required"));
        }
        if reference_laws.len() != patch.elements.len() || cards.len() != reference_laws.len() {
            return Err(invalid(
                "one original viscosity law and source card per triangle required",
            ));
        }
        let mut laws = Vec::with_capacity(reference_laws.len());
        for (index, ((element, sheet), reference)) in patch
            .elements
            .iter()
            .zip(&patch.sheets)
            .zip(&reference_laws)
            .enumerate()
        {
            poll(&mut cancelled)?;
            if reference.source() != sheet.material() || reference.plane() != MaterialPlane::Axes12
            {
                return Err(invalid(
                    "thermal viscosity needs the exact reference sheet bundle and axes",
                ));
            }
            let card = cards[index];
            if card.content_hash() != reference.source().card_identity() {
                return Err(invalid(
                    "thermal viscosity needs the original immutable material card",
                ));
            }
            let law = match &element.law {
                MembraneLaw::Hyperelastic(elastic)
                    if elastic.mean_thermal_expansion().is_some() =>
                {
                    endpoint_law(
                        reference,
                        card,
                        elastic.model().context(),
                        index,
                        &mut cancelled,
                    )?
                }
                _ => reference.clone(),
            };
            laws.push(law);
        }
        poll(&mut cancelled)?;
        Self::assemble(
            patch,
            laws,
            Some(reference_laws),
            resolved_band_rad_s,
            cancelled,
        )
    }
}

fn endpoint_law(
    reference: &PlaneStressViscosity,
    card: &MaterialCard,
    endpoint: &QueryPoint,
    element: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PlaneStressViscosity> {
    let mut requirements = Vec::with_capacity(6);
    let mut pins = Vec::with_capacity(6);
    for name in reference.properties() {
        let property = reference
            .source()
            .property(name)
            .ok_or_else(|| invalid("reference viscosity lost an original property"))?;
        requirements.push(property.requirement().clone());
        pins.push((name.clone(), property.answer().receipt.selected));
    }
    let source_error = |source| MembraneViscosityError::Source { element, source };
    // Select the original viscosity subset before refreshing: reference mass
    // and thickness may legitimately have no support at the endpoint.
    poll(cancelled)?;
    let original = resolve_material_state_point(
        card,
        &reference.source().query(),
        &requirements,
        MaterialPropertySelection::PinnedByProperty(pins),
    )
    .map_err(source_error)?;
    poll(cancelled)?;
    let actual = original
        .at_query_point(card, endpoint)
        .map_err(source_error)?;
    if actual.query() != *endpoint {
        return Err(invalid(
            "viscosity query differs from the exact elastic endpoint context",
        ));
    }
    poll(cancelled)?;
    PlaneStressViscosity::from_material(
        actual,
        reference.plane(),
        reference.properties().each_ref().map(String::as_str),
    )
    .map_err(|source| MembraneViscosityError::Material { element, source })
}
