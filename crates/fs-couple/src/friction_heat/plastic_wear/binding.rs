//! Original caloric and normal material source admission shared by wear owners.

use fs_contact::interface_binding::BoundNormalInterfaceState;
use fs_exec::Cx;
use fs_matdb::{MaterialCard, PropertyKey};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::Density;

use super::{PlasticSlidingHeatError, checkpoint};

/// Retain each original ordered card, complete fixed context and density claim.
/// The existing thermal owner separately checks its registered caloric object
/// and accepted temperature; no chart, query or material history is rebuilt.
pub(super) fn admit_caloric_sources<E>(
    normal: &BoundNormalInterfaceState,
    material_cards: [&MaterialCard; 2],
    caloric: [&ResolvedHeatCapacityCurve; 2],
    cx: &Cx<'_>,
) -> Result<(), PlasticSlidingHeatError<E>> {
    let elastic = normal.elastic();
    let bulk = [elastic.surface_a_state(), elastic.surface_b_state()];
    for i in 0..2 {
        let resolved = bulk[i].resolved();
        let density_name = bulk[i]
            .elastic_state()
            .derivation()
            .map_or(DENSITY_PROPERTY, |d| d.density_requirement().name());
        let density = resolved
            .property_by_key(&PropertyKey::new(density_name, Density::DIMS))
            .ok_or(PlasticSlidingHeatError::Invalid(
                "original normal density source is missing",
            ))?;
        if material_cards[i].content_hash() != resolved.card_identity()
            || caloric[i].curve().material_card_identity() != resolved.card_identity()
            || caloric[i].receipts().is_empty()
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "heat and hardness must retain the original ordered material cards and sourced calorics",
            ));
        }
        for (index, receipt) in caloric[i].receipts().iter().enumerate() {
            if index % 64 == 0 {
                checkpoint(cx)?;
            }
            if receipt.axis_quantities != *resolved.axis_quantities()
                || receipt
                    .query_point
                    .iter()
                    .filter(|(axis, _)| axis != "T")
                    .ne(resolved
                        .query_point()
                        .iter()
                        .filter(|(axis, _)| axis != "T"))
                || (index % 2 == 1
                    && (receipt.property != density_name
                        || receipt.selected != density.answer().receipt.selected))
            {
                return Err(PlasticSlidingHeatError::Invalid(
                    "caloric and normal sources require identical fixed conditions, quantity axes and density claims",
                ));
            }
        }
    }
    Ok(())
}
