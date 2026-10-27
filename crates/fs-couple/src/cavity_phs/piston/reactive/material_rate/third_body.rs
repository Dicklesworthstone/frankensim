//! Private concentration arithmetic shared by the declared collider regimes.
//!
//! `None` is exact absence of effective colliders, before any normalization
//! arithmetic. `Some(log_q)` is a positive represented concentration ratio.
//! Source admission and both kinetic laws remain in the original owners.

use super::collider::{ColliderEfficiency, REFERENCE_CONCENTRATION_PROPERTY};
use super::{
    ClosedGasParcel, GasReactionError, ResolvedMaterialStatePoint, invalid, poll, positive,
};
use crate::cavity_phs::piston::reactive::mechanism::Sum;

pub(super) fn normalized_concentration_log(
    efficiencies: &[ColliderEfficiency],
    resolved: &ResolvedMaterialStatePoint,
    gas: &ClosedGasParcel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, GasReactionError> {
    let property = |name| {
        resolved
            .property(name)
            .map(|p| p.value_si())
            .ok_or_else(|| invalid("resolved source omitted a required collider coefficient"))
    };
    let c0 = property(REFERENCE_CONCENTRATION_PROPERTY)?;
    let model = gas.model().thermodynamics();
    let mut weighted = Sum::default();
    for (species, fraction) in model
        .components()
        .iter()
        .zip(model.mole_composition().fractions())
    {
        poll(cancelled)?;
        let index = efficiencies
            .binary_search_by(|b| b.species.cmp(species.species()))
            .map_err(|_| invalid("live gas has no declared collider efficiency"))?;
        let efficiency = property(&efficiencies[index].property)?;
        let contribution = efficiency * fraction;
        if !contribution.is_finite() || (efficiency > 0.0 && *fraction > 0.0 && contribution == 0.0)
        {
            return Err(invalid("positive collider contribution is unrepresentable"));
        }
        weighted.add(contribution)?;
    }
    // Resolve all source claims before this path. Exact absence never attempts
    // a potentially unrepresentable concentration ratio or creates 0*infinity.
    if weighted.value == 0.0 {
        poll(cancelled)?;
        return Ok(None);
    }
    let molar_mass = gas
        .state()
        .thermodynamics()
        .properties()
        .molar_mass()
        .value();
    let concentration = (gas.mass_kg() / molar_mass) / gas.volume_m3();
    if !positive(concentration) || !positive(weighted.value) || !positive(c0) {
        return Err(invalid(
            "positive represented collider concentration required",
        ));
    }
    let log_q =
        fs_math::det::ln(weighted.value) + fs_math::det::ln(concentration) - fs_math::det::ln(c0);
    if !log_q.is_finite() {
        return Err(invalid(
            "collider concentration ratio logarithm is unrepresentable",
        ));
    }
    Ok(Some(log_q))
}
