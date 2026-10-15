//! Temperature samples lowered through the original material sources.
//!
//! This is the bridge from thermal solver output to existing free-expansion
//! mechanics. It restores the full reference query (including quantity kinds),
//! changes only absolute T, and delegates both alpha integration and mass
//! transport. It neither copies a thermal-expansion law nor replaces matter.

use fs_ga::Vec3;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_material::state_point::{MaterialPropertySelection, integrate_isotropic_thermal_expansion};

use super::{ExpandedBody, ThermalExpansionBinding, expand_body};
use crate::{MaterialInertiaError, ResolvedBody, poll};

/// One solver-supplied uniform region temperature and its original source.
#[derive(Clone, Debug)]
pub struct RegionTemperature<'a> {
    /// Exact reference region key; all regions must be supplied once.
    pub key: &'a str,
    /// The same immutable material card used to construct the reference body.
    pub card: &'a MaterialCard,
    /// Absolute temperature of the configured region [K].
    pub temperature_k: f64,
    /// Expansion-source selection at the reference state. The existing material
    /// integrator pins that selected claim throughout the temperature path.
    pub selection: MaterialPropertySelection,
    /// Prescribed free-expansion fixed point in the reference body frame [m].
    pub fixed_point_m: Vec3,
}

/// Reconfigure an existing body directly from material-bound temperature data.
///
/// Every non-temperature coordinate and its quantity descriptor is restored
/// from the original density receipt. A temperature sample cannot silently
/// change pressure, composition, material lot or the reference geometry. The
/// original material owner integrates alpha(T), and `expand_body` conserves
/// reference mass while updating geometry and the complete inertia tensor.
///
/// Samples can be permuted but must cover the whole body. Missing T at the
/// reference state, a foreign card, unsupported source paths, nonpositive
/// temperature or late cancellation return no configuration. The source card
/// is immutable and selection is always anchored at the original reference T,
/// so unchanged selection inputs retain the same claim and never compound free strain.
/// This retains the free-expansion (not constrained thermoelastic) scope of
/// `expand_body`. Uniform temperature per region is an explicit approximation.
pub fn expand_body_at_temperatures<'a>(
    reference: &'a ResolvedBody,
    temperatures: &[RegionTemperature<'_>],
    mut cancelled: impl FnMut() -> bool,
) -> Result<ExpandedBody<'a>, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if temperatures.len() != reference.regions().len() {
        return Err(MaterialInertiaError::Invalid("temperature samples must cover every reference region"));
    }
    let mut ordered: Vec<_> = temperatures.iter().collect();
    ordered.sort_by_key(|sample| sample.key);
    // Admit identities before performing any constitutive integrations.
    for (region, sample) in reference.regions().iter().zip(&ordered) {
        poll(&mut cancelled)?;
        if region.key() != sample.key {
            return Err(MaterialInertiaError::Invalid("temperature region keys are missing, repeated or foreign"));
        }
        if sample.card.content_hash() != region.material().card_identity()
            || sample.card.id() != region.material().material() {
            return Err(MaterialInertiaError::Invalid("temperature source must be the original reference material card"));
        }
        if !sample.temperature_k.is_finite() || sample.temperature_k <= 0.0 {
            return Err(MaterialInertiaError::Invalid("region temperature must be finite positive absolute kelvin"));
        }
    }
    let mut paths = Vec::with_capacity(ordered.len());
    for (region, sample) in reference.regions().iter().zip(&ordered) {
        poll(&mut cancelled)?;
        let mut point = QueryPoint::new();
        for (axis, value) in region.material().query_point() {
            point = match region.material().axis_quantities().get(axis) {
                Some(&quantity) => point.with_quantity(axis.clone(), quantity, *value),
                None => point.with(axis.clone(), *value),
            }.map_err(|_| MaterialInertiaError::Invalid("stored reference query could not be restored"))?;
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        if !point.axes().contains_key("T") || point.axis_quantities().get("T").is_some_and(|q| {
            *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)
        }) {
            return Err(MaterialInertiaError::Invalid("reference T must be an absolute-kelvin coordinate"));
        }
        let current = point.clone().with("T", sample.temperature_k)
            .map_err(|_| MaterialInertiaError::Invalid("temperature is incompatible with the reference quantity schema"))?;
        paths.push(integrate_isotropic_thermal_expansion(
            sample.card, &point, &current, sample.selection.clone(),
        )?);
    }
    let bindings: Vec<_> = ordered.iter().zip(&paths).map(|(sample, path)| {
        ThermalExpansionBinding { key: sample.key, expansion: path, fixed_point_m: sample.fixed_point_m }
    }).collect();
    expand_body(reference, &bindings, &mut cancelled)
}
