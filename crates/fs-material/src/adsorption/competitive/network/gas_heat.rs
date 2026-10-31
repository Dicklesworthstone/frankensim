//! Source-bound sensible heat between a finite gas and its finite adsorbent.
//!
//! Conductance includes the explicitly authored exposed geometry. It excludes
//! transported species enthalpy, adsorption binding energy and thermal storage;
//! those remain owned by the original gas and dry-plus-bound caloric laws.

use super::*;
use fs_matdb::PropertyKey;
use gas_network::SpeciesNetworkExternalHeat;

/// Thermal conductance [W/K], not conductivity [W/(m K)] or a film [W/(m² K)].
pub const GAS_SOLID_THERMAL_CONDUCTANCE_DIMS: Dims = Dims([2, 1, -3, -1, 0, 0]);

/// An explicitly calibrated, nonnegative conductance for this exposed path.
///
/// The complete typed query has exactly `gas_temperature`, `solid_temperature`
/// (absolute kelvin), and `pressure` (actual gas pressure). Scalar constants,
/// curves and surfaces must retain their original declared support. Other value
/// or applicability coordinates, including composition, coverage, legacy T/p,
/// frequency and rate, require an evolving owner and are refused here. The
/// selected immutable claim remains pinned through both actual endpoint queries.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbentGasHeatSource {
    material: MaterialCard,
    source: ResolvedMaterialStatePoint,
    property: String,
}

fn coordinates() -> [(&'static str, QuantitySpec); 3] {
    [
        ("gas_temperature", temp_quantity()),
        ("solid_temperature", temp_quantity()),
        (
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
        ),
    ]
}

impl AdsorbentGasHeatSource {
    pub fn try_new(
        material: MaterialCard,
        property: &PropertyKey,
        reference: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> core::result::Result<Self, AdsorptionError> {
        let owned = coordinates();
        if property.quantity() != QuantitySpec::dimensional(GAS_SOLID_THERMAL_CONDUCTANCE_DIMS)
            || reference.axes().len() != owned.len()
            || owned.iter().any(|(axis, quantity)| {
                reference.axis_quantities().get(*axis) != Some(quantity)
                    || !reference.axes().get(*axis).copied().is_some_and(positive)
            })
        {
            return Err(invalid(
                "gas/solid heat requires W/K and exactly typed actual gas/solid temperatures and pressure",
            ));
        }
        let requirement =
            ScalarPropertyRequirement::try_with_key(property, ScalarAdmissibility::NonNegative)?;
        let source = resolve_material_state_point(&material, reference, &[requirement], selection)?;
        let result = Self {
            material,
            source,
            property: property.name().to_owned(),
        };
        result.conductance(&result.source)?;
        Ok(result)
    }

    /// Original selected source, including its full reference query and receipt.
    pub fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }

    fn conductance(
        &self,
        source: &ResolvedMaterialStatePoint,
    ) -> core::result::Result<f64, AdsorptionError> {
        let property = source
            .property(&self.property)
            .ok_or(invalid("gas/solid heat lost its selected conductance"))?;
        let owned = coordinates();
        let admitted = |axis: &str| owned.iter().any(|(name, _)| *name == axis);
        if !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
            || property
                .source_axes()
                .iter()
                .any(|(axis, _)| !admitted(axis))
            || property
                .answer()
                .evidence
                .model
                .validity
                .bounds()
                .keys()
                .any(|axis| !admitted(axis))
        {
            return Err(invalid(
                "gas/solid conductance has an unowned physical coordinate or unsupported scalar response",
            ));
        }
        Ok(property.value_si())
    }

    fn at(
        &self,
        gas: &ClosedGasParcel,
        solid: &CompetitiveSolid,
    ) -> core::result::Result<(f64, ResolvedMaterialStatePoint), AdsorptionError> {
        let state = gas.state().state();
        let query = self
            .source
            .query()
            .with("gas_temperature", state.temperature)
            .and_then(|q| q.with("solid_temperature", solid.temperature_k))
            .and_then(|q| q.with("pressure", state.pressure))
            .map_err(|_| invalid("gas/solid heat endpoint coordinates are unrepresentable"))?;
        let source = self.source.at_query_point(&self.material, &query)?;
        let value = self.conductance(&source)?;
        Ok((value, source))
    }
}

/// One physical sensible-heat path to the gas already bound to this solid.
#[derive(Clone, Copy, Debug)]
pub struct AdsorbentGasHeatLink<'a> {
    pub id: u64,
    /// Original AdsorptionNetworkBinding identity (pore-cell index in a domain).
    pub solid: u64,
    pub source: &'a AdsorbentGasHeatSource,
}

/// Actual original-state forcing and both admitted material-source endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbentGasHeatSample {
    pub id: u64,
    pub solid: u64,
    pub parcel: u64,
    /// [gas, solid] at the accepted start, in kelvin.
    pub temperatures_k: [f64; 2],
    pub conductance_w_k: f64,
    pub power_gas_to_solid_w: f64,
    /// Signed heat, debited separately from gas and credited once to the solid.
    pub energy_gas_to_solid_j: f64,
    /// Original pinned claim queried at the actual [start, final] states.
    pub source: [ResolvedMaterialStatePoint; 2],
}

fn failure(id: u64, error: impl core::fmt::Display) -> Error {
    Error::Conduction {
        id,
        reason: error.to_string(),
    }
}

pub(super) struct GasHeat {
    pub into_solid_j: BTreeMap<u64, f64>,
    pub into_gas: Vec<SpeciesNetworkExternalHeat>,
    pub balance_residual_j: f64,
    pub samples: Vec<AdsorbentGasHeatSample>,
}

pub(super) fn sample(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    solids: &BTreeMap<u64, &AdsorptionNetworkBinding<'_>>,
    links: &[AdsorbentGasHeatLink<'_>],
    duration_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> NetworkResult<GasHeat> {
    if links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(Error::Invalid("too many physical gas/solid heat paths"));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(cancelled)?;
        if !solids.contains_key(&link.solid) || ordered.insert(link.id, link).is_some() {
            return Err(Error::Invalid(
                "gas/solid heat requires unique path identities and retained solid owners",
            ));
        }
    }
    let mut parts: BTreeMap<u64, Vec<f64>> = BTreeMap::new();
    let mut into_gas = Vec::with_capacity(links.len());
    let mut samples = Vec::with_capacity(links.len());
    for link in ordered.into_values() {
        poll(cancelled)?;
        let binding = solids[&link.solid];
        let gas = &parcels[&binding.parcel];
        let temperatures_k = [gas.state().state().temperature, binding.solid.temperature_k];
        let (conductance_w_k, source) = link
            .source
            .at(gas, binding.solid)
            .map_err(|e| failure(link.id, e))?;
        let difference = temperatures_k[0] - temperatures_k[1];
        let power = conductance_w_k * difference;
        let energy = power * duration_s;
        if !power.is_finite()
            || !energy.is_finite()
            || conductance_w_k != 0.0 && difference != 0.0 && power == 0.0
            || power != 0.0 && energy == 0.0
        {
            return Err(failure(
                link.id,
                "nonzero gas/solid sensible heat is unrepresentable",
            ));
        }
        parts.entry(link.solid).or_default().push(energy);
        into_gas.push(SpeciesNetworkExternalHeat {
            id: link.id,
            parcel: binding.parcel,
            heat_into_gas_j: -energy,
        });
        samples.push(AdsorbentGasHeatSample {
            id: link.id,
            solid: link.solid,
            parcel: binding.parcel,
            temperatures_k,
            conductance_w_k,
            power_gas_to_solid_w: power,
            energy_gas_to_solid_j: energy,
            source: [source.clone(), source],
        });
    }
    let into_solid_j: BTreeMap<_, _> = parts
        .into_iter()
        .map(|(id, parts)| sum(parts).map(|value| (id, value)))
        .collect::<core::result::Result<_, AdsorptionError>>()?;
    let balance_residual_j = sum(into_solid_j
        .values()
        .copied()
        .chain(into_gas.iter().map(|row| row.heat_into_gas_j)))?;
    Ok(GasHeat {
        into_solid_j,
        into_gas,
        balance_residual_j,
        samples,
    })
}

pub(super) fn admit_endpoint(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    prepared: &[(&AdsorptionNetworkBinding<'_>, PreparedAdsorption)],
    links: &[AdsorbentGasHeatLink<'_>],
    samples: &mut [AdsorbentGasHeatSample],
    cancelled: &mut impl FnMut() -> bool,
) -> NetworkResult<()> {
    let solids: BTreeMap<_, _> = prepared
        .iter()
        .map(|(binding, proposal)| (binding.id, &proposal.next))
        .collect();
    let sources: BTreeMap<_, _> = links.iter().map(|link| (link.id, link.source)).collect();
    for sample in samples {
        poll(cancelled)?;
        let (_, source) = sources[&sample.id]
            .at(&parcels[&sample.parcel], solids[&sample.solid])
            .map_err(|e| failure(sample.id, e))?;
        sample.source[1] = source;
    }
    Ok(())
}
