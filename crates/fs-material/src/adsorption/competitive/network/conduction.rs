//! Original integrated k(T) and ordered interfaces between finite adsorbents.
//!
//! Link power is held at the accepted starting state for one native adsorption
//! interval. Its equal/opposite energy enters the original dry+bound caloric
//! inversion with sorption heat. This adds no thermal storage, time integrator,
//! inferred skeleton conductivity, or second gas/solid heat account.

use super::*;
use crate::phase::{EquilibriumPhaseState, SolidLiquidPhase};
use crate::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use crate::thermal::{
    conduction::{
        LayerConductionConfig, LayerConductionResponse, LayerGeometry,
        source::ResolvedConductionLayer,
    },
    interface::{AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, ThermalInterfacePoint},
    source::ThermalMaterialError,
};
use fs_blake3::ContentHash;
use fs_matdb::PropertyLineSlice;
use fs_thermal::enthalpy::DepositionError;

/// One pinned original ordered interface, refreshed at both actual solid states.
#[derive(Clone, Debug)]
pub struct AdsorbentThermalContact<'a> {
    point: ThermalInterfacePoint<'a>,
    source: ResolvedInterfaceStatePoint,
}
impl<'a> AdsorbentThermalContact<'a> {
    pub fn try_new(
        point: ThermalInterfacePoint<'a>,
        states: [EquilibriumPhaseState; 2],
        selection: MaterialPropertySelection,
    ) -> core::result::Result<Self, ThermalMaterialError> {
        let source = point.resolve_resistance(states, selection)?;
        Ok(Self { point, source })
    }
    pub fn source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
    fn at(
        &self,
        states: [EquilibriumPhaseState; 2],
    ) -> core::result::Result<ResolvedInterfaceStatePoint, ThermalMaterialError> {
        let point = self.point.at(states)?;
        self.source
            .at_query_point(self.point.interface, &point)
            .map_err(ThermalMaterialError::Material)
    }
}

/// An explicitly authored physical heat path between two native site populations.
/// Layer caloric objects must be the SAME retained dry skeleton charts. Geometry
/// is independent of gas throats and does not infer pore size, dry mass or density.
#[derive(Clone, Copy, Debug)]
pub struct AdsorbentConductionLink<'a, 'm> {
    pub id: u64,
    /// Original AdsorptionNetworkBinding IDs, in physical layer/interface order.
    pub solids: [u64; 2],
    pub layers: [&'a ResolvedConductionLayer<'m>; 2],
    pub geometry: LayerGeometry,
    /// None explicitly declares perfect contact. A finite interface must retain
    /// its original ordered source through AdsorbentThermalContact.
    pub contact: Option<&'a AdsorbentThermalContact<'m>>,
    pub solve: LayerConductionConfig,
}

/// Physical power, applied heat and sources from the actual starting snapshot.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbentConductionSample {
    pub id: u64,
    pub solids: [u64; 2],
    pub temperatures_k: [f64; 2],
    pub response: LayerConductionResponse,
    pub geometry: LayerGeometry,
    /// Signed internal heat; positive leaves the first solid and enters the second.
    pub energy_a_to_b_j: f64,
    pub source_lines: [Arc<PropertyLineSlice>; 2],
    pub material_card_identities: [ContentHash; 2],
    pub caloric_curve_identities: [ContentHash; 2],
    pub interface: Option<ResolvedInterfaceStatePoint>,
}

fn failure(id: u64, error: impl core::fmt::Display) -> Error {
    Error::Conduction {
        id,
        reason: error.to_string(),
    }
}
fn flux_failure(id: u64, error: DepositionError) -> Error {
    if matches!(error, DepositionError::Cancelled) {
        Error::Cancelled
    } else {
        failure(id, error)
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> core::result::Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

impl AdsorbentConductionLink<'_, '_> {
    fn states(
        &self,
        solids: [&CompetitiveSolid; 2],
    ) -> NetworkResult<([f64; 2], Option<ResolvedInterfaceStatePoint>)> {
        let at = |side: usize| {
            let dry = solids[side].model.dry();
            if self.layers[side].phase() != SolidLiquidPhase::Solid
                || !std::ptr::eq(self.layers[side].caloric_curve(), dry.curve())
            {
                return Err(failure(
                    self.id,
                    "conduction must retain the original dry skeleton caloric object",
                ));
            }
            let state = dry
                .state_at_temperature(solids[side].temperature_k)
                .map_err(|error| failure(self.id, error))?;
            self.layers[side]
                .temperature_at_enthalpy(state.specific_enthalpy_j_kg())
                .map_err(|error| failure(self.id, error))?;
            Ok(state)
        };
        let states = [at(0)?, at(1)?];
        let temperatures = states.map(|state| state.temperature_k());
        // The original integrated layer law needs continuous support over the
        // complete interval, including at the final state of a completed window.
        for layer in self.layers {
            for temperature in temperatures {
                layer
                    .conductivity_curve()
                    .conductivity_w_m_k(temperature)
                    .map_err(|error| flux_failure(self.id, error))?;
            }
        }
        let interface = self
            .contact
            .map(|contact| contact.at(states))
            .transpose()
            .map_err(|error| failure(self.id, error))?;
        Ok((temperatures, interface))
    }
}

pub(super) struct ConductionHeat {
    pub into_solid_j: BTreeMap<u64, f64>,
    pub balance_residual_j: f64,
    pub samples: Vec<AdsorbentConductionSample>,
}

pub(super) fn sample(
    solids: &BTreeMap<u64, &AdsorptionNetworkBinding<'_>>,
    links: &[AdsorbentConductionLink<'_, '_>],
    duration_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> NetworkResult<ConductionHeat> {
    if links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(Error::Invalid(
            "too many physical adsorbent conduction links",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(cancelled)?;
        if link.solids[0] == link.solids[1]
            || link.solids.iter().any(|id| !solids.contains_key(id))
            || ordered.insert(link.id, link).is_some()
        {
            return Err(failure(
                link.id,
                "conduction requires unique link IDs and two existing distinct solids",
            ));
        }
    }
    let mut parts: BTreeMap<u64, Vec<f64>> = BTreeMap::new();
    let mut samples = Vec::with_capacity(links.len());
    for link in ordered.into_values() {
        poll(cancelled)?;
        let (temperatures_k, interface) = link.states(link.solids.map(|id| solids[&id].solid))?;
        let resistance = interface
            .as_ref()
            .map(|source| {
                source
                    .property(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
                    .map(|property| property.value_si())
                    .ok_or_else(|| {
                        failure(
                            link.id,
                            "original thermal interface lost its selected resistance",
                        )
                    })
            })
            .transpose()?
            .unwrap_or(0.0);
        let law = link
            .geometry
            .conduction(
                link.layers.map(|layer| layer.conductivity_curve()),
                resistance,
            )
            .map_err(|error| flux_failure(link.id, error))?;
        let response = law
            .evaluate(temperatures_k, link.solve, || checkpoint(cancelled))
            .map_err(|error| flux_failure(link.id, error))?;
        let energy = response.power_a_to_b_w * duration_s;
        if !energy.is_finite() || response.power_a_to_b_w != 0.0 && energy == 0.0 {
            return Err(failure(
                link.id,
                "nonzero conductive heat is unrepresentable",
            ));
        }
        parts.entry(link.solids[0]).or_default().push(-energy);
        parts.entry(link.solids[1]).or_default().push(energy);
        samples.push(AdsorbentConductionSample {
            id: link.id,
            solids: link.solids,
            temperatures_k,
            response,
            geometry: link.geometry,
            energy_a_to_b_j: energy,
            source_lines: link.layers.map(|layer| Arc::clone(layer.source_line())),
            material_card_identities: link
                .layers
                .map(|layer| layer.caloric_curve().material_card_identity()),
            caloric_curve_identities: link.layers.map(|layer| layer.caloric_curve().identity()),
            interface,
        });
    }
    let into_solid_j: BTreeMap<_, _> = parts
        .into_iter()
        .map(|(id, parts)| sum(parts).map(|value| (id, value)))
        .collect::<core::result::Result<_, AdsorptionError>>()?;
    let balance_residual_j = sum(into_solid_j.values().copied())?;
    Ok(ConductionHeat {
        into_solid_j,
        balance_residual_j,
        samples,
    })
}

pub(super) fn admit_endpoint(
    prepared: &[(&AdsorptionNetworkBinding<'_>, PreparedAdsorption)],
    links: &[AdsorbentConductionLink<'_, '_>],
    cancelled: &mut impl FnMut() -> bool,
) -> NetworkResult<()> {
    let solids: BTreeMap<_, _> = prepared
        .iter()
        .map(|(binding, proposal)| (binding.id, &proposal.next))
        .collect();
    for link in links {
        poll(cancelled)?;
        link.states(link.solids.map(|id| solids[&id]))?;
    }
    Ok(())
}
