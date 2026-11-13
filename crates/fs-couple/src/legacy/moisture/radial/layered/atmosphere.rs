//! Coated radial bodies exchanging water with an original finite atmosphere.
//!
//! One native layered solve spans every material at each candidate gas RH.
//! Actual shell inventories debit the original gas once; an interior material
//! interface never creates a gas port or an additional sorption boundary. The
//! finite-vapor owner reconstructs pressure/composition and admits water,
//! carrier and total mass before publishing all accepted bodies together.
//!
//! Gas free volume, temperature and reference geometry remain fixed. This does
//! not infer carrier permeation, condensation, hysteresis or sorption heat.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureCandidate, AtmosphereMoistureError, AtmosphereMoistureSolve,
        AtmosphereMoistureStep, FiniteMoistureAtmosphere, SharedAtmosphereMoistureCandidate,
        SharedMoistureSolid, admit_participants,
    },
    reservoir::water_sum,
    surface::admit_environment,
};

pub use crate::moisture::layered::atmosphere::{
    LayeredRadialAtmosphereSdirk2Candidate, SharedLayeredRadialAtmosphereSdirk2Candidate,
};

/// A complete independently sourced coated body in one common gas inventory.
pub struct LayeredRadialMoistureParticipant<'a> {
    /// Stable, unique physical-inventory ID; fixes summation and result order.
    pub id: u64,
    pub moisture: &'a MoistureLayeredRadial,
    /// Inner/outer physical surfaces; false seals. A solid center must be false.
    pub exposed: [bool; 2],
    pub transport: MoistureStepConfig,
    pub nonlinear: MoistureNonlinearConfig,
}

/// One original finite gas and complete heterogeneous radial body.
pub type LayeredRadialAtmosphereMoistureCandidate =
    AtmosphereMoistureCandidate<MoistureLayeredRadial, MoistureLayeredReport>;
/// Independently sourced coated radial bodies sharing one finite vapor account.
pub type SharedLayeredRadialAtmosphereMoistureCandidate =
    SharedAtmosphereMoistureCandidate<MoistureLayeredRadial, MoistureLayeredReport>;

pub(in crate::moisture) fn water_change(
    old: &MoistureLayeredRadial,
    next: &MoistureLayeredRadial,
) -> f64 {
    water_sum(
        old.layers.iter().zip(&next.layers).flat_map(|(before, after)| {
            before
                .radial
                .moisture
                .iter()
                .zip(&after.radial.moisture)
                .zip(&before.radial.dry_masses_kg)
                .map(|((&old, &new), &dry)| (new - old) * dry)
        }),
    )
}

/// A closed exterior law stays closed. Otherwise every material must support
/// the common gas trial and its current fixed temperature, without hidden gas
/// pressure/composition coordinates that would change during the exchange.
fn admit_exchange(
    moisture: &MoistureLayeredRadial,
    temperature: f64,
    exposed: [bool; 2],
    controls: AtmosphereMoistureSolve,
    cx: &Cx<'_>,
) -> Result<[bool; 2], AtmosphereMoistureError> {
    poll(cx)?;
    if exposed[0] && moisture.boundary_areas[0] == 0.0 {
        return Err(AtmosphereMoistureError::Invalid(
            "a solid radial center has zero area and cannot contact a finite atmosphere",
        ));
    }
    if moisture.layers[0].sorption.temperature_kelvin().to_bits() != temperature.to_bits() {
        return Err(AtmosphereMoistureError::Invalid(
            "layered radial bodies and gas need the same prescribed temperature",
        ));
    }
    let exposed = std::array::from_fn(|side| {
        exposed[side] && moisture.exterior_transfer[side] != MoistureActivityTransfer::Sealed
    });
    if exposed.into_iter().any(|face| face) {
        let supported = moisture.activity_bounds;
        if controls.relative_humidity_bracket[0] < supported[0]
            || controls.relative_humidity_bracket[1] > supported[1]
        {
            return Err(AtmosphereMoistureError::Invalid(
                "finite humidity bracket exceeds a layered radial source domain",
            ));
        }
        for layer in &moisture.layers {
            poll(cx)?;
            admit_environment(layer.sorption.source_point())?;
        }
    }
    Ok(exposed)
}

/// Re-query the selected end-material claim at the actual solid-side surface
/// activity. An explicit exterior resistance permits a difference from gas RH.
pub(in crate::moisture) fn surface_sorption(
    moisture: &MoistureLayeredRadial,
    transport: &MoistureLayeredReport,
    cx: &Cx<'_>,
) -> Result<[Option<ResolvedMaterialStatePoint>; 2], MoistureError> {
    let mut sorption = [None, None];
    for (side, layer) in [0, moisture.layers.len() - 1].into_iter().enumerate() {
        poll(cx)?;
        if let Some(activity) = transport.surface_activities[side] {
            let layer = &moisture.layers[layer];
            let claim = layer
                .radial
                .card
                .claims()
                .claim(layer.sorption.selected_claim())
                .expect("admitted radial layer sorption source");
            let point = layer
                .sorption
                .source_point()
                .clone()
                .with("relative_humidity", activity)
                .map_err(|_| MoistureError::Invalid("invalid accepted radial surface activity"))?;
            let requirement = ScalarPropertyRequirement::try_with_quantity(
                claim.key.name(),
                quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
                ScalarAdmissibility::NonNegative,
            )?;
            sorption[side] = Some(resolve_material_state_point(
                &layer.radial.card,
                &point,
                &[requirement],
                MaterialPropertySelection::PinnedByProperty(vec![(
                    claim.key.name().to_owned(),
                    layer.sorption.selected_claim(),
                )]),
            )?);
        }
    }
    Ok(sorption)
}

impl FiniteMoistureAtmosphere {
    /// Couple one entire coated radial body to this finite vapor inventory.
    /// Its original end-layer isotherms and explicit exterior/interface laws
    /// survive every trial; no additional boundary sorption model overrides them.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial(
        &self,
        moisture: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        let participant = LayeredRadialMoistureParticipant {
            id: 0,
            moisture,
            exposed,
            transport,
            nonlinear,
        };
        let (atmosphere, mut solids, shared) = self
            .advance_shared_layered_radial(
                &[participant],
                dt_s,
                moisture.activity.len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted layered radial body");
        let sorption = surface_sorption(moisture, &solid.transport, cx)?;
        let step = AtmosphereMoistureStep {
            sorption,
            transport: solid.transport,
            relative_humidity: shared.relative_humidity,
            vapor_mass_kg: shared.vapor_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: shared.maximum_carrier_mass_residual_kg,
            gas_internal_energy_change_j: shared.gas_internal_energy_change_j,
            iterations: shared.iterations,
        };
        poll(cx)?;
        Ok(AtmosphereMoistureCandidate::from_accepted(
            atmosphere,
            solid.moisture,
            step,
        ))
    }

    /// Simultaneous uptake and release of coated spheres or cylinders sharing
    /// one well-mixed finite gas. The total budget includes all layers/shells,
    /// including sealed bodies that still redistribute water internally. Stable
    /// IDs determine canonical order; every body starts from its accepted state
    /// on every original gas-root trial. Failure leaves all inputs unchanged.
    pub fn advance_shared_layered_radial(
        &self,
        participants: &[LayeredRadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredRadialAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(AtmosphereMoistureError::Invalid(
                "layered radial specimen IDs must be unique",
            ));
        }
        let time_s = ordered[0].moisture.time_s();
        let temperature = self.gas().state().state().temperature;
        let mut cells = 0usize;
        let mut exposure = Vec::with_capacity(ordered.len());
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.moisture.activity.len())
                .filter(|&n| n <= maximum_total_cells)
                .ok_or(AtmosphereMoistureError::Invalid(
                    "shared layered radial cell budget exhausted",
                ))?;
            if body.moisture.time_s().to_bits() != time_s.to_bits() {
                return Err(AtmosphereMoistureError::Invalid(
                    "layered radial specimens need the same accepted clock",
                ));
            }
            exposure.push(admit_exchange(
                body.moisture,
                temperature,
                body.exposed,
                controls,
                cx,
            )?);
        }
        let all_sealed = exposure.iter().all(|faces| *faces == [false; 2]);
        let (atmosphere, solids, step) =
            self.advance_exchange(time_s, all_sealed, controls, cx, |rh| {
                let mut solids = Vec::with_capacity(ordered.len());
                for (body, exposed) in ordered.iter().zip(&exposure) {
                    poll(cx)?;
                    let boundary = exposed.map(|exposed| {
                        if exposed {
                            MoistureLayeredFace::PrescribedRelativeHumidity(rh)
                        } else {
                            MoistureLayeredFace::Sealed
                        }
                    });
                    let (moisture, transport) = body.moisture.advance(
                        dt_s,
                        &boundary,
                        body.transport,
                        body.nonlinear,
                        cx,
                    )?;
                    solids.push(SharedMoistureSolid {
                        id: body.id,
                        moisture,
                        transport,
                    });
                }
                let change = water_sum(
                    ordered
                        .iter()
                        .zip(&solids)
                        .map(|(old, next)| water_change(old.moisture, &next.moisture)),
                );
                let end_time = solids[0].moisture.time_s();
                Ok((solids, change, end_time))
            })?;
        poll(cx)?;
        Ok(SharedAtmosphereMoistureCandidate::from_accepted(
            atmosphere, solids, step,
        ))
    }
}
