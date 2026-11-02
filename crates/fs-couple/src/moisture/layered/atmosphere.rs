//! Layered transport against the original finite, source-qualified gas owner.
//!
//! All exposed exterior ends use one candidate endpoint RH. Each trial starts
//! from the accepted layered fields; internal material interfaces retain their
//! own source partition. The original atmosphere owner alone reconstructs gas
//! composition/pressure and accepts vapor, carrier and total mass conservation.
//! This is isothermal, fixed-volume transport, without sorption heat or
//! condensation. Source dependence on changing gas pressure/composition refuses.
//! Each stack's explicit exterior transfer laws are retained during every trial.

use super::*;
use crate::moisture::atmosphere::{
    AtmosphereMoistureCandidate, AtmosphereMoistureError, AtmosphereMoistureSolve,
    AtmosphereMoistureStep, FiniteMoistureAtmosphere, SharedAtmosphereMoistureCandidate,
    SharedMoistureSolid, admit_participants,
};
use crate::moisture::reservoir::water_sum;
use crate::moisture::surface::admit_environment;

/// One heterogeneous stack sharing the same finite atmospheric water inventory.
pub struct LayeredMoistureParticipant<'a> {
    pub id: u64,
    pub moisture: &'a MoistureLayeredSlab,
    /// Left/right exterior ends: true exchanges with the atmosphere, false seals.
    pub exposed: [bool; 2],
    pub transport: MoistureStepConfig,
    pub nonlinear: MoistureNonlinearConfig,
}

pub type LayeredAtmosphereMoistureCandidate =
    AtmosphereMoistureCandidate<MoistureLayeredSlab, MoistureLayeredReport>;
pub type SharedLayeredAtmosphereMoistureCandidate =
    SharedAtmosphereMoistureCandidate<MoistureLayeredSlab, MoistureLayeredReport>;

fn water_change(old: &MoistureLayeredSlab, next: &MoistureLayeredSlab) -> f64 {
    water_sum(
        old.layers
            .iter()
            .zip(&next.layers)
            .flat_map(|(before, after)| {
                before
                    .slab
                    .moisture
                    .iter()
                    .zip(&after.slab.moisture)
                    .map(move |(old, new)| (new - old) * before.slab.dry_mass_per_cell_kg)
            }),
    )
}

/// Resolve the same pinned source at actual solid surface RH for the gas-step
/// report. This does not replace or refit the layer's prepared isotherm.
fn endpoint_sorption(
    layer: &MoistureLayer,
    activity: f64,
) -> Result<ResolvedMaterialStatePoint, MoistureError> {
    let claim = layer
        .slab
        .card
        .claims()
        .claim(layer.sorption.selected_claim())
        .expect("admitted layer sorption source");
    let point = layer
        .sorption
        .source_point()
        .clone()
        .with("relative_humidity", activity)
        .map_err(|_| MoistureError::Invalid("invalid accepted layer humidity"))?;
    let requirement = ScalarPropertyRequirement::try_with_quantity(
        claim.key.name(),
        quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
        ScalarAdmissibility::NonNegative,
    )?;
    Ok(resolve_material_state_point(
        &layer.slab.card,
        &point,
        &[requirement],
        MaterialPropertySelection::PinnedByProperty(vec![(
            claim.key.name().to_owned(),
            layer.sorption.selected_claim(),
        )]),
    )?)
}

impl FiniteMoistureAtmosphere {
    /// Couple one whole layered stack to this finite vapor inventory. The
    /// exterior source law comes from each end layer; no separate sorption law
    /// may silently override the material's interface partition relation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        let participant = LayeredMoistureParticipant {
            id: 0,
            moisture,
            exposed,
            transport,
            nonlinear,
        };
        let candidate = self.advance_shared_layered(
            &[participant],
            dt_s,
            moisture.activity.len(),
            controls,
            cx,
        )?;
        let (atmosphere, mut solids, shared) = candidate.into_parts();
        let solid = solids.pop().expect("one admitted layered stack");
        let mut sorption = [None, None];
        for (side, layer) in [0, moisture.layers.len() - 1].into_iter().enumerate() {
            poll(cx)?;
            if let Some(activity) = solid.transport.surface_activities[side] {
                sorption[side] = Some(endpoint_sorption(&moisture.layers[layer], activity)?);
            }
        }
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

    /// Simultaneous uptake/desorption of separately identified layered stacks
    /// debits one vapor inventory once. Native layered solves, source limits and
    /// interface reports survive; stable IDs determine deterministic summation.
    pub fn advance_shared_layered(
        &self,
        participants: &[LayeredMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(AtmosphereMoistureError::Invalid(
                "shared-atmosphere stack IDs must be unique",
            ));
        }
        let time_s = ordered[0].moisture.time_s();
        let temperature = self.gas().state().state().temperature;
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.moisture.activity.len())
                .filter(|&n| n <= maximum_total_cells)
                .ok_or(AtmosphereMoistureError::Invalid(
                    "shared-atmosphere total layered cell budget exhausted",
                ))?;
            if body.moisture.time_s().to_bits() != time_s.to_bits()
                || body.moisture.layers[0]
                    .sorption
                    .temperature_kelvin()
                    .to_bits()
                    != temperature.to_bits()
            {
                return Err(AtmosphereMoistureError::Invalid(
                    "layered stacks need the same accepted clock and gas temperature",
                ));
            }
            if effective_exposure(body).into_iter().any(|face| face) {
                let supported = body.moisture.activity_bounds;
                if controls.relative_humidity_bracket[0] < supported[0]
                    || controls.relative_humidity_bracket[1] > supported[1]
                {
                    return Err(AtmosphereMoistureError::Invalid(
                        "finite humidity bracket exceeds a layered source domain",
                    ));
                }
                for layer in &body.moisture.layers {
                    poll(cx)?;
                    // Pressure/composition cannot remain hidden fixed inputs
                    // while this very gas owner changes them during exchange.
                    admit_environment(layer.sorption.source_point())?;
                }
            }
        }
        let all_sealed = ordered
            .iter()
            .all(|body| effective_exposure(body) == [false; 2]);
        let evaluate = |rh| -> Result<_, AtmosphereMoistureError> {
            let mut solids = Vec::with_capacity(ordered.len());
            for body in &ordered {
                poll(cx)?;
                let boundary = effective_exposure(body).map(|exposed| {
                    if exposed {
                        MoistureLayeredFace::PrescribedRelativeHumidity(rh)
                    } else {
                        MoistureLayeredFace::Sealed
                    }
                });
                let (moisture, transport) =
                    body.moisture
                        .advance(dt_s, &boundary, body.transport, body.nonlinear, cx)?;
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
        };
        let (atmosphere, solids, step) =
            self.advance_exchange(time_s, all_sealed, controls, cx, evaluate)?;
        Ok(SharedAtmosphereMoistureCandidate::from_accepted(
            atmosphere, solids, step,
        ))
    }
}

fn effective_exposure(body: &LayeredMoistureParticipant<'_>) -> [bool; 2] {
    std::array::from_fn(|side| {
        body.exposed[side]
            && body.moisture.exterior_transfer[side] != MoistureActivityTransfer::Sealed
    })
}
