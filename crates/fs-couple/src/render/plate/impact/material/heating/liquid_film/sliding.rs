//! Original sliding layers bind to the existing finite thermal destinations.
//! The enclosing liquid-film owner resolves source states and accepts heat;
//! ShearedFilm alone recompiles pressure, entrainment, traction and layer work.

use std::ops::Range;

use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_tribo::resistive_film::shear::ShearedFilm;

use super::{HeatedImpactSystem, ImpactError, invalid};

/// One original sliding-film layer at its finite liquid store's temperature.
///
/// Every layer of every channel in a selected sliding island is bound once.
/// Uniform channels have one implicit layer numbered zero. Sharing a store
/// between rows gives uniform temperature without replacing different material
/// laws or averaging their viscosities. The source and caloric material must
/// agree, and the store must be that layer's original explicit heat destination.
#[derive(Clone, Debug)]
pub struct LiquidSlidingFilmLayerBinding<'m> {
    /// Original SLIDING-film index, independent of normal-pressure film indices.
    pub film: usize,
    /// Original channel insertion index.
    pub channel: usize,
    /// Original lower-to-upper layer index, zero for an implicit uniform layer.
    pub layer: usize,
    /// Existing finite heat destination setting the layer's temperature.
    pub reservoir: String,
    /// Original source-backed liquid law, retaining fixed nonthermal context.
    pub law: LiquidTemperatureLaw<'m>,
    /// Maximum `abs(mu_next-mu_old)/mu_old` for this layer in one mechanical
    /// sample, in (0,1]. Exceeding it requests refinement and never clips mu.
    pub maximum_relative_viscosity_change: f64,
}

pub(super) struct BindingGroup {
    pub(super) film: usize,
    rows: Range<usize>,
}

impl BindingGroup {
    pub(super) fn compile(
        &self,
        film: &ShearedFilm,
        sources: &[ResolvedLiquidState],
    ) -> Result<ShearedFilm, ImpactError> {
        let viscosities: Vec<_> = sources[self.rows.clone()]
            .iter()
            .map(ResolvedLiquidState::dynamic_viscosity_pa_s)
            .collect();
        film.with_layer_viscosities(&viscosities)
            .map_err(|error| invalid(error.0))
    }
}

pub(super) fn groups(
    inner: &HeatedImpactSystem<'_>,
    bindings: &[LiquidSlidingFilmLayerBinding<'_>],
    offset: usize,
) -> Result<Vec<BindingGroup>, ImpactError> {
    let mut result = Vec::new();
    let mut start = 0;
    while start < bindings.len() {
        let film_index = bindings[start].film;
        let end = start + bindings[start..].partition_point(|binding| binding.film == film_index);
        let film = inner
            .mechanics()
            .sliding_film(film_index)
            .ok_or_else(|| invalid("liquid feedback names an absent sliding pressure island"))?;
        if end - start != film.layer_count() {
            return Err(invalid(
                "sliding feedback requires every original channel and layer exactly once",
            ));
        }
        let mut next = start;
        for channel in 0..film.pressure_operator().channel_count() {
            let range = film
                .channel_layer_range(channel)
                .expect("original sliding channel");
            let destinations = inner
                .film_routes
                .iter()
                .find(|route| {
                    route.sliding_layers.is_some()
                        && route.film == film_index
                        && route.channel == Some(channel)
                })
                .and_then(|route| route.sliding_layers.as_ref())
                .ok_or_else(|| {
                    invalid(
                        "sliding feedback requires an explicit layer heat route for every channel",
                    )
                })?;
            if destinations.len() != range.len() {
                return Err(invalid(
                    "sliding feedback must retain original layer heat coverage",
                ));
            }
            for (layer, &destination) in destinations.iter().enumerate() {
                let binding = &bindings[next];
                if binding.channel != channel
                    || binding.layer != layer
                    || binding.reservoir.as_str() != inner.reservoirs()[destination].key()
                {
                    return Err(invalid(
                        "sliding feedback needs unique original layers at their existing heat destinations",
                    ));
                }
                next += 1;
            }
        }
        debug_assert_eq!(next, end);
        result.push(BindingGroup {
            film: film_index,
            rows: offset + start..offset + end,
        });
        start = end;
    }
    Ok(result)
}
