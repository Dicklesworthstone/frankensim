//! Source-temperature feedback for retained, explicitly heated viscosity layers.
//!
//! The existing liquid-film owner supplies the single physical clock and all
//! material/enthalpy admission. This module binds fixed original interfaces,
//! compiles each island once, and stages replacement heat fractions alongside
//! its hydraulic coefficient. No extra heat, fluid inventory or mixing is added.

use std::ops::Range;

use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;

use super::super::layered_film::{GapViscosityProfile, LayerHeatRoute, ViscosityLayer};
use super::{
    CompiledFilm, FilmBindingGroup, FilmBindingKind, HeatedImpactSystem, ImpactError,
    ResistiveFilm, invalid,
};

/// One original viscosity layer at its explicit finite liquid store's temperature.
///
/// Every layer in every channel of a selected island must appear exactly once.
/// The original layer heat destination is also its temperature reservoir. The
/// original normalized interfaces remain fixed as viscosity changes; source
/// transport data never creates or rescales a layer's mass or enthalpy.
#[derive(Clone, Debug)]
pub struct LiquidFilmLayerBinding<'m> {
    /// Existing normal-pressure film insertion index.
    pub film: usize,
    /// Existing hydraulic channel insertion index.
    pub channel: usize,
    /// Original lower-wall to upper-wall layer insertion index.
    pub layer: usize,
    /// Stable key of that layer's existing caloric heat destination.
    pub reservoir: String,
    /// Original pinned liquid source and unchanged nonthermal context.
    pub law: LiquidTemperatureLaw<'m>,
    /// Maximum `abs(mu_next-mu_old)/mu_old` for this individual layer in one
    /// sample, in (0,1]. This requests time refinement, not a viscosity clip.
    pub maximum_relative_viscosity_change: f64,
}

/// Immutable interface/destination template and canonical source locations.
/// Only its boundaries and destination indices are reused after construction;
/// accepted viscosity and work fractions live in the original heated owner.
pub(super) struct ChannelBinding {
    channel: usize,
    route: usize,
    rows: Range<usize>,
    template: LayerHeatRoute,
}

pub(super) fn append_groups(
    inner: &HeatedImpactSystem<'_>,
    bindings: &[LiquidFilmLayerBinding<'_>],
    offset: usize,
    groups: &mut Vec<FilmBindingGroup>,
) -> Result<(), ImpactError> {
    let mut start = 0;
    while start < bindings.len() {
        let film_index = bindings[start].film;
        let end = start + bindings[start..].partition_point(|row| row.film == film_index);
        let film = inner
            .mechanics()
            .squeeze_film(film_index)
            .ok_or_else(|| invalid("layer feedback names an absent normal-pressure film"))?;
        if groups.iter().any(|group| group.film == film_index) {
            return Err(invalid(
                "layer feedback cannot overlap uniform or channel assignments",
            ));
        }
        let mut channels = Vec::with_capacity(film.channel_count());
        let mut next = start;
        for channel in 0..film.channel_count() {
            if next == end || bindings[next].channel != channel {
                return Err(invalid(
                    "layer feedback requires every original pressure channel",
                ));
            }
            let next_end = next + bindings[next..end].partition_point(|row| row.channel == channel);
            let route = inner.film_routes.iter().position(|route|
                route.is_squeeze() && route.film == film_index
                    && route.channel == Some(channel) && route.layers.is_some()
            ).ok_or_else(|| invalid(
                "layer feedback requires an existing retained profile and heat route for every channel",
            ))?;
            let template = inner.film_routes[route]
                .layers
                .as_ref()
                .expect("selected layer route");
            if next_end - next != template.profile.layers().len()
                || bindings[next..next_end]
                    .iter()
                    .enumerate()
                    .any(|(layer, binding)| {
                        binding.layer != layer
                            || inner.reservoirs()[template.destinations[layer]].key()
                                != binding.reservoir.as_str()
                    })
            {
                return Err(invalid(
                    "layer feedback needs complete unique layer coverage at each original heat destination",
                ));
            }
            channels.push(ChannelBinding {
                channel,
                route,
                rows: offset + next..offset + next_end,
                template: template.clone(),
            });
            next = next_end;
        }
        if next != end {
            return Err(invalid(
                "layer feedback names an unknown or duplicate pressure channel",
            ));
        }
        groups.push(FilmBindingGroup {
            film: film_index,
            rows: offset + start..offset + end,
            kind: FilmBindingKind::Layers(channels),
        });
        start = end;
    }
    Ok(())
}

pub(super) fn compile(
    geometry: &ResistiveFilm,
    sources: &[ResolvedLiquidState],
    channels: &[ChannelBinding],
) -> Result<CompiledFilm, ImpactError> {
    let mut profiles = Vec::with_capacity(channels.len());
    let mut viscosity = Vec::with_capacity(channels.len());
    for channel in channels {
        let layers = channel
            .template
            .profile
            .layers()
            .iter()
            .zip(&sources[channel.rows.clone()])
            .map(|(old, source)| ViscosityLayer {
                upper_fraction: old.upper_fraction,
                dynamic_viscosity_pa_s: source.dynamic_viscosity_pa_s(),
            })
            .collect();
        let profile = GapViscosityProfile::new(layers).map_err(|error| invalid(error.0))?;
        viscosity.push(profile.hydraulic_equivalent_viscosity_pa_s());
        profiles.push(profile);
    }
    let operator = geometry
        .with_channel_viscosities(&viscosity)
        .map_err(|error| invalid(error.0))?;
    let mut layer_routes = Vec::with_capacity(channels.len());
    for (channel, profile) in channels.iter().zip(profiles) {
        layer_routes.push((
            channel.route,
            channel
                .template
                .with_profile(&operator, channel.channel, profile)?,
        ));
    }
    Ok(CompiledFilm {
        operator,
        layer_routes,
    })
}

/// All source/profile/operator and native time/heat gates already succeeded.
/// Fixed route indices cannot change through the public callback surfaces.
pub(super) fn publish(inner: &mut HeatedImpactSystem<'_>, updates: &mut [(usize, LayerHeatRoute)]) {
    for (index, update) in updates {
        std::mem::swap(
            inner.film_routes[*index]
                .layers
                .as_mut()
                .expect("bound original layer route"),
            update,
        );
    }
}
