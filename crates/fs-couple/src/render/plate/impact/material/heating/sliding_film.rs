//! Complete sliding-film layer work deposited on the original impact clock.
//!
//! Native passive and pressure-fed sliding films report integral(tau^2/mu)
//! for every retained layer,
//! including the signed pressure/Couette cross term. Its original pressure,
//! force, tangent and phase-space stepping remain the mechanical owners.
//! Finite thermal destinations are explicit; no wall partition, layer mass,
//! material exchange or source-temperature feedback is inferred here.

use fs_tribo::resistive_film::shear::{MAX_LAYER_REPORTS, ShearChannelReport};
use fs_tribo::resistive_film::{MAX_CELLS, MAX_CHANNELS};

use super::{
    CaloricReservoir, FilmRoute, FilmRouteKind, ImpactError, ImpactSystem, MAX_IMPACT_MODES, PreparedImpactSystem,
    invalid,
};

/// Complete finite destinations for one original sliding passage's layers.
///
/// Passive routes use `sliding_film_count`; the driven-route argument of
/// `HeatedImpactSystem::new_with_driven_sliding_losses` instead uses
/// `driven_sliding_film_count`. These and squeeze-film indices are separate.
/// Each selected channel is routed once; unselected channels remain
/// mechanical losses outside this thermal boundary. A uniform channel has
/// one implicit layer. Explicit profile layers retain lower-to-upper order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SlidingFilmLayerHeatRoute {
    /// Original insertion index within the selected passive or driven namespace.
    pub film: usize,
    /// Original channel insertion index inside that sliding island.
    pub channel: usize,
    /// One existing finite reservoir key per original layer; keys may repeat.
    pub layer_destinations: Vec<String>,
}

pub(super) fn compile(
    mechanics: &PreparedImpactSystem,
    route: SlidingFilmLayerHeatRoute,
    stores: &[CaloricReservoir<'_>],
    compiled: &[FilmRoute],
) -> Result<FilmRoute, ImpactError> {
    compile_kind(mechanics, route, stores, compiled, FilmRouteKind::Sliding)
}

pub(super) fn compile_driven(
    mechanics: &PreparedImpactSystem,
    route: SlidingFilmLayerHeatRoute,
    stores: &[CaloricReservoir<'_>],
    compiled: &[FilmRoute],
) -> Result<FilmRoute, ImpactError> {
    compile_kind(mechanics, route, stores, compiled, FilmRouteKind::DrivenSliding)
}

fn compile_kind(
    mechanics: &PreparedImpactSystem,
    route: SlidingFilmLayerHeatRoute,
    stores: &[CaloricReservoir<'_>],
    compiled: &[FilmRoute],
    kind: FilmRouteKind,
) -> Result<FilmRoute, ImpactError> {
    let range = match kind {
        FilmRouteKind::Sliding => mechanics
            .sliding_film(route.film)
            .ok_or_else(|| invalid("sliding heat names an absent sliding pressure island"))?
            .channel_layer_range(route.channel),
        FilmRouteKind::DrivenSliding => mechanics
            .driven_sliding_film(route.film)
            .ok_or_else(|| invalid("driven sliding heat names an absent pressure island"))?
            .channel_layer_range(route.channel),
        _ => return Err(invalid("sliding heat requires its original film namespace")),
    }
    .ok_or_else(|| invalid("sliding heat names an absent original passage"))?;
    if route.layer_destinations.len() != range.len() {
        return Err(invalid(
            "sliding heat needs exactly one destination per original layer",
        ));
    }
    if compiled.iter().any(|old| {
        old.kind == kind && old.film == route.film && old.channel == Some(route.channel)
    }) {
        return Err(invalid("sliding channel heat must be routed at most once"));
    }
    let mut destinations = Vec::with_capacity(route.layer_destinations.len());
    for key in &route.layer_destinations {
        destinations.push(
            stores
                .iter()
                .position(|store| store.key() == key.as_str())
                .ok_or_else(|| invalid("sliding layer heat names an absent finite reservoir"))?,
        );
    }
    let first = destinations[0];
    Ok(FilmRoute {
        kind,
        film: route.film,
        channel: Some(route.channel),
        a: first,
        b: first,
        fraction: 1.0,
        layers: None,
        sliding_layers: Some(destinations),
    })
}

/// Keep the bounded large layer report off the legacy normal-film/contact stack.
/// Only the original selected sliding islands invoke this helper.
#[inline(never)]
pub(super) fn stage(
    system: &ImpactSystem,
    midpoint: &[f64],
    velocity: &[f64],
    routes: &[FilmRoute],
    heat_j: &mut [f64],
) -> Result<f64, ImpactError> {
    let n = system.modes;
    let mut configuration = [0.0; MAX_IMPACT_MODES];
    for i in 0..n {
        configuration[i] = midpoint[2 * i];
    }
    let mut force = [0.0; MAX_IMPACT_MODES];
    let mut pressure = [0.0; MAX_CELLS];
    let mut channels = [ShearChannelReport::default(); MAX_CHANNELS];
    let mut layers = [0.0; MAX_LAYER_REPORTS];
    let mut current_film = None;
    let mut deposited = 0.0;
    for route in routes.iter().filter(|route| route.sliding_layers.is_some()) {
        let channel = route
            .channel
            .ok_or_else(|| invalid("sliding heat requires an original channel"))?;
        let identity = (route.kind, route.film);
        let range = match route.kind {
            FilmRouteKind::Sliding => {
                let film = system
                    .sliding_film(route.film)
                    .ok_or_else(|| invalid("sliding heat lost its original pressure island"))?;
                let native = film.pressure_operator();
                if current_film != Some(identity) {
                    film.evaluate_with_layers_into(
                        &configuration[..n],
                        velocity,
                        &mut force[..n],
                        &mut pressure[..native.cell_count()],
                        &mut channels[..native.channel_count()],
                        &mut layers[..film.layer_count()],
                    )
                    .map_err(|error| invalid(error.0))?;
                }
                film.channel_layer_range(channel)
            }
            FilmRouteKind::DrivenSliding => {
                let film = system
                    .driven_sliding_film(route.film)
                    .ok_or_else(|| invalid("driven sliding heat lost its original pressure island"))?;
                if current_film != Some(identity) {
                    // Only native viscous layer work is deposited. The same
                    // solve's signed reservoir work already belongs to the
                    // accepted mechanical supply account.
                    film.evaluate_with_layers_into(
                        &configuration[..n],
                        velocity,
                        &mut force[..n],
                        &mut pressure[..film.cell_count()],
                        &mut channels[..film.channel_count()],
                        &mut layers[..film.layer_count()],
                    )
                    .map_err(|error| invalid(error.0))?;
                }
                film.channel_layer_range(channel)
            }
            _ => return Err(invalid("sliding heat lost its original film namespace")),
        }
        .ok_or_else(|| invalid("sliding heat lost its original layer range"))?;
        current_film = Some(identity);
        let destinations = route
            .sliding_layers
            .as_ref()
            .expect("selected sliding route");
        if destinations.len() != range.len() {
            return Err(invalid("sliding layer heat changed its original coverage"));
        }
        let power = channels[channel].dissipated_power_w();
        let heat = system.config.dt_s * power;
        if !heat.is_finite() || heat < 0.0 || (power > 0.0 && heat == 0.0) {
            return Err(invalid("sliding film heat is unrepresentable"));
        }

        // Integrate the actual combined layer powers. Place only the floating-
        // point closure remainder in the largest contributor; this conserves
        // the native channel heat without erasing a tiny positive layer.
        let mut largest = 0;
        for local in 1..range.len() {
            if layers[range.start + local] > layers[range.start + largest] {
                largest = local;
            }
        }
        let mut assigned = 0.0;
        let mut correction = 0.0;
        for (local, &destination) in destinations.iter().enumerate() {
            if local == largest {
                continue;
            }
            let power = layers[range.start + local];
            let amount = system.config.dt_s * power;
            if !amount.is_finite() || amount < 0.0 || (power > 0.0 && amount == 0.0) {
                return Err(invalid("positive sliding layer heat is unrepresentable"));
            }
            heat_j[destination] += amount;
            let increment = amount - correction;
            let next = assigned + increment;
            correction = (next - assigned) - increment;
            assigned = next;
        }
        let remainder = heat - assigned;
        if !remainder.is_finite() || remainder < 0.0 || (heat > 0.0 && remainder == 0.0) {
            return Err(invalid(
                "conservative sliding layer heat is unrepresentable",
            ));
        }
        heat_j[destinations[largest]] += remainder;
        // Count physical work once per channel, independently of layer count.
        deposited += heat;
    }
    Ok(deposited)
}
