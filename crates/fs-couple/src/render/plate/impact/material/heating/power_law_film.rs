//! Original generalized-Newtonian film work becomes finite material heat.
//!
//! The native monotone stress law and pressure solve already supply physical
//! channel Q*delta_p and its separate numerical work defect. Only the former
//! is thermalized here, at the same midpoint configuration and actual
//! discrete-gradient velocity as mechanics. Laws remain frozen; the existing
//! enthalpy/transport owners supply the thermal clock and all publication.

use fs_exec::CancelGate;
use fs_tribo::resistive_film::{FilmChannelReport, FilmError, MAX_CELLS, MAX_CHANNELS};

use super::{
    CaloricReservoir, FilmRoute, FilmRouteKind, ImpactError, ImpactSystem, MAX_IMPACT_MODES,
    PreparedImpactSystem, invalid,
};

/// Finite thermal destinations for one complete power-law island or passage.
///
/// Indices follow `power_law_film_count`, independently of squeeze/sliding
/// indices. A whole island and its channels cannot both be selected; repeated
/// channels also refuse. Unselected work stays in the original mechanical loss
/// account outside this thermal boundary. Wall/fluid heat shares are explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct PowerLawFilmHeatRoute {
    /// Original generalized-Newtonian pressure-island insertion index.
    pub film: usize,
    /// Original passage index, or None for the complete pressure island.
    pub channel: Option<usize>,
    /// First existing finite reservoir key.
    pub to_a: String,
    /// Second existing finite reservoir key; may equal the first.
    pub to_b: String,
    /// Fraction delivered to A, in [0,1]; B receives the remainder.
    pub fraction_to_a: f64,
}

pub(super) fn compile(
    mechanics: &PreparedImpactSystem,
    route: PowerLawFilmHeatRoute,
    stores: &[CaloricReservoir<'_>],
    compiled: &[FilmRoute],
) -> Result<FilmRoute, ImpactError> {
    let film = mechanics
        .power_law_film(route.film)
        .ok_or_else(|| invalid("power-law heat names an absent constitutive island"))?;
    if route
        .channel
        .is_some_and(|channel| channel >= film.pressure_operator().channel_count())
        || compiled.iter().any(|old| {
            old.kind == FilmRouteKind::PowerLaw
                && old.film == route.film
                && (old.channel.is_none() || route.channel.is_none() || old.channel == route.channel)
        })
    {
        return Err(invalid(
            "power-law heat needs original channels without duplicate or whole-island overlap",
        ));
    }
    if !route.fraction_to_a.is_finite() || !(0.0..=1.0).contains(&route.fraction_to_a) {
        return Err(invalid("power-law heat partition must be finite and in [0,1]"));
    }
    let destination = |key: &str| {
        stores
            .iter()
            .position(|store| store.key() == key)
            .ok_or_else(|| invalid("power-law heat names an absent finite reservoir"))
    };
    Ok(FilmRoute {
        kind: FilmRouteKind::PowerLaw,
        film: route.film,
        channel: route.channel,
        a: destination(&route.to_a)?,
        b: destination(&route.to_b)?,
        fraction: route.fraction_to_a,
        layers: None,
        sliding_layers: None,
    })
}

/// Each selected constitutive island is solved once, regardless of route count.
/// The original bounded nonlinear solver polls cancellation during its work.
#[inline(never)]
pub(super) fn stage(
    system: &ImpactSystem,
    midpoint: &[f64],
    velocity: &[f64],
    routes: &[FilmRoute],
    heat_j: &mut [f64],
    gate: &CancelGate,
) -> Result<f64, ImpactError> {
    let n = system.modes;
    let mut configuration = [0.0; MAX_IMPACT_MODES];
    for i in 0..n {
        configuration[i] = midpoint[2 * i];
    }
    let mut force = [0.0; MAX_IMPACT_MODES];
    let mut pressure = [0.0; MAX_CELLS];
    let mut channels = [FilmChannelReport::default(); MAX_CHANNELS];
    let mut current_film = None;
    let mut total_power = 0.0;
    let mut deposited = 0.0;
    for route in routes.iter().filter(|route| route.kind == FilmRouteKind::PowerLaw) {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if current_film != Some(route.film) {
            let film = system
                .power_law_film(route.film)
                .ok_or_else(|| invalid("power-law heat lost its original constitutive island"))?;
            let native = film.pressure_operator();
            let report = film
                .evaluate_with_channels_and_check_into(
                    &configuration[..n],
                    velocity,
                    &mut force[..n],
                    &mut pressure[..native.cell_count()],
                    &mut channels[..native.channel_count()],
                    || {
                        if gate.is_requested() {
                            Err(FilmError("power-law thermal staging cancelled"))
                        } else {
                            Ok(())
                        }
                    },
                )
                .map_err(|error| {
                    if gate.is_requested() {
                        ImpactError::Cancelled
                    } else {
                        invalid(error.0)
                    }
                })?;
            total_power = report.film.dissipated_power_w;
            current_film = Some(route.film);
        }
        let power = route
            .channel
            .map_or(total_power, |channel| channels[channel].dissipated_power_w);
        let amount = system.config.dt_s * power;
        if !amount.is_finite() || amount < 0.0 || (power > 0.0 && amount == 0.0) {
            return Err(invalid("positive power-law film heat is unrepresentable"));
        }
        let first = amount * route.fraction;
        heat_j[route.a] += first;
        heat_j[route.b] += amount - first;
        deposited += amount;
    }
    Ok(deposited)
}
