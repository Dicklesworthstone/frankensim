//! Finite retarded observations integrated on a detector's proper clock.
//!
//! The observable supplies nonnegative rates per receiver proper second at the
//! original null ray. Four-point Gauss--Legendre quadrature integrates each
//! authored proper-time panel. The opt-in motion-knot path also splits at
//! detector turns and the received source turns. Material-history changes and
//! smooth-panel refinement remain the caller's accuracy controls.
//! This is a deterministic numerical estimate, not an interval error bound or
//! a guarantee about unsampled source behavior, visibility or coincidences.

use super::{
    PiecewiseVacuumWorldline, RetardedVacuumError, RetardedVacuumRay, UniformVacuumWorldline,
    checked, clock,
};
use fs_math::dd::Dd;

pub const MAX_EXPOSURE_CHANNELS: usize = 32;

/// Nodes and weights on [0,1], shared with the thermal spectral integration.
pub(crate) const GAUSS4: [(f64, f64); 4] = [
    (0.069_431_844_202_973_71, 0.173_927_422_568_726_92),
    (0.330_009_478_207_571_87, 0.326_072_577_431_273_05),
    (0.669_990_521_792_428_1, 0.326_072_577_431_273_05),
    (0.930_568_155_797_026_2, 0.173_927_422_568_726_92),
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedExposureSettings {
    /// Includes history endpoints and four nodes per panel. Automatic alignment
    /// also charges its two planning endpoints and each source-knot arrival.
    pub max_ray_queries: usize,
    /// Original absolute geometry/clock admission tolerance in metres.
    pub max_residual_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RetardedExposureError<E> {
    InvalidPanels,
    InvalidChannels,
    RayBudget,
    InvalidRate,
    NumericalRange,
    Cancelled,
    Geometry(RetardedVacuumError),
    Observable(E),
}

impl<E: core::fmt::Debug> core::fmt::Display for RetardedExposureError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "retarded exposure: {self:?}")
    }
}
impl<E: core::fmt::Debug> core::error::Error for RetardedExposureError<E> {}

/// Rate integrals, in each original observable's unit times proper seconds.
/// No extra frequency, time-dilation, solid-angle or area factor is applied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedExposureIntegral<const CHANNELS: usize> {
    pub integrated: [f64; CHANNELS],
    pub receiver_proper_time_s: [Dd; 2],
    pub ray_queries: usize,
    pub observable_evaluations: usize,
}

/// Integrate an actual source/receiver pair without allocation or time stepping.
///
/// Both endpoint rays are admitted before the first observable callback.
/// Subluminal endpoint worldlines give a monotone emission clock, so this also
/// admits the complete interval against their finite authored time histories.
/// A callback may use the supplied checkpoint to bound its own work. It must
/// read source history without advancing it: failure publishes no integral,
/// but cannot undo side effects performed by a caller's callback.
pub fn integrate_uniform<const C: usize, E>(
    source: &UniformVacuumWorldline,
    receiver: &UniformVacuumWorldline,
    panel_edges_s: &[Dd],
    settings: RetardedExposureSettings,
    at_ray: impl FnMut(&RetardedVacuumRay, &mut dyn FnMut() -> bool) -> Result<[f64; C], E>,
    cancelled: impl FnMut() -> bool,
) -> Result<RetardedExposureIntegral<C>, RetardedExposureError<E>> {
    integrate(
        panel_edges_s,
        settings,
        |proper| {
            source.observe(
                receiver.reception_at_proper_time(proper, settings.max_residual_m)?,
                settings.max_residual_m,
            )
        },
        at_ray,
        cancelled,
    )
}

/// The same proper-clock exposure over finite piecewise-inertial histories.
/// Authored panel edges, rather than hidden adaptive subdivisions, control the
/// temporal approximation across both endpoint histories' velocity changes.
pub fn integrate_piecewise<const C: usize, E>(
    source: &PiecewiseVacuumWorldline,
    receiver: &PiecewiseVacuumWorldline,
    panel_edges_s: &[Dd],
    settings: RetardedExposureSettings,
    at_ray: impl FnMut(&RetardedVacuumRay, &mut dyn FnMut() -> bool) -> Result<[f64; C], E>,
    cancelled: impl FnMut() -> bool,
) -> Result<RetardedExposureIntegral<C>, RetardedExposureError<E>> {
    integrate(
        panel_edges_s,
        settings,
        |proper| {
            source.observe(
                receiver.reception_at_proper_time(proper, settings.max_residual_m)?,
                settings.max_residual_m,
            )
        },
        at_ray,
        cancelled,
    )
}

/// Add the physical motion boundaries to caller-authored proper-time panels.
///
/// Receiver knots split at their accumulated proper times; source knots split
/// when their emitted light actually reaches this moving receiver. Only source
/// knots strictly inside the admitted endpoint emission interval are queried.
/// This resolves velocity jumps, not temperature, visibility or other material
/// discontinuities. Those still require caller-authored edges and refinement.
///
/// One fallible edge allocation is bounded by `max_ray_queries`. Ordered streams
/// are merged without sorting; exact cross-stream duplicates share one edge.
/// Planning, history, allocation and total-ray-budget admission finish before
/// any observable callback. The original authored-only API remains available.
pub fn integrate_piecewise_with_knots<const C: usize, E>(
    source: &PiecewiseVacuumWorldline,
    receiver: &PiecewiseVacuumWorldline,
    panel_edges_s: &[Dd],
    settings: RetardedExposureSettings,
    at_ray: impl FnMut(&RetardedVacuumRay, &mut dyn FnMut() -> bool) -> Result<[f64; C], E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RetardedExposureIntegral<C>, RetardedExposureError<E>> {
    let authored_queries =
        admit_panels::<C, E>(panel_edges_s, settings.max_ray_queries, &mut cancelled)?;
    if authored_queries
        .checked_add(2)
        .is_none_or(|n| n > settings.max_ray_queries)
    {
        return Err(RetardedExposureError::RayBudget);
    }
    let interval = [
        clock(panel_edges_s[0]),
        clock(panel_edges_s[panel_edges_s.len() - 1]),
    ];
    let interval = [
        interval[0].map_err(RetardedExposureError::Geometry)?,
        interval[1].map_err(RetardedExposureError::Geometry)?,
    ];
    let mut endpoint = |proper| {
        if cancelled() {
            return Err(RetardedExposureError::Cancelled);
        }
        source
            .observe(
                receiver
                    .reception_at_proper_time(proper, settings.max_residual_m)
                    .map_err(RetardedExposureError::Geometry)?,
                settings.max_residual_m,
            )
            .map_err(RetardedExposureError::Geometry)
    };
    let emission_interval = [
        endpoint(interval[0])?.emission_coordinate_time_s(),
        endpoint(interval[1])?.emission_coordinate_time_s(),
    ];
    if !emission_interval[0].lt(emission_interval[1]) {
        return Err(RetardedExposureError::NumericalRange);
    }
    let source_range = interior_knots(
        source.segment_count(),
        emission_interval,
        |i| Ok(source.knot(i).0),
        &mut cancelled,
    )?;
    let receiver_range = interior_knots(
        receiver.segment_count(),
        interval,
        |i| receiver.proper_knot_time(i),
        &mut cancelled,
    )?;
    let planning_queries = source_range
        .len()
        .checked_add(2)
        .ok_or(RetardedExposureError::RayBudget)?;
    let budget = settings
        .max_ray_queries
        .checked_sub(planning_queries)
        .filter(|&n| n >= authored_queries)
        .ok_or(RetardedExposureError::RayBudget)?;
    let edge_cap = (budget - 2) / GAUSS4.len() + 1;
    let worst_edges = panel_edges_s
        .len()
        .checked_add(source_range.len())
        .and_then(|n| n.checked_add(receiver_range.len()))
        .ok_or(RetardedExposureError::RayBudget)?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(worst_edges.min(edge_cap))
        .map_err(|_| RetardedExposureError::Geometry(RetardedVacuumError::AllocationRefused))?;
    let (mut authored, mut source_index, mut receiver_index) =
        (0, source_range.start, receiver_range.start);
    let mut source_edge = None;
    let mut previous_source_edge = interval[0];
    loop {
        if cancelled() {
            return Err(RetardedExposureError::Cancelled);
        }
        if source_edge.is_none() && source_index < source_range.end {
            let (time, position) = source.knot(source_index);
            let arrival = receiver
                .receive_emission(time, position, settings.max_residual_m)
                .map_err(RetardedExposureError::Geometry)?;
            let proper = arrival
                .proper_time_s()
                .ok_or(RetardedExposureError::NumericalRange)?;
            if !previous_source_edge.lt(proper) || !proper.lt(interval[1]) {
                return Err(RetardedExposureError::NumericalRange);
            }
            previous_source_edge = proper;
            source_edge = Some(proper);
        }
        let receiver_edge = if receiver_index < receiver_range.end {
            let proper = receiver
                .proper_knot_time(receiver_index)
                .map_err(RetardedExposureError::Geometry)?;
            if !receiver
                .proper_knot_time(receiver_index - 1)
                .map_err(RetardedExposureError::Geometry)?
                .lt(proper)
                || !proper.lt(receiver
                    .proper_knot_time(receiver_index + 1)
                    .map_err(RetardedExposureError::Geometry)?)
            {
                return Err(RetardedExposureError::NumericalRange);
            }
            Some(proper)
        } else {
            None
        };
        let authored_edge = panel_edges_s
            .get(authored)
            .copied()
            .map(clock)
            .transpose()
            .map_err(RetardedExposureError::Geometry)?;
        let mut next = authored_edge;
        for candidate in [source_edge, receiver_edge].into_iter().flatten() {
            if next.is_none_or(|t| candidate.lt(t)) {
                next = Some(candidate);
            }
        }
        let Some(next) = next else {
            break;
        };
        if edges.last() != Some(&next) {
            if edges.len() == edge_cap {
                return Err(RetardedExposureError::RayBudget);
            }
            edges.push(next);
        }
        if authored_edge == Some(next) {
            authored += 1;
        }
        if receiver_edge == Some(next) {
            receiver_index += 1;
        }
        if source_edge == Some(next) {
            source_index += 1;
            source_edge = None;
        }
    }
    let mut result = integrate_piecewise(
        source,
        receiver,
        &edges,
        RetardedExposureSettings {
            max_ray_queries: budget,
            ..settings
        },
        at_ray,
        cancelled,
    )?;
    result.ray_queries += planning_queries;
    Ok(result)
}

/// Binary bounds exclude endpoint knots before any future-cone query is made.
fn interior_knots<E>(
    count: usize,
    interval: [Dd; 2],
    mut time_at: impl FnMut(usize) -> Result<Dd, RetardedVacuumError>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<core::ops::Range<usize>, RetardedExposureError<E>> {
    let mut bounds = [1, count];
    for (side, edge) in interval.into_iter().enumerate() {
        let (mut lower, mut upper) = (1, count);
        while lower < upper {
            if cancelled() {
                return Err(RetardedExposureError::Cancelled);
            }
            let middle = lower + (upper - lower) / 2;
            let time = time_at(middle).map_err(RetardedExposureError::Geometry)?;
            if time.lt(edge) || (side == 0 && time == edge) {
                lower = middle + 1;
            } else {
                upper = middle;
            }
        }
        bounds[side] = lower;
    }
    Ok(bounds[0]..bounds[1])
}

fn integrate<const C: usize, E>(
    edges: &[Dd],
    settings: RetardedExposureSettings,
    mut ray_at: impl FnMut(Dd) -> Result<RetardedVacuumRay, RetardedVacuumError>,
    mut at_ray: impl FnMut(&RetardedVacuumRay, &mut dyn FnMut() -> bool) -> Result<[f64; C], E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RetardedExposureIntegral<C>, RetardedExposureError<E>> {
    let queries = admit_panels::<C, E>(edges, settings.max_ray_queries, &mut cancelled)?;
    let endpoints = [clock(edges[0]), clock(edges[edges.len() - 1])];
    let endpoints = [
        endpoints[0].map_err(RetardedExposureError::Geometry)?,
        endpoints[1].map_err(RetardedExposureError::Geometry)?,
    ];
    for endpoint in endpoints {
        if cancelled() {
            return Err(RetardedExposureError::Cancelled);
        }
        ray_at(endpoint).map_err(RetardedExposureError::Geometry)?;
    }
    let mut totals = [Dd::ZERO; C];
    for pair in edges.windows(2) {
        let lower = clock(pair[0]).map_err(RetardedExposureError::Geometry)?;
        let upper = clock(pair[1]).map_err(RetardedExposureError::Geometry)?;
        let width = upper - lower;
        let mut previous = lower;
        for (fraction, weight) in GAUSS4 {
            if cancelled() {
                return Err(RetardedExposureError::Cancelled);
            }
            let proper = checked(lower + width * Dd::from_f64(fraction))
                .map_err(RetardedExposureError::Geometry)?;
            if !previous.lt(proper) || !proper.lt(upper) {
                return Err(RetardedExposureError::NumericalRange);
            }
            previous = proper;
            let ray = ray_at(proper).map_err(RetardedExposureError::Geometry)?;
            let rates = at_ray(&ray, &mut cancelled).map_err(RetardedExposureError::Observable)?;
            for (total, rate) in totals.iter_mut().zip(rates) {
                if !rate.is_finite() || rate < 0.0 {
                    return Err(RetardedExposureError::InvalidRate);
                }
                let term = Dd::from_f64(rate) * (width * Dd::from_f64(weight));
                if rate > 0.0 && term == Dd::ZERO {
                    return Err(RetardedExposureError::NumericalRange);
                }
                *total =
                    checked(*total + term).map_err(|_| RetardedExposureError::NumericalRange)?;
            }
        }
    }
    if cancelled() {
        return Err(RetardedExposureError::Cancelled);
    }
    let integrated = totals.map(Dd::to_f64);
    if integrated
        .iter()
        .zip(totals)
        .any(|(&x, dd)| !x.is_finite() || x < 0.0 || (x == 0.0 && dd != Dd::ZERO))
    {
        return Err(RetardedExposureError::NumericalRange);
    }
    Ok(RetardedExposureIntegral {
        integrated,
        receiver_proper_time_s: endpoints,
        ray_queries: queries,
        observable_evaluations: queries - 2,
    })
}

fn admit_panels<const C: usize, E>(
    edges: &[Dd],
    budget: usize,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<usize, RetardedExposureError<E>> {
    if C == 0 || C > MAX_EXPOSURE_CHANNELS {
        return Err(RetardedExposureError::InvalidChannels);
    }
    let queries = edges
        .len()
        .checked_sub(1)
        .filter(|&n| n > 0)
        .and_then(|n| n.checked_mul(GAUSS4.len()))
        .and_then(|n| n.checked_add(2))
        .ok_or(RetardedExposureError::InvalidPanels)?;
    if queries > budget {
        return Err(RetardedExposureError::RayBudget);
    }
    for pair in edges.windows(2) {
        if cancelled() {
            return Err(RetardedExposureError::Cancelled);
        }
        let lower = clock(pair[0]).map_err(RetardedExposureError::Geometry)?;
        let upper = clock(pair[1]).map_err(RetardedExposureError::Geometry)?;
        if !lower.lt(upper) {
            return Err(RetardedExposureError::InvalidPanels);
        }
        checked(upper - lower).map_err(RetardedExposureError::Geometry)?;
    }
    Ok(queries)
}

#[cfg(test)]
mod tests;
