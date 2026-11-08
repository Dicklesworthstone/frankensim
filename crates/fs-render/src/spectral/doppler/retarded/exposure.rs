//! Finite retarded observations integrated on a detector's proper clock.
//!
//! The observable supplies nonnegative rates per receiver proper second at the
//! original null ray. Four-point Gauss--Legendre quadrature integrates each
//! authored proper-time panel. Panel edges must resolve source discontinuities
//! and retarded velocity changes; refinement is the caller's accuracy control.
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
    /// Includes both history-preflight endpoints and four nodes per panel.
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

fn integrate<const C: usize, E>(
    edges: &[Dd],
    settings: RetardedExposureSettings,
    mut ray_at: impl FnMut(Dd) -> Result<RetardedVacuumRay, RetardedVacuumError>,
    mut at_ray: impl FnMut(&RetardedVacuumRay, &mut dyn FnMut() -> bool) -> Result<[f64; C], E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RetardedExposureIntegral<C>, RetardedExposureError<E>> {
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
    if queries > settings.max_ray_queries {
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

#[cfg(test)]
mod tests;
