//! Resolve gas-film transport on every original physical deposition interval.
//!
//! A source resolver sees the current finite gas and both condensed states.
//! Growing owners call it after rebinding the mass-owned exposed area. The
//! existing adaptive controller, molecular streams and gas reconstruction stay
//! authoritative; this module supplies no second evolution or storage model.

use super::advance::{
    GasResistiveAdvanceCandidate, GasResistiveAdvanceError, GasResistiveAdvanceSolve,
};
use super::driven::GasResistivePower;
use super::{
    ClosedGasParcel, GasResistiveSurface, SublimationError, ThermalAccommodation,
    ThermochemicalGasModel,
};
use std::collections::BTreeMap;

/// Preserve the original source error independently of the evolution error.
#[derive(Clone, Debug, PartialEq)]
pub enum GasResistiveFilmError<E> {
    Source(E),
    Advance(GasResistiveAdvanceError),
}
impl<E: core::fmt::Debug> core::fmt::Display for GasResistiveFilmError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved resistive gas film: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for GasResistiveFilmError<E> {}

/// Internal interval seam shared by fixed and source-resolved film inputs.
/// The returned coefficients are borrowed only during the original kernel.
pub(super) trait TrialFilms {
    fn resolve(
        &mut self,
        surfaces: &[GasResistiveSurface<'_>],
        gas: &ClosedGasParcel,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Option<&BTreeMap<u64, f64>>, GasResistiveAdvanceError>;
}

impl TrialFilms for Option<&BTreeMap<u64, f64>> {
    fn resolve(
        &mut self,
        _surfaces: &[GasResistiveSurface<'_>],
        _gas: &ClosedGasParcel,
        _cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Option<&BTreeMap<u64, f64>>, GasResistiveAdvanceError> {
        Ok(*self)
    }
}

struct ResolvedFilms<'a, F, E> {
    resolver: &'a mut F,
    coefficients: BTreeMap<u64, f64>,
    failure: Option<E>,
}
impl<F, E> TrialFilms for ResolvedFilms<'_, F, E>
where
    F: FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        &mut BTreeMap<u64, f64>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(), E>,
{
    fn resolve(
        &mut self,
        surfaces: &[GasResistiveSurface<'_>],
        gas: &ClosedGasParcel,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Option<&BTreeMap<u64, f64>>, GasResistiveAdvanceError> {
        if cancelled() {
            return Err(SublimationError::Cancelled.into());
        }
        // Keep stable map nodes across trials, but never silently reuse a
        // coefficient from an earlier gas/geometry state. A sparse resolver
        // must overwrite each retained key or explicitly remove that binding.
        for coefficient in self.coefficients.values_mut() {
            if cancelled() {
                return Err(SublimationError::Cancelled.into());
            }
            *coefficient = f64::NAN;
        }
        if let Err(error) =
            (self.resolver)(surfaces, gas, &mut self.coefficients, cancelled)
        {
            self.failure = Some(error);
            // This non-retryable bridge error exits the original controller.
            // with_resolved_films restores the complete typed source error.
            return Err(GasResistiveAdvanceError::Invalid(
                "original gas-film source resolution refused",
            ));
        }
        if cancelled() {
            return Err(SublimationError::Cancelled.into());
        }
        if self.coefficients.len() > surfaces.len() {
            return Err(GasResistiveAdvanceError::Invalid(
                "resolved gas-film roster exceeds the surface work bound",
            ));
        }
        // The original chamber validates IDs and finite nonnegative values,
        // including any poisoned coefficient the resolver failed to refresh.
        Ok(Some(&self.coefficients))
    }
}

pub(super) fn with_resolved_films<F, E, T>(
    resolver: &mut F,
    run: impl FnOnce(&mut dyn TrialFilms) -> Result<T, GasResistiveAdvanceError>,
) -> Result<T, GasResistiveFilmError<E>>
where
    F: FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        &mut BTreeMap<u64, f64>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(), E>,
{
    let mut films = ResolvedFilms {
        resolver,
        coefficients: BTreeMap::new(),
        failure: None,
    };
    let result = run(&mut films);
    match films.failure {
        Some(error) => Err(GasResistiveFilmError::Source(error)),
        None => result.map_err(GasResistiveFilmError::Advance),
    }
}

/// Evolve with film coefficients derived from each current physical state.
///
/// The resolver receives canonical surface IDs, the current original gas and
/// a reusable output map. It supplies mass-transfer coefficients [m/s] for
/// selected IDs; omitted IDs retain free-molecular transfer. Existing output
/// values are invalidated before each call. Overwrite all retained entries,
/// or remove a key to return that surface to free-molecular transfer. The
/// original chamber refuses unknown IDs, negative/nonfinite coefficients and
/// incompletely refreshed outputs.
///
/// Source resolution runs for every coarse and half interval, including the
/// zero-duration admission. It must be deterministic for its physical inputs,
/// bounded and cancellation-aware; it must not commit caller state on a trial.
/// Source errors retain their original type and abort the whole history.
/// Prescribed powers and nonlinear contacts retain the driven owner's meaning.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_adaptive_resolved<'catalog, F, E>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    mut resolver: F,
    powers: &BTreeMap<u64, GasResistivePower<'_>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveFilmError<E>>
where
    F: FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        &mut BTreeMap<u64, f64>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(), E>,
{
    with_resolved_films(&mut resolver, |films| {
        super::driven::advance_with_films(
            surfaces,
            accommodation,
            films,
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    })
}

#[cfg(test)]
mod tests;
