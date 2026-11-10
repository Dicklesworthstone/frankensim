//! Evolving source-backed droplets over one accuracy-controlled time window.
//!
//! The original material controller compares one coarse interval with two
//! conservative fine intervals. This adapter reruns the original sourced
//! chamber at EVERY trial state, refreshing each actual radius, molecular D,
//! selected Sherwood model and forced mean-film properties from that trial's
//! droplets and shared gas. Only accepted fine intervals contribute physical
//! mass/heat accounts; only the last accepted fine source receipts escape.
//!
//! Initial zero-duration admission also runs for a positive window, so neither
//! zero elapsed time nor adaptive rejection can hide an unsupported original
//! source or volatile identity. Only an original phase StepSizeLimit permits
//! physical refinement. Source, correlation, dilute, conservation and
//! cancellation failures retain their original owner and remain terminal.
//!
//! Local coarse/fine discrepancies are numerical estimates, not certified
//! global error bounds. Each original interval still freezes its initial film
//! and geometry, while subsequent intervals rebuild them. Gas volume remains
//! fixed; original positive-liquid support, no dry-out/nucleation, and no
//! droplet motion, drag, shielding or breakup assumptions remain unchanged.

use super::{
    ClosedGasParcel, GasFilmError, GasFilmHeatSolve, GasLiquidChamberHeatSolve,
    SourceDropletFilmStep, SourceGasDropletHeatSurface, SourceGasDropletSurface,
    SphericalDropletFilmError, ThermochemicalGasModel, advance_droplets_with_gas as advance_source,
    advance_droplets_with_gas_and_heat as advance_source_with_heat,
};
use fs_material::evaporation::gas::{
    GasEvaporationError,
    droplet::chamber::{
        GasDropletChamberCandidate,
        adaptive::{
            AdaptiveGasDropletCandidate, AdaptiveGasDropletError, AdaptiveGasDropletReport,
            GasDropletSeed, advance_droplets_with_gas_adaptive_using,
        },
    },
};

pub use fs_material::evaporation::gas::droplet::chamber::adaptive::AdaptiveGasDropletControls;

/// Original source refusals retain their precise error and attempted interval.
pub type SourceAdaptiveGasDropletError = AdaptiveGasDropletError<SphericalDropletFilmError>;

/// Accepted physical trajectory, numerical accounts and final source evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceAdaptiveGasDropletCandidate<'catalog> {
    original: AdaptiveGasDropletCandidate<'catalog, Vec<SourceDropletFilmStep>>,
}
impl<'catalog> SourceAdaptiveGasDropletCandidate<'catalog> {
    pub const fn original(
        &self,
    ) -> &AdaptiveGasDropletCandidate<'catalog, Vec<SourceDropletFilmStep>> {
        &self.original
    }
    /// Last accepted fine interval's used and refreshed source samples.
    /// A zero-duration request retains its original admission samples instead.
    pub fn films(&self) -> &[SourceDropletFilmStep] {
        self.original.payload()
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletChamberCandidate<'catalog>,
        Vec<SourceDropletFilmStep>,
        AdaptiveGasDropletReport,
    ) {
        self.original.into_parts()
    }
}

/// Advance changing droplet geometry and sourced films with shared-gas feedback.
/// Source claims and explicit transport selections stay fixed throughout the
/// request; actual coefficients refresh on every original trial interval.
pub fn advance_droplets_with_gas_adaptive<'catalog>(
    surfaces: &[SourceGasDropletSurface<'_, '_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceAdaptiveGasDropletCandidate<'catalog>, SourceAdaptiveGasDropletError> {
    advance(
        surfaces.iter().copied().map(|surface| (surface, None)),
        gas,
        catalog,
        duration_s,
        controls,
        None,
        cancelled,
    )
}

/// The same adaptive trajectory with original reciprocal collective heat.
/// Explicit heat coefficients remain independent of phase-carried enthalpy;
/// each trial's current area supplies its physical heat-transfer area.
pub fn advance_droplets_with_gas_and_heat_adaptive<'catalog>(
    surfaces: &[SourceGasDropletHeatSurface<'_, '_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceAdaptiveGasDropletCandidate<'catalog>, SourceAdaptiveGasDropletError> {
    advance(
        surfaces.iter().map(|surface| {
            (
                surface.surface,
                Some((
                    surface.heat_transfer_coefficient_w_m2_k,
                    surface.heat_controls,
                )),
            )
        }),
        gas,
        catalog,
        duration_s,
        controls,
        Some(heat_controls),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'surface, 'm: 'surface, 'catalog>(
    surfaces: impl ExactSizeIterator<
        Item = (
            SourceGasDropletSurface<'surface, 'm>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceAdaptiveGasDropletCandidate<'catalog>, SourceAdaptiveGasDropletError> {
    // Bound storage before the one canonical source-roster allocation. The
    // original controller retains all other numerical and identity admission.
    if controls.solve.maximum_surfaces == 0 || surfaces.len() > controls.solve.maximum_surfaces {
        return Err(AdaptiveGasDropletError::Invalid(
            "sourced adaptive droplet count exceeds its positive work ceiling",
        ));
    }
    let mut ordered: Vec<_> = surfaces.collect();
    ordered.sort_unstable_by_key(|(surface, _)| surface.id);
    let seeds: Vec<_> = ordered
        .iter()
        .map(|(surface, _)| GasDropletSeed {
            id: surface.id,
            law: surface.law,
            binding: surface.binding,
            droplet: surface.droplet,
            temperature_from_enthalpy: surface.temperature_from_enthalpy,
        })
        .collect();
    let original = advance_droplets_with_gas_adaptive_using(
        &seeds,
        gas,
        catalog,
        duration_s,
        controls,
        |current, gas, duration, check| {
            if current.len() != ordered.len()
                || current
                    .iter()
                    .zip(&ordered)
                    .any(|(seed, (surface, _))| seed.id != surface.id)
            {
                return Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(
                    "original adaptive droplet controller changed the source roster",
                )));
            }
            let next = if let Some(heat_controls) = heat_controls {
                let rows: Vec<_> = current
                    .iter()
                    .zip(&ordered)
                    .map(|(seed, (surface, heat))| {
                        let (coefficient, heat_controls) = heat.ok_or(GasFilmError::Invalid(
                            "sourced adaptive droplet heat row is missing",
                        ))?;
                        Ok(SourceGasDropletHeatSurface {
                            surface: SourceGasDropletSurface {
                                droplet: seed.droplet,
                                ..*surface
                            },
                            heat_transfer_coefficient_w_m2_k: coefficient,
                            heat_controls,
                        })
                    })
                    .collect::<Result<_, GasFilmError>>()?;
                advance_source_with_heat(
                    &rows,
                    gas,
                    catalog,
                    duration,
                    controls.solve,
                    heat_controls,
                    check,
                )?
            } else {
                let rows: Vec<_> = current
                    .iter()
                    .zip(&ordered)
                    .map(|(seed, (surface, _))| SourceGasDropletSurface {
                        droplet: seed.droplet,
                        ..*surface
                    })
                    .collect();
                advance_source(&rows, gas, catalog, duration, controls.solve, check)?
            };
            Ok(next.into_parts())
        },
        |error| match error {
            SphericalDropletFilmError::Phase(GasEvaporationError::StepSizeLimit {
                maximum_duration_s,
                ..
            }) => Some(*maximum_duration_s),
            _ => None,
        },
        cancelled,
    )?;
    Ok(SourceAdaptiveGasDropletCandidate { original })
}
