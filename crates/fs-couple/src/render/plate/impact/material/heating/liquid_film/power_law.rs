//! Explicit separable temperature shifts of the original nonlinear film law.
//!
//! Native geometry, flow index and physical Newtonian-core shear stay fixed.
//! A caller-supplied reference consistency is scaled by the original pinned
//! liquid viscosity's dimensionless temperature ratio. The liquid source does
//! not independently validate that rheological shape or its reference K.

use std::ops::Range;

use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_math::det;
use fs_tribo::resistive_film::power_law::PowerLawFilm;

use super::{HeatedImpactSystem, ImpactError, invalid};

/// One original power-law channel equilibrated to a finite liquid store.
///
/// Explicitly selects the separable thermo-rheological approximation
/// `K(T) = K_ref * mu_source(T) / mu_source(T_ref)`. `T_ref` and the viscosity
/// claims belong to `law.reference_state()`; `K_ref` is supplied independently
/// in Pa s^n. The original native exponent n and physical core shear gamma0
/// stay fixed, so the SAME ratio scales viscosity at every fixed shear rate.
/// Reference liquid viscosity is not reinterpreted as consistency, and source
/// viscosity data alone do not establish temperature-independent n or gamma0.
/// The caller must justify this approximation over the admitted temperature
/// and shear range; it is not a general temperature-dependent rheology model.
///
/// Bind every original channel of a selected island exactly once. Distinct
/// passages may retain distinct source laws and local reservoir temperatures;
/// sharing a reservoir does not average their constitutive coefficients.
/// Temperature bindings never deposit heat: explicit power-law heat routes
/// continue to own the original Q*delta_p work and its selected destinations.
#[derive(Clone, Debug)]
pub struct LiquidPowerLawFilmChannelBinding<'m> {
    /// Original POWER-LAW island index, independent of normal/sliding indices.
    pub film: usize,
    /// Original hydraulic channel index, including currently closed channels.
    pub channel: usize,
    /// Existing finite liquid reservoir setting this passage's temperature.
    pub reservoir: String,
    /// Pinned source law whose viscosity ratio shifts the reference consistency.
    pub law: LiquidTemperatureLaw<'m>,
    /// Positive finite consistency K_ref [Pa s^n] at `law`'s reference T.
    /// Retained across accepted steps and extraction/rebinding; never replaced
    /// by the most recently shifted consistency.
    pub reference_consistency_pa_s_n: f64,
    /// Maximum `abs(mu_source_next-mu_source_old)/mu_source_old` per sample,
    /// in (0,1]. This also bounds relative K and viscosity at fixed shear. It
    /// requests refinement when exceeded; no clipping or error estimate.
    pub maximum_relative_viscosity_change: f64,
}

pub(super) struct BindingGroup {
    pub(super) film: usize,
    rows: Range<usize>,
    shifts: Vec<ConsistencyShift>,
}

struct ConsistencyShift {
    reference_consistency: f64,
    reference_viscosity: f64,
}

impl ConsistencyShift {
    fn at_viscosity(&self, viscosity: f64) -> Result<f64, ImpactError> {
        // Preserve the exact source-reference image. Ordinary magnitudes use
        // one ratio; logarithms only rescue representable combined products
        // whose intermediate ratio overflowed or underflowed.
        let direct = self.reference_consistency * (viscosity / self.reference_viscosity);
        let consistency = if direct.is_finite() && direct > 0.0 {
            direct
        } else {
            det::exp(
                det::ln(self.reference_consistency) + det::ln(viscosity)
                    - det::ln(self.reference_viscosity),
            )
        };
        if !consistency.is_finite() || consistency <= 0.0 {
            return Err(invalid(
                "temperature-shifted film consistency is unrepresentable",
            ));
        }
        Ok(consistency)
    }
}

impl BindingGroup {
    pub(super) fn compile(
        &self,
        film: &PowerLawFilm,
        sources: &[ResolvedLiquidState],
    ) -> Result<PowerLawFilm, ImpactError> {
        let mut laws = film.laws().to_vec();
        for ((law, shift), source) in laws
            .iter_mut()
            .zip(&self.shifts)
            .zip(&sources[self.rows.clone()])
        {
            law.consistency_pa_s_n = shift.at_viscosity(source.dynamic_viscosity_pa_s())?;
        }
        if laws == film.laws() {
            return Ok(film.clone());
        }
        // The native constructor alone admits K, its core viscosity/stress,
        // and the original differential pressure operator. No force or heat
        // is approximated by the core-viscosity geometry view.
        PowerLawFilm::new(film.pressure_operator(), laws, film.controls())
            .map_err(|error| invalid(error.0))
    }
}

pub(super) fn groups(
    inner: &HeatedImpactSystem<'_>,
    bindings: &[LiquidPowerLawFilmChannelBinding<'_>],
    offset: usize,
) -> Result<Vec<BindingGroup>, ImpactError> {
    let mut groups = Vec::new();
    let mut start = 0;
    while start < bindings.len() {
        let film_index = bindings[start].film;
        let end = start + bindings[start..].partition_point(|binding| binding.film == film_index);
        let film = inner
            .mechanics()
            .power_law_film(film_index)
            .ok_or_else(|| invalid("liquid feedback names an absent power-law pressure island"))?;
        if end - start != film.laws().len()
            || bindings[start..end]
                .iter()
                .enumerate()
                .any(|(channel, binding)| binding.channel != channel)
        {
            return Err(invalid(
                "power-law feedback requires every original channel exactly once",
            ));
        }
        let shifts = bindings[start..end]
            .iter()
            .map(|binding| {
                let reference_consistency = binding.reference_consistency_pa_s_n;
                let reference_viscosity = binding.law.reference_state().dynamic_viscosity_pa_s();
                if !reference_consistency.is_finite()
                    || reference_consistency <= 0.0
                    || !reference_viscosity.is_finite()
                    || reference_viscosity <= 0.0
                {
                    return Err(invalid(
                        "power-law feedback needs positive finite reference K and viscosity",
                    ));
                }
                Ok(ConsistencyShift {
                    reference_consistency,
                    reference_viscosity,
                })
            })
            .collect::<Result<Vec<_>, ImpactError>>()?;
        groups.push(BindingGroup {
            film: film_index,
            rows: offset + start..offset + end,
            shifts,
        });
        start = end;
    }
    Ok(groups)
}
