//! A source-owned nonlinear film in series with the original material integral.

use super::{ConductivityCurve, DepositionError, LayerConductionConfig, ReservoirLayerGeometry};

mod coupled;
pub use coupled::CoupledLayerResponse;

/// Preserve the original thermal/source failure, or the actual unclosed power.
#[derive(Debug)]
pub enum ConvectiveLayerError<E> {
    Thermal(DepositionError),
    Film(E),
    SurfaceBalance { residual_w: f64, iterations: usize },
}
impl<E> From<DepositionError> for ConvectiveLayerError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ConvectiveLayerError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "conducting convective surface: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ConvectiveLayerError<E> {}

/// One supported surface and its original film evaluation. This stores no heat.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvectiveLayerResponse<F> {
    pub surface_temperature_k: f64,
    /// Signed material integral-coordinate power, positive toward the bath.
    /// Reconstructing it from rounded absolute surface T would lose accuracy
    /// for a thin, highly conducting layer.
    pub solid_power_to_bath_w: f64,
    /// Actual callback film heat. Apply this once to the reciprocal heat owners.
    pub power_to_bath_w: f64,
    /// Actual film Q/(Tcell-Tbulk). At equilibrium this uses the local material
    /// slope in series with the explicitly supplied film G; G=0 stays zero.
    pub conductance_w_k: f64,
    /// Integral-coordinate solid power minus actual film power [W]. This is
    /// neither a rounded-T reconstruction nor a global flux-error certificate.
    pub heat_residual_w: f64,
    /// Interior root trials; at most two additional film calls form the bracket.
    pub iterations: usize,
    pub film: F,
}

struct Film<F> {
    power: f64,
    conductance: f64,
    payload: F,
}

fn evaluate<F, E>(
    surface: f64,
    bulk: f64,
    film: &mut impl FnMut(f64) -> Result<(f64, f64, F), E>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Film<F>, ConvectiveLayerError<E>> {
    checkpoint()?;
    let (power, conductance, payload) = film(surface).map_err(ConvectiveLayerError::Film)?;
    checkpoint()?;
    if !power.is_finite() || !conductance.is_finite() || conductance < 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    let difference = surface - bulk;
    if (difference == 0.0 && power != 0.0)
        || (difference != 0.0
            && (power == 0.0 || conductance == 0.0 || power.signum() != difference.signum()))
    {
        return Err(DepositionError::InvalidInput(
            "convective film must be passive, with positive driven heat and zero equilibrium heat",
        )
        .into());
    }
    Ok(Film {
        power,
        conductance,
        payload,
    })
}

fn closed<E>(
    solid: f64,
    film: f64,
    config: LayerConductionConfig,
) -> Result<bool, ConvectiveLayerError<E>> {
    residual_closed(solid - film, solid.abs().max(film.abs()), config)
}

fn residual_closed<E>(
    residual: f64,
    scale: f64,
    config: LayerConductionConfig,
) -> Result<bool, ConvectiveLayerError<E>> {
    let tolerance = config.absolute_heat_tolerance_w + config.relative_heat_tolerance * scale;
    if !residual.is_finite() || !tolerance.is_finite() {
        return Err(DepositionError::Arithmetic.into());
    }
    Ok(residual.abs() <= tolerance)
}

struct PowerTrial<F> {
    residual_w: f64,
    scale_w: f64,
    admissible: bool,
    state: F,
}

// The one- and two-bath owners use the SAME safeguarded power iteration.
// Their original source callbacks/inverse supply each trial and cancellation
// checks; no source is queried or rounded temperature flux reconstructed here.
#[allow(clippy::too_many_arguments)]
fn power_root<F, E>(
    [mut low, mut high]: [f64; 2],
    [mut low_residual, mut high_residual]: [f64; 2],
    guess: f64,
    orientation: f64,
    config: LayerConductionConfig,
    mut evaluate_trial: impl FnMut(f64) -> Result<PowerTrial<F>, ConvectiveLayerError<E>>,
) -> Result<(f64, F, usize), ConvectiveLayerError<E>> {
    let mut power = if guess.is_finite() && guess > low && guess < high {
        guess
    } else {
        low + 0.5 * (high - low)
    };
    let mut interpolated_width = None;
    let mut residual = orientation * high_residual;
    for iteration in 1..=config.max_iterations {
        let trial = evaluate_trial(power)?;
        residual = trial.residual_w;
        if trial.admissible && residual_closed(residual, trial.scale_w, config)? {
            return Ok((power, trial.state, iteration));
        }
        if orientation * residual > 0.0 {
            high = power;
            high_residual = orientation * residual;
        } else {
            low = power;
            low_residual = orientation * residual;
        }
        let width = high - low;
        let midpoint = low + 0.5 * width;
        let next = if interpolated_width.is_none_or(|previous| width <= 0.5 * previous) {
            let a = low_residual.abs();
            let b = high_residual.abs();
            let weight = if a >= b {
                1.0 / (1.0 + b / a)
            } else {
                let ratio = a / b;
                ratio / (1.0 + ratio)
            };
            let secant = low + width * weight;
            if secant > low && secant < high {
                interpolated_width = Some(width);
                secant
            } else {
                interpolated_width = None;
                midpoint
            }
        } else {
            interpolated_width = None;
            midpoint
        };
        if next <= low || next >= high || next == power {
            return Err(ConvectiveLayerError::SurfaceBalance {
                residual_w: residual,
                iterations: iteration,
            });
        }
        power = next;
    }
    Err(ConvectiveLayerError::SurfaceBalance {
        residual_w: residual,
        iterations: config.max_iterations,
    })
}

fn finish<F, E>(
    surface: f64,
    solid: f64,
    difference: f64,
    film: Film<F>,
    iterations: usize,
) -> Result<ConvectiveLayerResponse<F>, ConvectiveLayerError<E>> {
    let conductance = film.power / difference;
    if !conductance.is_finite() || conductance <= 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    Ok(ConvectiveLayerResponse {
        surface_temperature_k: surface,
        solid_power_to_bath_w: solid,
        power_to_bath_w: film.power,
        conductance_w_k: conductance,
        heat_residual_w: solid - film.power,
        iterations,
        film: film.payload,
    })
}

impl ReservoirLayerGeometry {
    /// Solve an original nonlinear film against the exact supported k(T) path.
    /// Endpoints are [Tcell,Tbulk]. The callback receives actual surface kelvin
    /// and returns (signed surface-to-bulk watts, film G in W/K, source payload).
    /// Its own material/regime/geometry admission remains authoritative. Film G
    /// supplies only an initial series guess and the explicit equilibrium limit.
    /// The caller must retain the same physical bulk state throughout this solve.
    ///
    /// Trial power is nonnegative along the cell-to-bulk direction. Inverting
    /// the original material integral avoids amplifying rounded surface T into
    /// a thin-metal flux error. Only the material trace and solved surface need
    /// conductivity support; the bulk may lie outside it. No extrapolation,
    /// temperature-bracket acceptance, numerical derivative or nested root is
    /// used. At most max_iterations+2 film evaluations occur, each surrounded
    /// by cancellation checks. A safeguarded secant reuses bracket residuals
    /// and forces bisection when interpolation fails to halve the bracket.
    ///
    /// This accepts a supported passive root with the declared power residual;
    /// it does not establish uniqueness or monotonicity of an arbitrary source
    /// film. The layer is homogeneous, pure-phase, source-free and quasi-steady;
    /// the consuming solver still owns storage, time and reciprocal heat.
    ///
    /// A face drop below Kelvin resolution is never floored; the original
    /// inverse refuses it. Only when the first power trial is refused that way
    /// and the cell face itself (zero solid power, its film already evaluated)
    /// closes within the declared power tolerance is that face returned, with
    /// zero iterations and its actual residual.
    pub fn evaluate_convective_material_to_bath<F, E>(
        self,
        curve: &ConductivityCurve,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        mut evaluate_film: impl FnMut(f64) -> Result<(f64, f64, F), E>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConvectiveLayerResponse<F>, ConvectiveLayerError<E>> {
        checkpoint()?;
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(
                DepositionError::InvalidInput("invalid convective layer solve controls").into(),
            );
        }
        let [cell, bulk] = temperatures_k;
        if !cell.is_finite() || cell <= 0.0 || !bulk.is_finite() || bulk <= 0.0 {
            return Err(DepositionError::MaterialDomain.into());
        }
        let conductivity = curve.conductivity_w_m_k(cell)?;
        let difference = cell - bulk;
        let first = evaluate(cell, bulk, &mut evaluate_film, &mut checkpoint)?;
        if difference == 0.0 {
            let conductance = if first.conductance == 0.0 {
                0.0
            } else {
                self.slope(conductivity, 1.0 / first.conductance)?
            };
            checkpoint()?;
            return Ok(ConvectiveLayerResponse {
                surface_temperature_k: cell,
                solid_power_to_bath_w: 0.0,
                power_to_bath_w: first.power,
                conductance_w_k: conductance,
                heat_residual_w: 0.0,
                iterations: 0,
                film: first.payload,
            });
        }
        let sign = difference.signum();
        let span = curve.temperature_range_k();
        // Only the supported search endpoint is clipped, never a source query.
        let limit = bulk.clamp(span[0], span[1]);
        if limit == cell {
            return Err(DepositionError::MaterialDomain.into());
        }
        let maximum = self.multiplier_m * curve.integral_w_m(cell, limit)?.abs();
        if !maximum.is_finite() || maximum <= 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        let low_residual = -sign * first.power;
        let guess = self.slope(conductivity, 1.0 / first.conductance)? * difference.abs();
        // Use the exact source-limit face, avoiding a multiply/divide round trip
        // that could place an otherwise valid boundary query outside the chart.
        let end = evaluate(limit, bulk, &mut evaluate_film, &mut checkpoint)?;
        let high_residual = maximum - sign * end.power;
        if !high_residual.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        if high_residual < 0.0 {
            return Err(DepositionError::MaterialDomain.into());
        }
        if end.power != 0.0 && closed(sign * maximum, end.power, config)? {
            return finish(limit, sign * maximum, difference, end, 0);
        }
        drop(end);
        // The cell face is the power-zero end of this bracket; its film was
        // evaluated above. A few ulp from the bath, the series guess can move
        // the face by less than half an ulp, and the original inverse refuses
        // that unrepresentable drop. Return the cell face only when that first
        // trial (the same start and inverse call as power_root's) would be
        // refused AND the face closes under the declared power tolerance, the
        // residual rule every accepted trial passes. A drop whose omission
        // leaves an unclosed power residual stays refused.
        if guess.is_finite()
            && guess > 0.0
            && guess < maximum
            && closed(0.0, first.power, config)?
            && matches!(
                curve.temperature_after_integral_w_m(cell, -sign * (guess / self.multiplier_m)),
                Err(DepositionError::Arithmetic)
            )
        {
            return finish(cell, 0.0, difference, first, 0);
        }
        drop(first);
        let (power, (surface, film), iterations) = power_root(
            [0.0, maximum],
            [low_residual, high_residual],
            guess,
            sign,
            config,
            |power| {
                checkpoint()?;
                let surface = curve
                    .temperature_after_integral_w_m(cell, -sign * (power / self.multiplier_m))?;
                if !surface.is_finite() || surface < cell.min(limit) || surface > cell.max(limit) {
                    return Err(DepositionError::Arithmetic.into());
                }
                let film = evaluate(surface, bulk, &mut evaluate_film, &mut checkpoint)?;
                Ok(PowerTrial {
                    residual_w: sign * power - film.power,
                    scale_w: power.abs().max(film.power.abs()),
                    admissible: film.power != 0.0,
                    state: (surface, film),
                })
            },
        )?;
        finish(surface, sign * power, difference, film, iterations)
    }
}

#[cfg(test)]
mod tests;
