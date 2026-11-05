//! Two original passive baths on one source-integrated material surface.

use super::*;

/// One zero-storage surface, two distinct bath powers, and frozen secants.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledLayerResponse<F> {
    pub surface_temperature_k: f64,
    /// Signed integral-coordinate heat from the stored material to its surface.
    /// Retained directly; do not reconstruct thin-layer heat from rounded Ts.
    pub solid_power_to_surface_w: f64,
    /// Original callback heat, positive from the surface into each bath.
    /// The material receives the opposite SUM, even when these oppose exactly.
    pub powers_to_baths_w: [f64; 2],
    /// Positive material/bath pair secants after eliminating the shared surface.
    /// These are held-coefficient decay scales, not nonlinear derivatives.
    pub material_bath_conductances_w_k: [f64; 2],
    /// The bath-to-bath path through the SAME surface, including a neutral body.
    pub bath_pair_conductance_w_k: f64,
    /// Integral-coordinate material power minus the sum of actual bath powers.
    pub heat_residual_w: f64,
    /// Interior trials; at most two additional callback calls admit the bracket.
    pub iterations: usize,
    pub film: F,
}

struct Boundary<F> {
    powers: [f64; 2],
    conductances: [f64; 2],
    total: f64,
    scale: f64,
    payload: F,
}

fn boundary<F, E>(
    surface: f64,
    baths: [f64; 2],
    evaluate: &mut impl FnMut(f64) -> Result<([f64; 2], [f64; 2], F), E>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Boundary<F>, ConvectiveLayerError<E>> {
    checkpoint()?;
    let (powers, conductances, payload) = evaluate(surface).map_err(ConvectiveLayerError::Film)?;
    checkpoint()?;
    for ((power, g), bath) in powers.into_iter().zip(conductances).zip(baths) {
        if !power.is_finite() || !g.is_finite() || g < 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        let difference = surface - bath;
        let represented = g * difference;
        if !represented.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        if (g == 0.0 && power != 0.0)
            || (difference == 0.0 && power != 0.0)
            || (g > 0.0
                && difference != 0.0
                && (power == 0.0 || represented == 0.0 || power.signum() != difference.signum()))
            || (power - represented).abs() > 8.0 * f64::EPSILON * power.abs().max(represented.abs())
        {
            return Err(DepositionError::InvalidInput(
                "each surface bath needs its passive heat/secant or explicit zero insulation",
            )
            .into());
        }
    }
    let total = powers[0] + powers[1];
    if !total.is_finite() {
        return Err(DepositionError::Arithmetic.into());
    }
    Ok(Boundary {
        powers,
        conductances,
        total,
        scale: powers[0].abs().max(powers[1].abs()),
        payload,
    })
}

fn finish<F, E>(
    geometry: ReservoirLayerGeometry,
    curve: &ConductivityCurve,
    cell: f64,
    surface: f64,
    solid: f64,
    boundary: Boundary<F>,
    iterations: usize,
) -> Result<CoupledLayerResponse<F>, ConvectiveLayerError<E>> {
    // Use the integral MEAN k, not Q/(Tc-Ts): the latter would amplify the
    // rounded thin-wall temperature drop into an erroneous decay coefficient.
    let mean_k = if cell == surface {
        curve.conductivity_w_m_k(cell)?
    } else {
        curve.integral_w_m(surface, cell)? / (cell - surface)
    };
    let gs = geometry.multiplier_m * mean_k;
    let [g0, g1] = boundary.conductances;
    if !gs.is_finite() || gs <= 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    // Normalize before summation/products; the individual conductances may be
    // representable even when their unscaled sum or product would overflow.
    let scale = gs.max(g0).max(g1);
    let denominator = gs / scale + g0 / scale + g1 / scale;
    let pair = |a: f64, b: f64| -> Result<f64, DepositionError> {
        if a == 0.0 || b == 0.0 {
            return Ok(0.0);
        }
        let value = a.min(b) * ((a.max(b) / scale) / denominator);
        if !value.is_finite() || value <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok(value)
    };
    Ok(CoupledLayerResponse {
        surface_temperature_k: surface,
        solid_power_to_surface_w: solid,
        powers_to_baths_w: boundary.powers,
        material_bath_conductances_w_k: [pair(gs, g0)?, pair(gs, g1)?],
        bath_pair_conductance_w_k: pair(g0, g1)?,
        heat_residual_w: solid - boundary.total,
        iterations,
        film: boundary.payload,
    })
}

impl ReservoirLayerGeometry {
    /// Close one exact supported k(T) integral against TWO original bath laws.
    /// The callback receives actual surface kelvin and returns
    /// `(powers_to_baths_w, conductances_w_k, payload)`, in the same bath order.
    /// Each power is G*(Ts-Tbath), with only floating-point roundoff admitted;
    /// G=Q=0 explicitly insulates that bath. Their sum may vanish while heat
    /// passes between baths. No weighted ambient or second root is introduced.
    ///
    /// `surface_temperature_range_k` is an inclusive ORIGINAL source domain for
    /// the boundary evaluation (e.g. gas film and surface emissivity). It need
    /// not contain the material cell. Only the cell-to-surface path must lie
    /// on the conductivity chart. Bath temperatures themselves need no material
    /// or boundary-source coverage. The physical three-temperature hull is
    /// intersected with these supports; no source query is extrapolated.
    ///
    /// The common safeguarded power root retains signed material heat and
    /// accepts only the declared residual, scaled by the largest individual
    /// material/bath power. Separate callback powers preserve accounting when
    /// drives oppose or the cell and one bath have identical temperatures.
    /// Pair conductances are frozen secants for the caller's capacity rows;
    /// they do not differentiate changing source coefficients. Original source
    /// and cancellation errors propagate. No uniqueness or model-error bound
    /// is claimed for arbitrary nonlinear passive callbacks.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_material_to_two_baths<F, E>(
        self,
        curve: &ConductivityCurve,
        cell_temperature_k: f64,
        bath_temperatures_k: [f64; 2],
        surface_temperature_range_k: [f64; 2],
        config: LayerConductionConfig,
        mut evaluate_boundary: impl FnMut(f64) -> Result<([f64; 2], [f64; 2], F), E>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledLayerResponse<F>, ConvectiveLayerError<E>> {
        checkpoint()?;
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(
                DepositionError::InvalidInput("invalid coupled layer solve controls").into(),
            );
        }
        let cell = cell_temperature_k;
        let baths = bath_temperatures_k;
        let [source_lo, source_hi] = surface_temperature_range_k;
        if !cell.is_finite()
            || cell <= 0.0
            || baths.into_iter().any(|t| !t.is_finite() || t <= 0.0)
            || !source_lo.is_finite()
            || source_lo <= 0.0
            || !source_hi.is_finite()
            || source_hi < source_lo
        {
            return Err(DepositionError::MaterialDomain.into());
        }
        curve.conductivity_w_m_k(cell)?;
        let [material_lo, material_hi] = curve.temperature_range_k();
        let lo = cell
            .min(baths[0])
            .min(baths[1])
            .max(source_lo)
            .max(material_lo);
        let hi = cell
            .max(baths[0])
            .max(baths[1])
            .min(source_hi)
            .min(material_hi);
        if lo > hi {
            return Err(DepositionError::MaterialDomain.into());
        }
        let low = self.multiplier_m * curve.integral_w_m(hi, cell)?;
        let high = self.multiplier_m * curve.integral_w_m(lo, cell)?;
        if !low.is_finite() || !high.is_finite() || !(high - low).is_finite() || low > high {
            return Err(DepositionError::Arithmetic.into());
        }
        let first = boundary(hi, baths, &mut evaluate_boundary, &mut checkpoint)?;
        let low_residual = low - first.total;
        if (low == 0.0 || first.total != 0.0)
            && residual_closed(low_residual, low.abs().max(first.scale), config)?
        {
            return finish(self, curve, cell, hi, low, first, 0);
        }
        if lo == hi || low_residual > 0.0 {
            return Err(DepositionError::MaterialDomain.into());
        }
        drop(first);
        let last = boundary(lo, baths, &mut evaluate_boundary, &mut checkpoint)?;
        let high_residual = high - last.total;
        if (high == 0.0 || last.total != 0.0)
            && residual_closed(high_residual, high.abs().max(last.scale), config)?
        {
            return finish(self, curve, cell, lo, high, last, 0);
        }
        if high_residual < 0.0 {
            return Err(DepositionError::MaterialDomain.into());
        }
        drop(last);
        let (power, (surface, boundary), iterations) = power_root(
            [low, high],
            [low_residual, high_residual],
            0.0,
            1.0,
            config,
            |power| {
                checkpoint()?;
                let surface =
                    curve.temperature_after_integral_w_m(cell, -power / self.multiplier_m)?;
                if !surface.is_finite() || surface < lo || surface > hi {
                    return Err(DepositionError::Arithmetic.into());
                }
                let state = boundary(surface, baths, &mut evaluate_boundary, &mut checkpoint)?;
                Ok(PowerTrial {
                    residual_w: power - state.total,
                    scale_w: power.abs().max(state.scale),
                    admissible: power == 0.0 || state.total != 0.0,
                    state: (surface, state),
                })
            },
        )?;
        finish(self, curve, cell, surface, power, boundary, iterations)
    }
}

#[cfg(test)]
mod tests;
