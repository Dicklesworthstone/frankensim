//! Integrated molecular streams through coating birth and complete extinction.
//!
//! At frozen surface/gas temperatures and exposed area, the original kinetic
//! coefficients give `dm_v/dt = c - lambda*m_v`. The finite-vapor solution and
//! BOTH integrated streams are analytical; previously exchanged molecules may
//! return during the same interval. Limiting each gross stream to its original
//! inventory would incorrectly prohibit this turnover and prevent film birth.
//!
//! The caller explicitly admits a barrier-free supported coating. At zero film
//! mass, the outgoing stream is the smaller of its saturation capacity and the
//! incoming stream. An undersaturated bare support therefore reemits arriving
//! sticking molecules immediately at the surface temperature. It supplies no
//! volatile substrate mass. A supersaturated interface starts accumulating a
//! film immediately. This is a zero-residence, complete-thermalization boundary
//! closure, not a nucleation, adsorption-coverage or finite-residence-time law.
//!
//! A first extinction event is located before integrating the subsequent bare
//! interval. Its equal gross streams still transport heat when temperatures
//! differ. The energy owner must use the original molecular caloric law for
//! each stream, debit/credit both bodies, and check the reported mass residuals.
//! No heat, phase datum, morphology or gas energy evolution is inferred here.
//! The optional film retains the original dilute series-resistance semantics.
//! Scalar work and storage are bounded independently of duration; no allocation
//! or nonlinear iteration is performed. This is exact for the declared frozen rates,
//! with floating-point residuals, not an exact nonisothermal time integrator.

use crate::sublimation::{
    MOLAR_GAS_CONSTANT, SolidSurfaceState, SublimationError, SublimationLaw, SurfaceKinetics,
    VaporReservoir, coefficients, finite, nonnegative, positive, properties,
};

/// One simultaneous stream/inventory candidate in the original mass units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularTurnover {
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    /// Positive into vapor. Evaluated independently of subtracting large gross streams.
    pub transferred_mass_kg: f64,
    pub deposit_mass_kg: f64,
    pub reservoir: VaporReservoir,
    /// First complete coating loss within the interval, measured from its start.
    /// Zero means an initially bare surface with a positive evaporation demand.
    /// A sealed interface or an initially bare equilibrium has no extinction event.
    pub extinction_time_s: Option<f64>,
    /// `(outgoing-incoming)-transferred`; consumers own an absolute admission budget.
    pub stream_mass_residual_kg: f64,
    /// `(deposit_final-deposit_initial)+transferred`, including representation loss.
    pub deposit_mass_residual_kg: f64,
    /// `(vapor_final-vapor_initial)-transferred`; zero for a prescribed reservoir.
    pub vapor_mass_residual_kg: f64,
}

fn invalid(message: &'static str) -> SublimationError {
    SublimationError::Invalid(message)
}

// Capture a represented-inventory/stream defect even when its small first
// term disappears in the following subtraction of a much larger datum.
fn balance(a: f64, b: f64, c: f64) -> Result<f64, SublimationError> {
    let (ab, low_ab) = fs_math::eft::two_sum(a, b);
    let (sum, low_c) = fs_math::eft::two_sum(ab, c);
    finite(sum + (low_ab + low_c))
}

// (z-1+exp(-z))/z^2. The first ten terms avoid subtracting the two nearly
// equal O(z) terms. At z<=1/8 the omitted term is below 2e-18 absolutely.
fn second_exponential_ratio(z: f64) -> f64 {
    0.5 + z
        * (-1.0 / 6.0
            + z * (1.0 / 24.0
                + z * (-1.0 / 120.0
                    + z * (1.0 / 720.0
                        + z * (-1.0 / 5040.0
                            + z * (1.0 / 40320.0
                                + z * (-1.0 / 362880.0
                                    + z * (1.0 / 3628800.0 - z / 39916800.0))))))))
}

// -ln(1-x)/x. The positive series keeps the linear-extinction limit when
// 1-x rounds to one. Its first 18 terms have <4e-18 truncation at x<=1/8.
fn extinction_ratio(x: f64) -> f64 {
    if x <= 0.125 {
        let mut ratio = 1.0 / 18.0;
        for denominator in (1..18).rev() {
            ratio = 1.0 / f64::from(denominator) + x * ratio;
        }
        ratio
    } else {
        -fs_math::det::ln(1.0 - x) / x
    }
}

// Saturated departure c and initial arrival i0 give
// incoming = i0*(1-exp(-lambda*t))/lambda
//          + c*t*[1-(1-exp(-lambda*t))/(lambda*t)].
// Both terms are nonnegative; this retains the O(t^2) return stream from an
// initially empty vapor without subtracting two O(t) gross quantities.
fn coated_streams(
    outgoing_rate: f64,
    incoming_rate: f64,
    net_rate: f64,
    decay_rate: f64,
    duration: f64,
) -> Result<(f64, f64, f64, f64), SublimationError> {
    let exposure = decay_rate * duration;
    if exposure.is_nan() || (duration > 0.0 && exposure == 0.0) {
        return Err(invalid("molecular turnover exposure is unrepresentable"));
    }
    let fraction = fs_math::decay::decayed_fraction(exposure);
    let effective_time = if exposure <= 0.125 {
        if exposure == 0.0 {
            duration
        } else {
            duration * (fraction / exposure)
        }
    } else {
        fraction / decay_rate
    };
    let outgoing = finite(outgoing_rate * duration)?;
    let returning = if exposure <= 0.125 {
        (outgoing * exposure) * second_exponential_ratio(exposure)
    } else {
        outgoing * (1.0 - fraction / exposure)
    };
    let incoming = finite(incoming_rate * effective_time + returning)?;
    Ok((
        outgoing,
        incoming,
        finite(net_rate * effective_time)?,
        exposure,
    ))
}

/// Integrate the original molecular law with a nonvolatile supporting surface.
///
/// `deposit.mass_kg` may be zero. Its specific enthalpy is admitted but never
/// changed or used as a substitute caloric law. Frozen area/temperature,
/// heterogeneous deposition and instantaneous reemission at a bare interface
/// are explicit caller-owned approximations. No original input is mutated.
///
/// A finite reservoir retains its original volume and temperature; the caller
/// must close gas sensible energy separately. No fractional gross-withdrawal
/// limiter applies: positivity follows from integrated backpressure and the
/// exact zero-film boundary. Positive asymptotic mass tails are preserved when
/// independently representable, even if subtracting the signed transfer from
/// the old mass rounds to zero. The separate residuals expose that encoding.
/// Zero duration, zero area and a sealed interface preserve both inventories.
#[allow(clippy::too_many_arguments)]
pub fn integrate_molecular_turnover<L: SublimationLaw + ?Sized>(
    law: &L,
    surface_temperature_k: f64,
    kinetics: SurfaceKinetics,
    deposit: SolidSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
) -> Result<MolecularTurnover, SublimationError> {
    if !nonnegative(deposit.mass_kg)
        || !deposit.specific_enthalpy_j_kg.is_finite()
        || !nonnegative(deposit.area_m2)
        || !nonnegative(duration_s)
    {
        return Err(invalid(
            "invalid supported molecular turnover state or duration",
        ));
    }
    let p = properties(law, surface_temperature_k)?;
    let vapor_pressure = reservoir.partial_pressure_pa(p.molar_mass_kg_mol)?;
    let gas_temperature = reservoir.temperature_k();
    let (a, b) = coefficients(p, surface_temperature_k, gas_temperature, kinetics)?;
    let mut result = MolecularTurnover {
        outgoing_mass_kg: 0.0,
        incoming_mass_kg: 0.0,
        transferred_mass_kg: 0.0,
        deposit_mass_kg: deposit.mass_kg,
        reservoir,
        extinction_time_s: None,
        stream_mass_residual_kg: 0.0,
        deposit_mass_residual_kg: 0.0,
        vapor_mass_residual_kg: 0.0,
    };
    if duration_s == 0.0 || deposit.area_m2 == 0.0 || (a == 0.0 && b == 0.0) {
        return Ok(result);
    }
    let outgoing_flux = finite(a * p.saturation_pressure_pa)?;
    let incoming_flux = finite(b * vapor_pressure)?;
    let net_flux = if surface_temperature_k == gas_temperature {
        a * (p.saturation_pressure_pa - vapor_pressure)
    } else {
        outgoing_flux - incoming_flux
    };
    let outgoing_rate = finite(deposit.area_m2 * outgoing_flux)?;
    let incoming_rate = finite(deposit.area_m2 * incoming_flux)?;
    let net_rate = finite(deposit.area_m2 * net_flux)?;
    if (outgoing_flux > 0.0 && outgoing_rate == 0.0)
        || (incoming_flux > 0.0 && incoming_rate == 0.0)
        || (net_flux != 0.0 && net_rate == 0.0)
    {
        return Err(invalid("positive molecular turnover rate underflows"));
    }
    if net_flux == 0.0 && outgoing_flux == incoming_flux {
        result.outgoing_mass_kg = finite(outgoing_rate * duration_s)?;
        result.incoming_mass_kg = result.outgoing_mass_kg;
        return Ok(result);
    }
    // Once bare and undersaturated, the gas and deposit inventories are fixed;
    // sticking molecules can still transfer heat by immediate thermalized return.
    if deposit.mass_kg == 0.0 && net_rate >= 0.0 {
        result.incoming_mass_kg = finite(incoming_rate * duration_s)?;
        result.outgoing_mass_kg = result.incoming_mass_kg;
        if net_rate > 0.0 {
            result.extinction_time_s = Some(0.0);
        }
        return Ok(result);
    }

    match reservoir {
        VaporReservoir::Prescribed { .. } => {
            let event = if net_rate > 0.0 {
                let time = deposit.mass_kg / net_rate;
                if time == 0.0 {
                    return Err(invalid("positive coating extinction time underflows"));
                }
                (time <= duration_s).then_some(time)
            } else {
                None
            };
            let coated_duration = event.unwrap_or(duration_s);
            let bare_duration = duration_s - coated_duration;
            result.outgoing_mass_kg =
                finite(outgoing_rate * coated_duration + incoming_rate * bare_duration)?;
            result.incoming_mass_kg = finite(incoming_rate * duration_s)?;
            result.transferred_mass_kg = match event {
                Some(_) => deposit.mass_kg,
                None => finite(net_rate * duration_s)?,
            };
            result.deposit_mass_kg = deposit.mass_kg - result.transferred_mass_kg;
            result.extinction_time_s = event;
        }
        VaporReservoir::Finite {
            vapor_mass_kg,
            volume_m3,
            temperature_k,
        } => {
            let pressure_per_mass =
                finite((MOLAR_GAS_CONSTANT / p.molar_mass_kg_mol) * temperature_k / volume_m3)?;
            let decay_rate = finite(deposit.area_m2 * b * pressure_per_mass)?;
            if !positive(decay_rate) {
                return Err(invalid(
                    "finite-vapor molecular turnover rate is unrepresentable",
                ));
            }
            let equilibrium_change = net_rate / decay_rate;
            let mut future_event = None;
            if net_rate > 0.0 && deposit.mass_kg < equilibrium_change {
                let linear_time = deposit.mass_kg / net_rate;
                if linear_time == 0.0 {
                    return Err(invalid("positive coating extinction time underflows"));
                }
                if linear_time <= duration_s {
                    let fraction = if equilibrium_change.is_finite() {
                        deposit.mass_kg / equilibrium_change
                    } else {
                        decay_rate * linear_time
                    };
                    if !(0.0..1.0).contains(&fraction) {
                        return Err(invalid("coating extinction fraction is unrepresentable"));
                    }
                    future_event = Some(linear_time * extinction_ratio(fraction));
                }
            }
            let event = future_event.filter(|time| *time <= duration_s);
            let coated_duration = event.unwrap_or(duration_s);
            let (outgoing, incoming, dm, exposure) = coated_streams(
                outgoing_rate,
                incoming_rate,
                net_rate,
                decay_rate,
                coated_duration,
            )?;
            result.outgoing_mass_kg = outgoing;
            result.incoming_mass_kg = incoming;
            result.extinction_time_s = event;
            let vapor_final;
            if event.is_some() {
                result.transferred_mass_kg = deposit.mass_kg;
                result.deposit_mass_kg = 0.0;
                vapor_final = finite(vapor_mass_kg + deposit.mass_kg)?;
                let final_pressure = VaporReservoir::Finite {
                    vapor_mass_kg: vapor_final,
                    volume_m3,
                    temperature_k,
                }
                .partial_pressure_pa(p.molar_mass_kg_mol)?;
                let bare_stream =
                    finite(deposit.area_m2 * b * final_pressure * (duration_s - coated_duration))?;
                result.outgoing_mass_kg = finite(outgoing + bare_stream)?;
                result.incoming_mass_kg = finite(incoming + bare_stream)?;
            } else {
                result.transferred_mass_kg = dm;
                result.deposit_mass_kg = deposit.mass_kg - dm;
                // Keep positive late tails. Formation-energy consumers must use
                // the same stable signed transfer and admit the encoding residual.
                if net_rate > 0.0 && exposure > core::f64::consts::LN_2 {
                    if equilibrium_change <= deposit.mass_kg {
                        result.deposit_mass_kg = finite(
                            (deposit.mass_kg - equilibrium_change)
                                + fs_math::decay::attenuate(equilibrium_change, exposure),
                        )?;
                    } else if result.deposit_mass_kg <= 0.0 {
                        let time = future_event.ok_or_else(|| {
                            invalid("coating remainder cannot be resolved before extinction")
                        })?;
                        result.deposit_mass_kg = finite(
                            (equilibrium_change - deposit.mass_kg)
                                * fs_math::det::expm1(decay_rate * (time - duration_s)),
                        )?;
                    }
                }
                vapor_final = if exposure <= core::f64::consts::LN_2 {
                    finite(vapor_mass_kg + dm)?
                } else {
                    finite(
                        fs_math::decay::attenuate(vapor_mass_kg, exposure)
                            + outgoing_rate
                                * (fs_math::decay::decayed_fraction(exposure) / decay_rate),
                    )?
                };
            }
            result.reservoir = VaporReservoir::Finite {
                vapor_mass_kg: vapor_final,
                volume_m3,
                temperature_k,
            };
            result.vapor_mass_residual_kg =
                balance(vapor_final, -vapor_mass_kg, -result.transferred_mass_kg)?;
        }
    }
    if !nonnegative(result.deposit_mass_kg)
        || !nonnegative(result.incoming_mass_kg)
        || !nonnegative(result.outgoing_mass_kg)
    {
        return Err(invalid("molecular turnover endpoint is unrepresentable"));
    }
    result.stream_mass_residual_kg = balance(
        result.outgoing_mass_kg,
        -result.incoming_mass_kg,
        -result.transferred_mass_kg,
    )?;
    result.deposit_mass_residual_kg = balance(
        result.deposit_mass_kg,
        -deposit.mass_kg,
        result.transferred_mass_kg,
    )?;
    Ok(result)
}

#[cfg(test)]
mod tests;
