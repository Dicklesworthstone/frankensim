//! One conservative face correction bounded by original caloric entropy.
use super::*;

pub(super) fn limit_face(
    thermal: &DepositionSimulation<'_>,
    up: CellId,
    down: CellId,
    mut increment: FaceIncrement,
    differentiate: bool,
    calorics: Option<&ReconstructionCalorics<'_>>,
) -> Result<FaceIncrement, DepositionError> {
    let donor = &thermal.cells()[&up];
    let receiver = &thermal.cells()[&down];
    if donor.material_id != receiver.material_id {
        return Err(invalid(
            "thermal entropy face limiting requires one original caloric law",
        ));
    }
    let law = thermal.materials[&donor.material_id];
    let state = |id, h| -> Result<(f64, f64), DepositionError> {
        let (temperature, slope) = if let Some(calorics) = calorics {
            calorics(id)?
        } else {
            (law.temperature_k(h)?, 0.0)
        };
        if !positive(temperature) || !slope.is_finite() || slope < 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((temperature, slope))
    };
    let (ta, slope_a) = state(up, donor.specific_h_j_kg)?;
    let (tb, slope_b) = state(down, receiver.specific_h_j_kg)?;
    let dh = finite(donor.specific_h_j_kg - receiver.specific_h_j_kg)?;
    let ds = finite(
        law.specific_entropy_change_j_kg_k(receiver.specific_h_j_kg, donor.specific_h_j_kg)?
            .ok_or(invalid(
                "thermal entropy face limiting requires original entropy differences",
            ))?,
    )?;
    let donor_gap = finite(dh / tb - ds)?;
    let receiver_gap = finite(ds - dh / ta)?;
    let reciprocal_jump = finite((ta - tb) / ta / tb)?;
    let correction = finite(increment.value * reciprocal_jump)?;
    // Recovering the concavity gap subtracts O(dh) quantities to obtain an
    // O(dh^2) term. Its uncertainty must use these finite increments, not an
    // arbitrary formation enthalpy or absolute entropy reference.
    let roundoff = finite(
        64.0 * f64::EPSILON * (dh / tb).abs()
            + 64.0 * f64::EPSILON * (dh / ta).abs()
            + 128.0 * f64::EPSILON * ds.abs()
            + 64.0 * f64::EPSILON * correction.abs(),
    )?;
    if donor_gap < -roundoff || receiver_gap < -roundoff {
        return Err(invalid(
            "original caloric entropy is not concave across the reconstructed face",
        ));
    }
    // At coexistence both entropy variables coincide. All bulk latent traces
    // share the same thermal effort and no quotient or apparent Cp is needed.
    if ta == tb {
        return Ok(increment);
    }
    increment.guard_power_roundoff = true;
    if finite(donor_gap + correction)? >= -roundoff {
        return Ok(increment);
    }
    if donor_gap <= roundoff
        || receiver_gap <= roundoff
        || (ta - tb).abs() <= 64.0 * f64::EPSILON * ta.max(tb)
    {
        // A cancellation-dominated active bound cannot supply a trustworthy
        // quotient tangent. Retain the original donor flux, not guessed data.
        return Ok(FaceIncrement {
            value: 0.0,
            partials: None,
            guard_power_roundoff: false,
        });
    }
    let factor = finite(donor_gap / -correction)?;
    if !(0.0..1.0).contains(&factor) {
        return Err(DepositionError::Arithmetic);
    }
    let value = finite(increment.value * factor)?;
    if !positive(law.temperature_k(finite(donor.specific_h_j_kg + value)?)?) {
        return Err(DepositionError::MaterialDomain);
    }
    let partials = if differentiate {
        if calorics.is_none() {
            return Err(invalid(
                "entropy-limited tangent requires original endpoint dT/dh",
            ));
        }
        // P=D_A+c*(1/T_B-1/T_A)=0. Differentiation uses ds/dh=1/T
        // and the original endpoint dT/dh, including one-sided phase tangents.
        // These are derivatives of the FULL face h; subtract the donor's unit
        // derivative because the native upwind tangent remains authoritative.
        let a = face_derivative(donor_gap, slope_a, tb / (ta - tb))?;
        let b = face_derivative(receiver_gap, slope_b, ta / (ta - tb))?;
        Some([(up, finite(a - 1.0)?), (down, b)])
    } else {
        None
    };
    Ok(FaceIncrement {
        value,
        partials,
        guard_power_roundoff: true,
    })
}

fn face_derivative(gap: f64, slope: f64, ratio: f64) -> Result<f64, DepositionError> {
    let value = if ratio.abs() <= 1.0 {
        finite(finite(gap * ratio)? * finite(slope * ratio)?)?
    } else {
        finite(finite(finite(gap * slope)? * ratio)? * ratio)?
    };
    if value < 0.0 || (gap != 0.0 && slope != 0.0 && value == 0.0) {
        return Err(DepositionError::Arithmetic);
    }
    Ok(value)
}

/// Do not let rounding against a distant physical component reference apply a
/// larger antidiffusive correction than the admitted face increment. Moving one
/// representable power toward the donor normally suffices; otherwise retain
/// donor flux. This changes the shared power, never its later entropy account.
pub(super) fn apply_correction(old_power: f64, correction: f64) -> Result<f64, DepositionError> {
    let mut power = finite(old_power + correction)?;
    let exceeds = |delta: f64| {
        (correction > 0.0 && delta > correction) || (correction < 0.0 && delta < correction)
    };
    if exceeds(finite(power - old_power)?) {
        power = if correction > 0.0 {
            power.next_down()
        } else {
            power.next_up()
        };
    }
    let applied = finite(power - old_power)?;
    if exceeds(applied)
        || (correction > 0.0 && applied < 0.0)
        || (correction < 0.0 && applied > 0.0)
    {
        return Ok(old_power);
    }
    Ok(power)
}
