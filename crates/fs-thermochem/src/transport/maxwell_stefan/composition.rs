//! Exact reduction by the physical composition, without a trace cutoff.
//!
//! A zero-composition row is diagonal: its flux is still driven by its own
//! force, and contributes to the active rows and the reference-frame constraint.
//! Only the positive-composition block needs a coupled solve. The caller checks
//! every original row again after primary flux rounding.

use super::{
    Frame, MAX_MAXWELL_STEFAN_SPECIES_V1, MaxwellStefanErrorV1, MaxwellStefanWorkspaceV1, poll,
    positive, scaled_product_ratio, solve, sum,
};

pub(super) fn evaluate(
    fractions: &[f64],
    scaled_masses: &[f64],
    gauge: usize,
    frame: Frame,
    drive_scale: f64,
    work: &mut MaxwellStefanWorkspaceV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaxwellStefanErrorV1> {
    let n = fractions.len();
    if n == 2 {
        // There is one pair, so D_min/D_01 = 1. Use the original row not
        // replaced by the frame constraint; the caller also checks the other.
        poll(cancelled)?;
        let i = 1 - gauge;
        let wi = weight(frame, scaled_masses, i);
        let wg = weight(frame, scaled_masses, gauge);
        let denominator = fractions[i] * wi + fractions[gauge] * wg;
        if !positive(denominator) {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
        let rhs = -work.driving[i] / drive_scale;
        work.right[i] = scaled_product_ratio(&[rhs, wg], &[denominator])?;
        work.right[gauge] = scaled_product_ratio(&[-rhs, wi], &[denominator])?;
        return Ok(());
    }

    // Equal active pair diffusivities reduce the ORIGINAL friction matrix to
    // S I - x 1^T, S=sum(x). Do not approximate unequal coefficients, or let
    // an absent/absent pair choose the model: its friction is identically zero.
    // This avoids the cubic positive-composition solve. Complete admission,
    // original-row checks and cancellation-safe pair rescue remain O(n^2).
    let mut common_diffusivity = true;
    for i in 0..n {
        poll(cancelled)?;
        for j in i + 1..n {
            if (fractions[i] > 0.0 || fractions[j] > 0.0) && work.ratios[i * n + j] != 1.0 {
                common_diffusivity = false;
                break;
            }
        }
        if !common_diffusivity {
            break;
        }
    }
    if common_diffusivity {
        evaluate_common(
            fractions,
            scaled_masses,
            gauge,
            frame,
            drive_scale,
            work,
            cancelled,
        )?;
        return Ok(());
    }

    let mut active = [0_usize; MAX_MAXWELL_STEFAN_SPECIES_V1];
    let mut count = 0;
    for (i, &x) in fractions.iter().enumerate() {
        if x > 0.0 {
            active[count] = i;
            count += 1;
        }
    }
    debug_assert!(count > 0);
    work.right.fill(0.0);
    for i in 0..n {
        poll(cancelled)?;
        if fractions[i] != 0.0 {
            continue;
        }
        let diagonal = sum(active[..count]
            .iter()
            .map(|&j| fractions[j] * work.ratios[i * n + j]));
        if !positive(diagonal) {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
        work.right[i] = (-work.driving[i] / drive_scale) / diagonal;
        if !work.right[i].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "zero-composition diffusion response",
            ));
        }
    }
    let frame_rhs = -sum((0..n)
        .filter(|&i| fractions[i] == 0.0)
        .map(|i| weight(frame, scaled_masses, i) * work.right[i]));
    if !frame_rhs.is_finite() {
        return Err(MaxwellStefanErrorV1::Arithmetic(
            "zero-composition frame contribution",
        ));
    }
    if count == 1 {
        // No coupled matrix is needed: the single carrier closes the frame.
        poll(cancelled)?;
        work.right[gauge] = frame_rhs / weight(frame, scaled_masses, gauge);
        if !work.right[gauge].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("single-carrier diffusion"));
        }
        return Ok(());
    }

    let mut rhs = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
    let mut terms = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
    for row in 0..count {
        poll(cancelled)?;
        let i = active[row];
        if i == gauge {
            for col in 0..count {
                work.matrix[row * count + col] = weight(frame, scaled_masses, active[col]);
            }
            rhs[row] = frame_rhs;
            continue;
        }
        for col in 0..count {
            let j = active[col];
            work.matrix[row * count + col] = if i == j {
                sum(active[..count]
                    .iter()
                    .copied()
                    .filter(|&k| k != i)
                    .map(|k| fractions[k] * work.ratios[i * n + k]))
            } else {
                -fractions[i] * work.ratios[i * n + j]
            };
        }
        terms[..n].fill(0.0);
        for j in 0..n {
            if fractions[j] == 0.0 && work.right[j] != 0.0 {
                // Keep a small fraction/rate from erasing a representable
                // contribution carried by a large dimensionless response.
                terms[j] = scaled_product_ratio(
                    &[fractions[i], work.ratios[i * n + j], work.right[j]],
                    &[],
                )?;
            }
        }
        rhs[row] = -work.driving[i] / drive_scale + sum(terms[..n].iter().copied());
        if !rhs[row].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "reduced diffusion driving force",
            ));
        }
    }
    solve(
        &mut work.matrix[..count * count],
        &mut rhs[..count],
        &mut work.row_scales[..count],
        count,
        cancelled,
    )?;
    for row in 0..count {
        work.right[active[row]] = rhs[row];
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn evaluate_common(
    fractions: &[f64],
    scaled_masses: &[f64],
    gauge: usize,
    frame: Frame,
    drive_scale: f64,
    work: &mut MaxwellStefanWorkspaceV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaxwellStefanErrorV1> {
    let n = fractions.len();
    let mut forcing = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
    for i in 0..n {
        forcing[i] = -work.driving[i] / drive_scale;
    }
    // Retain the same omitted original row as the general solve. A caller may
    // admit a small Gibbs-Duhem defect; do not distribute it over other rows.
    // The gauge equation is still checked against its original force later.
    forcing[gauge] = -sum((0..n).filter(|&i| i != gauge).map(|i| forcing[i]));
    let fraction_sum = sum(fractions.iter().copied());
    let frame_weight = sum((0..n).map(|i| weight(frame, scaled_masses, i) * fractions[i]));
    if !positive(fraction_sum) || !positive(frame_weight) {
        return Err(MaxwellStefanErrorV1::SingularSystem);
    }
    let frame_force = sum((0..n).map(|i| weight(frame, scaled_masses, i) * forcing[i]));
    for i in 0..n {
        poll(cancelled)?;
        if i == gauge {
            continue;
        }
        // z=(f-x*(w.f)/(w.x))/S closes w.z=0. Keep multiplication and
        // division together so a small positive fraction is not lost just
        // because (w.f)/(w.x) overflows before it is multiplied by x_i.
        let correction = scaled_product_ratio(&[fractions[i], frame_force], &[frame_weight])?;
        let difference = forcing[i] - correction;
        let scale = forcing[i].abs() + correction.abs();
        work.right[i] = if scale > 0.0 && difference.abs() <= 1.490_116_119_384_765_6e-8 * scale {
            // Near cancellation, expand the same numerator into off-diagonal
            // pairs: sum_j w_j*(x_j*f_i-x_i*f_j)/(S*w.x). The self-pair
            // cancels algebraically, retaining a tiny heavy-species flux that
            // subtracting two O(1) terms (or dense elimination) can erase.
            let mut pairs = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
            for j in 0..n {
                if j != i {
                    let wj = weight(frame, scaled_masses, j);
                    pairs[j] = scaled_product_ratio(
                        &[wj, fractions[j], forcing[i]],
                        &[fraction_sum, frame_weight],
                    )? - scaled_product_ratio(
                        &[wj, fractions[i], forcing[j]],
                        &[fraction_sum, frame_weight],
                    )?;
                }
            }
            sum(pairs[..n].iter().copied())
        } else {
            difference / fraction_sum
        };
        if !work.right[i].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "common-diffusivity response",
            ));
        }
    }
    // Close the original frame directly; subtracting two order-one terms
    // for the dominant composition can erase its tiny compensating flux.
    work.right[gauge] = -scaled_product_ratio(
        &[sum((0..n)
            .filter(|&i| i != gauge)
            .map(|i| weight(frame, scaled_masses, i) * work.right[i]))],
        &[weight(frame, scaled_masses, gauge)],
    )?;
    Ok(())
}

fn weight(frame: Frame, masses: &[f64], i: usize) -> f64 {
    match frame {
        Frame::MassAveraged => masses[i],
        Frame::MolarAveraged => 1.0,
    }
}
