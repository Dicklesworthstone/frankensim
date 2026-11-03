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

fn weight(frame: Frame, masses: &[f64], i: usize) -> f64 {
    match frame {
        Frame::MassAveraged => masses[i],
        Frame::MolarAveraged => 1.0,
    }
}
