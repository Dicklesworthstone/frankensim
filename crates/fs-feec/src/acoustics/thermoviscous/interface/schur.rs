//! Complete work-dual harmonic connection, independent of geometric binding.
use super::{
    C64, Dd, LuComplex, ThermoviscousError, finite, invalid, lu_complex_with_cancel, poll,
};

pub(super) struct ConformingSchur {
    factor: LuComplex,
    rows: Vec<f64>,
    columns: Vec<f64>,
    tolerance: f64,
}

impl ConformingSchur {
    pub(super) fn prepare<C: FnMut() -> bool>(
        loads: [&[f64]; 2],
        sizes: [usize; 2],
        count: usize,
        tolerance: f64,
        solve: &mut impl FnMut(usize, &[C64], &mut C) -> Result<Vec<C64>, ThermoviscousError>,
        cancelled: &mut C,
    ) -> Result<Self, ThermoviscousError> {
        let mut schur = vec![C64::ZERO; count * count];
        // One original response column is live. The complete volume response
        // is never fitted or truncated, and no n-by-count cache is retained.
        for side in 0..2 {
            let mut rhs = vec![C64::ZERO; sizes[side]];
            for column in 0..count {
                poll(cancelled)?;
                for (row, value) in rhs.iter_mut().enumerate() {
                    *value = C64::new(loads[side][row * count + column], 0.0);
                }
                let response = solve(side, &rhs, cancelled)?;
                for row in 0..count {
                    schur[row * count + column] =
                        schur[row * count + column] + dual(loads[side], count, row, &response);
                }
            }
        }
        let (factor, rows, columns) = factor_schur(&mut schur, count, cancelled)?;
        Ok(Self {
            factor,
            rows,
            columns,
            tolerance,
        })
    }

    /// Correct against the physical trace of freshly corrected ORIGINAL
    /// volume solves, not just the rounded Schur matrix. At most three small
    /// corrections; a failed constraint never returns a best iterate.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn response<C: FnMut() -> bool>(
        &self,
        loads: [&[f64]; 2],
        mechanical: usize,
        velocity_end: [usize; 2],
        temperature_start: [usize; 2],
        rhs: [&[C64]; 2],
        solve: &mut impl FnMut(usize, &[C64], &mut C) -> Result<Vec<C64>, ThermoviscousError>,
        cancelled: &mut C,
    ) -> Result<([Vec<C64>; 2], Vec<C64>), ThermoviscousError> {
        poll(cancelled)?;
        let count = self.rows.len();
        let mut fields = [solve(0, rhs[0], cancelled)?, solve(1, rhs[1], cancelled)?];
        let mut reference_scale = [0.0_f64; 2];
        let mut reaction = vec![C64::ZERO; count];
        for iteration in 0..4 {
            poll(cancelled)?;
            // Preserve the free-field scale when a constraint arrests a trace.
            // Velocity and theta use separate units; pressure cannot mask shear.
            for side in 0..2 {
                for value in &fields[side][..velocity_end[side]] {
                    reference_scale[0] = reference_scale[0].max(value.abs());
                }
                for value in &fields[side][temperature_start[side]..] {
                    reference_scale[1] = reference_scale[1].max(value.abs());
                }
            }
            let mut residual = vec![C64::ZERO; count];
            let mut defect: f64 = 0.0;
            for a in 0..count {
                residual[a] = C64::ZERO
                    - (dual(loads[0], count, a, &fields[0]) + dual(loads[1], count, a, &fields[1]));
                let physical_scale = reference_scale[usize::from(a >= mechanical)];
                let scale: f64 = loads
                    .iter()
                    .flat_map(|side| side.chunks_exact(count))
                    .map(|row| row[a].abs() * physical_scale)
                    .sum();
                let error = residual[a].abs();
                if !scale.is_finite() || !error.is_finite() {
                    return Err(ThermoviscousError::NonFinite);
                }
                if error > 0.0 {
                    defect = defect.max(if scale > 0.0 {
                        error / scale
                    } else {
                        f64::INFINITY
                    });
                }
            }
            if defect <= self.tolerance {
                return Ok((fields, reaction));
            }
            if iteration == 3 {
                break;
            }
            for (entry, &scale) in residual.iter_mut().zip(&self.rows) {
                *entry = entry.scale(scale);
            }
            self.factor
                .solve_with_cancel(&mut residual, &mut *cancelled)?;
            for a in 0..count {
                reaction[a] = reaction[a] + residual[a].scale(self.columns[a]);
                if !finite(reaction[a]) {
                    return Err(ThermoviscousError::NonFinite);
                }
            }
            for side in 0..2 {
                let mut coupled = rhs[side].to_vec();
                for (row, value) in coupled.iter_mut().enumerate() {
                    poll(cancelled)?;
                    for a in 0..count {
                        *value = *value + reaction[a].scale(loads[side][row * count + a]);
                    }
                }
                fields[side] = solve(side, &coupled, cancelled)?;
            }
        }
        Err(invalid("conforming harmonic continuity residual exceeded"))
    }
}

pub(super) fn dual(loads: &[f64], count: usize, column: usize, field: &[C64]) -> C64 {
    let mut re = Dd::ZERO;
    let mut im = Dd::ZERO;
    for (row, value) in field.iter().enumerate() {
        let coefficient = Dd::from_f64(loads[row * count + column]);
        re = re + coefficient * Dd::from_f64(value.re);
        im = im + coefficient * Dd::from_f64(value.im);
    }
    C64::new(re.to_f64(), im.to_f64())
}

fn factor_schur(
    a: &mut [C64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(LuComplex, Vec<f64>, Vec<f64>), ThermoviscousError> {
    let mut rows = vec![0.0; n];
    let mut columns = vec![0.0; n];
    for i in 0..n {
        poll(cancelled)?;
        let maximum = a[i * n..(i + 1) * n]
            .iter()
            .map(|z| z.abs())
            .fold(0.0, f64::max);
        if !maximum.is_finite() || maximum == 0.0 {
            return Err(invalid("nonzero finite interface rank required"));
        }
        rows[i] = 1.0 / maximum;
        for j in 0..n {
            a[i * n + j] = a[i * n + j].scale(rows[i]);
        }
    }
    for j in 0..n {
        poll(cancelled)?;
        let maximum = (0..n).map(|i| a[i * n + j].abs()).fold(0.0, f64::max);
        if !maximum.is_finite() || maximum == 0.0 {
            return Err(invalid("independent finite interface loads required"));
        }
        columns[j] = 1.0 / maximum;
        for i in 0..n {
            a[i * n + j] = a[i * n + j].scale(columns[j]);
        }
    }
    if rows.iter().chain(&columns).any(|v| !v.is_finite()) || a.iter().any(|&z| !finite(z)) {
        return Err(ThermoviscousError::NonFinite);
    }
    Ok((
        lu_complex_with_cancel(a, n, &mut *cancelled)?,
        rows,
        columns,
    ))
}

#[cfg(test)]
mod tests;
