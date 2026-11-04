//! Bordered elimination for a diagonal structural pencil plus physical ports.
//!
//! A = D + L R^T, R = B_hat and L = -i*mu*(B_hat Z_hat).
//! Small/zero bare diagonals remain in an explicitly bounded pivot block.
//! Neither the impedance nor a retained bare diagonal is inverted.

use super::impedance::{PortScaling, PreparedPortImpedance};
use super::{HarmonicResponseError, finite, invalid, lu_error, poll, scale, scale_complex};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_la::eigen::admit_jacobi_eigh;
use fs_la::eigen_complex::{
    LuComplex, lu_complex_blocked_with_cancel, lu_complex_blocked_workspace_bytes,
};
use fs_math::c64::C64;

const ELIMINATED: usize = usize::MAX;

#[derive(Clone, Debug)]
pub(super) struct CompactPortFactor {
    n: usize,
    ports: usize,
    retained: Vec<usize>,
    slots: Vec<usize>,
    b: Vec<f64>,
    bz: Vec<C64>,
    eliminated_response: Vec<C64>,
    frequency: f64,
    physical_frequency: f64,
    row_exponents: Vec<i32>,
    column_exponents: Vec<i32>,
    factor: LuComplex,
}

impl CompactPortFactor {
    #[allow(clippy::too_many_arguments)] // original physical scaling and residual policy
    pub(super) fn new(
        diagonal: &[C64],
        load: &PreparedPortImpedance,
        coordinates: &[i32],
        time_exponent: i32,
        frequency: f64,
        physical_frequency: f64,
        scaling: &PortScaling,
        maximum_relative_residual: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        let (n, p) = (diagonal.len(), load.ports);
        let (b, bz) = load.scaled_factors(coordinates, time_exponent, scaling, cx)?;
        let mut column_bounds = vec![0.0; p];
        for i in 0..n {
            poll(cx)?;
            for a in 0..p {
                column_bounds[a] += b[i * p + a].abs();
            }
        }
        // This selects an algebraic pivot block, not a physical approximation.
        // The final audit still uses every original loaded matrix entry.
        let relative_cutoff = (64.0 * f64::EPSILON / maximum_relative_residual).min(0.25);
        let mut retained = Vec::with_capacity(n);
        let mut slots = vec![ELIMINATED; n];
        let mut exact_zeros = 0;
        for i in 0..n {
            poll(cx)?;
            let mut load_bound = 0.0;
            for a in 0..p {
                load_bound += frequency * bz[i * p + a].abs() * column_bounds[a];
            }
            if !load_bound.is_finite() {
                return Err(invalid("compact port row bound is not representable"));
            }
            if diagonal[i] == C64::ZERO {
                exact_zeros += 1;
                // A rank-at-most-p update cannot remove more than p exact
                // null directions. Do not let rounded core pivots hide this.
                if exact_zeros > p {
                    return Err(HarmonicResponseError::Singular {
                        angular_frequency_rad_s: physical_frequency,
                        pivot: i,
                    });
                }
            }
            if diagonal[i] == C64::ZERO || diagonal[i].abs() <= relative_cutoff * load_bound {
                slots[i] = retained.len();
                retained.push(i);
            }
        }
        let k = p
            .checked_add(retained.len())
            .ok_or(invalid("compact port core extent overflow"))?;
        admit(n, p, k, load, cx)?;
        let mut eliminated_response = vec![C64::ZERO; n * p];
        let mut core = vec![C64::ZERO; k * k];
        let mut correction = vec![C64::ZERO; k * k];
        for a in 0..p {
            core[a * k + a] = C64::ONE;
        }
        for i in 0..n {
            poll(cx)?;
            if slots[i] != ELIMINATED {
                continue;
            }
            for a in 0..p {
                let load_entry = dynamic(bz[i * p + a], frequency);
                let value = if load_entry == C64::ZERO {
                    C64::ZERO
                } else {
                    load_entry / diagonal[i]
                };
                if !finite(value) {
                    return Err(invalid("compact modal elimination is not representable"));
                }
                eliminated_response[i * p + a] = value;
            }
            for a in 0..p {
                for c in 0..p {
                    let at = a * k + c;
                    add(
                        &mut core[at],
                        &mut correction[at],
                        eliminated_response[i * p + c].scale(b[i * p + a]),
                    );
                }
            }
        }
        for (slot, &i) in retained.iter().enumerate() {
            poll(cx)?;
            for a in 0..p {
                core[a * k + p + slot] = C64::from_re(-b[i * p + a]);
                core[(p + slot) * k + a] = dynamic(bz[i * p + a], frequency);
            }
            core[(p + slot) * k + p + slot] = diagonal[i];
        }
        if core.iter().any(|&v| !finite(v)) {
            return Err(invalid("compact port core is not representable"));
        }
        // Independent binary row/column equilibration changes no coefficient
        // or mode. It also keeps port-flow and modal-coordinate units separate.
        let (row_exponents, column_exponents) = equilibrate(&mut core, k, physical_frequency, cx)?;
        let factor = lu_complex_blocked_with_cancel(&core, k, 1, || cx.checkpoint().is_err())
            .map_err(|error| lu_error(error, physical_frequency))?;
        poll(cx)?;
        Ok(Self {
            n,
            ports: p,
            retained,
            slots,
            b,
            bz,
            eliminated_response,
            frequency,
            physical_frequency,
            row_exponents,
            column_exponents,
            factor,
        })
    }

    pub(super) fn dimension(&self) -> usize {
        self.ports + self.retained.len()
    }

    pub(super) fn solve(
        &self,
        diagonal: &[C64],
        rhs: &mut [C64],
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        let p = self.ports;
        let k = self.dimension();
        for i in 0..self.n {
            poll(cx)?;
            if self.slots[i] == ELIMINATED {
                rhs[i] = rhs[i] / diagonal[i];
            }
        }
        if rhs.iter().any(|&v| !finite(v)) {
            return Err(invalid("compact port source is not representable"));
        }
        // Contract each actual port source with its own exponent. An
        // uncoupled large modal response must not erase a small retained
        // force before row equilibration can restore its natural scale.
        let mut source = vec![C64::ZERO; k];
        let mut source_exponents = vec![0_i32; k];
        for a in 0..p {
            poll(cx)?;
            let mut exponent = i32::MIN;
            for i in 0..self.n {
                let coefficient = self.b[i * p + a];
                if self.slots[i] != ELIMINATED || coefficient == 0.0 {
                    continue;
                }
                let coefficient_exponent = quadratic_fraction(coefficient).1;
                for component in [rhs[i].re, rhs[i].im] {
                    if component != 0.0 {
                        exponent =
                            exponent.max(quadratic_fraction(component).1 + coefficient_exponent);
                    }
                }
            }
            if exponent == i32::MIN {
                continue;
            }
            source_exponents[a] = exponent;
            let mut correction = C64::ZERO;
            for i in 0..self.n {
                if self.slots[i] == ELIMINATED && self.b[i * p + a] != 0.0 {
                    add(
                        &mut source[a],
                        &mut correction,
                        scaled_real_product(rhs[i], self.b[i * p + a], -exponent)?,
                    );
                }
            }
        }
        for (slot, &i) in self.retained.iter().enumerate() {
            source[p + slot] = rhs[i];
        }
        if source.iter().any(|&v| !finite(v)) {
            return Err(invalid("compact port contraction is not representable"));
        }
        let mut amplitude = i32::MIN;
        for (i, value) in source.iter().enumerate() {
            for component in [value.re, value.im] {
                if component != 0.0 {
                    amplitude = amplitude.max(
                        quadratic_fraction(component).1 + source_exponents[i]
                            - self.row_exponents[i],
                    );
                }
            }
        }
        if amplitude == i32::MIN {
            amplitude = 0;
        }
        for (i, value) in source.iter_mut().enumerate() {
            *value = scale_complex(
                *value,
                source_exponents[i] - self.row_exponents[i] - amplitude,
            )?;
        }
        self.factor
            .solve_with_cancel(&mut source, || cx.checkpoint().is_err())
            .map_err(|error| lu_error(error, self.physical_frequency))?;
        for (i, value) in source.iter_mut().enumerate() {
            *value = scale_complex(*value, amplitude - self.column_exponents[i])?;
        }
        for i in 0..self.n {
            poll(cx)?;
            if self.slots[i] != ELIMINATED {
                rhs[i] = source[p + self.slots[i]];
                continue;
            }
            let mut reaction = C64::ZERO;
            let mut error = C64::ZERO;
            for a in 0..p {
                add(
                    &mut reaction,
                    &mut error,
                    self.eliminated_response[i * p + a] * source[a],
                );
            }
            rhs[i] = rhs[i] - reaction;
            if !finite(rhs[i]) {
                return Err(invalid("compact modal reconstruction is not representable"));
            }
        }
        Ok(())
    }

    // Stream the original loaded coefficient for the unchanged componentwise
    // audit. No bound on |L||R^T| replaces the actual |A_ij| denominator.
    pub(super) fn entry(
        &self,
        diagonal: &[C64],
        row: usize,
        column: usize,
    ) -> Result<C64, HarmonicResponseError> {
        let mut value = C64::ZERO;
        for a in 0..self.ports {
            value = value + self.bz[row * self.ports + a].scale(self.b[column * self.ports + a]);
        }
        let mut coefficient = dynamic(value, self.frequency);
        if row == column {
            coefficient = diagonal[row] + coefficient;
        }
        if !finite(coefficient) {
            return Err(invalid("loaded dynamic stiffness is not representable"));
        }
        Ok(coefficient)
    }
}

fn dynamic(value: C64, frequency: f64) -> C64 {
    C64::new(value.im * frequency, -value.re * frequency)
}

fn add(sum: &mut C64, correction: &mut C64, term: C64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

// Combine product and row-normalization exponents before materializing a
// term. Multiplying tiny/huge dimensional factors first could lose a finite
// equilibrated source. Nonzero underflow continues to refuse explicitly.
fn scaled_real_product(
    value: C64,
    coefficient: f64,
    exponent: i32,
) -> Result<C64, HarmonicResponseError> {
    if coefficient == 0.0 || value == C64::ZERO {
        return Ok(C64::ZERO);
    }
    let coefficient_exponent = quadratic_fraction(coefficient).1;
    let normalized_coefficient = scale(coefficient, -coefficient_exponent)?;
    let component = |v: f64| -> Result<f64, HarmonicResponseError> {
        if v == 0.0 {
            return Ok(v);
        }
        let value_exponent = quadratic_fraction(v).1;
        scale(
            scale(v, -value_exponent)? * normalized_coefficient,
            value_exponent + coefficient_exponent + exponent,
        )
    };
    Ok(C64::new(component(value.re)?, component(value.im)?))
}

fn vector_exponent(values: &[C64]) -> i32 {
    values
        .iter()
        .flat_map(|v| [v.re, v.im])
        .filter(|v| *v != 0.0)
        .map(|v| quadratic_fraction(v).1)
        .max()
        .unwrap_or(0)
}

fn equilibrate(
    core: &mut [C64],
    k: usize,
    physical_frequency: f64,
    cx: &Cx<'_>,
) -> Result<(Vec<i32>, Vec<i32>), HarmonicResponseError> {
    let mut rows = Vec::with_capacity(k);
    for i in 0..k {
        poll(cx)?;
        let row = &core[i * k..(i + 1) * k];
        if row.iter().all(|&v| v == C64::ZERO) {
            return Err(HarmonicResponseError::Singular {
                angular_frequency_rad_s: physical_frequency,
                pivot: i,
            });
        }
        let exponent = vector_exponent(row);
        rows.push(exponent);
        for value in &mut core[i * k..(i + 1) * k] {
            *value = scale_complex(*value, -exponent)?;
        }
    }
    let mut columns = Vec::with_capacity(k);
    for j in 0..k {
        poll(cx)?;
        let mut exponent = i32::MIN;
        for i in 0..k {
            for value in [core[i * k + j].re, core[i * k + j].im] {
                if value != 0.0 {
                    exponent = exponent.max(quadratic_fraction(value).1);
                }
            }
        }
        if exponent == i32::MIN {
            return Err(HarmonicResponseError::Singular {
                angular_frequency_rad_s: physical_frequency,
                pivot: j,
            });
        }
        columns.push(exponent);
        for i in 0..k {
            core[i * k + j] = scale_complex(core[i * k + j], -exponent)?;
        }
    }
    Ok((rows, columns))
}

fn admit(
    n: usize,
    p: usize,
    k: usize,
    load: &PreparedPortImpedance,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    poll(cx)?;
    let np = n
        .checked_mul(p)
        .ok_or(invalid("compact port extent overflow"))?;
    let pp = p
        .checked_mul(p)
        .ok_or(invalid("compact port extent overflow"))?;
    let kk = k
        .checked_mul(k)
        .ok_or(invalid("compact port core extent overflow"))?;
    let workspace = lu_complex_blocked_workspace_bytes(k, 1)
        .ok_or(invalid("compact port factor workspace overflow"))?;
    let bytes = 64_usize
        .checked_mul(np)
        .and_then(|x| pp.checked_mul(256).and_then(|y| x.checked_add(y)))
        .and_then(|x| p.checked_mul(128).and_then(|y| x.checked_add(y)))
        .and_then(|x| kk.checked_mul(144).and_then(|y| x.checked_add(y)))
        .and_then(|x| k.checked_mul(256).and_then(|y| x.checked_add(y)))
        .and_then(|x| x.checked_add(workspace))
        .ok_or(invalid("compact port payload overflow"))?;
    if bytes > load.controls.maximum_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "compact load numerical payload bytes",
            requested: bytes,
            cap: load.controls.maximum_payload_bytes,
        });
    }
    let work = np
        .checked_mul(n)
        .and_then(|x| {
            np.checked_mul(p)
                .and_then(|y| y.checked_mul(2))
                .and_then(|y| x.checked_add(y))
        })
        .and_then(|x| np.checked_mul(8).and_then(|y| x.checked_add(y)))
        .and_then(|x| kk.checked_mul(k).and_then(|y| x.checked_add(y)))
        .and_then(|x| u64::try_from(x).ok())
        .ok_or(invalid("compact port work overflow"))?;
    if work > load.controls.maximum_projection_terms {
        return Err(HarmonicResponseError::WorkBudget {
            requested: work,
            cap: load.controls.maximum_projection_terms,
        });
    }
    admit_jacobi_eigh(k)
        .map_err(|_| invalid("compact port core exceeds fs-la bounded dense work cap"))?;
    Ok(())
}
