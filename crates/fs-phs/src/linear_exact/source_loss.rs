//! Positive source-square observation of an existing held-input trajectory.
use super::{
    DenseExponentialBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite, positive_heat, refuse,
};
use crate::{PhsError, PreparedStepError};
use fs_math::{
    det,
    eft::{two_prod, two_sum},
};

/// Cached integral of `||S*x(t)+D*u||²` for the unchanged `x'=A*x+B*u`.
///
/// This is a read-only channel observation. It owns no physical state, advances
/// no clock, and makes no passivity claim for A. Its caller retains the physical
/// meaning of the ORIGINAL source rows S and the full trajectory's admission.
/// No symmetric decomposition, replacement resistance, endpoint quadrature or
/// total-energy deficit supplies this channel. Positive Gauss/QR doubling is
/// the same kernel as the existing exact port-Hamiltonian heat owner.
#[derive(Clone, Debug)]
pub struct PreparedSourceLoss {
    n: usize,
    m: usize,
    source_count: usize,
    flow: Vec<f64>,
    input: Vec<f64>,
    source: Vec<f64>,
    feedthrough: Vec<f64>,
    root: Vec<f64>,
    rate: f64,
    duration: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl PreparedSourceLoss {
    /// Prepare a source channel once for its exact generator, input map and dt.
    /// A is n×n, B is n×m and source rows have n entries each (at most n rows).
    /// The 2n observation coordinates must fit the existing 256-coordinate cap.
    /// Explicit byte/work allowances cover caches and conservative preparation
    /// before allocation. One bounded positive integral is the largest atomic
    /// operation; cancellation is checked on either side and on every warm row.
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        n: usize,
        m: usize,
        flow: &[f64],
        input: &[f64],
        source_rows: &[f64],
        duration: f64,
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            n,
            m,
            flow,
            input,
            source_rows,
            None,
            duration,
            budget,
            cancelled,
        )
    }

    /// Integrate physical source squares with direct held-input feedthrough.
    ///
    /// Each row is `S*x+D*u`, in units of sqrt(W). For `r` source rows,
    /// `source_rows` is row-major r×n and `feedthrough_rows` is r×m, with
    /// `r <= n`. This includes a real source resistor's
    /// `sqrt(G)*(V_bus-V_source)` and driven acoustic boundary losses.
    /// Squaring the state contribution alone would omit physical source work.
    ///
    /// No extra state or clock is introduced. With held u, the existing
    /// observation coordinates obey `y'=S*v`, `v'=A*v`; only their initial
    /// values change to `y(0)=S*x(0)+D*u`, `v(0)=A*x(0)+B*u`.
    /// Preparation therefore retains the original 2n dimension and positive
    /// Gauss/QR factor. Inputs must remain constant over this entire interval.
    /// Warm evaluation uses compensated products/sums to retain small physical
    /// differences between large source and state values, without allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_feedthrough_controlled(
        n: usize,
        m: usize,
        flow: &[f64],
        input: &[f64],
        source_rows: &[f64],
        feedthrough_rows: &[f64],
        duration: f64,
        budget: DenseExponentialBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            n,
            m,
            flow,
            input,
            source_rows,
            Some(feedthrough_rows),
            duration,
            budget,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        n: usize,
        m: usize,
        flow: &[f64],
        input: &[f64],
        source_rows: &[f64],
        feedthrough_rows: Option<&[f64]>,
        duration: f64,
        budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let square = n
            .checked_mul(n)
            .ok_or_else(|| refuse("source loss extent overflow"))?;
        let inputs = n
            .checked_mul(m)
            .ok_or_else(|| refuse("source loss input extent overflow"))?;
        let d = n
            .checked_mul(2)
            .ok_or_else(|| refuse("source loss doubled extent overflow"))?;
        let bytes = d
            .checked_mul(d)
            .and_then(|v| v.checked_mul(32))
            .and_then(|v| v.checked_add(square.checked_mul(3)?))
            .and_then(|v| v.checked_add(inputs))
            .and_then(|v| {
                v.checked_add(if feedthrough_rows.is_some() {
                    inputs
                } else {
                    0
                })
            })
            .and_then(|v| v.checked_add(source_rows.len()))
            .and_then(|v| v.checked_mul(8));
        let work = (d as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(2048))
            .and_then(|v| v.checked_add(inputs as u64));
        if n == 0
            || d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || d > budget.max_dimension
            || bytes.is_none_or(|v| v > budget.max_dense_bytes)
            || work.is_none_or(|v| v > budget.max_multiply_terms)
        {
            return Err(
                refuse("source loss exceeds complete dimension/byte/work allowance").into(),
            );
        }
        if flow.len() != square
            || input.len() != inputs
            || source_rows.len() % n != 0
            || source_rows.len() > square
            || !finite(flow)
            || !finite(input)
            || !finite(source_rows)
            || !duration.is_finite()
            || duration <= 0.0
        {
            return Err(refuse(
                "source loss needs finite complete generator, ports, source rows and clock",
            )
            .into());
        }
        if feedthrough_rows
            .is_some_and(|rows| rows.len() != (source_rows.len() / n) * m || !finite(rows))
        {
            return Err(
                refuse("source loss needs one finite feedthrough row per source row").into(),
            );
        }
        let mut source = source_rows.to_vec();
        source.resize(square, 0.0);
        let feedthrough = match feedthrough_rows {
            Some(rows) if rows.iter().any(|&value| value != 0.0) => {
                let mut values = rows.to_vec();
                values.resize(inputs, 0.0);
                values
            }
            _ => Vec::new(),
        };
        let (rate, root) = if source.iter().all(|&value| value == 0.0) && feedthrough.is_empty() {
            (1.0, Vec::new())
        } else {
            positive_heat::source_integral(flow, &source, n, duration)?
        };
        poll(&mut cancelled)?;
        Ok(Self {
            n,
            m,
            source_count: source_rows.len() / n,
            flow: flow.to_vec(),
            input: input.to_vec(),
            source,
            feedthrough,
            root,
            rate,
            duration,
        })
    }

    /// Exact held-input interval represented by this immutable cache.
    pub const fn duration(&self) -> f64 {
        self.duration
    }

    /// Integrate the original positive source squares without advancing state.
    /// Warm work is O(n²+nm); no matrix exponential or residual solve occurs.
    pub fn evaluate_controlled(
        &self,
        initial: &[f64],
        held_inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        let mut coordinates = if self.root.is_empty() {
            Vec::new()
        } else {
            vec![0.0; 2 * self.n]
        };
        self.evaluate_into_controlled(initial, held_inputs, &mut coordinates, cancelled)
    }

    /// The identical positive source integral with caller-owned `2*n` scratch
    /// (unused and permitted empty when the source is identically zero).
    /// No allocation or exponential occurs during a successful warm call.
    /// Initial source and acceleration contractions use compensated products
    /// and sums, including state-only sources and explicit zero feedthrough.
    /// Scratch can change on refusal; no physical state is retained here.
    pub fn evaluate_into_controlled(
        &self,
        initial: &[f64],
        held_inputs: &[f64],
        coordinates: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        self.initialize_coordinates(initial, held_inputs, coordinates, &mut cancelled)?;
        if self.root.is_empty() {
            return Ok(0.0);
        }
        let mut loss = 0.0;
        for row in self.root.chunks_exact(2 * self.n) {
            poll(&mut cancelled)?;
            let output = row
                .iter()
                .zip(coordinates.iter())
                .map(|(a, x)| a * x)
                .sum::<f64>();
            loss += output * output;
        }
        if !loss.is_finite() {
            return Err(PhsError::Dimension {
                what: "source loss integral is unrepresentable",
            }
            .into());
        }
        Ok(loss)
    }
    /// Signed integral of one ORIGINAL source row `S*x(t)+D*u`, without
    /// squaring it or subtracting large state and source averages.
    ///
    /// `source_row` follows the rows supplied at preparation. The unchanged
    /// derivative coordinates retain a small source difference about a driven
    /// equilibrium. Their first n basis trajectories have constant physical
    /// output sqrt(rate), so the existing positive Gramian factor R gives
    /// `integral source_i dt = (R^T R z)_i / sqrt(rate)` directly. Initialization
    /// and both factor contractions use compensated products/sums, including
    /// state-only source differences. This is an independent linear observation,
    /// never a heat or energy-balance remainder.
    ///
    /// No new factor, cache, dimension or allocation is introduced. Caller
    /// scratch is the same `2*n` coordinates as [`Self::evaluate_into_controlled`].
    /// Its contents may change on refusal; no physical state is retained.
    pub fn signed_integral_into_controlled(
        &self,
        source_row: usize,
        initial: &[f64],
        held_inputs: &[f64],
        coordinates: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        poll(&mut cancelled)?;
        if source_row >= self.source_count {
            return Err(refuse("signed source integral requires an original source row").into());
        }
        self.initialize_coordinates(initial, held_inputs, coordinates, &mut cancelled)?;
        if self.root.is_empty() {
            return Ok(0.0);
        }
        let root_rate = det::sqrt(self.rate);
        let (mut high, mut low) = (0.0, 0.0);
        for row in self.root.chunks_exact(2 * self.n) {
            poll(&mut cancelled)?;
            // An unrelated source can overflow even though this selected
            // linear observation is finite. Its exact zero coefficient must
            // exclude it before evaluating the unrelated factor row.
            if row[source_row] == 0.0 {
                continue;
            }
            let output = compensated_dot(row, coordinates, &[], &[]);
            let (product, product_error) = two_prod(row[source_row] / root_rate, output);
            let (next, sum_error) = two_sum(high, product);
            high = next;
            low += product_error + sum_error;
        }
        let integral = high + low;
        if !integral.is_finite() {
            return Err(refuse("signed source integral is unrepresentable").into());
        }
        Ok(integral)
    }

    // A conservative exact witness for a constant held-input trajectory.
    // A rounded compensated zero is insufficient: an underflowing derivative
    // can accumulate a representable change over a long interval. Only accept
    // products and additions proved exact in the ORIGINAL A*x+B*u reduction.
    pub(crate) fn has_exactly_zero_derivative_controlled(
        &self,
        initial: &[f64],
        held_inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<bool, PreparedStepError> {
        poll(&mut cancelled)?;
        if initial.len() != self.n
            || held_inputs.len() != self.m
            || !finite(initial)
            || !finite(held_inputs)
        {
            return Err(
                refuse("stationarity requires complete finite state and held inputs").into(),
            );
        }
        // This conservative range keeps all 106 significand bits of an exact
        // f64 product representable, so a zero TwoProd error cannot hide an
        // underflowed remainder. It is an arithmetic witness, not a power floor.
        let exact_product_floor = f64::MIN_POSITIVE / (f64::EPSILON * f64::EPSILON);
        for i in 0..self.n {
            poll(&mut cancelled)?;
            let state_terms = self.flow[i * self.n..(i + 1) * self.n].iter().zip(initial);
            let input_terms = self.input[i * self.m..(i + 1) * self.m]
                .iter()
                .zip(held_inputs);
            let mut sum = 0.0;
            for (&coefficient, &value) in state_terms.chain(input_terms) {
                if coefficient == 0.0 || value == 0.0 {
                    continue;
                }
                let (product, product_error) = two_prod(coefficient, value);
                if !product.is_finite()
                    || product.abs() < exact_product_floor
                    || product_error != 0.0
                {
                    return Ok(false);
                }
                let (next, sum_error) = two_sum(sum, product);
                if !next.is_finite() || sum_error != 0.0 {
                    return Ok(false);
                }
                sum = next;
            }
            if sum != 0.0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn initialize_coordinates(
        &self,
        initial: &[f64],
        held_inputs: &[f64],
        coordinates: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PreparedStepError> {
        poll(cancelled)?;
        let n = self.n;
        if initial.len() != n
            || held_inputs.len() != self.m
            || (!self.root.is_empty() && coordinates.len() != 2 * n)
            || !finite(initial)
            || !finite(held_inputs)
        {
            return Err(refuse(
                "source loss requires complete finite initial state and held inputs",
            )
            .into());
        }
        if self.root.is_empty() {
            return Ok(());
        }
        let root_rate = det::sqrt(self.rate);
        for i in 0..n {
            poll(cancelled)?;
            // A small dissipative strain or acceleration can lie between
            // large, exactly cancelling lossless coordinates even with D=0.
            // Use the same compensated contraction for every source and its
            // complete derivative BEFORE the positive integral squares it.
            let direct = if self.feedthrough.is_empty() {
                &[][..]
            } else {
                &self.feedthrough[i * self.m..(i + 1) * self.m]
            };
            coordinates[i] = compensated_dot(
                &self.source[i * n..(i + 1) * n],
                initial,
                direct,
                held_inputs,
            ) / root_rate;
            coordinates[n + i] = compensated_dot(
                &self.flow[i * n..(i + 1) * n],
                initial,
                &self.input[i * self.m..(i + 1) * self.m],
                held_inputs,
            ) / self.rate;
        }
        if !finite(coordinates) {
            return Err(
                refuse("source loss initial output or derivative is unrepresentable").into(),
            );
        }
        Ok(())
    }
}

// The same linear form spans state and source coordinates. Separate rounded
// dot products can erase a weak loss or acceleration about a driven equilibrium.
fn compensated_dot(a: &[f64], x: &[f64], b: &[f64], u: &[f64]) -> f64 {
    let mut high = 0.0;
    let mut low = 0.0;
    for (&coefficient, &value) in a.iter().zip(x).chain(b.iter().zip(u)) {
        let (product, product_error) = two_prod(coefficient, value);
        let (next, sum_error) = two_sum(high, product);
        high = next;
        low += product_error + sum_error;
    }
    high + low
}
