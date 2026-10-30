//! Positive source-square observation of an existing held-input trajectory.
use super::{
    DenseExponentialBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite, positive_heat, refuse,
};
use crate::{PhsError, PreparedStepError};
use fs_math::det;

/// Cached integral of `||S*x(t)||²` for the unchanged `x'=A*x+B*u`.
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
    flow: Vec<f64>,
    input: Vec<f64>,
    source: Vec<f64>,
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
        let mut source = source_rows.to_vec();
        source.resize(square, 0.0);
        let (rate, root) = if source.iter().all(|&value| value == 0.0) {
            (1.0, Vec::new())
        } else {
            positive_heat::source_integral(flow, &source, n, duration)?
        };
        poll(&mut cancelled)?;
        Ok(Self {
            n,
            m,
            flow: flow.to_vec(),
            input: input.to_vec(),
            source,
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        poll(&mut cancelled)?;
        let n = self.n;
        if initial.len() != n
            || held_inputs.len() != self.m
            || !finite(initial)
            || !finite(held_inputs)
        {
            return Err(refuse(
                "source loss requires complete finite initial state and held inputs",
            )
            .into());
        }
        if self.root.is_empty() {
            return Ok(0.0);
        }
        let mut coordinates = vec![0.0; 2 * n];
        let root_rate = det::sqrt(self.rate);
        for i in 0..n {
            poll(&mut cancelled)?;
            coordinates[i] = self.source[i * n..(i + 1) * n]
                .iter()
                .zip(initial)
                .map(|(a, x)| a * x)
                .sum::<f64>()
                / root_rate;
            let free = self.flow[i * n..(i + 1) * n]
                .iter()
                .zip(initial)
                .map(|(a, x)| a * x)
                .sum::<f64>();
            let forcing = self.input[i * self.m..(i + 1) * self.m]
                .iter()
                .zip(held_inputs)
                .map(|(b, u)| b * u)
                .sum::<f64>();
            coordinates[n + i] = (free + forcing) / self.rate;
        }
        if !finite(&coordinates) {
            return Err(
                refuse("source loss initial output or derivative is unrepresentable").into(),
            );
        }
        let mut loss = 0.0;
        for row in self.root.chunks_exact(2 * n) {
            poll(&mut cancelled)?;
            let output = row
                .iter()
                .zip(&coordinates)
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
}
