//! The Maxwell arms' chains, eight arms at a time.
//!
//! An arm's strain is the chain `sqrt(k) b_0 q_0 + sqrt(k) b_1 q_1 + ..` over its projection's displacement entries,
//! from +0.0 in mode order; the prepared operator's mismatch is the chain `-z + c_0 q_0 + c_1 q_1 + ..`. Each chain
//! keeps its own order, so it has the bits of that arm's own fold. Eight arms' chains advance together, so one arm's
//! additions do not wait for another's; the arms after the last whole eight are each a fold of their own.
//!
//! This file names nothing of its crate. `tests/material_arms.rs` compiles it a second time and holds each value to
//! the plain loop it replaced, bit for bit.

/// Arms whose chains advance together.
const ARMS: usize = 8;

/// For each `(start, scale, row)` of `arms`, the chain `start + scale * row[STRIDE * i] * q[i] + ..` added in index
/// order over the indices that both `row` and `q` reach: the bits of that arm's own fold.
pub(crate) fn arm_chains<const STRIDE: usize>(arms: &[(f64, f64, &[f64])], q: &[f64]) -> Vec<f64> {
    let fold = |&(start, scale, row): &(f64, f64, &[f64])| {
        let terms = row.iter().step_by(STRIDE).zip(q);
        terms.fold(start, |sum, (b, x)| sum + scale * b * x)
    };
    let reach = |row: &[f64]| row.len().div_ceil(STRIDE);
    let (blocks, rest) = arms.as_chunks::<ARMS>();
    let mut chains = Vec::with_capacity(arms.len());
    for block in blocks {
        // A whole block's rows and scales are held in arrays of eight, so the loop over arms has a fixed length.
        let together = block
            .iter()
            .fold(q.len(), |n, (_, _, row)| n.min(reach(row)));
        let rows: [&[f64]; ARMS] = core::array::from_fn(|a| block[a].2);
        let scales: [f64; ARMS] = core::array::from_fn(|a| block[a].1);
        let mut sums: [f64; ARMS] = core::array::from_fn(|a| block[a].0);
        for (i, x) in q[..together].iter().enumerate() {
            for ((sum, scale), row) in sums.iter_mut().zip(&scales).zip(&rows) {
                *sum += scale * row[STRIDE * i] * x;
            }
        }
        for ((sum, scale), row) in sums.iter_mut().zip(&scales).zip(&rows) {
            for (b, x) in row.iter().step_by(STRIDE).zip(q).skip(together) {
                *sum += scale * b * x;
            }
        }
        chains.extend_from_slice(&sums);
    }
    chains.extend(rest.iter().map(fold));
    chains
}

/// Stage each arm `(sqrt(k), b)` of `arms` with kept strain `memory` at the mechanical state `endpoint`, laid out
/// `[q_0, v_0, q_1, v_1, ..]`: its mismatch `sqrt(k) b.endpoint - z`; the arms' storage, the sum of `m^2 / 2` in arm
/// order, which is returned; and its restoring load `sqrt(k) b_i m`, taken from `forces[i]` in arm order. `finite`
/// checks each mismatch, each partial storage and each force, as the loop this replaced did.
///
/// The strain adds only the displacement terms, and nothing checks its terms. That is the loop's arithmetic for
/// projections whose velocity entries are +-0 and an `endpoint` whose velocities are finite:
/// - a velocity term `sqrt(k) * (+-0) * v` is then +-0, and a chain from +0.0 is never -0, so adding it changes no bit;
/// - a chain that leaves the finite range never returns to it, and its mismatch is then not finite either, so the
///   mismatch's check refuses each arm that a check of each term refused, with the forces of the arms before it
///   staged and none of its own.
pub(crate) fn staged_arms<E>(
    arms: &[(f64, &[f64])],
    memory: &[f64],
    endpoint: &[f64],
    forces: &mut [f64],
    finite: impl Fn(f64) -> Result<f64, E>,
) -> Result<f64, E> {
    let q: Vec<_> = endpoint.iter().step_by(2).copied().collect();
    let chains: Vec<_> = arms
        .iter()
        .map(|&(root, projection)| (0.0, root, projection))
        .collect();
    let mut energy = 0.0;
    for ((&(root, projection), z), strain) in
        arms.iter().zip(memory).zip(arm_chains::<2>(&chains, &q))
    {
        let mismatch = finite(strain - z)?;
        energy = finite(energy + 0.5 * mismatch * mismatch)?;
        for (force, &b) in forces.iter_mut().zip(projection.iter().step_by(2)) {
            *force = finite(*force - root * b * mismatch)?;
        }
    }
    Ok(energy)
}
