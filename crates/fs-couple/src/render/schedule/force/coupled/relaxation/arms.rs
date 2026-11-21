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

/// [`arm_chains`] on `threads` threads (the calling thread among them): each takes a contiguous run of whole blocks of
/// eight arms (the last run also the arms after the last whole block), and the chains are returned in arm order. Each
/// chain is computed by one thread exactly as [`arm_chains`] computes it, so it has its bits for every thread count.
pub(crate) fn arm_chains_on<const STRIDE: usize>(
    arms: &[(f64, f64, &[f64])],
    q: &[f64],
    threads: usize,
) -> Vec<f64> {
    let per = (arms.len().div_ceil(ARMS).div_ceil(threads.max(1)) * ARMS).max(1);
    if threads <= 1 || per >= arms.len() {
        return arm_chains::<STRIDE>(arms, q);
    }
    let mut parts = arms.chunks(per);
    let first = parts.next().unwrap_or_default();
    // The scope joins every thread it started before it returns, and panics if one of them did.
    std::thread::scope(|scope| {
        let rest: Vec<_> = parts
            .map(|part| scope.spawn(move || arm_chains::<STRIDE>(part, q)))
            .collect();
        let mut chains = arm_chains::<STRIDE>(first, q);
        for part in rest {
            chains.extend(joined(part));
        }
        chains
    })
}

/// The value of a scoped thread, or its panic raised again on the joining thread.
fn joined<T>(handle: std::thread::ScopedJoinHandle<'_, T>) -> T {
    handle
        .join()
        .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
}

/// Stage each arm `(sqrt(k), b)` of `arms` with kept strain `memory` at the mechanical state `endpoint`, laid out
/// `[q_0, v_0, q_1, v_1, ..]`, on `threads` threads: its mismatch `sqrt(k) b.endpoint - z`; the arms' storage, the
/// sum of `m^2 / 2` in arm order, which is returned; and its restoring load `sqrt(k) b_i m`, taken from `forces[i]`
/// in arm order. `finite`, which refuses exactly the values that are not finite, checks each mismatch and each
/// partial storage, as the loop this replaced did; the forces are checked as the last paragraph says.
///
/// The strain adds only the displacement terms, and nothing checks its terms. That is the loop's arithmetic for
/// projections whose velocity entries are +-0 and an `endpoint` whose velocities are finite:
/// - a velocity term `sqrt(k) * (+-0) * v` is then +-0, and a chain from +0.0 is never -0, so adding it changes no bit;
/// - a chain that leaves the finite range never returns to it, and its mismatch is then not finite either, so the
///   mismatch's check refuses each arm that a check of each term refused, with the forces of the arms before it
///   staged and none of its own.
///
/// The strains are shared out by blocks of arms ([`arm_chains_on`]); the mismatches and the storage follow in arm
/// order on the calling thread, up to the first that refuses; then the forces of the arms before it are shared out by
/// runs of force entries, each entry taking every arm's load in arm order. A force that leaves the finite range never
/// returns to it, so the entries are checked once, at the end; if one is not finite the forces are put back and staged
/// again in the loop's order, which refuses where the loop refused, with the same forces staged. So every value, every
/// refusal and every staged force is the loop's for every thread count.
pub(crate) fn staged_arms<E>(
    arms: &[(f64, &[f64])],
    memory: &[f64],
    endpoint: &[f64],
    forces: &mut [f64],
    threads: usize,
    finite: impl Fn(f64) -> Result<f64, E>,
) -> Result<f64, E> {
    let q: Vec<_> = endpoint.iter().step_by(2).copied().collect();
    let chains: Vec<_> = arms
        .iter()
        .map(|&(root, projection)| (0.0, root, projection))
        .collect();
    let (mut energy, mut mismatches, mut refusal) = (0.0, Vec::with_capacity(arms.len()), None);
    let strains = arm_chains_on::<2>(&chains, &q, threads);
    for (strain, z) in strains.into_iter().zip(memory) {
        let staged = finite(strain - z).and_then(|m| finite(energy + 0.5 * m * m).map(|e| (m, e)));
        match staged {
            Ok((mismatch, partial)) => {
                mismatches.push(mismatch);
                energy = partial;
            }
            Err(error) => {
                refusal = Some(error);
                break;
            }
        }
    }
    let staging = &arms[..mismatches.len()];
    let before = forces.to_vec();
    if !loads_on(staging, &mismatches, forces, threads) {
        forces.copy_from_slice(&before);
        for (&(root, projection), mismatch) in staging.iter().zip(&mismatches) {
            for (force, &b) in forces.iter_mut().zip(projection.iter().step_by(2)) {
                *force = finite(*force - root * b * mismatch)?;
            }
        }
    }
    refusal.map_or(Ok(energy), Err)
}

/// Take each arm's load `sqrt(k) b_i m` from `forces[i]` in arm order, the entries shared out by contiguous runs to
/// `threads` threads, with no check; whether every force is still finite.
fn loads_on(
    arms: &[(f64, &[f64])],
    mismatches: &[f64],
    forces: &mut [f64],
    threads: usize,
) -> bool {
    let run = |start: usize, part: &mut [f64]| {
        for (&(root, projection), mismatch) in arms.iter().zip(mismatches) {
            let entries = projection.iter().step_by(2).skip(start);
            for (force, &b) in part.iter_mut().zip(entries) {
                *force -= root * b * mismatch;
            }
        }
        part.iter().all(|force| force.is_finite())
    };
    let per = forces.len().div_ceil(threads.max(1)).max(1);
    if threads <= 1 || per >= forces.len() {
        return run(0, forces);
    }
    let mut parts = (0..).step_by(per).zip(forces.chunks_mut(per));
    let first = parts.next();
    std::thread::scope(|scope| {
        let rest: Vec<_> = parts
            .map(|(at, part)| scope.spawn(move || run(at, part)))
            .collect();
        let mine = first.is_none_or(|(start, part)| run(start, part));
        rest.into_iter().map(joined).fold(mine, |all, ok| all & ok)
    })
}
