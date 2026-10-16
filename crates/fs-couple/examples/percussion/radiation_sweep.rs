//! Cold, ordered frequency parallelism. No solver arithmetic is parallelized.
//! Each worker owns its dynamic quadrature and dense solve; only immutable
//! geometry is shared. The highest frequency always runs alone, first.
use super::Error;

const THREAD_ENV: &str = "FRANKENSIM_BEM_THREADS";
const MAX_THREADS: usize = 64;

fn parse_threads(value: Option<&str>) -> Result<usize, Error> {
    let Some(value) = value else { return Ok(1); };
    let threads = value.parse::<usize>().ok().filter(|n| (1..=MAX_THREADS).contains(n));
    threads.ok_or_else(|| format!("{THREAD_ENV} must be an integer in 1..={MAX_THREADS}").into())
}

pub(super) fn threads_from_env() -> Result<usize, Error> {
    match std::env::var(THREAD_ENV) {
        Ok(value) => parse_threads(Some(&value)),
        Err(std::env::VarError::NotPresent) => parse_threads(None),
        Err(std::env::VarError::NotUnicode(_)) => Err(format!("{THREAD_ENV} must be Unicode integer text").into()),
    }
}

/// Visit samples in descending index order. `solve` can run concurrently, but
/// `accept` runs only on the caller, in serial order. With one worker neither
/// threads nor string conversions of the example's non-Send Error are used.
///
/// On failure the caller must discard its cold output buffers. Every started
/// worker is joined, including on a refusal, a panic, or an OS spawn failure.
/// Concurrent scratch storage is bounded by `threads` frequency workspaces;
/// the caller's panel/source limits still apply to EACH workspace.
pub(super) fn run<T, F, A>(samples: usize, threads: usize, solve: F, mut accept: A)
    -> Result<(), Error>
where
    T: Send,
    F: Fn(usize) -> Result<T, Error> + Sync,
    A: FnMut(usize, T),
{
    if samples == 0 || !(1..=MAX_THREADS).contains(&threads) {
        return Err("radiation sweep needs samples and 1..=64 threads".into());
    }
    // The owner's highest-frequency wavelength refusal must precede ALL lower
    // samples, even when a caller explicitly requests parallel preparation.
    let highest = samples - 1;
    accept(highest, solve(highest)?);
    if threads == 1 || highest == 0 {
        for index in (0..highest).rev() { accept(index, solve(index)?); }
        return Ok(());
    }
    let indices: Vec<_> = (0..highest).rev().collect();
    let workers = threads.min(indices.len());
    let chunk_size = indices.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        let mut spawn_error = None;
        let solve = &solve;
        for chunk in indices.chunks(chunk_size) {
            let worker = std::thread::Builder::new().spawn_scoped(scope, move || {
                let mut rows = Vec::with_capacity(chunk.len());
                for &index in chunk {
                    // Box<dyn Error> is intentionally not Send. Preserve its
                    // message, rather than imposing Send on every example error.
                    let row = solve(index).map_err(|error| error.to_string())?;
                    rows.push((index, row));
                }
                Ok::<_, String>(rows)
            });
            match worker {
                Ok(handle) => handles.push(handle),
                Err(error) => { spawn_error = Some(error); break; }
            }
        }
        let mut failure: Option<Error> = None;
        // Handles follow contiguous descending chunks, not completion order.
        // Do not early-return: an unjoined panicking worker would make scope
        // unwind, and no borrowed geometry may escape the cold preparation.
        for handle in handles {
            let result = handle.join().unwrap_or_else(|_| {
                Err("radiation frequency worker panicked".to_owned())
            });
            if failure.is_none() {
                match result {
                    Ok(rows) => for (index, row) in rows { accept(index, row); },
                    Err(error) => failure = Some(error.into()),
                }
            }
        }
        if let Some(error) = failure { return Err(error); }
        if let Some(error) = spawn_error {
            return Err(format!("could not start radiation frequency worker: {error}").into());
        }
        Ok(())
    })
}

#[cfg(test)]
#[path = "radiation_sweep/tests.rs"]
mod tests;
