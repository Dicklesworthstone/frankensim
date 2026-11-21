//! A fixed set of worker threads that share one pass over the rows of a matrix at a time.
//!
//! Each row's value is computed by exactly one thread, with the kernel the serial pass uses, so a pass has the same
//! bits for every worker count. The threads start once and wait between passes by spinning on a generation counter,
//! then yielding, then sleeping: a pass of a few hundred microseconds is too short for a thread start or a condition
//! variable's wake (`Instrument::share_head_load_workers`, the heads' radiation loads).

use std::ops::Range;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::thread::JoinHandle;
use std::time::Duration;

/// What a pass computes for row `row` of its matrix and its vector `x`.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Kind {
    /// `crate::dot(row, x)`.
    Dot,
    /// `init` plus each `row[i] * x[i]`, added in ascending `i` (`ordered_sums::sums_of_products_from`).
    Ordered(f64),
}

#[derive(Clone)]
struct Job {
    matrix: Arc<[f64]>,
    len: usize,
    x: Arc<[f64]>,
    kind: Kind,
}

struct Shared {
    threads: usize,
    job: Mutex<Option<Job>>,
    generation: AtomicU64,
    finished: AtomicUsize,
    stop: AtomicBool,
    failed: AtomicBool,
    /// One per thread; the calling thread (index 0) writes its rows straight into the output.
    slots: Vec<Mutex<Vec<f64>>>,
}

/// `threads` threads in all: the calling thread and `threads - 1` workers.
pub(crate) struct Pool {
    shared: Arc<Shared>,
    handles: Vec<JoinHandle<()>>,
    /// One pass at a time, whichever thread calls.
    dispatch: Mutex<()>,
}

impl std::fmt::Debug for Pool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Pool").field("threads", &self.shared.threads).finish()
    }
}

/// Spins this many times, then yields this many times, then sleeps between looks.
const SPINS: u32 = 1 << 14;
const YIELDS: u32 = 1 << 12;
const NAP: Duration = Duration::from_micros(50);

fn wait_until(done: impl Fn() -> bool) {
    let mut looks = 0u32;
    while !done() {
        if looks < SPINS {
            std::hint::spin_loop();
        } else if looks < SPINS + YIELDS {
            std::thread::yield_now();
        } else {
            std::thread::sleep(NAP);
        }
        looks = looks.saturating_add(1);
    }
}

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// Thread `t`'s rows of `rows`: whole blocks of eight, the same share for each thread, the last share shorter.
fn share(rows: usize, threads: usize, t: usize) -> Range<usize> {
    let per = rows.div_ceil(8).div_ceil(threads.max(1)) * 8;
    (t * per).min(rows)..((t + 1) * per).min(rows)
}

fn compute(job: &Job, rows: Range<usize>, out: &mut [f64]) {
    let (m, len, x) = (&job.matrix, job.len, &job.x[..]);
    let start = rows.start;
    match job.kind {
        Kind::Dot => {
            for (o, r) in out.iter_mut().zip(rows) {
                *o = crate::dot(&m[r * len..(r + 1) * len], x);
            }
        }
        Kind::Ordered(init) => {
            crate::ordered_sums::sums_of_products_from(init, |i| &m[(start + i) * len..(start + i + 1) * len], x, out);
        }
    }
}

fn work(shared: &Shared, t: usize) {
    let mut seen = 0u64;
    loop {
        wait_until(|| shared.generation.load(Ordering::Acquire) != seen);
        seen = shared.generation.load(Ordering::Acquire);
        if shared.stop.load(Ordering::Acquire) {
            return;
        }
        let job = lock(&shared.job).clone();
        if let Some(job) = job {
            let ran = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let rows = share(job.matrix.len() / job.len, shared.threads, t);
                let mut slot = lock(&shared.slots[t]);
                slot.clear();
                slot.resize(rows.len(), 0.0);
                compute(&job, rows, &mut slot);
            }));
            if ran.is_err() {
                shared.failed.store(true, Ordering::Release);
            }
        }
        shared.finished.fetch_add(1, Ordering::AcqRel);
    }
}

impl Pool {
    /// A pool of `threads` threads, the caller among them: `None` below two, or when a worker cannot be started.
    pub(crate) fn new(threads: usize) -> Option<Arc<Pool>> {
        if threads < 2 {
            return None;
        }
        let shared = Arc::new(Shared {
            threads,
            job: Mutex::new(None),
            generation: AtomicU64::new(0),
            finished: AtomicUsize::new(0),
            stop: AtomicBool::new(false),
            failed: AtomicBool::new(false),
            slots: (0..threads).map(|_| Mutex::new(Vec::new())).collect(),
        });
        let mut pool = Pool { shared, handles: Vec::with_capacity(threads - 1), dispatch: Mutex::new(()) };
        for t in 1..threads {
            let shared = Arc::clone(&pool.shared);
            let handle = std::thread::Builder::new()
                .name(format!("fs-vkplate-pass-{t}"))
                .spawn(move || work(&shared, t))
                .ok()?; // dropping `pool` stops the workers already started
            pool.handles.push(handle);
        }
        Some(Arc::new(pool))
    }

    pub(crate) fn threads(&self) -> usize {
        self.shared.threads
    }

    /// `out[r]` for each row `r` of `matrix` (rows of `len` entries, `out.len()` of them) and `x`, by `kind`.
    pub(crate) fn run(&self, matrix: &Arc<[f64]>, len: usize, x: &[f64], kind: Kind, out: &mut [f64]) {
        assert_eq!(matrix.len(), out.len() * len, "a pass's matrix has out.len() rows of len entries");
        if len == 0 {
            // rows of no entries: each kernel's value of an empty row, on this thread
            out.fill(match kind {
                Kind::Dot => crate::dot(&[], x),
                Kind::Ordered(init) => init,
            });
            return;
        }
        let _one_pass = lock(&self.dispatch);
        let shared = &*self.shared;
        let job = Job { matrix: Arc::clone(matrix), len, x: Arc::from(x), kind };
        *lock(&shared.job) = Some(job.clone());
        shared.finished.store(0, Ordering::Release);
        shared.generation.fetch_add(1, Ordering::AcqRel);
        let rows = out.len();
        let mine = share(rows, shared.threads, 0);
        compute(&job, mine.clone(), &mut out[mine]);
        wait_until(|| shared.finished.load(Ordering::Acquire) == shared.threads - 1);
        *lock(&shared.job) = None;
        assert!(!shared.failed.load(Ordering::Acquire), "a pass worker panicked");
        for t in 1..shared.threads {
            let theirs = share(rows, shared.threads, t);
            out[theirs].copy_from_slice(&lock(&shared.slots[t]));
        }
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.shared.stop.store(true, Ordering::Release);
        self.shared.generation.fetch_add(1, Ordering::AcqRel);
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: usize, len: usize, seed: u64) -> Vec<f64> {
        let mut s = seed | 1;
        (0..rows * len)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                let v = (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
                if s.is_multiple_of(17) { -0.0 } else { v * 10f64.powi((s % 7) as i32 - 3) }
            })
            .collect()
    }

    fn serial(m: &[f64], len: usize, x: &[f64], kind: Kind) -> Vec<f64> {
        let rows = m.len() / len;
        let mut out = vec![0.0; rows];
        compute(&Job { matrix: Arc::from(m), len, x: Arc::from(x), kind }, 0..rows, &mut out);
        out
    }

    #[test]
    fn every_thread_count_has_the_serial_bits() {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        for (rows, len) in [(1, 5), (3, 8), (7, 13), (8, 64), (9, 1), (37, 29), (347, 61)] {
            let m: Arc<[f64]> = Arc::from(matrix(rows, len, (rows * 131 + len) as u64));
            let x = matrix(1, len, 977 + len as u64);
            for kind in [Kind::Dot, Kind::Ordered(0.0), Kind::Ordered(-0.0), Kind::Ordered(1.25)] {
                let want = serial(&m, len, &x, kind);
                for threads in 2..=7 {
                    let pool = Pool::new(threads).expect("workers start");
                    let mut out = vec![f64::NAN; rows];
                    for _ in 0..3 {
                        pool.run(&m, len, &x, kind, &mut out);
                        assert_eq!(bits(&out), bits(&want), "{rows} x {len}, {kind:?}, {threads} threads");
                    }
                }
            }
        }
    }

    #[test]
    fn ordered_rows_add_from_init_in_ascending_order() {
        let m = matrix(5, 9, 41);
        let x = matrix(1, 9, 43);
        for init in [0.0, -0.0, 3.5] {
            let got = serial(&m, 9, &x, Kind::Ordered(init));
            for (r, g) in got.iter().enumerate() {
                let plain = m[r * 9..(r + 1) * 9].iter().zip(&x).fold(init, |acc, (a, b)| acc + a * b);
                assert_eq!(g.to_bits(), plain.to_bits(), "row {r} from {init}");
            }
        }
    }

    #[test]
    fn rows_of_no_entries_give_the_kernels_empty_values() {
        let pool = Pool::new(3).expect("workers start");
        let m: Arc<[f64]> = Arc::from(Vec::new());
        let mut out = vec![f64::NAN; 4];
        pool.run(&m, 0, &[], Kind::Ordered(2.0), &mut out);
        assert_eq!(out, vec![2.0; 4]);
        pool.run(&m, 0, &[], Kind::Dot, &mut out);
        assert_eq!(out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), vec![crate::dot(&[], &[]).to_bits(); 4]);
    }
}
