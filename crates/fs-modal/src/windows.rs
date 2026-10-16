//! Ordered parallel windows on the existing certified slicing owner.
//! Endpoints are factored once, windows share one immutable symbolic pencil
//! and mass factor, and every window retains its serial arithmetic/seed.

use super::*;
use fs_exec::Cx;
use std::time::Instant;

/// All resource decisions are explicit. Adjacent or disjoint ordered windows
/// are accepted; overlaps and repeated windows are refused.
#[derive(Clone, Debug)]
pub struct WindowPlan {
    /// Windows in the pencil's eigenvalue units; exact endpoints may refuse.
    pub windows: Vec<(f64, f64)>,
    /// Requested scoped worker count (wasm uses one).
    pub threads: usize,
    /// Maximum admitted number of windows.
    pub maximum_windows: usize,
    /// Numerical payload admission, not an OS RSS/allocator-overhead claim.
    /// Includes retained inputs/symbolic/mass, simultaneous factors, Krylov
    /// vectors, certificates and outputs. Symbolic construction is transient
    /// and not yet cancellation/memory bounded inside fs-sparse.
    pub maximum_bytes: usize,
    /// Required residual / |lambda| for every returned mode.
    pub maximum_relative_halfwidth: f64,
    /// Optional xyz per row, used by geometric ND only.
    pub coordinates: Option<Vec<[f64; 3]>>,
}

/// Work accounting; elapsed time is excluded from deterministic result bits.
#[derive(Clone, Debug)]
pub struct SliceTimingReceipt {
    /// Wall seconds for the whole admitted call.
    pub wall_s: f64,
    /// Effective workers (one on wasm).
    pub threads: usize,
    /// Shared shifted-pencil analyses, always one on success.
    pub pencil_symbolic_analyses: usize,
    /// Separate mass analyses for the SPD gate/certificate norm.
    pub mass_symbolic_analyses: usize,
    /// Actual unique endpoint factorization calls.
    pub endpoint_factorizations: usize,
    /// Interior attempts, including singular-pole fallbacks.
    pub interior_factorizations: usize,
    /// Conservative admitted numerical payload, excluding process overhead.
    pub admitted_peak_bytes: usize,
    /// Largest residual / |lambda| among returned modes.
    pub maximum_relative_halfwidth: f64,
}

/// Complete, window-ordered results; failure returns no partial report.
#[derive(Clone, Debug)]
pub struct MultiSliceReport {
    /// One report per requested window, in caller order.
    pub reports: Vec<SliceReport>,
    /// Shared symbolic pencil's L nonzero count.
    pub symbolic_nnz_l: usize,
    /// Non-numerical work/timing accounting.
    pub receipt: SliceTimingReceipt,
}

fn check(cx: &Cx<'_>) -> Result<(), ModalError> {
    cx.checkpoint().map_err(|_| ModalError::Cancelled)
}
fn admit(bytes: usize, cap: usize) -> Result<(), ModalError> {
    if bytes > cap {
        Err(ModalError::MemoryBudget {
            estimated_bytes: bytes,
            cap,
        })
    } else {
        Ok(())
    }
}
fn csr_bytes(a: &Csr) -> usize {
    a.nnz()
        .saturating_mul(8 + core::mem::size_of::<usize>())
        .saturating_add((a.nrows() + 1).saturating_mul(core::mem::size_of::<usize>()))
}

fn endpoint(v: f64) -> f64 {
    if v == 0.0 { 0.0 } else { v }
}

/// Return every requested slice or a refusal after executor drain. The shared
/// endpoints use the serial owner's nonsingular-endpoint rule: exact spectral
/// endpoints refuse, never get nudged into a made-up count. A factorization
/// and symbolic analysis are indivisible cancellation units; polls occur
/// around them, at Lanczos restarts and before publishing any report.
pub fn slice_windows(
    k: &Csr,
    m: &Csr,
    plan: &WindowPlan,
    options: &SliceOptions,
    cx: &Cx<'_>,
) -> Result<MultiSliceReport, ModalError> {
    let start = Instant::now();
    check(cx)?;
    let n = check_shapes(k, m)?;
    if plan.windows.len() > plan.maximum_windows {
        return Err(ModalError::WindowBudget {
            requested: plan.windows.len(),
            cap: plan.maximum_windows,
        });
    }
    if plan.windows.is_empty()
        || plan.threads == 0
        || plan.maximum_bytes == 0
        || options.ldlt.max_factor_bytes.is_none_or(|cap| cap == 0)
        || !plan.maximum_relative_halfwidth.is_finite()
        || plan.maximum_relative_halfwidth <= 0.0
        || !options.ritz_tol.is_finite()
        || options.ritz_tol <= 0.0
    {
        return Err(ModalError::WindowPlan {
            field: "nonempty windows, positive threads/bytes/explicit factor cap/finite tolerances",
        });
    }
    for (i, &(lo, hi)) in plan.windows.iter().enumerate() {
        if !(lo.is_finite() && hi.is_finite() && lo < hi) {
            return Err(ModalError::InvalidWindow { low: lo, high: hi });
        }
        if i > 0 && lo < plan.windows[i - 1].1 {
            return Err(ModalError::WindowPlan {
                field: "ordered nonoverlapping windows",
            });
        }
    }
    let threads = if cfg!(target_arch = "wasm32") {
        1
    } else {
        plan.threads
    };
    let mut opts = *options;
    opts.ldlt.max_factor_bytes = Some(
        opts.ldlt
            .max_factor_bytes
            .unwrap_or(plan.maximum_bytes)
            .min(plan.maximum_bytes),
    );
    let input = csr_bytes(k).saturating_add(csr_bytes(m)).saturating_add(
        plan.coordinates
            .as_ref()
            .map_or(0, |p| p.len().saturating_mul(24)),
    );
    admit(input, plan.maximum_bytes)?;
    admit(
        input.saturating_add(opts.ldlt.max_factor_bytes.expect("explicit factor cap")),
        plan.maximum_bytes,
    )?;
    let prepared = PreparedPencil::new(
        k,
        m,
        plan.windows[0].0,
        &opts,
        plan.coordinates.as_deref(),
        &mut || check(cx),
    )?;
    check(cx)?;
    let (wk, wm) = prepared.work();
    let shared = input
        .saturating_add(if prepared.scales.is_some() {
            n.saturating_mul(8)
        } else {
            0
        })
        .saturating_add(prepared.p_sym.symbolic_bytes())
        .saturating_add(prepared.m_fact.stats().peak_factor_bytes)
        .saturating_add(if prepared.scaled.is_some() {
            csr_bytes(wk).saturating_add(csr_bytes(wm))
        } else {
            0
        });
    let factor_cap = opts
        .ldlt
        .max_factor_bytes
        .expect("validated explicit factor cap");
    // Endpoint count pass has one numeric factor and one assembled shift live.
    let endpoint_peak = shared.saturating_add(factor_cap).saturating_add(
        csr_bytes(wk)
            .saturating_add(csr_bytes(wm))
            .saturating_mul(3),
    );
    admit(endpoint_peak, plan.maximum_bytes)?;
    let mut endpoints: Vec<_> = plan
        .windows
        .iter()
        .flat_map(|&(a, b)| [endpoint(a), endpoint(b)])
        .collect();
    endpoints.sort_by(f64::total_cmp);
    endpoints.dedup_by(|a, b| a.to_bits() == b.to_bits());
    let mut counts = Vec::with_capacity(endpoints.len());
    for &shift in &endpoints {
        check(cx)?;
        let a = shifted_pencil(wk, wm, shift);
        let f = prepared
            .p_sym
            .factor(&a, &opts.ldlt)
            .map_err(|source| ModalError::Factor { shift, source })?;
        admit(
            shared
                .saturating_add(f.stats().peak_factor_bytes)
                .saturating_add(csr_bytes(&a)),
            plan.maximum_bytes,
        )?;
        counts.push(f.inertia().negative);
    }
    check(cx)?;
    let window_counts: Vec<_> = plan
        .windows
        .iter()
        .map(|&(lo, hi)| {
            let l = endpoints
                .binary_search_by(|x| x.total_cmp(&endpoint(lo)))
                .expect("declared endpoint");
            let h = endpoints
                .binary_search_by(|x| x.total_cmp(&endpoint(hi)))
                .expect("declared endpoint");
            (counts[l], counts[h])
        })
        .collect();
    let mut work = Vec::new();
    let mut output = 0usize;
    for &(lo, hi) in &window_counts {
        let expected = hi.checked_sub(lo).ok_or(ModalError::WindowPlan {
            field: "inertia counts decreased",
        })?;
        let sweep = (20usize.saturating_add(6usize.saturating_mul(expected))).min(
            if opts.max_lanczos == 0 {
                n
            } else {
                opts.max_lanczos.min(n)
            },
        );
        let vectors = 2usize
            .saturating_mul(sweep)
            .saturating_add(3usize.saturating_mul(expected))
            .saturating_add(16);
        let bytes = factor_cap
            .saturating_add(n.saturating_mul(vectors).saturating_mul(8))
            .saturating_add(sweep.saturating_mul(sweep).saturating_mul(32))
            .saturating_add(
                csr_bytes(wk)
                    .saturating_add(csr_bytes(wm))
                    .saturating_mul(3),
            );
        work.push(bytes);
        output = output.saturating_add(n.saturating_mul(expected).saturating_mul(8));
    }
    work.sort_unstable_by(|a, b| b.cmp(a));
    let peak = shared
        .saturating_add(output)
        .saturating_add(
            work.iter()
                .take(threads)
                .fold(0usize, |a, &b| a.saturating_add(b)),
        )
        .max(endpoint_peak);
    admit(peak, plan.maximum_bytes)?;
    let reports = execute(&prepared, plan, &opts, &window_counts, cx, threads)?;
    check(cx)?;
    let mut max_width = 0.0f64;
    for report in &reports {
        for mode in &report.modes {
            let relative = mode.residual / mode.lambda.abs().max(f64::MIN_POSITIVE);
            max_width = max_width.max(relative);
            if !relative.is_finite() || relative > plan.maximum_relative_halfwidth {
                return Err(ModalError::WindowUnresolved {
                    expected: report.expected,
                    converged: report
                        .modes
                        .iter()
                        .filter(|m| {
                            m.residual / m.lambda.abs().max(f64::MIN_POSITIVE)
                                <= plan.maximum_relative_halfwidth
                        })
                        .count(),
                    window: report.window,
                });
            }
        }
    }
    let interior = reports.iter().map(|r| r.stats.factorizations).sum();
    Ok(MultiSliceReport {
        symbolic_nnz_l: prepared.p_sym.nnz_l(),
        reports,
        receipt: SliceTimingReceipt {
            wall_s: start.elapsed().as_secs_f64(),
            threads,
            pencil_symbolic_analyses: 1,
            mass_symbolic_analyses: 1,
            endpoint_factorizations: endpoints.len(),
            interior_factorizations: interior,
            admitted_peak_bytes: peak,
            maximum_relative_halfwidth: max_width,
        },
    })
}

#[cfg(target_arch = "wasm32")]
fn execute(
    p: &PreparedPencil<'_>,
    plan: &WindowPlan,
    opts: &SliceOptions,
    counts: &[(usize, usize)],
    cx: &Cx<'_>,
    _: usize,
) -> Result<Vec<SliceReport>, ModalError> {
    plan.windows
        .iter()
        .zip(counts)
        .map(|(&w, &c)| p.slice(w, opts, Some(c), &mut || check(cx)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cancel_after_interior_factor_before_lanczos_returns_no_slice() {
        let mut k = Coo::new(8, 8);
        for i in 0..8 {
            k.push(i, i, i as f64 + 1.0);
        }
        let k = k.assemble();
        let m = Csr::identity(8);
        let opts = SliceOptions::default();
        let p = PreparedPencil::new(&k, &m, 0.1, &opts, None, &mut || Ok(())).unwrap();
        let mut polls = 0;
        let result = p.slice((0.1, 3.5), &opts, Some((0, 3)), &mut || {
            polls += 1;
            if polls == 2 {
                Err(ModalError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(result, Err(ModalError::Cancelled)));
        assert_eq!(polls, 2);
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn execute(
    p: &PreparedPencil<'_>,
    plan: &WindowPlan,
    opts: &SliceOptions,
    counts: &[(usize, usize)],
    cx: &Cx<'_>,
    threads: usize,
) -> Result<Vec<SliceReport>, ModalError> {
    use core::ops::ControlFlow;
    use fs_exec::{Cancelled, RunId, TileKernel, TilePlan, TilePool};
    use std::sync::Mutex;
    struct Windows<'a, 'b> {
        prepared: &'a PreparedPencil<'b>,
        plan: &'a WindowPlan,
        opts: &'a SliceOptions,
        counts: &'a [(usize, usize)],
        slots: Vec<Mutex<Option<Result<SliceReport, ModalError>>>>,
    }
    impl TileKernel for Windows<'_, '_> {
        type Out = u64;
        fn tiles(&self) -> TilePlan {
            TilePlan::new("fs-modal/spectrum-windows", self.plan.windows.len() as u64)
        }
        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let i = tile as usize;
            let result = self.prepared.slice(
                self.plan.windows[i],
                self.opts,
                Some(self.counts[i]),
                &mut || check(cx),
            );
            if matches!(result, Err(ModalError::Cancelled)) {
                return ControlFlow::Break(Cancelled);
            }
            match self.slots[i].lock() {
                Ok(mut slot) => *slot = Some(result),
                Err(_) => return ControlFlow::Break(Cancelled),
            }
            ControlFlow::Continue(1)
        }
    }
    let kernel = Windows {
        prepared: p,
        plan,
        opts,
        counts,
        slots: (0..plan.windows.len()).map(|_| Mutex::new(None)).collect(),
    };
    let pool = TilePool::for_host(threads, cx.stream_key().seed);
    // Existing executor bridge carries the ORIGINAL Cx cancellation authority
    // into scoped workers, and drains before returning. No polling thread or
    // shared arena/Cx is invented. Numeric Vec payload has separate admission.
    let lease = fs_alloc::OperationMemoryLease::bounded(plan.maximum_bytes as u64);
    let (run, _) =
        pool.run_declared_leased_with_cx(cx, &kernel, RunId(cx.stream_key().kernel_id), &lease);
    check(cx)?;
    run.map_err(|e| ModalError::Executor {
        detail: format!("{e:?}"),
    })?;
    kernel
        .slots
        .into_iter()
        .map(|slot| {
            slot.into_inner()
                .map_err(|_| ModalError::Executor {
                    detail: "poisoned window slot".into(),
                })?
                .ok_or(ModalError::Executor {
                    detail: "missing drained window".into(),
                })?
        })
        .collect()
}
