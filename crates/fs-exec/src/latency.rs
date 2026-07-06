//! The latency lane: asupersync's native async scheduling for orchestration
//! (plan §5.2) — ledger I/O, IR interpretation, progress streaming,
//! watchdogs. Milliseconds matter, throughput doesn't; this lane is what
//! keeps HELM conversational while the tile pool saturates the cores.
//!
//! fs-exec deliberately adds NO scheduling policy of its own here: the lane
//! is a thin, configured handle on the asupersync runtime, so its
//! cancellation semantics (request → drain → finalize through region state
//! machines) apply unmodified.

use core::fmt;

/// Structured lane-construction failure.
#[derive(Debug)]
pub struct LaneError {
    /// What the runtime builder reported.
    pub detail: String,
}

impl fmt::Display for LaneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "latency lane construction failed: {}; orchestration cannot start without its \
             runtime — check thread limits and runtime configuration",
            self.detail
        )
    }
}

impl core::error::Error for LaneError {}

/// The orchestration runtime handle.
pub struct LatencyLane {
    runtime: asupersync::runtime::Runtime,
}

impl LatencyLane {
    /// Build a lane with `threads` worker threads (clamped to at least 1).
    /// Keep this SMALL — one or two threads; compute belongs to the tile
    /// pool.
    ///
    /// # Errors
    /// [`LaneError`] when the asupersync runtime cannot be built.
    pub fn new(threads: usize) -> Result<Self, LaneError> {
        asupersync::runtime::RuntimeBuilder::new()
            .worker_threads(threads.max(1))
            .build()
            .map(|runtime| LatencyLane { runtime })
            .map_err(|e| LaneError {
                detail: format!("{e:?}"),
            })
    }

    /// Drive a future to completion on the lane (installs the ambient
    /// asupersync `Cx` for structured spawning inside).
    pub fn block_on<F: Future>(&self, future: F) -> F::Output {
        self.runtime.block_on(future)
    }

    /// The underlying asupersync runtime, for scope/spawn composition.
    #[must_use]
    pub fn runtime(&self) -> &asupersync::runtime::Runtime {
        &self.runtime
    }
}

impl fmt::Debug for LatencyLane {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LatencyLane").finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lane_builds_and_drives_futures() {
        let lane = LatencyLane::new(1).expect("lane");
        let out = lane.block_on(async { 41 + 1 });
        assert_eq!(out, 42);
    }

    #[test]
    fn lane_supports_structured_spawn_and_join() {
        // The canonical asupersync shape: spawn a ROOT task through the
        // runtime handle; inside it, the ambient Cx spawns owned children.
        let lane = LatencyLane::new(1).expect("lane");
        let root = lane.runtime().handle().spawn(async {
            let cx = asupersync::Cx::current().expect("task Cx");
            let mut task = cx.spawn(|_child| async move { 7u64 }).expect("spawn child");
            task.join(&cx).await.expect("child joins")
        });
        let joined = lane.block_on(root);
        assert_eq!(joined, 7);
    }
}
