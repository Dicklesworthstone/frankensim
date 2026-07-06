//! The tile-kernel contract (plan §5.2, Appendix B): every hot kernel is a
//! tile program. Tiles are the unit of scheduling, cancellation,
//! determinism, and NUMA placement.

use crate::cx::{Cancelled, Cx};
use core::ops::ControlFlow;

/// A value that can be merged up a FIXED-SHAPE reduction tree. The executor
/// stores each tile's output in a slot keyed by tile index and folds slots
/// in tile order — so the reduction shape (and therefore every rounding
/// decision) is a function of the PLAN, never of scheduling. `merge` need
/// not be commutative; it is always applied in ascending tile order.
pub trait Reduce: Send {
    /// The fold identity (`identity().merge(x) == x`).
    fn identity() -> Self;
    /// Fold `other` (the NEXT tile in order) into `self`.
    #[must_use]
    fn merge(self, other: Self) -> Self;
}

impl Reduce for () {
    fn identity() {}

    fn merge(self, (): Self) {}
}

impl Reduce for u64 {
    fn identity() -> Self {
        0
    }

    fn merge(self, other: Self) -> Self {
        self.wrapping_add(other)
    }
}

impl Reduce for f64 {
    fn identity() -> Self {
        0.0
    }

    fn merge(self, other: Self) -> Self {
        self + other
    }
}

impl<T: Send> Reduce for Vec<T> {
    fn identity() -> Self {
        Vec::new()
    }

    fn merge(mut self, mut other: Self) -> Self {
        self.append(&mut other);
        self
    }
}

/// The tile geometry of one kernel invocation. The autotuner will eventually
/// choose tile counts/shapes; kernels state them explicitly until then.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TilePlan {
    /// Number of tiles (tile indices are `0..tiles`).
    pub tiles: u64,
    /// Stable kernel name; hashed into every tile's [`crate::StreamKey`]
    /// and stamped on events.
    pub kernel: &'static str,
}

impl TilePlan {
    /// Plan `tiles` tiles for the named kernel.
    #[must_use]
    pub const fn new(kernel: &'static str, tiles: u64) -> Self {
        TilePlan { tiles, kernel }
    }

    /// Stable kernel id (FNV-1a of the kernel name).
    #[must_use]
    pub fn kernel_id(&self) -> u64 {
        fs_obs::fnv1a64(self.kernel.as_bytes())
    }
}

/// Every hot kernel is a `TileKernel` (plan Appendix B). `run` must poll
/// `cx.checkpoint()` at least once per tile (the executor also polls at
/// every tile boundary) and return `Break(Cancelled)` promptly when it
/// observes a request.
pub trait TileKernel: Sync {
    /// Per-tile output, merged up the fixed-shape tree.
    type Out: Reduce;

    /// The tile geometry of this invocation.
    fn tiles(&self) -> TilePlan;

    /// Execute one tile under its context.
    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduce_identities_hold() {
        assert_eq!(u64::identity().merge(7), 7);
        assert_eq!(7u64.merge(u64::identity()), 7);
        assert!((f64::identity().merge(1.5) - 1.5).abs() < f64::EPSILON);
        let v = Vec::<u32>::identity().merge(vec![1, 2]);
        assert_eq!(v, vec![1, 2]);
        let _: () = <() as Reduce>::identity();
    }

    #[test]
    fn plan_kernel_ids_are_stable_and_distinct() {
        let a = TilePlan::new("kernel-a", 4);
        let b = TilePlan::new("kernel-b", 4);
        assert_eq!(a.kernel_id(), TilePlan::new("kernel-a", 9).kernel_id());
        assert_ne!(a.kernel_id(), b.kernel_id());
    }
}
