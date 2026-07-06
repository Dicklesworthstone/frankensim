//! Built-in registry kernels: the fs-simd primitive set (the only hot
//! kernels that exist yet) plus a deliberately de-optimized kernel for
//! harness meta-tests.
//!
//! Intensity models (the hand-calculation basis for `rf_002`):
//! - axpy `y = a·x + y`: reads x and y, writes y → 24 B/elem, 2 flop/elem.
//! - dot `Σ x·y`: reads x and y → 16 B/elem, 2 flop/elem.
//! - sum `Σ x`: reads x → 8 B/elem, 1 flop/elem.
//!
//! Targets here are report-only in v0 (`target_fraction: None` except where
//! a band is deliberately claimed for meta-testing): CI runners are shared
//! machines and §14 bands belong to fingerprinted reference machines.

use crate::{KernelSpec, RooflineKernel, Threading};

/// `y = a·x + y` through the dispatched fs-simd table.
pub struct AxpyKernel {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl AxpyKernel {
    /// Buffers of `n` elements each (pick `n` large enough to stream past
    /// the last-level cache when measuring the bandwidth roof).
    #[must_use]
    pub fn new(n: usize) -> AxpyKernel {
        AxpyKernel {
            x: vec![1.5; n],
            y: vec![0.5; n],
        }
    }
}

impl RooflineKernel for AxpyKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-axpy-f64",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        (fs_simd::ops().axpy)(1.000_000_1, &self.x, &mut self.y);
        std::hint::black_box(self.y[self.y.len() / 2]);
    }
}

/// `Σ x·y` through the dispatched fs-simd table.
pub struct DotKernel {
    x: Vec<f64>,
    y: Vec<f64>,
    out: f64,
}

impl DotKernel {
    /// Buffers of `n` elements each.
    #[must_use]
    pub fn new(n: usize) -> DotKernel {
        DotKernel {
            x: vec![1.5; n],
            y: vec![0.5; n],
            out: 0.0,
        }
    }
}

impl RooflineKernel for DotKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-dot-f64",
            version: "1",
            bytes_per_elem: 16.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        self.out = (fs_simd::ops().dot)(&self.x, &self.y);
        std::hint::black_box(self.out);
    }
}

/// `Σ x` through the dispatched fs-simd table.
pub struct SumKernel {
    x: Vec<f64>,
    out: f64,
}

impl SumKernel {
    /// A buffer of `n` elements.
    #[must_use]
    pub fn new(n: usize) -> SumKernel {
        SumKernel {
            x: vec![0.25; n],
            out: 0.0,
        }
    }
}

impl RooflineKernel for SumKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-sum-f64",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        self.out = (fs_simd::ops().sum)(&self.x);
        std::hint::black_box(self.out);
    }
}

/// Deliberately de-optimized kernel with a band it cannot meet: a serial
/// dependency chain strided across a buffer, claiming 90% of the bandwidth
/// roof. The harness meta-test asserts it reports `BelowBand` — proof that
/// a slow kernel is caught, not absorbed (bead acceptance criterion).
pub struct SeededSlowKernel {
    x: Vec<f64>,
    out: f64,
}

impl SeededSlowKernel {
    /// A buffer of `n` elements.
    #[must_use]
    pub fn new(n: usize) -> SeededSlowKernel {
        SeededSlowKernel {
            x: (0..n).map(|i| (i % 7) as f64).collect(),
            out: 0.0,
        }
    }
}

impl RooflineKernel for SeededSlowKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "seeded-slow",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.9),
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        // Serial chain + division: nowhere near any roof, by construction.
        let mut acc = 1.0f64;
        for &v in &self.x {
            acc = (acc / 1.000_000_01) + v.sqrt().sin();
        }
        self.out = acc;
        std::hint::black_box(self.out);
    }
}

/// The default registry: everything that exists today.
#[must_use]
pub fn default_registry(n: usize) -> Vec<Box<dyn RooflineKernel>> {
    vec![
        Box::new(AxpyKernel::new(n)),
        Box::new(DotKernel::new(n)),
        Box::new(SumKernel::new(n)),
    ]
}
