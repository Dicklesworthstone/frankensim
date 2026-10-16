//! Claim-integrity regressions for the fixture topology hints and the
//! sampled-SDF reconstruction radius (`docs/CLAIM_INTEGRITY.md`).
//!
//! - `frankensim-extreal-program-f85xj.2.18` — all three fixture charts
//!   published `BettiBounds::exact` with `b2 = 1`, i.e. "exactly one
//!   enclosed void", for SOLID regions that enclose none (and the torus
//!   additionally claimed the wrong `b1`). `BettiBounds` describes the
//!   region where the signed distance is NEGATIVE, so the honest values
//!   are ball `(1, 0, 0)`, box `(1, 0, 0)`, solid torus `(1, 1, 0)` —
//!   the convention `fs_topo::cubical::verify_topology` measures and its
//!   conformance suite pins (`tests/conformance.rs`, topo-003).
//! - `.2.19` — `SampledSdf::bound()` was documented as a reconstruction
//!   bound "useful even when the source makes no abstract-distance
//!   claim", but for any non-`ExactDistance` source it is built from a
//!   maximum of per-node LOCAL Lipschitz values, which bound nothing
//!   across a cell. It now carries its authority.

use asupersync::types::Budget;
use fs_evidence::NumericalKind;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart, TorusChart};
use fs_geom::{Aabb, BettiBounds, Chart, Convert, ErrBudget, Point3, SampledSdf, TraceStepClaim};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x2_18,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// ---------------------------------------------------------------- .2.18

/// Independent evidence for the `b2 = 0` half of the claim, and for
/// `b0 = 1`: voxelize the chart on a coarse grid over an inflated box and
/// check (a) every EMPTY cell is 6-connected to the grid boundary, so no
/// void is enclosed, and (b) every FILLED cell is 6-connected to every
/// other, so there is one component.
///
/// This deliberately does not import `fs-topo` (which depends on
/// `fs-geom`); it is a small independent flood fill, not the oracle.
fn flood_check(chart: &dyn Chart, cx: &Cx<'_>, n: usize) -> (usize, usize) {
    let support = chart.support();
    let pad = 0.35
        * (support.max.x - support.min.x)
            .max(support.max.y - support.min.y)
            .max(support.max.z - support.min.z);
    let lo = Point3::new(
        support.min.x - pad,
        support.min.y - pad,
        support.min.z - pad,
    );
    let hi = Point3::new(
        support.max.x + pad,
        support.max.y + pad,
        support.max.z + pad,
    );
    let coord = |lo: f64, hi: f64, i: usize| lo + (hi - lo) * (i as f64 + 0.5) / (n as f64);
    let mut filled = vec![false; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let p = Point3::new(
                    coord(lo.x, hi.x, i),
                    coord(lo.y, hi.y, j),
                    coord(lo.z, hi.z, k),
                );
                filled[(k * n + j) * n + i] = chart.eval(p, cx).signed_distance < 0.0;
            }
        }
    }
    let flood = |want: bool, seeds: Vec<usize>| -> usize {
        let mut seen = vec![false; n * n * n];
        let mut stack = Vec::new();
        for s in seeds {
            if filled[s] == want && !seen[s] {
                seen[s] = true;
                stack.push(s);
            }
        }
        let mut count = 0;
        while let Some(c) = stack.pop() {
            count += 1;
            let (i, j, k) = (c % n, (c / n) % n, c / (n * n));
            let mut push = |i: usize, j: usize, k: usize, stack: &mut Vec<usize>| {
                let idx = (k * n + j) * n + i;
                if filled[idx] == want && !seen[idx] {
                    seen[idx] = true;
                    stack.push(idx);
                }
            };
            if i > 0 {
                push(i - 1, j, k, &mut stack);
            }
            if i + 1 < n {
                push(i + 1, j, k, &mut stack);
            }
            if j > 0 {
                push(i, j - 1, k, &mut stack);
            }
            if j + 1 < n {
                push(i, j + 1, k, &mut stack);
            }
            if k > 0 {
                push(i, j, k - 1, &mut stack);
            }
            if k + 1 < n {
                push(i, j, k + 1, &mut stack);
            }
        }
        count
    };
    // Empty phase, seeded from the whole grid boundary.
    let mut boundary = Vec::new();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                if i == 0 || j == 0 || k == 0 || i + 1 == n || j + 1 == n || k + 1 == n {
                    boundary.push((k * n + j) * n + i);
                }
            }
        }
    }
    let empty_total = filled.iter().filter(|f| !**f).count();
    let empty_reached = flood(false, boundary);
    // Filled phase, seeded from the first filled cell.
    let first_filled = filled
        .iter()
        .position(|f| *f)
        .expect("the fixture has a non-empty solid region");
    let filled_total = filled.iter().filter(|f| **f).count();
    let filled_reached = flood(true, vec![first_filled]);
    assert_eq!(
        empty_total, empty_reached,
        "an unreachable empty cell would be an ENCLOSED VOID (b2 > 0)"
    );
    assert_eq!(
        filled_total, filled_reached,
        "the solid region must be one connected component (b0 = 1)"
    );
    (filled_total, empty_total)
}

#[test]
fn solid_fixture_charts_publish_solid_region_betti_numbers() {
    let sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let boxed = BoxChart {
        aabb: Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 0.5, 2.0)),
    };
    let torus = TorusChart {
        center: Point3::new(0.0, 0.0, 0.0),
        major: 1.0,
        minor: 0.35,
    };

    // A solid ball, a solid box and a SOLID torus enclose no void.
    assert_eq!(sphere.topology_hint(), BettiBounds::exact(1, 0, 0));
    assert_eq!(boxed.topology_hint(), BettiBounds::exact(1, 0, 0));
    assert_eq!(torus.topology_hint(), BettiBounds::exact(1, 1, 0));

    with_cx(|cx| {
        flood_check(&sphere, cx, 24);
        flood_check(&boxed, cx, 24);
        flood_check(&torus, cx, 32);
    });
}

#[test]
fn degenerate_fixture_parameters_publish_unknown_topology() {
    let flat = BoxChart {
        aabb: Aabb::new(Point3::new(-1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
    };
    assert_eq!(flat.topology_hint(), BettiBounds::unknown());

    // Horn/spindle parameters are not a solid torus; no claim is made.
    let spindle = TorusChart {
        center: Point3::new(0.0, 0.0, 0.0),
        major: 0.5,
        minor: 1.0,
    };
    assert_eq!(spindle.topology_hint(), BettiBounds::unknown());

    let empty = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.0,
    };
    assert_eq!(
        empty.topology_hint(),
        BettiBounds::unknown(),
        "a zero-radius sphere has an empty negative region, not one component"
    );
}

// ---------------------------------------------------------------- .2.19

#[test]
fn the_reconstruction_radius_carries_its_authority() {
    with_cx(|cx| {
        // ExactDistance source: the unit-Lipschitz theorem is GLOBAL, so
        // `L · cell diagonal` really does bound the trilinear error.
        let exact_source = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        assert_eq!(
            exact_source.trace_step_claim(),
            TraceStepClaim::ExactDistance
        );
        let converted: fs_evidence::Certified<SampledSdf> = exact_source
            .convert(ErrBudget { abs_sd_error: 0.5 }, cx)
            .expect("a unit sphere converts inside a 0.5 budget");
        assert_eq!(
            converted.value.nominal_field_bound_kind(),
            NumericalKind::Enclosure,
            "a global Lipschitz theorem makes the reconstruction radius rigorous"
        );

        // Weaker source: `l_max` is a maximum of per-node LOCAL Lipschitz
        // values with no stated radius, so a sub-cell slope spike can
        // exceed the published radius. It must NOT read as a bound.
        let weak_source = TorusChart {
            center: Point3::new(0.0, 0.0, 0.0),
            major: 0.5,
            minor: 1.0,
        };
        assert_ne!(
            weak_source.trace_step_claim(),
            TraceStepClaim::ExactDistance,
            "the spindle torus is the weak-source fixture"
        );
        let weak: fs_evidence::Evidence<SampledSdf> = weak_source
            .convert_clipped(
                ErrBudget { abs_sd_error: 0.5 },
                Aabb::new(Point3::new(-2.0, -2.0, -1.5), Point3::new(2.0, 2.0, 1.5)),
                cx,
            )
            .expect("the spindle torus converts inside a 0.5 budget");
        assert_eq!(
            weak.value.nominal_field_bound_kind(),
            NumericalKind::Estimate,
            "sample-local Lipschitz maxima parameterize an ESTIMATE, not a bound"
        );
        assert_eq!(
            weak.value.bound().to_bits(),
            weak.value.nominal_field_bound().to_bits(),
            "`bound()` is an alias and inherits the same authority"
        );
    });
}
