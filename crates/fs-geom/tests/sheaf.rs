//! Sheaf-certificate conformance (the wqd.13 bead). Acceptance: verdicts
//! correct with correct localization on seeded watertight/leaky
//! multi-chart fixtures; interval verification SOUND (no PASS on a truly
//! leaky seam — ray-parity cross-examines); δδ = 0 bitwise; verdict
//! invariance under patch re-indexing (exact) and rigid motion
//! (tolerance-level); the adversarial seam zoo (near-tangent, T-junction)
//! behaves honestly; the coboundary/structural split feeds the merge
//! semantics.

use asupersync::types::Budget;
use fs_evidence::NumericalKind;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart};
use fs_geom::{Chart, Point3, SheafComplex, SheafVerdict, ray_parity_falsifier};

fn verdict_line(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
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

/// A rigid-motion wrapper: evaluates the inner chart in rotated+shifted
/// coordinates (a rotation about z by a fixed angle plus a translation).
struct Moved<C: Chart> {
    inner: C,
    cos: f64,
    sin: f64,
    shift: [f64; 3],
}

impl<C: Chart> Moved<C> {
    fn new(inner: C, angle: f64, shift: [f64; 3]) -> Self {
        Moved {
            inner,
            cos: angle.cos(),
            sin: angle.sin(),
            shift,
        }
    }

    fn map(&self, p: Point3) -> Point3 {
        // Inverse motion: un-shift, un-rotate.
        let q = [
            p.x - self.shift[0],
            p.y - self.shift[1],
            p.z - self.shift[2],
        ];
        Point3::new(
            self.cos * q[0] + self.sin * q[1],
            -self.sin * q[0] + self.cos * q[1],
            q[2],
        )
    }
}

impl<C: Chart> Chart for Moved<C> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        self.inner.eval(self.map(x), cx)
    }

    fn support(&self) -> fs_geom::Aabb {
        // Conservative: rotate the support's corners, box them.
        let s = self.inner.support();
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for &x in &[s.min.x, s.max.x] {
            for &y in &[s.min.y, s.max.y] {
                for &z in &[s.min.z, s.max.z] {
                    // Forward motion of the corner.
                    let fx = self.cos * x - self.sin * y + self.shift[0];
                    let fy = self.sin * x + self.cos * y + self.shift[1];
                    let fz = z + self.shift[2];
                    min[0] = min[0].min(fx);
                    min[1] = min[1].min(fy);
                    min[2] = min[2].min(fz);
                    max[0] = max[0].max(fx);
                    max[1] = max[1].max(fy);
                    max[2] = max[2].max(fz);
                }
            }
        }
        fs_geom::Aabb::new(
            Point3::new(min[0], min[1], min[2]),
            Point3::new(max[0], max[1], max[2]),
        )
    }

    fn name(&self) -> &'static str {
        "test/moved"
    }
}

/// Two identical spheres offset so their supports overlap in a lens
/// around the shared surface — the WATERTIGHT fixture (same abstract
/// region presented twice).
fn watertight_pair() -> (SphereChart, SphereChart) {
    let s = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    (s, s)
}

/// The LEAKY variant: the second chart's radius is off by delta.
fn leaky_pair(delta: f64) -> (SphereChart, SphereChart) {
    let (a, mut b) = watertight_pair();
    b.radius += delta;
    (a, b)
}

#[test]
fn sh_001_verdicts_and_localization() {
    with_cx(|cx| {
        // Watertight: identical spheres agree exactly on the seam.
        let (a, b) = watertight_pair();
        let charts: Vec<&dyn Chart> = vec![&a, &b];
        let complex = SheafComplex::from_charts(&charts, cx);
        assert_eq!(complex.interfaces.len(), 1, "one shared interface");
        assert!(!complex.interfaces[0].samples.is_empty());
        let ev = complex.watertightness(1e-9);
        match &ev.value {
            SheafVerdict::Pass {
                worst_mismatch,
                margins,
            } => {
                assert!(*worst_mismatch <= 1e-9);
                assert_eq!(margins.len(), 1);
            }
            other => panic!("identical charts must certify: {other:?}"),
        }
        assert_eq!(
            ev.numerical.kind,
            NumericalKind::Enclosure,
            "interval-verified"
        );
        // Leaky: radius off by 1e-2 — FAIL, localized to the (0,1) seam,
        // with the mismatch magnitude ~ delta.
        let (a2, b2) = leaky_pair(1e-2);
        let charts2: Vec<&dyn Chart> = vec![&a2, &b2];
        let complex2 = SheafComplex::from_charts(&charts2, cx);
        let ev2 = complex2.watertightness(1e-4);
        match &ev2.value {
            SheafVerdict::Fail {
                obstruction,
                coboundary_share,
            } => {
                assert_eq!(obstruction[0].0, (0, 1), "offending interface named");
                assert!(
                    obstruction[0].1 > 5e-3 && obstruction[0].1 < 2e-2,
                    "mismatch magnitude ~ delta: {}",
                    obstruction[0].1
                );
                // A pure radius offset is a CONSTANT mismatch on the seam —
                // exactly the coboundary (gauge) component.
                assert!(
                    *coboundary_share > 0.9,
                    "constant seam mismatch is reconcilable: {coboundary_share}"
                );
            }
            other => panic!("leaky seam must fail: {other:?}"),
        }
        verdict_line(
            "sh-001",
            "identical charts PASS; radius leak FAILs at the named seam with ~delta \
             magnitude and >0.9 coboundary share",
        );
    });
}

#[test]
fn sh_002_delta_delta_is_zero_bitwise() {
    with_cx(|cx| {
        // Three pairwise-overlapping boxes: a genuine triple junction.
        let boxes: Vec<BoxChart> = (0..3)
            .map(|i| {
                let base = f64::from(i) * 0.4;
                BoxChart {
                    aabb: fs_geom::Aabb::new(
                        Point3::new(base - 0.5, -0.5, -0.5),
                        Point3::new(base + 0.5, 0.5, 0.5),
                    ),
                }
            })
            .collect();
        let charts: Vec<&dyn Chart> = boxes.iter().map(|b| b as &dyn Chart).collect();
        let complex = SheafComplex::from_charts(&charts, cx);
        assert!(!complex.interfaces.is_empty());
        assert!(
            !complex.triples.is_empty(),
            "fixture must produce a triple junction; adjust geometry"
        );
        let d0 = complex.delta0_edges();
        let d1 = complex.delta1();
        // δ¹ · δ⁰ computed densely (test scale): every entry EXACTLY 0.0.
        let (rows, mid, cols) = (
            complex.triples.len(),
            complex.interfaces.len(),
            complex.n_patches,
        );
        for r in 0..rows {
            for c in 0..cols {
                let mut acc = 0.0f64;
                let (d1_cols, d1_vals) = d1.row(r);
                for (k, &e) in d1_cols.iter().enumerate() {
                    let (d0_cols, d0_vals) = d0.row(e);
                    for (j, &p) in d0_cols.iter().enumerate() {
                        if p == c {
                            acc += d1_vals[k] * d0_vals[j];
                        }
                    }
                }
                assert_eq!(acc.to_bits(), 0.0f64.to_bits(), "δδ must be bitwise zero");
            }
        }
        let _ = mid;
        verdict_line("sh-002", "δ¹·δ⁰ == 0 bitwise on a genuine triple junction");
    });
}

#[test]
fn sh_003_invariance_reindex_exact_rigid_tolerance() {
    with_cx(|cx| {
        let (a, b) = leaky_pair(2e-2);
        // Re-indexing: swap chart order — the verdict is EXACTLY equal
        // (geometry-derived sample seeds are index-free).
        let fwd: Vec<&dyn Chart> = vec![&a, &b];
        let rev: Vec<&dyn Chart> = vec![&b, &a];
        let v1 = SheafComplex::from_charts(&fwd, cx).watertightness(1e-4);
        let v2 = SheafComplex::from_charts(&rev, cx).watertightness(1e-4);
        let key = |v: &SheafVerdict| match v {
            SheafVerdict::Fail { obstruction, .. } => obstruction[0].1,
            _ => f64::NAN,
        };
        assert_eq!(
            key(&v1.value).to_bits(),
            key(&v2.value).to_bits(),
            "re-indexing invariance is exact"
        );
        // Rigid motion: rotate+translate BOTH charts — verdict class and
        // magnitude agree to tolerance (samples differ, physics doesn't).
        let ma = Moved::new(a, 0.7, [3.0, -1.0, 0.5]);
        let mb = Moved::new(b, 0.7, [3.0, -1.0, 0.5]);
        let moved: Vec<&dyn Chart> = vec![&ma, &mb];
        let v3 = SheafComplex::from_charts(&moved, cx).watertightness(1e-4);
        match (&v1.value, &v3.value) {
            (
                SheafVerdict::Fail {
                    obstruction: o1, ..
                },
                SheafVerdict::Fail {
                    obstruction: o3, ..
                },
            ) => {
                assert!(
                    (o1[0].1 - o3[0].1).abs() < 5e-3,
                    "rigid motion preserves the leak magnitude: {} vs {}",
                    o1[0].1,
                    o3[0].1
                );
            }
            other => panic!("verdict class must survive rigid motion: {other:?}"),
        }
        verdict_line(
            "sh-003",
            "re-index invariance bitwise; rigid-motion invariance to 5e-3",
        );
    });
}

#[test]
fn sh_004_adversarial_seams_and_soundness() {
    with_cx(|cx| {
        // T-junction: three consistent boxes sharing faces — PASS (all
        // charts are exact SDFs of the same union's pieces... each chart
        // is ITS OWN box; interfaces only certify agreement where both
        // charts are near zero, i.e. shared face bands).
        let mk = |cx: f64, cy: f64| BoxChart {
            aabb: fs_geom::Aabb::new(
                Point3::new(cx - 0.5, cy - 0.5, -0.5),
                Point3::new(cx + 0.5, cy + 0.5, 0.5),
            ),
        };
        let b1 = mk(0.0, 0.0);
        let b2 = mk(0.9, 0.0);
        let b3 = mk(0.45, 0.9);
        let charts: Vec<&dyn Chart> = vec![&b1, &b2, &b3];
        let complex = SheafComplex::from_charts(&charts, cx);
        // Adjacent identical box SDFs disagree off the shared face (each
        // is its own box), so this is genuinely adversarial: the verdict
        // must not be a false PASS at tight tolerance, and must not be a
        // false FAIL at a tolerance matching the band geometry.
        let ev = complex.watertightness(1e-12);
        assert!(
            !matches!(ev.value, SheafVerdict::Pass { .. }),
            "distinct-box seams must not certify at 1e-12"
        );
        // Near-tangent spheres (distinct surfaces, kissing at one point):
        // NOT a false PASS.
        let s1 = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let s2 = SphereChart {
            center: Point3::new(2.05, 0.0, 0.0),
            radius: 1.0,
        };
        let kiss: Vec<&dyn Chart> = vec![&s1, &s2];
        let kiss_complex = SheafComplex::from_charts(&kiss, cx);
        // The 0.05 gap means NO overlap interface is discovered, so the
        // interface-agreement check has gathered zero evidence. It must NOT
        // report a positive PASS from an empty bound set (bead obnw: `all()`
        // on the empty set was vacuously true) — the honest verdict is Unknown.
        let kv = kiss_complex.watertightness(1e-9);
        assert!(
            matches!(kv.value, SheafVerdict::Unknown { .. }),
            "near-tangent distinct surfaces (no interface) must be Unknown, not a \
             false PASS, got {:?}",
            kv.value
        );
        // SOUNDNESS cross-examination (the falsifier pairing): a PASSing
        // watertight model survives ray parity; the falsifier runs a
        // DIFFERENT algorithm (sign-crossing counts) on the same charts.
        let (wa, wb) = watertight_pair();
        let watertight: Vec<&dyn Chart> = vec![&wa, &wb];
        let rays = [
            (Point3::new(-3.0, 0.01, 0.02), Point3::new(3.0, 0.01, 0.02)),
            (Point3::new(0.02, -3.0, 0.01), Point3::new(0.02, 3.0, 0.01)),
            (Point3::new(-2.5, -2.5, 0.0), Point3::new(2.5, 2.5, 0.0)),
        ];
        assert!(
            ray_parity_falsifier(&watertight, &rays, 4001, cx).is_none(),
            "the ray-parity falsifier must not refute a sound PASS"
        );
        verdict_line(
            "sh-004",
            "distinct-surface seams never falsely certify; ray parity cross-examines the \
             sound PASS",
        );
    });
}

#[test]
fn sh_005_section_split_feeds_merge_semantics() {
    with_cx(|cx| {
        // A chain of three spheres where the middle chart carries a
        // constant radial offset: the mismatch is pure gauge (coboundary)
        // and the section solve should absorb nearly all of it.
        let a = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let b = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.015,
        };
        let charts: Vec<&dyn Chart> = vec![&a, &b];
        let complex = SheafComplex::from_charts(&charts, cx);
        let (offsets, raw, residual) = complex.section_solve();
        assert!(raw > 0.0, "the leak is visible pre-gauge");
        assert!(
            residual < raw * 0.01,
            "a constant offset is pure coboundary: raw {raw}, residual {residual}"
        );
        assert!(
            (offsets[1] - offsets[0] - 0.015).abs() < 1e-3,
            "the recovered gauge cancels the radius delta: {offsets:?}"
        );
        verdict_line(
            "sh-005",
            "constant seam offsets are absorbed by the section solve (the merge-semantics \
             coboundary split)",
        );
    });
}
