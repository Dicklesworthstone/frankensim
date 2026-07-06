//! Dense tiled SDF grids (plan §7.2): f32 STORAGE / f64 EVALUATION on
//! fs-substrate's Morton/tile-major fields, with C¹ triquadratic B-spline
//! reconstruction — continuous gradients matter because shape optimization
//! DIFFERENTIATES through samples (plan §7.6).
//!
//! Error honesty: the chart carries a CONSTRUCTED bound (sampling `L·h` +
//! B-spline smoothing `L·h` + f32 quantization) as its declared enclosure,
//! and MEASURED eikonal statistics ("how much is this NOT a distance
//! field") as separate, clearly-labeled evidence for sphere-tracing
//! safety margins.

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3, Vec3};

/// A dense signed-distance grid over a box, built from any source chart
/// with a certified Lipschitz bound.
#[derive(Debug)]
pub struct TiledSdf {
    field: fs_substrate::field::TiledField<f32>,
    box_: Aabb,
    /// Grid step per axis.
    h: [f64; 3],
    /// Samples per axis.
    n: [u32; 3],
    /// Rigorous |chart - source| bound inside the box (constructed).
    bound: f64,
    /// The source's certified Lipschitz constant.
    source_lipschitz: f64,
    /// Measured max |,∇φ| − 1, over the construction probe set (an
    /// ESTIMATE, clearly labeled — see [`TiledSdf::eikonal_stats`]).
    eikonal_dev: f64,
}

/// Construction failure (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum SdfBuildError {
    /// The source chart certifies no finite Lipschitz bound.
    NoLipschitzBound,
    /// The requested step would need more samples per axis than the cap.
    ResolutionTooFine {
        /// Samples/axis the request needs.
        need: u64,
        /// The cap.
        cap: u64,
    },
}

impl core::fmt::Display for SdfBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SdfBuildError::NoLipschitzBound => write!(
                f,
                "dense SDF build refused: the source chart certifies no Lipschitz bound, so \
                 no rigorous sampling error exists; use a certified source"
            ),
            SdfBuildError::ResolutionTooFine { need, cap } => write!(
                f,
                "dense SDF build refused: {need} samples/axis exceed the {cap} cap; coarsen \
                 the step, shrink the box, or use the sparse VDB/adaptive charts"
            ),
        }
    }
}

impl core::error::Error for SdfBuildError {}

/// Per-axis sample cap (beyond this, dense is the wrong tool — the
/// refusal says which chart to use instead).
pub const DENSE_MAX_SAMPLES_PER_AXIS: u64 = 512;

/// Measured eikonal statistics (`|∇φ| − 1` over a seeded probe set):
/// evidence for sphere-tracing safety, LABELED as measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EikonalStats {
    /// Mean absolute deviation.
    pub mean_abs_dev: f64,
    /// Maximum absolute deviation observed.
    pub max_abs_dev: f64,
    /// Probe count.
    pub probes: u64,
}

impl TiledSdf {
    /// Sample `source` over its inflated support at step `target_h`.
    /// The declared error bound is `2·L·h_max + q` (sampling + quadratic
    /// B-spline smoothing of an L-Lipschitz field, plus f32 quantization),
    /// conservative by construction.
    ///
    /// # Errors
    /// [`SdfBuildError`] (teaching refusals; nothing runs before checks).
    pub fn build(
        source: &dyn Chart,
        target_h: f64,
        cx: &Cx<'_>,
    ) -> Result<TiledSdf, SdfBuildError> {
        let box_ = source.support().inflate(3.0 * target_h.max(1e-9));
        let probe = source.eval(
            Point3::new(
                f64::midpoint(box_.min.x, box_.max.x),
                f64::midpoint(box_.min.y, box_.max.y),
                f64::midpoint(box_.min.z, box_.max.z),
            ),
            cx,
        );
        let lipschitz = match probe.lipschitz {
            Some(l) if l.is_finite() => l,
            _ => return Err(SdfBuildError::NoLipschitzBound),
        };
        let edges = [
            box_.max.x - box_.min.x,
            box_.max.y - box_.min.y,
            box_.max.z - box_.min.z,
        ];
        let mut n = [0u32; 3];
        let mut h = [0.0f64; 3];
        for a in 0..3 {
            let need = (edges[a] / target_h).ceil() as u64 + 1;
            if need > DENSE_MAX_SAMPLES_PER_AXIS {
                return Err(SdfBuildError::ResolutionTooFine {
                    need,
                    cap: DENSE_MAX_SAMPLES_PER_AXIS,
                });
            }
            n[a] = (need as u32).max(4);
            h[a] = edges[a] / f64::from(n[a] - 1);
        }
        let grid = fs_substrate::tile::TileGrid::new(n, fs_substrate::tile::TileEdge::E8)
            .expect("caps keep the grid within Morton bounds");
        let mut field = fs_substrate::field::TiledField::new(grid, 0.0f32);
        let mut max_abs = 0.0f64;
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let p = Point3::new(
                        box_.min.x + f64::from(i) * h[0],
                        box_.min.y + f64::from(j) * h[1],
                        box_.min.z + f64::from(k) * h[2],
                    );
                    let sd = source.eval(p, cx).signed_distance;
                    max_abs = max_abs.max(sd.abs());
                    field.set([i, j, k], sd as f32);
                }
            }
        }
        let h_max = h[0].max(h[1]).max(h[2]);
        // f32 quantization of values up to max_abs.
        let quant = max_abs * f64::from(f32::EPSILON);
        let bound = 2.0 * lipschitz * h_max + quant;
        let mut sdf = TiledSdf {
            field,
            box_,
            h,
            n,
            bound,
            source_lipschitz: lipschitz,
            eikonal_dev: 0.0,
        };
        sdf.eikonal_dev = sdf.measure_eikonal(0x51DF_0001, 2_000, cx).max_abs_dev;
        Ok(sdf)
    }

    /// The declared in-box error bound.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Grid steps per axis.
    #[must_use]
    pub fn steps(&self) -> [f64; 3] {
        self.h
    }

    /// Re-sample every grid sample whose position lies inside `region`
    /// (inflated by one spline-support cell so reconstruction around the
    /// edit seam sees fresh data) — the incremental path for locally
    /// edited sources. Samples are recomputed at EXACTLY the original
    /// positions, so an incremental update is bit-identical to a full
    /// rebuild of the same source (the converter beads' G5 law). Returns
    /// the number of samples refreshed.
    ///
    /// # Errors
    /// [`fs_exec::Cancelled`] mid-update (the field stays valid: every
    /// written sample is complete).
    pub fn resample_box(
        &mut self,
        source: &dyn Chart,
        region: fs_geom::Aabb,
        cx: &Cx<'_>,
    ) -> Result<u64, fs_exec::Cancelled> {
        let pad = 2.0 * self.h[0].max(self.h[1]).max(self.h[2]);
        let r = region.inflate(pad);
        let lo = |axis: usize, min: f64| -> u32 {
            let base = match axis {
                0 => self.box_.min.x,
                1 => self.box_.min.y,
                _ => self.box_.min.z,
            };
            (((min - base) / self.h[axis]).floor().max(0.0)) as u32
        };
        let hi = |axis: usize, max: f64, n: u32| -> u32 {
            let base = match axis {
                0 => self.box_.min.x,
                1 => self.box_.min.y,
                _ => self.box_.min.z,
            };
            (((max - base) / self.h[axis]).ceil() as u32).min(n - 1)
        };
        let (i0, i1) = (lo(0, r.min.x), hi(0, r.max.x, self.n[0]));
        let (j0, j1) = (lo(1, r.min.y), hi(1, r.max.y, self.n[1]));
        let (k0, k1) = (lo(2, r.min.z), hi(2, r.max.z, self.n[2]));
        let mut updated = 0u64;
        for k in k0..=k1 {
            for j in j0..=j1 {
                cx.checkpoint()?;
                for i in i0..=i1 {
                    let p = Point3::new(
                        self.box_.min.x + f64::from(i) * self.h[0],
                        self.box_.min.y + f64::from(j) * self.h[1],
                        self.box_.min.z + f64::from(k) * self.h[2],
                    );
                    let sd = source.eval(p, cx).signed_distance;
                    self.field.set([i, j, k], sd as f32);
                    updated += 1;
                }
            }
        }
        Ok(updated)
    }

    /// Quadratic B-spline basis at fractional offset `t ∈ [0,1)` for the
    /// three support samples (standard uniform quadratic B-spline).
    fn bspline_w(t: f64) -> [f64; 3] {
        [
            0.5 * (1.0 - t) * (1.0 - t),
            0.5 + t * (1.0 - t),
            0.5 * t * t,
        ]
    }

    /// Derivative of [`Self::bspline_w`] with respect to t.
    fn bspline_dw(t: f64) -> [f64; 3] {
        [t - 1.0, 1.0 - 2.0 * t, t]
    }

    fn sample_raw(&self, i: i64, j: i64, k: i64) -> f64 {
        let c = [
            i.clamp(0, i64::from(self.n[0]) - 1) as u32,
            j.clamp(0, i64::from(self.n[1]) - 1) as u32,
            k.clamp(0, i64::from(self.n[2]) - 1) as u32,
        ];
        f64::from(self.field.get(c))
    }

    /// Evaluate value and gradient of the C¹ triquadratic reconstruction.
    fn spline_eval(&self, x: Point3) -> (f64, Vec3) {
        // Cell-centered quadratic B-spline: base sample nearest the point,
        // fractional offset in [-0.5, 0.5) mapped to t = frac + 0.5.
        let mut base = [0i64; 3];
        let mut t = [0.0f64; 3];
        let coords = [
            (x.x - self.box_.min.x) / self.h[0],
            (x.y - self.box_.min.y) / self.h[1],
            (x.z - self.box_.min.z) / self.h[2],
        ];
        for a in 0..3 {
            let c = coords[a].clamp(0.0, f64::from(self.n[a] - 1));
            let b = c.round();
            base[a] = b as i64;
            t[a] = (c - b) + 0.5;
        }
        let (wx, wy, wz) = (
            Self::bspline_w(t[0]),
            Self::bspline_w(t[1]),
            Self::bspline_w(t[2]),
        );
        let (dwx, dwy, dwz) = (
            Self::bspline_dw(t[0]),
            Self::bspline_dw(t[1]),
            Self::bspline_dw(t[2]),
        );
        let mut v = 0.0f64;
        let mut g = [0.0f64; 3];
        for (dk, wzk) in wz.iter().enumerate() {
            for (dj, wyj) in wy.iter().enumerate() {
                for (di, wxi) in wx.iter().enumerate() {
                    let s = self.sample_raw(
                        base[0] + i64::try_from(di).expect("di<3") - 1,
                        base[1] + i64::try_from(dj).expect("dj<3") - 1,
                        base[2] + i64::try_from(dk).expect("dk<3") - 1,
                    );
                    v += s * wxi * wyj * wzk;
                    g[0] += s * dwx[di] * wyj * wzk;
                    g[1] += s * wxi * dwy[dj] * wzk;
                    g[2] += s * wxi * wyj * dwz[dk];
                }
            }
        }
        (
            v,
            Vec3::new(g[0] / self.h[0], g[1] / self.h[1], g[2] / self.h[2]),
        )
    }

    /// Measured eikonal statistics over a seeded probe set (evidence for
    /// sphere-tracing safety; a MEASUREMENT, not a certificate).
    #[must_use]
    pub fn measure_eikonal(&self, seed: u64, probes: u64, _cx: &Cx<'_>) -> EikonalStats {
        let mut state = seed | 1;
        let mut unit = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / (1u64 << 53) as f64
        };
        let (mut sum, mut max) = (0.0f64, 0.0f64);
        for _ in 0..probes {
            let p = Point3::new(
                self.box_.min.x + (self.box_.max.x - self.box_.min.x) * unit(),
                self.box_.min.y + (self.box_.max.y - self.box_.min.y) * unit(),
                self.box_.min.z + (self.box_.max.z - self.box_.min.z) * unit(),
            );
            let (_, g) = self.spline_eval(p);
            let dev = (g.norm() - 1.0).abs();
            sum += dev;
            max = max.max(dev);
        }
        EikonalStats {
            mean_abs_dev: sum / probes.max(1) as f64,
            max_abs_dev: max,
            probes,
        }
    }

    /// Mean curvature via central second differences at step `h` (an
    /// ESTIMATE — certified curvature stencils are fs-ivl-integration
    /// follow-up work; see CONTRACT no-claims).
    #[must_use]
    pub fn mean_curvature_estimate(&self, x: Point3) -> f64 {
        let h = self.h[0].min(self.h[1]).min(self.h[2]);
        let f = |p: Point3| self.spline_eval(p).0;
        let lap = (f(x.offset(Vec3::new(h, 0.0, 0.0)))
            + f(x.offset(Vec3::new(-h, 0.0, 0.0)))
            + f(x.offset(Vec3::new(0.0, h, 0.0)))
            + f(x.offset(Vec3::new(0.0, -h, 0.0)))
            + f(x.offset(Vec3::new(0.0, 0.0, h)))
            + f(x.offset(Vec3::new(0.0, 0.0, -h)))
            - 6.0 * f(x))
            / (h * h);
        0.5 * lap
    }

    /// Sphere-trace a ray against this chart using its own declared bound
    /// and Lipschitz claim (steps shrink by the safety factor; the
    /// certified-never-tunnel property of plan §10.2). Returns the hit
    /// parameter `t` when the surface is crossed within `t_max`.
    #[must_use]
    pub fn raycast(&self, origin: Point3, dir: Vec3, t_max: f64, cx: &Cx<'_>) -> Option<f64> {
        let n = dir.norm();
        if n <= 0.0 {
            return None;
        }
        let d = dir.scale(1.0 / n);
        let lip = self.chart_lipschitz();
        let mut t = 0.0f64;
        for _ in 0..10_000 {
            if cx.is_cancel_requested() || t > t_max {
                return None;
            }
            let p = origin.offset(d.scale(t));
            let sd = self.spline_eval(p).0;
            if sd <= self.bound {
                return Some(t);
            }
            // Safe step: the field moves at most `lip` per unit distance,
            // and lies within `bound` of the true SDF.
            t += ((sd - self.bound) / lip).max(1e-7);
        }
        None
    }

    /// The conservative Lipschitz claim for the reconstruction: the spline
    /// slope is bounded by the samples' slope, which is at most the
    /// source's constant plus the quantization slack across one cell.
    fn chart_lipschitz(&self) -> f64 {
        let h_min = self.h[0].min(self.h[1]).min(self.h[2]);
        self.source_lipschitz + 2.0 * self.bound / h_min
    }
}

impl Chart for TiledSdf {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let clamped = Point3::new(
            x.x.clamp(self.box_.min.x, self.box_.max.x),
            x.y.clamp(self.box_.min.y, self.box_.max.y),
            x.z.clamp(self.box_.min.z, self.box_.max.z),
        );
        let dist_out = x.delta_from(clamped).norm();
        let (v, g) = self.spline_eval(clamped);
        if dist_out == 0.0 {
            ChartSample {
                signed_distance: v,
                gradient: Some(g),
                lipschitz: Some(self.chart_lipschitz()),
                error: NumericalCertificate::enclosure(v - self.bound, v + self.bound),
            }
        } else {
            let sd = v + dist_out;
            ChartSample {
                signed_distance: sd,
                gradient: None,
                lipschitz: Some(self.chart_lipschitz()),
                error: NumericalCertificate::enclosure(
                    sd - self.bound - (1.0 + self.source_lipschitz) * dist_out,
                    sd + self.bound,
                ),
            }
        }
    }

    fn support(&self) -> Aabb {
        self.box_
    }

    fn name(&self) -> &'static str {
        "rep-sdf/dense"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C1
    }
}
