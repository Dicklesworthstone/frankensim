//! Narrow-band level sets on FrankenVDB (plan §7.2): the storage mode for
//! level-set evolution (topology optimization's workhorse, §9.5). v1
//! provides the storage + maintenance + transport primitives — band
//! construction from any chart, semi-Lagrangian advection with band
//! rebuild, and eikonal reinitialization sweeps — with drift measured and
//! ledgered on fixtures (WENO/fast-iterative upgrades are the
//! topo-levelset bead's, per CONTRACT no-claims).

use crate::dense::{SdfBuildError, finite_positive};
use crate::vdb::VdbGrid;
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ClippedChart, Point3, SamplingDomain, Vec3};
use std::fmt::Write as _;

/// Maximum dense lattice work admitted before narrow-band sparsification.
pub const NARROW_BAND_MAX_SCAN_SAMPLES: u64 = 16_777_216;

/// Largest per-axis sample count whose non-negative VDB coordinates and
/// one-cell maintenance halo remain representable as `i32`.
pub const NARROW_BAND_MAX_SAMPLES_PER_AXIS: u64 = 2_147_483_647;

/// A narrow-band signed-distance function: φ stored on active voxels
/// within `half_width` cells of the zero crossing.
pub struct NarrowBand {
    grid: VdbGrid<f32>,
    /// Voxel edge length (world units).
    h: f64,
    /// World position of voxel (0,0,0)'s center.
    origin: Point3,
    /// Band half-width in cells.
    half_width_cells: u32,
}

/// Band statistics (ledgered evidence for maintenance/drift audits).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandStats {
    /// Active voxels in the band.
    pub active: u64,
    /// Mean |,∇φ| − 1, over interior band voxels (reinit quality).
    pub mean_eikonal_dev: f64,
    /// Max |φ| in the band (bandwidth sanity).
    pub max_abs_phi: f64,
}

impl BandStats {
    /// Canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(80);
        let _ = write!(
            s,
            "{{\"active\":{},\"mean_eikonal_dev\":{:.6},\"max_abs_phi\":{:.6}}}",
            self.active, self.mean_eikonal_dev, self.max_abs_phi
        );
        s
    }
}

fn checked_axis_samples(span: f64, h: f64, axis: usize) -> Result<u64, SdfBuildError> {
    // The sparse grid stores one isotropic scalar `h`, so keep exactly uniform
    // lattice spacing. Flooring plus a max-anchored origin covers the original
    // support because this scan domain was padded by 2h, while never placing a
    // final node beyond the admitted finite AABB.
    let cells = (span / h).floor();
    if !cells.is_finite() || cells < 0.0 || cells > (NARROW_BAND_MAX_SAMPLES_PER_AXIS - 1) as f64 {
        let need = if cells.is_finite() && cells >= 0.0 {
            (cells as u64).saturating_add(1)
        } else {
            u64::MAX
        };
        return Err(SdfBuildError::CoordinateRange {
            axis,
            need,
            cap: NARROW_BAND_MAX_SAMPLES_PER_AXIS,
        });
    }
    (cells as u64)
        .checked_add(1)
        .ok_or(SdfBuildError::CoordinateRange {
            axis,
            need: u64::MAX,
            cap: NARROW_BAND_MAX_SAMPLES_PER_AXIS,
        })
}

fn checked_scan_count(dims: [u64; 3]) -> Result<u64, SdfBuildError> {
    let need = dims
        .iter()
        .try_fold(1u128, |product, &dim| product.checked_mul(u128::from(dim)))
        .unwrap_or(u128::MAX);
    if need > u128::from(NARROW_BAND_MAX_SCAN_SAMPLES) {
        return Err(SdfBuildError::BandScanLimit {
            dims,
            need,
            cap: NARROW_BAND_MAX_SCAN_SAMPLES,
        });
    }
    u64::try_from(need).map_err(|_| SdfBuildError::BandScanLimit {
        dims,
        need,
        cap: NARROW_BAND_MAX_SCAN_SAMPLES,
    })
}

impl NarrowBand {
    /// Build a band of `half_width_cells` around `source`'s zero level set
    /// at voxel size `h`, scanning the source's inflated support.
    /// Polls cancellation at most every 256 source evaluations.
    ///
    /// # Errors
    /// [`SdfBuildError`] when spacing, the finite sampling domain, VDB
    /// coordinates, or dense scan work is inadmissible, or when cancellation
    /// is observed mid-build.
    pub fn from_chart(
        source: &dyn Chart,
        h: f64,
        half_width_cells: u32,
        cx: &Cx<'_>,
    ) -> Result<NarrowBand, SdfBuildError> {
        let h = finite_positive(h, "h")?;
        let cutoff = f64::from(half_width_cells) * h;
        if !cutoff.is_finite() {
            return Err(SdfBuildError::InvalidSpacing {
                field: "half_width_cells * h",
                value: cutoff,
            });
        }
        let padding = finite_positive(2.0 * h, "2 * h")?;
        let raw_support = SamplingDomain::admit(source.support(), None)?.bounds();
        let support_domain = SamplingDomain::admit(raw_support.inflate(padding), None)?;
        let support = support_domain.bounds();
        let spans = support_domain.spans();
        let dims = [
            checked_axis_samples(spans.x, h, 0)?,
            checked_axis_samples(spans.y, h, 1)?,
            checked_axis_samples(spans.z, h, 2)?,
        ];
        let _scan_count = checked_scan_count(dims)?;
        let n = dims.map(|dim| i32::try_from(dim).expect("coordinate cap checked"));
        let origin = Point3::new(
            support.max.x - f64::from(n[0] - 1) * h,
            support.max.y - f64::from(n[1] - 1) * h,
            support.max.z - f64::from(n[2] - 1) * h,
        );
        let mut grid = VdbGrid::new(f32::MAX);
        let mut completed_samples = 0u64;
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    if completed_samples.is_multiple_of(256) {
                        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
                    }
                    let p = Point3::new(
                        origin.x + f64::from(i) * h,
                        origin.y + f64::from(j) * h,
                        origin.z + f64::from(k) * h,
                    );
                    let sd = source.eval(p, cx).signed_distance;
                    let stored = sd as f32;
                    if !sd.is_finite() || !stored.is_finite() {
                        return Err(SdfBuildError::InvalidSample {
                            point: p,
                            value_bits: sd.to_bits(),
                        });
                    }
                    completed_samples += 1;
                    if sd.abs() <= cutoff {
                        grid.set([i, j, k], stored);
                    }
                }
            }
        }
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        Ok(NarrowBand {
            grid,
            h,
            origin,
            half_width_cells,
        })
    }

    /// Build a narrow band of the geometric intersection `source ∩ clip`.
    /// The clip and dense scan work are validated before any source
    /// evaluation or sparse-grid population.
    ///
    /// # Errors
    /// [`SdfBuildError`] under the same conditions as [`Self::from_chart`],
    /// plus an invalid, empty, or degenerate explicit clip.
    pub fn from_chart_clipped(
        source: &dyn Chart,
        h: f64,
        half_width_cells: u32,
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<NarrowBand, SdfBuildError> {
        let clipped = ClippedChart::new(source, clip)?;
        Self::from_chart(&clipped, h, half_width_cells, cx)
    }

    /// Voxel size.
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// The underlying sparse grid (read access for consumers/tests).
    #[must_use]
    pub fn grid(&self) -> &VdbGrid<f32> {
        &self.grid
    }

    /// Mutable grid access (test fixtures and the level-set evolution bead
    /// distort/edit φ directly; band invariants are re-established by
    /// `reinitialize`/`rebuild`).
    pub fn grid_mut(&mut self) -> &mut VdbGrid<f32> {
        &mut self.grid
    }

    fn world(&self, c: [i32; 3]) -> Point3 {
        Point3::new(
            self.origin.x + f64::from(c[0]) * self.h,
            self.origin.y + f64::from(c[1]) * self.h,
            self.origin.z + f64::from(c[2]) * self.h,
        )
    }

    /// Trilinear φ at a world point; `None` when any stencil corner is
    /// outside the band (callers treat that as "far").
    #[must_use]
    pub fn interpolate(&self, p: Point3) -> Option<f64> {
        let fx = (p.x - self.origin.x) / self.h;
        let fy = (p.y - self.origin.y) / self.h;
        let fz = (p.z - self.origin.z) / self.h;
        let (i0, j0, k0) = (fx.floor() as i32, fy.floor() as i32, fz.floor() as i32);
        let (tx, ty, tz) = (fx - f64::from(i0), fy - f64::from(j0), fz - f64::from(k0));
        let mut c = [0.0f64; 8];
        for (idx, slot) in c.iter_mut().enumerate() {
            let d = [
                i0 + i32::from(idx & 1 == 1),
                j0 + i32::from((idx >> 1) & 1 == 1),
                k0 + i32::from((idx >> 2) & 1 == 1),
            ];
            if !self.grid.is_active(d) {
                return None;
            }
            *slot = f64::from(self.grid.get(d));
        }
        let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
        let c00 = lerp(c[0], c[1], tx);
        let c10 = lerp(c[2], c[3], tx);
        let c01 = lerp(c[4], c[5], tx);
        let c11 = lerp(c[6], c[7], tx);
        Some(lerp(lerp(c00, c10, ty), lerp(c01, c11, ty), tz))
    }

    /// One semi-Lagrangian advection step: φⁿ⁺¹(x) = φⁿ(x − v(x)·dt),
    /// then band rebuild (dilate + trim). First-order, unconditionally
    /// stable; drift is measured, not assumed (rsdf-005).
    pub fn advect(&mut self, velocity: impl Fn(Point3) -> Vec3, dt: f64) {
        let snapshot: Vec<([i32; 3], f32)> = self.grid.iter_active().collect();
        // Dilate once so the band can follow the interface.
        self.grid.dilate();
        let targets: Vec<[i32; 3]> = self.grid.iter_active().map(|(c, _)| c).collect();
        let old = {
            let mut g = VdbGrid::new(f32::MAX);
            for (c, v) in snapshot {
                g.set(c, v);
            }
            g
        };
        let old_band = NarrowBand {
            grid: old,
            h: self.h,
            origin: self.origin,
            half_width_cells: self.half_width_cells,
        };
        let cutoff = (f64::from(self.half_width_cells) * self.h) as f32;
        for c in targets {
            let p = self.world(c);
            let v = velocity(p);
            let back = p.offset(v.scale(-dt));
            match old_band.interpolate(back) {
                Some(phi) if phi.abs() as f32 <= cutoff => self.grid.set(c, phi as f32),
                _ => self.grid.deactivate(c),
            }
        }
    }

    /// Eikonal reinitialization sweeps: relax φ toward |∇φ| = 1 with the
    /// upwind Godunov gradient and sub-CFL pseudo-time steps, keeping the
    /// zero crossing pinned by the smoothed sign function.
    pub fn reinitialize(&mut self, sweeps: u32) {
        let dtau = 0.3 * self.h;
        for _ in 0..sweeps {
            let snapshot: Vec<([i32; 3], f32)> = self.grid.iter_active().collect();
            let read = |c: [i32; 3], fallback: f32| -> f64 {
                if self.grid.is_active(c) {
                    f64::from(self.grid.get(c))
                } else {
                    f64::from(fallback)
                }
            };
            let mut updates = Vec::with_capacity(snapshot.len());
            for (c, v0) in &snapshot {
                let phi = f64::from(*v0);
                let g = self.upwind_gradient_norm(*c, *v0, &read);
                let sign = phi / (phi * phi + self.h * self.h).sqrt();
                updates.push((*c, (phi - dtau * sign * (g - 1.0)) as f32));
            }
            for (c, v) in updates {
                self.grid.set(c, v);
            }
        }
    }

    fn upwind_gradient_norm(
        &self,
        c: [i32; 3],
        v0: f32,
        read: &impl Fn([i32; 3], f32) -> f64,
    ) -> f64 {
        let phi = f64::from(v0);
        let mut acc = 0.0f64;
        for axis in 0..3 {
            let mut plus = c;
            plus[axis] += 1;
            let mut minus = c;
            minus[axis] -= 1;
            let dp = (read(plus, v0) - phi) / self.h;
            let dm = (phi - read(minus, v0)) / self.h;
            // Godunov upwinding keyed by the sign of φ.
            let d = if phi >= 0.0 {
                dm.max(0.0).powi(2).max(dp.min(0.0).powi(2))
            } else {
                dp.max(0.0).powi(2).max(dm.min(0.0).powi(2))
            };
            acc += d;
        }
        acc.sqrt()
    }

    /// Band statistics (measured; ledgered by callers/tests).
    #[must_use]
    pub fn stats(&self) -> BandStats {
        let mut n = 0u64;
        let mut dev_sum = 0.0f64;
        let mut interior = 0u64;
        let mut max_abs = 0.0f64;
        let read = |c: [i32; 3], fallback: f32| -> f64 {
            if self.grid.is_active(c) {
                f64::from(self.grid.get(c))
            } else {
                f64::from(fallback)
            }
        };
        for (c, v) in self.grid.iter_active() {
            n += 1;
            max_abs = max_abs.max(f64::from(v).abs());
            // Only fully-interior voxels give meaningful eikonal residuals.
            let all_in = (0..3).all(|axis| {
                let mut p = c;
                p[axis] += 1;
                let mut m = c;
                m[axis] -= 1;
                self.grid.is_active(p) && self.grid.is_active(m)
            });
            if all_in {
                interior += 1;
                dev_sum += (self.upwind_gradient_norm(c, v, &read) - 1.0).abs();
            }
        }
        BandStats {
            active: n,
            mean_eikonal_dev: dev_sum / interior.max(1) as f64,
            max_abs_phi: max_abs,
        }
    }
}
