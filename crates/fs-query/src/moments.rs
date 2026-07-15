//! Certified geometric moments (bead rjnd, E1 query upgrades, part 1).
//!
//! Rigorous enclosures of volume, first moments, and origin-referenced
//! second moments of the region presented by a chart, computed by
//! certified cell classification over the chart's exact-distance field
//! (the [`TraceStepClaim::ExactDistance`] capability). The pattern
//! generalizes fs-geocon's `volume_certified`: cells whose inflated
//! distance enclosure proves them inside contribute exact closed-form
//! box integrals with outward rounding; cells straddling the boundary
//! contribute conservative brackets; everything else contributes
//! nothing. A chart that cannot certify its field REFUSES instead of
//! guessing, per the capability-routing doctrine.
//!
//! Geometry never knows materials: these are unit-density geometric
//! integrals. Rigid-body mass properties are minted downstream by
//! combining [`GeometricMoments`] with a density receipt (fs-matdb),
//! and inertia about other frames follows from [`GeometricMoments::
//! translated`], whose covariance laws are outward-rounded exact
//! identities.

use crate::{ContactInflation, QueryError};
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, TraceStepClaim};

/// Hard bound on quadrature cells per moments call: refusal, not
/// truncation, beyond this (the caller owns coarsening or splitting).
pub const MAX_MOMENT_CELLS: u64 = 1 << 22;

/// Cancellation-poll stride in cells.
const CHECKPOINT_STRIDE: u64 = 256;

/// One outward-rounded scalar enclosure `[lo, hi]` of a real integral.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MomentEnclosure {
    /// Certified lower bound.
    pub lo: f64,
    /// Certified upper bound.
    pub hi: f64,
}

impl MomentEnclosure {
    const ZERO: MomentEnclosure = MomentEnclosure { lo: 0.0, hi: 0.0 };

    /// Whether the enclosure contains `value`.
    #[must_use]
    pub fn contains(&self, value: f64) -> bool {
        self.lo <= value && value <= self.hi
    }

    /// Enclosure width (`hi - lo`, rounded up).
    #[must_use]
    pub fn width(&self) -> f64 {
        (self.hi - self.lo).next_up()
    }

    /// Whether two enclosures of the same real value are consistent
    /// (a sound metamorphic: both contain the truth, so they overlap).
    #[must_use]
    pub fn overlaps(&self, other: &MomentEnclosure) -> bool {
        self.lo <= other.hi && other.lo <= self.hi
    }

    fn add(self, other: MomentEnclosure) -> MomentEnclosure {
        MomentEnclosure {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }

    /// Interval product, outward rounded (general signs).
    fn mul(self, other: MomentEnclosure) -> MomentEnclosure {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for p in products {
            lo = lo.min(p);
            hi = hi.max(p);
        }
        MomentEnclosure {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }

    fn scalar(value: f64) -> MomentEnclosure {
        MomentEnclosure {
            lo: value.next_down().next_down(),
            hi: value.next_up().next_up(),
        }
    }
}

/// Origin-referenced second-moment enclosures (the inertia integrand
/// family `∫ x_i x_j dV`, NOT the inertia tensor itself).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondMoments {
    /// `∫ x² dV`.
    pub xx: MomentEnclosure,
    /// `∫ y² dV`.
    pub yy: MomentEnclosure,
    /// `∫ z² dV`.
    pub zz: MomentEnclosure,
    /// `∫ x·y dV`.
    pub xy: MomentEnclosure,
    /// `∫ x·z dV`.
    pub xz: MomentEnclosure,
    /// `∫ y·z dV`.
    pub yz: MomentEnclosure,
}

/// Certified unit-density geometric moments of a chart's region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometricMoments {
    /// `∫ 1 dV` enclosure.
    pub volume: MomentEnclosure,
    /// `[∫x dV, ∫y dV, ∫z dV]` enclosures.
    pub first: [MomentEnclosure; 3],
    /// Origin-referenced second-moment enclosures.
    pub second: SecondMoments,
    /// Requested cell-spacing bound the grid honored.
    pub h: f64,
    /// Cells proven fully inside.
    pub sure_cells: u64,
    /// Boundary-straddling cells bracketed conservatively.
    pub band_cells: u64,
}

impl GeometricMoments {
    /// Center-of-mass enclosure `first / volume` per axis.
    ///
    /// # Errors
    /// [`QueryError::MomentsVolumeUnproven`] unless the certified
    /// volume lower bound is strictly positive (a region that may be
    /// empty has no certifiable center).
    pub fn com_enclosure(&self) -> Result<[MomentEnclosure; 3], QueryError> {
        if !(self.volume.lo.is_finite() && self.volume.lo > 0.0) {
            return Err(QueryError::MomentsVolumeUnproven {
                volume_lo: self.volume.lo,
            });
        }
        let inv = MomentEnclosure {
            lo: (1.0 / self.volume.hi).next_down(),
            hi: (1.0 / self.volume.lo).next_up(),
        };
        Ok([
            self.first[0].mul(inv),
            self.first[1].mul(inv),
            self.first[2].mul(inv),
        ])
    }

    /// Exact translation covariance, outward rounded: the moments of
    /// the same region shifted by `d` (equivalently, of the original
    /// region measured in a frame whose origin moved by `-d`). The
    /// parallel-axis law is the diagonal of this identity:
    /// `∫(x+dx)² = ∫x² + 2·dx·∫x + dx²·V`.
    #[must_use]
    pub fn translated(&self, d: [f64; 3]) -> GeometricMoments {
        let dx = MomentEnclosure::scalar(d[0]);
        let dy = MomentEnclosure::scalar(d[1]);
        let dz = MomentEnclosure::scalar(d[2]);
        let two = MomentEnclosure { lo: 2.0, hi: 2.0 };
        let first = [
            self.first[0].add(dx.mul(self.volume)),
            self.first[1].add(dy.mul(self.volume)),
            self.first[2].add(dz.mul(self.volume)),
        ];
        let diag = |m: MomentEnclosure, m1: MomentEnclosure, dd: MomentEnclosure| {
            m.add(two.mul(dd).mul(m1)).add(dd.mul(dd).mul(self.volume))
        };
        let cross = |m: MomentEnclosure,
                     ma: MomentEnclosure,
                     mb: MomentEnclosure,
                     da: MomentEnclosure,
                     db: MomentEnclosure| {
            m.add(da.mul(mb))
                .add(db.mul(ma))
                .add(da.mul(db).mul(self.volume))
        };
        GeometricMoments {
            volume: self.volume,
            first,
            second: SecondMoments {
                xx: diag(self.second.xx, self.first[0], dx),
                yy: diag(self.second.yy, self.first[1], dy),
                zz: diag(self.second.zz, self.first[2], dz),
                xy: cross(self.second.xy, self.first[0], self.first[1], dx, dy),
                xz: cross(self.second.xz, self.first[0], self.first[2], dx, dz),
                yz: cross(self.second.yz, self.first[1], self.first[2], dy, dz),
            },
            h: self.h,
            sure_cells: self.sure_cells,
            band_cells: self.band_cells,
        }
    }
}

struct MomentAccumulator {
    volume: MomentEnclosure,
    first: [MomentEnclosure; 3],
    diag: [MomentEnclosure; 3],
    cross: [MomentEnclosure; 3],
}

/// Certified geometric moments of `chart`'s region over `domain`.
///
/// `domain` must contain the chart's declared support box: moments are
/// claims about the WHOLE region, so a domain that could exclude part
/// of it refuses instead of silently integrating a subset. `h` bounds
/// the cell spacing; the grid uses `ceil(extent/h)` cells per axis.
///
/// # Errors
/// [`QueryError::MomentsUncertifiedChart`] unless the chart claims
/// `ExactDistance`; [`QueryError::MomentsInvalidDomain`] for
/// non-finite/inverted domains or a domain that does not contain the
/// chart support; [`QueryError::MomentsInvalidSpacing`] for a
/// non-finite or non-positive `h`; [`QueryError::MomentsExcessiveWork`]
/// beyond [`MAX_MOMENT_CELLS`]; [`QueryError::MomentsInvalidSample`]
/// when a sample's enclosure is missing, non-finite, inverted, or only
/// Estimate/NoClaim class; [`QueryError::Cancelled`] on cancellation.
pub fn geometric_moments(
    chart: &dyn Chart,
    domain: &Aabb,
    h: f64,
    cx: &Cx<'_>,
) -> Result<GeometricMoments, QueryError> {
    geometric_moments_with_inflation(chart, domain, h, cx, ContactInflation::exact_zero())
}

/// Certified moments with a proof-bearing absolute representation error added
/// to every cell's exact-distance classification radius.
///
/// The uncertainty can only turn sure inside/outside cells into conservative
/// boundary-band cells. It cannot tighten an integral enclosure. Exact zero
/// delegates with the original classification-radius bits unchanged.
///
/// # Errors
/// The refusals from [`geometric_moments`], plus
/// [`QueryError::InvalidContactInflation`] if the cell radius and retained
/// uncertainty cannot be combined in the finite domain.
pub fn geometric_moments_with_inflation(
    chart: &dyn Chart,
    domain: &Aabb,
    h: f64,
    cx: &Cx<'_>,
    inflation: ContactInflation,
) -> Result<GeometricMoments, QueryError> {
    geometric_moments_impl(chart, domain, h, cx, inflation)
}

fn geometric_moments_impl(
    chart: &dyn Chart,
    domain: &Aabb,
    h: f64,
    cx: &Cx<'_>,
    inflation: ContactInflation,
) -> Result<GeometricMoments, QueryError> {
    let claim = chart.trace_step_claim();
    if claim != TraceStepClaim::ExactDistance {
        return Err(QueryError::MomentsUncertifiedChart { claim });
    }
    if !h.is_finite() || h <= 0.0 {
        return Err(QueryError::MomentsInvalidSpacing {
            spacing_bits: h.to_bits(),
        });
    }
    let min = [domain.min.x, domain.min.y, domain.min.z];
    let max = [domain.max.x, domain.max.y, domain.max.z];
    for a in 0..3 {
        if !(min[a].is_finite() && max[a].is_finite() && min[a] < max[a]) {
            return Err(QueryError::MomentsInvalidDomain {
                detail: "domain must be finite with min strictly below max on every axis",
            });
        }
    }
    let support = chart.support();
    let smin = [support.min.x, support.min.y, support.min.z];
    let smax = [support.max.x, support.max.y, support.max.z];
    for a in 0..3 {
        if !(smin[a] >= min[a] && smax[a] <= max[a]) {
            return Err(QueryError::MomentsInvalidDomain {
                detail: "domain must contain the chart support box \
                         (moments are whole-region claims)",
            });
        }
    }
    let (dims, width, radius, cell_vol) = moment_grid(&min, &max, h)?;
    let classification_radius = inflation.inflate_nonnegative(radius)?;
    let mut acc = MomentAccumulator {
        volume: MomentEnclosure::ZERO,
        first: [MomentEnclosure::ZERO; 3],
        diag: [MomentEnclosure::ZERO; 3],
        cross: [MomentEnclosure::ZERO; 3],
    };
    let mut sure_cells = 0u64;
    let mut band_cells = 0u64;
    let mut completed = 0u64;
    for i in 0..dims[0] {
        for j in 0..dims[1] {
            for k in 0..dims[2] {
                if completed.is_multiple_of(CHECKPOINT_STRIDE) && cx.checkpoint().is_err() {
                    return Err(QueryError::Cancelled);
                }
                completed += 1;
                let index = [i, j, k];
                let mut center = [0f64; 3];
                let mut cell_lo = [0f64; 3];
                let mut cell_hi = [0f64; 3];
                for a in 0..3 {
                    #[allow(clippy::cast_precision_loss)]
                    let base = index[a] as f64;
                    cell_lo[a] = (min[a] + base * width[a]).next_down().next_down();
                    cell_hi[a] = (min[a] + (base + 1.0) * width[a]).next_up().next_up();
                    center[a] = min[a] + (base + 0.5) * width[a];
                }
                let p = fs_geom::Point3::new(center[0], center[1], center[2]);
                let sample = chart.eval(p, cx);
                let enclosure = chart.trace_value_enclosure(p, &sample, cx);
                validate_moment_enclosure(&enclosure, center)?;
                let field_hi = (enclosure.hi + classification_radius).next_up();
                let field_lo = (enclosure.lo - classification_radius).next_down();
                if field_hi <= 0.0 {
                    sure_cells += 1;
                    accumulate_sure(&mut acc, cell_vol, center, width);
                } else if field_lo < 0.0 {
                    band_cells += 1;
                    accumulate_band(&mut acc, cell_vol.hi, cell_lo, cell_hi);
                }
            }
        }
    }
    if cx.checkpoint().is_err() {
        return Err(QueryError::Cancelled);
    }
    Ok(GeometricMoments {
        volume: acc.volume,
        first: acc.first,
        second: SecondMoments {
            xx: acc.diag[0],
            yy: acc.diag[1],
            zz: acc.diag[2],
            xy: acc.cross[0],
            xz: acc.cross[1],
            yz: acc.cross[2],
        },
        h,
        sure_cells,
        band_cells,
    })
}

/// Grid sizing under the deterministic work ceiling, plus the L = 1
/// cell inflation radius (circumradius, rounded up) and outward-rounded
/// cell volume bounds.
#[allow(clippy::type_complexity)]
fn moment_grid(
    min: &[f64; 3],
    max: &[f64; 3],
    h: f64,
) -> Result<([u64; 3], [f64; 3], f64, MomentEnclosure), QueryError> {
    let mut dims = [0u64; 3];
    let mut width = [0f64; 3];
    let mut cells: u64 = 1;
    for a in 0..3 {
        let extent = max[a] - min[a];
        let n = (extent / h).ceil();
        if !n.is_finite() || n < 0.0 {
            return Err(QueryError::MomentsInvalidSpacing {
                spacing_bits: h.to_bits(),
            });
        }
        // Refuse before sizing anything: the total below can only grow,
        // so an early per-axis explosion is final.
        #[allow(clippy::cast_precision_loss)]
        if n > MAX_MOMENT_CELLS as f64 {
            return Err(QueryError::MomentsExcessiveWork {
                max_cells: MAX_MOMENT_CELLS,
            });
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = (n as u64).max(1);
        dims[a] = n;
        #[allow(clippy::cast_precision_loss)]
        let w = extent / n as f64;
        width[a] = w;
        cells = cells.saturating_mul(n);
    }
    if cells > MAX_MOMENT_CELLS {
        return Err(QueryError::MomentsExcessiveWork {
            max_cells: MAX_MOMENT_CELLS,
        });
    }
    let radius = {
        let q = |w: f64| {
            let half = 0.5 * w;
            (half * half).next_up()
        };
        (q(width[0]) + q(width[1]) + q(width[2]))
            .next_up()
            .sqrt()
            .next_up()
    };
    let cell_vol = MomentEnclosure {
        lo: ((width[0] * width[1]).next_down() * width[2]).next_down(),
        hi: ((width[0] * width[1]).next_up() * width[2]).next_up(),
    };
    Ok((dims, width, radius, cell_vol))
}

fn validate_moment_enclosure(
    enclosure: &NumericalCertificate,
    at: [f64; 3],
) -> Result<(), QueryError> {
    let sound = matches!(
        enclosure.kind,
        NumericalKind::Exact | NumericalKind::Enclosure
    ) && enclosure.lo.is_finite()
        && enclosure.hi.is_finite()
        && enclosure.lo <= enclosure.hi;
    if sound {
        Ok(())
    } else {
        Err(QueryError::MomentsInvalidSample { at })
    }
}

/// Exact closed-form box integrals for a cell proven fully inside,
/// outward rounded: `∫1 = V`, `∫x = V·c_x`, `∫x² = V·(c_x² + w_x²/12)`,
/// `∫xy = V·c_x·c_y`.
fn accumulate_sure(
    acc: &mut MomentAccumulator,
    cell_vol: MomentEnclosure,
    center: [f64; 3],
    width: [f64; 3],
) {
    acc.volume = acc.volume.add(cell_vol);
    let c = [
        MomentEnclosure::scalar(center[0]),
        MomentEnclosure::scalar(center[1]),
        MomentEnclosure::scalar(center[2]),
    ];
    for a in 0..3 {
        acc.first[a] = acc.first[a].add(cell_vol.mul(c[a]));
        let w2_12 = MomentEnclosure::scalar(width[a] * width[a] / 12.0);
        acc.diag[a] = acc.diag[a].add(cell_vol.mul(c[a].mul(c[a]).add(w2_12)));
    }
    let pairs = [(0usize, 1usize), (0, 2), (1, 2)];
    for (slot, (a, b)) in pairs.into_iter().enumerate() {
        acc.cross[slot] = acc.cross[slot].add(cell_vol.mul(c[a].mul(c[b])));
    }
}

/// Conservative brackets for a boundary-straddling cell: the region
/// mass inside is between 0 and the full cell volume, and each
/// integrand ranges over its cell extremes, so every contribution is
/// bracketed by the worst signed combination (zero mass always being
/// admissible).
fn accumulate_band(acc: &mut MomentAccumulator, vol_hi: f64, cell_lo: [f64; 3], cell_hi: [f64; 3]) {
    acc.volume = acc.volume.add(MomentEnclosure {
        lo: 0.0,
        hi: vol_hi,
    });
    for a in 0..3 {
        acc.first[a] = acc.first[a].add(MomentEnclosure {
            lo: (vol_hi * cell_lo[a].min(0.0)).next_down(),
            hi: (vol_hi * cell_hi[a].max(0.0)).next_up(),
        });
        let q_max = (cell_lo[a] * cell_lo[a])
            .max(cell_hi[a] * cell_hi[a])
            .next_up();
        acc.diag[a] = acc.diag[a].add(MomentEnclosure {
            lo: 0.0,
            hi: (vol_hi * q_max).next_up(),
        });
    }
    let pairs = [(0usize, 1usize), (0, 2), (1, 2)];
    for (slot, (a, b)) in pairs.into_iter().enumerate() {
        let corners = [
            cell_lo[a] * cell_lo[b],
            cell_lo[a] * cell_hi[b],
            cell_hi[a] * cell_lo[b],
            cell_hi[a] * cell_hi[b],
        ];
        let mut p_min = 0.0f64;
        let mut p_max = 0.0f64;
        for p in corners {
            p_min = p_min.min(p);
            p_max = p_max.max(p);
        }
        acc.cross[slot] = acc.cross[slot].add(MomentEnclosure {
            lo: (vol_hi * p_min).next_down(),
            hi: (vol_hi * p_max).next_up(),
        });
    }
}
