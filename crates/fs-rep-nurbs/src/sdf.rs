//! CONVERTER NURBS → SDF (plan §7.3 edge 3, bead wqd.11; [F] — behind
//! the `nurbs-sdf` feature until its Gauntlet tier is green): measured distance
//! estimates to trimmed NURBS shells. Bézier control hulls drive a useful
//! branch-and-bound bracket and damped Gauss–Newton improves an evaluated-point
//! estimate, but Cartesian division, hull inflation, surface evaluation, and
//! distance arithmetic are ordinary f64 operations rather than outward-rounded
//! enclosures. The chart therefore emits `Estimate`, no Lipschitz authority,
//! and an explicitly estimated name. Trim classification can further widen the
//! estimate. A successor interval/Taylor path owns certified distance and sign.

use crate::NurbsError;
use crate::closest::{
    CLOSEST_MAX_BASE_WORK_UNITS, CLOSEST_MAX_SPLITS, norm3, preflight_surface_subdivision,
    surface_base_work_units, surface_subdivision_work_per_split,
};
use crate::rat::Rat;
use crate::surface::{AdmittedNurbsSurface, NurbsSurface};
use crate::trim::{Classification, TRIM_CLASSIFY_MAX_WORK_UNITS, TrimmedPatch};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{Aabb, BettiBounds, Chart, ChartSample, Differentiability, Point3, Vec3};
use fs_math::{next_down, next_up};

/// Sign policy for the generated field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Surface normals (du × dv) point OUTWARD (B-rep topology says so):
    /// the lower-authority distance estimate uses an orientation-based sign.
    Outward,
    /// No orientation claim: the field is UNSIGNED (all non-negative)
    /// and the chart name says so.
    Unknown,
}

/// One measured distance-bracket query answer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SdfQuery {
    /// Convex-hull lower estimate with heuristic f64 inflation.
    pub lower: f64,
    /// Evaluated-point distance estimate.
    pub upper: f64,
    /// The best parameter found (u, v).
    pub param: [f64; 2],
    /// Which shell surface owned the minimum.
    pub surface: usize,
    /// The closest point fell outside the kept trim region (or in the
    /// boundary band): `upper` is infinite because no kept-surface witness was
    /// established. `param` and `surface` remain diagnostic only.
    pub trim_downgrade: bool,
    /// Branch-and-bound splits spent (throughput evidence).
    pub splits: u32,
}

/// A NURBS shell presented as a measured distance-field approximation.
#[derive(Debug)]
pub struct ShellSdf {
    surfaces: Vec<NurbsSurface<f64>>,
    trims: Vec<Option<TrimmedPatch>>,
    orientation: Orientation,
    base_work_units: u128,
    trimmed_surface_count: u128,
    split_work_units_per_round: u128,
    polish_work_units: u128,
    sign_work_units: u128,
}

/// Gauss–Newton polish iterations (evaluated-distance improvement only).
const POLISH_STEPS: usize = 8;

/// Trim-classification grid. The entire closed cell containing the floating
/// parameter is classified; no point sample stands in for the original value.
const TRIM_SCALE: i128 = 1 << 20;
const TRIM_SCALE_F64: f64 = 1_048_576.0;
const MAX_EXACT_GRID_INDEX: f64 = 9_007_199_254_740_991.0;

/// Defensive sample ceiling for one legacy tile allocation.
const SDF_TILE_MAX_SAMPLES: usize = 16_777_216;

/// Defensive worst-case total-work ceiling for one legacy tile request.
const SDF_TILE_MAX_WORST_CASE_WORK: u128 = 1_073_741_824;

/// Defensive bound for one shell query including immutable Bézier seed work,
/// requested splits, polish/sign fallback, and exact trim classification.
const SDF_QUERY_MAX_WORK_UNITS: u128 = 67_108_864;

/// Fixed cap on deterministic regular-witness candidates used when the closest
/// parameter lies on a singular surface chart.
const SIGN_FALLBACK_MAX_CANDIDATES: u128 = 80;

fn checked_work_product(values: &[u128], stage: &str) -> Result<u128, NurbsError> {
    values.iter().try_fold(1u128, |acc, value| {
        acc.checked_mul(*value).ok_or_else(|| NurbsError::Domain {
            what: format!("NURBS shell {stage} work accounting overflows u128"),
        })
    })
}

fn validate_distance_request(q: [f64; 3], tol: f64) -> Result<(), NurbsError> {
    if q.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(NurbsError::Domain {
            what: "NURBS SDF query point must be finite".to_string(),
        });
    }
    if !tol.is_finite() || tol < 0.0 {
        return Err(NurbsError::Domain {
            what: "NURBS SDF tolerance must be finite and non-negative".to_string(),
        });
    }
    Ok(())
}

/// Conservative scalar-operation coefficients for one surface after its
/// immutable Bézier-conversion cost. The legacy algorithms operate on whole
/// isocurve nets during partial evaluation, so charging only one unit per
/// split/polish/candidate would leave degree and structure size unbounded.
fn surface_evaluation_work(surface: &NurbsSurface<f64>) -> Result<u128, NurbsError> {
    let nu = surface.cpw.len() as u128;
    let nv = surface.cpw.first().map_or(0, Vec::len) as u128;
    let input_controls = checked_work_product(&[nu, nv], "input-control")?;
    let evaluation_work = input_controls
        .checked_add(surface.knots_u.knots.len() as u128)
        .and_then(|work| work.checked_add(surface.knots_v.knots.len() as u128))
        .and_then(|work| work.checked_mul(16))
        .ok_or_else(|| NurbsError::Domain {
            what: "NURBS shell differential-evaluation work accounting overflows u128".to_string(),
        })?;
    Ok(evaluation_work)
}

impl ShellSdf {
    /// A shell from surfaces + optional trims (parallel arrays).
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on length mismatch or an empty shell.
    pub fn new(
        surfaces: Vec<NurbsSurface<f64>>,
        trims: Vec<Option<TrimmedPatch>>,
        orientation: Orientation,
    ) -> Result<ShellSdf, NurbsError> {
        if surfaces.is_empty() {
            return Err(NurbsError::Structure {
                what: "a shell needs at least one surface".to_string(),
            });
        }
        if surfaces.len() != trims.len() {
            return Err(NurbsError::Structure {
                what: format!(
                    "{} surfaces but {} trim slots (parallel arrays)",
                    surfaces.len(),
                    trims.len()
                ),
            });
        }
        let mut trim_validation_work_remaining = SDF_QUERY_MAX_WORK_UNITS;
        if (trims.len() as u128) > trim_validation_work_remaining {
            return Err(NurbsError::Domain {
                what: format!(
                    "NURBS shell has {} trim slots above the defensive construction scan ceiling {trim_validation_work_remaining}",
                    trims.len()
                ),
            });
        }
        for trim in trims.iter().flatten() {
            trim.validate_live_with_budget(&mut trim_validation_work_remaining)?;
        }
        let mut base_work_units = 0u128;
        let mut split_work_units_per_round = 0u128;
        let mut polish_work_units = 0u128;
        let mut max_evaluation_work = 0u128;
        for surface in &surfaces {
            let work = surface_base_work_units(surface)?;
            if work > CLOSEST_MAX_BASE_WORK_UNITS {
                return Err(NurbsError::Domain {
                    what: format!(
                        "NURBS shell surface base work {work} exceeds defensive closest-point ceiling {CLOSEST_MAX_BASE_WORK_UNITS}"
                    ),
                });
            }
            base_work_units =
                base_work_units
                    .checked_add(work)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "NURBS shell base-work accounting overflows u128".to_string(),
                    })?;
            let split_work = surface_subdivision_work_per_split(surface)?;
            let evaluation_work = surface_evaluation_work(surface)?;
            split_work_units_per_round = split_work_units_per_round
                .checked_add(split_work)
                .ok_or_else(|| NurbsError::Domain {
                    what: "NURBS shell split-work coefficient overflows u128".to_string(),
                })?;
            // One partial plus up to four damped candidate evaluations per
            // polishing step, for every surface in the shell.
            polish_work_units = polish_work_units
                .checked_add(checked_work_product(
                    &[evaluation_work, POLISH_STEPS as u128, 5],
                    "Gauss-Newton polish",
                )?)
                .ok_or_else(|| NurbsError::Domain {
                    what: "NURBS shell polish-work accounting overflows u128".to_string(),
                })?;
            max_evaluation_work = max_evaluation_work.max(evaluation_work);
        }
        let trimmed_surface_count = trims.iter().filter(|trim| trim.is_some()).count() as u128;
        // The winning surface incurs one ordinary gradient/sign evaluation and
        // may evaluate every fixed regular-witness candidate twice (position,
        // then partials). Charge the largest surface coefficient.
        let sign_multiplier = match orientation {
            Orientation::Unknown => 1,
            Orientation::Outward => SIGN_FALLBACK_MAX_CANDIDATES
                .checked_mul(2)
                .and_then(|work| work.checked_add(1))
                .ok_or_else(|| NurbsError::Domain {
                    what: "NURBS shell sign-candidate accounting overflows u128".to_string(),
                })?,
        };
        let sign_work_units =
            checked_work_product(&[max_evaluation_work, sign_multiplier], "sign fallback")?;
        let shell = ShellSdf {
            surfaces,
            trims,
            orientation,
            base_work_units,
            trimmed_surface_count,
            split_work_units_per_round,
            polish_work_units,
            sign_work_units,
        };
        // Do not successfully construct a shell for which even a zero-split
        // query is deterministically inadmissible.
        let _ = shell.query_work_units(0)?;
        Ok(shell)
    }

    fn query_work_units(&self, max_splits: u32) -> Result<u128, NurbsError> {
        if max_splits > CLOSEST_MAX_SPLITS {
            return Err(NurbsError::Domain {
                what: "NURBS shell split request exceeds the defensive legacy-path ceiling"
                    .to_string(),
            });
        }
        // Keep shell admission exactly aligned with the delegated closest
        // solver's per-surface split-work and retained-frontier envelope. The
        // surfaces execute sequentially, so each must fit the per-query cap;
        // summing their retained bytes would overstate simultaneous memory.
        for surface in &self.surfaces {
            preflight_surface_subdivision(surface, max_splits)?;
        }
        let split_work = u128::from(max_splits)
            .checked_mul(self.split_work_units_per_round)
            .ok_or_else(|| NurbsError::Domain {
                what: "NURBS shell split-work accounting overflows u128".to_string(),
            })?;
        // A trimmed query may have to classify both the retained B&B witness
        // and a distinct optional polished witness. Charge both before work.
        let trim_work = checked_work_product(
            &[self.trimmed_surface_count, TRIM_CLASSIFY_MAX_WORK_UNITS, 2],
            "trim classification",
        )?;
        let requested = self
            .base_work_units
            .checked_add(split_work)
            .and_then(|work| work.checked_add(trim_work))
            .and_then(|work| work.checked_add(self.polish_work_units))
            .and_then(|work| work.checked_add(self.sign_work_units))
            .ok_or_else(|| NurbsError::Domain {
                what: "NURBS shell query-work accounting overflows u128".to_string(),
            })?;
        if requested > SDF_QUERY_MAX_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "NURBS shell query requests {requested} work units above defensive ceiling {SDF_QUERY_MAX_WORK_UNITS}"
                ),
            });
        }
        Ok(requested)
    }

    /// One-ULP-outward control-net support, padded outward by one further ULP.
    /// For the exact-real interpretation of the stored f64 homogeneous
    /// controls, correctly rounded division plus this expansion contains each
    /// Cartesian control point and hence the rational surface.
    ///
    /// # Errors
    /// Returns [`NurbsError::Domain`] when `pad` is non-finite or negative.
    pub fn control_aabb(&self, pad: f64) -> Result<Aabb, NurbsError> {
        if !pad.is_finite() || pad < 0.0 {
            return Err(NurbsError::Domain {
                what: "NURBS SDF support padding must be finite and non-negative".to_string(),
            });
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for s in &self.surfaces {
            for row in &s.cpw {
                for h in row {
                    let c = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
                    for k in 0..3 {
                        min[k] = min[k].min(next_down(c[k]));
                        max[k] = max[k].max(next_up(c[k]));
                    }
                }
            }
        }
        let support = Aabb::new(
            Point3::new(
                next_down(min[0] - pad),
                next_down(min[1] - pad),
                next_down(min[2] - pad),
            ),
            Point3::new(
                next_up(max[0] + pad),
                next_up(max[1] + pad),
                next_up(max[2] + pad),
            ),
        );
        if !support.is_finite() {
            return Err(NurbsError::Domain {
                what: "NURBS SDF support overflowed the finite AABB domain after outward padding"
                    .to_string(),
            });
        }
        Ok(support)
    }

    /// Measured unsigned-distance bracket per surface, Gauss–Newton polish,
    /// then trim classification of the winning parameter.
    ///
    /// # Errors
    /// Returns a domain refusal for malformed query settings or excessive
    /// admitted work, and propagates surface/trim evaluation errors.
    pub fn distance(&self, q: [f64; 3], tol: f64, max_splits: u32) -> Result<SdfQuery, NurbsError> {
        self.distance_with_admission(q, tol, max_splits)
            .map(|(query, _)| query)
    }

    fn distance_with_admission(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
    ) -> Result<(SdfQuery, AdmittedNurbsSurface<'_, f64>), NurbsError> {
        validate_distance_request(q, tol)?;
        let _ = self.query_work_units(max_splits)?;
        let surface_count = u64::try_from(self.surfaces.len()).map_err(|_| NurbsError::Domain {
            what: "NURBS shell surface count cannot be represented as u64".to_string(),
        })?;
        let requested_splits = u64::from(max_splits)
            .checked_mul(surface_count)
            .ok_or_else(|| NurbsError::Domain {
                what: "NURBS shell split accounting overflows u64".to_string(),
            })?;
        if requested_splits > u64::from(u32::MAX) {
            return Err(NurbsError::Domain {
                what: "NURBS shell worst-case split accounting exceeds u32".to_string(),
            });
        }
        let mut best: Option<(SdfQuery, AdmittedNurbsSurface<'_, f64>)> = None;
        let mut global_lower = f64::INFINITY;
        let mut total_splits = 0u32;
        for (idx, s) in self.surfaces.iter().enumerate() {
            let admitted = s.admit()?;
            let cd = admitted.closest_point(q, tol, max_splits)?;
            global_lower = global_lower.min(cd.lower);
            total_splits =
                total_splits
                    .checked_add(cd.iterations)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "NURBS shell split accounting exceeds u32".to_string(),
                    })?;
            let polished = polish_upper(admitted, q, cd.param, cd.upper);
            let (upper, param, trim_downgrade) = if let Some(trim) = &self.trims[idx] {
                select_trimmed_witness(trim, (cd.upper, cd.param), polished)?
            } else {
                (polished.0, polished.1, false)
            };
            let cand = SdfQuery {
                lower: cd.lower.min(upper),
                upper,
                param,
                surface: idx,
                trim_downgrade,
                splits: 0,
            };
            best = Some(match best {
                None => (cand, admitted),
                // A point on a kept surface is a usable upper witness and
                // always outranks a numerically closer point that was trimmed
                // away. Within the same trim class, retain the smaller
                // evaluated distance deterministically.
                Some((b, _)) if b.trim_downgrade && !cand.trim_downgrade => (cand, admitted),
                Some((b, b_surface)) if !b.trim_downgrade && cand.trim_downgrade => (b, b_surface),
                Some((b, _)) if cand.upper < b.upper => (cand, admitted),
                Some(existing) => existing,
            });
        }
        let (mut out, admitted) = best.expect("non-empty shell");
        out.lower = global_lower.min(out.upper);
        out.splits = total_splits;
        if out.trim_downgrade {
            out.upper = f64::INFINITY;
        }
        Ok((out, admitted))
    }
}

/// Exact rational endpoints of the 2^-20 cell containing a finite f64.
/// Multiplication by this power of two is exact while representable. Values
/// outside the exactly integral f64 grid-index range fail closed so no
/// saturating float-to-integer cast can manufacture trim authority.
fn trim_parameter_cell(value: f64) -> Option<(Rat, Rat)> {
    let scaled = value * TRIM_SCALE_F64;
    if !scaled.is_finite() || scaled.abs() > MAX_EXACT_GRID_INDEX {
        return None;
    }
    let lo_f = scaled.floor();
    let hi_f = scaled.ceil();
    if lo_f.abs() > MAX_EXACT_GRID_INDEX || hi_f.abs() > MAX_EXACT_GRID_INDEX {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    let lo = i128::from(lo_f as i64);
    #[allow(clippy::cast_possible_truncation)]
    let hi = i128::from(hi_f as i64);
    Some((Rat::new(lo, TRIM_SCALE), Rat::new(hi, TRIM_SCALE)))
}

fn trim_witness_is_kept(trim: &TrimmedPatch, parameter: [f64; 2]) -> Result<bool, NurbsError> {
    match (
        trim_parameter_cell(parameter[0]),
        trim_parameter_cell(parameter[1]),
    ) {
        (Some((umin, umax)), Some((vmin, vmax))) => {
            Ok(trim.classify_box([umin, vmin], [umax, vmax])? == Classification::Inside)
        }
        // A parameter outside the checked rationalization domain cannot
        // acquire trim authority from a saturating cast.
        _ => Ok(false),
    }
}

fn select_trimmed_witness(
    trim: &TrimmedPatch,
    base: (f64, [f64; 2]),
    polished: (f64, [f64; 2]),
) -> Result<(f64, [f64; 2], bool), NurbsError> {
    let base_kept = trim_witness_is_kept(trim, base.1)?;
    let polished_kept = if polished.1.map(f64::to_bits) == base.1.map(f64::to_bits) {
        base_kept
    } else {
        trim_witness_is_kept(trim, polished.1)?
    };
    Ok(match (base_kept, polished_kept) {
        (true, true) if polished.0 < base.0 => (polished.0, polished.1, false),
        (true, _) => (base.0, base.1, false),
        (false, true) => (polished.0, polished.1, false),
        (false, false) => (polished.0, polished.1, true),
    })
}

/// Damped Gauss–Newton on `min ‖S(u,v) − q‖²`: only an ACCEPTED
/// improvement (a strictly smaller evaluated distance inside the domain)
/// updates the answer, so the returned upper estimate comes from a genuinely
/// evaluated surface point (still ordinary, non-directed f64 arithmetic).
fn polish_upper(
    s: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    start: [f64; 2],
    upper0: f64,
) -> (f64, [f64; 2]) {
    let (ulo, uhi) = s.knots_u().domain();
    let (vlo, vhi) = s.knots_v().domain();
    let mut best = (upper0, start);
    let mut uv = start;
    for _ in 0..POLISH_STEPS {
        let Ok((pos, du, dv)) = s.partials(uv[0], uv[1]) else {
            break;
        };
        let r = [pos[0] - q[0], pos[1] - q[1], pos[2] - q[2]];
        // Normal equations of the 2x3 Jacobian.
        let (a, b, c) = (dot3(du, du), dot3(du, dv), dot3(dv, dv));
        let (g0, g1) = (dot3(du, r), dot3(dv, r));
        let detm = a * c - b * b;
        if detm.abs() < 1e-300 {
            break;
        }
        let step = [-(c * g0 - b * g1) / detm, -(a * g1 - b * g0) / detm];
        let mut damp = 1.0f64;
        let mut improved = false;
        for _ in 0..4 {
            let cand = [
                (uv[0] + damp * step[0]).clamp(ulo, uhi),
                (uv[1] + damp * step[1]).clamp(vlo, vhi),
            ];
            if let Ok(p) = s.eval(cand[0], cand[1]) {
                let d = norm3([p[0] - q[0], p[1] - q[1], p[2] - q[2]]);
                if d < best.0 {
                    best = (d, cand);
                    uv = cand;
                    improved = true;
                    break;
                }
            }
            damp *= 0.25;
        }
        if !improved {
            break;
        }
    }
    best
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The Chart presentation of a [`ShellSdf`] (the router-visible form).
#[derive(Debug)]
pub struct ShellSdfChart {
    shell: ShellSdf,
    tol: f64,
    max_splits: u32,
    support_pad: f64,
}

impl ShellSdfChart {
    /// Wrap a shell with query effort settings.
    ///
    /// # Errors
    /// Returns a structured domain error for non-finite/negative settings or a
    /// configured query whose statically bounded work exceeds the legacy-path
    /// ceiling.
    pub fn new(
        shell: ShellSdf,
        tol: f64,
        max_splits: u32,
        support_pad: f64,
    ) -> Result<ShellSdfChart, NurbsError> {
        if !tol.is_finite() || tol < 0.0 {
            return Err(NurbsError::Domain {
                what: "NURBS SDF tolerance must be finite and non-negative".to_string(),
            });
        }
        if !support_pad.is_finite() || support_pad < 0.0 {
            return Err(NurbsError::Domain {
                what: "NURBS SDF support padding must be finite and non-negative".to_string(),
            });
        }
        let _ = shell.query_work_units(max_splits)?;
        let _ = shell.control_aabb(support_pad)?;
        Ok(ShellSdfChart {
            shell,
            tol,
            max_splits,
            support_pad,
        })
    }

    /// Signed-or-unsigned distance estimate for one point.
    fn sample(&self, x: Point3) -> Result<ChartSample, NurbsError> {
        let q = [x.x, x.y, x.z];
        let (query, admitted_surface) =
            self.shell
                .distance_with_admission(q, self.tol, self.max_splits)?;
        if query.trim_downgrade {
            // The found point is trimmed away, so it is not a witness on the
            // kept surface. Do not retain its finite value as the chart's
            // nominal distance: consumers that forget to inspect `error` must
            // fail closed too.
            return Ok(ChartSample {
                signed_distance: f64::INFINITY,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            });
        }
        let (mut lo, mut hi) = (query.lower, query.upper);
        let Some((sign, gradient)) = self.sign_and_gradient(q, &query, admitted_surface) else {
            // A declared orientation cannot assign sign through a singular or
            // contradictory local chart. Preserve no nominal finite value that
            // a consumer could misread while ignoring the certificate.
            return Ok(ChartSample {
                signed_distance: f64::INFINITY,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            });
        };
        let signed = sign * query.upper;
        if sign < 0.0 {
            (lo, hi) = (-hi, -lo);
        }
        Ok(ChartSample {
            signed_distance: signed,
            gradient,
            lipschitz: None,
            error: NumericalCertificate::estimate(lo, hi),
        })
    }

    /// Sign from declared orientation (unit normal · unit offset); gradient from
    /// the offset direction when it is well-defined. A singular closest
    /// parameter triggers a deterministic nearby regular-witness search on an
    /// untrimmed surface. If no consistent regular witness stays within the
    /// measured upper-distance tolerance, oriented sign fails closed.
    fn sign_and_gradient(
        &self,
        q: [f64; 3],
        query: &SdfQuery,
        surface: AdmittedNurbsSurface<'_, f64>,
    ) -> Option<(f64, Option<Vec3>)> {
        match self.shell.orientation {
            Orientation::Unknown => Some((1.0, unsigned_gradient(surface, q, query.param))),
            Orientation::Outward => oriented_sign_at(surface, q, query.param)
                .or_else(|| self.regular_witness_sign(q, query, surface)),
        }
    }

    fn regular_witness_sign(
        &self,
        q: [f64; 3],
        query: &SdfQuery,
        surface: AdmittedNurbsSurface<'_, f64>,
    ) -> Option<(f64, Option<Vec3>)> {
        // A nearby parameter could cross a trim curve. Until the trim path can
        // certify a connected kept neighborhood, do not use it for sign repair.
        if self.shell.trims[query.surface].is_some() {
            return None;
        }
        let (ulo, uhi) = surface.knots_u().domain();
        let (vlo, vhi) = surface.knots_v().domain();
        let spans = [uhi - ulo, vhi - vlo];
        if spans.iter().any(|span| !span.is_finite() || *span <= 0.0) {
            return None;
        }
        let slack = self
            .tol
            .max(256.0 * f64::EPSILON * query.upper.abs().max(1.0));
        let admitted_upper = query.upper + slack;
        if !admitted_upper.is_finite() {
            return None;
        }
        let mut best: Option<(f64, [f64; 2], f64, Option<Vec3>)> = None;
        let mut observed_sign = None;
        for exponent in [8, 12, 16, 20, 24, 28, 32, 36, 40, 44] {
            let scale = fs_math::det::powi(2.0, -exponent);
            let du = spans[0] * scale;
            let dv = spans[1] * scale;
            for candidate in [
                [(query.param[0] + du).min(uhi), query.param[1]],
                [(query.param[0] - du).max(ulo), query.param[1]],
                [query.param[0], (query.param[1] + dv).min(vhi)],
                [query.param[0], (query.param[1] - dv).max(vlo)],
                [
                    (query.param[0] + du).min(uhi),
                    (query.param[1] + dv).min(vhi),
                ],
                [
                    (query.param[0] + du).min(uhi),
                    (query.param[1] - dv).max(vlo),
                ],
                [
                    (query.param[0] - du).max(ulo),
                    (query.param[1] + dv).min(vhi),
                ],
                [
                    (query.param[0] - du).max(ulo),
                    (query.param[1] - dv).max(vlo),
                ],
            ] {
                if candidate == query.param {
                    continue;
                }
                let Ok(position) = surface.eval(candidate[0], candidate[1]) else {
                    continue;
                };
                let distance = norm3([position[0] - q[0], position[1] - q[1], position[2] - q[2]]);
                if !distance.is_finite() || distance > admitted_upper {
                    continue;
                }
                let Some((sign, gradient)) = oriented_sign_at(surface, q, candidate) else {
                    continue;
                };
                if observed_sign.is_some_and(|seen| seen != sign) {
                    return None;
                }
                observed_sign = Some(sign);
                let replace = best.as_ref().is_none_or(|current| {
                    distance.total_cmp(&current.0).is_lt()
                        || (distance.to_bits() == current.0.to_bits()
                            && candidate
                                .into_iter()
                                .map(f64::to_bits)
                                .cmp(current.1.into_iter().map(f64::to_bits))
                                .is_lt())
                });
                if replace {
                    best = Some((distance, candidate, sign, gradient));
                }
            }
        }
        best.map(|(_, _, sign, gradient)| (sign, gradient))
    }
}

fn unsigned_gradient(
    surface: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    param: [f64; 2],
) -> Option<Vec3> {
    let (position, _, _) = surface.partials(param[0], param[1]).ok()?;
    let offset = [q[0] - position[0], q[1] - position[1], q[2] - position[2]];
    let norm = norm3(offset);
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    Some(Vec3::new(
        offset[0] / norm,
        offset[1] / norm,
        offset[2] / norm,
    ))
}

fn oriented_sign_at(
    surface: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    param: [f64; 2],
) -> Option<(f64, Option<Vec3>)> {
    let (position, du, dv) = surface.partials(param[0], param[1]).ok()?;
    if position
        .iter()
        .chain(du.iter())
        .chain(dv.iter())
        .any(|value| !value.is_finite())
    {
        return None;
    }
    let du_norm = norm3(du);
    let dv_norm = norm3(dv);
    if !du_norm.is_finite() || !dv_norm.is_finite() || du_norm == 0.0 || dv_norm == 0.0 {
        return None;
    }
    let unit_u = du.map(|value| value / du_norm);
    let unit_v = dv.map(|value| value / dv_norm);
    let normal = [
        unit_u[1] * unit_v[2] - unit_u[2] * unit_v[1],
        unit_u[2] * unit_v[0] - unit_u[0] * unit_v[2],
        unit_u[0] * unit_v[1] - unit_u[1] * unit_v[0],
    ];
    let normal_norm = norm3(normal);
    if !normal_norm.is_finite() || normal_norm <= 256.0 * f64::EPSILON {
        return None;
    }
    let unit_normal = normal.map(|value| value / normal_norm);
    let offset = [q[0] - position[0], q[1] - position[1], q[2] - position[2]];
    let offset_norm = norm3(offset);
    if !offset_norm.is_finite() {
        return None;
    }
    if offset_norm == 0.0 {
        return Some((
            1.0,
            Some(Vec3::new(unit_normal[0], unit_normal[1], unit_normal[2])),
        ));
    }
    let unit_offset = offset.map(|value| value / offset_norm);
    let alignment = dot3(unit_normal, unit_offset);
    if !alignment.is_finite() || alignment.abs() <= 256.0 * f64::EPSILON {
        return None;
    }
    let sign = if alignment < 0.0 { -1.0 } else { 1.0 };
    let gradient = Vec3::new(unit_offset[0], unit_offset[1], unit_offset[2]);
    Some((
        sign,
        Some(if sign < 0.0 {
            gradient.scale(-1.0)
        } else {
            gradient
        }),
    ))
}

impl Chart for ShellSdfChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        self.sample(x).unwrap_or(ChartSample {
            signed_distance: f64::INFINITY,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        })
    }

    fn support(&self) -> Aabb {
        self.shell
            .control_aabb(self.support_pad)
            .expect("ShellSdfChart constructor validates its private support padding")
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    fn name(&self) -> &'static str {
        match self.shell.orientation {
            Orientation::Outward => "nurbs-sdf/estimated-signed",
            Orientation::Unknown => "nurbs-sdf/estimated-unsigned",
        }
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::Unknown
    }
}

#[cfg(test)]
mod tests {
    use super::{select_trimmed_witness, store_f32_with_error};
    use crate::basis::KnotVector;
    use crate::curve::NurbsCurve;
    use crate::rat::Rat;
    use crate::trim::{TrimLoop, TrimmedPatch};

    fn unit_square_trim() -> TrimmedPatch {
        let knots = KnotVector::new(
            vec![
                Rat::int(0),
                Rat::int(0),
                Rat::new(1, 4),
                Rat::new(1, 2),
                Rat::new(3, 4),
                Rat::int(1),
                Rat::int(1),
            ],
            1,
        )
        .expect("square knots");
        let curve = NurbsCurve::new(
            knots,
            &[
                [Rat::int(0), Rat::int(0)],
                [Rat::int(1), Rat::int(0)],
                [Rat::int(1), Rat::int(1)],
                [Rat::int(0), Rat::int(1)],
                [Rat::int(0), Rat::int(0)],
            ],
            &[Rat::int(1); 5],
        )
        .expect("square curve");
        TrimmedPatch::new(vec![TrimLoop::new(curve).expect("closed square")])
    }

    #[test]
    fn optional_polish_cannot_erase_a_kept_base_witness() {
        let trim = unit_square_trim();
        let selected = select_trimmed_witness(&trim, (2.0, [0.5, 0.5]), (1.0, [2.0, 2.0]))
            .expect("classification");
        assert_eq!(selected, (2.0, [0.5, 0.5], false));

        let repaired = select_trimmed_witness(&trim, (2.0, [2.0, 2.0]), (1.0, [0.5, 0.5]))
            .expect("classification");
        assert_eq!(repaired, (1.0, [0.5, 0.5], false));

        let absent = select_trimmed_witness(&trim, (2.0, [2.0, 2.0]), (1.0, [3.0, 3.0]))
            .expect("classification");
        assert_eq!(absent, (1.0, [3.0, 3.0], true));
    }

    #[test]
    fn f32_storage_error_covers_rounding_and_signed_underflow() {
        let ordinary = 1.0 + f64::EPSILON;
        let (stored, error) = store_f32_with_error(ordinary);
        assert_eq!(stored, 1.0);
        assert!(error >= (ordinary - f64::from(stored)).abs());

        let tiny_negative = -f64::from(f32::from_bits(1)) / 4.0;
        let (stored_zero, underflow_error) = store_f32_with_error(tiny_negative);
        assert_eq!(stored_zero.to_bits(), (-0.0f32).to_bits());
        assert!(underflow_error >= tiny_negative.abs());
        assert!(underflow_error.is_finite());
    }
}

/// One generated tile of measured samples on a regular grid.
#[derive(Debug, Clone)]
pub struct SdfTile {
    /// Grid resolution per axis.
    pub n: usize,
    /// Field values (x-fastest), signed per the shell orientation; a
    /// downgraded cell stores positive infinity as an unusable sentinel.
    /// Finite values use ordinary IEEE-754 f64→f32 rounding; a nonzero value
    /// may therefore become a signed zero, with that loss reported separately.
    pub values: Vec<f32>,
    /// Worst pre-storage measured bracket width among near-surface cells.
    pub worst_near_width: f64,
    /// Worst pre-storage width among far cells (cheap bounds are allowed there).
    pub worst_far_width: f64,
    /// Worst outward-expanded absolute f64→f32 storage error among usable
    /// near-surface cells; infinity if any such cell is unusable.
    pub worst_near_storage_error: f64,
    /// Worst outward-expanded absolute f64→f32 storage error among usable far
    /// cells; infinity if any such cell is unusable.
    pub worst_far_storage_error: f64,
    /// Total branch-and-bound splits (the throughput ledger line).
    pub total_splits: u64,
    /// Total unusable cells represented by the positive-infinity sentinel.
    pub downgraded: usize,
    /// Unusable because trim classification could not retain a kept-surface
    /// witness.
    pub trim_downgraded: usize,
    /// Unusable because declared orientation could not establish a consistent
    /// sign at a regular witness.
    pub sign_downgraded: usize,
    /// Unusable because a finite f64 value is outside finite f32 range.
    pub range_downgraded: usize,
}

fn store_f32_with_error(value: f64) -> (f32, f64) {
    #[allow(clippy::cast_possible_truncation)]
    let stored = value as f32;
    let error = (value - f64::from(stored)).abs();
    let outward_error = if error == 0.0 { 0.0 } else { next_up(error) };
    (stored, outward_error)
}

/// Tiled generation with adaptive effort under defensive static ceilings:
/// a tighter requested tolerance inside the near band (|d| ≤ 2 cell
/// diagonals), and cheaper measured bounds elsewhere.
///
/// # Errors
/// Returns a structured domain/structure error for malformed tile settings,
/// overflow, requests above the defensive static ceilings, allocation refusal,
/// or a structural surface error.
pub fn generate_tile(
    chart: &ShellSdfChart,
    aabb: &Aabb,
    n: usize,
    tol_near: f64,
    max_splits: u32,
) -> Result<SdfTile, NurbsError> {
    if n < 2 {
        return Err(NurbsError::Domain {
            what: "an SDF tile needs at least 2 samples per axis".to_string(),
        });
    }
    if !aabb.is_finite() {
        return Err(NurbsError::Domain {
            what: "an SDF tile requires a finite well-formed AABB".to_string(),
        });
    }
    if !tol_near.is_finite() || tol_near < 0.0 {
        return Err(NurbsError::Domain {
            what: "SDF tile tolerance must be finite and non-negative".to_string(),
        });
    }
    if max_splits > CLOSEST_MAX_SPLITS {
        return Err(NurbsError::Domain {
            what: "SDF tile split request exceeds the defensive legacy-path ceiling".to_string(),
        });
    }
    let sample_count = n
        .checked_mul(n)
        .and_then(|square| square.checked_mul(n))
        .ok_or_else(|| NurbsError::Domain {
            what: "SDF tile sample count overflows usize".to_string(),
        })?;
    if sample_count > SDF_TILE_MAX_SAMPLES {
        return Err(NurbsError::Domain {
            what: format!(
                "SDF tile sample count {sample_count} exceeds defensive ceiling {SDF_TILE_MAX_SAMPLES}"
            ),
        });
    }
    let sample_count_u128 = u128::try_from(sample_count).map_err(|_| NurbsError::Domain {
        what: "SDF tile sample count cannot be represented as u128".to_string(),
    })?;
    let surface_count =
        u128::try_from(chart.shell.surfaces.len()).map_err(|_| NurbsError::Domain {
            what: "SDF shell surface count cannot be represented as u128".to_string(),
        })?;
    let splits_per_query = u128::from(max_splits)
        .checked_mul(surface_count)
        .ok_or_else(|| NurbsError::Domain {
            what: "SDF shell per-query split count overflows u128".to_string(),
        })?;
    if splits_per_query > u128::from(u32::MAX) {
        return Err(NurbsError::Domain {
            what: "SDF shell per-query split accounting exceeds u32".to_string(),
        });
    }
    let coarse_work = chart.shell.query_work_units(max_splits / 4)?;
    let fine_work = chart.shell.query_work_units(max_splits)?;
    let worst_case_work = sample_count_u128
        .checked_mul(
            coarse_work
                .checked_add(fine_work)
                .ok_or_else(|| NurbsError::Domain {
                    what: "SDF tile per-sample work accounting overflows u128".to_string(),
                })?,
        )
        .ok_or_else(|| NurbsError::Domain {
            what: "SDF tile worst-case work count overflows u128".to_string(),
        })?;
    if worst_case_work > SDF_TILE_MAX_WORST_CASE_WORK {
        return Err(NurbsError::Domain {
            what: format!(
                "SDF tile worst-case work request {worst_case_work} exceeds defensive ceiling {SDF_TILE_MAX_WORST_CASE_WORK}"
            ),
        });
    }
    #[allow(clippy::cast_precision_loss)]
    let step = [
        (aabb.max.x - aabb.min.x) / (n - 1) as f64,
        (aabb.max.y - aabb.min.y) / (n - 1) as f64,
        (aabb.max.z - aabb.min.z) / (n - 1) as f64,
    ];
    let diag = norm3(step);
    if step.iter().any(|value| !value.is_finite()) || !diag.is_finite() {
        return Err(NurbsError::Domain {
            what: "SDF tile step or diagonal is not representable as finite f64".to_string(),
        });
    }
    // Refinement fires within two diagonals of the surface; the tighter
    // measured-effort lane is used only for cells adjacent by this heuristic
    // (medial-axis
    // cells are equidistant from many patches — hull pruning stalls
    // there, and the estimate honestly stays wider).
    let refine_band = 2.0 * diag;
    let near_band = 0.75 * diag;
    let tol_far = (refine_band * 0.5).max(tol_near);
    if !refine_band.is_finite() || !near_band.is_finite() || !tol_far.is_finite() {
        return Err(NurbsError::Domain {
            what: "SDF tile refinement scale is not representable as finite f64".to_string(),
        });
    }
    let mut values = Vec::new();
    values
        .try_reserve_exact(sample_count)
        .map_err(|_| NurbsError::Structure {
            what: format!("SDF tile allocation refused for {sample_count} samples"),
        })?;
    let mut tile = SdfTile {
        n,
        values,
        worst_near_width: 0.0,
        worst_far_width: 0.0,
        worst_near_storage_error: 0.0,
        worst_far_storage_error: 0.0,
        total_splits: 0,
        downgraded: 0,
        trim_downgraded: 0,
        sign_downgraded: 0,
        range_downgraded: 0,
    };
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                #[allow(clippy::cast_precision_loss)]
                let q = [
                    aabb.min.x + i as f64 * step[0],
                    aabb.min.y + j as f64 * step[1],
                    aabb.min.z + k as f64 * step[2],
                ];
                // Cheap pass first; refine only inside the near band.
                let coarse = chart
                    .shell
                    .distance_with_admission(q, tol_far, max_splits / 4)?;
                let (query, admitted_surface) = if coarse.0.lower <= refine_band {
                    let fine = chart
                        .shell
                        .distance_with_admission(q, tol_near, max_splits)?;
                    tile.total_splits += u64::from(coarse.0.splits);
                    fine
                } else {
                    coarse
                };
                tile.total_splits += u64::from(query.splits);
                let width = query.upper - query.lower;
                let near_cell = query.upper <= near_band;
                if near_cell {
                    tile.worst_near_width = tile.worst_near_width.max(width);
                } else {
                    tile.worst_far_width = tile.worst_far_width.max(width);
                }
                if query.trim_downgrade {
                    tile.downgraded += 1;
                    tile.trim_downgraded += 1;
                    if near_cell {
                        tile.worst_near_width = f64::INFINITY;
                        tile.worst_near_storage_error = f64::INFINITY;
                    } else {
                        tile.worst_far_width = f64::INFINITY;
                        tile.worst_far_storage_error = f64::INFINITY;
                    }
                    tile.values.push(f32::INFINITY);
                    continue;
                }
                let Some((sign, _)) = chart.sign_and_gradient(q, &query, admitted_surface) else {
                    tile.downgraded += 1;
                    tile.sign_downgraded += 1;
                    if near_cell {
                        tile.worst_near_width = f64::INFINITY;
                        tile.worst_near_storage_error = f64::INFINITY;
                    } else {
                        tile.worst_far_width = f64::INFINITY;
                        tile.worst_far_storage_error = f64::INFINITY;
                    }
                    tile.values.push(f32::INFINITY);
                    continue;
                };
                let value = sign * query.upper;
                if !value.is_finite() || value.abs() > f64::from(f32::MAX) {
                    tile.downgraded += 1;
                    tile.range_downgraded += 1;
                    if near_cell {
                        tile.worst_near_width = f64::INFINITY;
                        tile.worst_near_storage_error = f64::INFINITY;
                    } else {
                        tile.worst_far_width = f64::INFINITY;
                        tile.worst_far_storage_error = f64::INFINITY;
                    }
                    tile.values.push(f32::INFINITY);
                    continue;
                }
                let (stored, storage_error) = store_f32_with_error(value);
                if near_cell {
                    tile.worst_near_storage_error =
                        tile.worst_near_storage_error.max(storage_error);
                } else {
                    tile.worst_far_storage_error = tile.worst_far_storage_error.max(storage_error);
                }
                tile.values.push(stored);
            }
        }
    }
    Ok(tile)
}
