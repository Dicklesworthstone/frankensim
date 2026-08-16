//! Mesh -> A(x) bore extraction: centerline + per-station cross-section
//! areas for TUBE lumens (music program bead
//! `frankensim-music-v8-root-3ez8g.3.1`; ingest law).
//!
//! Pipeline: [`crate::medial_poles`] pole cloud -> deterministic thinning ->
//! Euclidean MST -> longest path (tube spine) with an explicit BRANCH
//! REFUSAL -> smoothing + arc-length resampling -> per-station polar
//! plane-sections (bisected boundary rays) -> cross-checks against the
//! medial-pole radii and [`crate::thickness_at`] -> axial volume vs the
//! CERTIFIED whole-lumen [`crate::geometric_moments`] enclosure.
//!
//! AUTHORITY (binding): the centerline, areas, and cross-checks are
//! `Estimate` BY CONSTRUCTION — the underlying queries carry no
//! no-tunneling theorem and this module adds none. Only the volume-closure
//! enclosure is rigorous, and only when the chart itself publishes
//! `TraceStepClaim::ExactDistance`; every weaker chart yields
//! [`VolumeClosure::Unavailable`] with the claim recorded, never a
//! downgraded number. The receipt is one sentence of honesty: "Estimate,
//! cross-checked by medial-pole radii and thickness, closed by a certified
//! volume enclosure (when the chart can certify one)."
//!
//! v1 topology contract: single unbranched tubes, open or closed-loop
//! (a full torus chains into a loop; the MST drops exactly one loop edge
//! and endpoint proximity restores it). Branched lumens (valve clusters)
//! REFUSE with [`BoreError::BranchedLumen`] — segment-wise extraction
//! between authored cut planes is the caller's workflow, never a guess
//! here.

use fs_evidence::NumericalKind;
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, Point3, TraceStepClaim, Vec3};
use fs_rep_mesh::Soup;

use crate::{QueryError, geometric_moments, thickness_at};

/// Configuration for one extraction. Every knob is an authored, receipt-
/// visible choice; none is derived silently.
#[derive(Debug, Clone)]
pub struct BoreConfig {
    /// Medial-pole sliver filter multiplier (see [`crate::medial_poles`]).
    pub lambda: f64,
    /// Number of stations along the resampled centerline.
    pub stations: usize,
    /// Polar boundary rays per station (area quadrature resolution).
    pub rays_per_station: usize,
    /// Bisection tolerance for boundary rays and end extension [m].
    pub root_tolerance_m: f64,
    /// Deterministic pole-cloud thinning cap (MST is O(n^2)).
    pub max_poles: usize,
    /// An off-spine MST subtree REACHING farther from the spine than this
    /// multiple of the local pole radius is a BRANCH -> refusal (medial
    /// noise hugs the tube; a genuine side arm extends away from it).
    pub branch_reach_factor: f64,
    /// Centerline moving-average half-width (0 = no smoothing).
    pub smoothing_half_width: usize,
    /// Boundary rays give up (refuse) beyond this multiple of the local
    /// pole radius.
    pub ray_cap_factor: f64,
    /// Grid spacing for the certified volume closure; `None` skips the
    /// closure entirely (recorded as skipped, not silently absent).
    pub closure_h_m: Option<f64>,
    /// Run the per-station thickness cross-check (one `thickness_at` per
    /// station; local failures are counted, not fatal).
    pub thickness_cross_check: bool,
}

impl Default for BoreConfig {
    fn default() -> Self {
        BoreConfig {
            lambda: 1.2,
            stations: 33,
            rays_per_station: 64,
            root_tolerance_m: 1.0e-6,
            max_poles: 1500,
            branch_reach_factor: 2.0,
            smoothing_half_width: 2,
            ray_cap_factor: 4.0,
            closure_h_m: None,
            thickness_cross_check: true,
        }
    }
}

/// One station of the extracted bore.
#[derive(Debug, Clone)]
pub struct BoreStation {
    /// Arc length from the first station [m].
    pub arc_length_m: f64,
    /// Centerline point.
    pub center_m: Point3,
    /// Unit tangent (direction of increasing arc length).
    pub tangent: Vec3,
    /// Cross-section area from the polar plane-section [m^2].
    pub area_m2: f64,
    /// `sqrt(area/pi)` — the equivalent circular radius [m].
    pub equivalent_radius_m: f64,
    /// Interpolated medial-pole radius at this station [m].
    pub pole_radius_m: f64,
    /// `|area - pi*pole_radius^2| / area` — the medial cross-check.
    pub pole_area_deviation: f64,
    /// Thickness at the station's first boundary point (None = skipped).
    pub thickness_m: Option<f64>,
    /// `|pi*(t/2)^2 - area| / area` (None = thickness skipped).
    pub thickness_area_deviation: Option<f64>,
    /// Shortest boundary ray [m] (eccentricity diagnostic).
    pub min_ray_m: f64,
    /// Longest boundary ray [m].
    pub max_ray_m: f64,
}

/// Outcome of the whole-lumen volume-closure comparison.
#[derive(Debug, Clone, PartialEq)]
pub enum VolumeClosure {
    /// The chart certified a rigorous volume enclosure.
    Certified {
        /// Enclosure lower bound [m^3].
        lo_m3: f64,
        /// Enclosure upper bound [m^3].
        hi_m3: f64,
        /// Grid spacing used [m].
        h_m: f64,
        /// `|axial - mid| / mid` where mid is the enclosure midpoint.
        axial_vs_mid_rel: f64,
    },
    /// The chart's trace claim cannot support a certificate — recorded,
    /// never downgraded into a number.
    Unavailable {
        /// The claim the chart actually published.
        claim: TraceStepClaim,
    },
    /// The caller skipped the closure (`closure_h_m: None`).
    Skipped,
}

/// The A(x) receipt: samples, centerline, cross-check deviations, closure
/// result, source identity, and the authority label.
#[derive(Debug, Clone)]
pub struct BoreExtraction {
    /// Ordered stations along the centerline.
    pub stations: Vec<BoreStation>,
    /// The chained tube closes on itself (torus-like lumen).
    pub closed_loop: bool,
    /// Total centerline arc length [m] (includes the closing segment for
    /// loops).
    pub total_length_m: f64,
    /// Trapezoidal `integral A ds` over the stations [m^3].
    pub axial_volume_m3: f64,
    /// The certified-closure comparison.
    pub volume_closure: VolumeClosure,
    /// Caller-supplied immutable source identity for provenance.
    pub source_geometry_id: String,
    /// FNV-1a digest of the boundary soup bytes (positions + triangles).
    pub boundary_digest: u64,
    /// Stations whose thickness cross-check locally failed (skipped).
    pub thickness_skipped: u32,
    /// The extraction's authority: always `Estimate` (see module doc).
    pub authority: NumericalKind,
}

impl BoreExtraction {
    /// Largest medial cross-check deviation across all stations.
    #[must_use]
    pub fn worst_pole_deviation(&self) -> f64 {
        self.stations
            .iter()
            .map(|s| s.pole_area_deviation)
            .fold(0.0, f64::max)
    }

    /// Whether the certified closure is within `band` of the enclosure
    /// (relative to the midpoint, band widened by the enclosure's own
    /// half-width). `None` when no certificate exists.
    #[must_use]
    pub fn closure_within(&self, band: f64) -> Option<bool> {
        match &self.volume_closure {
            VolumeClosure::Certified {
                lo_m3,
                hi_m3,
                axial_vs_mid_rel,
                ..
            } => {
                let mid = 0.5 * (lo_m3 + hi_m3);
                let half_width_rel = if mid > 0.0 {
                    0.5 * (hi_m3 - lo_m3) / mid
                } else {
                    f64::INFINITY
                };
                Some(*axial_vs_mid_rel <= band + half_width_rel)
            }
            _ => None,
        }
    }

    /// Per-station JSON-lines debug records — enough to diagnose a bad
    /// extraction from logs alone (s, A, pole radius, deviations, rays).
    #[must_use]
    pub fn debug_lines(&self) -> Vec<String> {
        self.stations
            .iter()
            .enumerate()
            .map(|(i, s)| {
                format!(
                    "{{\"station\":{i},\"s_m\":{:.6e},\"area_m2\":{:.6e},\
                     \"eq_radius_m\":{:.6e},\"pole_radius_m\":{:.6e},\
                     \"pole_dev\":{:.3e},\"thickness_dev\":{},\
                     \"min_ray_m\":{:.6e},\"max_ray_m\":{:.6e}}}",
                    s.arc_length_m,
                    s.area_m2,
                    s.equivalent_radius_m,
                    s.pole_radius_m,
                    s.pole_area_deviation,
                    s.thickness_area_deviation
                        .map_or_else(|| "null".to_string(), |d| format!("{d:.3e}")),
                    s.min_ray_m,
                    s.max_ray_m,
                )
            })
            .collect()
    }
}

/// Typed refusals from bore extraction.
#[derive(Debug, Clone, PartialEq)]
pub enum BoreError {
    /// An underlying query refused.
    Query(QueryError),
    /// A configuration knob is unusable as stated.
    InvalidConfig {
        /// Which knob and why.
        what: &'static str,
    },
    /// The filtered pole cloud is too small to chain.
    TooFewPoles {
        /// Poles surviving the filter and thinning.
        found: usize,
        /// Minimum this pipeline needs.
        needed: usize,
    },
    /// The lumen branches: an off-spine subtree is longer than the branch
    /// tolerance. Segment-wise extraction between authored cut planes is
    /// the supported workflow.
    BranchedLumen {
        /// Spine point where the branch attaches.
        at: [f64; 3],
        /// Arc length of the offending subtree [m].
        side_length_m: f64,
    },
    /// The pole chain is degenerate (zero-length spine, coincident poles).
    DegeneratePoleChain {
        /// Diagnosis.
        reason: &'static str,
    },
    /// A resampled centerline point is not inside the lumen — the chain
    /// left the region and no area can be measured there.
    CenterOutsideLumen {
        /// Station index.
        station: usize,
        /// Arc length of the bad station [m].
        arc_length_m: f64,
    },
    /// A boundary ray never exited the region within the ray cap.
    RayNoExit {
        /// Station index.
        station: usize,
        /// Ray index.
        ray: usize,
        /// Station center.
        center: [f64; 3],
        /// Station unit tangent.
        tangent: [f64; 3],
        /// Ray cap that was exhausted [m].
        cap_m: f64,
        /// Local pole radius the cap was scaled from [m].
        pole_radius_m: f64,
    },
}

impl core::fmt::Display for BoreError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BoreError::Query(e) => write!(f, "bore extraction query failure: {e:?}"),
            BoreError::InvalidConfig { what } => write!(f, "bad bore config: {what}"),
            BoreError::TooFewPoles { found, needed } => {
                write!(f, "too few medial poles to chain: {found} < {needed}")
            }
            BoreError::BranchedLumen { at, side_length_m } => write!(
                f,
                "branched lumen at ({}, {}, {}): side arm {side_length_m} m; \
                 extract segment-wise between authored cut planes",
                at[0], at[1], at[2]
            ),
            BoreError::DegeneratePoleChain { reason } => {
                write!(f, "degenerate pole chain: {reason}")
            }
            BoreError::CenterOutsideLumen {
                station,
                arc_length_m,
            } => write!(
                f,
                "centerline station {station} (s = {arc_length_m} m) is outside the lumen"
            ),
            BoreError::RayNoExit {
                station,
                ray,
                cap_m,
                ..
            } => write!(
                f,
                "boundary ray {ray} at station {station} never exited the region within {cap_m} m"
            ),
        }
    }
}

impl core::error::Error for BoreError {}

impl From<QueryError> for BoreError {
    fn from(e: QueryError) -> Self {
        BoreError::Query(e)
    }
}

/// Extract the bore of a tube lumen: centerline, A(s) samples, cross-checks,
/// and the volume-closure comparison. See the module doc for the authority
/// statement and topology contract.
///
/// # Errors
/// [`BoreError`] on refusal — branched lumens, degenerate chains, stations
/// outside the region, unusable configuration, or any underlying
/// [`QueryError`].
pub fn extract_bore(
    chart: &dyn Chart,
    boundary: &Soup,
    config: &BoreConfig,
    source_geometry_id: &str,
    cx: &Cx<'_>,
) -> Result<BoreExtraction, BoreError> {
    validate_config(config)?;
    let digest = soup_digest(boundary);

    // 1. Medial pole cloud (Estimate; the existing query owns the sliver
    //    filter), thinned, then reduced to SPINAL poles: keep a pole only
    //    when its radius is near-maximal among the poles within its own
    //    radius. Sheet poles (the disc-shaped medial set near a flat tube
    //    end) always sit next to a strictly larger axis pole, so this
    //    local-maximality filter deletes sheets and keeps the 1-D spine.
    let poles = crate::medial_poles(chart, boundary, config.lambda, cx)?;
    let poles = thin_poles(poles, config.max_poles);
    let poles = spinal_filter(&poles, cx)?;
    if poles.len() < 4 {
        return Err(BoreError::TooFewPoles {
            found: poles.len(),
            needed: 4,
        });
    }

    // 2. MST + longest path = the tube spine; branches refuse.
    let edges = mst_edges(&poles, cx)?;
    let spine = longest_path_poles(&poles, &edges);
    if spine.len() < 3 {
        return Err(BoreError::DegeneratePoleChain {
            reason: "spine shorter than three poles",
        });
    }
    audit_branches(&poles, &edges, &spine, config.branch_reach_factor)?;

    // 3. Loop detection: MST breaks a closed loop at exactly one edge, so
    //    a loop's spine endpoints sit within hop distance of each other.
    let first = poles[spine[0]];
    let last = poles[spine[spine.len() - 1]];
    let end_gap = last.0.delta_from(first.0).norm();
    let spine_len = polyline_length(&poles, &spine);
    if spine_len <= 0.0 {
        return Err(BoreError::DegeneratePoleChain {
            reason: "zero-length spine",
        });
    }
    let closed_loop = end_gap < 1.5 * (first.1 + last.1) && spine_len > 4.0 * end_gap;

    // 4. Snap the centerline to the CLOUD, not the spine: the medial cloud
    //    is a 3-D thicket and the MST diameter path threads its extremes
    //    (its endpoints are off-axis rim poles). Every pole is parameterized
    //    by the arc length of its nearest spine node, binned along the
    //    spine, and each bin's radius^2-weighted centroid becomes a chain
    //    point — the cloud's own center of medial mass at that station.
    let mut chain = bin_centroids(&poles, &spine, config.stations, closed_loop, cx)?;
    chain = smooth_chain(&chain, config.smoothing_half_width, closed_loop);
    let (start_slab_m, end_slab_m) = if closed_loop {
        (0.0, 0.0)
    } else {
        extend_to_boundary(chart, &mut chain, config, cx)?
    };

    // 5. Resample by arc length into stations.
    let sampled = resample(&chain, config.stations, closed_loop);

    // 6. Per-station plane sections + cross-checks.
    let mut stations = Vec::with_capacity(sampled.len());
    let mut thickness_skipped = 0u32;
    let mut extra_start_slab_m = 0.0f64;
    let mut extra_end_slab_m = 0.0f64;
    for (index, (s_len, center, tangent, pole_r)) in sampled.iter().enumerate() {
        crate::query_checkpoint(cx)?;
        let inside = chart.eval(*center, cx).signed_distance < 0.0;
        if !inside {
            return Err(BoreError::CenterOutsideLumen {
                station: index,
                arc_length_m: *s_len,
            });
        }
        let (area, min_ray, max_ray, boundary_points, retreat_m) = {
            // An END station whose shortest ray is far below the local
            // radius has clipped its section plane against the tube's cut
            // face (any tangent tilt does this at a face); retreat the
            // station inward along the tangent and remeasure, keeping the
            // face-to-station slab for the axial volume.
            let is_end = !closed_loop && (index == 0 || index + 1 == sampled.len());
            let mut center_now = *center;
            let mut retreat = 0.0f64;
            let mut result =
                plane_section(chart, center_now, *tangent, *pole_r, index, config, cx)?;
            if is_end {
                let inward = if index == 0 { 1.0 } else { -1.0 };
                let mut tries = 0;
                while result.1 < 0.4 * *pole_r && tries < 3 {
                    let step = 0.4 * *pole_r;
                    retreat += step;
                    center_now = center_now.offset(tangent.scale(inward * step));
                    if chart.eval(center_now, cx).signed_distance >= 0.0 {
                        return Err(BoreError::CenterOutsideLumen {
                            station: index,
                            arc_length_m: *s_len,
                        });
                    }
                    result =
                        plane_section(chart, center_now, *tangent, *pole_r, index, config, cx)?;
                    tries += 1;
                }
            }
            (result.0, result.1, result.2, result.3.clone(), retreat)
        };
        if retreat_m > 0.0 {
            if index == 0 {
                extra_start_slab_m = retreat_m;
            } else {
                extra_end_slab_m = retreat_m;
            }
        }
        let pole_area = core::f64::consts::PI * pole_r * pole_r;
        let pole_dev = (area - pole_area).abs() / area.max(f64::MIN_POSITIVE);
        let (thickness, t_dev) = if config.thickness_cross_check {
            let mut outcome = (None, None);
            for bp in &boundary_points {
                match thickness_at(chart, *bp, cx) {
                    Ok(t) => {
                        let t_area = core::f64::consts::PI * (0.5 * t.value) * (0.5 * t.value);
                        outcome = (
                            Some(t.value),
                            Some((t_area - area).abs() / area.max(f64::MIN_POSITIVE)),
                        );
                        break;
                    }
                    Err(QueryError::Cancelled) => return Err(QueryError::Cancelled.into()),
                    Err(_) => {}
                }
            }
            if outcome.0.is_none() {
                thickness_skipped += 1;
            }
            outcome
        } else {
            (None, None)
        };
        stations.push(BoreStation {
            arc_length_m: *s_len,
            center_m: *center,
            tangent: *tangent,
            area_m2: area,
            equivalent_radius_m: (area / core::f64::consts::PI).sqrt(),
            pole_radius_m: *pole_r,
            pole_area_deviation: pole_dev,
            thickness_m: thickness,
            thickness_area_deviation: t_dev,
            min_ray_m: min_ray,
            max_ray_m: max_ray,
        });
    }

    // 7. Axial volume (trapezoid; loops add the wrap segment; open tubes
    //    add the constant-area end slabs between each face and its inset
    //    end station).
    let axial_volume_m3 = axial_volume(&stations, closed_loop)
        + stations[0].area_m2 * (start_slab_m + extra_start_slab_m)
        + stations[stations.len() - 1].area_m2 * (end_slab_m + extra_end_slab_m);
    let total_length_m = if closed_loop {
        let wrap = stations[0]
            .center_m
            .delta_from(stations[stations.len() - 1].center_m)
            .norm();
        stations[stations.len() - 1].arc_length_m + wrap
    } else {
        stations[stations.len() - 1].arc_length_m
    };

    // 8. Certified volume closure (or the honest reason there is none).
    let volume_closure = match config.closure_h_m {
        None => VolumeClosure::Skipped,
        Some(h) => {
            let claim = chart.trace_step_claim();
            if claim == TraceStepClaim::ExactDistance {
                let support = chart.support();
                let pad = 2.0 * h;
                let domain = Aabb::new(
                    support.min.offset(Vec3::new(-pad, -pad, -pad)),
                    support.max.offset(Vec3::new(pad, pad, pad)),
                );
                let moments = geometric_moments(chart, &domain, h, cx)?;
                let (lo, hi) = (moments.volume.lo, moments.volume.hi);
                let mid = 0.5 * (lo + hi);
                let rel = if mid > 0.0 {
                    (axial_volume_m3 - mid).abs() / mid
                } else {
                    f64::INFINITY
                };
                VolumeClosure::Certified {
                    lo_m3: lo,
                    hi_m3: hi,
                    h_m: h,
                    axial_vs_mid_rel: rel,
                }
            } else {
                VolumeClosure::Unavailable { claim }
            }
        }
    };

    Ok(BoreExtraction {
        stations,
        closed_loop,
        total_length_m,
        axial_volume_m3,
        volume_closure,
        source_geometry_id: source_geometry_id.to_string(),
        boundary_digest: digest,
        thickness_skipped,
        authority: NumericalKind::Estimate,
    })
}

fn validate_config(config: &BoreConfig) -> Result<(), BoreError> {
    if !(config.lambda.is_finite() && config.lambda >= 0.0) {
        return Err(BoreError::InvalidConfig {
            what: "lambda must be finite and non-negative",
        });
    }
    if config.stations < 3 {
        return Err(BoreError::InvalidConfig {
            what: "at least three stations",
        });
    }
    if config.rays_per_station < 8 {
        return Err(BoreError::InvalidConfig {
            what: "at least eight rays per station",
        });
    }
    if !(config.root_tolerance_m.is_finite() && config.root_tolerance_m > 0.0) {
        return Err(BoreError::InvalidConfig {
            what: "root tolerance must be finite and positive",
        });
    }
    if config.max_poles < 8 {
        return Err(BoreError::InvalidConfig {
            what: "pole thinning cap below eight",
        });
    }
    if !(config.branch_reach_factor.is_finite() && config.branch_reach_factor > 0.0) {
        return Err(BoreError::InvalidConfig {
            what: "branch tolerance must be finite and positive",
        });
    }
    if !(config.ray_cap_factor.is_finite() && config.ray_cap_factor > 1.0) {
        return Err(BoreError::InvalidConfig {
            what: "ray cap factor must exceed one",
        });
    }
    if let Some(h) = config.closure_h_m
        && !(h.is_finite() && h > 0.0)
    {
        return Err(BoreError::InvalidConfig {
            what: "closure grid spacing must be finite and positive",
        });
    }
    Ok(())
}

/// FNV-1a over the soup's position bits and triangle indices.
fn soup_digest(soup: &Soup) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            hash ^= u64::from(b);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for p in &soup.positions {
        eat(&p.x.to_bits().to_le_bytes());
        eat(&p.y.to_bits().to_le_bytes());
        eat(&p.z.to_bits().to_le_bytes());
    }
    for t in &soup.triangles {
        for &i in t {
            eat(&i.to_le_bytes());
        }
    }
    hash
}

/// Keep only poles whose radius is within 10% of the largest radius in
/// their own-radius neighborhood — the local-maximality property that
/// distinguishes 1-D spine poles from medial SHEET poles (end discs).
fn spinal_filter(poles: &[(Point3, f64)], cx: &Cx<'_>) -> Result<Vec<(Point3, f64)>, BoreError> {
    let mut kept = Vec::with_capacity(poles.len());
    for &(p, r) in poles {
        crate::query_checkpoint(cx)?;
        let mut local_max = r;
        for &(q, rq) in poles {
            if q.delta_from(p).norm() < r {
                local_max = local_max.max(rq);
            }
        }
        if r >= 0.9 * local_max {
            kept.push((p, r));
        }
    }
    Ok(kept)
}

/// Deterministic thinning: total-order sort, then stride sampling.
fn thin_poles(mut poles: Vec<(Point3, f64)>, cap: usize) -> Vec<(Point3, f64)> {
    poles.sort_by(|a, b| {
        a.0.x
            .total_cmp(&b.0.x)
            .then(a.0.y.total_cmp(&b.0.y))
            .then(a.0.z.total_cmp(&b.0.z))
    });
    if poles.len() <= cap {
        return poles;
    }
    let stride = poles.len().div_ceil(cap);
    poles.into_iter().step_by(stride).collect()
}

/// Prim MST over the full Euclidean graph (O(n^2); n is capped upstream).
fn mst_edges(poles: &[(Point3, f64)], cx: &Cx<'_>) -> Result<Vec<(usize, usize)>, BoreError> {
    let n = poles.len();
    let mut in_tree = vec![false; n];
    let mut best = vec![f64::INFINITY; n];
    let mut best_from = vec![0usize; n];
    let mut edges = Vec::with_capacity(n - 1);
    in_tree[0] = true;
    for j in 1..n {
        best[j] = poles[j].0.delta_from(poles[0].0).norm();
    }
    for _ in 1..n {
        crate::query_checkpoint(cx)?;
        let mut pick = usize::MAX;
        let mut pick_d = f64::INFINITY;
        for j in 0..n {
            if !in_tree[j] && best[j] < pick_d {
                pick = j;
                pick_d = best[j];
            }
        }
        if pick == usize::MAX {
            return Err(BoreError::DegeneratePoleChain {
                reason: "disconnected pole cloud under finite distances",
            });
        }
        in_tree[pick] = true;
        edges.push((best_from[pick], pick));
        for j in 0..n {
            if !in_tree[j] {
                let d = poles[j].0.delta_from(poles[pick].0).norm();
                if d < best[j] {
                    best[j] = d;
                    best_from[j] = pick;
                }
            }
        }
    }
    Ok(edges)
}

fn adjacency(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        adj[a].push(b);
        adj[b].push(a);
    }
    adj
}

/// Farthest node from `start` by tree distance, with parents for path
/// recovery.
fn farthest(
    adj: &[Vec<usize>],
    poles: &[(Point3, f64)],
    start: usize,
) -> (usize, Vec<Option<usize>>) {
    let n = adj.len();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut dist = vec![f64::NEG_INFINITY; n];
    let mut stack = vec![start];
    dist[start] = 0.0;
    parent[start] = Some(start);
    while let Some(u) = stack.pop() {
        for &v in &adj[u] {
            if parent[v].is_none() {
                parent[v] = Some(u);
                dist[v] = dist[u] + poles[v].0.delta_from(poles[u].0).norm();
                stack.push(v);
            }
        }
    }
    let mut far = start;
    for (j, &d) in dist.iter().enumerate() {
        if d > dist[far] {
            far = j;
        }
    }
    parent[start] = None;
    (far, parent)
}

fn longest_path_poles(poles: &[(Point3, f64)], edges: &[(usize, usize)]) -> Vec<usize> {
    let adj = adjacency(poles.len(), edges);
    let (u, _) = farthest(&adj, poles, 0);
    let (v, parent) = farthest(&adj, poles, u);
    let mut path = vec![v];
    let mut cur = v;
    while let Some(p) = parent[cur] {
        path.push(p);
        cur = p;
    }
    path.reverse();
    path
}

/// Refuse when an off-spine subtree REACHES farther from the spine than
/// `factor` times the local pole radius: that is a BRANCH, not medial
/// noise (noise poles hug the tube interior; a genuine side arm extends
/// multiple radii away).
fn audit_branches(
    poles: &[(Point3, f64)],
    edges: &[(usize, usize)],
    spine: &[usize],
    factor: f64,
) -> Result<(), BoreError> {
    let n = poles.len();
    let adj = adjacency(n, edges);
    let mut on_spine = vec![false; n];
    for &i in spine {
        on_spine[i] = true;
    }
    let dist_to_spine = |p: Point3| -> f64 {
        spine
            .iter()
            .map(|&s| poles[s].0.delta_from(p).norm())
            .fold(f64::INFINITY, f64::min)
    };
    for &s in spine {
        for &nb in &adj[s] {
            if on_spine[nb] {
                continue;
            }
            // Farthest reach of the subtree hanging off the spine at s.
            let mut reach = 0.0f64;
            let mut reach_len = 0.0f64;
            let mut seen = vec![false; n];
            seen[s] = true;
            seen[nb] = true;
            let mut stack = vec![nb];
            while let Some(u) = stack.pop() {
                reach = reach.max(dist_to_spine(poles[u].0));
                for &v in &adj[u] {
                    if !seen[v] && !on_spine[v] {
                        seen[v] = true;
                        reach_len += poles[v].0.delta_from(poles[u].0).norm();
                        stack.push(v);
                    }
                }
            }
            if reach > factor * poles[s].1 {
                return Err(BoreError::BranchedLumen {
                    at: [poles[s].0.x, poles[s].0.y, poles[s].0.z],
                    side_length_m: reach_len.max(reach),
                });
            }
        }
    }
    Ok(())
}

fn polyline_length(poles: &[(Point3, f64)], spine: &[usize]) -> f64 {
    spine
        .windows(2)
        .map(|w| poles[w[1]].0.delta_from(poles[w[0]].0).norm())
        .sum()
}

/// Radius^2-weighted pole centroids binned by arc length along the
/// SMOOTHED spine polyline (segment projection, not nearest node — the raw
/// MST diameter path zigzags through the 3-D pole thicket and nearest-node
/// parameterization scrambles the bins). The result is the cloud-snapped
/// centerline chain. Empty bins are skipped.
fn bin_centroids(
    poles: &[(Point3, f64)],
    spine: &[usize],
    bins: usize,
    closed: bool,
    cx: &Cx<'_>,
) -> Result<Vec<(Point3, f64)>, BoreError> {
    // Heavily smoothed spine as the parameterization skeleton.
    let spine_pts: Vec<(Point3, f64)> = spine.iter().map(|&i| poles[i]).collect();
    let mut skeleton = smooth_chain(&spine_pts, 3, closed);
    // Near a flat tube end the medial set is a disc-shaped SHEET reaching
    // about one local radius in from the cap; the MST diameter path curls
    // through it and scrambles the tail ordering. Trim one end radius of
    // skeleton arc from each OPEN end — the end-extension step rebuilds
    // the physical ends from the clean interior tangents.
    if !closed {
        let arc_at = |sk: &[(Point3, f64)]| -> Vec<f64> {
            let mut c = vec![0.0f64; sk.len()];
            for i in 1..sk.len() {
                c[i] = c[i - 1] + sk[i].0.delta_from(sk[i - 1].0).norm();
            }
            c
        };
        let c = arc_at(&skeleton);
        let total_raw = c[skeleton.len() - 1];
        let r_start = skeleton[0].1;
        let r_end = skeleton[skeleton.len() - 1].1;
        if total_raw > 3.0 * (r_start + r_end) {
            let keep: Vec<(Point3, f64)> = skeleton
                .iter()
                .zip(&c)
                .filter(|&(_, &a)| a >= r_start && a <= total_raw - r_end)
                .map(|(&node, _)| node)
                .collect();
            if keep.len() >= 3 {
                skeleton = keep;
            }
        }
    }
    let m = skeleton.len();
    let mut cum = vec![0.0f64; m];
    for i in 1..m {
        cum[i] = cum[i - 1] + skeleton[i].0.delta_from(skeleton[i - 1].0).norm();
    }
    let total = cum[m - 1];
    if !(total > 0.0) {
        return Err(BoreError::DegeneratePoleChain {
            reason: "zero-length spine",
        });
    }
    // Pass 1: arc-length parameter of every pole (segment projection) and
    // the per-bin maximum radius.
    let mut params = Vec::with_capacity(poles.len());
    let mut bin_max = vec![0.0f64; bins];
    for &(p, r) in poles {
        crate::query_checkpoint(cx)?;
        // Closest point on the skeleton polyline (segment projection).
        let mut best_d2 = f64::INFINITY;
        let mut best_s = 0.0f64;
        for i in 0..m - 1 {
            let (a, _) = skeleton[i];
            let (b, _) = skeleton[i + 1];
            let ab = b.delta_from(a);
            let len2 = ab.dot(ab);
            let t = if len2 > 0.0 {
                (p.delta_from(a).dot(ab) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let q = a.offset(ab.scale(t));
            let d = p.delta_from(q);
            let d2 = d.dot(d);
            if d2 < best_d2 {
                best_d2 = d2;
                best_s = cum[i] + t * (cum[i + 1] - cum[i]);
            }
        }
        let bin = ((best_s / total * bins as f64) as usize).min(bins - 1);
        params.push(bin);
        bin_max[bin] = bin_max[bin].max(r);
    }
    // Pass 2: accumulate only the locally-most-medial poles of each bin
    // (radius within 20% of the bin maximum). Near a flat tube end the
    // medial set is a disc-shaped SHEET, not a curve; its off-axis poles
    // have smaller inscribed radii, so this filter recovers the 1-D spine
    // — radius-local-maximality is exactly what makes a pole spinal.
    let mut acc: Vec<(f64, f64, f64, f64, f64)> = vec![(0.0, 0.0, 0.0, 0.0, 0.0); bins];
    for (&(p, r), &bin) in poles.iter().zip(&params) {
        if r < 0.8 * bin_max[bin] {
            continue;
        }
        let w = r * r;
        let slot = &mut acc[bin];
        slot.0 += w * p.x;
        slot.1 += w * p.y;
        slot.2 += w * p.z;
        slot.3 += w * r;
        slot.4 += w;
    }
    let chain: Vec<(Point3, f64)> = acc
        .into_iter()
        .filter(|slot| slot.4 > 0.0)
        .map(|(x, y, z, r, w)| (Point3::new(x / w, y / w, z / w), r / w))
        .collect();
    if chain.len() < 3 {
        return Err(BoreError::DegeneratePoleChain {
            reason: "fewer than three populated centerline bins",
        });
    }
    Ok(chain)
}

/// Moving-average smoothing; endpoints held for open chains, wrapped for
/// loops. Radii smoothed with the same window.
fn smooth_chain(chain: &[(Point3, f64)], half_width: usize, closed: bool) -> Vec<(Point3, f64)> {
    if half_width == 0 || chain.len() < 3 {
        return chain.to_vec();
    }
    let n = chain.len();
    let idx = |i: isize| -> usize {
        if closed {
            (i.rem_euclid(n as isize)) as usize
        } else {
            i.clamp(0, n as isize - 1) as usize
        }
    };
    (0..n)
        .map(|i| {
            if !closed && (i == 0 || i == n - 1) {
                return chain[i];
            }
            let mut px = 0.0;
            let mut py = 0.0;
            let mut pz = 0.0;
            let mut pr = 0.0;
            let mut count = 0.0;
            let hw = half_width as isize;
            for o in -hw..=hw {
                let (p, r) = chain[idx(i as isize + o)];
                px += p.x;
                py += p.y;
                pz += p.z;
                pr += r;
                count += 1.0;
            }
            (Point3::new(px / count, py / count, pz / count), pr / count)
        })
        .collect()
}

/// March each open end along the averaged interior end tangent to the
/// region boundary, then place the end station INSET one tenth of the
/// local radius inside the face (a station sitting on the face clips its
/// own cross-section plane against the cap under any tangent tilt). The
/// returned pair is the (start, end) slab length between each face and
/// its inset station, so the axial volume can cover the full tube.
fn extend_to_boundary(
    chart: &dyn Chart,
    chain: &mut Vec<(Point3, f64)>,
    config: &BoreConfig,
    cx: &Cx<'_>,
) -> Result<(f64, f64), BoreError> {
    let mut slabs = [0.0f64; 2];
    for (slot, end) in [(0usize, false), (1usize, true)] {
        let n = chain.len();
        let (tip, r) = if end {
            (chain[n - 1].0, chain[n - 1].1)
        } else {
            (chain[0].0, chain[0].1)
        };
        // Averaged direction of the last few interior segments (a single
        // noisy bin-to-bin hop tilts the end station's section plane).
        let span = 3.min(n - 1);
        let mut dir = Vec3::new(0.0, 0.0, 0.0);
        for k in 0..span {
            let (a, b) = if end {
                (chain[n - 2 - k].0, chain[n - 1 - k].0)
            } else {
                (chain[1 + k].0, chain[k].0)
            };
            let d = b.delta_from(a);
            let l = d.norm();
            if l > 0.0 {
                let du = d.scale(1.0 / l);
                dir = Vec3::new(dir.x + du.x, dir.y + du.y, dir.z + du.z);
            }
        }
        let len = dir.norm();
        if !(len > 0.0) {
            continue;
        }
        let dir = dir.scale(1.0 / len);
        if chart.eval(tip, cx).signed_distance >= 0.0 {
            continue; // tip already at/outside the boundary; nothing to add
        }
        // Expand to an outside bracket.
        let cap = 4.0 * r.max(config.root_tolerance_m);
        let mut hi = 0.25 * r.max(config.root_tolerance_m);
        let mut found = false;
        while hi <= cap {
            crate::query_checkpoint(cx)?;
            if chart.eval(tip.offset(dir.scale(hi)), cx).signed_distance >= 0.0 {
                found = true;
                break;
            }
            hi *= 2.0;
        }
        if !found {
            continue; // no nearby end face; leave the chain as chained
        }
        let mut lo = 0.0f64;
        while hi - lo > config.root_tolerance_m {
            crate::query_checkpoint(cx)?;
            let mid = 0.5 * (lo + hi);
            if chart.eval(tip.offset(dir.scale(mid)), cx).signed_distance < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        // Station inset from the face; the face-to-station slab length is
        // recorded for the axial volume.
        let inset = 0.1 * r;
        let t = (lo - inset).max(0.0);
        let face = tip.offset(dir.scale(t));
        slabs[slot] = (lo - t).max(0.0);
        if end {
            chain.push((face, r));
        } else {
            chain.insert(0, (face, r));
        }
    }
    Ok((slabs[0], slabs[1]))
}

/// Resample the chain into `stations` points by arc length; returns
/// (arc length, center, unit tangent, interpolated pole radius).
fn resample(
    chain: &[(Point3, f64)],
    stations: usize,
    closed: bool,
) -> Vec<(f64, Point3, Vec3, f64)> {
    let n = chain.len();
    // Cumulative arc length (loops append the wrap segment).
    let mut cum = vec![0.0f64; n];
    for i in 1..n {
        cum[i] = cum[i - 1] + chain[i].0.delta_from(chain[i - 1].0).norm();
    }
    let total = if closed {
        cum[n - 1] + chain[0].0.delta_from(chain[n - 1].0).norm()
    } else {
        cum[n - 1]
    };
    let at = |s: f64| -> (Point3, f64) {
        // Locate the segment containing s (wrap for loops).
        let s = if closed { s % total } else { s.min(cum[n - 1]) };
        let mut seg = 0;
        while seg + 1 < n && cum[seg + 1] < s {
            seg += 1;
        }
        let (a, ra) = chain[seg];
        let (b, rb, sa, sb) = if seg + 1 < n {
            (chain[seg + 1].0, chain[seg + 1].1, cum[seg], cum[seg + 1])
        } else {
            (chain[0].0, chain[0].1, cum[n - 1], total)
        };
        let span = (sb - sa).max(f64::MIN_POSITIVE);
        let f = ((s - sa) / span).clamp(0.0, 1.0);
        let d = b.delta_from(a);
        (a.offset(d.scale(f)), ra + f * (rb - ra))
    };
    let count = stations;
    (0..count)
        .map(|k| {
            let s = if closed {
                total * k as f64 / count as f64
            } else {
                total * k as f64 / (count - 1) as f64
            };
            let (p, r) = at(s);
            // Tangent by symmetric difference at a small arc offset.
            let ds = (total / count as f64).max(f64::MIN_POSITIVE);
            let (pf, _) = at(if closed { s + ds } else { (s + ds).min(total) });
            let (pb, _) = at(if closed {
                (s - ds).rem_euclid(total)
            } else {
                (s - ds).max(0.0)
            });
            let d = pf.delta_from(pb);
            let nrm = d.norm().max(f64::MIN_POSITIVE);
            (s, p, d.scale(1.0 / nrm), r)
        })
        .collect()
}

/// Polar plane-section: bisected boundary rays in the plane normal to the
/// tangent; area by the exact polygonal polar quadrature `sum r^2/2 dtheta`.
/// Returns (area, min ray, max ray, candidate boundary points for the
/// thickness cross-check).
fn plane_section(
    chart: &dyn Chart,
    center: Point3,
    tangent: Vec3,
    pole_radius: f64,
    station: usize,
    config: &BoreConfig,
    cx: &Cx<'_>,
) -> Result<(f64, f64, f64, Vec<Point3>), BoreError> {
    // Orthonormal in-plane frame.
    let seed = if tangent.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = {
        let proj = tangent.scale(seed.dot(tangent));
        let raw = Vec3::new(seed.x - proj.x, seed.y - proj.y, seed.z - proj.z);
        raw.scale(1.0 / raw.norm().max(f64::MIN_POSITIVE))
    };
    let v = cross(tangent, u);
    let rays = config.rays_per_station;
    let cap = config.ray_cap_factor * pole_radius.max(config.root_tolerance_m);
    let mut area = 0.0f64;
    let mut min_ray = f64::INFINITY;
    let mut max_ray = 0.0f64;
    let mut boundary_points = Vec::new();
    for k in 0..rays {
        crate::query_checkpoint(cx)?;
        let theta = core::f64::consts::TAU * k as f64 / rays as f64;
        let dir = {
            let (s, c) = theta.sin_cos();
            Vec3::new(c * u.x + s * v.x, c * u.y + s * v.y, c * u.z + s * v.z)
        };
        // Expand to an outside bracket.
        let mut hi = 0.5 * pole_radius.max(config.root_tolerance_m);
        let mut found = false;
        while hi <= cap {
            crate::query_checkpoint(cx)?;
            if chart.eval(center.offset(dir.scale(hi)), cx).signed_distance >= 0.0 {
                found = true;
                break;
            }
            hi *= 2.0;
        }
        if !found {
            return Err(BoreError::RayNoExit {
                station,
                ray: k,
                center: [center.x, center.y, center.z],
                tangent: [tangent.x, tangent.y, tangent.z],
                cap_m: cap,
                pole_radius_m: pole_radius,
            });
        }
        let mut lo = 0.0f64;
        while hi - lo > config.root_tolerance_m {
            let mid = 0.5 * (lo + hi);
            if chart
                .eval(center.offset(dir.scale(mid)), cx)
                .signed_distance
                < 0.0
            {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let r = 0.5 * (lo + hi);
        if k % (rays / 8).max(1) == 0 && boundary_points.len() < 4 {
            // The thickness cross-check needs points within the boundary
            // tolerance of downstream queries: refine a few rays much
            // deeper and keep the strictly-INSIDE endpoints. Several
            // angles are kept because a diametral chord of a shape whose
            // support box is tight ends ON the sampling-domain boundary,
            // where `thickness_at` correctly refuses (NoOppositeWall);
            // an oblique chord exits in the domain interior.
            while hi - lo > config.root_tolerance_m / 64.0 {
                let mid = 0.5 * (lo + hi);
                if chart
                    .eval(center.offset(dir.scale(mid)), cx)
                    .signed_distance
                    < 0.0
                {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            boundary_points.push(center.offset(dir.scale(lo)));
        }
        min_ray = min_ray.min(r);
        max_ray = max_ray.max(r);
        area += 0.5 * r * r * (core::f64::consts::TAU / rays as f64);
    }
    Ok((area, min_ray, max_ray, boundary_points))
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

/// Trapezoidal `integral A ds` over the stations.
fn axial_volume(stations: &[BoreStation], closed: bool) -> f64 {
    let mut vol = 0.0f64;
    for w in stations.windows(2) {
        let ds = w[1].arc_length_m - w[0].arc_length_m;
        vol += 0.5 * (w[0].area_m2 + w[1].area_m2) * ds;
    }
    if closed {
        let first = &stations[0];
        let last = &stations[stations.len() - 1];
        let ds = last.center_m.delta_from(first.center_m).norm();
        vol += 0.5 * (first.area_m2 + last.area_m2) * ds;
    }
    vol
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> (Point3, f64) {
        (Point3::new(x, y, z), 0.5)
    }

    #[test]
    fn chaining_orders_a_noisy_line() {
        // Poles along x with slight jitter, shuffled by construction of
        // the sort (thin_poles sorts by x anyway) — MST + diameter must
        // recover the ordered spine.
        let poles: Vec<(Point3, f64)> = (0..20)
            .map(|i| {
                let x = i as f64 * 0.3;
                let j = if i % 2 == 0 { 0.02 } else { -0.02 };
                (Point3::new(x, j, 0.0), 0.5)
            })
            .collect();
        let gate = fs_exec::CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0,
                    kernel_id: 900,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let edges = mst_edges(&poles, &cx).expect("mst");
            let spine = longest_path_poles(&poles, &edges);
            assert_eq!(spine.len(), 20, "the spine must cover every pole");
            let xs: Vec<f64> = spine.iter().map(|&i| poles[i].0.x).collect();
            let increasing = xs.windows(2).all(|w| w[1] > w[0]);
            let decreasing = xs.windows(2).all(|w| w[1] < w[0]);
            assert!(increasing || decreasing, "spine must be x-monotone: {xs:?}");
        });
    }

    #[test]
    fn a_y_junction_refuses_as_branched() {
        // A main line plus a genuine side arm several radii long.
        let mut poles: Vec<(Point3, f64)> = (0..15).map(|i| p(i as f64 * 0.3, 0.0, 0.0)).collect();
        for j in 1..=6 {
            poles.push(p(2.1, j as f64 * 0.3, 0.0));
        }
        let gate = fs_exec::CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0,
                    kernel_id: 901,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let edges = mst_edges(&poles, &cx).expect("mst");
            let spine = longest_path_poles(&poles, &edges);
            let refusal = audit_branches(&poles, &edges, &spine, 2.0);
            assert!(
                matches!(refusal, Err(BoreError::BranchedLumen { .. })),
                "a six-pole side arm must refuse: {refusal:?}"
            );
        });
    }

    #[test]
    fn smoothing_holds_open_endpoints_fixed() {
        let chain: Vec<(Point3, f64)> = (0..9)
            .map(|i| p(i as f64, (i % 2) as f64 * 0.2, 0.0))
            .collect();
        let out = smooth_chain(&chain, 2, false);
        assert_eq!(out[0].0.x, chain[0].0.x);
        assert_eq!(out[8].0.x, chain[8].0.x);
        // Interior jitter must shrink.
        let jitter_in: f64 = chain[1..8].iter().map(|c| c.0.y).sum::<f64>();
        let jitter_out: f64 = out[1..8].iter().map(|c| c.0.y).sum::<f64>();
        assert!(
            (jitter_out - jitter_in).abs() < 1.0e-12 || jitter_out < jitter_in,
            "smoothing must not amplify jitter"
        );
    }
}
