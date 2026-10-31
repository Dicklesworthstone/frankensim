//! Numerical enclosures for the original direct point-fixture calculation.
//!
//! The enclosed real function uses the retained binary64 flux, normalized
//! angular coefficients/cutoffs, inverse-rotation quaternion, source/probe
//! positions, receiver normal and gray absorption/scattering coefficients.
//! Native opaque triangles and their retained forward rigid placements determine
//! visibility on the source-to-receiver segment, excluding the source endpoint
//! and the final `ray_epsilon_m` at the receiver, as in the nominal calculation.
//! A visibility predicate which cannot decide an edge/endpoint/parallel case
//! retains both possibilities; it never silently becomes a clear ray.
//!
//! Constructors' normalization, filtering and input-conversion errors are NOT
//! included. Neither are physical calibration, geometry fidelity or the physical
//! effect of the excluded receiver segment, spectral/model discrepancy, or
//! omitted reflected/scattered light. The total error therefore stays unknown.
//!
//! Physical-grid transport has no enclosure here. Supported source laws are
//! isotropic, uniform cone, integer cosine powers 0..=64 and two-column
//! axisymmetric tables with at most 64 polar knots. Cutoff uncertainty encloses
//! both branches. The original native nominal calculation is unchanged.

use super::{AngularLaw, PhotometricSource};
use crate::photometry::{
    CalculationGrid, IlluminanceEstimate, IlluminanceProbe, PhotometricScene,
    PhotometryError, Receivers,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_ivl::{Interval, decay::exp_neg};
use std::io::{self, Write};

mod visibility;
use visibility::{Visibility, VisibilityIndex};
pub use visibility::MAX_VISIBILITY_TRIANGLES;

/// Version of the represented-function and numerical-enclosure interpretation.
/// Version 2 adds bounded opaque-mesh visibility and geometry-bound identities.
pub const DIRECT_ENCLOSURE_SEMANTICS_VERSION: u32 = 2;

/// Why the unchanged nominal result has no numerical enclosure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnclosureUnavailable {
    /// Unsupported geometry or a nonrepresentable world-coordinate bound.
    GeometryVisibility,
    /// Optional visibility metadata exceeds the documented admission ceiling.
    GeometryBudgetExceeded,
    /// Heterogeneous spectral attenuation/cell crossings are not enclosed.
    PhysicalGridTransport,
    /// This source needs an unsupported angular law or a larger table.
    SourceLaw {
        /// Source index in the original caller order.
        source: usize,
    },
    /// An intermediate or absolute error bound could not be enclosed finitely.
    ArithmeticRange,
}

/// Direct lux with an optional enclosure of the retained numerical model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnclosedIlluminance {
    /// Original nominal lux; its numerical field is promoted only on success.
    /// Its total-error bound is always still unknown.
    pub estimate: IlluminanceEstimate,
    /// Outward bound on |nominal - exact retained-model lux|, not physical error.
    /// This includes any unresolved visibility branches, not just small rounding.
    pub numerical_error_bound_lux: Option<f64>,
    /// Content identity of retained source, referenced geometry and probe inputs.
    /// Spectra/phase are irrelevant to gray unscattered lux and are excluded,
    /// as are inactive source details and constructor history. This is not a
    /// complete scene hash or authenticated calibration provenance.
    pub numerical_model_identity: Option<ContentHash>,
    /// Sources with nonzero possible contribution and unresolved visibility.
    /// `Some(0)` means all relevant visibility predicates were decided; `None`
    /// means no enclosure was produced. Unknown branches widen the lux interval.
    pub uncertain_visibility_sources: Option<usize>,
    /// Explicit absence reason, mutually exclusive with a numerical enclosure.
    pub unavailable: Option<EnclosureUnavailable>,
}

impl EnclosedIlluminance {
    /// Header for the direct numerical-accounting CSV, distinct from path CSV.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,direct_lux,numerical_lower_lux,numerical_upper_lux,numerical_error_bound_lux,numerical_status,numerical_model_identity,source_construction_error_bound_lux,calibration_error_bound_lux,physical_model_error_bound_lux,total_error_bound_lux,enclosure_version,uncertain_visibility_sources\n")
    }

    /// Stream round-trip numeric values with literal unknown physical bounds.
    /// An uncertain visibility branch can give a wide, but valid, enclosure.
    /// Failure propagates; atomic publication belongs to the caller's writer.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        write!(writer, "{index},{},", self.estimate.lux)?;
        match self.numerical_error_bound_lux {
            Some(error) => write!(writer, "{},{},{error},enclosure,",
                self.estimate.numerical.lo, self.estimate.numerical.hi)?,
            None => {
                let reason = match self.unavailable {
                    Some(EnclosureUnavailable::GeometryVisibility) => "geometry-visibility",
                    Some(EnclosureUnavailable::GeometryBudgetExceeded) => "geometry-budget",
                    Some(EnclosureUnavailable::PhysicalGridTransport) => "physical-grid-transport",
                    Some(EnclosureUnavailable::SourceLaw { .. }) => "source-law",
                    _ => "arithmetic-range",
                };
                write!(writer, "unknown,unknown,unknown,{reason},")?;
            }
        }
        match self.numerical_model_identity {
            Some(identity) => write!(writer, "{identity}")?,
            None => writer.write_all(b"unknown")?,
        }
        write!(writer, ",unknown,unknown,unknown,unknown,{DIRECT_ENCLOSURE_SEMANTICS_VERSION},")?;
        match self.uncertain_visibility_sources {
            Some(count) => writeln!(writer, "{count}"),
            None => writer.write_all(b"unknown\n"),
        }
    }
}

enum Prepared {
    Available { identity: ContentHash, visibility: VisibilityIndex },
    Unavailable(EnclosureUnavailable),
}

impl PhotometricScene<'_> {
    /// Evaluate the original direct lux plus a supported arithmetic enclosure.
    /// Unsupported enclosure physics does not reject an otherwise valid nominal
    /// calculation: inspect `unavailable`. Cancellation and nominal transport
    /// refusals still return errors. No total physical-accuracy claim is made.
    /// Use [`EnclosedIlluminanceJob`] to amortize visibility preparation over
    /// multiple probes or an entire calculation grid.
    pub fn illuminance_with_enclosure(
        &self, cx: &Cx<'_>, probe: &IlluminanceProbe,
    ) -> Result<EnclosedIlluminance, PhotometryError> {
        evaluate(self, cx, probe, &prepare(self, cx)?)
    }
}

/// Progressive direct calculation with enclosure and provenance on each value.
/// The original scene/probes stay borrowed. The numerical-model identity and
/// bounded, stackless visibility hierarchy are prepared once, not per probe.
/// Each receiver publishes atomically; failed/cancelled receivers are absent.
pub struct EnclosedIlluminanceJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    receivers: Receivers<'a>,
    prepared: Prepared,
    count: usize,
    values: Vec<EnclosedIlluminance>,
}

impl<'a, 'scene> EnclosedIlluminanceJob<'a, 'scene> {
    /// Admit explicit probe positions without copying them.
    pub fn new(scene: &'a PhotometricScene<'scene>, probes: &'a [IlluminanceProbe], cx: &Cx<'_>) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Probes(probes), probes.len(), cx)
    }

    /// Admit the native point-sampled calculation plane without a probe copy.
    pub fn for_grid(scene: &'a PhotometricScene<'scene>, grid: CalculationGrid, cx: &Cx<'_>) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Grid(grid), grid.len(), cx)
    }

    fn begin(scene: &'a PhotometricScene<'scene>, receivers: Receivers<'a>, count: usize, cx: &Cx<'_>) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if count > scene.options.max_outputs { return Err(PhotometryError::OutputBudgetExceeded); }
        let prepared = prepare(scene, cx)?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| PhotometryError::AllocationRefused)?;
        cx.checkpoint()?;
        Ok(Self { scene, receivers, prepared, count, values })
    }

    /// Append only complete receivers. Chunk sizes do not enter the arithmetic.
    pub fn advance(&mut self, cx: &Cx<'_>, max_points: usize) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_points).min(self.count);
        for index in start..end {
            let probe = match &self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index)?,
            };
            let value = evaluate(self.scene, cx, &probe, &self.prepared)?;
            cx.checkpoint()?;
            self.values.push(value);
        }
        Ok(self.values.len() - start)
    }

    /// Successfully completed prefix only.
    #[must_use]
    pub fn results(&self) -> &[EnclosedIlluminance] { &self.values }

    /// Whether every admitted receiver has a complete result.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.values.len() == self.count }
}

fn prepare(scene: &PhotometricScene<'_>, cx: &Cx<'_>) -> Result<Prepared, PhotometryError> {
    cx.checkpoint()?;
    let Some(medium) = scene.medium.homogeneous() else {
        return Ok(Prepared::Unavailable(EnclosureUnavailable::PhysicalGridTransport));
    };
    let mut hash = DomainHasher::new("org.frankensim.photometry.direct-enclosure.model.v2");
    hash.update(&DIRECT_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&fs_ivl::INTERVAL_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&fs_ivl::decay::DECAY_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&crate::photometry::PHOTOMETRY_SEMANTICS_VERSION.to_le_bytes());
    for value in [medium.absorption_per_m(), medium.scattering_per_m(), scene.options.ray_epsilon_m, scene.options.max_distance_m] {
        hash_float(&mut hash, value);
    }
    hash.update(&(scene.sources.len() as u64).to_le_bytes());
    for (index, source) in scene.sources.iter().enumerate() {
        cx.checkpoint()?;
        hash_float(&mut hash, source.emitted_flux_lm);
        if source.emitted_flux_lm == 0.0 { continue; }
        for value in source.transform.translation_m().into_iter().chain(source.world_to_local_rotation.rotation_xyzw()) {
            hash_float(&mut hash, value);
        }
        hash_float(&mut hash, source.distribution.peak_per_sr);
        match &source.distribution.law {
            AngularLaw::Isotropic => { hash.update(&[0]); },
            AngularLaw::Cone { cosine_cutoff } => {
                hash.update(&[1]); hash_float(&mut hash, *cosine_cutoff);
            }
            AngularLaw::CosinePower { exponent, cosine_cutoff } if *exponent <= 64.0 && exponent.fract() == 0.0 => {
                hash.update(&[2]); hash_float(&mut hash, *exponent); hash_float(&mut hash, *cosine_cutoff);
            }
            AngularLaw::Table(table) if table.phi.len() == 2 && table.cosine.len() <= 64 => {
                hash.update(&[3]);
                hash.update(&(table.cosine.len() as u64).to_le_bytes());
                for (row, mu) in table.cosine.iter().enumerate() {
                    hash_float(&mut hash, *mu); hash_float(&mut hash, table.values[2 * row]);
                }
            }
            _ => return Ok(Prepared::Unavailable(EnclosureUnavailable::SourceLaw { source: index })),
        }
    }
    let visibility = match VisibilityIndex::prepare(scene.primitives, &mut hash, cx)? {
        Ok(index) => index,
        Err(reason) => return Ok(Prepared::Unavailable(reason)),
    };
    cx.checkpoint()?;
    Ok(Prepared::Available { identity: hash.finalize(), visibility })
}

fn evaluate(scene: &PhotometricScene<'_>, cx: &Cx<'_>, probe: &IlluminanceProbe, prepared: &Prepared) -> Result<EnclosedIlluminance, PhotometryError> {
    let mut estimate = scene.illuminance(cx, probe)?;
    cx.checkpoint()?;
    let no_claim = |reason| EnclosedIlluminance {
        estimate, numerical_error_bound_lux: None, numerical_model_identity: None,
        uncertain_visibility_sources: None, unavailable: Some(reason),
    };
    let (model, visibility) = match prepared {
        Prepared::Unavailable(reason) => return Ok(no_claim(*reason)),
        Prepared::Available { identity, visibility } => (identity, visibility),
    };
    let Some((enclosure, uncertain)) = direct_bound(scene, cx, probe, visibility)? else {
        return Ok(no_claim(EnclosureUnavailable::ArithmeticRange));
    };
    // Include the native nominal value without using it to establish the real
    // model enclosure. A nominal visibility disagreement widens, never narrows.
    let enclosure = enclosure.hull(Interval::point(estimate.lux));
    let error = if enclosure.lo() == enclosure.hi() { 0.0 } else {
        (Interval::point(estimate.lux) - Interval::point(enclosure.lo())).hi()
            .max((Interval::point(enclosure.hi()) - Interval::point(estimate.lux)).hi())
    };
    if !error.is_finite() { return Ok(no_claim(EnclosureUnavailable::ArithmeticRange)); }
    estimate.numerical = NumericalCertificate::enclosure(enclosure.lo(), enclosure.hi());
    let mut hash = DomainHasher::new("org.frankensim.photometry.direct-enclosure.probe.v2");
    hash.update(&model.0);
    for value in [probe.position_m.x, probe.position_m.y, probe.position_m.z, probe.normal.x, probe.normal.y, probe.normal.z] { hash_float(&mut hash, value); }
    cx.checkpoint()?;
    Ok(EnclosedIlluminance {
        estimate, numerical_error_bound_lux: Some(error), numerical_model_identity: Some(hash.finalize()),
        uncertain_visibility_sources: Some(uncertain), unavailable: None,
    })
}

fn direct_bound(scene: &PhotometricScene<'_>, cx: &Cx<'_>, probe: &IlluminanceProbe, visibility: &VisibilityIndex) -> Result<Option<(Interval, usize)>, PhotometryError> {
    let Some(medium) = scene.medium.homogeneous() else { return Ok(None); };
    let extinction = nonnegative(Interval::point(medium.absorption_per_m()) + Interval::point(medium.scattering_per_m()));
    let vacuum = medium.extinction_per_m() == 0.0;
    let mut sum = Interval::point(0.0);
    let mut uncertain = 0;
    for source in scene.sources {
        cx.checkpoint()?;
        if source.emitted_flux_lm == 0.0 { continue; }
        let origin = source.transform.translation_m();
        let receiver = [probe.position_m.x, probe.position_m.y, probe.position_m.z];
        let delta = std::array::from_fn(|i| Interval::point(origin[i]) - Interval::point(receiver[i]));
        let Some(distance) = norm(delta) else { return Ok(None); };
        if distance.lo() <= 0.0 { return Ok(None); }
        let direction = delta.map(|v| v / distance);
        let normal = [probe.normal.x, probe.normal.y, probe.normal.z].map(Interval::point);
        let cosine = nonnegative(dot(normal, direction));
        if cosine.hi() == 0.0 { continue; }
        let Some(intensity) = intensity_bound(source, direction.map(|v| -v)) else { return Ok(None); };
        if intensity.hi() == 0.0 { continue; }
        // Parameterize from source to receiver, exactly as native visibility.
        // The receiver exclusion changes only this segment, not distance² or
        // extinction. Interval end points preserve uncertain boundary crossings.
        let end = Interval::point(1.0) - Interval::point(scene.options.ray_epsilon_m) / distance;
        let end = Interval::new(end.lo().clamp(0.0, 1.0), end.hi().clamp(0.0, 1.0));
        let state = visibility.visibility(scene.primitives, origin.map(Interval::point), delta.map(|v| -v), end, cx)?;
        if state == Visibility::Blocked { continue; }
        let transmission = if vacuum { Interval::point(1.0) } else {
            let Some(value) = exp_neg(nonnegative(extinction * distance)) else { return Ok(None); };
            value
        };
        let mut contribution = nonnegative(intensity / distance / distance * cosine * transmission);
        if !finite(contribution) { return Ok(None); }
        if state == Visibility::Uncertain {
            contribution = Interval::new(0.0, contribution.hi());
            uncertain += 1;
        }
        sum = nonnegative(sum + contribution);
        if !finite(sum) { return Ok(None); }
    }
    cx.checkpoint()?;
    Ok(Some((sum, uncertain)))
}

fn intensity_bound(source: &PhotometricSource, outgoing: [Interval; 3]) -> Option<Interval> {
    let peak = Interval::point(source.distribution.peak_per_sr);
    if matches!(&source.distribution.law, AngularLaw::Isotropic) {
        return Some(nonnegative(Interval::point(source.emitted_flux_lm) * peak));
    }
    let [x, y, z, w] = source.world_to_local_rotation.rotation_xyzw().map(Interval::point);
    let q = [x, y, z];
    let twice = cross(q, outgoing).map(|v| Interval::point(2.0) * v);
    let double = cross(q, twice);
    let local = std::array::from_fn(|i| outgoing[i] + w * twice[i] + double[i]);
    let length = norm(local)?;
    if length.lo() <= 0.0 { return None; }
    let cosine = local[2] / length;
    if !finite(cosine) { return None; }
    let mu = Interval::new(cosine.lo().clamp(-1.0, 1.0), cosine.hi().clamp(-1.0, 1.0));
    let angular = match &source.distribution.law {
        AngularLaw::Isotropic => peak,
        AngularLaw::Cone { cosine_cutoff } => cutoff(peak, mu, *cosine_cutoff),
        AngularLaw::CosinePower { exponent, cosine_cutoff } => {
            let mut power = Interval::point(1.0);
            let positive_mu = nonnegative(mu);
            for _ in 0..(*exponent as u32) { power = nonnegative(power * positive_mu); }
            cutoff(nonnegative(peak * power), mu, *cosine_cutoff)
        }
        AngularLaw::Table(table) => {
            let mut result: Option<Interval> = None;
            for (row, cell) in table.cosine.windows(2).enumerate() {
                let lo = mu.lo().max(cell[1]);
                let hi = mu.hi().min(cell[0]);
                if lo > hi { continue; }
                let denominator = Interval::point(cell[0]) - Interval::point(cell[1]);
                if denominator.lo() <= 0.0 { return None; }
                let t = (Interval::point(cell[0]) - Interval::new(lo, hi)) / denominator;
                let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
                let a = Interval::point(table.values[2 * row]);
                let b = Interval::point(table.values[2 * row + 2]);
                let value = nonnegative((Interval::point(1.0) - t) * a + t * b);
                result = Some(result.map_or(value, |old| old.hull(value)));
            }
            result?
        }
    };
    let result = nonnegative(Interval::point(source.emitted_flux_lm) * angular);
    finite(result).then_some(result)
}

fn cutoff(value: Interval, mu: Interval, boundary: f64) -> Interval {
    if mu.hi() < boundary { Interval::point(0.0) }
    else if mu.lo() >= boundary { value }
    else { Interval::new(0.0, value.hi()) }
}

fn norm(vector: [Interval; 3]) -> Option<Interval> {
    if vector.iter().any(|v| !finite(*v)) { return None; }
    let scale = vector.iter().map(|v| v.lo().abs().max(v.hi().abs())).fold(0.0_f64, f64::max);
    if scale == 0.0 { return Some(Interval::point(0.0)); }
    let scaled = vector.map(|v| v / Interval::point(scale));
    let squares = scaled.map(|v| nonnegative(v * v));
    let length = nonnegative(squares[0] + squares[1] + squares[2]).sqrt() * Interval::point(scale);
    finite(length).then_some(nonnegative(length))
}

fn dot(a: [Interval; 3], b: [Interval; 3]) -> Interval { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn cross(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn finite(value: Interval) -> bool { value.lo().is_finite() && value.hi().is_finite() }
fn nonnegative(value: Interval) -> Interval { Interval::new(value.lo().max(0.0), value.hi().max(0.0)) }
fn hash_float(hash: &mut DomainHasher, value: f64) { hash.update(&value.to_bits().to_le_bytes()); }

#[cfg(test)]
mod tests;
