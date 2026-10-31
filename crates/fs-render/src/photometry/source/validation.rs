//! Numerical enclosures for the original direct point-fixture calculation.
//!
//! The enclosed real function uses the retained binary64 flux, normalized
//! angular coefficients/cutoffs, inverse-rotation quaternion, source/probe
//! positions, receiver normal and gray absorption/scattering coefficients.
//! Constructors' normalization, filtering and input-conversion errors are NOT
//! included. Neither are physical calibration, spectral/model discrepancy or
//! omitted reflected/scattered light. The total error therefore stays unknown.
//!
//! Geometry visibility and physical-grid transport have no enclosure here.
//! Supported source laws are isotropic, uniform cone, integer cosine powers
//! 0..=64 and two-column axisymmetric tables with at most 64 polar knots.
//! Cutoff uncertainty encloses both branches; it never chooses a nominal side
//! and labels it exact. The native nominal calculation is unchanged.

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

/// Version of the represented-function and numerical-enclosure interpretation.
pub const DIRECT_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

/// Why the unchanged nominal result has no numerical enclosure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnclosureUnavailable {
    /// Even an apparently clear native mesh ray lacks certified visibility.
    GeometryVisibility,
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
    pub numerical_error_bound_lux: Option<f64>,
    /// Content identity of inputs actually used by this enclosed real function.
    /// It excludes spectra/phase (irrelevant to gray unscattered lux), inactive
    /// source details and all constructor history. It is not a full scene hash
    /// or authenticated calibration provenance.
    pub numerical_model_identity: Option<ContentHash>,
    /// Explicit absence reason, mutually exclusive with a numerical enclosure.
    pub unavailable: Option<EnclosureUnavailable>,
}

impl EnclosedIlluminance {
    /// Header for the direct numerical-accounting CSV, distinct from path CSV.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(b"index,direct_lux,numerical_lower_lux,numerical_upper_lux,numerical_error_bound_lux,numerical_status,numerical_model_identity,source_construction_error_bound_lux,calibration_error_bound_lux,physical_model_error_bound_lux,total_error_bound_lux,enclosure_version\n")
    }

    /// Stream round-trip numeric values with literal unknown physical bounds.
    /// No I/O side effect occurs until the caller supplies a writer. Failure
    /// propagates; atomic file publication remains the caller's responsibility.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        write!(writer, "{index},{},", self.estimate.lux)?;
        match self.numerical_error_bound_lux {
            Some(error) => write!(writer, "{},{},{error},enclosure,",
                self.estimate.numerical.lo, self.estimate.numerical.hi)?,
            None => {
                let reason = match self.unavailable {
                    Some(EnclosureUnavailable::GeometryVisibility) => "geometry-visibility",
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
        writeln!(writer, ",unknown,unknown,unknown,unknown,{DIRECT_ENCLOSURE_SEMANTICS_VERSION}")
    }
}

#[derive(Clone, Copy)]
enum Prepared {
    Available(ContentHash),
    Unavailable(EnclosureUnavailable),
}

impl PhotometricScene<'_> {
    /// Evaluate the original direct lux plus a supported arithmetic enclosure.
    /// Unsupported enclosure physics does not reject an otherwise valid nominal
    /// calculation: inspect `unavailable`. Cancellation and nominal transport
    /// refusals still return errors. No total physical-accuracy claim is made.
    pub fn illuminance_with_enclosure(
        &self,
        cx: &Cx<'_>,
        probe: &IlluminanceProbe,
    ) -> Result<EnclosedIlluminance, PhotometryError> {
        evaluate(self, cx, probe, prepare(self, cx)?)
    }
}

/// Progressive direct calculation with enclosure and provenance on each value.
/// The original scene/probes stay borrowed. Static numerical-model identity is
/// prepared once; no source-table hashing is repeated at every grid point.
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
    pub fn new(
        scene: &'a PhotometricScene<'scene>,
        probes: &'a [IlluminanceProbe],
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Probes(probes), probes.len(), cx)
    }

    /// Admit the native point-sampled calculation plane without a probe copy.
    pub fn for_grid(
        scene: &'a PhotometricScene<'scene>,
        grid: CalculationGrid,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Grid(grid), grid.len(), cx)
    }

    fn begin(
        scene: &'a PhotometricScene<'scene>,
        receivers: Receivers<'a>,
        count: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if count > scene.options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        let prepared = prepare(scene, cx)?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| PhotometryError::AllocationRefused)?;
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
            let value = evaluate(self.scene, cx, &probe, self.prepared)?;
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
    if !scene.primitives.is_empty() {
        return Ok(Prepared::Unavailable(EnclosureUnavailable::GeometryVisibility));
    }
    let Some(medium) = scene.medium.homogeneous() else {
        return Ok(Prepared::Unavailable(EnclosureUnavailable::PhysicalGridTransport));
    };
    let mut hash = DomainHasher::new("org.frankensim.photometry.direct-enclosure.model.v1");
    hash.update(&DIRECT_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&fs_ivl::INTERVAL_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&fs_ivl::decay::DECAY_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&crate::photometry::PHOTOMETRY_SEMANTICS_VERSION.to_le_bytes());
    for value in [medium.absorption_per_m(), medium.scattering_per_m(),
        scene.options.ray_epsilon_m, scene.options.max_distance_m] {
        hash_float(&mut hash, value);
    }
    hash.update(&(scene.sources.len() as u64).to_le_bytes());
    for (index, source) in scene.sources.iter().enumerate() {
        cx.checkpoint()?;
        hash_float(&mut hash, source.emitted_flux_lm);
        if source.emitted_flux_lm == 0.0 { continue; }
        for value in source.transform.translation_m().into_iter()
            .chain(source.world_to_local_rotation.rotation_xyzw()) {
            hash_float(&mut hash, value);
        }
        hash_float(&mut hash, source.distribution.peak_per_sr);
        match &source.distribution.law {
            AngularLaw::Isotropic => { hash.update(&[0]); },
            AngularLaw::Cone { cosine_cutoff } => {
                hash.update(&[1]);
                hash_float(&mut hash, *cosine_cutoff);
            }
            AngularLaw::CosinePower { exponent, cosine_cutoff }
                if *exponent <= 64.0 && exponent.fract() == 0.0 => {
                hash.update(&[2]);
                hash_float(&mut hash, *exponent);
                hash_float(&mut hash, *cosine_cutoff);
            }
            AngularLaw::Table(table) if table.phi.len() == 2 && table.cosine.len() <= 64 => {
                hash.update(&[3]);
                hash.update(&(table.cosine.len() as u64).to_le_bytes());
                for (row, mu) in table.cosine.iter().enumerate() {
                    hash_float(&mut hash, *mu);
                    hash_float(&mut hash, table.values[2 * row]);
                }
            }
            _ => return Ok(Prepared::Unavailable(EnclosureUnavailable::SourceLaw { source: index })),
        }
    }
    cx.checkpoint()?;
    Ok(Prepared::Available(hash.finalize()))
}

fn evaluate(
    scene: &PhotometricScene<'_>, cx: &Cx<'_>, probe: &IlluminanceProbe, prepared: Prepared,
) -> Result<EnclosedIlluminance, PhotometryError> {
    let mut estimate = scene.illuminance(cx, probe)?;
    cx.checkpoint()?;
    let no_claim = |reason| EnclosedIlluminance {
        estimate, numerical_error_bound_lux: None,
        numerical_model_identity: None, unavailable: Some(reason),
    };
    let model = match prepared {
        Prepared::Unavailable(reason) => return Ok(no_claim(reason)),
        Prepared::Available(model) => model,
    };
    let Some(enclosure) = direct_bound(scene, cx, probe)? else {
        return Ok(no_claim(EnclosureUnavailable::ArithmeticRange));
    };
    // Also contain the native nominal value, without using it as the oracle.
    let enclosure = enclosure.hull(Interval::point(estimate.lux));
    let error = if enclosure.lo() == enclosure.hi() { 0.0 } else {
        (Interval::point(estimate.lux) - Interval::point(enclosure.lo())).hi()
            .max((Interval::point(enclosure.hi()) - Interval::point(estimate.lux)).hi())
    };
    if !error.is_finite() { return Ok(no_claim(EnclosureUnavailable::ArithmeticRange)); }
    estimate.numerical = NumericalCertificate::enclosure(enclosure.lo(), enclosure.hi());
    let mut hash = DomainHasher::new("org.frankensim.photometry.direct-enclosure.probe.v1");
    hash.update(&model.0);
    for value in [probe.position_m.x, probe.position_m.y, probe.position_m.z,
        probe.normal.x, probe.normal.y, probe.normal.z] { hash_float(&mut hash, value); }
    cx.checkpoint()?;
    Ok(EnclosedIlluminance {
        estimate, numerical_error_bound_lux: Some(error),
        numerical_model_identity: Some(hash.finalize()), unavailable: None,
    })
}

fn direct_bound(
    scene: &PhotometricScene<'_>, cx: &Cx<'_>, probe: &IlluminanceProbe,
) -> Result<Option<Interval>, PhotometryError> {
    let Some(medium) = scene.medium.homogeneous() else { return Ok(None); };
    let extinction = nonnegative(Interval::point(medium.absorption_per_m())
        + Interval::point(medium.scattering_per_m()));
    let vacuum = medium.extinction_per_m() == 0.0;
    let mut sum = Interval::point(0.0);
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
        let transmission = if vacuum { Interval::point(1.0) } else {
            let Some(value) = exp_neg(nonnegative(extinction * distance)) else { return Ok(None); };
            value
        };
        let contribution = nonnegative(intensity / distance / distance * cosine * transmission);
        if !finite(contribution) { return Ok(None); }
        sum = nonnegative(sum + contribution);
        if !finite(sum) { return Ok(None); }
    }
    cx.checkpoint()?;
    Ok(Some(sum))
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
fn hash_float(hash: &mut DomainHasher, value: f64) {
    hash.update(&value.to_bits().to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{AngularDistribution, SourceSpectrum};
    use crate::instances::RigidTransform;
    use crate::photometry::TransportOptions;
    use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
    use fs_evidence::NumericalKind;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use std::f64::consts::{FRAC_PI_2, PI};

    fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena,
                StreamKey { seed: 173, kernel_id: 0x5048_454e, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            operation(&cx, &gate)
        })
    }

    fn probe(x: f64, y: f64) -> IlluminanceProbe {
        IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
    }

    fn source(flux: f64, height: f64, law: AngularDistribution) -> PhotometricSource {
        PhotometricSource::try_new(
            RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 0.0, height]).unwrap(),
            flux, 1.0, law, SourceSpectrum::equal_energy()).unwrap()
    }

    fn contains(report: EnclosedIlluminance, reference: f64) {
        assert_eq!(report.unavailable, None);
        assert_eq!(report.estimate.numerical.kind, NumericalKind::Enclosure);
        assert!(report.estimate.numerical.lo <= reference, "{report:?}, reference={reference}");
        assert!(report.estimate.numerical.hi >= reference, "{report:?}, reference={reference}");
        assert!((report.estimate.lux-reference).abs() <= report.numerical_error_bound_lux.unwrap());
        assert!(report.estimate.total_error_bound_lux.is_none());
        assert!(report.numerical_model_identity.is_some());
    }

    #[test]
    fn cie171_diffuse_point_case_has_an_enclosure_not_a_total_accuracy_claim() {
        // Same publicly described CIE 171 section 5.2 analytical subcase as
        // tests/photometric_path/cie171_analytical.rs. No reference table copied.
        // I=1000*cos(theta), h=3: E=9000/(9+x*x+y*y)^2.
        with_cx(|cx, _| {
            let sources = [source(1000.0*PI, 3.0,
                AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap())];
            let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
                TransportOptions::default(), cx).unwrap();
            for [x,y] in [[0.0,0.0],[0.5,0.0],[1.0,0.0],[1.5,0.0],[0.5,0.5],
                [1.0,0.5],[1.5,0.5],[1.0,1.0],[1.5,1.0],[1.5,1.5]] {
                let p = probe(x,y);
                let report = scene.illuminance_with_enclosure(cx, &p).unwrap();
                let squared = 9.0+x*x+y*y;
                contains(report, 9000.0/(squared*squared));
                assert_eq!(report.estimate.lux.to_bits(), scene.illuminance(cx, &p).unwrap().lux.to_bits());
                assert!(report.numerical_error_bound_lux.unwrap() < 1e-9);
            }
        });
    }

    #[test]
    fn gray_haze_encloses_absorption_plus_out_scattering_and_positive_underflow() {
        with_cx(|cx, _| {
            let sources = [source(1000.0*PI, 5.0, AngularDistribution::isotropic())];
            for (a,s) in [(0.2,0.0),(0.0,0.2),(0.07,0.13)] {
                let medium = HomogeneousMedium::try_new(a,s,PhaseFunction::HenyeyGreenstein {g:0.7}).unwrap();
                let scene = PhotometricScene::try_new(&[], &sources, medium,
                    TransportOptions::default(), cx).unwrap();
                let report = scene.illuminance_with_enclosure(cx, &probe(0.0,0.0)).unwrap();
                contains(report, 3.678_794_411_714_423_3); // independent 10/e
                assert!(report.numerical_error_bound_lux.unwrap() < 1e-10);
            }
            let scene = PhotometricScene::try_new(&[], &sources,
                HomogeneousMedium::try_new(200.0,0.0,PhaseFunction::Isotropic).unwrap(),
                TransportOptions::default(), cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx, &probe(0.0,0.0)).unwrap();
            assert_eq!(report.unavailable, None);
            assert!(report.estimate.numerical.hi > 0.0);
            assert!(report.estimate.numerical.hi < f64::MIN_POSITIVE);
            assert!(report.numerical_error_bound_lux.unwrap() > 0.0);
            assert!(report.estimate.total_error_bound_lux.is_none());
        });
    }

    #[test]
    fn axisymmetric_beam_tables_and_uncertain_cutoffs_keep_their_support() {
        with_cx(|cx, _| {
            let sources = [source(1000.0, 3.0, AngularDistribution::beam_field(0.2,0.4,0.6).unwrap())];
            let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
                TransportOptions::default(), cx).unwrap();
            for x in [0.0,0.1,0.3,0.5,0.7,1.0] {
                let p = probe(x,0.0);
                let report = scene.illuminance_with_enclosure(cx, &p).unwrap();
                // A consistency check of retained table interpolation, separate
                // from the independent CIE and airlight reference checks.
                contains(report, scene.illuminance(cx, &p).unwrap().lux);
                assert!(report.numerical_error_bound_lux.unwrap() < 1e-8);
            }
        });
        let peak = Interval::point(20.0);
        assert_eq!(cutoff(peak, Interval::new(0.5_f64.next_down(),0.5_f64.next_up()),0.5),
            Interval::new(0.0,20.0));
        assert_eq!(cutoff(peak,Interval::point(0.5),0.5),peak);
        assert_eq!(cutoff(peak,Interval::point(0.5_f64.next_down()),0.5),Interval::point(0.0));
    }

    #[test]
    fn unsupported_visibility_and_angular_laws_retain_nominal_no_claim() {
        use crate::charts::TriMesh;
        use crate::spectral::SpectralReflectance;
        use crate::tracer::{Material, Primitive, Shape};
        with_cx(|cx, _| {
            let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
            let geometry = [Primitive {
                shape: Shape::Mesh(TriMesh::new(vec![[-1.0,-1.0,1.0],[1.0,-1.0,1.0],[0.0,1.0,1.0]],vec![[0,1,2]])),
                material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() },
                emission: None,
            }];
            let scene = PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            assert_eq!(report.estimate.lux,0.0);
            assert_eq!(report.unavailable,Some(EnclosureUnavailable::GeometryVisibility));
            assert_eq!(report.estimate.numerical.kind,NumericalKind::NoClaim);
            assert!(report.numerical_error_bound_lux.is_none());
            assert!(report.numerical_model_identity.is_none());
            let sources = [source(1000.0,3.0,AngularDistribution::cosine_power(0.5,FRAC_PI_2).unwrap())];
            let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            assert!(report.estimate.lux > 0.0);
            assert_eq!(report.unavailable,Some(EnclosureUnavailable::SourceLaw {source:0}));
            assert!(report.estimate.total_error_bound_lux.is_none());
        });
    }

    #[test]
    fn progressive_grid_cancel_resume_preserves_values_bounds_and_identities() {
        with_cx(|cx,gate| {
            let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
            let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
            let grid = CalculationGrid::try_new(Point3::new(-1.0,-1.0,0.0),Vec3::new(0.5,0.0,0.0),Vec3::new(0.0,0.5,0.0),5,5).unwrap();
            let mut full = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
            full.advance(cx,usize::MAX).unwrap();
            let mut partial = EnclosedIlluminanceJob::for_grid(&scene,grid,cx).unwrap();
            partial.advance(cx,3).unwrap();
            assert_eq!(partial.results(),&full.results()[..3]);
            gate.request();
            assert!(matches!(partial.advance(cx,5),Err(PhotometryError::Cancelled)));
            assert_eq!(partial.results(),&full.results()[..3]);
            with_cx(|resume,_| {
                while !partial.is_complete() { partial.advance(resume,7).unwrap(); }
            });
            assert_eq!(partial.results(),full.results());
        });
    }

    #[test]
    fn source_and_receiver_changes_change_the_numerical_model_identity() {
        with_cx(|cx,_| {
            let a = [source(1000.0,3.0,AngularDistribution::isotropic())];
            let b = [source(2000.0,3.0,AngularDistribution::isotropic())];
            let sa = PhotometricScene::try_new(&[], &a, HomogeneousMedium::vacuum(),
                TransportOptions::default(), cx).unwrap();
            let sb = PhotometricScene::try_new(&[], &b, HomogeneousMedium::vacuum(),
                TransportOptions::default(), cx).unwrap();
            let ea = sa.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            let eb = sb.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            let ec = sa.illuminance_with_enclosure(cx,&probe(1.0,0.0)).unwrap();
            assert_ne!(ea.numerical_model_identity,eb.numerical_model_identity);
            assert_ne!(ea.numerical_model_identity,ec.numerical_model_identity);
            assert_eq!(ea,sa.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap());
        });
    }

    #[test]
    fn csv_keeps_numerical_and_unknown_physical_errors_separate() {
        with_cx(|cx,_| {
            let sources = [source(1000.0,3.0,AngularDistribution::isotropic())];
            let scene = PhotometricScene::try_new(&[],&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
            let report = scene.illuminance_with_enclosure(cx,&probe(0.0,0.0)).unwrap();
            let mut bytes = Vec::new();
            EnclosedIlluminance::write_csv_header(&mut bytes).unwrap();
            report.write_csv_row(0,&mut bytes).unwrap();
            let csv = String::from_utf8(bytes).unwrap();
            let rows: Vec<_> = csv.lines().collect();
            let columns: Vec<_> = rows[1].split(',').collect();
            assert_eq!(columns.len(),12);
            assert_eq!(rows[0].split(',').count(),12);
            assert_eq!(columns[5],"enclosure");
            assert_eq!(columns[6].len(),64);
            assert_eq!(&columns[7..11],&["unknown";4]);
            assert_eq!(columns[1].parse::<f64>().unwrap().to_bits(),report.estimate.lux.to_bits());
            struct Refuse;
            impl Write for Refuse {
                fn write(&mut self,_: &[u8])->io::Result<usize> { Err(io::Error::other("test sink")) }
                fn flush(&mut self)->io::Result<()> { Ok(()) }
            }
            assert!(report.write_csv_row(0,&mut Refuse).is_err());
        });
    }
}
