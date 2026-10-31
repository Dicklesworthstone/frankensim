//! Differential ray tubes through smooth, isotropic refractive field segments.
//!
//! Two derivatives with respect to orthonormal launch angles travel alongside
//! the central ray. The field explicitly supplies its analytic spatial Hessian;
//! a zero local gradient never implies a zero Hessian. No finite-difference
//! fallback or representation-dependent smoothing is inferred.
//!
//! Writing `g = grad(n)/n` and `K = Hess(n)/n - g g^T`, the first variations
//! of the parent's arc-length equations are `X' = V` and
//! `V' = (I-tt^T) K X - (t.g) V - t (g.V)`. Internal and final unit-direction
//! projections are differentiated too. An augmented Dormand–Prince 5(4) step
//! controls the central optics and both derivative columns together.
//!
//! Endpoint-plane projection accounts for the neighboring rays' different
//! arrival lengths. The resulting `dA/dOmega` converts point-source intensity
//! to irradiance, including the central ray's absorption exactly once.
//! It describes a plane THROUGH the integrated endpoint; it does not search
//! for an intersection or a connecting ray. Local estimates are not global
//! derivative certificates. Caustics and poorly conditioned/grazing projections
//! refuse an irradiance result. Validated continuation preserves the original
//! angular coordinates across separately differentiated interface maps. This
//! geometric owner does not supply Fresnel power weights; a physical interface
//! consumer must retain those separately. Visibility changes, field-parameter
//! gradients, diffraction and scalar-camera/BDPT MIS remain separate problems.
//!
//! Continuous-ray differentiation and the need for analytic field Hessians:
//! Teh et al., *Adjoint Nonlinear Ray Tracing* (2022),
//! <https://doi.org/10.1145/3528223.3530077>. This owner uses forward variations
//! in physical arc length, rather than that paper's canonical-parameter adjoint.

use super::*;
use fs_ivl::Interval;

pub mod continuation;
pub mod interface;
pub use continuation::{
    BundleContinuation, propagate_seeded_bundle, propagate_seeded_bundle_with_control,
};

/// Numerical semantics of smooth-field angular ray bundles.
pub const RAY_BUNDLE_SEMANTICS_VERSION: u32 = 2;

/// Analytic, symmetric spatial Hessian of phase index, in inverse square metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexHessian {
    diagonal: [f64; 3],
    off_diagonal: [f64; 3],
}

impl IndexHessian {
    /// Exact zero for a spatially homogeneous or affine phase-index law.
    pub const ZERO: Self = Self {
        diagonal: [0.0; 3],
        off_diagonal: [0.0; 3],
    };

    /// Diagonal entries are `(xx,yy,zz)`; off-diagonal entries are `(xy,xz,yz)`.
    /// Symmetry is structural; finite values are required. The source remains
    /// responsible for these being derivatives of its declared phase law.
    pub fn try_new(diagonal: [f64; 3], off_diagonal: [f64; 3]) -> Result<Self, BundleError> {
        if diagonal
            .into_iter()
            .chain(off_diagonal)
            .any(|x| !x.is_finite())
        {
            return Err(BundleError::InvalidDifferential);
        }
        Ok(Self {
            diagonal,
            off_diagonal,
        })
    }

    pub const fn diagonal(self) -> [f64; 3] {
        self.diagonal
    }
    pub const fn off_diagonal(self) -> [f64; 3] {
        self.off_diagonal
    }

    fn apply(self, v: Vec3) -> Vec3 {
        let [xx, yy, zz] = self.diagonal;
        let [xy, xz, yz] = self.off_diagonal;
        Vec3::new(
            xx * v.x + xy * v.y + xz * v.z,
            xy * v.x + yy * v.y + yz * v.z,
            xz * v.x + yz * v.y + zz * v.z,
        )
    }
}

/// One consistent optical sample and its analytic phase-index Hessian.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexJet {
    optical: GradientIndexSample,
    hessian: IndexHessian,
}

impl IndexJet {
    /// Both inputs already carry their finite-value admission.
    pub const fn new(optical: GradientIndexSample, hessian: IndexHessian) -> Self {
        Self { optical, hessian }
    }
    pub const fn optical(self) -> GradientIndexSample {
        self.optical
    }
    pub const fn hessian(self) -> IndexHessian {
        self.hessian
    }
}

/// Opt-in second spatial derivative capability for a frozen smooth field.
/// The jet's optical sample must equal the ordinary field sample at the same
/// position/wavelength. One query returns both, avoiding duplicate evaluation.
pub trait SecondOrderIndexField: GradientIndexField {
    fn sample_second_order(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError>;
}

impl SecondOrderIndexField for HomogeneousIndexField<'_> {
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        Ok(IndexJet::new(
            self.sample(p, wavelength_nm)?,
            IndexHessian::ZERO,
        ))
    }
}
impl SecondOrderIndexField for AffineIndexField<'_> {
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        Ok(IndexJet::new(
            self.sample(p, wavelength_nm)?,
            IndexHessian::ZERO,
        ))
    }
}

/// Bounded quadratic phase index on a native spectral material:
/// `n(r,lambda) = n_reference(lambda) + g.d + 0.5*d.H.d`, `d=r-reference`.
/// The authored spatial offset is independent of wavelength, so it contributes
/// the same offset to group index. Absorption retains the material's law.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticIndexField<'a> {
    glass: &'a DielectricGlass,
    reference: Point3,
    gradient: Vec3,
    hessian: IndexHessian,
    bounds: Aabb,
}

impl<'a> QuadraticIndexField<'a> {
    /// An outward interval enclosure admits the complete box, including its
    /// interior extrema, over every supported wavelength. Conservative refusal
    /// is possible for a box whose enclosure is wider than its exact range.
    pub fn try_new(
        glass: &'a DielectricGlass,
        reference: Point3,
        gradient: Vec3,
        hessian: IndexHessian,
        bounds: Aabb,
    ) -> Result<Self, BundleError> {
        if !finite_point(reference)
            || !finite_vec(gradient)
            || !finite_point(bounds.min)
            || !finite_point(bounds.max)
            || bounds.min.x >= bounds.max.x
            || bounds.min.y >= bounds.max.y
            || bounds.min.z >= bounds.max.z
        {
            return Err(BundleError::InvalidDifferential);
        }
        let candidate = Self {
            glass,
            reference,
            gradient,
            hessian,
            bounds,
        };
        let offset = candidate.offset_enclosure();
        let check = |wavelength| -> Result<(), BundleError> {
            let n = glass
                .ior()
                .eval(wavelength)
                .map_err(GradientIndexError::from)?;
            let range = if gradient == Vec3::new(0.0, 0.0, 0.0) && hessian == IndexHessian::ZERO {
                Interval::point(n)
            } else {
                Interval::point(n) + offset
            };
            if !range.lo().is_finite()
                || !range.hi().is_finite()
                || range.lo() < 1.0
                || range.hi() > 3.0
            {
                return Err(BundleError::InvalidDifferential);
            }
            Ok(())
        };
        check(LAMBDA_MIN)?;
        check(LAMBDA_MAX)?;
        if let DielectricIor::Tabulated(table) = glass.ior() {
            for sample in table.samples() {
                check(sample.wavelength_nm())?;
            }
        }
        Ok(candidate)
    }

    fn offset_enclosure(&self) -> Interval {
        let coordinates = [
            (
                self.bounds.min.x,
                self.bounds.max.x,
                self.reference.x,
                self.gradient.x,
            ),
            (
                self.bounds.min.y,
                self.bounds.max.y,
                self.reference.y,
                self.gradient.y,
            ),
            (
                self.bounds.min.z,
                self.bounds.max.z,
                self.reference.z,
                self.gradient.z,
            ),
        ];
        let d = coordinates.map(|(lo, hi, r, _)| Interval::new(lo, hi) - Interval::point(r));
        let mut offset = Interval::point(0.0);
        for i in 0..3 {
            if coordinates[i].3 != 0.0 {
                offset = offset + Interval::point(coordinates[i].3) * d[i];
            }
            if self.hessian.diagonal[i] != 0.0 {
                let a = d[i].abs();
                offset = offset
                    + Interval::point(0.5) * Interval::point(self.hessian.diagonal[i]) * (a * a);
            }
        }
        for (index, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            if self.hessian.off_diagonal[index] != 0.0 {
                offset = offset + Interval::point(self.hessian.off_diagonal[index]) * d[i] * d[j];
            }
        }
        offset
    }
}

impl GradientIndexField for QuadraticIndexField<'_> {
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if !finite_point(position) || !self.bounds.contains(position) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let d = position.delta_from(self.reference);
        let hd = self.hessian.apply(d);
        let gradient = Vec3::new(
            self.gradient.x + hd.x,
            self.gradient.y + hd.y,
            self.gradient.z + hd.z,
        );
        let offset = self.gradient.dot(d) + 0.5 * d.dot(hd);
        material_sample(self.glass, wavelength_nm, offset, gradient)
    }
}
impl SecondOrderIndexField for QuadraticIndexField<'_> {
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        Ok(IndexJet::new(self.sample(p, wavelength_nm)?, self.hessian))
    }
}

/// Refusal from ray-bundle propagation or a local receiving-plane projection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BundleError {
    Propagation(GradientIndexError),
    InvalidDifferential,
    InvalidOptions,
    InvalidProjection,
    GrazingProjection,
    SingularProjection,
    InvalidIntensity,
    /// A mapped endpoint changed position or violated differential invariants.
    InvalidContinuation,
    /// The next field does not match the optical state retained at the endpoint.
    InitialSampleMismatch,
    NumericalRange,
}
impl From<GradientIndexError> for BundleError {
    fn from(error: GradientIndexError) -> Self {
        Self::Propagation(error)
    }
}
impl core::fmt::Display for BundleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ray bundle: {self:?}")
    }
}
impl core::error::Error for BundleError {}

/// Central optical controls and additional absolute local derivative tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleOptions {
    pub integration: GradientIndexOptions,
    /// Metres per radian of orthonormal launch-angle perturbation.
    pub position_derivative_tolerance_m: f64,
    /// Unit-direction change per radian of launch-angle perturbation.
    pub direction_derivative_tolerance: f64,
}
impl Default for BundleOptions {
    fn default() -> Self {
        Self {
            integration: GradientIndexOptions::default(),
            position_derivative_tolerance_m: 1e-9,
            direction_derivative_tolerance: 1e-10,
        }
    }
}
impl BundleOptions {
    fn validate(self) -> Result<(), BundleError> {
        self.integration.validate()?;
        if [
            self.position_derivative_tolerance_m,
            self.direction_derivative_tolerance,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
        {
            return Err(BundleError::InvalidOptions);
        }
        Ok(())
    }
}

/// One first variation with respect to a dimensionless launch angle in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularDerivative {
    pub position_m_per_rad: Vec3,
    pub direction_per_rad: Vec3,
}

/// Smooth-field angular bundle, retaining its immutable physical launch.
/// The two local coordinates are `normalize(t0 + alpha*e0 + beta*e1)`;
/// at zero their area measure is exactly `dOmega = d_alpha d_beta`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayBundlePropagation {
    initial_ray: Ray,
    wavelength_nm: f64,
    launch_basis: [Vec3; 2],
    propagation: GradientIndexPropagation,
    derivatives: [AngularDerivative; 2],
    interface_events: u32,
}
impl RayBundlePropagation {
    pub const fn initial_ray(&self) -> Ray {
        self.initial_ray
    }
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    pub const fn launch_basis(&self) -> [Vec3; 2] {
        self.launch_basis
    }
    pub const fn propagation(&self) -> &GradientIndexPropagation {
        &self.propagation
    }
    pub const fn angular_derivatives(&self) -> [AngularDerivative; 2] {
        self.derivatives
    }
    /// Number of explicitly mapped discrete interfaces in this geometric path.
    /// Their Fresnel power/polarization weights are owned by the caller.
    pub const fn interface_events(&self) -> u32 {
        self.interface_events
    }

    /// Derivatives of arrival at the plane through the EXISTING endpoint.
    /// No position is advanced/projected and no intersection is searched.
    pub fn project_at_endpoint(
        &self,
        normal: Vec3,
        options: ProjectionOptions,
    ) -> Result<PlanarBundleProjection, BundleError> {
        options.validate()?;
        let normal = unit(normal).ok_or(BundleError::InvalidProjection)?;
        let tangent = self.propagation.ray.dir;
        let cosine = normal.dot(tangent);
        if cosine.abs() < options.min_incidence_cosine {
            return Err(BundleError::GrazingProjection);
        }
        let sample = self.propagation.final_sample;
        let g = sample.phase_gradient_per_m.scale(1.0 / sample.phase_index);
        let curvature = sub(g, tangent.scale(tangent.dot(g)));
        let mut arrival = self.derivatives;
        let mut arc_derivatives = [0.0; 2];
        for i in 0..2 {
            let ds = -normal.dot(self.derivatives[i].position_m_per_rad) / cosine;
            arc_derivatives[i] = ds;
            let position = self.derivatives[i].position_m_per_rad;
            let direction = self.derivatives[i].direction_per_rad;
            arrival[i] = AngularDerivative {
                position_m_per_rad: sub(position, tangent.scale(-ds)),
                direction_per_rad: sub(direction, curvature.scale(-ds)),
            };
            if !ds.is_finite()
                || !finite_vec(arrival[i].position_m_per_rad)
                || !finite_vec(arrival[i].direction_per_rad)
            {
                return Err(BundleError::NumericalRange);
            }
        }
        let a = arrival[0].position_m_per_rad;
        let b = arrival[1].position_m_per_rad;
        let scale = vec_norm(a).max(vec_norm(b));
        if !scale.is_finite() {
            return Err(BundleError::NumericalRange);
        }
        if scale == 0.0 {
            return Err(BundleError::SingularProjection);
        }
        // Division keeps a subnormal footprint scale from overflowing an
        // intermediate reciprocal before the explicit caustic admission.
        let a = Vec3::new(a.x / scale, a.y / scale, a.z / scale);
        let b = Vec3::new(b.x / scale, b.y / scale, b.z / scale);
        let area_scaled = normal.dot(cross(a, b));
        let aa = a.dot(a);
        let bb = b.dot(b);
        let ab = a.dot(b);
        let largest_eigenvalue = 0.5 * (aa + bb + det::hypot(aa - bb, 2.0 * ab));
        let condition = largest_eigenvalue / area_scaled.abs();
        let signed_area = (scale * area_scaled) * scale;
        if !signed_area.is_finite() {
            return Err(BundleError::NumericalRange);
        }
        if signed_area.abs() < options.min_area_m2_per_sr
            || !condition.is_finite()
            || condition > options.max_condition_number
        {
            return Err(BundleError::SingularProjection);
        }
        Ok(PlanarBundleProjection {
            bundle: *self,
            normal,
            derivatives: arrival,
            arc_derivatives,
            signed_area_m2_per_sr: signed_area,
            condition_number: condition,
        })
    }
}

/// Explicit resolution limits on a local receiving-plane footprint.
/// These are admission thresholds, not inferred global integration-error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProjectionOptions {
    pub min_incidence_cosine: f64,
    pub min_area_m2_per_sr: f64,
    pub max_condition_number: f64,
}
impl Default for ProjectionOptions {
    fn default() -> Self {
        Self {
            min_incidence_cosine: 1e-6,
            min_area_m2_per_sr: 1e-18,
            max_condition_number: 1e10,
        }
    }
}
impl ProjectionOptions {
    fn validate(self) -> Result<(), BundleError> {
        if !self.min_incidence_cosine.is_finite()
            || self.min_incidence_cosine <= 0.0
            || self.min_incidence_cosine > 1.0
            || !self.min_area_m2_per_sr.is_finite()
            || self.min_area_m2_per_sr <= 0.0
            || !self.max_condition_number.is_finite()
            || self.max_condition_number < 1.0
        {
            return Err(BundleError::InvalidProjection);
        }
        Ok(())
    }
}

/// A local differential receiving-plane footprint with nonzero admitted area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarBundleProjection {
    bundle: RayBundlePropagation,
    normal: Vec3,
    derivatives: [AngularDerivative; 2],
    arc_derivatives: [f64; 2],
    signed_area_m2_per_sr: f64,
    condition_number: f64,
}
impl PlanarBundleProjection {
    pub const fn bundle(&self) -> &RayBundlePropagation {
        &self.bundle
    }
    pub const fn normal(&self) -> Vec3 {
        self.normal
    }
    pub const fn angular_derivatives(&self) -> [AngularDerivative; 2] {
        self.derivatives
    }
    pub const fn arc_length_derivatives_m_per_rad(&self) -> [f64; 2] {
        self.arc_derivatives
    }
    pub fn area_per_solid_angle_m2(&self) -> f64 {
        self.signed_area_m2_per_sr.abs()
    }
    pub const fn signed_area_per_solid_angle_m2(&self) -> f64 {
        self.signed_area_m2_per_sr
    }
    pub const fn condition_number(&self) -> f64 {
        self.condition_number
    }

    /// Incident point-source irradiance from radiant intensity [W/sr]. Spectral
    /// intensity [W/(sr nm)] instead produces spectral irradiance [W/(m² nm)].
    /// Source intensity must belong to the retained launch and wavelength.
    /// Includes absorption once; no extra cosine or n² factor is required.
    /// This geometric utility excludes Fresnel weights for any retained
    /// interface events; physical interface consumers must apply those once.
    /// Multiple arrival branches add separately; this is a local ray tube,
    /// not a finite-pixel average or a diffraction-regularized caustic.
    pub fn point_source_irradiance(&self, intensity: f64) -> Result<f64, BundleError> {
        if !intensity.is_finite() || intensity < 0.0 {
            return Err(BundleError::InvalidIntensity);
        }
        let result = (intensity * self.bundle.propagation.power_transmittance())
            / self.area_per_solid_angle_m2();
        if result.is_finite() && (result > 0.0 || intensity == 0.0) {
            return Ok(result);
        }
        // An underflowed transmission/product can still give representable
        // irradiance after source intensity and footprint area are combined.
        // Preserve the ordinary arithmetic above for the common case.
        let logarithm = det::ln(intensity)
            - self.bundle.propagation.optical_depth
            - det::ln(self.area_per_solid_angle_m2());
        let result = det::exp(logarithm);
        if !result.is_finite() {
            return Err(BundleError::NumericalRange);
        }
        Ok(result)
    }
}

/// Propagate the central optics and two launch-angle derivatives to fixed arc length.
pub fn propagate_bundle<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    distance_m: f64,
    options: BundleOptions,
) -> Result<RayBundlePropagation, BundleError> {
    propagate_bundle_with_control(
        field,
        ray,
        transverse_axis,
        wavelength_nm,
        distance_m,
        options,
        || true,
    )
}

/// Cancellation is polled before the first field query and every attempted
/// step, including consecutive rejections. Failure publishes no partial bundle.
pub fn propagate_bundle_with_control<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    distance_m: f64,
    options: BundleOptions,
    mut control: impl FnMut() -> bool,
) -> Result<RayBundlePropagation, BundleError> {
    options.validate()?;
    if !distance_m.is_finite() || distance_m < 0.0 {
        return Err(BundleError::InvalidOptions);
    }
    if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
        return Err(GradientIndexError::Dielectric(DielectricError::InvalidWavelength).into());
    }
    if !finite_point(ray.origin) {
        return Err(GradientIndexError::InvalidRay.into());
    }
    let direction = unit(ray.dir).ok_or(GradientIndexError::InvalidRay)?;
    let axis =
        transverse(direction, transverse_axis).ok_or(GradientIndexError::InvalidTransverseAxis)?;
    let vertical = unit(cross(direction, axis)).ok_or(GradientIndexError::InvalidTransverseAxis)?;
    let launch = Ray {
        origin: ray.origin,
        dir: direction,
    };
    let basis = [axis, vertical];
    let columns = basis.map(|e| AngularDerivative {
        position_m_per_rad: Vec3::new(0.0, 0.0, 0.0),
        direction_per_rad: e,
    });
    let (propagation, derivatives) = propagate_segment(
        field,
        BundleSegmentStart {
            ray: launch,
            axis,
            columns,
            expected_sample: None,
        },
        wavelength_nm,
        distance_m,
        options,
        &mut control,
    )?;
    Ok(RayBundlePropagation {
        initial_ray: launch,
        wavelength_nm,
        launch_basis: basis,
        propagation,
        derivatives,
        interface_events: 0,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct BundleSegmentStart {
    ray: Ray,
    axis: Vec3,
    columns: [AngularDerivative; 2],
    expected_sample: Option<GradientIndexSample>,
}

// One solver handles both fresh launch angles and arbitrary validated
// endpoint columns. Its displacement and optical integrals remain local to
// this segment; the continuation owner then accumulates the physical path.
fn propagate_segment<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    start: BundleSegmentStart,
    wavelength_nm: f64,
    distance_m: f64,
    options: BundleOptions,
    control: &mut impl FnMut() -> bool,
) -> Result<(GradientIndexPropagation, [AngularDerivative; 2]), BundleError> {
    options.validate()?;
    if !distance_m.is_finite() || distance_m < 0.0 {
        return Err(BundleError::InvalidOptions);
    }
    let ray = start.ray;
    let direction = ray.dir;
    let axis = start.axis;
    if !control() {
        return Err(GradientIndexError::Cancelled.into());
    }
    let initial = field.sample_second_order(ray.origin, wavelength_nm)?;
    if start
        .expected_sample
        .is_some_and(|sample| sample != initial.optical)
    {
        return Err(BundleError::InitialSampleMismatch);
    }
    if field.is_spatially_uniform() {
        if initial.optical.phase_gradient_per_m != Vec3::new(0.0, 0.0, 0.0)
            || initial.hessian != IndexHessian::ZERO
        {
            return Err(GradientIndexError::InvalidSample.into());
        }
        // Reuse the exact native uniform optics. Only the query above touches
        // the authored field; the internal adapter serves that constant sample.
        struct ConstantSample(GradientIndexSample);
        impl GradientIndexField for ConstantSample {
            fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
                Ok(self.0)
            }
            fn is_spatially_uniform(&self) -> bool {
                true
            }
        }
        let propagation = super::propagate(
            &ConstantSample(initial.optical),
            ray,
            axis,
            wavelength_nm,
            distance_m,
            options.integration,
        )?;
        let derivatives = start.columns.map(|column| AngularDerivative {
            position_m_per_rad: sub(
                column.position_m_per_rad,
                column.direction_per_rad.scale(-distance_m),
            ),
            direction_per_rad: column.direction_per_rad,
        });
        if derivatives
            .iter()
            .any(|column| !finite_vec(column.position_m_per_rad))
        {
            return Err(BundleError::NumericalRange);
        }
        return Ok((propagation, derivatives));
    }
    let track_group = initial.optical.group_index.is_some();
    let mut state = [0.0; 24];
    write_vec(&mut state, 3, direction);
    write_vec(&mut state, 6, axis);
    for (i, column) in start.columns.iter().enumerate() {
        write_vec(&mut state, 12 + 6 * i, column.position_m_per_rad);
        write_vec(&mut state, 15 + 6 * i, column.direction_per_rad);
    }
    let mut jet = initial;
    let mut derivative = bundle_rhs(&state, jet, track_group)?;
    let mut distance = 0.0;
    let mut step = options.integration.initial_step_m.min(distance_m);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut evaluations = 1;
    let mut max_error = 0.0_f64;
    while distance < distance_m {
        if !control() {
            return Err(GradientIndexError::Cancelled.into());
        }
        if accepted + rejected >= options.integration.max_steps {
            return Err(GradientIndexError::StepBudgetExceeded.into());
        }
        let remaining = distance_m - distance;
        step = step.min(remaining);
        if !step.is_finite() || step <= 0.0 || distance + step == distance {
            return Err(BundleError::NumericalRange);
        }
        let trial = bundle_step(
            field,
            ray.origin,
            wavelength_nm,
            &state,
            derivative,
            step,
            track_group,
        )?;
        evaluations += 6;
        let error = bundle_error(&state, &trial, step, options)?;
        let factor = if error == 0.0 {
            5.0
        } else {
            (0.9 * det::pow(error, -0.2)).clamp(0.2, 5.0)
        };
        if error <= 1.0 {
            state = trial.state;
            jet = trial.jet;
            derivative = bundle_rhs(&state, jet, track_group)?;
            distance = if step == remaining {
                distance_m
            } else {
                distance + step
            };
            accepted += 1;
            max_error = max_error.max(error);
            step = (step * factor)
                .min(options.integration.max_step_m)
                .max(options.integration.min_step_m);
        } else {
            rejected += 1;
            if step <= options.integration.min_step_m {
                return Err(GradientIndexError::AccuracyNotReached.into());
            }
            step = (step * factor.min(0.9)).max(options.integration.min_step_m);
        }
    }
    let endpoint = ray.origin.offset(read_vec(&state, 0));
    if !finite_point(endpoint) {
        return Err(BundleError::NumericalRange);
    }
    let propagation = GradientIndexPropagation {
        ray: Ray {
            origin: endpoint,
            dir: read_vec(&state, 3),
        },
        transverse_axis: read_vec(&state, 6),
        initial_sample: initial.optical,
        final_sample: jet.optical,
        distance_m,
        optical_path_length_m: state[9],
        group_optical_path_length_m: track_group.then_some(state[10]),
        optical_depth: state[11],
        accepted_steps: accepted,
        rejected_steps: rejected,
        field_evaluations: evaluations,
        max_normalized_error_estimate: max_error,
    };
    let derivatives = core::array::from_fn(|i| AngularDerivative {
        position_m_per_rad: read_vec(&state, 12 + 6 * i),
        direction_per_rad: read_vec(&state, 15 + 6 * i),
    });
    Ok((propagation, derivatives))
}

type BundleState = [f64; 24];
struct BundleTrial {
    state: BundleState,
    error: BundleState,
    projection_error: f64,
    column_projection_error: [f64; 2],
    jet: IndexJet,
}
fn read_vec(state: &BundleState, i: usize) -> Vec3 {
    Vec3::new(state[i], state[i + 1], state[i + 2])
}
fn write_vec(state: &mut BundleState, i: usize, v: Vec3) {
    state[i] = v.x;
    state[i + 1] = v.y;
    state[i + 2] = v.z;
}
fn optical_state(state: &BundleState) -> State {
    core::array::from_fn(|i| state[i])
}
fn vec_norm(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}

fn bundle_rhs(
    state: &BundleState,
    jet: IndexJet,
    track_group: bool,
) -> Result<BundleState, GradientIndexError> {
    let base = rhs(&optical_state(state), jet.optical, track_group)?;
    let direction = get_vec(&base, 0);
    let raw_norm = vec_norm(read_vec(state, 3));
    if !raw_norm.is_finite() || raw_norm <= 0.0 {
        return Err(GradientIndexError::NumericalRange);
    }
    let inverse_n = 1.0 / jet.optical.phase_index;
    let g = jet.optical.phase_gradient_per_m.scale(inverse_n);
    let mut result = [0.0; 24];
    result[..12].copy_from_slice(&base);
    for i in 0..2 {
        let x = read_vec(state, 12 + 6 * i);
        let raw_v = read_vec(state, 15 + 6 * i);
        let v = sub(raw_v, direction.scale(direction.dot(raw_v))).scale(1.0 / raw_norm);
        let kx = sub(jet.hessian.apply(x).scale(inverse_n), g.scale(g.dot(x)));
        let dv = sub(
            sub(
                sub(kx, direction.scale(direction.dot(kx))),
                v.scale(direction.dot(g)),
            ),
            direction.scale(g.dot(v)),
        );
        write_vec(&mut result, 12 + 6 * i, v);
        write_vec(&mut result, 15 + 6 * i, dv);
    }
    if result.iter().any(|x| !x.is_finite()) {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn bundle_step<F: SecondOrderIndexField + ?Sized>(
    field: &F,
    origin: Point3,
    wavelength: f64,
    state: &BundleState,
    first: BundleState,
    step: f64,
    track_group: bool,
) -> Result<BundleTrial, GradientIndexError> {
    let mut k = [[0.0; 24]; 7];
    k[0] = first;
    let mut stage = *state;
    let mut last = None;
    for (i, row) in A.iter().enumerate() {
        for j in 0..24 {
            let mut increment = 0.0;
            for a in 0..=i {
                increment += row[a] * k[a][j];
            }
            stage[j] = state[j] + step * increment;
        }
        if stage.iter().any(|x| !x.is_finite()) {
            return Err(GradientIndexError::NumericalRange);
        }
        let point = origin.offset(read_vec(&stage, 0));
        if !finite_point(point) {
            return Err(GradientIndexError::NumericalRange);
        }
        let jet = field.sample_second_order(point, wavelength)?;
        k[i + 1] = bundle_rhs(&stage, jet, track_group)?;
        last = Some(jet);
    }
    let mut error = [0.0; 24];
    for j in 0..24 {
        for i in 0..7 {
            error[j] += ERROR[i] * k[i][j];
        }
        error[j] *= step;
    }
    let raw = read_vec(&stage, 3);
    let raw_axis = read_vec(&stage, 6);
    let direction = unit(raw).ok_or(GradientIndexError::NumericalRange)?;
    let axis = transverse(direction, raw_axis).ok_or(GradientIndexError::NumericalRange)?;
    let projection_error = norm(sub(direction, raw)).max(norm(sub(axis, raw_axis)));
    let raw_norm = vec_norm(raw);
    let mut column_projection_error = [0.0; 2];
    for (i, correction) in column_projection_error.iter_mut().enumerate() {
        let v = read_vec(&stage, 15 + 6 * i);
        let projected = sub(v, direction.scale(direction.dot(v))).scale(1.0 / raw_norm);
        *correction = vec_norm(sub(projected, v));
        write_vec(&mut stage, 15 + 6 * i, projected);
    }
    write_vec(&mut stage, 3, direction);
    write_vec(&mut stage, 6, axis);
    if stage[9] < state[9]
        || stage[10] < state[10]
        || stage[11] < state[11]
        || stage.iter().any(|x| !x.is_finite())
    {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(BundleTrial {
        state: stage,
        error,
        projection_error,
        column_projection_error,
        jet: last.ok_or(GradientIndexError::NumericalRange)?,
    })
}

fn bundle_error(
    state: &BundleState,
    trial: &BundleTrial,
    step: f64,
    options: BundleOptions,
) -> Result<f64, GradientIndexError> {
    let base = Trial {
        state: optical_state(&trial.state),
        error: optical_state(&trial.error),
        projection_error: trial.projection_error,
        sample: trial.jet.optical,
    };
    let mut error = normalized_error(&optical_state(state), &base, step, options.integration)?;
    for i in 12..24 {
        let absolute = if (i - 12) % 6 < 3 {
            options.position_derivative_tolerance_m
        } else {
            options.direction_derivative_tolerance
        };
        let scale = absolute
            + options.integration.relative_tolerance * state[i].abs().max(trial.state[i].abs());
        let normalized = trial.error[i].abs() / scale;
        if !scale.is_finite() || scale <= 0.0 || !normalized.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        error = error.max(normalized);
    }
    for i in 0..2 {
        let scale = options.direction_derivative_tolerance
            + options.integration.relative_tolerance * vec_norm(read_vec(&trial.state, 15 + 6 * i));
        let normalized = trial.column_projection_error[i] / scale;
        if !scale.is_finite() || scale <= 0.0 || !normalized.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        error = error.max(normalized);
    }
    Ok(error)
}

#[cfg(test)]
mod tests;
