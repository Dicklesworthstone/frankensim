//! Datum-priority (3-2-1) registration against declared datum features.
//!
//! Engineering drawings define tolerances from DATUM features with a strict
//! priority order: the primary datum A (a plane) constrains orientation and
//! one translation, the secondary datum B (a line/direction) constrains the
//! remaining in-plane rotation and one translation, and the tertiary datum C
//! (a point) fixes the last translation. A global best fit can hide exactly
//! the local deviation a datum-based scheme is designed to expose, because it
//! spreads a single defective feature's error across every fiducial. This
//! module therefore reports residuals PER DATUM and publishes the
//! datum-versus-global delta as a first-class diagnostic instead of choosing
//! one answer silently.
//!
//! Constraint sequencing is structural: datum B's out-of-plane information is
//! discarded by construction (only its component in the A plane is used), and
//! datum C contributes only the final translation direction. Perturbing B
//! along the A normal or C transverse to its constraint direction provably
//! cannot move the datum pose.
//!
//! The registration direction matches the crate convention: design →
//! measured. Orientation pairing assumes the as-built pose is within 90
//! degrees of nominal for the A normal and the projected B direction (parts
//! are fixtured approximately right before fine registration); an exactly
//! perpendicular pairing refuses rather than guessing a sign.

#![allow(clippy::needless_range_loop)] // Fixed 3x3 indices expose the axis ordering.
#![allow(clippy::float_cmp)] // Exact zeros distinguish structural refusals from IEEE noise.

use crate::rigid3::{
    Fiducial3, Point3, Rigid3Error, Rigid3Registration, add3, cross3, dot3, mat3_mul,
    mat3_transpose, mat3_vec, norm3, register3, rotation_angle_from_matrix, scale3, sub3,
};

type Vec3 = [f64; 3];
type Mat3 = [[f64; 3]; 3];

const POLL_STRIDE: usize = 256;
/// Relative spectral-rank gate for datum feature fits, matching the module
/// rank-gate convention in `rigid3`.
const RANK_RELATIVE_TOLERANCE: f64 = 1e-12;

/// Datum feature label, in priority order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatumLabel {
    /// Primary plane datum.
    A,
    /// Secondary line/direction datum.
    B,
    /// Tertiary point datum.
    C,
}

impl core::fmt::Display for DatumLabel {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::A => write!(formatter, "A"),
            Self::B => write!(formatter, "B"),
            Self::C => write!(formatter, "C"),
        }
    }
}

/// Which correspondence side a datum feature fit failed on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitSide {
    /// The design (nominal) point set.
    Design,
    /// The measured (as-built) point set.
    Measured,
}

impl core::fmt::Display for FitSide {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Design => write!(formatter, "design"),
            Self::Measured => write!(formatter, "measured"),
        }
    }
}

/// Geometric diagnosis of a degenerate datum feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatumDegeneracy {
    /// The datum's points coincide.
    Coincident,
    /// The datum's points are collinear where a plane is required.
    Collinear,
}

impl core::fmt::Display for DatumDegeneracy {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coincident => write!(formatter, "points coincide"),
            Self::Collinear => write!(formatter, "points are collinear"),
        }
    }
}

/// A structured datum-registration failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatumError {
    /// The datum system shape is invalid.
    System {
        /// Stable reason.
        reason: &'static str,
    },
    /// A datum index does not address a supplied fiducial.
    IndexOutOfRange {
        /// Offending index.
        index: usize,
        /// Supplied fiducial count.
        count: usize,
    },
    /// A fiducial index appears in more than one datum role.
    DuplicateIndex {
        /// Offending index.
        index: usize,
    },
    /// A datum feature fit is geometrically degenerate.
    DegenerateDatum {
        /// Which datum.
        datum: DatumLabel,
        /// Which side.
        side: FitSide,
        /// Geometric diagnosis.
        diagnosis: DatumDegeneracy,
    },
    /// The B direction is parallel to the A normal, so it constrains no
    /// in-plane rotation.
    BParallelToA {
        /// Which side collapsed under projection.
        side: FitSide,
    },
    /// A design/measured feature pairing is exactly perpendicular, so the
    /// sign convention cannot orient it. The precondition is an as-built
    /// pose within 90 degrees of nominal.
    OrientationAmbiguous {
        /// Which datum pairing was ambiguous.
        datum: DatumLabel,
    },
    /// The embedded global fit or a shared numeric refusal.
    Registration(Rigid3Error),
    /// A finite input produced an unrepresentable aggregate.
    ArithmeticOverflow {
        /// Stable aggregate name.
        field: &'static str,
    },
    /// A bounded output allocation could not be reserved.
    AllocationFailed,
    /// Cancellation was observed at a bounded scan boundary.
    Cancelled {
        /// Stable phase name.
        phase: &'static str,
    },
}

impl core::fmt::Display for DatumError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::System { reason } => write!(formatter, "invalid datum system: {reason}"),
            Self::IndexOutOfRange { index, count } => {
                write!(
                    formatter,
                    "datum index {index} out of range for {count} fiducials"
                )
            }
            Self::DuplicateIndex { index } => {
                write!(
                    formatter,
                    "fiducial {index} appears in more than one datum role"
                )
            }
            Self::DegenerateDatum {
                datum,
                side,
                diagnosis,
            } => write!(
                formatter,
                "datum {datum} {side} feature is degenerate: {diagnosis}"
            ),
            Self::BParallelToA { side } => write!(
                formatter,
                "datum B {side} direction is parallel to the datum A normal"
            ),
            Self::OrientationAmbiguous { datum } => write!(
                formatter,
                "datum {datum} design/measured pairing is exactly perpendicular; orientation sign is ambiguous"
            ),
            Self::Registration(error) => write!(formatter, "embedded registration failed: {error}"),
            Self::ArithmeticOverflow { field } => {
                write!(formatter, "{field} overflowed the finite binary64 range")
            }
            Self::AllocationFailed => write!(formatter, "bounded output allocation failed"),
            Self::Cancelled { phase } => write!(formatter, "cancelled during {phase}"),
        }
    }
}

impl core::error::Error for DatumError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Registration(error) => Some(error),
            _ => None,
        }
    }
}

impl From<Rigid3Error> for DatumError {
    fn from(error: Rigid3Error) -> Self {
        Self::Registration(error)
    }
}

/// A declared datum system: fiducial indices for the primary plane (A, at
/// least three points), the secondary direction (B, at least two points), and
/// the tertiary point (C). The three roles must be disjoint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatumSystem {
    a: Vec<usize>,
    b: Vec<usize>,
    c: usize,
}

impl DatumSystem {
    /// Declare a datum system over fiducial indices.
    ///
    /// # Errors
    /// Too few A or B targets, or an index appearing in more than one role.
    pub fn new(a: Vec<usize>, b: Vec<usize>, c: usize) -> Result<Self, DatumError> {
        if a.len() < 3 {
            return Err(DatumError::System {
                reason: "datum A needs at least three target indices",
            });
        }
        if b.len() < 2 {
            return Err(DatumError::System {
                reason: "datum B needs at least two target indices",
            });
        }
        let mut all = Vec::new();
        all.try_reserve_exact(a.len() + b.len() + 1)
            .map_err(|_| DatumError::AllocationFailed)?;
        all.extend_from_slice(&a);
        all.extend_from_slice(&b);
        all.push(c);
        let mut sorted = all.clone();
        sorted.sort_unstable();
        for pair in sorted.windows(2) {
            if pair[0] == pair[1] {
                return Err(DatumError::DuplicateIndex { index: pair[0] });
            }
        }
        Ok(Self { a, b, c })
    }

    /// Primary plane target indices.
    #[must_use]
    pub fn a(&self) -> &[usize] {
        &self.a
    }

    /// Secondary direction target indices.
    #[must_use]
    pub fn b(&self) -> &[usize] {
        &self.b
    }

    /// Tertiary point target index.
    #[must_use]
    pub const fn c(&self) -> usize {
        self.c
    }
}

/// The orthonormal measured-side constraint frame the datum translation was
/// solved in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DatumFrame {
    plane_normal: Vec3,
    in_plane_transverse: Vec3,
    line_direction: Vec3,
}

impl DatumFrame {
    /// Measured datum A plane unit normal (the A constraint direction).
    #[must_use]
    pub const fn plane_normal(&self) -> Vec3 {
        self.plane_normal
    }

    /// In-plane unit direction perpendicular to the B line (the B constraint
    /// direction).
    #[must_use]
    pub const fn in_plane_transverse(&self) -> Vec3 {
        self.in_plane_transverse
    }

    /// In-plane unit direction along the B line (the C constraint
    /// direction).
    #[must_use]
    pub const fn line_direction(&self) -> Vec3 {
        self.line_direction
    }
}

/// The datum-versus-global diagnostic delta.
#[derive(Debug, Clone, PartialEq)]
pub struct DatumGlobalComparison {
    rotation_delta_rad: f64,
    translation_delta: Vec3,
    residual_norm_deltas: Vec<f64>,
    global: Rigid3Registration,
}

impl DatumGlobalComparison {
    /// Rotation angle of the relative transform between the datum pose and
    /// the global best-fit pose.
    #[must_use]
    pub const fn rotation_delta_rad(&self) -> f64 {
        self.rotation_delta_rad
    }

    /// Translation of the relative transform between the datum pose and the
    /// global best-fit pose.
    #[must_use]
    pub const fn translation_delta(&self) -> Vec3 {
        self.translation_delta
    }

    /// Per-fiducial `||r_datum|| - ||r_global||`. A large-magnitude entry
    /// marks a fiducial the two frames disagree about; its sign is
    /// geometry-dependent (each fit absorbs deviation in proportion to its
    /// own leverage structure), so investigate large entries in the datum
    /// frame the drawing tolerances rather than reading the sign as a
    /// verdict.
    #[must_use]
    pub fn residual_norm_deltas(&self) -> &[f64] {
        &self.residual_norm_deltas
    }

    /// The embedded global Kabsch registration.
    #[must_use]
    pub const fn global(&self) -> &Rigid3Registration {
        &self.global
    }
}

/// A datum-priority registration (design → measured) with per-datum
/// residuals and the datum-versus-global diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub struct DatumRegistration {
    rotation: Mat3,
    translation: Vec3,
    frame: DatumFrame,
    a_out_of_plane: Vec<f64>,
    b_off_line: Vec<f64>,
    c_along_line: f64,
    residual_norms: Vec<f64>,
    comparison: DatumGlobalComparison,
}

impl DatumRegistration {
    /// Row-major rotation matrix (design → measured).
    #[must_use]
    pub const fn rotation(&self) -> &Mat3 {
        &self.rotation
    }

    /// Translation vector (design → measured).
    #[must_use]
    pub const fn translation(&self) -> Vec3 {
        self.translation
    }

    /// The measured-side constraint frame used for the translation solve and
    /// the per-datum residual decomposition.
    #[must_use]
    pub const fn frame(&self) -> &DatumFrame {
        &self.frame
    }

    /// Signed out-of-plane residual of each datum A target, in A-index
    /// order, measured along the plane normal.
    #[must_use]
    pub fn a_out_of_plane(&self) -> &[f64] {
        &self.a_out_of_plane
    }

    /// Signed off-line residual of each datum B target, in B-index order,
    /// measured along the in-plane transverse direction.
    #[must_use]
    pub fn b_off_line(&self) -> &[f64] {
        &self.b_off_line
    }

    /// Signed residual of the datum C target along the line direction.
    #[must_use]
    pub const fn c_along_line(&self) -> f64 {
        self.c_along_line
    }

    /// Residual norm of every fiducial under the datum pose, in input order.
    #[must_use]
    pub fn residual_norms(&self) -> &[f64] {
        &self.residual_norms
    }

    /// The datum-versus-global diagnostic delta.
    #[must_use]
    pub const fn comparison(&self) -> &DatumGlobalComparison {
        &self.comparison
    }

    /// Rotation angle in radians, recovered from the matrix trace.
    #[must_use]
    pub fn rotation_angle_rad(&self) -> f64 {
        rotation_angle_from_matrix(&self.rotation)
    }

    /// Apply the datum registration to a design point.
    ///
    /// # Errors
    /// Refuses non-finite arithmetic overflow.
    pub fn apply(&self, point: Point3) -> Result<Point3, DatumError> {
        let mapped = add3(mat3_vec(&self.rotation, point.coords()), self.translation);
        Point3::new(mapped[0], mapped[1], mapped[2]).map_err(|_| DatumError::ArithmeticOverflow {
            field: "applied point",
        })
    }
}

fn checkpoint(cx: &fs_exec::Cx<'_>, ordinal: usize, phase: &'static str) -> Result<(), DatumError> {
    if ordinal.is_multiple_of(POLL_STRIDE) {
        cx.checkpoint()
            .map_err(|_| DatumError::Cancelled { phase })?;
    }
    Ok(())
}

struct FeatureFit {
    centroid: Vec3,
    /// Eigenvalues of the deviation scatter, ascending.
    eigenvalues: Vec3,
    /// Matching eigenvectors as columns, ascending eigenvalue order.
    eigenvectors: Mat3,
}

fn fit_feature(
    fiducials: &[Fiducial3],
    indices: &[usize],
    select_measured: bool,
    cx: &fs_exec::Cx<'_>,
    phase: &'static str,
) -> Result<FeatureFit, DatumError> {
    let mut centroid = [0.0f64; 3];
    for (ordinal, index) in indices.iter().enumerate() {
        checkpoint(cx, ordinal, phase)?;
        let point = if select_measured {
            fiducials[*index].measured().coords()
        } else {
            fiducials[*index].design().coords()
        };
        let gain = 1.0 / (ordinal as f64 + 1.0);
        for axis in 0..3 {
            centroid[axis] += gain * (point[axis] - centroid[axis]);
        }
    }
    let mut scatter = [[0.0f64; 3]; 3];
    for (ordinal, index) in indices.iter().enumerate() {
        checkpoint(cx, ordinal, phase)?;
        let point = if select_measured {
            fiducials[*index].measured().coords()
        } else {
            fiducials[*index].design().coords()
        };
        let deviation = sub3(point, centroid);
        for row in 0..3 {
            for column in 0..3 {
                scatter[row][column] =
                    deviation[row].mul_add(deviation[column], scatter[row][column]);
            }
        }
    }
    if !scatter.iter().flatten().all(|value| value.is_finite()) {
        return Err(DatumError::ArithmeticOverflow {
            field: "datum scatter",
        });
    }
    let flat = [
        scatter[0][0],
        scatter[0][1],
        scatter[0][2],
        scatter[1][0],
        scatter[1][1],
        scatter[1][2],
        scatter[2][0],
        scatter[2][1],
        scatter[2][2],
    ];
    let (eigenvalues, eigenvectors) = fs_la::eigen::jacobi_eigh(&flat, 3);
    let mut vectors = [[0.0f64; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            vectors[row][column] = eigenvectors[row * 3 + column];
        }
    }
    Ok(FeatureFit {
        centroid,
        eigenvalues: [eigenvalues[0], eigenvalues[1], eigenvalues[2]],
        eigenvectors: vectors,
    })
}

fn eigencolumn(fit: &FeatureFit, column: usize) -> Vec3 {
    [
        fit.eigenvectors[0][column],
        fit.eigenvectors[1][column],
        fit.eigenvectors[2][column],
    ]
}

/// Deterministic sign canonicalization: the first component whose magnitude
/// is nonzero is made positive.
fn canonical_sign(vector: Vec3) -> Vec3 {
    for component in vector {
        if component != 0.0 {
            return if component < 0.0 {
                scale3(vector, -1.0)
            } else {
                vector
            };
        }
    }
    vector
}

/// Plane fit: unit normal (smallest-eigenvalue direction) plus degeneracy
/// classification. A plane needs rank two in the deviations.
fn plane_of(fit: &FeatureFit, datum: DatumLabel, side: FitSide) -> Result<Vec3, DatumError> {
    let largest = fit.eigenvalues[2];
    if largest <= 0.0 {
        return Err(DatumError::DegenerateDatum {
            datum,
            side,
            diagnosis: DatumDegeneracy::Coincident,
        });
    }
    if fit.eigenvalues[1] <= RANK_RELATIVE_TOLERANCE * largest {
        return Err(DatumError::DegenerateDatum {
            datum,
            side,
            diagnosis: DatumDegeneracy::Collinear,
        });
    }
    Ok(canonical_sign(eigencolumn(fit, 0)))
}

/// Direction fit: unit direction (largest-eigenvalue direction) plus
/// degeneracy classification. A direction needs rank one in the deviations.
fn direction_of(fit: &FeatureFit, datum: DatumLabel, side: FitSide) -> Result<Vec3, DatumError> {
    if fit.eigenvalues[2] <= 0.0 {
        return Err(DatumError::DegenerateDatum {
            datum,
            side,
            diagnosis: DatumDegeneracy::Coincident,
        });
    }
    Ok(canonical_sign(eigencolumn(fit, 2)))
}

/// Orient `vector` so its inner product with `reference` is positive;
/// refuses an exactly perpendicular pairing.
fn orient_toward(vector: Vec3, reference: Vec3, datum: DatumLabel) -> Result<Vec3, DatumError> {
    let alignment = dot3(vector, reference);
    if alignment == 0.0 {
        return Err(DatumError::OrientationAmbiguous { datum });
    }
    Ok(if alignment < 0.0 {
        scale3(vector, -1.0)
    } else {
        vector
    })
}

/// Rodrigues rotation about a unit axis with the given cosine/sine.
fn rotation_about(axis: Vec3, cosine: f64, sine: f64) -> Mat3 {
    let mut rotation = [[0.0f64; 3]; 3];
    let versine = 1.0 - cosine;
    let cross = [
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ];
    for row in 0..3 {
        for column in 0..3 {
            let identity = if row == column { cosine } else { 0.0 };
            rotation[row][column] = (versine * axis[row])
                .mul_add(axis[column], sine.mul_add(cross[row][column], identity));
        }
    }
    rotation
}

/// Minimal rotation taking unit vector `from` onto unit vector `to`. The
/// caller guarantees `dot(from, to) > 0`.
fn minimal_rotation(from: Vec3, to: Vec3) -> Mat3 {
    let axis = cross3(from, to);
    let sine = norm3(axis);
    let cosine = dot3(from, to);
    if sine == 0.0 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    rotation_about(scale3(axis, 1.0 / sine), cosine, sine)
}

fn project_onto_plane(vector: Vec3, unit_normal: Vec3, side: FitSide) -> Result<Vec3, DatumError> {
    let out_of_plane = dot3(vector, unit_normal);
    let projected = sub3(vector, scale3(unit_normal, out_of_plane));
    let magnitude = norm3(projected);
    if magnitude <= RANK_RELATIVE_TOLERANCE * norm3(vector) || magnitude == 0.0 {
        return Err(DatumError::BParallelToA { side });
    }
    Ok(scale3(projected, 1.0 / magnitude))
}

/// Register the design onto the measured data through the declared datum
/// hierarchy: A (plane) then B (direction) then C (point), each constraint
/// consuming only the degrees of freedom its priority allows. Residuals are
/// reported per datum, and the delta against the global Kabsch best fit is a
/// first-class output.
///
/// # Errors
/// Invalid indices, degenerate datum features with a geometric diagnosis,
/// a B direction parallel to the A normal, exactly perpendicular orientation
/// pairings, embedded-registration refusals, non-finite aggregates, or a
/// structured [`DatumError::Cancelled`].
#[allow(clippy::too_many_lines)]
pub fn register3_datum(
    fiducials: &[Fiducial3],
    system: &DatumSystem,
    cx: &fs_exec::Cx<'_>,
) -> Result<DatumRegistration, DatumError> {
    for index in system
        .a
        .iter()
        .chain(system.b.iter())
        .chain(core::iter::once(&system.c))
    {
        if *index >= fiducials.len() {
            return Err(DatumError::IndexOutOfRange {
                index: *index,
                count: fiducials.len(),
            });
        }
    }

    // Datum A: plane orientation.
    let design_a = fit_feature(fiducials, &system.a, false, cx, "datum.a-design")?;
    let measured_a = fit_feature(fiducials, &system.a, true, cx, "datum.a-measured")?;
    let design_normal = plane_of(&design_a, DatumLabel::A, FitSide::Design)?;
    let measured_normal_raw = plane_of(&measured_a, DatumLabel::A, FitSide::Measured)?;
    let measured_normal = orient_toward(measured_normal_raw, design_normal, DatumLabel::A)?;
    let primary_rotation = minimal_rotation(design_normal, measured_normal);

    // Datum B: in-plane direction about the measured normal.
    let design_b = fit_feature(fiducials, &system.b, false, cx, "datum.b-design")?;
    let measured_b = fit_feature(fiducials, &system.b, true, cx, "datum.b-measured")?;
    let design_direction = direction_of(&design_b, DatumLabel::B, FitSide::Design)?;
    let measured_direction_raw = direction_of(&measured_b, DatumLabel::B, FitSide::Measured)?;
    let rotated_design_direction = mat3_vec(&primary_rotation, design_direction);
    let design_in_plane =
        project_onto_plane(rotated_design_direction, measured_normal, FitSide::Design)?;
    let measured_direction = orient_toward(measured_direction_raw, design_in_plane, DatumLabel::B)?;
    let measured_in_plane =
        project_onto_plane(measured_direction, measured_normal, FitSide::Measured)?;
    let in_plane_sine = dot3(measured_normal, cross3(design_in_plane, measured_in_plane));
    let in_plane_cosine = dot3(design_in_plane, measured_in_plane);
    let angle = in_plane_sine.atan2(in_plane_cosine);
    let secondary_rotation = rotation_about(measured_normal, angle.cos(), angle.sin());
    let rotation = mat3_mul(&secondary_rotation, &primary_rotation);

    // Constraint frame in the measured space: f1 out of plane (A), f3 along
    // the B line, f2 = f3 x f1 in-plane transverse (right-handed).
    let f1 = measured_normal;
    let f3 = measured_in_plane;
    let f2 = cross3(f3, f1);

    // Hierarchical translation: A fixes the out-of-plane component at the
    // fitted plane, B fixes the in-plane transverse component at the fitted
    // line, C fixes the along-line component at the tertiary target.
    let c_design = fiducials[system.c].design().coords();
    let c_measured = fiducials[system.c].measured().coords();
    let alpha = dot3(
        f1,
        sub3(measured_a.centroid, mat3_vec(&rotation, design_a.centroid)),
    );
    let beta = dot3(
        f2,
        sub3(measured_b.centroid, mat3_vec(&rotation, design_b.centroid)),
    );
    let gamma = dot3(f3, sub3(c_measured, mat3_vec(&rotation, c_design)));
    let translation = add3(add3(scale3(f1, alpha), scale3(f2, beta)), scale3(f3, gamma));
    if !translation.iter().all(|value| value.is_finite()) {
        return Err(DatumError::ArithmeticOverflow {
            field: "datum translation",
        });
    }

    // Residuals of every fiducial under the datum pose.
    let mut residual_vectors = Vec::new();
    residual_vectors
        .try_reserve_exact(fiducials.len())
        .map_err(|_| DatumError::AllocationFailed)?;
    let mut residual_norms = Vec::new();
    residual_norms
        .try_reserve_exact(fiducials.len())
        .map_err(|_| DatumError::AllocationFailed)?;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "datum.residuals")?;
        let mapped = add3(mat3_vec(&rotation, fiducial.design().coords()), translation);
        let residual = sub3(fiducial.measured().coords(), mapped);
        let norm = norm3(residual);
        if !norm.is_finite() {
            return Err(DatumError::ArithmeticOverflow {
                field: "datum residual",
            });
        }
        residual_vectors.push(residual);
        residual_norms.push(norm);
    }

    let mut a_out_of_plane = Vec::new();
    a_out_of_plane
        .try_reserve_exact(system.a.len())
        .map_err(|_| DatumError::AllocationFailed)?;
    for index in &system.a {
        a_out_of_plane.push(dot3(f1, residual_vectors[*index]));
    }
    let mut b_off_line = Vec::new();
    b_off_line
        .try_reserve_exact(system.b.len())
        .map_err(|_| DatumError::AllocationFailed)?;
    for index in &system.b {
        b_off_line.push(dot3(f2, residual_vectors[*index]));
    }
    let c_along_line = dot3(f3, residual_vectors[system.c]);

    // Global comparison: the datum-versus-global delta is the diagnostic.
    let global = register3(fiducials, cx)?;
    let relative_rotation = mat3_mul(&rotation, &mat3_transpose(global.rotation()));
    let rotation_delta_rad = rotation_angle_from_matrix(&relative_rotation);
    let translation_delta = sub3(
        translation,
        mat3_vec(&relative_rotation, global.translation()),
    );
    let mut residual_norm_deltas = Vec::new();
    residual_norm_deltas
        .try_reserve_exact(fiducials.len())
        .map_err(|_| DatumError::AllocationFailed)?;
    for (index, fiducial) in fiducials.iter().enumerate() {
        checkpoint(cx, index, "datum.comparison")?;
        let mapped = add3(
            mat3_vec(global.rotation(), fiducial.design().coords()),
            global.translation(),
        );
        let global_norm = norm3(sub3(fiducial.measured().coords(), mapped));
        residual_norm_deltas.push(residual_norms[index] - global_norm);
    }

    cx.checkpoint().map_err(|_| DatumError::Cancelled {
        phase: "datum.publish",
    })?;
    Ok(DatumRegistration {
        rotation,
        translation,
        frame: DatumFrame {
            plane_normal: f1,
            in_plane_transverse: f2,
            line_direction: f3,
        },
        a_out_of_plane,
        b_off_line,
        c_along_line,
        residual_norms,
        comparison: DatumGlobalComparison {
            rotation_delta_rad,
            translation_delta,
            residual_norm_deltas,
            global,
        },
    })
}
