//! Whole-element harmonic strain observation in the original MITC trial space.
//!
//! Membrane strain is constant, curvature is quadratic, and assumed transverse
//! shear is affine in triangle coordinates. Six degree-two Bernstein controls
//! therefore represent the complete strain polynomial, including the original
//! statically condensed director bubble. They are coefficients, not six sample
//! points at which a maximum is presumed to occur. Nonnegative Bernstein
//! weights sum to one everywhere in the triangle; the norm of a convex
//! combination cannot exceed the largest control norm. Strain is affine in
//! thickness, so the same argument covers every point between both layer faces.
//!
//! For each complex control tensor, the cycle maximum of the Frobenius norm is
//! the square root of the largest eigenvalue of its real/imaginary Gram matrix.
//! This retains coherent phase between components. The resulting spatial bound
//! is conservative and may exceed the actual maximum. Arithmetic is ordinary
//! deterministic floating point, not an outward-rounded interval certificate;
//! condensation conditioning and roundoff remain those of the original model.
//! The cold observation map changes summation order relative to direct point
//! evaluation; equivalence is to floating-point roundoff, not bit identity.
//! This is an incremental, linear, equivalent-single-layer shell observation,
//! not a three-dimensional interlaminar recovery or a measured damage limit.

use super::{NodalDirectors, PlateError, ShearShellSection, ShellMesh, dot, element};
use fs_math::det;

const POINTS: [(f64, f64); 6] = [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (0.5, 0.0),
    (0.5, 0.5),
    (0.0, 0.5),
];
const EDGES: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];
type Controls = [[f64; 8]; 6];

fn invalid(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}

/// Fixed Bernstein observation of the original condensed element trial space.
/// The 48-by-18 real map is prepared once from the original point evaluator.
/// Its size is independent of the number of retained modes and frequencies.
/// No geometry, duplicate kinematics or elastic factorization is retained.
#[derive(Debug, Clone)]
pub struct PreparedElementStrain {
    observation: [[[f64; 18]; 8]; 6],
}

impl PreparedElementStrain {
    /// Coarse fixed scalar-work ceiling for one preparation, including the
    /// original 20-by-20, seven-point elastic form, bubble condensation and
    /// observation of all 18 physical nodal basis vectors at six controls.
    /// Admission units are deliberately conservative, not CPU instructions.
    pub const PREPARATION_WORK: usize = 1_048_576;
    /// Coarse fixed scalar-work ceiling for the two 48-by-18 map products.
    /// No original kinematic basis is evaluated on this response path.
    pub const HARMONIC_FIELD_WORK: usize = 8_192;
    /// One real 48-by-18 product for an actual displacement snapshot.
    pub const REAL_FIELD_WORK: usize = 4_096;
    /// Whole-layer real tensor norms, with no harmonic phase maximization.
    pub const REAL_LAYER_BOUND_WORK: usize = 2_048;
    /// Coarse fixed scalar-work ceiling for one complete layer bound.
    pub const LAYER_BOUND_WORK: usize = 4_096;
    /// Work ceiling for the convenience field-plus-layer observation.
    pub const EVALUATION_WORK: usize = Self::HARMONIC_FIELD_WORK + Self::LAYER_BOUND_WORK;
    /// Exact retained inline storage; preparation has no retained heap data.
    pub const RETAINED_BYTES: usize = core::mem::size_of::<Self>();

    /// Prepare one facet from the exact mesh, directors and equilibrium section
    /// used by the structural pencil. Only this facet is inspected: constructing
    /// observers for E facets does not repeat an O(E) global mesh validation.
    /// The original [`NodalDirectors::new`] owns the whole-mesh director audit.
    ///
    /// The section must be the original elastic section that defined static
    /// condensation, including its A/B/D coupling, rather than a new modulus
    /// evaluated at the requested response frequency.
    pub fn new(
        mesh: &ShellMesh,
        element: usize,
        directors: &NodalDirectors,
        section: &ShearShellSection,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        poll(cancelled)?;
        let geometry = mesh.facet(element)?;
        if directors.vectors.len() != mesh.nodes.len() {
            return Err(PlateError::BadDirector {
                node: directors.vectors.len().min(mesh.nodes.len()),
            });
        }
        let tri = &mesh.tris[element];
        let normals = tri.map(|node| directors.vectors[node]);
        for (&node, normal) in tri.iter().zip(normals) {
            if normal.iter().any(|v| !v.is_finite())
                || (dot(normal, normal) - 1.0).abs() > 1e-12
                || dot(geometry.frame[2], normal) < directors.minimum_facet_dot
            {
                return Err(PlateError::BadDirector { node });
            }
        }
        section.validate()?;
        let field =
            element::PreparedStrainField::new(mesh, tri, geometry, normals, section, cancelled)?;
        let mut observation = [[[0.0; 18]; 8]; 6];
        for column in 0..18 {
            poll(cancelled)?;
            let mut basis = [0.0; 18];
            basis[column] = 1.0;
            let controls = original_controls(&field, &basis, cancelled)?;
            for control in 0..6 {
                for component in 0..8 {
                    observation[control][component][column] = controls[control][component];
                }
            }
        }
        poll(cancelled)?;
        Ok(Self { observation })
    }

    /// Reconstruct coherent Bernstein coefficients from physical global nodal
    /// phasors in triangle-node order, six DOFs per node. Translations are metres
    /// and rotations radians. Local support coordinates must already be expanded
    /// by their original shell owner. Both signs of a harmonic time convention
    /// give the same reported cycle envelope.
    pub fn harmonic_field(
        &self,
        real: &[f64; 18],
        imaginary: &[f64; 18],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicElementStrain, PlateError> {
        poll(cancelled)?;
        if real.iter().chain(imaginary).any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite physical harmonic element displacement"));
        }
        let mut field = HarmonicElementStrain {
            real: [[0.0; 8]; 6],
            imaginary: [[0.0; 8]; 6],
        };
        for control in 0..6 {
            poll(cancelled)?;
            for component in 0..8 {
                let row = &self.observation[control][component];
                let mut re = 0.0;
                let mut im = 0.0;
                for j in 0..18 {
                    re += row[j] * real[j];
                    im += row[j] * imaginary[j];
                }
                if !re.is_finite() || !im.is_finite() {
                    return Err(invalid(
                        "unrepresentable harmonic Bernstein strain controls",
                    ));
                }
                field.real[control][component] = re;
                field.imaginary[control][component] = im;
            }
        }
        Ok(field)
    }

    /// Recover one actual displacement snapshot using the SAME condensed map.
    /// The caller supplies all retained motion, including any retained static
    /// offset. Removed prestress, cure strain and omitted modes are not inferred.
    /// This needs one real map product and no complex temporary or cycle solve.
    pub fn real_field(
        &self,
        displacement: &[f64; 18],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RealElementStrain, PlateError> {
        poll(cancelled)?;
        if displacement.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite physical element displacement"));
        }
        let mut controls = [[0.0; 8]; 6];
        for (control, values) in controls.iter_mut().enumerate() {
            poll(cancelled)?;
            for (component, value) in values.iter_mut().enumerate() {
                for (coefficient, q) in self.observation[control][component]
                    .iter()
                    .zip(displacement)
                {
                    *value += coefficient * q;
                }
                if !value.is_finite() {
                    return Err(invalid("unrepresentable real Bernstein strain controls"));
                }
            }
        }
        Ok(RealElementStrain { controls })
    }

    /// Convenience for one film. Reuse [`Self::harmonic_field`] when multiple
    /// films share this facet and displacement, to avoid repeated observation.
    pub fn harmonic_layer_bound(
        &self,
        real: &[f64; 18],
        imaginary: &[f64; 18],
        faces_m: [f64; 2],
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TensorStrainBound, PlateError> {
        self.harmonic_field(real, imaginary, cancelled)?
            .isotropic_layer_bound(faces_m, poisson_ratio, cancelled)
    }
}

// Cold map construction samples the original condensed evaluator. Keeping this
// conversion here avoids a second set of strain, tying or bubble formulas.
fn original_controls(
    field: &element::PreparedStrainField,
    displacement: &[f64; 18],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Controls, PlateError> {
    let expanded = field.expand(displacement);
    expanded_controls(field, &expanded, cancelled)
}

// Thermal condensation adds an affine director offset to these same twenty
// coordinates. Both observers use the original point evaluator and Bernstein
// conversion, including that offset, rather than reconstructing kinematics.
pub(super) fn real_field_from_expanded(
    field: &element::PreparedStrainField,
    expanded: &[f64; 20],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealElementStrain, PlateError> {
    Ok(RealElementStrain {
        controls: expanded_controls(field, expanded, cancelled)?,
    })
}

fn expanded_controls(
    field: &element::PreparedStrainField,
    expanded: &[f64; 20],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Controls, PlateError> {
    if expanded.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "unrepresentable condensed harmonic element displacement",
        ));
    }
    let mut values = [[0.0; 8]; 6];
    for (i, &(r, s)) in POINTS.iter().enumerate() {
        poll(cancelled)?;
        let (membrane, curvature, shear) = field.eval(expanded, r, s);
        values[i][..3].copy_from_slice(&membrane);
        values[i][3..6].copy_from_slice(&curvature);
        values[i][6..].copy_from_slice(&shear);
    }
    // B_ij = 2*p(mid_ij) - (p(vertex_i) + p(vertex_j))/2.
    // Vertices remain untouched while all three edge controls are formed.
    for (edge, &(i, j)) in EDGES.iter().enumerate() {
        for component in 0..8 {
            values[3 + edge][component] = 2.0 * values[3 + edge][component]
                - 0.5 * values[i][component]
                - 0.5 * values[j][component];
        }
    }
    if values.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid(
            "unrepresentable harmonic Bernstein strain controls",
        ));
    }
    Ok(values)
}

/// One coherent complex strain polynomial in the original element trial space.
/// It can be observed on several signed physical layer intervals without
/// reconstructing the nodal field or the condensed director rotations again.
#[derive(Debug, Clone)]
pub struct HarmonicElementStrain {
    real: Controls,
    imaginary: Controls,
}

/// One actual real strain polynomial. Whole-layer bounds cover space at this
/// instant only; they do not bound unsampled motion between time endpoints.
#[derive(Debug, Clone)]
pub struct RealElementStrain {
    controls: Controls,
}

/// Engineering strains in the original facet frame. At signed thickness z,
/// in-plane total strain is `membrane + z * gradient_per_m`; transverse shear
/// is the native assumed field, constant through thickness. No thermal strain
/// is subtracted here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellStrainPoint {
    pub membrane: [f64; 3],
    pub gradient_per_m: [f64; 3],
    pub transverse_shear: [f64; 2],
}

impl ShellStrainPoint {
    fn from_coefficients(c: [f64; 8]) -> Self {
        Self {
            membrane: [c[0], c[1], c[2]],
            gradient_per_m: [c[3], c[4], c[5]],
            transverse_shear: [c[6], c[7]],
        }
    }
}

/// The SAME positive degree-five triangle rule used by native elastic and
/// thermal recovery. Barycentric points accompany area-normalized weights.
/// A consumer integrating material work should multiply by the original facet
/// area and retain its own thickness integration. Preparing this seven-entry
/// rule once avoids duplicating quadrature constants in coupled owners.
pub fn integration_rule() -> Vec<([f64; 3], f64)> {
    element::stiffness_rule()
}

impl RealElementStrain {
    pub const RETAINED_BYTES: usize = core::mem::size_of::<Self>();

    /// Evaluate the complete native strain polynomial at one point in the
    /// original triangle. This also retains an affine thermal bubble whenever
    /// this field came from `PreparedThermalShell::strain_field`.
    pub fn at_barycentric(&self, lambda: [f64; 3]) -> Result<ShellStrainPoint, PlateError> {
        if lambda.iter().any(|v| !v.is_finite() || *v < 0.0 || *v > 1.0)
            || (lambda.iter().sum::<f64>() - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(invalid("finite barycentric point inside the original facet required"));
        }
        let [a, b, c] = lambda;
        let weights = [a * a, b * b, c * c, 2.0 * a * b, 2.0 * b * c, 2.0 * c * a];
        let values = core::array::from_fn(|component| {
            (0..6).map(|i| weights[i] * self.controls[i][component]).sum()
        });
        if values.iter().any(|v: &f64| !v.is_finite()) {
            return Err(invalid("unrepresentable native point strain"));
        }
        Ok(ShellStrainPoint::from_coefficients(values))
    }

    /// Six Bernstein coefficients, not six physical samples. Their convex
    /// hull encloses the complete represented field. Constitutive consumers
    /// can use the hull and both physical layer faces to admit a rectangular
    /// total/elastic strain domain before evaluating integration points.
    pub fn bernstein_controls(&self) -> [ShellStrainPoint; 6] {
        self.controls.map(ShellStrainPoint::from_coefficients)
    }

    /// Bound the represented elastic strain over a whole facet and layer.
    /// `free_at_faces` is an affine engineering eigenstrain in this facet's
    /// axes. The norm includes xx, yy, xy, xz and yz, with tensor shear weights;
    /// no unmeasured out-of-plane normal strain or material failure criterion
    /// is inferred. Bernstein convexity covers the complete spatial polynomial,
    /// including an affine thermal bubble when supplied by the thermal owner.
    pub fn reduced_elastic_layer_bound(
        &self,
        faces_m: [f64; 2],
        free_at_faces: [[f64; 3]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, PlateError> {
        poll(cancelled)?;
        if faces_m.iter().any(|z| !z.is_finite())
            || faces_m[0] > faces_m[1]
            || free_at_faces.iter().flatten().any(|e| !e.is_finite())
        {
            return Err(invalid(
                "finite ordered faces and integrated free strains required",
            ));
        }
        let mut bound = 0.0_f64;
        for (face, &z) in faces_m.iter().enumerate() {
            for control in &self.controls {
                poll(cancelled)?;
                let values = [
                    control[0] + z * control[3] - free_at_faces[face][0],
                    control[1] + z * control[4] - free_at_faces[face][1],
                    control[2] + z * control[5] - free_at_faces[face][2],
                    control[6],
                    control[7],
                ];
                if values.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("unrepresentable thermal layer strain"));
                }
                let scale = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                let norm = if scale == 0.0 {
                    0.0
                } else {
                    scale
                        * det::sqrt(
                            values
                                .iter()
                                .enumerate()
                                .map(|(i, v)| {
                                    let x = v / scale;
                                    if i < 2 { x * x } else { 0.5 * x * x }
                                })
                                .sum(),
                        )
                };
                if !norm.is_finite() {
                    return Err(invalid("unrepresentable thermal layer strain bound"));
                }
                bound = bound.max(norm);
            }
        }
        Ok(bound)
    }

    /// Bound the full tensor of an isotropic plane-stress film using the same
    /// layer faces, tensor shear convention and Bernstein hull as harmonic
    /// observation. No harmonic phase envelope or temporal interpolation is
    /// introduced. Floating-point and equivalent-single-layer limits remain.
    pub fn isotropic_layer_bound(
        &self,
        faces_m: [f64; 2],
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TensorStrainBound, PlateError> {
        poll(cancelled)?;
        validate_layer(faces_m, poisson_ratio)?;
        let mut bound = empty_bound();
        for (face, &z) in faces_m.iter().enumerate() {
            for (control, coefficients) in self.controls.iter().enumerate() {
                poll(cancelled)?;
                let tensor = plane_stress_tensor(coefficients, z, poisson_ratio)?;
                let scale = tensor.iter().fold(0.0_f64, |m, value| m.max(value.abs()));
                let mut squared = 0.0;
                if scale > 0.0 {
                    for (i, value) in tensor.iter().enumerate() {
                        let normalized = value / scale;
                        squared += if i < 3 {
                            normalized * normalized
                        } else {
                            2.0 * normalized * normalized
                        };
                    }
                }
                let norm = scale * det::sqrt(squared);
                if !norm.is_finite() {
                    return Err(invalid("unrepresentable real tensor norm"));
                }
                if norm > bound.frobenius_upper_bound {
                    bound.frobenius_upper_bound = norm;
                    bound.ruling_control_point = control;
                    bound.ruling_z_face = face;
                }
                for (amplitude, value) in bound.component_amplitudes.iter_mut().zip(tensor) {
                    *amplitude = amplitude.max(value.abs());
                }
            }
        }
        Ok(bound)
    }
}

/// Conservative whole-triangle, whole-layer tensor envelope. Harmonic queries
/// additionally maximize over one cycle; real queries describe one instant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorStrainBound {
    /// Upper bound for max ||epsilon(t,x,z)||_F over this triangle and layer.
    /// Tensor shear entries carry weight two in the squared Frobenius norm.
    pub frobenius_upper_bound: f64,
    /// Individual component amplitude bounds in local facet coordinates,
    /// tensor-shear order [xx, yy, zz, xy, yz, zx]. These independent component
    /// maxima need not occur at the same phase or location.
    pub component_amplitudes: [f64; 6],
    /// Bernstein coefficient controlling the reported Frobenius bound:
    /// 0/1/2 are vertices; 3/4/5 are edge controls 01/12/20. An edge control is
    /// not the physical midpoint value and this is not a maximum-location claim.
    pub ruling_control_point: usize,
    /// Index 0 or 1 in the caller's ordered `faces_m` interval.
    pub ruling_z_face: usize,
}

impl HarmonicElementStrain {
    /// Exact retained inline storage; no heap data.
    pub const RETAINED_BYTES: usize = core::mem::size_of::<Self>();

    /// Bound an isotropic, real constant-Poisson-ratio plane-stress film.
    /// Reconstruct xx/yy/gamma_xy as membrane + z*curvature; use
    /// epsilon_zz = -nu/(1-nu)*(epsilon_xx+epsilon_yy). All engineering shears
    /// are divided by two before forming the three-dimensional tensor.
    /// The section shear correction belongs to energy, and never scales this
    /// kinematic strain. Transverse shear remains the original assumed MITC
    /// equivalent-single-layer field, not a recovered interlaminar profile.
    ///
    /// Faces must be finite and ordered, and -1 < nu < 0.5. Compare the result
    /// to an explicitly authored Frobenius small-strain envelope; it is not a
    /// source-derived failure threshold and excludes any mean/prestrain field.
    pub fn isotropic_layer_bound(
        &self,
        faces_m: [f64; 2],
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TensorStrainBound, PlateError> {
        poll(cancelled)?;
        validate_layer(faces_m, poisson_ratio)?;
        let mut bound = empty_bound();
        for (face, &z) in faces_m.iter().enumerate() {
            for control in 0..6 {
                poll(cancelled)?;
                let real = plane_stress_tensor(&self.real[control], z, poisson_ratio)?;
                let imaginary = plane_stress_tensor(&self.imaginary[control], z, poisson_ratio)?;
                let peak = tensor_cycle_norm(&real, &imaginary)?;
                if peak > bound.frobenius_upper_bound {
                    bound.frobenius_upper_bound = peak;
                    bound.ruling_control_point = control;
                    bound.ruling_z_face = face;
                }
                for i in 0..6 {
                    bound.component_amplitudes[i] =
                        bound.component_amplitudes[i].max(hypot(real[i], imaginary[i])?);
                }
            }
        }
        Ok(bound)
    }
}

fn validate_layer(faces_m: [f64; 2], poisson_ratio: f64) -> Result<(), PlateError> {
    if faces_m.iter().any(|z| !z.is_finite())
        || faces_m[0] > faces_m[1]
        || !poisson_ratio.is_finite()
        || poisson_ratio <= -1.0
        || poisson_ratio >= 0.5
    {
        return Err(invalid(
            "finite ordered layer faces and -1 < Poisson ratio < 0.5 required",
        ));
    }
    Ok(())
}

fn empty_bound() -> TensorStrainBound {
    TensorStrainBound {
        frobenius_upper_bound: 0.0,
        component_amplitudes: [0.0; 6],
        ruling_control_point: 0,
        ruling_z_face: 0,
    }
}

fn plane_stress_tensor(strain: &[f64; 8], z: f64, nu: f64) -> Result<[f64; 6], PlateError> {
    let xx = strain[0] + z * strain[3];
    let yy = strain[1] + z * strain[4];
    let tensor = [
        xx,
        yy,
        -nu / (1.0 - nu) * (xx + yy),
        0.5 * (strain[2] + z * strain[5]),
        0.5 * strain[7],
        0.5 * strain[6],
    ];
    if tensor.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable harmonic layer tensor"));
    }
    Ok(tensor)
}

fn hypot(a: f64, b: f64) -> Result<f64, PlateError> {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    let value = scale * det::sqrt((a / scale).powi(2) + (b / scale).powi(2));
    if !value.is_finite() {
        return Err(invalid("unrepresentable harmonic strain amplitude"));
    }
    Ok(value)
}

fn tensor_cycle_norm(real: &[f64; 6], imaginary: &[f64; 6]) -> Result<f64, PlateError> {
    // Scale before squaring: very small/large but representable strains should
    // not silently underflow to zero or overflow their Gram matrix.
    let scale = real
        .iter()
        .chain(imaginary)
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 {
        return Ok(0.0);
    }
    let mut aa = 0.0;
    let mut bb = 0.0;
    let mut ab = 0.0;
    for i in 0..6 {
        let a = real[i] / scale;
        let b = imaginary[i] / scale;
        let weight = if i < 3 { 1.0 } else { 2.0 };
        aa += weight * a * a;
        bb += weight * b * b;
        ab += weight * a * b;
    }
    let eigenvalue = 0.5 * (aa + bb + hypot(aa - bb, 2.0 * ab)?);
    let peak = scale * det::sqrt(eigenvalue);
    if !peak.is_finite() {
        return Err(invalid("unrepresentable harmonic tensor norm"));
    }
    Ok(peak)
}

#[cfg(test)]
#[path = "strain/tests.rs"]
mod tests;
