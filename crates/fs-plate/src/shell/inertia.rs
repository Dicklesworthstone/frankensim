//! Kinetic interpolation for the existing flat-facet CST/DKT shell owner.
//!
//! P1 midsurface velocity and P1 physical director rotation are integrated
//! through the section thickness. This is not a new DKT displacement field.

use super::FacetGeometry;
use crate::{PlateError, PlateSection};

mod retained;
pub use retained::{
    ReferenceKineticIncrement, ShellReferenceMass, prepare_reference_mass,
    reference_kinetic_bilinear, reference_kinetic_increment,
};

/// Mass discretization of the flat-facet CST/DKT shell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CstDktMass {
    /// Original diagonal translation/rotation mass, retained for frozen pencils.
    Lumped,
    /// Exact P1 midsurface/director kinetic integral and a relative drilling gauge.
    /// The gauge retains the original `rho*h^3/12` scale but vanishes for rigid
    /// motion; it is numerical stabilization, not material normal-spin inertia.
    Consistent,
}

/// Reference density moments in the normalized thickness coordinate `s = z/h`.
/// `q0 = integral rho ds`, `q1 = integral rho*s ds`, and
/// `q2 = integral rho*s² ds` all have units kg/m³. The physical section moments
/// are `I0 = h*q0`, `I1 = h²*q1`, and `I2 = h³*q2` about the mesh surface.
/// Interpolating these values preserves layer fractions as measured thickness
/// varies; it does not average physical areal moments before applying thickness.
///
/// The mesh normal sets positive `z`. Reversing that normal for the same
/// specimen requires reversing `q1`. A different reference surface likewise
/// changes these moments and the translational coordinates together.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceDensityMoments {
    zeroth_kg_m3: f64,
    first_kg_m3: f64,
    second_kg_m3: f64,
}

impl ReferenceDensityMoments {
    /// Admit finite moments of an actual positive, nonzero-thickness inventory.
    /// `q0*q2 > q1²` requires strictly positive mass variance; `q1` is signed.
    /// # Errors
    /// Nonfinite moments, nonpositive mass/second moment, or degenerate variance.
    pub fn new(zeroth_kg_m3: f64, first_kg_m3: f64, second_kg_m3: f64) -> Result<Self, PlateError> {
        if !zeroth_kg_m3.is_finite()
            || zeroth_kg_m3 <= 0.0
            || !first_kg_m3.is_finite()
            || !second_kg_m3.is_finite()
            || second_kg_m3 <= 0.0
        {
            return Err(PlateError::BadSection {
                what: "shell density moments require finite positive mass and strictly positive variance",
            });
        }
        // Eliminate the larger diagonal: a valid cross moment cannot overflow
        // this division/product. A sqrt-normalized correlation can round an
        // exactly rank-one block such as [2,1;1,0.5] below one and admit it.
        let variance = (-first_kg_m3).mul_add(
            first_kg_m3 / zeroth_kg_m3.max(second_kg_m3),
            zeroth_kg_m3.min(second_kg_m3),
        );
        if !variance.is_finite() || variance <= 0.0 {
            return Err(PlateError::BadSection {
                what: "shell density moments require strictly positive representable variance",
            });
        }
        Ok(Self {
            zeroth_kg_m3,
            first_kg_m3,
            second_kg_m3,
        })
    }

    /// Uniform density symmetric about the geometric midsurface, `-1/2 <= s <= 1/2`.
    /// This limit retains the original scalar/P1-density mass kernel exactly.
    /// # Errors
    /// Nonpositive, nonfinite or unrepresentable density moments.
    pub fn from_uniform(density_kg_m3: f64) -> Result<Self, PlateError> {
        Self::new(density_kg_m3, 0.0, density_kg_m3 / 12.0)
    }

    /// Zeroth normalized moment [kg/m³]; physical areal mass is `h*q0`.
    #[must_use]
    pub const fn zeroth_kg_m3(self) -> f64 {
        self.zeroth_kg_m3
    }
    /// Signed first normalized moment [kg/m³]; physical first moment is `h²*q1`.
    #[must_use]
    pub const fn first_kg_m3(self) -> f64 {
        self.first_kg_m3
    }
    /// Second normalized moment [kg/m³]; physical rotary inertia is `h³*q2`.
    #[must_use]
    pub const fn second_kg_m3(self) -> f64 {
        self.second_kg_m3
    }
}

/// The existing assembly has one optional authored reference inventory.
#[derive(Clone, Copy)]
pub(super) enum ReferenceInertia<'a> {
    Density(&'a [[f64; 3]]),
    Moments(&'a [[ReferenceDensityMoments; 3]]),
}

impl ReferenceInertia<'_> {
    pub(super) fn validate(self, facets: usize) -> Result<(), PlateError> {
        match self {
            Self::Density(density) => validate_reference_density(facets, density),
            Self::Moments(moments) if moments.len() != facets => Err(PlateError::SectionCount {
                expected: facets,
                actual: moments.len(),
            }),
            Self::Moments(_) => Ok(()), // Each immutable value was admitted by its constructor.
        }
    }

    pub(super) fn zeroth_at(self, facet: usize) -> [f64; 3] {
        match self {
            Self::Density(density) => density[facet],
            Self::Moments(moments) => moments[facet].map(ReferenceDensityMoments::zeroth_kg_m3),
        }
    }

    pub(super) fn consistent_mass_at(
        self,
        facet: usize,
        g: &FacetGeometry,
        section: &PlateSection,
        thickness: Option<[f64; 3]>,
    ) -> Result<[f64; 324], PlateError> {
        match self {
            Self::Density(density) => {
                consistent_mass_with_reference_density(g, section, thickness, density[facet])
            }
            Self::Moments(moments) => consistent_mass_with_reference_density_moments(
                g,
                section,
                thickness,
                moments[facet],
            ),
        }
    }
}

/// Local `[u,v,w,theta_x,theta_y,theta_z]` coordinates at each vertex.
pub(super) fn consistent_mass(
    g: &FacetGeometry,
    s: &PlateSection,
) -> Result<[f64; 324], PlateError> {
    let total = s.density * s.thickness * g.area_m2;
    let rotary = total * s.thickness * s.thickness / 12.0;
    if !total.is_finite() || total <= 0.0 || !rotary.is_finite() || rotary <= 0.0 {
        return Err(PlateError::BadSection {
            what: "consistent shell inertia is unrepresentable",
        });
    }
    let mut m = [0.0; 324];
    // Integral Ni*Nj / area = (1 + delta_ij)/12. Rotation about the
    // facet normal does not move the through-thickness material director.
    for i in 0..3 {
        for j in 0..3 {
            let gram = if i == j { 1.0 / 6.0 } else { 1.0 / 12.0 };
            for c in 0..5 {
                m[18 * (6 * i + c) + 6 * j + c] = gram * if c < 3 { total } else { rotary };
            }
        }
    }
    // The original drilling stiffness penalizes theta_z - (v_x-u_y)/2.
    // Stabilize its RATE, so an actual rigid normal spin acquires no invented
    // section inertia. Positive degree-two edge quadrature is exact for P1.
    let slots = [0, 1, 5, 6, 7, 11, 12, 13, 17];
    for h in [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] {
        let spin: [f64; 9] = std::array::from_fn(|i| match i % 3 {
            0 => 0.5 * g.gradient[i / 3][1],
            1 => -0.5 * g.gradient[i / 3][0],
            _ => h[i / 3],
        });
        for i in 0..9 {
            for j in i..9 {
                let value = rotary / 3.0 * spin[i] * spin[j];
                m[18 * slots[i] + slots[j]] += value;
                if i != j {
                    m[18 * slots[j] + slots[i]] += value;
                }
            }
        }
    }
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "consistent shell kinetic integral overflow",
        });
    }
    Ok(m)
}

/// The same kinetic field with thickness interpolated from the measured vertices.
/// The section supplies density; its facet-average thickness still owns stiffness.
/// Positive degree-five triangle quadrature integrates h*Ni*Nj and h^3*Ni*Nj
/// exactly (up to roundoff), including the relative drilling rate. The material
/// occupies -h(x,y)/2 <= z <= h(x,y)/2 about this flat reference facet.
pub(super) fn consistent_mass_with_thickness(
    g: &FacetGeometry,
    s: &PlateSection,
    thickness: [f64; 3],
) -> Result<[f64; 324], PlateError> {
    if thickness[0] == thickness[1] && thickness[1] == thickness[2] {
        // Preserve the original section-average arithmetic, including rounding,
        // for every homogeneous facet and its existing modal pencils.
        return consistent_mass(g, s);
    }
    // Symmetric seven-point degree-five rule, weights normalized to unit area.
    // The two orbits have repeated coordinate (6 +/- sqrt(15))/21 and
    // weight (155 +/- sqrt(15))/1200. All points and weights are positive.
    const A: f64 = 0.4701420641051151;
    const B: f64 = 0.05971587178976982;
    const C: f64 = 0.10128650732345634;
    const D: f64 = 0.7974269853530873;
    const W: f64 = 0.1323941527885062;
    const V: f64 = 0.12593918054482714;
    let rule = [
        ([1.0 / 3.0; 3], 0.225),
        ([B, A, A], W),
        ([A, B, A], W),
        ([A, A, B], W),
        ([D, C, C], V),
        ([C, D, C], V),
        ([C, C, D], V),
    ];
    let mut m = [0.0; 324];
    for (n, weight) in rule {
        let h = n[0] * thickness[0] + n[1] * thickness[1] + n[2] * thickness[2];
        let translation = s.density * h * g.area_m2 * weight;
        let rotary = translation * h * h / 12.0;
        if !translation.is_finite() || translation <= 0.0 || !rotary.is_finite() || rotary <= 0.0 {
            return Err(PlateError::BadSection {
                what: "variable-thickness shell inertia is unrepresentable",
            });
        }
        accumulate_point(&mut m, g, n, translation, rotary);
    }
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "variable-thickness shell kinetic integral overflow",
        });
    }
    Ok(m)
}

/// Admission is shared by numerical assembly and measured-shell mass reporting.
pub(super) fn validate_reference_density(
    facets: usize,
    density: &[[f64; 3]],
) -> Result<(), PlateError> {
    if density.len() != facets {
        return Err(PlateError::SectionCount {
            expected: facets,
            actual: density.len(),
        });
    }
    if density
        .iter()
        .flatten()
        .any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return Err(PlateError::BadSection {
            what: "shell reference density needs three positive finite corner values per facet",
        });
    }
    Ok(())
}

/// Exact original-volume mass for two P1 fields on a reference triangle.
/// The cross terms retain correlation between measured thickness and density.
pub(super) fn reference_mass_kg(
    area_m2: f64,
    thickness: [f64; 3],
    density: [f64; 3],
) -> Result<f64, PlateError> {
    let mut mass = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            mass +=
                area_m2 * density[i] * thickness[j] * if i == j { 1.0 / 6.0 } else { 1.0 / 12.0 };
        }
    }
    if !mass.is_finite() || mass <= 0.0 {
        return Err(PlateError::BadSection {
            what: "graded shell reference mass is unrepresentable",
        });
    }
    Ok(mass)
}

/// P1 density and optional P1 thickness on the original P1/director kinetic field.
/// The density is mass per original volume, constant through each local thickness.
/// Keeping all sixteen positive weights makes every local kinetic contribution
/// nonnegative, including the original relative drilling-rate gauge.
pub(super) fn consistent_mass_with_reference_density(
    g: &FacetGeometry,
    section: &PlateSection,
    thickness: Option<[f64; 3]>,
    density: [f64; 3],
) -> Result<[f64; 324], PlateError> {
    if density[0] == density[1] && density[1] == density[2] {
        // Retain the existing arithmetic, even for a measured thickness field.
        let mut uniform = *section;
        uniform.density = density[0];
        return match thickness {
            Some(h) => consistent_mass_with_thickness(g, &uniform, h),
            None => consistent_mass(g, &uniform),
        };
    }
    let thickness = thickness.unwrap_or([section.thickness; 3]);
    let mut m = [0.0; 324];
    for_each_reference_point(|n, weight| {
        let h = n[0] * thickness[0] + n[1] * thickness[1] + n[2] * thickness[2];
        let rho = n[0] * density[0] + n[1] * density[1] + n[2] * density[2];
        let translation = rho * h * g.area_m2 * weight;
        let rotary = translation * h * h / 12.0;
        if !translation.is_finite() || translation <= 0.0 || !rotary.is_finite() || rotary <= 0.0 {
            return Err(PlateError::BadSection {
                what: "graded shell inertia is unrepresentable",
            });
        }
        accumulate_point(&mut m, g, n, translation, rotary);
        Ok(())
    })?;
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "graded shell kinetic integral overflow",
        });
    }
    Ok(m)
}

/// Same physical P1 velocity/director basis with full through-thickness moments.
/// `u_dot + z*(theta_dot cross normal)` gives `+I1*u_dot*theta_y_dot`
/// and `-I1*v_dot*theta_x_dot`; these are signed mass cross blocks, not springs.
pub(super) fn consistent_mass_with_reference_density_moments(
    g: &FacetGeometry,
    section: &PlateSection,
    thickness: Option<[f64; 3]>,
    moments: [ReferenceDensityMoments; 3],
) -> Result<[f64; 324], PlateError> {
    if moments
        .iter()
        .all(|q| q.first_kg_m3 == 0.0 && q.second_kg_m3 == q.zeroth_kg_m3 / 12.0)
    {
        return consistent_mass_with_reference_density(
            g,
            section,
            thickness,
            moments.map(ReferenceDensityMoments::zeroth_kg_m3),
        );
    }
    let thickness = thickness.unwrap_or([section.thickness; 3]);
    let mut m = [0.0; 324];
    for_each_reference_point(|n, weight| {
        let (translation, first, rotary) =
            reference_point_moments(g, thickness, moments, n, weight)?;
        accumulate_point(&mut m, g, n, translation, rotary);
        accumulate_first_moment(&mut m, n, first);
        Ok(())
    })?;
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "layered shell kinetic integral overflow",
        });
    }
    Ok(m)
}

// Same admitted normalized density moments and physical thickness at a native
// reference quadrature point. Keep the original multiplication order.
fn reference_point_moments(
    g: &FacetGeometry,
    thickness: [f64; 3],
    moments: [ReferenceDensityMoments; 3],
    n: [f64; 3],
    weight: f64,
) -> Result<(f64, f64, f64), PlateError> {
    let h = n[0] * thickness[0] + n[1] * thickness[1] + n[2] * thickness[2];
    let at = |value: fn(ReferenceDensityMoments) -> f64| {
        n[0] * value(moments[0]) + n[1] * value(moments[1]) + n[2] * value(moments[2])
    };
    let translation = at(ReferenceDensityMoments::zeroth_kg_m3) * h * g.area_m2 * weight;
    let first = at(ReferenceDensityMoments::first_kg_m3) * h * h * g.area_m2 * weight;
    let rotary = at(ReferenceDensityMoments::second_kg_m3) * h * h * h * g.area_m2 * weight;
    if !translation.is_finite()
        || translation <= 0.0
        || !first.is_finite()
        || !rotary.is_finite()
        || rotary <= 0.0
    {
        return Err(PlateError::BadSection {
            what: "layered shell reference inertia is unrepresentable",
        });
    }
    Ok((translation, first, rotary))
}

// Shared positive degree-six rule for both reference-density representations.
// A positive P1 interpolation of admitted moment matrices remains positive.
fn for_each_reference_point(
    mut point: impl FnMut([f64; 3], f64) -> Result<(), PlateError>,
) -> Result<(), PlateError> {
    // Four-point Gauss-Legendre on [0,1], followed by the Duffy map
    // (N0,N1,N2)=((1-u)(1-v),u,(1-u)v), Jacobian 2*A*(1-u).
    // q2*h^3*Ni*Nj has degree six in barycentric coordinates. After the
    // map its degree is at most seven in u and six in v, so 4x4 is exact.
    const POINTS: [f64; 4] = [
        0.06943184420297371,
        0.33000947820757187,
        0.6699905217924281,
        0.9305681557970262,
    ];
    const WEIGHTS: [f64; 4] = [
        0.17392742256872693,
        0.32607257743127305,
        0.32607257743127305,
        0.17392742256872693,
    ];
    for (i, u) in POINTS.into_iter().enumerate() {
        for (j, v) in POINTS.into_iter().enumerate() {
            let n = [(1.0 - u) * (1.0 - v), u, (1.0 - u) * v];
            let weight = 2.0 * (1.0 - u) * WEIGHTS[i] * WEIGHTS[j];
            point(n, weight)?;
        }
    }
    Ok(())
}

// The same physical velocity/director and drilling-rate basis for both measured
// thickness and graded reference density; no independent inertia model is added.
fn accumulate_point(
    m: &mut [f64; 324],
    g: &FacetGeometry,
    n: [f64; 3],
    translation: f64,
    rotary: f64,
) {
    accumulate_physical_point(m, n, translation, rotary);
    accumulate_drilling_point(m, g, n, rotary);
}

// Physical P1 midsurface velocity and director rotation. Normal drilling
// rotation does not move the through-thickness material director.
fn accumulate_physical_point(m: &mut [f64; 324], n: [f64; 3], translation: f64, rotary: f64) {
    for i in 0..3 {
        for j in i..3 {
            for c in 0..5 {
                let value = n[i] * n[j] * if c < 3 { translation } else { rotary };
                m[18 * (6 * i + c) + 6 * j + c] += value;
                if i != j {
                    m[18 * (6 * j + c) + 6 * i + c] += value;
                }
            }
        }
    }
}

// Signed first mass moment from u_dot + z*(theta_dot cross normal).
fn accumulate_first_moment(m: &mut [f64; 324], n: [f64; 3], first: f64) {
    if first != 0.0 {
        for i in 0..3 {
            for j in 0..3 {
                let value = n[i] * n[j] * first;
                for (translation_axis, rotation_axis, sign) in [(0, 4, 1.0), (1, 3, -1.0)] {
                    let row = 6 * i + translation_axis;
                    let col = 6 * j + rotation_axis;
                    m[18 * row + col] += sign * value;
                    m[18 * col + row] += sign * value;
                }
            }
        }
    }
}

// Original numerical relative-spin gauge, kept separate by retained motion
// owners so its storage and integration losses cannot be labeled material heat.
fn accumulate_drilling_point(m: &mut [f64; 324], g: &FacetGeometry, n: [f64; 3], rotary: f64) {
    let slots = [0, 1, 5, 6, 7, 11, 12, 13, 17];
    let spin: [f64; 9] = std::array::from_fn(|i| match i % 3 {
        0 => 0.5 * g.gradient[i / 3][1],
        1 => -0.5 * g.gradient[i / 3][0],
        _ => n[i / 3],
    });
    for i in 0..9 {
        for j in i..9 {
            let value = rotary * spin[i] * spin[j];
            m[18 * slots[i] + slots[j]] += value;
            if i != j {
                m[18 * slots[j] + slots[i]] += value;
            }
        }
    }
}
