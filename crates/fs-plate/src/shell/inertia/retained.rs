//! The original shell kinetic field without an elastic constitutive placeholder.

use super::{
    FacetGeometry, PlateError, ReferenceDensityMoments, accumulate_drilling_point,
    accumulate_first_moment, accumulate_physical_point, for_each_reference_point,
    reference_mass_kg, reference_point_moments,
};

/// Original local shell reference inertia, with numerical drilling kept explicit.
///
/// Both row-major 18-by-18 matrices use `[u,v,w,theta_x,theta_y,theta_z]` at
/// each vertex in the facet frame. Translational coordinates are metres and
/// physical axial rotations are radians. A caller rotates both vector blocks
/// with the original facet frame and applies any numerical coordinate scaling
/// by congruence. The physical and drilling matrices sum to the native
/// consistent mass; neither matrix is an elastic material property.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellReferenceMass {
    /// P1 midsurface/director kinetic integral, including signed first-moment
    /// translation/rotation coupling. Normal director spin has no physical mass.
    pub physical: [f64; 324],
    /// Original numerical relative drilling-rate gauge. Its kinetic storage
    /// and time-discretization losses must remain separate from material heat.
    pub drilling: [f64; 324],
    /// Actual reference inventory in kg; numerical drilling adds no material.
    pub reference_mass_kg: f64,
}

/// Prepare the original reference kinetic integral independently of stiffness.
///
/// `geometry` is the native geometry returned by `ShellMesh::facet`. Positive
/// thickness and the admitted normalized density moments are P1 vertex fields.
/// At each point their physical moments are `h*q0`, `h²*q1`, and `h³*q2` about
/// that same reference surface. Signed first moments therefore retain eccentric
/// and layered inventories. No effective uniform density or surrogate elastic
/// section is introduced. The existing positive degree-six quadrature and
/// physical/director/gauge kernels own every entry.
///
/// Work and storage are fixed-size for one facet. A mesh owner supplies its own
/// aggregate budget and cancellation checkpoints between facets. Retained
/// motion must keep the original inventory; changing it needs an independent
/// matter/momentum-transfer law.
///
/// # Errors
/// Nonfinite native geometry, nonpositive/nonfinite thickness, or unrepresentable
/// physical mass, rotary inertia or matrix entries.
pub fn prepare_reference_mass(
    geometry: &FacetGeometry,
    thickness: [f64; 3],
    moments: [ReferenceDensityMoments; 3],
) -> Result<ShellReferenceMass, PlateError> {
    let reference_mass_kg = validate(geometry, thickness, moments)?;
    let mut physical = [0.0; 324];
    let mut drilling = [0.0; 324];
    for_each_reference_point(|n, weight| {
        let (translation, first, rotary) =
            reference_point_moments(geometry, thickness, moments, n, weight)?;
        accumulate_physical_point(&mut physical, n, translation, rotary);
        accumulate_first_moment(&mut physical, n, first);
        accumulate_drilling_point(&mut drilling, geometry, n, rotary);
        Ok(())
    })?;
    if physical.iter().chain(&drilling).any(|x| !x.is_finite()) {
        return Err(PlateError::BadSection {
            what: "shell reference inertia matrix overflow",
        });
    }
    Ok(ShellReferenceMass {
        physical,
        drilling,
        reference_mass_kg,
    })
}

fn validate(
    geometry: &FacetGeometry,
    thickness: [f64; 3],
    moments: [ReferenceDensityMoments; 3],
) -> Result<f64, PlateError> {
    if !geometry.area_m2.is_finite()
        || geometry.area_m2 <= 0.0
        || geometry.frame.iter().flatten().any(|x| !x.is_finite())
        || geometry.gradient.iter().flatten().any(|x| !x.is_finite())
        || geometry.x.iter().chain(&geometry.y).any(|x| !x.is_finite())
        || thickness.iter().any(|h| !h.is_finite() || *h <= 0.0)
    {
        return Err(PlateError::BadSection {
            what: "shell reference inertia needs finite native geometry and positive thickness",
        });
    }
    reference_mass_kg(
        geometry.area_m2,
        thickness,
        moments.map(ReferenceDensityMoments::zeroth_kg_m3),
    )
}

/// Evaluate `a^T M_physical b` and `a^T M_drilling b` without building matrices.
///
/// Inputs have the same local physical ordering as [`ShellReferenceMass`].
/// Neither returned value includes a factor of one half. This is the same
/// reference quadrature, density moments and kinetic basis used by
/// [`prepare_reference_mass`], with no subtraction of total and gauge energies.
/// Each coupled physical mass block is evaluated by square completion about its
/// larger diagonal. A quadratic query therefore retains nonnegative physical
/// storage even for eccentric layers and almost pure drilling motion.
///
/// # Errors
/// The same reference-inventory failures as [`prepare_reference_mass`],
/// nonfinite input vectors or unrepresentable kinetic products/positive variance.
pub fn reference_kinetic_bilinear(
    geometry: &FacetGeometry,
    thickness: [f64; 3],
    moments: [ReferenceDensityMoments; 3],
    a: &[f64; 18],
    b: &[f64; 18],
) -> Result<(f64, f64), PlateError> {
    validate(geometry, thickness, moments)?;
    if a.iter().chain(b).any(|x| !x.is_finite()) {
        return Err(PlateError::BadSection {
            what: "shell kinetic products require finite physical coordinates",
        });
    }
    let spin = |v: &[f64; 18]| {
        let mut gradient = 0.0;
        for i in 1..3 {
            gradient += geometry.gradient[i][0] * (v[6 * i + 1] - v[1])
                - geometry.gradient[i][1] * (v[6 * i] - v[0]);
        }
        0.5 * gradient
    };
    let spin_a = spin(a);
    let spin_b = spin(b);
    let mut physical = 0.0;
    let mut drilling = 0.0;
    for_each_reference_point(|n, weight| {
        let (translation, first, rotary) =
            reference_point_moments(geometry, thickness, moments, n, weight)?;
        let at = |v: &[f64; 18]| -> [f64; 6] {
            core::array::from_fn(|c| v[c] + n[1] * (v[6 + c] - v[c]) + n[2] * (v[12 + c] - v[c]))
        };
        let av = at(a);
        let bv = at(b);
        let coupling = first / translation.max(rotary);
        let variance = (-first).mul_add(coupling, translation.min(rotary));
        if !variance.is_finite() || variance <= 0.0 {
            return Err(PlateError::BadSection {
                what: "shell kinetic mass variance must remain positive and representable",
            });
        }
        physical += translation * av[2] * bv[2];
        if translation >= rotary {
            physical += translation
                * ((av[0] + coupling * av[4]) * (bv[0] + coupling * bv[4])
                    + (av[1] - coupling * av[3]) * (bv[1] - coupling * bv[3]))
                + variance * (av[3] * bv[3] + av[4] * bv[4]);
        } else {
            physical += rotary
                * ((av[4] + coupling * av[0]) * (bv[4] + coupling * bv[0])
                    + (av[3] - coupling * av[1]) * (bv[3] - coupling * bv[1]))
                + variance * (av[0] * bv[0] + av[1] * bv[1]);
        }
        drilling += rotary * (av[5] - spin_a) * (bv[5] - spin_b);
        Ok(())
    })?;
    if !physical.is_finite() || !drilling.is_finite() {
        return Err(PlateError::BadSection {
            what: "shell kinetic product overflow",
        });
    }
    Ok((physical, drilling))
}
