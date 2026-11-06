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

fn membrane_spin(geometry: &FacetGeometry, v: &[f64; 18]) -> f64 {
    let mut gradient = 0.0;
    for i in 1..3 {
        gradient += geometry.gradient[i][0] * (v[6 * i + 1] - v[1])
            - geometry.gradient[i][1] * (v[6 * i] - v[0]);
    }
    0.5 * gradient
}

fn point_velocity(n: [f64; 3], v: &[f64; 18]) -> [f64; 6] {
    core::array::from_fn(|c| v[c] + n[1] * (v[6 + c] - v[c]) + n[2] * (v[12 + c] - v[c]))
}

// Pivoted square completion of the same physical I0/I1/I2 matrix. Sharing
// amplitudes preserves small relative motion beside large common velocities.
struct PhysicalBasis {
    weights: [f64; 5],
    coupling: f64,
    translation_dominant: bool,
}
impl PhysicalBasis {
    fn new(translation: f64, first: f64, rotary: f64) -> Result<Self, PlateError> {
        let coupling = first / translation.max(rotary);
        let variance = (-first).mul_add(coupling, translation.min(rotary));
        if !variance.is_finite() || variance <= 0.0 {
            return Err(PlateError::BadSection {
                what: "shell kinetic mass variance must remain positive and representable",
            });
        }
        let pivot = translation.max(rotary);
        Ok(Self {
            weights: [translation, pivot, pivot, variance, variance],
            coupling,
            translation_dominant: translation >= rotary,
        })
    }

    fn amplitudes(&self, v: [f64; 6]) -> [f64; 5] {
        if self.translation_dominant {
            [
                v[2],
                v[0] + self.coupling * v[4],
                v[1] - self.coupling * v[3],
                v[3],
                v[4],
            ]
        } else {
            [
                v[2],
                v[4] + self.coupling * v[0],
                v[3] - self.coupling * v[1],
                v[0],
                v[1],
            ]
        }
    }
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
    let spin_a = membrane_spin(geometry, a);
    let spin_b = membrane_spin(geometry, b);
    let mut physical = 0.0;
    let mut drilling = 0.0;
    for_each_reference_point(|n, weight| {
        let (translation, first, rotary) =
            reference_point_moments(geometry, thickness, moments, n, weight)?;
        let av = point_velocity(n, a);
        let bv = point_velocity(n, b);
        let basis = PhysicalBasis::new(translation, first, rotary)?;
        let aa = basis.amplitudes(av);
        let ba = basis.amplitudes(bv);
        for i in 0..5 {
            physical += basis.weights[i] * aa[i] * ba[i];
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

/// One physical or numerical-gauge contribution to an inertial step, in joules.
/// Storage and work include the conventional factor of one half where required.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReferenceKineticIncrement {
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub kinetic_change_j: f64,
    /// Backward-Euler loss `0.5 * delta_velocity^T M delta_velocity`.
    pub numerical_loss_j: f64,
    /// Endpoint momentum work `velocity_after^T M delta_velocity`.
    pub inertial_work_j: f64,
}
impl ReferenceKineticIncrement {
    fn accumulate(&mut self, weight: f64, before: f64, after: f64) {
        let delta = after - before;
        let mean = 0.5 * after + 0.5 * before;
        let impulse = weight * delta;
        self.kinetic_before_j += 0.5 * before * (weight * before);
        self.kinetic_after_j += 0.5 * after * (weight * after);
        self.kinetic_change_j += mean * impulse;
        self.numerical_loss_j += 0.5 * delta * impulse;
        self.inertial_work_j += after * impulse;
    }

    fn validate(self) -> Result<Self, PlateError> {
        let scale = self
            .kinetic_before_j
            .max(self.kinetic_after_j)
            .max(self.numerical_loss_j)
            .max(self.inertial_work_j.abs());
        if [
            self.kinetic_before_j,
            self.kinetic_after_j,
            self.numerical_loss_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || !self.kinetic_change_j.is_finite()
            || !self.inertial_work_j.is_finite()
            || (self.inertial_work_j - self.kinetic_change_j - self.numerical_loss_j).abs()
                > 2048.0 * f64::EPSILON * scale
        {
            return Err(PlateError::BadSection {
                what: "shell native kinetic increment misses its work identity",
            });
        }
        Ok(self)
    }
}

/// Integrate one backward-Euler kinetic-work account in one native traversal.
///
/// The first result belongs to physical midsurface/director motion; the second
/// belongs to the original numerical relative-drilling gauge. They use the
/// same inventory and local physical velocity ordering as
/// [`reference_kinetic_bilinear`]. Endpoint, mean and difference products are
/// formed from the same original point amplitudes. Nearly rigid velocities
/// therefore preserve the tiny gauge account's work identity without imposing
/// a large absolute tolerance or relabeling a roundoff deficit as heat.
///
/// No matrices or heap allocations are needed. These are reference kinetic
/// integrals only; material dissipation and external work retain their owners.
///
/// # Errors
/// Invalid reference data or velocities, unrepresentable positive mass
/// variance, nonfinite energies, or failure of the native kinetic-work identity.
pub fn reference_kinetic_increment(
    geometry: &FacetGeometry,
    thickness: [f64; 3],
    moments: [ReferenceDensityMoments; 3],
    before: &[f64; 18],
    after: &[f64; 18],
) -> Result<(ReferenceKineticIncrement, ReferenceKineticIncrement), PlateError> {
    validate(geometry, thickness, moments)?;
    if before.iter().chain(after).any(|x| !x.is_finite()) {
        return Err(PlateError::BadSection {
            what: "shell kinetic increments require finite physical velocities",
        });
    }
    let spin_before = membrane_spin(geometry, before);
    let spin_after = membrane_spin(geometry, after);
    let mut physical = ReferenceKineticIncrement::default();
    let mut drilling = ReferenceKineticIncrement::default();
    for_each_reference_point(|n, weight| {
        let (translation, first, rotary) =
            reference_point_moments(geometry, thickness, moments, n, weight)?;
        let av = point_velocity(n, before);
        let bv = point_velocity(n, after);
        let basis = PhysicalBasis::new(translation, first, rotary)?;
        let aa = basis.amplitudes(av);
        let ba = basis.amplitudes(bv);
        for i in 0..5 {
            physical.accumulate(basis.weights[i], aa[i], ba[i]);
        }
        drilling.accumulate(rotary, av[5] - spin_before, bv[5] - spin_after);
        Ok(())
    })?;
    Ok((physical.validate()?, drilling.validate()?))
}
