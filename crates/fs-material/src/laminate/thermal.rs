//! Prescribed thermal expansion of the original, ordered elastic laminate.
//!
//! Each layer has its own integrated free strain, affine between its faces.
//! Thermal force, thermal moment and elastic storage are thickness integrals
//! of the retained layer stiffness, not an averaged expansion coefficient.
//! Two-point Gauss integration is exact for these piecewise affine strains.
//! Energy is evaluated from elastic strain before squaring, preserving the
//! residual self-stress energy of a freely bending incompatible laminate.
//!
//! This is infinitesimal-strain, prescribed-temperature elasticity with fixed
//! section stiffness. It does not evolve heat, layer interfaces or viscoelastic
//! history. Source integration and temperature-domain admission remain with
//! the material state-point owners; callers retain those original receipts.

use super::{LaminateSection, LayerLaw};
use crate::state_point::{
    IntegratedIsotropicThermalExpansion, IntegratedOrthotropicThermalExpansion,
};
use fs_math::det;
use std::sync::Arc;

/// Refusal of one complete preparation or response; no partial state is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminateThermalError {
    /// Exactly one profile is required per original layer.
    LayerCount { expected: usize, actual: usize },
    /// The requested recovery layer does not exist.
    LayerIndex { index: usize, count: usize },
    /// Finite inputs and representable physical integrals are required.
    Nonfinite { quantity: &'static str },
    /// A thermal history law is required for this viscoelastic layer.
    ViscoelasticLayer { layer: usize },
    /// The caller cancelled before the complete result was available.
    Cancelled,
}

impl core::fmt::Display for LaminateThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-THERMAL: {self:?}")
    }
}
impl std::error::Error for LaminateThermalError {}

/// Integrated engineering free strain `[epsilon_11, epsilon_22, gamma_12]`
/// at the two faces of one layer, in that layer's original principal frame.
/// Values are strains, never CTEs. Axis 3 is explicitly normal to the section.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminaThermalStrain {
    bottom: [f64; 3],
    top: [f64; 3],
}

impl LaminaThermalStrain {
    /// A constant already integrated free strain throughout this layer.
    pub fn uniform(strain: [f64; 3]) -> Result<Self, LaminateThermalError> {
        Self::from_faces(strain, strain)
    }

    /// Declare an affine free-strain profile between the original layer faces.
    /// A nonlinear temperature/CTE profile must first be integrated into a
    /// suitable layer discretization; no nonlinear profile is inferred here.
    pub fn from_faces(bottom: [f64; 3], top: [f64; 3]) -> Result<Self, LaminateThermalError> {
        finite(bottom.iter().chain(&top), "layer free strain")?;
        Ok(Self { bottom, top })
    }

    /// Lower the material owner's admitted integrated isotropic strain.
    /// Retain `path` for its original card, temperature path and source receipt.
    pub fn from_isotropic(
        path: &IntegratedIsotropicThermalExpansion,
    ) -> Result<Self, LaminateThermalError> {
        let e = path.free_linear_strain();
        Self::uniform([e, e, 0.0])
    }

    /// Lower admitted principal strains with axes 1/2 in the layer plane and
    /// axis 3 normal. `LaminaLayer::orientation_rad` supplies the original
    /// in-plane orientation. Retain `path` for its source identity and domain.
    pub fn from_orthotropic(
        path: &IntegratedOrthotropicThermalExpansion,
    ) -> Result<Self, LaminateThermalError> {
        let [a, b, _] = *path.free_principal_strains();
        Self::uniform([a, b, 0.0])
    }

    /// Bottom and top free strains in the original layer principal frame.
    pub const fn faces(&self) -> [[f64; 3]; 2] {
        [self.bottom, self.top]
    }
}

/// In-plane response per original midsurface area, about the declared reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateThermalResponse {
    /// Physical `[Nxx, Nyy, Nxy]` [N/m], including thermal restraint stress.
    pub membrane_force_n_m: [f64; 3],
    /// Physical `[Mxx, Myy, Mxy]` [N], with signed reference-surface moments.
    pub bending_moment_n: [f64; 3],
    /// Physical membrane/bending elastic storage [J/m²], including self-stress.
    /// Transverse shear storage belongs to the surrounding shell operator.
    pub strain_energy_j_m2: f64,
}

/// Physical in-plane stress `[sigma_xx, sigma_yy, tau_xy]` in section axes.
/// The two stresses at an interface belong to their respective layers and
/// need not agree; this does not recover interlaminar traction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerThermalStress {
    /// Stress at this layer's bottom face [Pa].
    pub bottom_pa: [f64; 3],
    /// Stress at this layer's top face [Pa].
    pub top_pa: [f64; 3],
}

#[derive(Clone, Debug, PartialEq)]
struct LayerSamples {
    q: [f64; 9],
    z: [f64; 2],
    free: [[f64; 3]; 2],
    face_free: [[f64; 3]; 2],
    weight_m: f64,
}

/// Immutable thermal preparation sharing the exact original section and its
/// retained material sources. Geometry, mass, stiffness and layer order remain
/// unchanged. Preparation is O(layers); evaluation allocates nothing.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedLaminateThermal {
    section: Arc<LaminateSection>,
    profiles: Vec<LaminaThermalStrain>,
    samples: Vec<LayerSamples>,
    force: [f64; 3],
    moment: [f64; 3],
    restrained_energy: f64,
}

impl PreparedLaminateThermal {
    /// Integrate original per-layer free strains using original rotated elastic
    /// stiffness and thickness coordinates. A Prony layer refuses because its
    /// band-reference stiffness is not a thermal-history constitutive law.
    /// Cancellation is polled before each layer and after complete admission.
    pub fn try_new(
        section: Arc<LaminateSection>,
        profiles: &[LaminaThermalStrain],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateThermalError> {
        if cancelled() {
            return Err(LaminateThermalError::Cancelled);
        }
        if profiles.len() != section.layers().len() {
            return Err(LaminateThermalError::LayerCount {
                expected: section.layers().len(),
                actual: profiles.len(),
            });
        }
        let mut samples = Vec::with_capacity(profiles.len());
        let mut force = [0.0; 3];
        let mut moment = [0.0; 3];
        let mut restrained_energy = 0.0;
        // Weights 1/2 and positions (1 +/- 1/sqrt(3))/2 on each layer.
        let t = [0.211_324_865_405_187_13, 0.788_675_134_594_812_9];
        for (i, ((layer, integral), profile)) in section
            .layers()
            .iter()
            .zip(section.layer_integrals())
            .zip(profiles)
            .enumerate()
        {
            if cancelled() {
                return Err(LaminateThermalError::Cancelled);
            }
            if matches!(&layer.law, LayerLaw::Prony(_)) {
                return Err(LaminateThermalError::ViscoelasticLayer { layer: i });
            }
            let faces = profile
                .faces()
                .map(|e| rotate_strain(e, layer.orientation_rad));
            let h = integral.z_top_m - integral.z_bottom_m;
            if !h.is_finite() || h <= 0.0 {
                return Err(LaminateThermalError::Nonfinite {
                    quantity: "representable positive layer thickness",
                });
            }
            let sample = LayerSamples {
                q: integral.q_bar_pa,
                z: t.map(|u| (1.0 - u) * integral.z_bottom_m + u * integral.z_top_m),
                free: t
                    .map(|u| core::array::from_fn(|j| (1.0 - u) * faces[0][j] + u * faces[1][j])),
                face_free: faces,
                weight_m: 0.5 * h,
            };
            for p in 0..2 {
                let stress = matvec(&sample.q, sample.free[p]);
                for j in 0..3 {
                    force[j] += sample.weight_m * stress[j];
                    moment[j] += sample.weight_m * sample.z[p] * stress[j];
                }
                restrained_energy += 0.5 * sample.weight_m * dot(sample.free[p], stress);
            }
            finite(
                sample
                    .free
                    .iter()
                    .flatten()
                    .chain(sample.face_free.iter().flatten())
                    .chain(&sample.z),
                "rotated layer thermal profile",
            )?;
            samples.push(sample);
        }
        finite(
            force
                .iter()
                .chain(&moment)
                .chain(core::iter::once(&restrained_energy)),
            "thermal section integrals",
        )?;
        if cancelled() {
            return Err(LaminateThermalError::Cancelled);
        }
        Ok(Self {
            section,
            profiles: profiles.to_vec(),
            samples,
            force,
            moment,
            restrained_energy,
        })
    }

    /// Exact immutable section, including layer sources, orientation and inertia.
    pub fn section(&self) -> &LaminateSection {
        &self.section
    }
    /// Original integrated profiles in original layer order and principal axes.
    pub fn profiles(&self) -> &[LaminaThermalStrain] {
        &self.profiles
    }
    /// Equivalent thermal membrane loading `[Nxx, Nyy, Nxy]` [N/m].
    pub const fn membrane_force_n_m(&self) -> [f64; 3] {
        self.force
    }
    /// Equivalent thermal bending loading `[Mxx, Myy, Mxy]` [N].
    pub const fn bending_moment_n(&self) -> [f64; 3] {
        self.moment
    }
    /// Elastic energy per area at zero membrane strain and zero curvature.
    pub const fn restrained_energy_j_m2(&self) -> f64 {
        self.restrained_energy
    }

    /// Recover membrane/bending resultants and actual elastic storage for
    /// `epsilon(z) = membrane + z * gradient_per_m` in the section axes.
    /// Exact two-point integration retains incompatible layer self-stress.
    /// No allocation, source lookup or trigonometry occurs in this evaluation.
    pub fn evaluate(
        &self,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateThermalResponse, LaminateThermalError> {
        finite(membrane.iter().chain(&gradient_per_m), "section strain")?;
        let mut response = LaminateThermalResponse {
            membrane_force_n_m: [0.0; 3],
            bending_moment_n: [0.0; 3],
            strain_energy_j_m2: 0.0,
        };
        for sample in &self.samples {
            if cancelled() {
                return Err(LaminateThermalError::Cancelled);
            }
            for p in 0..2 {
                let elastic = core::array::from_fn(|j| {
                    membrane[j] + sample.z[p] * gradient_per_m[j] - sample.free[p][j]
                });
                let stress = matvec(&sample.q, elastic);
                for j in 0..3 {
                    response.membrane_force_n_m[j] += sample.weight_m * stress[j];
                    response.bending_moment_n[j] += sample.weight_m * sample.z[p] * stress[j];
                }
                response.strain_energy_j_m2 += 0.5 * sample.weight_m * dot(elastic, stress);
            }
        }
        finite(
            response
                .membrane_force_n_m
                .iter()
                .chain(&response.bending_moment_n)
                .chain(core::iter::once(&response.strain_energy_j_m2)),
            "thermal section response",
        )?;
        if cancelled() {
            return Err(LaminateThermalError::Cancelled);
        }
        Ok(response)
    }

    /// Recover physical face stress in section axes for one original layer.
    /// This fixed-size O(1) operation allocates nothing.
    pub fn layer_stress(
        &self,
        layer: usize,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerThermalStress, LaminateThermalError> {
        finite(membrane.iter().chain(&gradient_per_m), "section strain")?;
        let sample = self
            .samples
            .get(layer)
            .ok_or(LaminateThermalError::LayerIndex {
                index: layer,
                count: self.samples.len(),
            })?;
        let integral = &self.section.layer_integrals()[layer];
        let at = |p: usize, z: f64| {
            matvec(
                &sample.q,
                core::array::from_fn(|j| {
                    membrane[j] + z * gradient_per_m[j] - sample.face_free[p][j]
                }),
            )
        };
        let bottom_pa = at(0, integral.z_bottom_m);
        let top_pa = at(1, integral.z_top_m);
        finite(bottom_pa.iter().chain(&top_pa), "layer face stress")?;
        Ok(LayerThermalStress { bottom_pa, top_pa })
    }
}

fn finite<'a>(
    values: impl Iterator<Item = &'a f64>,
    quantity: &'static str,
) -> Result<(), LaminateThermalError> {
    if values.into_iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(LaminateThermalError::Nonfinite { quantity })
    }
}

fn matvec(q: &[f64; 9], strain: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| q[3 * i + j] * strain[j]).sum())
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn rotate_strain(e: [f64; 3], angle: f64) -> [f64; 3] {
    if angle == 0.0 {
        return e;
    }
    let (s, c) = (det::sin(angle), det::cos(angle));
    [
        c * c * e[0] + s * s * e[1] - c * s * e[2],
        s * s * e[0] + c * c * e[1] + c * s * e[2],
        2.0 * c * s * (e[0] - e[1]) + (c * c - s * s) * e[2],
    ]
}
