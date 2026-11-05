//! A held stress-free deformation, shared by elastic and history-bearing laws.
//!
//! For total deformation `F`, prescribed reference deformation `G` and
//! `H = G^-1`, the elastic deformation is `Fe = F H`. Constitutive energies are
//! per relaxed volume: the original-volume potential is `W0(F) = det(G) We(Fe)`.
//! Its first Piola stress and Hessian follow the same chain rule. The reference
//! tensors are prepared once; point evaluations allocate and invert nothing.
//!
//! This is a kinematic owner, not an expansion law. Source-derived engineering
//! length changes, a growth law or a prescribed prestrain must supply `G` in a
//! named material frame. Mass and actual volume remain owned by total `F` and
//! the original geometry. Changing `G` with retained viscous/plastic history
//! requires a separate, explicit history and external-work transformation.

use super::Hyperelastic;
use super::cards::HyperelasticMaterial;
use super::mixed::{
    MixedFiniteStrainPotential, PressurePotential, RemainderResponse, VolumetricTreatment, cofactor,
};
use crate::MaterialError;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn parameters(what: &str) -> MaterialError {
    MaterialError::Parameters {
        what: format!("reference deformation: {what}"),
    }
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("reference deformation: {what}"),
    }
}

/// Immutable, orientation-preserving stress-free deformation in material axes.
///
/// `G` may include anisotropic expansion and shear. The determinant, inverse
/// and cofactor must all be finite and representable. No projection, clamping
/// or inference of a new material frame occurs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReferenceDeformation {
    gradient: [f64; 9],
    inverse: [f64; 9],
    cofactor: [f64; 9],
    volume_ratio: f64,
    isotropic_stretch: Option<f64>,
}

impl Default for ReferenceDeformation {
    fn default() -> Self {
        Self::identity()
    }
}

impl ReferenceDeformation {
    /// Unchanged material reference; its constitutive operations are exact
    /// pass-throughs, including floating-point evaluation order.
    pub const fn identity() -> Self {
        Self {
            gradient: IDENTITY,
            inverse: IDENTITY,
            cofactor: IDENTITY,
            volume_ratio: 1.0,
            isotropic_stretch: Some(1.0),
        }
    }

    /// Prepare a row-major reference deformation. Inverted, singular and
    /// nonfinite maps are refused before any constitutive history is created.
    pub fn try_new(gradient: [f64; 9]) -> Result<Self, MaterialError> {
        if gradient.iter().any(|x| !x.is_finite()) {
            return Err(parameters("finite gradient required"));
        }
        if gradient == IDENTITY {
            return Ok(Self::identity());
        }
        let cofactor = cofactor(&gradient);
        let volume_ratio =
            gradient[0] * cofactor[0] + gradient[1] * cofactor[1] + gradient[2] * cofactor[2];
        if !(volume_ratio.is_finite() && volume_ratio > 0.0) {
            return Err(parameters("finite positive reference volume required"));
        }
        let inverse = core::array::from_fn(|ij| cofactor[3 * (ij % 3) + ij / 3] / volume_ratio);
        if inverse.iter().chain(&cofactor).any(|x| !x.is_finite()) {
            return Err(parameters("inverse or area ratio is unrepresentable"));
        }
        let isotropic_stretch = if gradient[0] > 0.0
            && gradient[0] == gradient[4]
            && gradient[0] == gradient[8]
            && [1, 2, 3, 5, 6, 7].iter().all(|i| gradient[*i] == 0.0)
        {
            Some(gradient[0])
        } else {
            None
        };
        Ok(Self {
            gradient,
            inverse,
            cofactor,
            volume_ratio,
            isotropic_stretch,
        })
    }

    /// Prepare positive principal length ratios in an explicitly supplied
    /// frame. Columns of `principal_to_reference` are the principal axes.
    /// The frame must be proper and orthonormal within 1e-12. The inverse and
    /// areas use the original principal lengths directly, without reinversion
    /// of a rotated tensor. Length ratios retain their source convention.
    pub fn from_principal_stretches(
        stretches: [f64; 3],
        principal_to_reference: [[f64; 3]; 3],
    ) -> Result<Self, MaterialError> {
        let q = principal_to_reference;
        if stretches.iter().any(|s| !s.is_finite() || *s <= 0.0)
            || q.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(parameters(
                "positive finite lengths and finite axes required",
            ));
        }
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| q[k][i] * q[k][j]).sum();
                if (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1.0e-12 {
                    return Err(parameters("orthonormal principal axes required"));
                }
            }
        }
        let flat = core::array::from_fn(|ij| q[ij / 3][ij % 3]);
        if (super::det3(&flat) - 1.0).abs() > 1.0e-12 {
            return Err(parameters("proper principal rotation required"));
        }
        if stretches[0] == stretches[1] && stretches[0] == stretches[2] {
            let s = stretches[0];
            return Self::try_new([s, 0.0, 0.0, 0.0, s, 0.0, 0.0, 0.0, s]);
        }
        let reciprocal = stretches.map(|s| 1.0 / s);
        let areas = [
            stretches[1] * stretches[2],
            stretches[0] * stretches[2],
            stretches[0] * stretches[1],
        ];
        let volume_ratio = stretches[0] * areas[0];
        if !(volume_ratio.is_finite() && volume_ratio > 0.0)
            || reciprocal
                .iter()
                .chain(&areas)
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(parameters(
                "reference volume, inverse or area ratio is unrepresentable",
            ));
        }
        let rotate = |values: [f64; 3]| -> [f64; 9] {
            core::array::from_fn(|ij| {
                let (i, j) = (ij / 3, ij % 3);
                (0..3).map(|a| (q[i][a] * values[a]) * q[j][a]).sum()
            })
        };
        let gradient = rotate(stretches);
        let inverse = rotate(reciprocal);
        let cofactor = rotate(areas);
        if gradient
            .iter()
            .chain(&inverse)
            .chain(&cofactor)
            .any(|x| !x.is_finite())
        {
            return Err(parameters("oriented reference tensors are unrepresentable"));
        }
        Ok(Self {
            gradient,
            inverse,
            cofactor,
            volume_ratio,
            isotropic_stretch: None,
        })
    }

    /// Row-major stress-free `G`, independent of actual compatible `F`.
    pub const fn deformation_gradient(&self) -> &[f64; 9] {
        &self.gradient
    }

    /// Relaxed volume per original reference volume, `det(G)`.
    pub const fn volume_ratio(&self) -> f64 {
        self.volume_ratio
    }

    /// Map total `F` to elastic `F G^-1`. Near the relaxed configuration,
    /// evaluate `I + (F-G) G^-1` to retain weak mechanical increments and
    /// return identity exactly at `F=G`. Away from it, direct multiplication
    /// avoids subtracting and restoring large reference entries.
    pub fn elastic_deformation(&self, total: &[f64; 9]) -> [f64; 9] {
        if let Some(s) = self.isotropic_stretch {
            return if s == 1.0 {
                *total
            } else {
                total.map(|x| x / s)
            };
        }
        let increment: [f64; 9] = core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            (0..3)
                .map(|a| (total[3 * i + a] - self.gradient[3 * i + a]) * self.inverse[3 * a + j])
                .sum()
        });
        if increment.iter().all(|x| x.abs() <= 0.25) {
            return core::array::from_fn(|ij| IDENTITY[ij] + increment[ij]);
        }
        core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            (0..3)
                .map(|a| total[3 * i + a] * self.inverse[3 * a + j])
                .sum()
        })
    }

    /// Original-volume first Piola stress, `det(G) Pe G^-T`.
    pub fn pullback_piola(&self, elastic: &[f64; 9]) -> [f64; 9] {
        if let Some(s) = self.isotropic_stretch {
            return if s == 1.0 {
                *elastic
            } else {
                elastic.map(|x| (s * s) * x)
            };
        }
        core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            (0..3)
                .map(|a| elastic[3 * i + a] * self.cofactor[3 * a + j])
                .sum()
        })
    }

    /// Original-volume constitutive Hessian. Two three-term contractions
    /// replace a nine-term contraction for every entry; no temporary heap
    /// storage, numerical differentiation or repeated inversions are needed.
    pub fn pullback_tangent(&self, elastic: &[[f64; 9]; 9]) -> [[f64; 9]; 9] {
        if let Some(s) = self.isotropic_stretch {
            return if s == 1.0 {
                *elastic
            } else {
                elastic.map(|row| row.map(|x| s * x))
            };
        }
        let right: [[f64; 9]; 9] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                let (k, l) = (column / 3, column % 3);
                (0..3)
                    .map(|b| elastic[row][3 * k + b] * self.inverse[3 * l + b])
                    .sum()
            })
        });
        core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                let (i, j) = (row / 3, row % 3);
                (0..3)
                    .map(|a| self.cofactor[3 * a + j] * right[3 * i + a][column])
                    .sum()
            })
        })
    }

    /// Transform the local pressure dual into the total-volume convention.
    /// Pressure remains physical tensile pressure: the mixed residual is
    /// `det(F) - det(G) Je(p)` and its cross derivative is `cof(F)`.
    /// Exact incompressibility therefore constrains elastic, not total, volume.
    pub fn pullback_pressure_potential(
        &self,
        elastic: PressurePotential,
        pressure_pa: f64,
    ) -> Result<PressurePotential, MaterialError> {
        let j = self.volume_ratio;
        let value = if j == 1.0 {
            elastic
        } else {
            PressurePotential {
                complementary_pa: j * elastic.complementary_pa + pressure_pa * (j - 1.0),
                bulk_stored_pa: j * elastic.bulk_stored_pa,
                volume_ratio: j * elastic.volume_ratio,
                volume_increment: (j - 1.0) + j * elastic.volume_increment,
                compliance_per_pa: j * elastic.compliance_per_pa,
            }
        };
        if !pressure_pa.is_finite()
            || ![
                value.complementary_pa,
                value.bulk_stored_pa,
                value.volume_ratio,
                value.volume_increment,
                value.compliance_per_pa,
            ]
            .iter()
            .all(|x| x.is_finite())
            || value.volume_ratio <= 0.0
            || value.compliance_per_pa < 0.0
        {
            return Err(state_error("pressure pullback is unrepresentable"));
        }
        Ok(value)
    }
}

/// A borrowed elastic material evaluated about a prescribed relaxed reference.
///
/// The original material, source card and constitutive limits are retained.
/// Limits apply to elastic `F G^-1`, while native geometry/recovery still see
/// total `F`. For a history-bearing law, select its own reference-map API so
/// all accepted states and work/dissipation ledgers retain the same reference.
#[derive(Clone, Copy, Debug)]
pub struct ReferencedHyperelastic<'a, L: HyperelasticMaterial + ?Sized> {
    material: &'a L,
    reference: ReferenceDeformation,
}

impl<'a, L: HyperelasticMaterial + ?Sized> ReferencedHyperelastic<'a, L> {
    /// Retain an admitted material and an immutable prescribed reference.
    pub const fn new(material: &'a L, reference: ReferenceDeformation) -> Self {
        Self {
            material,
            reference,
        }
    }

    /// Original material, including any source-backed constitutive evidence.
    pub const fn material(&self) -> &'a L {
        self.material
    }

    /// The held reference used by every constitutive operation.
    pub const fn reference_deformation(&self) -> &ReferenceDeformation {
        &self.reference
    }
}

impl<L: HyperelasticMaterial + ?Sized> HyperelasticMaterial for ReferencedHyperelastic<'_, L> {
    fn hyperelastic(&self) -> &Hyperelastic {
        self.material.hyperelastic()
    }

    fn energy_density(&self, total: &[f64; 9]) -> f64 {
        self.reference.volume_ratio()
            * self
                .material
                .energy_density(&self.reference.elastic_deformation(total))
    }

    fn first_piola(&self, total: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        let elastic = self
            .material
            .first_piola(&self.reference.elastic_deformation(total))?;
        let result = self.reference.pullback_piola(&elastic);
        if result.iter().any(|x| !x.is_finite()) {
            return Err(state_error("first Piola pullback is unrepresentable"));
        }
        Ok(result)
    }

    fn consistent_tangent(&self, total: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        let elastic = self
            .material
            .consistent_tangent(&self.reference.elastic_deformation(total))?;
        let result = self.reference.pullback_tangent(&elastic);
        if result.iter().flatten().any(|x| !x.is_finite()) {
            return Err(state_error("tangent pullback is unrepresentable"));
        }
        Ok(result)
    }
}

impl<L: HyperelasticMaterial + MixedFiniteStrainPotential + ?Sized> MixedFiniteStrainPotential
    for ReferencedHyperelastic<'_, L>
{
    fn equilibrium_law(&self) -> &Hyperelastic {
        self.material.equilibrium_law()
    }

    fn validate_mixed(&self) -> Result<(), MaterialError> {
        self.material.validate_mixed()
    }

    fn remainder_stress_and_potential(
        &self,
        total: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        let (piola, energy) = self
            .material
            .remainder_stress_and_potential(&self.reference.elastic_deformation(total))?;
        super::mixed::checked_stress(
            self.reference.volume_ratio() * energy,
            self.reference.pullback_piola(&piola),
        )
    }

    fn remainder_response(&self, total: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        let elastic = self
            .material
            .remainder_response(&self.reference.elastic_deformation(total))?;
        let result = RemainderResponse {
            potential_pa: self.reference.volume_ratio() * elastic.potential_pa,
            first_piola_pa: self.reference.pullback_piola(&elastic.first_piola_pa),
            tangent_pa: self.reference.pullback_tangent(&elastic.tangent_pa),
        };
        super::mixed::checked_stress(result.potential_pa, result.first_piola_pa)?;
        if result.tangent_pa.iter().flatten().any(|x| !x.is_finite()) {
            return Err(state_error("remainder tangent pullback is unrepresentable"));
        }
        Ok(result)
    }

    fn mixed_pressure_potential(
        &self,
        total: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        let elastic = self.material.mixed_pressure_potential(
            &self.reference.elastic_deformation(total),
            pressure_pa,
            treatment,
        )?;
        self.reference
            .pullback_pressure_potential(elastic, pressure_pa)
    }
}
