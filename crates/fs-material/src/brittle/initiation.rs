//! Tensile crack-plane discovery before handing an elastic increment to the
//! original fixed-plane crack band.
//!
//! The declared Rankine criterion is `lambda_max(stress) = tensile_strength`.
//! Stress is the ORIGINAL isotropic elastic stress, including Poisson coupling
//! and tensor shear. The maximum eigenvalue is convex on an affine strain path:
//! two strictly elastic endpoints cannot hide an interior tensile crossing.
//! A possible crossing is located with bounded bisection. Positive-definiteness
//! tests screen ordinary elastic points without allocation; unresolved points
//! and the final plane use the existing `fs_la` symmetric eigensolver.
//! [`TensileInitiator::activation_at`] instead observes one supplied strain under
//! its current coefficients, without inferring a path, crossing time or seed.
//!
//! Repeated principal tension has no unique physical crack direction. An
//! explicitly supplied material-direction hint is projected into that principal
//! eigenspace; absent or perpendicular hints refuse instead of choosing a world
//! axis. A simple principal plane is objective up to its immaterial sign.
//!
//! Handoff occurs at the ORIGINAL elastic increment start, followed by replay
//! through [`super::BrittleCrackBand::trial`]. This preserves the actual affine
//! path without rounding an onset strain into an already damaged state. Native
//! mixed-mode onset at that start, incompatible coefficients/plane, or a storage
//! discontinuity refuses. The virtual preload used to initialize native history
//! is not new external work. No accepted history is modified.
//!
//! This is a tensile initiation model, not an all-plane mixed-mode criterion,
//! shear/compressive fracture, crack-path solver, event-aware FE integrator,
//! topology change, or experimental calibration. After handoff the normal stays
//! fixed and the ORIGINAL mixed-mode law owns all subsequent evolution. Caller
//! owners must resolve reversals/non-affine paths and apply the detector only to
//! material that has remained elastic, never reset an existing damaged history.

mod activation;
pub use activation::TensileActivation;

use super::{
    BrittleCrackBand, BrittleError, BrittleState, finite, finite_tangent, finite_voigt, normalize,
    positive,
};
use crate::elastic::IsotropicElastic;
use crate::tensor::contract;
use crate::{Tangent6, Voigt};
use fs_math::det;

/// Numerical event resolution, separate from constitutive material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensileInitiationSettings {
    /// Frobenius norm of total tensor strain; at most 0.05 and the source limit.
    pub maximum_strain: f64,
    /// Absolute stress tolerance is this value times the tensile strength.
    /// Also resolves numerically repeated principal stresses. In [1024 eps, .01].
    pub relative_stress_tolerance: f64,
    /// Maximum width of the returned affine path-fraction bracket, in (0, 1].
    pub fraction_tolerance: f64,
    /// Bounded bisection work; in 1..=128. Exhaustion publishes no event.
    pub maximum_bisections: usize,
}

/// Immutable tensile criterion borrowing the original elastic material.
#[derive(Clone, Debug, PartialEq)]
pub struct TensileInitiator<'a> {
    elastic: &'a IsotropicElastic,
    strength_pa: f64,
    settings: TensileInitiationSettings,
    stiffness: Tangent6,
}

/// A numerical or physical refusal leaves the supplied elastic state untouched.
#[derive(Clone, Debug, PartialEq)]
pub enum TensileInitiationError {
    Material(BrittleError),
    Cancelled,
    StartNotElastic { principal_stress_pa: f64 },
    UnresolvedSpectrum,
    OnsetBudget,
    AmbiguousPlane { multiplicity: usize },
    DirectionOutsidePrincipalSpace,
    MaterialMismatch,
    PlaneMismatch,
    PriorMixedModeInitiation,
    StressDiscontinuity { maximum_pa: f64 },
    StorageDiscontinuity { defect_j_m3: f64 },
}
impl core::fmt::Display for TensileInitiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tensile crack initiation refused: {self:?}")
    }
}
impl std::error::Error for TensileInitiationError {}
impl From<BrittleError> for TensileInitiationError {
    fn from(error: BrittleError) -> Self {
        Self::Material(error)
    }
}

/// A disposable onset on the supplied affine increment. Fractions are numerical
/// resolution results, not interval-arithmetic certificates. The original
/// elastic source is borrowed through the result and remains inspectable.
#[derive(Clone, Debug, PartialEq)]
pub struct TensileOnset<'a> {
    elastic: &'a IsotropicElastic,
    strength_pa: f64,
    stress_tolerance_pa: f64,
    start: Voigt,
    end: Voigt,
    fraction: [f64; 2],
    normal: [f64; 3],
    principal_stresses_pa: [f64; 3],
    multiplicity: usize,
    eigen_residual_pa: f64,
    bisections: usize,
}

/// The existing elastic storage represented by an undamaged native band.
/// No work or dissipation is charged by this change of constitutive owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ElasticBrittleHandoff<'a> {
    pub state: BrittleState<'a>,
    pub elastic_storage_j_m3: f64,
    pub storage_defect_j_m3: f64,
    pub maximum_stress_defect_pa: f64,
}

impl<'a> TensileInitiator<'a> {
    /// Select a tensile criterion from explicit original coefficients. A source
    /// adapter should pass its exact admitted modulus, Poisson ratio and tensile
    /// strength; hardness or flexural modulus of rupture is not a substitute.
    pub fn new(
        elastic: &'a IsotropicElastic,
        tensile_strength_pa: f64,
        settings: TensileInitiationSettings,
    ) -> Result<Self, TensileInitiationError> {
        positive(elastic.youngs)?;
        positive(elastic.strain_limit)?;
        positive(tensile_strength_pa)?;
        positive(settings.maximum_strain)?;
        positive(settings.fraction_tolerance)?;
        if !(elastic.poisson > -1.0 && elastic.poisson < 0.5)
            || settings.maximum_strain > 0.05
            || settings.maximum_strain > elastic.strain_limit
            || !(1024.0 * f64::EPSILON..=0.01).contains(&settings.relative_stress_tolerance)
            || settings.fraction_tolerance > 1.0
            || !(1..=128).contains(&settings.maximum_bisections)
        {
            return Err(BrittleError::Invalid(
                "invalid tensile initiation domain or resolution budget",
            )
            .into());
        }
        positive(tensile_strength_pa * settings.relative_stress_tolerance)?;
        let stiffness = elastic.stiffness();
        finite_tangent(&stiffness)?;
        Ok(Self {
            elastic,
            strength_pa: tensile_strength_pa,
            settings,
            stiffness,
        })
    }

    #[must_use]
    pub const fn elastic(&self) -> &'a IsotropicElastic {
        self.elastic
    }
    #[must_use]
    pub const fn tensile_strength_pa(&self) -> f64 {
        self.strength_pa
    }
    #[must_use]
    pub const fn settings(&self) -> TensileInitiationSettings {
        self.settings
    }

    /// Locate the first tensile crossing, or return `None` for an entirely
    /// elastic affine increment. A start at/above the criterion refuses.
    /// `direction_hint` supplies a physical tie-break only when needed.
    pub fn between(
        &self,
        start: &Voigt,
        end: &Voigt,
        direction_hint: Option<[f64; 3]>,
    ) -> Result<Option<TensileOnset<'a>>, TensileInitiationError> {
        self.between_controlled(start, end, direction_hint, &mut || false)
    }

    /// Same operation with cancellation polled before spectral work and at
    /// every bisection. A scalar point performs at most 128 bisections; batch
    /// callers retain responsibility for their own work and memory admission.
    pub fn between_controlled<C: FnMut() -> bool>(
        &self,
        start: &Voigt,
        end: &Voigt,
        direction_hint: Option<[f64; 3]>,
        cancelled: &mut C,
    ) -> Result<Option<TensileOnset<'a>>, TensileInitiationError> {
        poll(cancelled)?;
        let hint = direction_hint.map(normalize).transpose()?;
        self.check_strain(start)?;
        self.check_strain(end)?;
        let a = self.stress(start)?;
        let b = self.stress(end)?;
        let tolerance = self.strength_pa * self.settings.relative_stress_tolerance;
        if !self.is_below(&a, tolerance, cancelled)? {
            let spectrum = spectrum(&a, tolerance, cancelled)?;
            return Err(TensileInitiationError::StartNotElastic {
                principal_stress_pa: spectrum.values[2],
            });
        }
        if self.is_below(&b, tolerance, cancelled)? {
            poll(cancelled)?;
            return Ok(None);
        }
        let delta = core::array::from_fn(|i| b[i] - a[i]);
        finite_voigt(&delta)?;
        // Symmetric matrix infinity norm bounds every eigenvalue's variation.
        let slope = finite(
            (delta[0].abs() + delta[3].abs() + delta[5].abs())
                .max(delta[1].abs() + delta[3].abs() + delta[4].abs())
                .max(delta[2].abs() + delta[4].abs() + delta[5].abs()),
        )?;
        let (mut low, mut high, mut count) = (0.0, 1.0, 0);
        while high - low > self.settings.fraction_tolerance || (high - low) * slope > tolerance {
            poll(cancelled)?;
            if count == self.settings.maximum_bisections {
                return Err(TensileInitiationError::OnsetBudget);
            }
            let mid = low + 0.5 * (high - low);
            if mid == low || mid == high {
                return Err(TensileInitiationError::OnsetBudget);
            }
            if self.is_below(&interpolate(&a, &b, mid), tolerance, cancelled)? {
                low = mid;
            } else {
                high = mid;
            }
            count += 1;
        }
        // Recompute from strain through the original stiffness instead of
        // carrying interpolation roundoff into the selected physical plane.
        let onset_stress = self.stress(&interpolate(start, end, high))?;
        let spectral = spectrum(&onset_stress, tolerance, cancelled)?;
        if (spectral.values[2] - self.strength_pa).abs() > 2.0 * tolerance {
            return Err(TensileInitiationError::UnresolvedSpectrum);
        }
        let (normal, multiplicity) = principal_plane(&spectral, tolerance, hint)?;
        poll(cancelled)?;
        Ok(Some(TensileOnset {
            elastic: self.elastic,
            strength_pa: self.strength_pa,
            stress_tolerance_pa: tolerance,
            start: *start,
            end: *end,
            fraction: [low, high],
            normal,
            principal_stresses_pa: spectral.values,
            multiplicity,
            eigen_residual_pa: spectral.residual_pa,
            bisections: count,
        }))
    }

    fn check_strain(&self, strain: &Voigt) -> Result<(), BrittleError> {
        finite_voigt(strain)?;
        let scale = strain.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        if scale == 0.0 {
            return Ok(());
        }
        let scaled = strain.map(|x| x / scale);
        let norm = finite(scale * det::sqrt(contract(&scaled, &scaled)))?;
        if norm > self.settings.maximum_strain {
            return Err(BrittleError::StrainLimit {
                norm,
                maximum: self.settings.maximum_strain,
            });
        }
        Ok(())
    }
    fn stress(&self, strain: &Voigt) -> Result<Voigt, BrittleError> {
        let stress =
            core::array::from_fn(|i| (0..6).map(|j| self.stiffness[i][j] * strain[j]).sum());
        finite_voigt(&stress)?;
        Ok(stress)
    }
    fn is_below<C: FnMut() -> bool>(
        &self,
        stress: &Voigt,
        tolerance: f64,
        cancelled: &mut C,
    ) -> Result<bool, TensileInitiationError> {
        match threshold_side(stress, self.strength_pa) {
            Some(below) => Ok(below),
            None => Ok(spectrum(stress, tolerance, cancelled)?.values[2] < self.strength_pa),
        }
    }
}

impl<'a> TensileOnset<'a> {
    #[must_use]
    pub const fn elastic(&self) -> &'a IsotropicElastic {
        self.elastic
    }
    /// Numerical before/after fractions in the original affine increment.
    #[must_use]
    pub const fn fraction_bracket(&self) -> [f64; 2] {
        self.fraction
    }
    #[must_use]
    pub fn strain_before_onset(&self) -> Voigt {
        interpolate(&self.start, &self.end, self.fraction[0])
    }
    #[must_use]
    pub fn strain_at_onset(&self) -> Voigt {
        interpolate(&self.start, &self.end, self.fraction[1])
    }
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }
    /// Ascending original principal elastic stresses at the located onset, Pa.
    #[must_use]
    pub const fn principal_stresses_pa(&self) -> [f64; 3] {
        self.principal_stresses_pa
    }
    #[must_use]
    pub const fn principal_multiplicity(&self) -> usize {
        self.multiplicity
    }
    #[must_use]
    pub const fn eigen_residual_pa(&self) -> f64 {
        self.eigen_residual_pa
    }
    #[must_use]
    pub const fn bisections(&self) -> usize {
        self.bisections
    }

    /// Initialize a matching native band at the increment START. The caller
    /// supplies directional geometry and all independently admitted fracture
    /// coefficients to that owner, then replays the increment through it.
    /// This method never accepts or resets an existing brittle history.
    pub fn seed<'b>(
        &self,
        band: &'b BrittleCrackBand,
    ) -> Result<ElasticBrittleHandoff<'b>, TensileInitiationError> {
        if band.elastic() != self.elastic
            || band.parameters().tensile_strength_pa != self.strength_pa
        {
            return Err(TensileInitiationError::MaterialMismatch);
        }
        let normal = band.normal();
        let sign = if (0..3).map(|i| normal[i] * self.normal[i]).sum::<f64>() < 0.0 {
            -1.0
        } else {
            1.0
        };
        if (0..3).any(|i| (normal[i] - sign * self.normal[i]).abs() > 128.0 * f64::EPSILON) {
            return Err(TensileInitiationError::PlaneMismatch);
        }
        band.check_strain(&self.start)?;
        band.check_strain(&self.end)?;
        let seed = band.trial(&self.start, &band.initial_state())?;
        if seed.state.damage() != 0.0
            || seed.state.dissipated_j_m3() != 0.0
            || seed.state.cohesive_point().mixed_mode_history().is_some()
        {
            return Err(TensileInitiationError::PriorMixedModeInitiation);
        }
        let stiffness = self.elastic.stiffness();
        let stress: Voigt =
            core::array::from_fn(|i| (0..6).map(|j| stiffness[i][j] * self.start[j]).sum());
        let maximum_stress_defect_pa = (0..6)
            .map(|i| (seed.stress_pa[i] - stress[i]).abs())
            .fold(0.0, f64::max);
        if maximum_stress_defect_pa > self.stress_tolerance_pa {
            return Err(TensileInitiationError::StressDiscontinuity {
                maximum_pa: maximum_stress_defect_pa,
            });
        }
        let elastic_storage_j_m3 = finite(0.5 * contract(&stress, &self.start))?;
        let storage_defect_j_m3 = finite(seed.state.stored_energy_j_m3() - elastic_storage_j_m3)?;
        if storage_defect_j_m3.abs() > band.settings().energy_tolerance_j_m3 {
            return Err(TensileInitiationError::StorageDiscontinuity {
                defect_j_m3: storage_defect_j_m3,
            });
        }
        Ok(ElasticBrittleHandoff {
            state: seed.state,
            elastic_storage_j_m3,
            storage_defect_j_m3,
            maximum_stress_defect_pa,
        })
    }
}

fn interpolate(a: &Voigt, b: &Voigt, fraction: f64) -> Voigt {
    core::array::from_fn(|i| (1.0 - fraction) * a[i] + fraction * b[i])
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), TensileInitiationError> {
    if cancelled() {
        Err(TensileInitiationError::Cancelled)
    } else {
        Ok(())
    }
}

// Sylvester's criterion on strength*I-stress, with scaled arithmetic and
// roundoff guards. Near-zero minors defer to the ORIGINAL spectral owner.
// These guards are numerical screens, not interval-certificate authority.
fn threshold_side(stress: &Voigt, strength: f64) -> Option<bool> {
    let scale = stress.iter().fold(strength, |a, b| a.max(b.abs()));
    let s = stress.map(|x| x / scale);
    let t = strength / scale;
    let guard = 32.0 * f64::EPSILON;
    let upper = (s[0] + s[3].abs() + s[5].abs())
        .max(s[1] + s[3].abs() + s[4].abs())
        .max(s[2] + s[4].abs() + s[5].abs());
    if upper < t - guard {
        return Some(true);
    }
    if s[0].max(s[1]).max(s[2]) > t + guard {
        return Some(false);
    }
    let [a, b, c, d, e, f] = [t - s[0], t - s[1], t - s[2], -s[3], -s[4], -s[5]];
    let second = a * b - d * d;
    let third = a * (b * c - e * e) - d * (d * c - e * f) + f * (d * e - b * f);
    let g1 = guard * (t.abs() + s[0].abs());
    let g2 = 2.0 * guard * ((a * b).abs() + d * d);
    let g3 = 4.0
        * guard
        * ((a * b * c).abs()
            + 2.0 * (d * e * f).abs()
            + (a * e * e).abs()
            + (b * f * f).abs()
            + (c * d * d).abs());
    if a > g1 && second > g2 && third > g3 {
        Some(true)
    } else if a < -g1 || second < -g2 || third < -g3 {
        Some(false)
    } else {
        None
    }
}

struct Spectrum {
    values: [f64; 3],
    vectors: [[f64; 3]; 3],
    residual_pa: f64,
}

// One objective plane-selection rule for a located affine crossing or an
// explicitly observed loaded endpoint. The optional hint is already normalized.
fn principal_plane(
    spectral: &Spectrum,
    tolerance: f64,
    hint: Option<[f64; 3]>,
) -> Result<([f64; 3], usize), TensileInitiationError> {
    let mut first = 2;
    while first > 0 && spectral.values[2] - spectral.values[first - 1] <= tolerance {
        first -= 1;
    }
    let multiplicity = 3 - first;
    let mut normal = if multiplicity == 1 {
        core::array::from_fn(|i| spectral.vectors[i][2])
    } else {
        let hint = hint.ok_or(TensileInitiationError::AmbiguousPlane { multiplicity })?;
        let mut projected = [0.0; 3];
        for j in first..3 {
            let coefficient = (0..3)
                .map(|i| spectral.vectors[i][j] * hint[i])
                .sum::<f64>();
            for i in 0..3 {
                projected[i] += coefficient * spectral.vectors[i][j];
            }
        }
        if projected.iter().fold(0.0_f64, |a, b| a.max(b.abs())) <= 128.0 * f64::EPSILON {
            return Err(TensileInitiationError::DirectionOutsidePrincipalSpace);
        }
        normalize(projected)?
    };
    normal = normalize(normal)?;
    // Fix only the sign, which does not change the physical plane. A hint
    // can carry orientation; otherwise use a deterministic largest entry.
    let sign = match hint {
        Some(hint) => (0..3).map(|i| normal[i] * hint[i]).sum::<f64>(),
        None => {
            let mut index = 0;
            for i in 1..3 {
                if normal[i].abs() > normal[index].abs() {
                    index = i;
                }
            }
            normal[index]
        }
    };
    if sign < 0.0 {
        normal = normal.map(|x| -x);
    }
    Ok((normal, multiplicity))
}

fn spectrum<C: FnMut() -> bool>(
    stress: &Voigt,
    tolerance: f64,
    cancelled: &mut C,
) -> Result<Spectrum, TensileInitiationError> {
    poll(cancelled)?;
    let scale = stress.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
    if scale == 0.0 {
        return Ok(Spectrum {
            values: [0.0; 3],
            vectors: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            residual_pa: 0.0,
        });
    }
    let s = stress.map(|x| x / scale);
    let matrix = [s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2]];
    let (values, vectors) = fs_la::eigen::jacobi_eigh(&matrix, 3);
    let v: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| vectors[3 * i + j]));
    let mut residual = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let av = (0..3).map(|k| matrix[3 * i + k] * v[k][j]).sum::<f64>();
            residual = residual.max((av - values[j] * v[i][j]).abs());
            let gram = (0..3).map(|k| v[k][i] * v[k][j]).sum::<f64>();
            if !gram.is_finite()
                || (gram - if i == j { 1.0 } else { 0.0 }).abs() > 256.0 * f64::EPSILON
            {
                return Err(TensileInitiationError::UnresolvedSpectrum);
            }
        }
    }
    let residual_pa = finite(scale * residual)?;
    if residual_pa > tolerance {
        return Err(TensileInitiationError::UnresolvedSpectrum);
    }
    let values = [
        finite(values[0] * scale)?,
        finite(values[1] * scale)?,
        finite(values[2] * scale)?,
    ];
    poll(cancelled)?;
    Ok(Spectrum {
        values,
        vectors: v,
        residual_pa,
    })
}
