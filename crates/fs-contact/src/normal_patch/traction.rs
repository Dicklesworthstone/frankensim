//! Spatial adhesive tractions from the original axisymmetric normal sources.
//!
//! JKR is the Hertz pressure minus a tensile flat-punch field; its integrable
//! edge singularity must not be sampled as a finite nodal pressure. Maugis
//! retains both its finite contact pressure and constant-traction cohesive
//! annulus (Leroy & Willert 2021, DOI 10.1002/zamm.201900313, equation 7).
//! Annular forces are exact integrals of these selected continuum estimates,
//! evaluated without numerical quadrature or subtraction of cumulative loads.
//!
//! A transfer preserves signed force and work for piecewise constant normal
//! motion on the supplied concentric rings. It does not solve nonaxisymmetric
//! deformation, add finite-body compliance, or supply a new contact tangent.

use super::{
    NormalLawSource,
    adhesive::{JkrError, SourceJkrSphereContact},
    maugis::{
        SourceMaugisSphereContact,
        continuation::{MaugisContinuationError, MaugisEquilibrium},
    },
};
use core::f64::consts::PI;
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_tribo::InputAuthority;

/// Fixed bound on one atomic transfer's work and caller-owned ring storage.
pub const MAX_ANNULAR_TRACTION_RINGS: usize = 4096;

/// Spatial-query or original-source refusal. No caller output changes on failure.
#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveTractionError {
    /// The original JKR source rejected the requested attached state.
    Jkr(JkrError),
    /// The original finite-zone source rejected the requested branch or state.
    Maugis(MaugisContinuationError),
    /// An equilibrium belongs to another original material/geometry owner.
    SourceMismatch,
    /// Invalid radii, partition, tolerance, or unrepresentable spatial data.
    Invalid(&'static str),
    /// Ideal JKR pressure is unbounded exactly at the attached contact edge.
    SingularEdge,
    /// The ring forces do not recover the original force within the given gate.
    ForceBalance {
        /// Sum of ring forces minus the original resultant, N.
        residual_n: f64,
        /// Caller-authored absolute arithmetic tolerance, N.
        tolerance_n: f64,
    },
}
impl core::fmt::Display for AdhesiveTractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive spatial traction: {self:?}")
    }
}
impl std::error::Error for AdhesiveTractionError {}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Distribution {
    Jkr {
        hertz_peak_pa: f64,
        tensile_pa: f64,
    },
    Maugis {
        hertz_peak_pa: f64,
        stress_pa: f64,
        annulus_root_m: f64,
    },
}

/// Immutable spatial profile, constructed from original source evaluations.
/// Setup retains source identity and authority; queries allocate no memory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTractionProfile {
    model_identity: ContentHash,
    normal_state_identity: ContentHash,
    authority: InputAuthority,
    contact_radius_m: f64,
    interaction_radius_m: f64,
    normal_force_n: f64,
    distribution: Distribution,
}

/// Integrated signed load on one ring, including any unsupported outer area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularAdhesiveTraction {
    /// Entire caller-specified ring area, m².
    pub area_m2: f64,
    /// Portion with real solid contact; a cohesive annulus is excluded, m².
    pub solid_contact_area_m2: f64,
    /// Integrated Hertz term, N; not the positive part of the total pressure.
    pub hertz_force_n: f64,
    /// Integrated tensile term, N; nonpositive.
    pub adhesive_force_n: f64,
    /// Signed sum of Hertz and adhesive contributions, N.
    pub normal_force_n: f64,
    /// Area mean over the entire ring, Pa; finite even at the JKR edge.
    pub mean_pressure_pa: f64,
}

impl NormalLawSource {
    /// Spatial traction for the existing JKR/Maugis mechanical source, retaining
    /// its branch admission (including opt-in Maugis noncontact continuation).
    /// Other bindings return `None`: a scalar force curve or DMT-M resultant
    /// does not uniquely determine a spatial traction field. No new footprint
    /// approximation is selected implicitly for those owners.
    pub fn adhesive_traction_profile(
        &self,
        indentation_m: f64,
    ) -> Result<Option<AdhesiveTractionProfile>, AdhesiveTractionError> {
        match self {
            Self::JkrAttached(source) => {
                AdhesiveTractionProfile::jkr_attached(source, indentation_m).map(Some)
            }
            Self::MaugisAttached(source) => {
                let equilibrium = source
                    .evaluate(indentation_m)
                    .map_err(AdhesiveTractionError::Maugis)?;
                AdhesiveTractionProfile::maugis(source.source(), &equilibrium).map(Some)
            }
            Self::HertzPoint(_) | Self::FiniteGap(_) | Self::DmtSphere(_) => Ok(None),
        }
    }
}

impl AdhesiveTractionProfile {
    /// Explicitly select the original attached JKR branch, including tensile
    /// net load and its displacement fold. No branch history is advanced.
    pub fn jkr_attached(
        source: &SourceJkrSphereContact,
        indentation_m: f64,
    ) -> Result<Self, AdhesiveTractionError> {
        let response = source
            .evaluate_attached(indentation_m)
            .map_err(AdhesiveTractionError::Jkr)?;
        let a = response.contact_radius_m;
        let e = source.normal_state().elastic().reduced_modulus_pa();
        let w = source.normal_state().adhesion_energy_j_per_m2();
        // p = p_H sqrt(1-r²/a²) - q/sqrt(1-r²/a²).
        // Its integrated attraction is -sqrt(8 pi E* W a³).
        let hertz_peak_pa = (2.0 / PI) * e * (a / source.effective_radius_m());
        let tensile_pa = ((2.0 / PI) * (e / a) * w).sqrt();
        finite(hertz_peak_pa)?;
        finite(tensile_pa)?;
        Ok(Self {
            model_identity: source.identity(),
            normal_state_identity: source.normal_state().identity(),
            authority: response.authority,
            contact_radius_m: a,
            interaction_radius_m: a,
            normal_force_n: response.normal_force_n,
            distribution: Distribution::Jkr {
                hertz_peak_pa,
                tensile_pa,
            },
        })
    }

    /// Retain an original immutable Maugis equilibrium. Noncontact attraction
    /// loads a finite disk but has zero solid area; detached states load nothing.
    /// At a contact fold pressure remains finite even though a tangent is absent.
    pub fn maugis(
        source: &SourceMaugisSphereContact,
        equilibrium: &MaugisEquilibrium,
    ) -> Result<Self, AdhesiveTractionError> {
        if source.identity() != equilibrium.model_identity() {
            return Err(AdhesiveTractionError::SourceMismatch);
        }
        let (a, l) = equilibrium
            .contact_response()
            .map_or((0.0, 0.0), |r| (r.contact_radius_m, r.annulus_root_m));
        let e = source.normal_state().elastic().reduced_modulus_pa();
        let hertz_peak_pa = (2.0 / PI) * e * (a / source.effective_radius_m());
        let stress_pa =
            source.normal_state().adhesion_energy_j_per_m2() / source.interaction_range_m();
        finite(hertz_peak_pa)?;
        finite(stress_pa)?;
        Ok(Self {
            model_identity: source.identity(),
            normal_state_identity: source.normal_state().identity(),
            authority: equilibrium.authority(),
            contact_radius_m: a,
            interaction_radius_m: equilibrium.cohesive_radius_m(),
            normal_force_n: equilibrium.normal_force_n(),
            distribution: Distribution::Maugis {
                hertz_peak_pa,
                stress_pa,
                annulus_root_m: l,
            },
        })
    }

    /// Original model identity, without changing material or geometry authority.
    #[must_use]
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }
    /// Original ordered-interface and bulk state identity.
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_state_identity
    }
    /// Original estimate/synthetic authority, never upgraded by integration.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }
    /// Radius of real solid contact, m; excludes the cohesive annulus.
    #[must_use]
    pub const fn contact_radius_m(&self) -> f64 {
        self.contact_radius_m
    }
    /// Outer radius carrying normal traction, m.
    #[must_use]
    pub const fn interaction_radius_m(&self) -> f64 {
        self.interaction_radius_m
    }
    /// Original signed source resultant, N.
    #[must_use]
    pub const fn normal_force_n(&self) -> f64 {
        self.normal_force_n
    }

    /// Pointwise continuum pressure, positive in compression. At an ideal JKR
    /// edge use `annulus` instead: clipping its singular tension would lose force.
    pub fn pressure_pa(&self, radius_m: f64) -> Result<f64, AdhesiveTractionError> {
        radius(radius_m)?;
        let a = self.contact_radius_m;
        if radius_m > self.interaction_radius_m || self.interaction_radius_m == 0.0 {
            return Ok(0.0);
        }
        let root = if radius_m < a {
            edge_root(a, radius_m)
        } else {
            0.0
        };
        finite(match self.distribution {
            Distribution::Jkr {
                hertz_peak_pa,
                tensile_pa,
            } => {
                if radius_m == a {
                    return Err(AdhesiveTractionError::SingularEdge);
                }
                hertz_peak_pa * root - tensile_pa / root
            }
            Distribution::Maugis {
                hertz_peak_pa,
                stress_pa,
                annulus_root_m,
            } => {
                if radius_m >= a {
                    -stress_pa
                } else {
                    hertz_peak_pa * root - (2.0 / PI) * stress_pa * annulus_root_m.atan2(a * root)
                }
            }
        })
    }

    /// Exact radial force integral over a concentric ring. Factored differences
    /// preserve narrow rings and finite JKR edge loads, without quadrature.
    /// Outside the interaction footprint the pressure is zero.
    pub fn annulus(
        &self,
        inner_m: f64,
        outer_m: f64,
    ) -> Result<AnnularAdhesiveTraction, AdhesiveTractionError> {
        radius(inner_m)?;
        radius(outer_m)?;
        if outer_m <= inner_m {
            return Err(AdhesiveTractionError::Invalid("ordered nonempty annulus"));
        }
        let area_m2 = ring_area(inner_m, outer_m)?;
        let a = self.contact_radius_m;
        let lo = inner_m.min(a);
        let hi = outer_m.min(a);
        let mut solid_contact_area_m2 = 0.0;
        let mut hertz_force_n = 0.0;
        let mut adhesive_force_n = 0.0;
        if hi > lo {
            solid_contact_area_m2 = ring_area(lo, hi)?;
            let x = edge_root(a, lo);
            let y = edge_root(a, hi);
            let sum = x + y;
            let hertz_mean = (2.0 / 3.0) * (x * x + x * y + y * y) / sum;
            let (peak, attraction_mean) = match self.distribution {
                Distribution::Jkr {
                    hertz_peak_pa,
                    tensile_pa,
                } => (hertz_peak_pa, -2.0 * (tensile_pa / sum)),
                Distribution::Maugis {
                    hertz_peak_pa,
                    stress_pa,
                    annulus_root_m,
                } => {
                    let angle = mean_cohesive_angle(a, annulus_root_m, lo, hi, x, y)?;
                    (hertz_peak_pa, -(2.0 / PI) * stress_pa * angle)
                }
            };
            hertz_force_n = finite(solid_contact_area_m2 * peak * hertz_mean)?;
            adhesive_force_n = finite(solid_contact_area_m2 * attraction_mean)?;
        }
        if let Distribution::Maugis { stress_pa, .. } = self.distribution {
            let left = inner_m.max(a).min(self.interaction_radius_m);
            let right = outer_m.min(self.interaction_radius_m);
            if right > left {
                adhesive_force_n = finite(adhesive_force_n - stress_pa * ring_area(left, right)?)?;
            }
        }
        let normal_force_n = finite(hertz_force_n + adhesive_force_n)?;
        Ok(AnnularAdhesiveTraction {
            area_m2,
            solid_contact_area_m2,
            hertz_force_n,
            adhesive_force_n,
            normal_force_n,
            mean_pressure_pa: finite(normal_force_n / area_m2)?,
        })
    }

    /// Fill caller-owned ring forces after checking the COMPLETE partition and
    /// original resultant. Radii start at zero and cover the full interaction
    /// footprint. Two linear passes use constant scratch memory; any refusal
    /// leaves `forces_n` untouched. Returns sum(forces)-original force, N.
    /// These are work-conjugate efforts for piecewise constant ring velocities.
    /// One atomic transfer is bounded by `MAX_ANNULAR_TRACTION_RINGS`.
    pub fn transfer_annular_forces(
        &self,
        radii_m: &[f64],
        forces_n: &mut [f64],
        absolute_tolerance_n: f64,
    ) -> Result<f64, AdhesiveTractionError> {
        if radii_m.len() < 2
            || forces_n.len() > MAX_ANNULAR_TRACTION_RINGS
            || forces_n.len() != radii_m.len() - 1
            || radii_m[0] != 0.0
            || !absolute_tolerance_n.is_finite()
            || absolute_tolerance_n < 0.0
        {
            return Err(AdhesiveTractionError::Invalid(
                "complete annular partition and force gate",
            ));
        }
        if *radii_m.last().unwrap() < self.interaction_radius_m {
            return Err(AdhesiveTractionError::Invalid(
                "partition does not cover adhesive footprint",
            ));
        }
        let mut sum = Dd::from_f64(0.0);
        for pair in radii_m.windows(2) {
            sum = sum + Dd::from_f64(self.annulus(pair[0], pair[1])?.normal_force_n);
        }
        let residual_n = finite((sum - Dd::from_f64(self.normal_force_n)).to_f64())?;
        if residual_n.abs() > absolute_tolerance_n {
            return Err(AdhesiveTractionError::ForceBalance {
                residual_n,
                tolerance_n: absolute_tolerance_n,
            });
        }
        for (pair, force) in radii_m.windows(2).zip(forces_n.iter_mut()) {
            *force = self.annulus(pair[0], pair[1])?.normal_force_n;
        }
        Ok(residual_n)
    }
}

// Dimensionless sqrt(a²-r²)/a, without forming a² or subtracting two squares.
fn edge_root(a: f64, r: f64) -> f64 {
    (((a - r) / a) * (1.0 + r / a)).sqrt()
}
fn radius(r: f64) -> Result<(), AdhesiveTractionError> {
    if !r.is_finite() || r < 0.0 {
        Err(AdhesiveTractionError::Invalid("finite nonnegative radius"))
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64, AdhesiveTractionError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(AdhesiveTractionError::Invalid(
            "unrepresentable spatial traction",
        ))
    }
}
fn ring_area(lo: f64, hi: f64) -> Result<f64, AdhesiveTractionError> {
    let area = finite(PI * (hi - lo) * (hi + lo))?;
    if area <= 0.0 {
        return Err(AdhesiveTractionError::Invalid(
            "unrepresentable nonzero ring area",
        ));
    }
    Ok(area)
}

// Average atan(l/sqrt(a²-r²)) with radial area measure. Integration by parts
// gives an exact primitive, but subtracting it loses thin-ring and thin-zone
// forces. This rearrangement has only nonnegative terms. Scaling by max(a,l)
// keeps the angle-difference denominator away from avoidable SI overflow.
fn mean_cohesive_angle(
    a: f64,
    l: f64,
    lo: f64,
    hi: f64,
    sx: f64,
    sy: f64,
) -> Result<f64, AdhesiveTractionError> {
    let scale = a.max(l);
    let x = sx * (a / scale);
    let y = sy * (a / scale);
    let l = l / scale;
    let sum = x + y;
    let dx = (((hi - lo) / scale) * (hi / scale + lo / scale)) / sum;
    let denominator = l * l + x * y;
    let numerator = l * dx;
    let (q, defect) = if numerator < 0.125 * denominator {
        let z = numerator / denominator;
        let z2 = z * z;
        let defect = z2
            * (1.0 / 3.0
                + z2 * (-1.0 / 5.0
                    + z2 * (1.0 / 7.0
                        + z2 * (-1.0 / 9.0
                            + z2 * (1.0 / 11.0
                                + z2 * (-1.0 / 13.0
                                    + z2 * (1.0 / 15.0 + z2 * (-1.0 / 17.0 + z2 / 19.0))))))));
        (1.0 - defect, defect)
    } else {
        let q = numerator.atan2(denominator) * (denominator / numerator);
        (q, 1.0 - q)
    };
    finite(l.atan2(x) + (l / sum) * (defect + ((y * dx) / denominator) * q))
}
