//! Two-dimensional transfer of the original axisymmetric adhesive field.
//!
//! A sector carries signed force AND first moments. Its resultant can vanish
//! while its couple remains nonzero, so no fictitious application point is
//! obtained by dividing moment by force. The local frame is right-handed,
//! with z along the applied normal traction. These numerical loads retain the
//! originating profile; they do not solve finite-body compliance or advance a
//! mechanical/thermal state. Counterpart bodies receive opposite wrenches.

mod centroid;

use super::{
    AdhesiveTractionError, AdhesiveTractionProfile, AnnularAdhesiveTraction, Distribution,
    edge_root, finite, mean_cohesive_angle,
};
use core::f64::consts::TAU;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::InputAuthority;

/// Bound on one atomic polar-grid publication.
pub const MAX_POLAR_TRACTION_SECTORS: usize = 4096;
/// Bound on all Maugis radial evaluations in one query or transfer.
pub const MAX_RADIAL_MOMENT_EVALUATIONS: usize = 65_536;

/// Accuracy for each radial first moment integral, ∫p(r)r²dr, in N m.
/// The relative scale is its nonnegative Maugis attractive contribution.
/// Exact Hertz/JKR and constant-annulus formulas consume no quadrature budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialMomentBudget {
    pub absolute_tolerance_n_m: f64,
    pub relative_tolerance: f64,
    /// Shared across every ring of a bulk transfer, not reset per sector.
    pub maximum_evaluations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarTransferBudget {
    pub radial_moment: RadialMomentBudget,
    pub force_absolute_tolerance_n: f64,
    pub moment_absolute_tolerance_n_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SectorTractionError {
    Source(AdhesiveTractionError),
    Invalid(&'static str),
    Cancelled,
    MomentBudget {
        evaluations: usize,
        maximum: usize,
    },
    /// Accuracy, geometric support or analytic force normalization failed.
    MomentUnresolved,
    Resultant {
        force_residual_n: f64,
        first_moment_residual_n_m: [f64; 2],
    },
}
impl From<AdhesiveTractionError> for SectorTractionError {
    fn from(value: AdhesiveTractionError) -> Self {
        Self::Source(value)
    }
}
impl core::fmt::Display for SectorTractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive polar traction: {self:?}")
    }
}
impl std::error::Error for SectorTractionError {}

/// Numerical sector load, to be retained with its originating source profile.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SectorAdhesiveTraction {
    pub area_m2: f64,
    pub solid_contact_area_m2: f64,
    pub hertz_force_n: f64,
    pub adhesive_force_n: f64,
    pub normal_force_n: f64,
    pub mean_pressure_pa: f64,
    /// [∫x p dA, ∫y p dA] about the original contact centre, N m.
    pub first_moment_n_m: [f64; 2],
    /// Numerical quadrature estimate for the Euclidean first-moment error.
    /// Excludes floating-point rounding, material and continuum-model errors;
    /// it is not a certified enclosure. Zero for the closed-form JKR path.
    pub moment_error_estimate_n_m: f64,
}

impl SectorAdhesiveTraction {
    /// Local normal force and torque about the original contact centre.
    #[must_use]
    pub fn wrench(&self) -> ([f64; 3], [f64; 3]) {
        (
            [0.0, 0.0, self.normal_force_n],
            [self.first_moment_n_m[1], -self.first_moment_n_m[0], 0.0],
        )
    }

    /// Work-conjugate power against v_z(x,y)=translation+gradient·[x,y].
    /// For rigid normal motion gradient=[-omega_y,omega_x]. Using relative
    /// velocity accounts for both opposing body efforts exactly once.
    /// Returns (power, moment-quadrature power-error estimate), W.
    pub fn affine_normal_power(
        &self,
        translation_m_per_s: f64,
        gradient_per_s: [f64; 2],
    ) -> Result<(f64, f64), SectorTractionError> {
        for value in [translation_m_per_s, gradient_per_s[0], gradient_per_s[1]] {
            finite(value)?;
        }
        let power = (Dd::from_f64(self.normal_force_n) * Dd::from_f64(translation_m_per_s)
            + Dd::from_f64(self.first_moment_n_m[0]) * Dd::from_f64(gradient_per_s[0])
            + Dd::from_f64(self.first_moment_n_m[1]) * Dd::from_f64(gradient_per_s[1]))
        .to_f64();
        Ok((
            finite(power)?,
            finite(self.moment_error_estimate_n_m * gradient_per_s[0].hypot(gradient_per_s[1]))?,
        ))
    }
}

/// Reusable caller-owned radial storage. Scratch may change on refusal;
/// accepted output never does. Fields are private to prevent stale reuse.
#[derive(Clone, Copy, Debug)]
pub struct PolarRingScratch {
    ring: AnnularAdhesiveTraction,
    first_moment_n_m: f64,
    error_estimate_n_m: f64,
}
impl Default for PolarRingScratch {
    fn default() -> Self {
        Self {
            ring: AnnularAdhesiveTraction {
                area_m2: 0.0,
                solid_contact_area_m2: 0.0,
                hertz_force_n: 0.0,
                adhesive_force_n: 0.0,
                normal_force_n: 0.0,
                mean_pressure_pa: 0.0,
            },
            first_moment_n_m: 0.0,
            error_estimate_n_m: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarTractionTransfer {
    pub model_identity: ContentHash,
    pub normal_state_identity: ContentHash,
    pub authority: InputAuthority,
    pub force_residual_n: f64,
    pub first_moment_residual_n_m: [f64; 2],
    pub maximum_sector_moment_error_estimate_n_m: f64,
    pub radial_evaluations: usize,
}

impl AdhesiveTractionProfile {
    /// Query a sector bounded by radii and increasing angles in [0,2π].
    /// Angular integration and force are analytic. Only the contacting Maugis
    /// radial first moment needs bounded numerical quadrature.
    #[allow(clippy::too_many_arguments)]
    pub fn sector(
        &self,
        inner_m: f64,
        outer_m: f64,
        start_rad: f64,
        end_rad: f64,
        budget: RadialMomentBudget,
        gate: &CancelGate,
    ) -> Result<SectorAdhesiveTraction, SectorTractionError> {
        check_budget(budget)?;
        angles(start_rad, end_rad)?;
        poll(gate)?;
        let mut evaluations = 0;
        let radial = if end_rad - start_rad == TAU {
            PolarRingScratch {
                ring: self.annulus(inner_m, outer_m)?,
                ..PolarRingScratch::default()
            }
        } else {
            self.radial_load(inner_m, outer_m, budget, &mut evaluations, gate)?
        };
        let result = sector_load(radial, start_rad, end_rad);
        check_sector(result)?;
        poll(gate)?;
        Ok(result)
    }

    /// Transfer a COMPLETE polar grid, radial-major then angular-minor.
    /// Rings cover the original interaction radius and angular edges cover
    /// exactly [0,2π]. Source force and zero axisymmetric first moment are
    /// checked before publication. One radial solve is reused across angles;
    /// the final bounded output copy has no fallible work or allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn transfer_polar_grid(
        &self,
        radii_m: &[f64],
        angles_rad: &[f64],
        output: &mut [SectorAdhesiveTraction],
        scratch: &mut [PolarRingScratch],
        budget: PolarTransferBudget,
        gate: &CancelGate,
    ) -> Result<PolarTractionTransfer, SectorTractionError> {
        check_budget(budget.radial_moment)?;
        poll(gate)?;
        if radii_m.len() < 2
            || angles_rad.len() < 2
            || radii_m[0] != 0.0
            || angles_rad[0] != 0.0
            || *angles_rad.last().unwrap() != TAU
            || *radii_m.last().unwrap() < self.interaction_radius_m()
            || scratch.len() != radii_m.len() - 1
            || (radii_m.len() - 1).checked_mul(angles_rad.len() - 1) != Some(output.len())
            || output.len() > MAX_POLAR_TRACTION_SECTORS
            || !budget.force_absolute_tolerance_n.is_finite()
            || budget.force_absolute_tolerance_n < 0.0
            || !budget.moment_absolute_tolerance_n_m.is_finite()
            || budget.moment_absolute_tolerance_n_m < 0.0
        {
            return Err(SectorTractionError::Invalid(
                "complete polar grid, storage and resultant gates",
            ));
        }
        for edge in angles_rad.windows(2) {
            angles(edge[0], edge[1])?;
        }
        let mut evaluations = 0;
        for (edge, staged) in radii_m.windows(2).zip(scratch.iter_mut()) {
            poll(gate)?;
            *staged = if angles_rad.len() == 2 {
                PolarRingScratch {
                    ring: self.annulus(edge[0], edge[1])?,
                    ..PolarRingScratch::default()
                }
            } else {
                self.radial_load(
                    edge[0],
                    edge[1],
                    budget.radial_moment,
                    &mut evaluations,
                    gate,
                )?
            };
        }
        let mut sums = [Dd::from_f64(0.0); 3];
        let mut largest_error = 0.0_f64;
        for radial in scratch.iter().copied() {
            poll(gate)?;
            for edge in angles_rad.windows(2) {
                let cell = sector_load(radial, edge[0], edge[1]);
                check_sector(cell)?;
                sums[0] = sums[0] + Dd::from_f64(cell.normal_force_n);
                sums[1] = sums[1] + Dd::from_f64(cell.first_moment_n_m[0]);
                sums[2] = sums[2] + Dd::from_f64(cell.first_moment_n_m[1]);
                largest_error = largest_error.max(cell.moment_error_estimate_n_m);
            }
        }
        let residual = finite((sums[0] - Dd::from_f64(self.normal_force_n())).to_f64())?;
        let moment = [finite(sums[1].to_f64())?, finite(sums[2].to_f64())?];
        if residual.abs() > budget.force_absolute_tolerance_n
            || moment[0].hypot(moment[1]) > budget.moment_absolute_tolerance_n_m
        {
            return Err(SectorTractionError::Resultant {
                force_residual_n: residual,
                first_moment_residual_n_m: moment,
            });
        }
        poll(gate)?;
        let sectors = angles_rad.len() - 1;
        for (radial, cells) in scratch.iter().copied().zip(output.chunks_mut(sectors)) {
            for (edge, cell) in angles_rad.windows(2).zip(cells) {
                *cell = sector_load(radial, edge[0], edge[1]);
            }
        }
        Ok(PolarTractionTransfer {
            model_identity: self.model_identity(),
            normal_state_identity: self.normal_state_identity(),
            authority: self.authority(),
            force_residual_n: residual,
            first_moment_residual_n_m: moment,
            maximum_sector_moment_error_estimate_n_m: largest_error,
            radial_evaluations: evaluations,
        })
    }

    fn radial_load(
        &self,
        inner: f64,
        outer: f64,
        budget: RadialMomentBudget,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<PolarRingScratch, SectorTractionError> {
        poll(gate)?;
        let ring = self.annulus(inner, outer)?;
        let a = self.contact_radius_m();
        let lo = inner.min(a);
        let hi = outer.min(a);
        let mut hertz = 0.0;
        let mut adhesive = Dd::from_f64(0.0);
        let mut error = 0.0;
        if hi > lo {
            let means = centroid::jkr(a, lo, hi).ok_or(SectorTractionError::MomentUnresolved)?;
            hertz = finite(ring.hertz_force_n * (means[0] / TAU))?;
            match self.distribution {
                Distribution::Jkr { .. } => {
                    adhesive = Dd::from_f64(finite(ring.adhesive_force_n * (means[1] / TAU))?);
                }
                Distribution::Maugis { annulus_root_m, .. } => {
                    let contact = self.annulus(lo, hi)?;
                    if contact.adhesive_force_n != 0.0 {
                        let scale = finite(contact.adhesive_force_n.abs() * (hi / TAU))?;
                        if scale == 0.0 {
                            return Err(SectorTractionError::MomentUnresolved);
                        }
                        let angle = mean_cohesive_angle(
                            a,
                            annulus_root_m,
                            lo,
                            hi,
                            edge_root(a, lo),
                            edge_root(a, hi),
                        )?;
                        let controls = centroid::Controls {
                            absolute_fraction_tolerance: (budget.absolute_tolerance_n_m / scale)
                                .min(1.0),
                            relative_tolerance: budget.relative_tolerance,
                            maximum_evaluations: budget
                                .maximum_evaluations
                                .saturating_sub(*evaluations),
                        };
                        let mean = centroid::maugis(
                            a,
                            annulus_root_m,
                            lo,
                            hi,
                            angle,
                            controls,
                            &mut || gate.is_requested(),
                        )
                        .map_err(|failure| match failure {
                            centroid::CentroidError::Cancelled => SectorTractionError::Cancelled,
                            centroid::CentroidError::Budget { evaluations: used } => {
                                SectorTractionError::MomentBudget {
                                    evaluations: *evaluations + used,
                                    maximum: budget.maximum_evaluations,
                                }
                            }
                            centroid::CentroidError::Invalid
                            | centroid::CentroidError::Unresolved => {
                                SectorTractionError::MomentUnresolved
                            }
                        })?;
                        *evaluations += mean.evaluations;
                        adhesive =
                            Dd::from_f64(finite(contact.adhesive_force_n * (mean.radius_m / TAU))?);
                        error =
                            finite(contact.adhesive_force_n.abs() * (mean.error_estimate_m / TAU))?;
                    }
                }
            }
        }
        if let Distribution::Maugis { .. } = self.distribution {
            let left = inner.max(a).min(self.interaction_radius_m());
            let right = outer.min(self.interaction_radius_m());
            if right > left {
                let annulus = self.annulus(left, right)?;
                adhesive = adhesive
                    + Dd::from_f64(finite(
                        annulus.adhesive_force_n * (centroid::uniform(left, right) / TAU),
                    )?);
            }
        }
        Ok(PolarRingScratch {
            ring,
            first_moment_n_m: finite((Dd::from_f64(hertz) + adhesive).to_f64())?,
            error_estimate_n_m: error,
        })
    }
}

fn sector_load(radial: PolarRingScratch, start: f64, end: f64) -> SectorAdhesiveTraction {
    let width = end - start;
    let fraction = width / TAU;
    let angular = if width == TAU {
        [0.0, 0.0]
    } else {
        let amplitude = 2.0 * (0.5 * width).sin();
        let mid = 0.5 * start + 0.5 * end;
        [amplitude * mid.cos(), amplitude * mid.sin()]
    };
    SectorAdhesiveTraction {
        area_m2: radial.ring.area_m2 * fraction,
        solid_contact_area_m2: radial.ring.solid_contact_area_m2 * fraction,
        hertz_force_n: radial.ring.hertz_force_n * fraction,
        adhesive_force_n: radial.ring.adhesive_force_n * fraction,
        normal_force_n: radial.ring.normal_force_n * fraction,
        mean_pressure_pa: radial.ring.mean_pressure_pa,
        first_moment_n_m: [
            radial.first_moment_n_m * angular[0],
            radial.first_moment_n_m * angular[1],
        ],
        moment_error_estimate_n_m: radial.error_estimate_n_m * angular[0].hypot(angular[1]),
    }
}
fn check_sector(value: SectorAdhesiveTraction) -> Result<(), SectorTractionError> {
    for x in [
        value.area_m2,
        value.solid_contact_area_m2,
        value.hertz_force_n,
        value.adhesive_force_n,
        value.normal_force_n,
        value.mean_pressure_pa,
        value.first_moment_n_m[0],
        value.first_moment_n_m[1],
        value.moment_error_estimate_n_m,
    ] {
        finite(x)?;
    }
    if value.area_m2 <= 0.0 {
        return Err(SectorTractionError::Invalid("unrepresentable sector area"));
    }
    Ok(())
}
fn angles(start: f64, end: f64) -> Result<(), SectorTractionError> {
    if !start.is_finite() || !end.is_finite() || start < 0.0 || end > TAU || end <= start {
        Err(SectorTractionError::Invalid(
            "increasing sector angles in [0,2pi]",
        ))
    } else {
        Ok(())
    }
}
fn check_budget(budget: RadialMomentBudget) -> Result<(), SectorTractionError> {
    if !budget.absolute_tolerance_n_m.is_finite()
        || budget.absolute_tolerance_n_m < 0.0
        || !budget.relative_tolerance.is_finite()
        || budget.relative_tolerance < 0.0
        || budget.maximum_evaluations > MAX_RADIAL_MOMENT_EVALUATIONS
    {
        Err(SectorTractionError::Invalid(
            "finite radial moment accuracy and bounded work",
        ))
    } else {
        Ok(())
    }
}
fn poll(gate: &CancelGate) -> Result<(), SectorTractionError> {
    if gate.is_requested() {
        Err(SectorTractionError::Cancelled)
    } else {
        Ok(())
    }
}
