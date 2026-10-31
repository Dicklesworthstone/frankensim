//! A physical creep/plastic interval on the native P1 tetrahedral Newton solve.
//!
//! This adapter keeps the original uniform/per-element constitutive field,
//! tensor conventions, sparse assembly, supports and globalization. It admits
//! exactly one load endpoint for one common physical duration. Repeated Newton
//! trials borrow the same accepted history; only the returned proposal contains
//! the next history. Failed equilibrium, material, work or cancellation checks
//! leave every caller-owned value unchanged.
//!
//! Elastic/hardening storage coefficients and the reference mesh are fixed
//! across an interval. Temperature may change kinetics and yield stress, but
//! changing stored-energy coefficients requires a parameter-work accounting law.
//! Physical creep/plastic dissipation is separate from backward-Euler loss;
//! only the former is available for thermal coupling. Inertia, finite strain,
//! primary/tertiary creep and service-life qualification are not implied.

use fs_material::creep::{CreepEnergy, CreepError, CreepLaw, CreepState, CreepTolerances};

use super::*;

/// Absolute local constitutive and whole-body work gates. These bound the
/// discrete numerical identities, not error of a measured creep prediction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepTetTolerances {
    pub material: CreepTolerances,
    pub energy_j: f64,
}

/// Work through actual displacement and the native constitutive energy split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepTetEnergy {
    pub initial_stored_j: f64,
    pub final_stored_j: f64,
    pub applied_work_j: f64,
    pub support_work_j: f64,
    pub physical_dissipation_j: f64,
    pub numerical_dissipation_j: f64,
    /// Applied plus support work minus independently integrated stress work.
    pub equilibrium_work_residual_j: f64,
    /// Applied plus support work minus storage change and both loss channels.
    pub balance_residual_j: f64,
}

/// One fully checked proposal, ready for atomic acceptance by the caller.
#[derive(Clone, Debug, PartialEq)]
pub struct CreepTetStep<S = CreepState> {
    pub solution: SmallStrainTetSolution<S>,
    pub duration_s: f64,
    /// Per-reference-volume accounts [J/m³], in tetrahedron order.
    pub element_energy: Vec<CreepEnergy>,
    /// The native reference element volumes; multiply physical dissipation
    /// densities by these volumes when transferring heat to a thermal owner.
    pub element_volume_m3: Vec<f64>,
    pub energy: CreepTetEnergy,
}

#[derive(Debug)]
pub enum CreepTetError {
    Invalid(&'static str),
    Mechanical(SmallStrainTetError),
    Material {
        element: usize,
        source: CreepError,
    },
    /// The constitutive field and its audited native owner disagree.
    ConstitutiveMismatch {
        element: usize,
    },
    EnergyResidual {
        residual_j: f64,
        tolerance_j: f64,
    },
}

impl core::fmt::Display for CreepTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "FS-SOLID-CREEP-INPUT: {reason}"),
            Self::Mechanical(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "FS-SOLID-CREEP-MATERIAL: element {element}: {source}")
            }
            Self::ConstitutiveMismatch { element } => {
                write!(
                    f,
                    "FS-SOLID-CREEP-STATE: inconsistent element {element} endpoint"
                )
            }
            Self::EnergyResidual {
                residual_j,
                tolerance_j,
            } => write!(
                f,
                "FS-SOLID-CREEP-WORK: residual {residual_j} J exceeds {tolerance_j} J"
            ),
        }
    }
}
impl std::error::Error for CreepTetError {}
impl From<SmallStrainTetError> for CreepTetError {
    fn from(value: SmallStrainTetError) -> Self {
        Self::Mechanical(value)
    }
}

fn finite(value: f64) -> Result<f64, CreepTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CreepTetError::Invalid(
            "unrepresentable creep work or storage",
        ))
    }
}

fn audit_residual(value: f64, tolerance: f64) -> Result<(), CreepTetError> {
    finite(value)?;
    if value.abs() > tolerance {
        Err(CreepTetError::EnergyResidual {
            residual_j: value,
            tolerance_j: tolerance,
        })
    } else {
        Ok(())
    }
}

fn admit_small_kinematics(
    element: &Element,
    displacement: &[f64],
    limit: f64,
) -> Result<(), CreepTetError> {
    // Symmetric strain alone cannot detect an arbitrarily large infinitesimal
    // rigid rotation. The complete displacement gradient must remain small.
    for component in 0..3 {
        for direction in 0..3 {
            let gradient = element
                .nodes
                .iter()
                .enumerate()
                .fold(0.0, |sum, (a, &node)| {
                    displacement[3 * node + component].mul_add(element.gradients[a][direction], sum)
                });
            if !gradient.is_finite() || gradient.abs() > limit {
                return Err(CreepTetError::Invalid(
                    "displacement gradient exceeds the small-strain/rotation envelope",
                ));
            }
        }
    }
    Ok(())
}

impl<L: CreepLaw> SmallStrainTetProblem<'_, L> {
    /// Solve and audit one physical creep interval from accepted displacement
    /// and internal variables. `load_path` must have exactly one absolute
    /// endpoint; every element law must represent the same positive duration.
    ///
    /// The previous history must have the same elastic coefficients and
    /// reference geometry. Temperature/rate changes preserve that history;
    /// callers must account separately for changing elastic storage or geometry.
    /// All Newton trials see the old state. Returned stresses are the stresses
    /// that balanced the endpoint loads, not a second constitutive advance.
    pub fn solve_creep_from(
        &self,
        displacement_m: &[f64],
        element_states: &[L::State],
        tolerances: CreepTetTolerances,
        cx: &Cx<'_>,
    ) -> Result<CreepTetStep<L::State>, CreepTetError> {
        checkpoint(cx)?;
        if self.load_path.len() != 1 {
            return Err(CreepTetError::Invalid(
                "one physical creep interval requires exactly one load endpoint",
            ));
        }
        if [
            tolerances.energy_j,
            tolerances.material.stress_pa,
            tolerances.material.energy_j_m3,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(CreepTetError::Invalid(
                "all work and flow tolerances must be finite and positive",
            ));
        }
        let elements = self.validate()?;
        self.check(displacement_m, element_states)?;
        let duration_s = self.laws.at(0).kinetics().duration_s();
        let mut old_strains = Vec::with_capacity(elements.len());
        for (index, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            if !element.volume.is_finite() || element.volume <= 0.0 {
                return Err(CreepTetError::Invalid(
                    "each element requires a finite positive reference volume",
                ));
            }
            let law = self.laws.at(index);
            if law.kinetics().duration_s() != duration_s {
                return Err(CreepTetError::Invalid(
                    "all element laws must represent the same physical duration",
                ));
            }
            admit_small_kinematics(
                element,
                displacement_m,
                law.kinetics().elastic().strain_limit,
            )?;
            let old_strain = Self::strain(element, displacement_m);
            law.validate_creep_history(&old_strain, &element_states[index])
                .map_err(|source| CreepTetError::Material {
                    element: index,
                    source,
                })?;
            old_strains.push(old_strain);
        }
        let solution = self.solve_from(displacement_m, element_states, cx)?;
        let next_u: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
        let mut element_energy = Vec::with_capacity(elements.len());
        let mut element_volume_m3 = Vec::with_capacity(elements.len());
        let mut energy = CreepTetEnergy {
            initial_stored_j: 0.0,
            final_stored_j: 0.0,
            applied_work_j: 0.0,
            support_work_j: 0.0,
            physical_dissipation_j: 0.0,
            numerical_dissipation_j: 0.0,
            equilibrium_work_residual_j: 0.0,
            balance_residual_j: 0.0,
        };
        let mut material_work = 0.0;
        let mut storage_change = 0.0;
        for (index, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let law = self.laws.at(index);
            admit_small_kinematics(element, &next_u, law.kinetics().elastic().strain_limit)?;
            let increment = law
                .creep_increment(
                    &old_strains[index],
                    &element_states[index],
                    &solution.element_strain[index],
                    tolerances.material,
                )
                .map_err(|source| CreepTetError::Material {
                    element: index,
                    source,
                })?;
            if increment.state != solution.element_states[index]
                || increment
                    .stress_pa
                    .iter()
                    .zip(&solution.element_stress[index])
                    .any(|(a, b)| !b.is_finite() || (a - b).abs() > tolerances.material.stress_pa)
            {
                return Err(CreepTetError::ConstitutiveMismatch { element: index });
            }
            let account = increment.energy;
            let volume = element.volume;
            energy.initial_stored_j =
                finite(energy.initial_stored_j + volume * account.initial_stored_j_m3)?;
            energy.final_stored_j =
                finite(energy.final_stored_j + volume * account.final_stored_j_m3)?;
            energy.physical_dissipation_j =
                finite(energy.physical_dissipation_j + volume * account.physical_dissipation_j_m3)?;
            energy.numerical_dissipation_j = finite(
                energy.numerical_dissipation_j + volume * account.numerical_dissipation_j_m3,
            )?;
            material_work = finite(material_work + volume * account.mechanical_work_j_m3)?;
            storage_change = finite(
                storage_change + volume * (account.final_stored_j_m3 - account.initial_stored_j_m3),
            )?;
            element_energy.push(account);
            element_volume_m3.push(volume);
        }
        for &(dof, force) in self.nodal_forces_n {
            energy.applied_work_j = finite(
                energy.applied_work_j
                    + self.load_path[0] * force * (next_u[dof] - displacement_m[dof]),
            )?;
        }
        for (&(dof, _), &reaction) in self
            .prescribed_m
            .iter()
            .zip(&solution.increments[0].reactions_n)
        {
            energy.support_work_j =
                finite(energy.support_work_j + reaction * (next_u[dof] - displacement_m[dof]))?;
        }
        let mechanical_work = finite(energy.applied_work_j + energy.support_work_j)?;
        energy.equilibrium_work_residual_j = mechanical_work - material_work;
        energy.balance_residual_j = mechanical_work
            - storage_change
            - energy.physical_dissipation_j
            - energy.numerical_dissipation_j;
        audit_residual(energy.equilibrium_work_residual_j, tolerances.energy_j)?;
        audit_residual(energy.balance_residual_j, tolerances.energy_j)?;
        checkpoint(cx)?;
        Ok(CreepTetStep {
            solution,
            duration_s,
            element_energy,
            element_volume_m3,
            energy,
        })
    }
}
