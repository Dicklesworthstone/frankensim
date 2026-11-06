//! Progressive ductile softening on the native backward-Euler momentum owner.
//!
//! The original material, geometry, crack-band lengths and consistent reference
//! mass survive every accepted fraction. Rejected trials advance neither motion
//! nor plastic/damage history. Affine endpoint loads are evaluated at each actual
//! accepted time, with the same duration driving inertia and the ductile law.
//!
//! Material numerical loss and kinetic numerical loss have separate budgets and
//! are never physical heating or fracture work. Small strain/rotation, held
//! temperature, fixed mass and prescribed crack-band direction remain the parent
//! model's assumptions. Adaptive work/damage limits do not certify trajectory
//! accuracy, resolve cracks, provide contact, or predict fragment motion.

use super::*;
use crate::inertia::{InertialEnergy, P1ReferenceMass};
use crate::plastic3::dynamics::BackwardEulerStep;

/// Independent material and inertial limits over one physical interval.
#[derive(Clone, Copy, Debug)]
pub struct DuctileMotionControls {
    /// Energy tolerance gates material, kinetic and combined work identities
    /// separately; their errors cannot cancel to conceal a failed account.
    pub continuation: DuctileEvolutionControls,
    /// Maximum sum of `1/2 delta_v M delta_v`, J, over the complete interval.
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Independent external, constitutive and inertial work accounts, in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileMotionEnergy {
    /// Applied endpoint forces and inertial support reactions dot actual du.
    pub external_work_j: f64,
    /// Its `external_work_j` is the work delivered to the constitutive response:
    /// complete external work minus the separately retained inertial work.
    pub material: DuctileTetEnergy,
    pub inertia: InertialEnergy,
    /// External work minus material storage/loss, kinetic change and kinetic loss.
    pub balance_error_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileMotionRecord {
    pub fraction: f64,
    pub load_factor: f64,
    pub duration_s: f64,
    pub maximum_damage: f64,
    pub energy: DuctileMotionEnergy,
}

/// Native endpoint and bounded scalar history. The caller's old inputs survive
/// whole-interval refusal unchanged, even after internally accepted fractions.
#[derive(Debug)]
pub struct DuctileMotionStep<'a, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Final native substep; the top-level energy aggregates the whole interval.
    pub mechanics: BackwardEulerStep<DuctileJ2State<'a, L>>,
    pub records: Vec<DuctileMotionRecord>,
    pub trials: usize,
    pub rejected_trials: usize,
    pub duration_s: f64,
    pub energy: DuctileMotionEnergy,
}

/// Immutable original specimen and mass. Density is admitted once; no endpoint
/// temperature, damage or reduced stiffness silently changes the mass inventory.
#[derive(Debug)]
pub struct DuctileTetMotion<'s, 'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    specimen: &'s DuctileTetEvolution<'g, 'm, L>,
    mass: P1ReferenceMass,
}

impl<'s, 'g, 'm, L> DuctileTetMotion<'s, 'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Bind positive reference densities in original tetrahedron order.
    /// The original specimen's admitted assembly/geometry bounds apply.
    pub fn new(
        specimen: &'s DuctileTetEvolution<'g, 'm, L>,
        reference_density_kg_m3: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Self, DuctileTetEvolutionError> {
        checkpoint(cx)?;
        let geometry = specimen.geometry;
        let ndofs =
            geometry
                .nodes_m()
                .len()
                .checked_mul(3)
                .ok_or(DuctileTetEvolutionError::Invalid(
                    "reference motion component count overflow",
                ))?;
        let mass = P1ReferenceMass::new(
            geometry
                .tetrahedra()
                .iter()
                .copied()
                .zip(geometry.element_volumes_m3().iter().copied()),
            reference_density_kg_m3,
            ndofs,
            cx,
        )
        .map_err(|error| match SmallStrainTetError::from(error) {
            SmallStrainTetError::Cancelled => DuctileTetEvolutionError::Cancelled,
            error => DuctileTetEvolutionError::Native(error),
        })?;
        Ok(Self { specimen, mass })
    }

    /// Advance one affine load interval. Old displacement, velocity and history
    /// are borrowed. Refine from the returned complete endpoint on success only.
    /// Forces and prescribed motion use the specimen's original load pattern.
    /// Initial prescribed displacements must equal that pattern times the
    /// starting load factor; a different support origin is not an affine path.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_from(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        states: &[DuctileJ2State<'m, L>],
        load_factors: [f64; 2],
        conditions: DuctileStepConditions<'_>,
        controls: DuctileMotionControls,
        cx: &Cx<'_>,
    ) -> Result<DuctileMotionStep<'m, L>, DuctileTetEvolutionError> {
        checkpoint(cx)?;
        let limits = controls.continuation;
        validate_controls(limits)?;
        let specimen = self.specimen;
        let count = specimen.materials.len();
        let cover = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(values) => values.len() == count,
            DuctileTemperatures::PerHexElement(_) => false,
        };
        if !cover
            || states.len() != count
            || specimen.geometry.nodes_m().len().checked_mul(3) != Some(displacement_m.len())
            || velocity_mps.len() != specimen.geometry.nodes_m().len()
            || displacement_m
                .iter()
                .chain(&load_factors)
                .chain(velocity_mps.as_flattened())
                .any(|v| !v.is_finite())
            || !positive(conditions.duration_s)
            || !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(DuctileTetEvolutionError::Invalid(
                "complete finite motion/history/temperature fields and nonnegative kinetic-loss budget required",
            ));
        }
        for (index, &(dof, value)) in specimen.prescribed_m.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            let expected = load_factors[0] * value;
            if !expected.is_finite() || displacement_m.get(dof) != Some(&expected) {
                return Err(DuctileTetEvolutionError::Invalid(
                    "initial prescribed displacement must match the starting load factor exactly",
                ));
            }
        }
        let mut laws = specimen.admit_history(displacement_m, states, conditions, cx)?;

        let mut current: Option<BackwardEulerStep<DuctileJ2State<'m, L>>> = None;
        let mut records: Vec<DuctileMotionRecord> =
            Vec::with_capacity(limits.maximum_accepted_steps.min(64));
        let mut prescribed = vec![(0, 0.0); specimen.prescribed_m.len()];
        let mut forces = vec![(0, 0.0); specimen.nodal_forces_n.len()];
        // Material's six original channels, external work, and three kinetic channels.
        let mut totals = [Sum::default(); 10];
        let mut progress = 0.0_f64;
        let mut fraction = limits.initial_fraction;
        let mut trials = 0;
        let mut rejected = 0;
        while progress < 1.0 {
            checkpoint(cx)?;
            if trials >= limits.maximum_trials || records.len() >= limits.maximum_accepted_steps {
                return Err(DuctileTetEvolutionError::WorkLimit {
                    accepted_steps: records.len(),
                    trials,
                });
            }
            let end = (progress + fraction).min(1.0);
            let span = end - progress;
            let duration_s = span * conditions.duration_s;
            if end <= progress || !positive(duration_s) {
                return Err(DuctileTetEvolutionError::Invalid(
                    "physical subinterval cannot be represented",
                ));
            }
            let target = (1.0 - end) * load_factors[0] + end * load_factors[1];
            let old_states = current
                .as_ref()
                .map_or(states, |v| v.solution.element_states.as_slice());
            let old_u = current
                .as_ref()
                .map_or(displacement_m, |v| v.solution.displacement_m.as_flattened());
            let old_v = current
                .as_ref()
                .map_or(velocity_mps, |v| v.velocity_mps.as_slice());
            laws.clear();
            for (element, model) in specimen.materials.iter().enumerate() {
                if element % 256 == 0 {
                    checkpoint(cx)?;
                }
                let model: &'m DuctileJ2<L> = *model;
                laws.push(
                    model
                        .small_strain_law(DuctileJ2Conditions {
                            duration_s,
                            temperature_k: conditions.temperatures_k.at(element),
                        })
                        .map_err(|source| DuctileTetEvolutionError::Material { element, source })?,
                );
            }
            // The native inertia API admits actual endpoint values at factor 1.
            for (index, (destination, &(dof, value))) in prescribed
                .iter_mut()
                .zip(specimen.prescribed_m)
                .chain(forces.iter_mut().zip(specimen.nodal_forces_n))
                .enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                *destination = (dof, target * value);
            }
            let problem = SmallStrainTetProblem {
                nodes_m: specimen.geometry.nodes_m(),
                tetrahedra: specimen.geometry.tetrahedra(),
                laws: SmallStrainLawField::PerElement(&laws),
                prescribed_m: &prescribed,
                nodal_forces_n: &forces,
                load_path: &[1.0],
                budget: specimen.budget,
                settings: specimen.settings,
            };
            trials += 1;
            let candidate = match problem.advance_backward_euler_with_mass(
                old_u, old_v, old_states, &self.mass, duration_s, cx,
            ) {
                Ok(value) => value,
                Err(SmallStrainTetError::Cancelled) => {
                    return Err(DuctileTetEvolutionError::Cancelled);
                }
                Err(error) if retryable(&error) => {
                    fraction = refine(
                        progress,
                        span,
                        limits.minimum_fraction,
                        DuctileRefinementCause::Native(error),
                    )?;
                    rejected += 1;
                    continue;
                }
                Err(error) => return Err(DuctileTetEvolutionError::Native(error)),
            };
            let mut material =
                specimen.energy(old_u, old_states, &candidate.solution, target, cx)?;
            let external_work_j = material.external_work_j;
            material.external_work_j -= candidate.energy.inertial_work_j;
            // Recompute from the displayed material channels after removing
            // inertia. Subtracting inertia from an already formed residual
            // could hide material work lost to cancellation against large K.
            material = DuctileTetEnergy::from_values(material.values());
            let energy = DuctileMotionEnergy {
                external_work_j,
                material,
                inertia: candidate.energy,
                balance_error_j: material.residual_j + candidate.energy.balance_error_j,
            };
            let allowances = [
                limits.maximum_numerical_loss_j * span,
                controls.maximum_numerical_kinetic_loss_j * span,
            ];
            let losses = [material.numerical_loss_j, energy.inertia.numerical_loss_j];
            if let Some((&observed_j, &allowance_j)) = losses
                .iter()
                .zip(&allowances)
                .find(|(loss, allowance)| loss > allowance)
            {
                fraction = refine(
                    progress,
                    span,
                    limits.minimum_fraction,
                    DuctileRefinementCause::NumericalLoss {
                        observed_j,
                        allowance_j,
                    },
                )?;
                rejected += 1;
                continue;
            }
            let tolerance = limits.energy_tolerance_j * span;
            for residual_j in [
                material.residual_j,
                energy.inertia.balance_error_j,
                energy.balance_error_j,
            ] {
                if !residual_j.is_finite() || residual_j.abs() > tolerance {
                    return Err(DuctileTetEvolutionError::Energy {
                        residual_j,
                        tolerance_j: tolerance,
                    });
                }
            }
            let mut growth = 2.0_f64;
            for (loss, allowance) in losses.into_iter().zip(allowances) {
                if loss > 0.0 {
                    growth = growth.min((0.9 * allowance / loss).clamp(0.5, 2.0));
                }
            }
            let mut maximum_damage = 0.0_f64;
            for (element, (old, next)) in old_states
                .iter()
                .zip(&candidate.solution.element_states)
                .enumerate()
            {
                if element % 256 == 0 {
                    checkpoint(cx)?;
                }
                maximum_damage = maximum_damage.max(next.damage());
                let increment = next.damage() - old.damage();
                if increment > 0.0 {
                    growth = growth.min(
                        (0.9 * specimen.materials[element]
                            .settings()
                            .maximum_damage_increment
                            / increment)
                            .clamp(0.5, 2.0),
                    );
                }
            }
            for (sum, value) in totals[..6].iter_mut().zip(material.values()) {
                sum.add(value);
            }
            totals[6].add(external_work_j);
            totals[7].add(energy.inertia.kinetic_change_j);
            totals[8].add(energy.inertia.numerical_loss_j);
            totals[9].add(energy.inertia.inertial_work_j);
            records.push(DuctileMotionRecord {
                fraction: end,
                load_factor: target,
                duration_s,
                maximum_damage,
                energy,
            });
            current = Some(candidate);
            progress = end;
            fraction = (span * growth).clamp(limits.minimum_fraction, limits.maximum_fraction);
        }
        let values = totals.map(Sum::value);
        if values.iter().any(|v| !v.is_finite()) {
            return Err(DuctileTetEvolutionError::Invalid(
                "complete motion interval energy cannot be represented",
            ));
        }
        let first = records.first().ok_or(DuctileTetEvolutionError::Invalid(
            "no motion increment was accepted",
        ))?;
        let last = records.last().ok_or(DuctileTetEvolutionError::Invalid(
            "no motion increment was accepted",
        ))?;
        let material = DuctileTetEnergy::from_values([
            values[0], values[1], values[2], values[3], values[4], values[5],
        ]);
        let inertia = InertialEnergy {
            reference_mass_kg: first.energy.inertia.reference_mass_kg,
            kinetic_before_j: first.energy.inertia.kinetic_before_j,
            kinetic_after_j: last.energy.inertia.kinetic_after_j,
            kinetic_change_j: values[7],
            numerical_loss_j: values[8],
            inertial_work_j: values[9],
            balance_error_j: (values[9] - values[7]) - values[8],
        };
        let energy = DuctileMotionEnergy {
            external_work_j: values[6],
            material,
            inertia,
            balance_error_j: material.residual_j + inertia.balance_error_j,
        };
        for (observed_j, allowance_j) in [
            (material.numerical_loss_j, limits.maximum_numerical_loss_j),
            (
                inertia.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            ),
        ] {
            if observed_j > allowance_j {
                return Err(DuctileTetEvolutionError::NumericalLoss {
                    observed_j,
                    allowance_j,
                });
            }
        }
        for residual_j in [
            material.residual_j,
            inertia.balance_error_j,
            energy.balance_error_j,
        ] {
            if !residual_j.is_finite() || residual_j.abs() > limits.energy_tolerance_j {
                return Err(DuctileTetEvolutionError::Energy {
                    residual_j,
                    tolerance_j: limits.energy_tolerance_j,
                });
            }
        }
        checkpoint(cx)?;
        Ok(DuctileMotionStep {
            mechanics: current.ok_or(DuctileTetEvolutionError::Invalid(
                "no motion increment was accepted",
            ))?,
            records,
            trials,
            rejected_trials: rejected,
            duration_s: conditions.duration_s,
            energy,
        })
    }
}
