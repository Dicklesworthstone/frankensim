//! Progressive ductile softening on the original tetrahedral equilibrium owner.
//!
//! The material supplies native J2 plasticity, a consistent damaged tangent,
//! unilateral volumetric compression and separate physical/numerical accounts.
//! This module binds its characteristic length to original directional element
//! geometry and adaptively subdivides one physical load interval. Rejected
//! attempts never advance the retained plastic, initiation or damage histories.
//! It adds no return map, sparse solver, element deletion or fragment topology.
//!
//! Subdivision controls material damage increments and the original endpoint
//! numerical work loss. Those gates are NOT a temporal-error estimator or proof
//! of a mesh-independent crack path. Load reversals must be separate intervals;
//! reference forces and prescribed displacements vary affinely with the same
//! load factor. Temperatures are held per element throughout each interval.

use core::fmt::Debug;

use fs_exec::Cx;
use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileEvolutionError, DuctileJ2, DuctileJ2Conditions,
    DuctileJ2State,
};

use super::{DuctileStepConditions, DuctileTemperatures};
use crate::TetAssemblyBudget;
use crate::brittle::{TetCrackBandGeometry, TetCrackBandGeometryError};
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};

/// Explicit subdivision, numerical-loss and accepted-history budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileEvolutionControls {
    /// First attempted fraction of the complete physical interval.
    pub initial_fraction: f64,
    /// Smallest allowed fraction, except for one final shorter remainder.
    pub minimum_fraction: f64,
    /// Largest attempted fraction; one means the complete interval.
    pub maximum_fraction: f64,
    /// Bound on published summary records and accepted native increments.
    pub maximum_accepted_steps: usize,
    /// Bound on all native solves, including rejected attempts.
    pub maximum_trials: usize,
    /// Whole-interval budget for endpoint numerical work loss, J. Each attempt
    /// receives this budget times its actual fraction of the physical interval.
    /// This numerical loss is never deposited as heat or fracture energy.
    pub maximum_numerical_loss_j: f64,
    /// Whole-interval absolute external-work accounting tolerance, J.
    pub energy_tolerance_j: f64,
}

/// Original external work and changes in the declared material accounts, J.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DuctileTetEnergy {
    pub external_work_j: f64,
    pub stored_change_j: f64,
    pub plastic_dissipation_j: f64,
    pub damage_release_j: f64,
    pub numerical_loss_j: f64,
    /// Independent diagnostic integral used by the softening calibration. This
    /// includes energetic hardening and is NOT another physical heat/loss term.
    pub degraded_equivalent_plastic_work_j: f64,
    /// External work minus storage change and the three separate loss accounts.
    pub residual_j: f64,
}

/// Bounded observation of one accepted native increment, without retaining a
/// copy of its entire displacement field or reaction vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileEvolutionRecord {
    pub fraction: f64,
    pub load_factor: f64,
    pub duration_s: f64,
    pub maximum_damage: f64,
    pub energy: DuctileTetEnergy,
}

/// One complete accepted physical interval. All returned material histories
/// continue to borrow their original immutable material and length owners.
#[derive(Debug)]
pub struct DuctileTetEvolutionStep<'a, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Final native endpoint; its native increment record describes the last
    /// substep, while `records` describes every accepted physical substep.
    pub mechanics: SmallStrainTetSolution<DuctileJ2State<'a, L>>,
    pub records: Vec<DuctileEvolutionRecord>,
    pub trials: usize,
    pub rejected_trials: usize,
    pub duration_s: f64,
    pub energy: DuctileTetEnergy,
}

/// Reason a physical subinterval needed refinement.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileRefinementCause {
    Native(SmallStrainTetError),
    NumericalLoss { observed_j: f64, allowance_j: f64 },
}

/// Refusal of a complete interval; no partly advanced solution is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileTetEvolutionError {
    Invalid(&'static str),
    Geometry(TetCrackBandGeometryError),
    Material {
        element: usize,
        source: DuctileEvolutionError,
    },
    Native(SmallStrainTetError),
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    NumericalLoss {
        observed_j: f64,
        allowance_j: f64,
    },
    WorkLimit {
        accepted_steps: usize,
        trials: usize,
    },
    MinimumInterval {
        accepted_fraction: f64,
        attempted_fraction: f64,
        cause: DuctileRefinementCause,
    },
    Cancelled,
}

impl core::fmt::Display for DuctileTetEvolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ductile tetrahedral evolution refused: {self:?}")
    }
}
impl std::error::Error for DuctileTetEvolutionError {}

/// Original geometry and material owners for a quasistatic softening specimen.
/// Model references permit source-backed material families to be retained
/// without cloning their calibration cards into an assembly field.
#[derive(Debug)]
pub struct DuctileTetEvolution<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    geometry: &'g TetLoadGeometry3<'g>,
    materials: &'m [&'m DuctileJ2<L>],
    bands: &'g [TetCrackBandGeometry],
    prescribed_m: &'g [(usize, f64)],
    nodal_forces_n: &'g [(usize, f64)],
    budget: TetAssemblyBudget,
    settings: IncrementSettings,
}

impl<'g, 'm, L> DuctileTetEvolution<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Check each band's original connectivity, volume and directional area
    /// against the admitted element, using the existing V/projected-area owner.
    /// The material retains the exact supplied band width. Geometric agreement
    /// does not authorize history transfer after changing the mesh or direction.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        geometry: &'g TetLoadGeometry3<'g>,
        materials: &'m [&'m DuctileJ2<L>],
        bands: &'g [TetCrackBandGeometry],
        prescribed_m: &'g [(usize, f64)],
        nodal_forces_n: &'g [(usize, f64)],
        budget: TetAssemblyBudget,
        settings: IncrementSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, DuctileTetEvolutionError> {
        checkpoint(cx)?;
        let count = geometry.tetrahedra().len();
        if count == 0
            || count > budget.maximum_tetrahedra
            || geometry.nodes_m().len() > budget.maximum_nodes
            || materials.len() != count
            || bands.len() != count
        {
            return Err(DuctileTetEvolutionError::Invalid(
                "geometry, band and material fields must cover every bounded tetrahedron",
            ));
        }
        for (element, (model, band)) in materials.iter().zip(bands).enumerate() {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            let actual = TetCrackBandGeometry::from_element(geometry, element, band.normal())
                .map_err(DuctileTetEvolutionError::Geometry)?;
            // Reconstructing from a unit normal renormalizes it; that operation
            // is not bit-idempotent. Only the derived area and width receive a
            // small relative roundoff allowance. Original topology/volume and
            // the model-to-record length binding remain exact.
            if actual.element() != band.element()
                || actual.nodes() != band.nodes()
                || actual.reference_volume_m3() != band.reference_volume_m3()
                || (actual.projected_area_m2() / band.projected_area_m2() - 1.0).abs()
                    > 64.0 * f64::EPSILON
                || (actual.band_width_m() / band.band_width_m() - 1.0).abs() > 64.0 * f64::EPSILON
                || model.settings().characteristic_length_m != band.band_width_m()
            {
                return Err(DuctileTetEvolutionError::Invalid(
                    "material length and band must belong to the actual original tetrahedron",
                ));
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            geometry,
            materials,
            bands,
            prescribed_m,
            nodal_forces_n,
            budget,
            settings,
        })
    }

    /// Original directional width/area records, in element order.
    pub fn bands(&self) -> &[TetCrackBandGeometry] {
        self.bands
    }

    /// Advance an affine load-factor interval from an accepted native state.
    /// Prescribed displacements, forces and the supplied old state must describe
    /// the same initial load factor. Each retry uses the same accepted history;
    /// only accepted fractions consume physical time. A later refusal discards
    /// the whole proposed interval and leaves all supplied inputs untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_from(
        &self,
        displacement_m: &[f64],
        states: &[DuctileJ2State<'m, L>],
        load_factors: [f64; 2],
        conditions: DuctileStepConditions<'_>,
        controls: DuctileEvolutionControls,
        cx: &Cx<'_>,
    ) -> Result<DuctileTetEvolutionStep<'m, L>, DuctileTetEvolutionError> {
        checkpoint(cx)?;
        validate_controls(controls)?;
        let count = self.materials.len();
        let cover = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(values) => values.len() == count,
            DuctileTemperatures::PerHexElement(_) => false,
        };
        if !cover
            || states.len() != count
            || self.geometry.nodes_m().len().checked_mul(3) != Some(displacement_m.len())
            || displacement_m
                .iter()
                .chain(&load_factors)
                .any(|v| !v.is_finite())
            || !positive(conditions.duration_s)
        {
            return Err(DuctileTetEvolutionError::Invalid(
                "finite bounded original state, positive physical duration and a complete tetrahedral temperature field required",
            ));
        }
        // Check original owner/history and the actual fixed temperature before
        // a refinement loop, so a permanent source mismatch is not retried.
        for (element, (model, state)) in self.materials.iter().zip(states).enumerate() {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            if !model.owns_state(state) || !state.is_valid() {
                return Err(DuctileTetEvolutionError::Invalid(
                    "retained history belongs to a different material or characteristic length",
                ));
            }
            model
                .small_strain_law(DuctileJ2Conditions {
                    duration_s: conditions.duration_s,
                    temperature_k: conditions.temperatures_k.at(element),
                })
                .map_err(|source| DuctileTetEvolutionError::Material { element, source })?;
        }

        let mut current: Option<SmallStrainTetSolution<DuctileJ2State<'m, L>>> = None;
        let mut records = Vec::with_capacity(controls.maximum_accepted_steps.min(64));
        let mut totals = [Sum::default(); 6];
        let mut progress = 0.0_f64;
        let mut fraction = controls.initial_fraction;
        let mut trials = 0;
        let mut rejected = 0;
        while progress < 1.0 {
            checkpoint(cx)?;
            if trials >= controls.maximum_trials || records.len() >= controls.maximum_accepted_steps
            {
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
                .map_or(states, |value| value.element_states.as_slice());
            let old_u = current
                .as_ref()
                .map_or(displacement_m, |value| value.displacement_m.as_flattened());
            let mut laws = Vec::with_capacity(count);
            for (element, model) in self.materials.iter().enumerate() {
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
            let targets = [target];
            let problem = SmallStrainTetProblem {
                nodes_m: self.geometry.nodes_m(),
                tetrahedra: self.geometry.tetrahedra(),
                laws: SmallStrainLawField::PerElement(&laws),
                prescribed_m: self.prescribed_m,
                nodal_forces_n: self.nodal_forces_n,
                load_path: &targets,
                budget: self.budget,
                settings: self.settings,
            };
            trials += 1;
            let proposed = problem.solve_from(old_u, old_states, cx);
            let mechanics = match proposed {
                Ok(candidate) => candidate,
                Err(SmallStrainTetError::Cancelled) => {
                    return Err(DuctileTetEvolutionError::Cancelled);
                }
                Err(error) if retryable(&error) => {
                    fraction = refine(
                        progress,
                        span,
                        controls.minimum_fraction,
                        DuctileRefinementCause::Native(error),
                    )?;
                    rejected += 1;
                    continue;
                }
                Err(error) => return Err(DuctileTetEvolutionError::Native(error)),
            };
            let energy = self.energy(old_u, old_states, &mechanics, target, cx)?;
            let allowance = controls.maximum_numerical_loss_j * span;
            if energy.numerical_loss_j > allowance {
                fraction = refine(
                    progress,
                    span,
                    controls.minimum_fraction,
                    DuctileRefinementCause::NumericalLoss {
                        observed_j: energy.numerical_loss_j,
                        allowance_j: allowance,
                    },
                )?;
                rejected += 1;
                continue;
            }
            let tolerance = controls.energy_tolerance_j * span;
            if energy.residual_j.abs() > tolerance {
                return Err(DuctileTetEvolutionError::Energy {
                    residual_j: energy.residual_j,
                    tolerance_j: tolerance,
                });
            }
            let mut growth = if energy.numerical_loss_j > 0.0 {
                (0.9 * allowance / energy.numerical_loss_j).clamp(0.5, 2.0)
            } else {
                2.0
            };
            let mut maximum_damage = 0.0_f64;
            for (element, (old, next)) in
                old_states.iter().zip(&mechanics.element_states).enumerate()
            {
                if element % 256 == 0 {
                    checkpoint(cx)?;
                }
                let increment = next.damage() - old.damage();
                maximum_damage = maximum_damage.max(next.damage());
                if increment > 0.0 {
                    growth = growth.min(
                        (0.9 * self.materials[element].settings().maximum_damage_increment
                            / increment)
                            .clamp(0.5, 2.0),
                    );
                }
            }
            for (sum, value) in totals.iter_mut().zip(energy.values()) {
                sum.add(value);
            }
            records.push(DuctileEvolutionRecord {
                fraction: end,
                load_factor: target,
                duration_s,
                maximum_damage,
                energy,
            });
            current = Some(mechanics);
            progress = end;
            fraction = (span * growth).clamp(controls.minimum_fraction, controls.maximum_fraction);
        }
        let energy = DuctileTetEnergy::from_values(totals.map(Sum::value));
        if energy
            .values()
            .iter()
            .chain([energy.residual_j].iter())
            .any(|v| !v.is_finite())
        {
            return Err(DuctileTetEvolutionError::Invalid(
                "complete interval energy cannot be represented",
            ));
        }
        if energy.numerical_loss_j > controls.maximum_numerical_loss_j {
            return Err(DuctileTetEvolutionError::NumericalLoss {
                observed_j: energy.numerical_loss_j,
                allowance_j: controls.maximum_numerical_loss_j,
            });
        }
        if energy.residual_j.abs() > controls.energy_tolerance_j {
            return Err(DuctileTetEvolutionError::Energy {
                residual_j: energy.residual_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        checkpoint(cx)?;
        Ok(DuctileTetEvolutionStep {
            mechanics: current.ok_or(DuctileTetEvolutionError::Invalid(
                "no physical increment was accepted",
            ))?,
            records,
            trials,
            rejected_trials: rejected,
            duration_s: conditions.duration_s,
            energy,
        })
    }

    fn energy(
        &self,
        old_u: &[f64],
        old_states: &[DuctileJ2State<'m, L>],
        candidate: &SmallStrainTetSolution<DuctileJ2State<'m, L>>,
        load_factor: f64,
        cx: &Cx<'_>,
    ) -> Result<DuctileTetEnergy, DuctileTetEvolutionError> {
        let mut sums = [Sum::default(); 6];
        let new_u = candidate.displacement_m.as_flattened();
        let record = candidate
            .increments
            .last()
            .ok_or(DuctileTetEvolutionError::Invalid(
                "native endpoint lacks its reaction record",
            ))?;
        if record.reactions_n.len() != self.prescribed_m.len() {
            return Err(DuctileTetEvolutionError::Invalid(
                "native reaction count changed",
            ));
        }
        for (index, &(dof, force)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            sums[0].add((load_factor * force) * (new_u[dof] - old_u[dof]));
        }
        for (index, (&(dof, _), reaction)) in self
            .prescribed_m
            .iter()
            .zip(&record.reactions_n)
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            sums[0].add(reaction * (new_u[dof] - old_u[dof]));
        }
        for (element, ((old, next), &volume)) in old_states
            .iter()
            .zip(&candidate.element_states)
            .zip(self.geometry.element_volumes_m3())
            .enumerate()
        {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            let plastic = next.dissipated_plastic_j_m3() - old.dissipated_plastic_j_m3();
            let damage = next.dissipated_damage_j_m3() - old.dissipated_damage_j_m3();
            let numerical = next.numerical_dissipation_j_m3() - old.numerical_dissipation_j_m3();
            if plastic < 0.0 || damage < 0.0 || numerical < 0.0 {
                return Err(DuctileTetEvolutionError::Invalid(
                    "native cumulative loss decreased",
                ));
            }
            sums[1].add(volume * (next.stored_energy_j_m3() - old.stored_energy_j_m3()));
            sums[2].add(volume * plastic);
            sums[3].add(volume * damage);
            sums[4].add(volume * numerical);
            sums[5].add(
                self.bands[element].projected_area_m2()
                    * (next.degraded_equivalent_plastic_work_j_m2()
                        - old.degraded_equivalent_plastic_work_j_m2()),
            );
        }
        let result = DuctileTetEnergy::from_values(sums.map(Sum::value));
        if result
            .values()
            .iter()
            .chain([result.residual_j].iter())
            .any(|v| !v.is_finite())
        {
            return Err(DuctileTetEvolutionError::Invalid(
                "whole-solid energy cannot be represented",
            ));
        }
        Ok(result)
    }
}

impl DuctileTetEnergy {
    fn values(self) -> [f64; 6] {
        [
            self.external_work_j,
            self.stored_change_j,
            self.plastic_dissipation_j,
            self.damage_release_j,
            self.numerical_loss_j,
            self.degraded_equivalent_plastic_work_j,
        ]
    }
    fn from_values(v: [f64; 6]) -> Self {
        Self {
            external_work_j: v[0],
            stored_change_j: v[1],
            plastic_dissipation_j: v[2],
            damage_release_j: v[3],
            numerical_loss_j: v[4],
            degraded_equivalent_plastic_work_j: v[5],
            residual_j: (((v[0] - v[1]) - v[2]) - v[3]) - v[4],
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), DuctileTetEvolutionError> {
    cx.checkpoint()
        .map_err(|_| DuctileTetEvolutionError::Cancelled)
}
fn validate_controls(c: DuctileEvolutionControls) -> Result<(), DuctileTetEvolutionError> {
    if !positive(c.minimum_fraction)
        || !positive(c.initial_fraction)
        || !positive(c.maximum_fraction)
        || c.minimum_fraction > c.initial_fraction
        || c.initial_fraction > c.maximum_fraction
        || c.maximum_fraction > 1.0
        || c.maximum_accepted_steps == 0
        || c.maximum_trials == 0
        || !c.maximum_numerical_loss_j.is_finite()
        || c.maximum_numerical_loss_j < 0.0
        || !positive(c.energy_tolerance_j)
    {
        return Err(DuctileTetEvolutionError::Invalid(
            "ordered positive interval fractions and finite nonnegative work budgets required",
        ));
    }
    Ok(())
}
fn retryable(error: &SmallStrainTetError) -> bool {
    match error {
        SmallStrainTetError::NewtonStalled { .. }
        | SmallStrainTetError::LinearSolveFailed { .. } => true,
        SmallStrainTetError::InvalidInput { what } => what.contains("constitutive"),
        _ => false,
    }
}
fn refine(
    progress: f64,
    span: f64,
    minimum: f64,
    cause: DuctileRefinementCause,
) -> Result<f64, DuctileTetEvolutionError> {
    let smaller = 0.5 * span;
    if smaller < minimum || progress + smaller <= progress {
        Err(DuctileTetEvolutionError::MinimumInterval {
            accepted_fraction: progress,
            attempted_fraction: span,
            cause,
        })
    } else {
        Ok(smaller)
    }
}
