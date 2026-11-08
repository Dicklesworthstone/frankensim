//! Progressive ductility on the original CST/DKT shell equilibrium owner.
//!
//! The existing layered section supplies the native J2 return, transverse
//! equilibrium and nonsymmetric damaged tangent at every surface/thickness
//! point. This module gates one accepted quasistatic or backward-Euler endpoint
//! using the material's own storage, plastic loss, damage release and numerical
//! loss. It introduces no scalar stiffness approximation or additional physical
//! update.
//!
//! Immutable reference drilling and consistent mass remain numerical/mechanical
//! owners of their own work. Neither their losses nor the damage work-scale
//! diagnostic becomes plastic heat. Every prepared layer must use the same
//! actual duration, equal to the inertial interval when motion is advanced.
//!
//! Large damage increments refuse through the original material gate; callers
//! must refine physical load/time intervals without changing retained history.
//! These balance gates do not estimate temporal error. Band lengths retain their
//! explicitly authored geometry meaning: thickness is not substituted for an
//! in-plane crack width. No finite forming, delamination, crack topology, erosion,
//! contact generation, heat partition or experimental qualification is supplied.

use core::fmt::Debug;
use std::cell::RefCell;

use fs_exec::Cx;
use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileJ2Law, DuctileJ2State,
};
use fs_material::plastic::reduction::PlaneStressWorkspace;

use super::drilling::{PlasticShellReferenceDrilling, ShellDrillingEnergy};
use super::dynamics::{PlasticShellReferenceInertia, ShellInertialEnergy};
use super::{
    Element, Operator, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, invalid,
};
use crate::plastic3::SmallStrainTetError;
use crate::plastic3::increment::IncrementOperator;

/// Independent absolute work-account allowances, not material calibration or
/// time-error estimates. Each original material retains its own native gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileShellEnergyTolerances {
    pub section_work_j_m2: f64,
    pub shell_work_j: f64,
}

/// Native mechanical work over the original shell, J. Every loss keeps its
/// constitutive meaning; the work-scale diagnostic is excluded from the sum.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DuctileShellEnergy {
    pub initial_stored_j: f64,
    pub final_stored_j: f64,
    /// Integrated pointwise changes, avoiding subtraction of large totals.
    pub stored_change_j: f64,
    pub plastic_dissipation_j: f64,
    pub damage_release_j: f64,
    pub material_numerical_dissipation_j: f64,
    pub material_work_j: f64,
    pub generalized_work_j: f64,
    pub transverse_residual_work_j: f64,
    /// Volume integral of `Gp * delta_damage / h`, not another physical loss.
    pub degraded_equivalent_plastic_work_j: f64,
    /// Original numerical spin gauge, kept outside material storage and heat.
    pub drilling: ShellDrillingEnergy,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Signed endpoint inertia work; zero for the quasistatic entry point.
    pub inertial_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    /// External plus transverse work minus material storage/loss, drilling
    /// storage/loss and inertia work. Motion reports kinetic accounts separately.
    pub balance_residual_j: f64,
}

/// One native quasistatic endpoint, published only after all work gates pass.
#[derive(Debug)]
pub struct DuctileShellEnergyStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub duration_s: f64,
    pub energy: DuctileShellEnergy,
}

/// One endpoint motion on unchanged reference mass and drilling stiffness.
#[derive(Debug)]
pub struct DuctileShellDynamicEnergyStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub velocity: Vec<[f64; 6]>,
    pub duration_s: f64,
    pub inertia: ShellInertialEnergy,
    pub energy: DuctileShellEnergy,
}

impl<S> DuctileShellDynamicEnergyStep<S> {
    /// Complete residual including the independently integrated kinetic work.
    pub fn balance_residual_j(&self) -> f64 {
        self.energy.balance_residual_j + self.inertia.total.balance_error_j
    }
}

fn gate(what: &str, residual: f64, tolerance: f64) -> Result<(), SmallStrainTetError> {
    if !residual.is_finite() || residual.abs() > tolerance {
        return Err(invalid(format!(
            "ductile shell {what}: {residual} J exceeds {tolerance} J"
        )));
    }
    Ok(())
}

fn validate_tolerances(value: DuctileShellEnergyTolerances) -> Result<(), SmallStrainTetError> {
    if [value.section_work_j_m2, value.shell_work_j]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid(
            "positive finite ductile shell work allowances required",
        ));
    }
    Ok(())
}

impl<'m, L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq>
    SmallStrainShellProblem<'_, DuctileJ2Law<'m, L>>
{
    /// Advance one physical quasistatic interval with the native damaged
    /// tangent and the original numerical drilling gauge. Exactly one load
    /// factor is required. Every section law represents the same duration;
    /// temperature may vary between original material layers.
    ///
    /// Incoming coordinates and histories remain unchanged on every refusal.
    /// Rebuild only the borrowed prepared laws when refining an interval;
    /// their original material owners and characteristic lengths must survive.
    pub fn solve_ductile_with_energy(
        &self,
        initial: Option<&SmallStrainShellSolution<DuctileJ2State<'m, L>>>,
        drilling: &PlasticShellReferenceDrilling<'_>,
        tolerances: DuctileShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<DuctileShellEnergyStep<DuctileJ2State<'m, L>>, SmallStrainTetError> {
        checkpoint(cx)?;
        validate_tolerances(tolerances)?;
        if self.load_path.len() != 1 {
            return Err(invalid(
                "ductile shell energy requires one physical load interval",
            ));
        }
        let elements = self.prepare_with_reference_drilling(Some(drilling), cx)?;
        let duration_s = self.ductile_duration(&elements, cx)?;
        let solution = self.solve_prepared(
            &elements,
            initial.map(|old| (old.coordinates.as_slice(), old.section_states.as_slice())),
            None,
            cx,
        )?;
        let energy =
            self.audit_ductile_solution(&elements, initial, &solution, 0.0, tolerances, cx)?;
        checkpoint(cx)?;
        Ok(DuctileShellEnergyStep {
            solution,
            duration_s,
            energy,
        })
    }

    /// Advance one backward-Euler momentum interval with exactly the same
    /// duration in every native ductile law. `load_path` must be `[1.0]`.
    /// Physical translation/director kinetics, numerical drilling-rate kinetics,
    /// material endpoint loss and damage release remain separate receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_ductile_backward_euler_with_energy(
        &self,
        initial: Option<&SmallStrainShellSolution<DuctileJ2State<'m, L>>>,
        velocity: &[[f64; 6]],
        inertia: &PlasticShellReferenceInertia<'_>,
        drilling: &PlasticShellReferenceDrilling<'_>,
        dt_s: f64,
        tolerances: DuctileShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<DuctileShellDynamicEnergyStep<DuctileJ2State<'m, L>>, SmallStrainTetError> {
        checkpoint(cx)?;
        validate_tolerances(tolerances)?;
        if self.load_path != [1.0] || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid(
                "one endpoint shell load and positive physical duration required",
            ));
        }
        let elements = self.prepare_with_reference_drilling(Some(drilling), cx)?;
        let duration_s = self.ductile_duration(&elements, cx)?;
        if duration_s != dt_s {
            return Err(invalid(
                "ductile shell material and inertia durations must agree",
            ));
        }
        let step =
            self.advance_backward_euler_prepared(&elements, initial, velocity, inertia, dt_s, cx)?;
        let energy = self.audit_ductile_solution(
            &elements,
            initial,
            &step.solution,
            step.inertia.total.inertial_work_j,
            tolerances,
            cx,
        )?;
        for (what, residual) in [
            ("inertial work", step.inertia.total.balance_error_j),
            (
                "physical kinetic work",
                step.inertia.physical.inertial_work_j
                    - step.inertia.physical.kinetic_change_j
                    - step.inertia.physical.numerical_loss_j,
            ),
            (
                "drilling kinetic work",
                step.inertia.drilling_gauge.inertial_work_j
                    - step.inertia.drilling_gauge.kinetic_change_j
                    - step.inertia.drilling_gauge.numerical_loss_j,
            ),
            (
                "inertial partition",
                step.inertia.total.inertial_work_j
                    - step.inertia.physical.inertial_work_j
                    - step.inertia.drilling_gauge.inertial_work_j,
            ),
            (
                "complete dynamic balance",
                energy.balance_residual_j + step.inertia.total.balance_error_j,
            ),
        ] {
            gate(what, residual, tolerances.shell_work_j)?;
        }
        checkpoint(cx)?;
        Ok(DuctileShellDynamicEnergyStep {
            solution: step.solution,
            velocity: step.velocity,
            duration_s,
            inertia: step.inertia,
            energy,
        })
    }

    fn ductile_duration(
        &self,
        elements: &[Element],
        cx: &Cx<'_>,
    ) -> Result<f64, SmallStrainTetError> {
        let mut duration = None;
        for e in 0..elements.len() {
            checkpoint(cx)?;
            for layer in self.sections.at(e).layers() {
                let value = layer.law.conditions().duration_s;
                if duration.is_some_and(|original| original != value) {
                    return Err(invalid(
                        "all ductile shell layers must share one physical interval",
                    ));
                }
                duration = Some(value);
            }
        }
        duration.ok_or_else(|| invalid("ductile shell requires original material sections"))
    }

    #[allow(clippy::too_many_arguments)]
    fn audit_ductile_solution(
        &self,
        elements: &[Element],
        initial: Option<&SmallStrainShellSolution<DuctileJ2State<'m, L>>>,
        solution: &SmallStrainShellSolution<DuctileJ2State<'m, L>>,
        inertial_work_j: f64,
        tolerances: DuctileShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<DuctileShellEnergy, SmallStrainTetError> {
        checkpoint(cx)?;
        if solution.increments.len() != 1
            || solution.coordinates.len() != self.mesh.nodes.len()
            || solution.section_states.len() != 3 * elements.len()
            || solution.section_resultants.len() != 3 * elements.len()
            || solution.increments[0].reactions.len() != self.prescribed.len()
            || initial.is_some_and(|old| {
                old.coordinates.len() != self.mesh.nodes.len()
                    || old.section_states.len() != 3 * elements.len()
            })
            || !inertial_work_j.is_finite()
        {
            return Err(invalid(
                "complete single ductile shell endpoints and finite inertia work required",
            ));
        }
        let operator = Operator {
            problem: self,
            elements,
            workspace: RefCell::new(
                PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?,
            ),
            inertia: None,
        };
        let virgin = if initial.is_none() {
            operator.initial_states()
        } else {
            Vec::new()
        };
        let old_states = initial.map_or(virgin.as_slice(), |old| old.section_states.as_slice());
        let length = self.coordinate_length_m;
        let old_q: Vec<_> = (0..6 * self.mesh.nodes.len())
            .map(|d| {
                initial.map_or(0.0, |old| old.coordinates[d / 6][d % 6])
                    * if d % 6 < 3 { 1.0 } else { length }
            })
            .collect();
        let new_q: Vec<_> = solution
            .coordinates
            .iter()
            .flat_map(|q| {
                q.iter()
                    .enumerate()
                    .map(move |(c, x)| x * if c < 3 { 1.0 } else { length })
            })
            .collect();
        let mut energy = DuctileShellEnergy {
            inertial_work_j,
            ..DuctileShellEnergy::default()
        };
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            for q in 0..3 {
                checkpoint(cx)?;
                let index = 3 * e + q;
                let account = self
                    .sections
                    .at(e)
                    .audit_ductile_energy(
                        &old_states[index],
                        &solution.section_states[index],
                        tolerances.section_work_j_m2,
                    )
                    .map_err(|error| {
                        invalid(format!("ductile shell element {e}, point {q}: {error}"))
                    })?;
                if account.resultants != solution.section_resultants[index] {
                    return Err(invalid(
                        "ductile work must retain the exact native shell resultants",
                    ));
                }
                let area = element.geometry.area_m2 / 3.0;
                energy.initial_stored_j += area * account.initial_stored_j_m2;
                energy.final_stored_j += area * account.final_stored_j_m2;
                energy.stored_change_j += area * account.stored_change_j_m2;
                energy.plastic_dissipation_j += area * account.plastic_dissipation_j_m2;
                energy.damage_release_j += area * account.damage_release_j_m2;
                energy.material_numerical_dissipation_j +=
                    area * account.numerical_dissipation_j_m2;
                energy.material_work_j += area * account.material_work_j_m2;
                energy.generalized_work_j += area * account.generalized_work_j_m2;
                energy.transverse_residual_work_j += area * account.transverse_residual_work_j_m2;
                energy.degraded_equivalent_plastic_work_j +=
                    area * account.degraded_equivalent_plastic_work_j_m2;
            }
            let old_local = operator.local_coordinates(element, &old_q)?;
            let new_local = operator.local_coordinates(element, &new_q)?;
            energy.drilling.add(element, &old_local, &new_local);
        }
        let increment = &solution.increments[0];
        let displacement = |d: usize| {
            solution.coordinates[d / 6][d % 6]
                - initial.map_or(0.0, |old| old.coordinates[d / 6][d % 6])
        };
        for &(d, load) in self.nodal_loads {
            energy.dead_load_endpoint_work_j += increment.load_factor * load * displacement(d);
        }
        for (&(d, _), &reaction) in self.prescribed.iter().zip(&increment.reactions) {
            energy.support_endpoint_work_j += reaction * displacement(d);
        }
        let external = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j;
        energy.equilibrium_work_residual_j =
            external - energy.generalized_work_j - energy.drilling.work_j - inertial_work_j;
        energy.balance_residual_j = external + energy.transverse_residual_work_j
            - energy.stored_change_j
            - energy.plastic_dissipation_j
            - energy.damage_release_j
            - energy.material_numerical_dissipation_j
            - (energy.drilling.final_j - energy.drilling.initial_j)
            - energy.drilling.numerical_loss_j
            - inertial_work_j;
        for (what, residual) in [
            ("equilibrium work", energy.equilibrium_work_residual_j),
            (
                "transverse work",
                energy.generalized_work_j + energy.transverse_residual_work_j
                    - energy.material_work_j,
            ),
            (
                "drilling work",
                energy.drilling.work_j
                    - (energy.drilling.final_j - energy.drilling.initial_j)
                    - energy.drilling.numerical_loss_j,
            ),
            ("complete energy balance", energy.balance_residual_j),
        ] {
            gate(what, residual, tolerances.shell_work_j)?;
        }
        if [
            energy.initial_stored_j,
            energy.final_stored_j,
            energy.plastic_dissipation_j,
            energy.damage_release_j,
            energy.material_numerical_dissipation_j,
            energy.degraded_equivalent_plastic_work_j,
            energy.drilling.initial_j,
            energy.drilling.final_j,
            energy.drilling.numerical_loss_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid(
                "ductile shell storage and losses must remain finite and nonnegative",
            ));
        }
        checkpoint(cx)?;
        Ok(energy)
    }
}
