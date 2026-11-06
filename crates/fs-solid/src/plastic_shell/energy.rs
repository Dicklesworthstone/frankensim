//! One-step native material work and numerical drilling energy on CST/DKT shells.
//!
//! Through-thickness energy comes from the original material laws and is
//! integrated on original facet area. Physical plastic loss is never inferred
//! from a global energy deficit. Numerical material loss, drilling storage and
//! drilling integration loss remain separate. Transverse residual work is
//! signed and explicitly included in the balance, rather than treated as heat.

use super::drilling::ShellDrillingEnergy;
use super::dynamics::{PlasticShellReferenceInertia, ShellInertialEnergy};
use super::{
    Element, Operator, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, invalid,
};
use crate::plastic3::SmallStrainTetError;
use crate::plastic3::increment::IncrementOperator;
use fs_exec::Cx;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::PlaneStressWorkspace;
use fs_material::plastic::section::energy::{NativeSectionEnergyLaw, SectionStoredEnergy};
use std::cell::RefCell;

/// Independent SI material, section and whole-shell allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellEnergyTolerances {
    pub material: EnergyTolerances,
    pub section_work_j_m2: f64,
    pub shell_work_j: f64,
}

/// Physical stored material energy over the whole reference shell, J.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellStoredEnergy {
    pub elastic_j: f64,
    pub hardening_j: f64,
    pub kinematic_hardening_j: f64,
}
impl ShellStoredEnergy {
    pub fn total_j(self) -> f64 {
        self.elastic_j + self.hardening_j + self.kinematic_hardening_j
    }
    fn add(&mut self, value: SectionStoredEnergy, area: f64) {
        self.elastic_j += area * value.elastic_j_m2;
        self.hardening_j += area * value.hardening_j_m2;
        self.kinematic_hardening_j += area * value.kinematic_hardening_j_m2;
    }
}

/// Actual mechanical endpoint work, storage and separately identified losses.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlasticShellEnergy {
    pub initial: ShellStoredEnergy,
    pub final_state: ShellStoredEnergy,
    /// Native yield/overstress/recovery heat under the selected energetic law.
    pub physical_dissipation_j: f64,
    pub material_numerical_dissipation_j: f64,
    pub material_work_j: f64,
    pub generalized_work_j: f64,
    pub transverse_residual_work_j: f64,
    /// Numerical elastic penalty storage, excluded from physical material storage.
    pub drilling_initial_j: f64,
    pub drilling_final_j: f64,
    pub drilling_work_j: f64,
    pub drilling_numerical_loss_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Signed endpoint inertia work, including any separately reported drilling
    /// rate gauge. Zero in a static solve; never physical material heat.
    pub inertial_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    /// External plus transverse work minus material/drilling storage and losses
    /// and the signed endpoint inertia work.
    pub balance_residual_j: f64,
}
impl PlasticShellEnergy {
    pub fn stored_change_j(&self) -> f64 {
        (self.final_state.elastic_j - self.initial.elastic_j)
            + (self.final_state.hardening_j - self.initial.hardening_j)
            + (self.final_state.kinematic_hardening_j - self.initial.kinematic_hardening_j)
    }
}

/// A complete solved endpoint published only after its work gates pass.
#[derive(Debug)]
pub struct PlasticShellEnergyStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub energy: PlasticShellEnergy,
}

/// Complete one-step motion, published only after material and inertial work
/// close. Physical director kinetics and numerical drilling-gauge kinetics
/// remain separate in `inertia`; neither backward-Euler loss heats the material.
#[derive(Debug)]
pub struct PlasticShellDynamicEnergyStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    /// Physical translation rates [m/s] followed by axial rotation rates [rad/s].
    pub velocity: Vec<[f64; 6]>,
    pub inertia: ShellInertialEnergy,
    pub energy: PlasticShellEnergy,
}
impl<S> PlasticShellDynamicEnergyStep<S> {
    /// External plus transverse work minus all physical and numerical storage
    /// changes and losses, including the independently computed kinetic work.
    pub fn balance_residual_j(&self) -> f64 {
        self.energy.balance_residual_j + self.inertia.total.balance_error_j
    }
}

fn gate(what: &str, value: f64, tolerance: f64) -> Result<(), SmallStrainTetError> {
    if !value.is_finite() || value.abs() > tolerance {
        return Err(invalid(format!(
            "shell {what}: {value} J exceeds {tolerance} J"
        )));
    }
    Ok(())
}

fn validate_tolerances(tolerances: ShellEnergyTolerances) -> Result<(), SmallStrainTetError> {
    if [
        tolerances.material.energy_j_m3,
        tolerances.material.yield_pa,
        tolerances.section_work_j_m2,
        tolerances.shell_work_j,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid("positive finite shell energy allowances required"));
    }
    Ok(())
}

impl<L: NativeSectionEnergyLaw> SmallStrainShellProblem<'_, L> {
    /// Solve one endpoint momentum equation and audit its exact native material
    /// return. The reference inertia remains bound to its original inventory;
    /// accepted coordinates, velocities and histories are never mutated.
    /// `load_path` must be `[1.0]`, and timed native laws must have exactly the
    /// same original physical duration as `dt_s`. This does not retime a law or
    /// transport material heat. Refine `dt_s` to control backward-Euler damping.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler_with_energy(
        &self,
        initial: Option<&SmallStrainShellSolution<L::State>>,
        velocity: &[[f64; 6]],
        reference: &PlasticShellReferenceInertia<'_>,
        dt_s: f64,
        tolerances: ShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellDynamicEnergyStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        validate_tolerances(tolerances)?;
        if self.load_path != [1.0] || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid(
                "one endpoint shell load and positive physical duration required",
            ));
        }
        let elements = self.prepare(cx)?;
        for e in 0..elements.len() {
            checkpoint(cx)?;
            for layer in self.sections.at(e).layers() {
                if layer
                    .law
                    .physical_duration_s()
                    .is_some_and(|duration| duration != dt_s)
                {
                    return Err(invalid(
                        "shell inertia and timed material durations must agree",
                    ));
                }
            }
        }
        let step = self
            .advance_backward_euler_prepared(&elements, initial, velocity, reference, dt_s, cx)?;
        let energy = self.audit_solution(
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
        Ok(PlasticShellDynamicEnergyStep {
            solution: step.solution,
            velocity: step.velocity,
            inertia: step.inertia,
            energy,
        })
    }

    /// Solve and audit one physical endpoint from virgin state or a complete
    /// accepted checkpoint. Exactly one `load_path` entry is required. The same
    /// law coefficients, reference geometry and energy interpretation must own
    /// both endpoints; changing stiffness is not an unaccounted energy source.
    ///
    /// The initial checkpoint remains unchanged even if a final energy gate
    /// refuses. Timed constitutive laws retain their original physical interval
    /// through every trial and the single audited return. No thermal transport
    /// or temperature evolution is performed by this mechanical account.
    pub fn solve_with_energy(
        &self,
        initial: Option<&SmallStrainShellSolution<L::State>>,
        tolerances: ShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellEnergyStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        validate_tolerances(tolerances)?;
        if self.load_path.len() != 1 {
            return Err(invalid("one shell physical interval required"));
        }
        let elements = self.prepare(cx)?;
        let solution = self.solve_prepared(
            &elements,
            initial.map(|old| (old.coordinates.as_slice(), old.section_states.as_slice())),
            None,
            cx,
        )?;
        let energy = self.audit_solution(&elements, initial, &solution, 0.0, tolerances, cx)?;
        Ok(PlasticShellEnergyStep { solution, energy })
    }

    /// Audit the same accepted native return for static and inertial endpoints.
    /// The motion owner supplies separately gated physical and gauge inertia;
    /// their endpoint work is removed here, never inferred as material heat.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn audit_solution(
        &self,
        elements: &[Element],
        initial: Option<&SmallStrainShellSolution<L::State>>,
        solution: &SmallStrainShellSolution<L::State>,
        inertial_work_j: f64,
        tolerances: ShellEnergyTolerances,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellEnergy, SmallStrainTetError> {
        checkpoint(cx)?;
        validate_tolerances(tolerances)?;
        if self.load_path.len() != 1
            || solution.increments.len() != 1
            || solution.coordinates.len() != self.mesh.nodes.len()
            || solution.section_states.len() != 3 * elements.len()
            || initial.is_some_and(|old| {
                old.coordinates.len() != self.mesh.nodes.len()
                    || old.section_states.len() != 3 * elements.len()
            })
            || !inertial_work_j.is_finite()
        {
            return Err(invalid(
                "complete single shell endpoints and finite inertia work required",
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
        let mut energy = PlasticShellEnergy {
            inertial_work_j,
            ..PlasticShellEnergy::default()
        };
        let mut drilling = ShellDrillingEnergy::default();
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            for q in 0..3 {
                checkpoint(cx)?;
                let account = self
                    .sections
                    .at(e)
                    .audit_energy(
                        &old_states[3 * e + q],
                        &solution.section_states[3 * e + q],
                        tolerances.material,
                        tolerances.section_work_j_m2,
                    )
                    .map_err(|error| {
                        invalid(format!("shell energy element {e}, point {q}: {error}"))
                    })?;
                let area = element.geometry.area_m2 / 3.0;
                energy.initial.add(account.initial, area);
                energy.final_state.add(account.final_state, area);
                energy.physical_dissipation_j += area * account.physical_dissipation_j_m2;
                energy.material_numerical_dissipation_j +=
                    area * account.numerical_dissipation_j_m2;
                energy.material_work_j += area * account.material_work_j_m2;
                energy.generalized_work_j += area * account.generalized_work_j_m2;
                energy.transverse_residual_work_j += area * account.transverse_residual_work_j_m2;
            }
            let old_local = operator.local_coordinates(element, &old_q)?;
            let new_local = operator.local_coordinates(element, &new_q)?;
            drilling.add(element, &old_local, &new_local);
        }
        energy.drilling_initial_j = drilling.initial_j;
        energy.drilling_final_j = drilling.final_j;
        energy.drilling_work_j = drilling.work_j;
        energy.drilling_numerical_loss_j = drilling.numerical_loss_j;
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
            external - energy.generalized_work_j - energy.drilling_work_j - energy.inertial_work_j;
        energy.balance_residual_j = external + energy.transverse_residual_work_j
            - energy.stored_change_j()
            - energy.physical_dissipation_j
            - energy.material_numerical_dissipation_j
            - (energy.drilling_final_j - energy.drilling_initial_j)
            - energy.drilling_numerical_loss_j
            - energy.inertial_work_j;
        for (what, residual) in [
            ("equilibrium work", energy.equilibrium_work_residual_j),
            (
                "transverse work",
                energy.generalized_work_j + energy.transverse_residual_work_j
                    - energy.material_work_j,
            ),
            (
                "drilling work",
                energy.drilling_work_j
                    - (energy.drilling_final_j - energy.drilling_initial_j)
                    - energy.drilling_numerical_loss_j,
            ),
            ("complete energy balance", energy.balance_residual_j),
        ] {
            gate(what, residual, tolerances.shell_work_j)?;
        }
        if [
            energy.initial.total_j(),
            energy.final_state.total_j(),
            energy.physical_dissipation_j,
            energy.material_numerical_dissipation_j,
            energy.drilling_initial_j,
            energy.drilling_final_j,
            energy.drilling_numerical_loss_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid(
                "shell energy storage and losses must be finite and nonnegative",
            ));
        }
        checkpoint(cx)?;
        Ok(energy)
    }
}
