//! One-step native material work and numerical drilling energy on CST/DKT shells.
//!
//! Through-thickness energy comes from the original material laws and is
//! integrated on original facet area. Physical plastic loss is never inferred
//! from a global energy deficit. Numerical material loss, drilling storage and
//! drilling integration loss remain separate. Transverse residual work is
//! signed and explicitly included in the balance, rather than treated as heat.

use super::{
    Operator, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, dot, invalid,
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
    pub equilibrium_work_residual_j: f64,
    /// External plus transverse work minus material and drilling storage/losses.
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

fn gate(what: &str, value: f64, tolerance: f64) -> Result<(), SmallStrainTetError> {
    if !value.is_finite() || value.abs() > tolerance {
        return Err(invalid(format!(
            "shell {what}: {value} J exceeds {tolerance} J"
        )));
    }
    Ok(())
}

impl<L: NativeSectionEnergyLaw> SmallStrainShellProblem<'_, L> {
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
        if self.load_path.len() != 1
            || [
                tolerances.material.energy_j_m3,
                tolerances.material.yield_pa,
                tolerances.section_work_j_m2,
                tolerances.shell_work_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "one shell physical interval and positive finite energy allowances required",
            ));
        }
        let elements = self.prepare(cx)?;
        let solution = self.solve_prepared(
            &elements,
            initial.map(|old| (old.coordinates.as_slice(), old.section_states.as_slice())),
            cx,
        )?;
        let operator = Operator {
            problem: self,
            elements: &elements,
            workspace: RefCell::new(
                PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?,
            ),
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
        let mut energy = PlasticShellEnergy::default();
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
            for row in &element.drilling {
                let before = dot(row, &old_local);
                let after = dot(row, &new_local);
                let delta = after - before;
                let weight = element.drilling_weight;
                energy.drilling_initial_j += 0.5 * weight * before * before;
                energy.drilling_final_j += 0.5 * weight * after * after;
                energy.drilling_work_j += weight * after * delta;
                energy.drilling_numerical_loss_j += 0.5 * weight * delta * delta;
            }
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
            external - energy.generalized_work_j - energy.drilling_work_j;
        energy.balance_residual_j = external + energy.transverse_residual_work_j
            - energy.stored_change_j()
            - energy.physical_dissipation_j
            - energy.material_numerical_dissipation_j
            - (energy.drilling_final_j - energy.drilling_initial_j)
            - energy.drilling_numerical_loss_j;
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
        Ok(PlasticShellEnergyStep { solution, energy })
    }
}
