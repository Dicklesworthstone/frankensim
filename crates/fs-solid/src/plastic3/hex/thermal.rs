//! Source-bound thermal history on the native Standard/B-bar Hex8 solid.
//!
//! Each original element family resolves eight independent endpoint laws from
//! pointwise specific enthalpy. Old history is admitted at its OLD temperature,
//! then the same frozen history drives every native equilibrium or inertia trial.
//! The original material owner independently audits each converged point before
//! any displacement, history, temperature, heat or velocity is returned.
//!
//! All work uses reference quadrature. Physical heat, reversible caloric transfer,
//! coefficient/free-energy storage, T times mechanical entropy, and numerical
//! material/kinetic loss remain distinct. This supplies a mechanical proposal to
//! a caloric owner, not a temperature predictor or a separate heat-transport law.
//! Families, reference geometry and formulation must retain their original
//! binding across continuation; no state remapping or finite strain is supplied.

use super::dynamics::{BackwardEulerHexStep, InertialEnergy, ReferenceDensity};
use super::{
    HexAssemblyBudget, SmallStrainHexProblem, SmallStrainHexSolution, checkpoint, tensor_norm,
};
use crate::{Formulation, IncrementSettings, SmallStrainLawField, SmallStrainTetError};
use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::section::thermal::source::{
    ThermalSectionIncrement, ThermalSectionSourceError, ThermalSectionSourceFamily,
};

type History<'a, F> = <<F as ThermalSectionSourceFamily<'a>>::Material as SmallStrainLaw>::State;

/// One immutable family per original cell; temperatures vary independently at
/// its eight native Gauss points. Forces/support values are actual endpoints.
#[derive(Clone, Debug)]
pub struct ThermalHexProblem<'p, F> {
    pub nodes_m: &'p [[f64; 3]],
    pub hexahedra: &'p [[usize; 8]],
    pub families: &'p [F],
    pub formulation: Formulation,
    pub strain_limit: f64,
    pub rotation_limit_rad: f64,
    pub prescribed_m: &'p [(usize, f64)],
    pub nodal_forces_n: &'p [(usize, f64)],
    pub budget: HexAssemblyBudget,
    pub settings: IncrementSettings,
}

/// Borrowed accepted state, in original element-major/Gauss-point order.
/// Refused proposals never modify these arrays.
#[derive(Clone, Debug)]
pub struct ThermalHexHistory<'h, S> {
    pub displacement_m: &'h [f64],
    pub quadrature_states: &'h [S],
    pub specific_enthalpy_j_kg: &'h [f64],
}

/// Absolute local constitutive and whole-body discrete work gates, not a
/// calibration or temporal-accuracy certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalHexTolerances {
    pub material: EnergyTolerances,
    pub energy_j: f64,
}

/// Integrated native point accounts [J]; the same carrier sums the whole solid.
/// Coefficient storage is INCLUDED in stored_change; creep heat is INCLUDED in
/// physical_heat. Neither subaccount is added to the balance a second time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalHexMaterialEnergy {
    pub mechanical_work_j: f64,
    pub stored_change_j: f64,
    pub coefficient_storage_change_j: f64,
    /// Change in T times mechanical entropy, distinct from free energy.
    pub thermoelastic_storage_change_j: f64,
    pub physical_heat_j: f64,
    pub creep_heat_j: f64,
    pub reversible_heat_j: f64,
    pub numerical_loss_j: f64,
}

impl ThermalHexMaterialEnergy {
    fn integrate<S>(increment: &ThermalSectionIncrement<S>, volume: f64) -> Self {
        Self {
            mechanical_work_j: volume * increment.mechanical_work_j_m3,
            stored_change_j: volume * increment.stored_change_j_m3,
            coefficient_storage_change_j: volume * increment.coefficient_storage_change_j_m3,
            thermoelastic_storage_change_j: volume * increment.thermoelastic_storage_change_j_m3,
            physical_heat_j: volume * increment.physical_heat_j_m3,
            creep_heat_j: volume * increment.creep_heat_j_m3,
            reversible_heat_j: volume * increment.reversible_heat_j_m3,
            numerical_loss_j: volume * increment.numerical_loss_j_m3,
        }
    }

    fn accumulate(&mut self, point: Self) -> Result<(), ThermalHexError> {
        for (sum, value) in [
            (&mut self.mechanical_work_j, point.mechanical_work_j),
            (&mut self.stored_change_j, point.stored_change_j),
            (
                &mut self.coefficient_storage_change_j,
                point.coefficient_storage_change_j,
            ),
            (
                &mut self.thermoelastic_storage_change_j,
                point.thermoelastic_storage_change_j,
            ),
            (&mut self.physical_heat_j, point.physical_heat_j),
            (&mut self.creep_heat_j, point.creep_heat_j),
            (&mut self.reversible_heat_j, point.reversible_heat_j),
            (&mut self.numerical_loss_j, point.numerical_loss_j),
        ] {
            *sum = finite(*sum + value)?;
        }
        if point.physical_heat_j < 0.0
            || point.creep_heat_j < 0.0
            || point.creep_heat_j > point.physical_heat_j
            || point.numerical_loss_j < 0.0
        {
            return Err(ThermalHexError::Invalid(
                "invalid native physical, creep or numerical loss",
            ));
        }
        Ok(())
    }
}

/// Independent external/support work and closed material/kinetic balance.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalHexEnergy {
    pub applied_work_j: f64,
    pub support_work_j: f64,
    pub material: ThermalHexMaterialEnergy,
    pub equilibrium_work_residual_j: f64,
    pub balance_residual_j: f64,
}

/// One complete candidate, ready for joint acceptance with its caloric owner.
#[derive(Clone, Debug)]
pub struct ThermalHexStep<S> {
    pub solution: SmallStrainHexSolution<S>,
    pub duration_s: f64,
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub quadrature_temperature_k: Vec<f64>,
    /// Original-volume energy at each point; no duplicate constitutive state.
    pub quadrature_energy: Vec<ThermalHexMaterialEnergy>,
    pub energy: ThermalHexEnergy,
    /// Present together only for the backward-Euler motion entry point.
    pub velocity_mps: Option<Vec<[f64; 3]>>,
    pub inertia: Option<InertialEnergy>,
}

impl<S> ThermalHexStep<S> {
    pub fn history(&self) -> ThermalHexHistory<'_, S> {
        ThermalHexHistory {
            displacement_m: self.solution.displacement_m.as_flattened(),
            quadrature_states: &self.solution.quadrature_states,
            specific_enthalpy_j_kg: &self.specific_enthalpy_j_kg,
        }
    }
}

#[derive(Debug)]
pub enum ThermalHexError {
    Invalid(&'static str),
    Mechanical(SmallStrainTetError),
    Material {
        point: usize,
        source: ThermalSectionSourceError,
    },
    ConstitutiveMismatch {
        point: usize,
    },
    EnergyResidual {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for ThermalHexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal Hex8 proposal refused: {self:?}")
    }
}
impl std::error::Error for ThermalHexError {}
impl From<SmallStrainTetError> for ThermalHexError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Mechanical(error)
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64) -> Result<f64, ThermalHexError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermalHexError::Invalid("nonfinite thermal Hex8 work"))
    }
}
fn gate(residual_j: f64, tolerance_j: f64) -> Result<(), ThermalHexError> {
    finite(residual_j)?;
    if residual_j.abs() > tolerance_j {
        Err(ThermalHexError::EnergyResidual {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

impl<F> ThermalHexProblem<'_, F> {
    /// One physical material interval at prescribed point enthalpies.
    /// Every endpoint law is resolved once with this exact duration.
    pub fn solve_from<'a>(
        &self,
        history: ThermalHexHistory<'_, History<'a, F>>,
        enthalpy_j_kg: &[f64],
        duration_s: f64,
        tolerances: ThermalHexTolerances,
        cx: &Cx<'_>,
    ) -> Result<ThermalHexStep<History<'a, F>>, ThermalHexError>
    where
        F: ThermalSectionSourceFamily<'a>,
    {
        self.advance(history, None, enthalpy_j_kg, duration_s, tolerances, cx)
    }

    /// Same audited thermal interval with original consistent reference mass.
    /// Density comes from the original families, never from the endpoint
    /// temperature or assumed B-bar strain. Kinetic loss is not caloric heat.
    pub fn advance_backward_euler<'a>(
        &self,
        history: ThermalHexHistory<'_, History<'a, F>>,
        velocity_mps: &[[f64; 3]],
        enthalpy_j_kg: &[f64],
        duration_s: f64,
        tolerances: ThermalHexTolerances,
        cx: &Cx<'_>,
    ) -> Result<ThermalHexStep<History<'a, F>>, ThermalHexError>
    where
        F: ThermalSectionSourceFamily<'a>,
    {
        self.advance(
            history,
            Some(velocity_mps),
            enthalpy_j_kg,
            duration_s,
            tolerances,
            cx,
        )
    }

    fn native<'p, L: SmallStrainLaw>(
        &'p self,
        laws: SmallStrainLawField<'p, L>,
    ) -> SmallStrainHexProblem<'p, L> {
        SmallStrainHexProblem {
            nodes_m: self.nodes_m,
            hexahedra: self.hexahedra,
            laws,
            formulation: self.formulation,
            strain_limit: self.strain_limit,
            rotation_limit_rad: self.rotation_limit_rad,
            prescribed_m: self.prescribed_m,
            nodal_forces_n: self.nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings: self.settings,
        }
    }

    fn advance<'a>(
        &self,
        history: ThermalHexHistory<'_, History<'a, F>>,
        velocity_mps: Option<&[[f64; 3]]>,
        enthalpy_j_kg: &[f64],
        duration_s: f64,
        tolerances: ThermalHexTolerances,
        cx: &Cx<'_>,
    ) -> Result<ThermalHexStep<History<'a, F>>, ThermalHexError>
    where
        F: ThermalSectionSourceFamily<'a>,
    {
        checkpoint(cx)?;
        let cells = self.hexahedra.len();
        let points = cells
            .checked_mul(8)
            .ok_or(ThermalHexError::Invalid("point count overflow"))?;
        let dofs = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(ThermalHexError::Invalid("DOF count overflow"))?;
        if cells == 0
            || self.families.len() != cells
            || history.displacement_m.len() != dofs
            || history.quadrature_states.len() != points
            || history.specific_enthalpy_j_kg.len() != points
            || enthalpy_j_kg.len() != points
            || cells > self.budget.maximum_elements
            || points > self.budget.maximum_quadrature_points
            || self.nodes_m.len() > self.budget.maximum_nodes
            || !positive(duration_s)
            || !positive(tolerances.energy_j)
            || !positive(tolerances.material.energy_j_m3)
            || !positive(tolerances.material.yield_pa)
            || history
                .displacement_m
                .iter()
                .chain(history.specific_enthalpy_j_kg)
                .chain(enthalpy_j_kg)
                .any(|x| !x.is_finite())
        {
            return Err(ThermalHexError::Invalid(
                "thermal Hex8 history, duration, budget or work gates",
            ));
        }
        let material_error = |point, source| ThermalHexError::Material { point, source };
        let mut laws = Vec::with_capacity(points);
        let mut temperatures = Vec::with_capacity(points);
        for point in 0..points {
            checkpoint(cx)?;
            let step = self.families[point / 8]
                .resolve_section_step(
                    history.specific_enthalpy_j_kg[point],
                    enthalpy_j_kg[point],
                    duration_s,
                )
                .map_err(|e| material_error(point, e))?;
            if !positive(step.temperature_k) {
                return Err(ThermalHexError::Invalid(
                    "nonpositive resolved absolute temperature",
                ));
            }
            temperatures.push(step.temperature_k);
            laws.push(step.law);
        }
        let fields: Vec<_> = laws.chunks_exact(8).collect();
        let native = self.native(SmallStrainLawField::PerQuadrature(&fields));
        let elements = native.prepare(cx)?;
        native.check(history.displacement_m, history.quadrature_states)?;
        let mut old_strains = Vec::with_capacity(points);
        for (element, geometry) in elements.iter().enumerate() {
            for (local, point) in geometry.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = 8 * element + local;
                let strain =
                    native.strain(geometry, point, history.displacement_m, element, local)?;
                self.families[element]
                    .validate_section_history(
                        history.specific_enthalpy_j_kg[index],
                        &strain,
                        &history.quadrature_states[index],
                        tolerances.material.yield_pa,
                    )
                    .map_err(|e| material_error(index, e))?;
                old_strains.push(strain);
            }
        }
        let (solution, velocity_mps, inertia) = if let Some(velocity) = velocity_mps {
            let densities: Vec<_> = self
                .families
                .iter()
                .map(F::reference_density_kg_m3)
                .collect();
            let BackwardEulerHexStep {
                solution,
                velocity_mps,
                energy,
            } = native.advance_prepared(
                &elements,
                history.displacement_m,
                velocity,
                history.quadrature_states,
                ReferenceDensity::Cell(&densities),
                duration_s,
                cx,
            )?;
            (solution, Some(velocity_mps), Some(energy))
        } else {
            (
                native.solve_prepared(
                    &elements,
                    Some((history.displacement_m, history.quadrature_states)),
                    None,
                    cx,
                )?,
                None,
                None,
            )
        };
        let mut energy = ThermalHexEnergy::default();
        let mut quadrature_energy = Vec::with_capacity(points);
        for point in 0..points {
            checkpoint(cx)?;
            let increment = self.families[point / 8]
                .section_increment(
                    history.specific_enthalpy_j_kg[point],
                    enthalpy_j_kg[point],
                    &old_strains[point],
                    &history.quadrature_states[point],
                    &solution.quadrature_strain[point],
                    tolerances.material,
                    duration_s,
                )
                .map_err(|e| material_error(point, e))?;
            let stress_error = core::array::from_fn(|i| {
                increment.stress_pa[i] - solution.quadrature_stress_pa[point][i]
            });
            if increment.state != solution.quadrature_states[point]
                || !tensor_norm(&stress_error).is_finite()
                || tensor_norm(&stress_error) > tolerances.material.yield_pa
            {
                return Err(ThermalHexError::ConstitutiveMismatch { point });
            }
            let account = ThermalHexMaterialEnergy::integrate(
                &increment,
                solution.quadrature_volume_m3[point],
            );
            energy.material.accumulate(account)?;
            quadrature_energy.push(account);
        }
        let next_u = solution.displacement_m.as_flattened();
        for &(dof, force) in self.nodal_forces_n {
            energy.applied_work_j = finite(
                energy.applied_work_j + force * (next_u[dof] - history.displacement_m[dof]),
            )?;
        }
        for (&(dof, _), &reaction) in self
            .prescribed_m
            .iter()
            .zip(&solution.increments[0].reactions_n)
        {
            energy.support_work_j = finite(
                energy.support_work_j + reaction * (next_u[dof] - history.displacement_m[dof]),
            )?;
        }
        let work = finite(energy.applied_work_j + energy.support_work_j)?;
        let material = energy.material;
        energy.equilibrium_work_residual_j = work - material.mechanical_work_j;
        energy.balance_residual_j = work
            - material.stored_change_j
            - material.thermoelastic_storage_change_j
            - material.physical_heat_j
            - material.reversible_heat_j
            - material.numerical_loss_j;
        if let Some(inertia) = inertia {
            gate(inertia.balance_error_j, tolerances.energy_j)?;
            energy.equilibrium_work_residual_j -= inertia.inertial_work_j;
            energy.balance_residual_j -= inertia.kinetic_change_j + inertia.numerical_loss_j;
        }
        gate(energy.equilibrium_work_residual_j, tolerances.energy_j)?;
        gate(energy.balance_residual_j, tolerances.energy_j)?;
        checkpoint(cx)?;
        Ok(ThermalHexStep {
            solution,
            duration_s,
            specific_enthalpy_j_kg: enthalpy_j_kg.to_vec(),
            quadrature_temperature_k: temperatures,
            quadrature_energy,
            energy,
            velocity_mps,
            inertia,
        })
    }
}
