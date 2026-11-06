//! Sourced thermal eigenstrain and tensile fracture on the native P1 solver.
//!
//! The original expansion owner integrates alpha(T) from one fixed reference.
//! Total displacement strain and thermal strain enter the SAME Newton trial;
//! only the resulting mechanical strain evolves the original crack-band law.
//! A free body can therefore expand without first acquiring artificial damage
//! from an intermediate fully restrained thermal update.
//!
//! Each endpoint first transfers coefficients at held mechanical history,
//! then solves displacement and eigenstrain together at the new coefficients.
//! Exact native path-mean stress separates mechanical and eigenstrain work.
//! Coefficient work, physical fracture loss and endpoint quadrature error have
//! distinct accounts. No partial state is published on refusal.
//!
//! This is small-strain, fixed-plane, quasistatic fracture with a resolved
//! temperature path. It does not supply heat capacity, entropy, finite thermal
//! expansion, automatic thermal substeps, new crack planes or fragment contact.
//! Large thermal jumps can cross a softening branch before Newton equilibrates;
//! callers must resolve them with physical temperature/load endpoints.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::BrittleSourceError;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::{BrittleEnergy, BrittleState};
use fs_material::tensor::contract;

use super::temperature::BrittleTetTemperatureElementEnergy;
use super::{BrittleTetControls, BrittleTetError, checkpoint, energy_gate};
use crate::brittle::TetCrackBandGeometry;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, SmallStrainTetSolution,
};

/// Immutable original geometry, expansion reference, sources and fracture history.
#[derive(Debug)]
pub struct ThermalBrittleTetSolid<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedThermalBrittle],
    bands: Arc<[TetCrackBandGeometry]>,
    budget: TetAssemblyBudget,
    solution: SmallStrainTetSolution<BrittleState<'a>>,
    stored_energy_j: f64,
    fracture_dissipation_j: f64,
}

/// Complete work split over the original reference mesh, in joules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalBrittleTetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed work of changing coefficients at held mechanical history.
    pub parameter_work_j: f64,
    /// Exact path-mean stress contracted with the TOTAL strain change.
    pub mechanical_work_j: f64,
    /// Minus path-mean stress contracted with the FREE thermal strain change.
    /// This is mechanical eigenstrain work, not an inferred caloric heat input.
    pub eigenstrain_work_j: f64,
    /// Original native work along the MECHANICAL strain increment.
    pub constitutive_strain_work_j: f64,
    /// Total new fracture loss, including coefficient-induced relaxation.
    pub fracture_dissipation_j: f64,
    pub parameter_fracture_dissipation_j: f64,
    /// Constitutive strain work minus mechanical plus eigenstrain work.
    pub strain_split_error_j: f64,
    /// Parameter plus constitutive work minus storage change and fracture loss.
    pub material_balance_error_j: f64,
    /// Endpoint stress contracted with the TOTAL strain change.
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus exact TOTAL-strain mechanical work.
    /// Signed numerical quadrature error; never counted as physical fracture.
    pub endpoint_work_defect_j: f64,
    /// External, parameter and eigenstrain work minus storage, fracture and
    /// the separately identified endpoint quadrature defect.
    pub total_balance_error_j: f64,
}

/// One accepted thermo-mechanical endpoint and original per-element accounts.
#[derive(Debug)]
pub struct ThermalBrittleTetStep<'a> {
    pub solid: ThermalBrittleTetSolid<'a>,
    /// Native mechanical-strain work/storage/fracture densities, J/m3.
    pub element_energy: Vec<BrittleEnergy>,
    /// Separate original coefficient-transfer densities, J/m3.
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub energy: ThermalBrittleTetEnergy,
}

/// Original source/material admission or native solid/work/cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalBrittleTetError {
    Source {
        element: usize,
        source: BrittleSourceError,
    },
    Solid(BrittleTetError),
}

impl core::fmt::Display for ThermalBrittleTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal brittle solid refused: {self:?}")
    }
}
impl std::error::Error for ThermalBrittleTetError {}
impl From<BrittleTetError> for ThermalBrittleTetError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}

impl<'a> ThermalBrittleTetSolid<'a> {
    /// Bind a virgin, stress-free reference state to the original source owners.
    ///
    /// Every free thermal strain must be exactly zero at construction. Existing
    /// prestress/history cannot be silently imported or reset. Each original
    /// plane selects its geometric width under the native assembly budget.
    pub fn try_new(
        geometry: &'a TetLoadGeometry3<'a>,
        materials: &'a [ResolvedThermalBrittle],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermalBrittleTetError> {
        checkpoint(cx)?;
        if materials.len() != geometry.tetrahedra().len() {
            return Err(BrittleTetError::MaterialCount {
                expected: geometry.tetrahedra().len(),
                got: materials.len(),
            }
            .into());
        }
        if geometry.nodes_m().len() > budget.maximum_nodes
            || materials.len() > budget.maximum_tetrahedra
            || geometry.nodes_m().len().checked_mul(3).is_none()
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
        {
            return Err(
                BrittleTetError::Invalid("invalid or exceeded reference assembly budget").into(),
            );
        }
        let mut bands = Vec::with_capacity(materials.len());
        let mut states = Vec::with_capacity(materials.len());
        let mut stresses = Vec::with_capacity(materials.len());
        for (element, material) in materials.iter().enumerate() {
            checkpoint(cx)?;
            if material.free_linear_strain() != 0.0 {
                return Err(BrittleTetError::Invalid(
                    "virgin thermal solid requires zero reference free strain",
                )
                .into());
            }
            let law = material.band().law();
            let band = TetCrackBandGeometry::from_element(geometry, element, law.normal())
                .map_err(BrittleTetError::Geometry)?;
            let required = band.band_width_m();
            let supplied = law.band_width_m();
            if (required - supplied).abs() > 64.0 * f64::EPSILON * required.max(supplied) {
                return Err(BrittleTetError::WidthMismatch {
                    element,
                    required_m: required,
                    material_m: supplied,
                }
                .into());
            }
            let virgin = law
                .trial(&[0.0; 6], &law.initial_state())
                .map_err(|source| BrittleTetError::Material { element, source })?;
            bands.push(band);
            states.push(virgin.state);
            stresses.push(virgin.stress_pa);
        }
        let solution = SmallStrainTetSolution {
            displacement_m: vec![[0.0; 3]; geometry.nodes_m().len()],
            element_strain: vec![[0.0; 6]; materials.len()],
            element_stress: stresses,
            element_states: states,
            increments: Vec::new(),
        };
        checkpoint(cx)?;
        Ok(Self {
            geometry,
            materials,
            bands: bands.into(),
            budget,
            solution,
            stored_energy_j: 0.0,
            fracture_dissipation_j: 0.0,
        })
    }

    #[must_use]
    pub const fn geometry(&self) -> &'a TetLoadGeometry3<'a> {
        self.geometry
    }

    #[must_use]
    pub const fn materials(&self) -> &'a [ResolvedThermalBrittle] {
        self.materials
    }

    #[must_use]
    pub fn band_geometry(&self) -> &[TetCrackBandGeometry] {
        &self.bands
    }

    /// Native FE strain is TOTAL; each native fracture state's strain is
    /// MECHANICAL (total minus the integrated isotropic thermal eigenstrain).
    #[must_use]
    pub const fn solution(&self) -> &SmallStrainTetSolution<BrittleState<'a>> {
        &self.solution
    }

    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }

    #[must_use]
    pub const fn fracture_dissipation_j(&self) -> f64 {
        self.fracture_dissipation_j
    }

    /// Transfer source coefficients and equilibrate thermal strain and motion.
    ///
    /// Targets are absolute prescribed displacements and dead nodal forces.
    /// Successors retain the original card, claims, reference, nonthermal
    /// conditions, plane, width and bounds. Every Newton proposal starts from
    /// the same transferred MECHANICAL history, so rejected proposals cannot
    /// predamage an otherwise freely expanding body.
    ///
    /// The material tolerance gates transfer, work split and full balance;
    /// the endpoint-defect tolerance applies only to TOTAL-strain mechanical
    /// work. Every refusal leaves the complete original state retryable.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetStep<'next>, ThermalBrittleTetError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if next_materials.len() != self.materials.len() {
            return Err(BrittleTetError::MaterialCount {
                expected: self.materials.len(),
                got: next_materials.len(),
            }
            .into());
        }
        if [
            controls.material_energy_tolerance_j,
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative work tolerances required").into(),
            );
        }
        let old_u = self.solution.displacement_m.as_flattened();
        if prescribed_m.len() > old_u.len() || nodal_forces_n.len() > old_u.len() {
            return Err(BrittleTetError::Invalid(
                "endpoint load entries exceed reference DOF count",
            )
            .into());
        }
        let mut states = Vec::with_capacity(next_materials.len());
        let mut laws = Vec::with_capacity(next_materials.len());
        let mut element_transfer = Vec::with_capacity(next_materials.len());
        let mut energy = ThermalBrittleTetEnergy {
            stored_before_j: self.stored_energy_j,
            ..ThermalBrittleTetEnergy::default()
        };
        let mut prepared_storage_j = 0.0;
        let mut prepared_loss_j = 0.0;
        for (element, (old_source, next_source)) in
            self.materials.iter().zip(next_materials).enumerate()
        {
            checkpoint(cx)?;
            let old_state = &self.solution.element_states[element];
            let step = old_source
                .transfer_coefficients(
                    old_state,
                    next_source,
                    next_source.band().law().settings().energy_tolerance_j_m3,
                )
                .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
            if !step.state.is_valid()
                || step.state.strain() != old_state.strain()
                || step.state.damage() < old_state.damage()
                || step.state.dissipated_j_m3() < old_state.dissipated_j_m3()
                || step.damage_dissipated_j_m3 < 0.0
            {
                return Err(BrittleTetError::Invalid(
                    "coefficient transfer changed strain or healed history",
                )
                .into());
            }
            let volume = self.bands[element].reference_volume_m3();
            prepared_storage_j =
                volume.mul_add(step.state.stored_energy_j_m3(), prepared_storage_j);
            prepared_loss_j = volume.mul_add(step.state.dissipated_j_m3(), prepared_loss_j);
            energy.parameter_work_j =
                volume.mul_add(step.parameter_work_j_m3, energy.parameter_work_j);
            energy.parameter_fracture_dissipation_j = volume.mul_add(
                step.damage_dissipated_j_m3,
                energy.parameter_fracture_dissipation_j,
            );
            element_transfer.push(BrittleTetTemperatureElementEnergy {
                parameter_work_j_m3: step.parameter_work_j_m3,
                damage_dissipated_j_m3: step.damage_dissipated_j_m3,
                stored_change_j_m3: step.stored_change_j_m3,
                defect_j_m3: step.defect_j_m3,
            });
            states.push(step.state);
            laws.push(next_source.small_strain_law());
        }
        energy_gate(
            "thermal coefficient work/storage/fracture",
            (energy.parameter_work_j - (prepared_storage_j - self.stored_energy_j))
                - energy.parameter_fracture_dissipation_j,
            controls.material_energy_tolerance_j,
        )?;
        energy_gate(
            "thermal coefficient cumulative fracture",
            (prepared_loss_j - self.fracture_dissipation_j)
                - energy.parameter_fracture_dissipation_j,
            controls.material_energy_tolerance_j,
        )?;
        let problem = SmallStrainTetProblem {
            nodes_m: self.geometry.nodes_m(),
            tetrahedra: self.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings,
        };
        let solution = problem
            .solve_from(old_u, &states, cx)
            .map_err(BrittleTetError::from)?;
        let mut element_energy = Vec::with_capacity(next_materials.len());
        let mut cumulative_loss = 0.0;
        energy.fracture_dissipation_j = energy.parameter_fracture_dissipation_j;
        for (element, material) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            let mechanical = material
                .mechanical_strain(&solution.element_strain[element])
                .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
            let path = material
                .band()
                .law()
                .trial_with_mean_stress(&mechanical, &states[element])
                .map_err(|source| BrittleTetError::Material { element, source })?;
            let step = path.endpoint;
            if step.state != solution.element_states[element]
                || step.stress_pa != solution.element_stress[element]
                || step.state.damage() < states[element].damage()
                || step.state.dissipated_j_m3() < states[element].dissipated_j_m3()
                || step.energy.dissipated < 0.0
            {
                return Err(BrittleTetError::Invalid(
                    "native thermal history/stress or irreversibility mismatch",
                )
                .into());
            }
            let change = core::array::from_fn(|i| {
                solution.element_strain[element][i] - self.solution.element_strain[element][i]
            });
            let thermal_change =
                material.free_linear_strain() - self.materials[element].free_linear_strain();
            let mechanical_work = contract(&path.mean_stress_pa, &change);
            let eigenstrain_work = -thermal_change
                * (path.mean_stress_pa[0] + path.mean_stress_pa[1] + path.mean_stress_pa[2]);
            let volume = self.bands[element].reference_volume_m3();
            energy.stored_after_j =
                volume.mul_add(step.state.stored_energy_j_m3(), energy.stored_after_j);
            energy.constitutive_strain_work_j =
                volume.mul_add(step.energy.work, energy.constitutive_strain_work_j);
            energy.mechanical_work_j = volume.mul_add(mechanical_work, energy.mechanical_work_j);
            energy.eigenstrain_work_j = volume.mul_add(eigenstrain_work, energy.eigenstrain_work_j);
            energy.internal_endpoint_work_j = volume.mul_add(
                contract(&step.stress_pa, &change),
                energy.internal_endpoint_work_j,
            );
            energy.fracture_dissipation_j =
                volume.mul_add(step.energy.dissipated, energy.fracture_dissipation_j);
            cumulative_loss = volume.mul_add(step.state.dissipated_j_m3(), cumulative_loss);
            element_energy.push(step.energy);
        }
        let new_u = solution.displacement_m.as_flattened();
        for &(dof, force) in nodal_forces_n {
            checkpoint(cx)?;
            energy.dead_load_endpoint_work_j =
                force.mul_add(new_u[dof] - old_u[dof], energy.dead_load_endpoint_work_j);
        }
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&solution.increments[0].reactions_n) {
            checkpoint(cx)?;
            energy.support_endpoint_work_j =
                reaction.mul_add(new_u[dof] - old_u[dof], energy.support_endpoint_work_j);
        }
        energy.strain_split_error_j = energy.constitutive_strain_work_j
            - (energy.mechanical_work_j + energy.eigenstrain_work_j);
        energy.material_balance_error_j = ((energy.parameter_work_j
            + energy.constitutive_strain_work_j)
            - (energy.stored_after_j - energy.stored_before_j))
            - energy.fracture_dissipation_j;
        let external_work = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j;
        energy.equilibrium_work_residual_j = external_work - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - energy.mechanical_work_j;
        energy.total_balance_error_j =
            (((external_work + energy.parameter_work_j + energy.eigenstrain_work_j)
                - (energy.stored_after_j - energy.stored_before_j))
                - energy.fracture_dissipation_j)
                - energy.endpoint_work_defect_j;
        for (account, defect) in [
            (
                "thermal mechanical/eigenstrain work split",
                energy.strain_split_error_j,
            ),
            (
                "thermal material work/storage/fracture",
                energy.material_balance_error_j,
            ),
            (
                "thermal cumulative fracture loss",
                (cumulative_loss - self.fracture_dissipation_j) - energy.fracture_dissipation_j,
            ),
            (
                "thermal complete work balance",
                energy.total_balance_error_j,
            ),
        ] {
            energy_gate(account, defect, controls.material_energy_tolerance_j)?;
        }
        energy_gate(
            "thermal equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "thermal mechanical endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if [
            prepared_storage_j,
            prepared_loss_j,
            energy.stored_after_j,
            energy.parameter_fracture_dissipation_j,
            energy.fracture_dissipation_j,
            cumulative_loss,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative physical energy required").into(),
            );
        }
        checkpoint(cx)?;
        Ok(ThermalBrittleTetStep {
            solid: ThermalBrittleTetSolid {
                geometry: self.geometry,
                materials: next_materials,
                bands: Arc::clone(&self.bands),
                budget: self.budget,
                solution,
                stored_energy_j: energy.stored_after_j,
                fracture_dissipation_j: cumulative_loss,
            },
            element_energy,
            element_transfer,
            energy,
        })
    }
}
