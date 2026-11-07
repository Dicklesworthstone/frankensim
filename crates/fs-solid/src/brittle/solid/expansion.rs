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
//! Coefficient work, physical fracture loss, surviving-bond viscous loss and
//! endpoint quadrature error have distinct accounts. The timed dynamic owner
//! uses total mean stress in BOTH strain work ports. No partial state is
//! published on refusal.
//!
//! The equilibrium owner is small-strain, fixed-plane, quasistatic fracture
//! with a resolved temperature path. The [`dynamics`] child retains the same
//! material preparation and work recovery during physical motion. Neither owner
//! supplies heat capacity, entropy, finite thermal expansion, automatic thermal
//! substeps or new crack planes. The [`fracture`] child releases explicitly
//! selected completed, traction-free cracks through the original topology owner.
//! Large thermal jumps can cross a softening branch before Newton equilibrates;
//! callers must resolve them with physical temperature/load endpoints.
//! The [`initiation`] child discovers tensile planes from retained elastic
//! thermal prehistory, then hands the old accepted state to this same owner.

/// Retained source-bound motion with joint displacement and thermal eigenstrain.
pub mod dynamics;

/// Completed native fracture release retaining thermal source and strain history.
pub mod fracture;

/// Thermal elastic prehistory, tensile-plane discovery and conservative handoff.
pub mod initiation;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::BrittleSourceError;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::{BrittleEnergy, BrittleState};
use fs_material::tensor::contract;

use super::fracture::ReleasedBrittleFace;
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
    released_faces: Arc<[ReleasedBrittleFace]>,
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
    /// Original native work along the MECHANICAL strain increment, including
    /// physical viscosity when recovered from an explicit timed path.
    pub constitutive_strain_work_j: f64,
    /// Total new fracture loss, including coefficient-induced relaxation.
    pub fracture_dissipation_j: f64,
    pub parameter_fracture_dissipation_j: f64,
    /// Exact physical surviving-bond loss; neither intrinsic fracture nor heat.
    pub viscous_dissipation_j: f64,
    /// Constitutive strain work minus mechanical plus eigenstrain work.
    pub strain_split_error_j: f64,
    /// Parameter plus constitutive work minus storage, fracture and viscous loss.
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
    /// Intrinsic mechanical-strain work/storage/fracture densities, J/m3.
    /// Timed total work and physical viscous loss live in the whole-solid account.
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

/// Checked coefficient-only state awaiting a complete native mechanical solve.
/// No public result exposes this potentially unequilibrated preparation.
#[derive(Debug)]
struct PreparedThermalBrittleTet<'a> {
    states: Vec<BrittleState<'a>>,
    element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    energy: ThermalBrittleTetEnergy,
}

/// Native solver stress semantics, distinct from retained physical endpoint stress.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ThermalBrittleResponse {
    Endpoint,
    Mean,
}

/// Original constitutive work recovered once for equilibrium and inertial schemes.
/// Applied/support work, inertia and the complete balance belong to the caller.
#[derive(Debug)]
struct ThermalBrittleTetRecovery {
    energy: ThermalBrittleTetEnergy,
    element_energy: Vec<BrittleEnergy>,
    cumulative_loss_j: f64,
    /// Actual solver response stress contracted with the TOTAL strain change.
    response_internal_work_j: f64,
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
            let virgin = law.initial_state();
            let stress = law
                .stress_at_rate(&virgin, &[0.0; 6])
                .map_err(|source| BrittleTetError::Material { element, source })?;
            bands.push(band);
            states.push(virgin);
            stresses.push(stress);
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
            released_faces: Arc::from([]),
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
        let prepared =
            self.prepare_materials(next_materials, controls.material_energy_tolerance_j, cx)?;
        let mut laws = Vec::with_capacity(next_materials.len());
        for (element, material) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            if material.band().law().viscosity().is_some() && !prepared.states[element].failed() {
                return Err(BrittleTetError::Invalid(
                    "surviving viscous thermal bonds require a timed dynamic strain path",
                ).into());
            }
            laws.push(material.small_strain_law());
        }
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
        let mut solution = problem
            .solve_from(old_u, &prepared.states, cx)
            .map_err(BrittleTetError::from)?;
        self.require_released_faces(&solution.displacement_m, cx)?;
        let ThermalBrittleTetRecovery {
            mut energy,
            element_energy,
            cumulative_loss_j,
            ..
        } = self.recover_material(
            next_materials,
            &prepared,
            &mut solution,
            ThermalBrittleResponse::Endpoint,
            None,
            controls.material_energy_tolerance_j,
            cx,
        )?;
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
        let external_work = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j;
        energy.equilibrium_work_residual_j = external_work - energy.internal_endpoint_work_j;
        energy.total_balance_error_j =
            ((((external_work + energy.parameter_work_j + energy.eigenstrain_work_j)
                - (energy.stored_after_j - energy.stored_before_j))
                - energy.fracture_dissipation_j)
                - energy.viscous_dissipation_j)
                - energy.endpoint_work_defect_j;
        energy_gate(
            "thermal complete work balance",
            energy.total_balance_error_j,
            controls.material_energy_tolerance_j,
        )?;
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
        checkpoint(cx)?;
        Ok(ThermalBrittleTetStep {
            solid: ThermalBrittleTetSolid {
                geometry: self.geometry,
                materials: next_materials,
                bands: Arc::clone(&self.bands),
                released_faces: Arc::clone(&self.released_faces),
                budget: self.budget,
                solution,
                stored_energy_j: energy.stored_after_j,
                fracture_dissipation_j: cumulative_loss_j,
            },
            element_energy,
            element_transfer: prepared.element_transfer,
            energy,
        })
    }

    /// Prepare original source coefficients while holding MECHANICAL history.
    /// The caller must complete its joint displacement/eigenstrain solve and
    /// work gates before exposing any successor. Preparation never applies the
    /// new thermal strain to a fictitiously restrained displacement state.
    fn prepare_materials<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        material_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<PreparedThermalBrittleTet<'next>, ThermalBrittleTetError> {
        checkpoint(cx)?;
        if next_materials.len() != self.materials.len() {
            return Err(BrittleTetError::MaterialCount {
                expected: self.materials.len(),
                got: next_materials.len(),
            }
            .into());
        }
        if !material_tolerance_j.is_finite() || material_tolerance_j < 0.0 {
            return Err(BrittleTetError::Invalid(
                "finite nonnegative material work tolerance required",
            )
            .into());
        }
        let mut states = Vec::with_capacity(next_materials.len());
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
        }
        energy_gate(
            "thermal coefficient work/storage/fracture",
            (energy.parameter_work_j - (prepared_storage_j - self.stored_energy_j))
                - energy.parameter_fracture_dissipation_j,
            material_tolerance_j,
        )?;
        energy_gate(
            "thermal coefficient cumulative fracture",
            (prepared_loss_j - self.fracture_dissipation_j)
                - energy.parameter_fracture_dissipation_j,
            material_tolerance_j,
        )?;
        if [
            prepared_storage_j,
            prepared_loss_j,
            energy.parameter_fracture_dissipation_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative physical energy required").into(),
            );
        }
        checkpoint(cx)?;
        Ok(PreparedThermalBrittleTet {
            states,
            element_transfer,
            energy,
        })
    }

    /// Recover native endpoint history and exact constitutive work from either
    /// endpoint-force or path-mean-force assembly. Endpoint stress is restored
    /// in the returned native solution; its reactions retain the solver's
    /// original force semantics. A timed response uses the interval secant rate;
    /// the dynamic caller must observe its actual endpoint velocity separately.
    /// No alternate history or fracture law is used.
    /// Public endpoints admit released-face opening separately; only unpublished
    /// native sustained-contact force probes may defer that geometric check.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn recover_material<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        prepared: &PreparedThermalBrittleTet<'next>,
        solution: &mut SmallStrainTetSolution<BrittleState<'next>>,
        response: ThermalBrittleResponse,
        dt_s: Option<f64>,
        material_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetRecovery, ThermalBrittleTetError> {
        checkpoint(cx)?;
        if !material_tolerance_j.is_finite() || material_tolerance_j < 0.0 {
            return Err(BrittleTetError::Invalid(
                "finite nonnegative material work tolerance required",
            )
            .into());
        }
        let mut energy = prepared.energy;
        let mut element_energy = Vec::with_capacity(next_materials.len());
        let mut cumulative_loss_j = 0.0;
        let mut response_internal_work_j = 0.0;
        energy.fracture_dissipation_j = energy.parameter_fracture_dissipation_j;
        for (element, material) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            let total = &solution.element_strain[element];
            let (step, endpoint_stress, mean_stress, strain_work, viscous_loss) =
                if let Some(dt) = dt_s {
                    let path = material
                        .trial_total_over_time_with_mean_stress(total, &prepared.states[element], dt)
                        .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
                    let work = path.endpoint.total_work_j_m3();
                    (
                        path.endpoint.equilibrium,
                        path.endpoint.stress_pa,
                        path.mean_stress_pa,
                        work,
                        path.endpoint.viscous_dissipated_j_m3,
                    )
                } else {
                    let mechanical = material
                        .mechanical_strain(total)
                        .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
                    let path = material
                        .band()
                        .law()
                        .trial_with_mean_stress(&mechanical, &prepared.states[element])
                        .map_err(|source| BrittleTetError::Material { element, source })?;
                    let stress = path.endpoint.stress_pa;
                    let work = path.endpoint.energy.work;
                    (path.endpoint, stress, path.mean_stress_pa, work, 0.0)
                };
            let expected_response = match response {
                ThermalBrittleResponse::Endpoint => endpoint_stress,
                ThermalBrittleResponse::Mean => mean_stress,
            };
            if step.state != solution.element_states[element]
                || expected_response != solution.element_stress[element]
                || step.state.damage() < prepared.states[element].damage()
                || step.state.dissipated_j_m3() < prepared.states[element].dissipated_j_m3()
                || step.energy.dissipated < 0.0
                || !viscous_loss.is_finite()
                || viscous_loss < 0.0
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
            let mechanical_work = contract(&mean_stress, &change);
            let eigenstrain_work = -thermal_change
                * (mean_stress[0] + mean_stress[1] + mean_stress[2]);
            let response_density = contract(&solution.element_stress[element], &change);
            let endpoint_density = match response {
                ThermalBrittleResponse::Endpoint => response_density,
                ThermalBrittleResponse::Mean => contract(&endpoint_stress, &change),
            };
            let volume = self.bands[element].reference_volume_m3();
            energy.stored_after_j =
                volume.mul_add(step.state.stored_energy_j_m3(), energy.stored_after_j);
            energy.constitutive_strain_work_j =
                volume.mul_add(strain_work, energy.constitutive_strain_work_j);
            energy.mechanical_work_j = volume.mul_add(mechanical_work, energy.mechanical_work_j);
            energy.eigenstrain_work_j = volume.mul_add(eigenstrain_work, energy.eigenstrain_work_j);
            energy.internal_endpoint_work_j =
                volume.mul_add(endpoint_density, energy.internal_endpoint_work_j);
            response_internal_work_j = volume.mul_add(response_density, response_internal_work_j);
            energy.fracture_dissipation_j =
                volume.mul_add(step.energy.dissipated, energy.fracture_dissipation_j);
            energy.viscous_dissipation_j =
                volume.mul_add(viscous_loss, energy.viscous_dissipation_j);
            cumulative_loss_j = volume.mul_add(step.state.dissipated_j_m3(), cumulative_loss_j);
            element_energy.push(step.energy);
            solution.element_stress[element] = endpoint_stress;
        }
        energy.strain_split_error_j = energy.constitutive_strain_work_j
            - (energy.mechanical_work_j + energy.eigenstrain_work_j);
        energy.material_balance_error_j = (((energy.parameter_work_j
            + energy.constitutive_strain_work_j)
            - (energy.stored_after_j - energy.stored_before_j))
            - energy.fracture_dissipation_j)
            - energy.viscous_dissipation_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - energy.mechanical_work_j;
        for (account, defect) in [
            (
                "thermal mechanical/eigenstrain work split",
                energy.strain_split_error_j,
            ),
            (
                "thermal material work/storage/fracture/viscosity",
                energy.material_balance_error_j,
            ),
            (
                "thermal cumulative fracture loss",
                (cumulative_loss_j - self.fracture_dissipation_j) - energy.fracture_dissipation_j,
            ),
        ] {
            energy_gate(account, defect, material_tolerance_j)?;
        }
        if [
            energy.stored_after_j,
            energy.fracture_dissipation_j,
            energy.viscous_dissipation_j,
            cumulative_loss_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative physical energy required").into(),
            );
        }
        if !response_internal_work_j.is_finite()
            || !energy.internal_endpoint_work_j.is_finite()
            || !energy.endpoint_work_defect_j.is_finite()
        {
            return Err(
                BrittleTetError::Invalid("finite thermal solver response work required").into(),
            );
        }
        checkpoint(cx)?;
        Ok(ThermalBrittleTetRecovery {
            energy,
            element_energy,
            cumulative_loss_j,
            response_internal_work_j,
        })
    }
}
