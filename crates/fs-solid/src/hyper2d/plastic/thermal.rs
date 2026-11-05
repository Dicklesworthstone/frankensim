//! Retained thermal sheet forming through the original material and sheet owners.
//!
//! Endpoint temperatures and the actual physical interval freeze the complete
//! source-bound family for every local/global trial. Accepted temperature,
//! deformation, history and time advance together after equilibrium and work
//! checks. The same P1/Q1 geometry, condensed tangent and increment solver serve
//! the fixed-temperature sheet; no parallel constitutive or Newton path exists.
//!
//! Native two-temperature accounts distinguish free/internal storage, physical
//! plastic/creep loss, signed reversible caloric exchange and endpoint work
//! defects. Only the native caloric transfer belongs in a separate heat solve.
//! Temperature is supplied by the caller; this module invents no caloric law,
//! heat solver, source data, shell bending, contact or forming-limit model.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::FiniteJ2State;
use fs_material::plastic::finite::plane_stress::FinitePlaneStressResponse;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteThermalJ2Material, FrozenThermalFiniteJ2,
};

use super::{
    Element, FinitePlasticSheetError, SheetAssemblyBudget, SheetOperator, SheetPointResponse,
    checkpoint, energy_gate, invalid, prepare_elements, validate_step_inputs,
};
use crate::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use crate::mesh2::Mesh2;
use crate::plastic3::increment::{self, IncrementControls};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

// These accounts are dimension-independent densities and whole-body energies.
// Reuse their established definitions instead of maintaining a second heat split.
pub use super::FinitePlasticSheetEvaluation as FiniteThermoplasticSheetEvaluation;
pub use crate::hyper3::plastic::thermal::{
    FiniteThermoplasticElementEnergy,
    FiniteThermoplasticTetEnergy as FiniteThermoplasticSheetEnergy,
};

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Native geometry/solver and thermal point errors retain their original owners.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteThermoplasticSheetError {
    Solid(FinitePlasticSheetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2ThermalError,
    },
}
impl From<SmallStrainTetError> for FiniteThermoplasticSheetError {
    fn from(value: SmallStrainTetError) -> Self {
        Self::Solid(value.into())
    }
}
impl From<FinitePlasticSheetError> for FiniteThermoplasticSheetError {
    fn from(value: FinitePlasticSheetError) -> Self {
        Self::Solid(value)
    }
}
impl core::fmt::Display for FiniteThermoplasticSheetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material {
                element,
                point,
                source,
            } => {
                write!(
                    f,
                    "finite thermal sheet element {element} point {point}: {source}"
                )
            }
        }
    }
}
impl std::error::Error for FiniteThermoplasticSheetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Material { source, .. } => Some(source),
        }
    }
}

/// Immutable complete accepted sheet; original source families remain borrowed.
#[derive(Debug)]
pub struct FiniteThermoplasticSheet<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    mesh: &'a Mesh2,
    families: FinitePlasticLawField<'a, L>,
    thickness_m: &'a [f64],
    elements: Arc<Vec<Element>>,
    local: PlaneStressSettings,
    budget: SheetAssemblyBudget,
    displacement_m: Vec<[f64; 2]>,
    states: Vec<FiniteJ2State>,
    deformation: Vec<[f64; 9]>,
    piola_pa: Vec<[f64; 9]>,
    temperatures_k: Vec<f64>,
    free_energy_j: f64,
    internal_energy_j: f64,
    physical_time_s: f64,
}

/// Complete physical endpoint; point accounts use original reference volume.
#[derive(Debug)]
pub struct FiniteThermoplasticSheetStep<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    pub sheet: FiniteThermoplasticSheet<'a, L>,
    pub equilibrium: IncrementRecord,
    /// Original element order, then the element's native quadrature order.
    pub quadrature_energy: Vec<FiniteThermoplasticElementEnergy>,
    pub energy: FiniteThermoplasticSheetEnergy,
}

struct FrozenSheet<'s, 'a, L: FiniteThermalJ2Material> {
    sheet: &'s FiniteThermoplasticSheet<'a, L>,
    laws: Vec<FrozenThermalFiniteJ2>,
}
impl<L: FiniteThermalJ2Material> SheetPointResponse for FrozenSheet<'_, '_, L> {
    type Error = FiniteJ2ThermalError;
    fn node_count(&self) -> usize {
        self.sheet.mesh.nodes.len()
    }
    fn elements(&self) -> &[Element] {
        &self.sheet.elements
    }
    fn states(&self) -> &[FiniteJ2State] {
        &self.sheet.states
    }
    fn response(
        &self,
        e: usize,
        q: usize,
        u: &[f64],
        states: &[FiniteJ2State],
    ) -> Result<FinitePlaneStressResponse, Self::Error> {
        let element = &self.sheet.elements[e];
        let index = element.offset + q;
        self.laws[e].plane_stress_response(
            &element.in_plane(&element.points[q], u),
            &states[index],
            [2, 5, 8].map(|i| self.sheet.deformation[index][i]),
            self.sheet.local,
        )
    }
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticSheet<'a, L> {
    /// Admit virgin identity geometry at each original family's reference T.
    /// This inspects thermodynamics without inventing an initialization interval.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        mesh: &'a Mesh2,
        families: FinitePlasticLawField<'a, L>,
        thickness_m: &'a [f64],
        initial_temperatures_k: &[f64],
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticSheetError> {
        let elements = prepare_elements(mesh, thickness_m, budget, cx)?;
        if let FinitePlasticLawField::PerElement(values) = families
            && values.len() != mesh.elems.len()
        {
            return Err(SmallStrainTetError::LawCountMismatch {
                expected: mesh.elems.len(),
                got: values.len(),
            }
            .into());
        }
        if initial_temperatures_k.len() != mesh.elems.len()
            || !local.relative_tolerance.is_finite()
            || local.relative_tolerance <= 0.0
            || local.relative_tolerance >= 1.0
            || local.max_iterations == 0
            || local.max_backtracks == 0
        {
            return Err(invalid("thermal sheet needs one initial temperature per element and valid plane-stress controls").into());
        }
        let mut states = Vec::new();
        let mut deformation = Vec::new();
        let mut piola_pa = Vec::new();
        let (mut free_energy_j, mut internal_energy_j) = (0.0, 0.0);
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let family = match families {
                FinitePlasticLawField::Uniform(law) => law,
                FinitePlasticLawField::PerElement(laws) => &laws[e],
            }
            .thermal_finite_j2();
            if initial_temperatures_k[e] != family.parameters().reference_temperature_k {
                return Err(invalid(format!("thermal sheet element {e} must start at its original family's reference temperature")).into());
            }
            let map = |source| FiniteThermoplasticSheetError::Material {
                element: e,
                point: 0,
                source,
            };
            let law = family
                .at_temperature_for_inspection(initial_temperatures_k[e])
                .map_err(map)?;
            let state = family.initial_state();
            law.validate_committed(&IDENTITY, &state, 0.0)
                .map_err(map)?;
            let thermodynamics = law.thermodynamics(&IDENTITY, &state).map_err(map)?;
            for point in &element.points {
                states.push(state.clone());
                deformation.push(IDENTITY);
                piola_pa.push(thermodynamics.first_piola_pa);
                free_energy_j += point.volume_m3 * thermodynamics.free_energy_j_m3;
                internal_energy_j += point.volume_m3 * thermodynamics.internal_energy_j_m3;
            }
        }
        if !free_energy_j.is_finite() || !internal_energy_j.is_finite() {
            return Err(invalid("thermal sheet initial energy is not representable").into());
        }
        checkpoint(cx)?;
        Ok(Self {
            mesh,
            families,
            thickness_m,
            elements,
            local,
            budget,
            displacement_m: vec![[0.0; 2]; mesh.nodes.len()],
            states,
            deformation,
            piola_pa,
            temperatures_k: initial_temperatures_k.to_vec(),
            free_energy_j,
            internal_energy_j,
            physical_time_s: 0.0,
        })
    }

    pub fn families(&self) -> FinitePlasticLawField<'a, L> {
        self.families
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    pub fn quadrature_states(&self) -> &[FiniteJ2State] {
        &self.states
    }
    pub fn quadrature_deformation(&self) -> &[[f64; 9]] {
        &self.deformation
    }
    pub fn quadrature_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.piola_pa
    }
    pub fn free_energy_j(&self) -> f64 {
        self.free_energy_j
    }
    /// Mechanical internal energy; the caller's caloric chart is separate.
    pub fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub fn physical_time_s(&self) -> f64 {
        self.physical_time_s
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements
            .get(element)
            .map(|element| element.points.iter().map(|p| p.volume_m3).sum())
    }
    pub fn quadrature_volume_m3(&self, element: usize, point: usize) -> Option<f64> {
        self.elements
            .get(element)?
            .points
            .get(point)
            .map(|point| point.volume_m3)
    }
    pub fn quadrature_thickness_m(&self) -> Vec<f64> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(e, element)| {
                (0..element.points.len())
                    .map(move |q| self.thickness_m[e] * self.deformation[element.offset + q][8])
            })
            .collect()
    }

    fn family(&self, element: usize) -> &AffineThermalFiniteJ2 {
        match self.families {
            FinitePlasticLawField::Uniform(law) => law,
            FinitePlasticLawField::PerElement(laws) => &laws[element],
        }
        .thermal_finite_j2()
    }

    fn freeze(
        &self,
        temperatures_k: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<FrozenSheet<'_, 'a, L>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        if temperatures_k.len() != self.elements.len() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("thermal sheet needs one endpoint temperature per element and a positive finite duration").into());
        }
        let mut laws = Vec::with_capacity(self.elements.len());
        for (e, &temperature) in temperatures_k.iter().enumerate() {
            checkpoint(cx)?;
            laws.push(
                self.family(e)
                    .at_temperature_for_step(temperature, dt_s)
                    .map_err(|source| FiniteThermoplasticSheetError::Material {
                        element: e,
                        point: 0,
                        source,
                    })?,
            );
        }
        Ok(FrozenSheet { sheet: self, laws })
    }

    /// Frozen-history force, potential and optional tangent at an endpoint T/dt.
    /// Neither accepted time nor temperature changes during this evaluation.
    pub fn evaluate_trial(
        &self,
        temperatures_k: &[f64],
        displacement_m: &[[f64; 2]],
        dt_s: f64,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetEvaluation, FiniteThermoplasticSheetError> {
        let frozen = self.freeze(temperatures_k, dt_s, cx)?;
        let assembled = SheetOperator { sheet: &frozen }
            .assemble(
                &self.states,
                displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )?
            .0;
        Ok(FiniteThermoplasticSheetEvaluation {
            internal_force_n: assembled
                .force
                .chunks_exact(2)
                .map(|v| [v[0], v[1]])
                .collect(),
            tangent_n_m: assembled.tangent,
            incremental_potential_j: assembled.potential,
        })
    }

    /// Solve one real interval at supplied element temperatures. Old history is
    /// inspected at its old T, while every trial freezes the new T and same dt.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetStep<'a, L>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        let dofs = 2 * self.mesh.nodes.len();
        let external = validate_step_inputs(
            dofs,
            self.budget,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
        )?;
        let frozen = self.freeze(new_temperatures_k, dt_s, cx)?;
        let physical_time_s = self.physical_time_s + dt_s;
        if !physical_time_s.is_finite() || physical_time_s <= self.physical_time_s {
            return Err(invalid("thermal sheet physical time cannot advance representably").into());
        }
        for (e, element) in self.elements.iter().enumerate() {
            let old_law = self
                .family(e)
                .at_temperature_for_inspection(self.temperatures_k[e])
                .map_err(|source| FiniteThermoplasticSheetError::Material {
                    element: e,
                    point: 0,
                    source,
                })?;
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                old_law
                    .validate_committed(
                        &self.deformation[index],
                        &self.states[index],
                        controls.yield_tolerance_pa,
                    )
                    .map_err(|source| FiniteThermoplasticSheetError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
            }
        }
        let operator = SheetOperator { sheet: &frozen };
        let solved = increment::solve(
            &operator,
            IncrementControls {
                dofs,
                prescribed_m,
                nodal_forces_n,
                load_path: &[1.0],
                settings,
            },
            Some((self.displacement_m.as_flattened(), &self.states)),
            None,
            cx,
        )?;
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .expect("one physical endpoint");
        let assembled = operator
            .assemble(&self.states, &solved.displacement_m, None, false, cx)?
            .0;
        let mut energy = FiniteThermoplasticSheetEnergy::default();
        let mut quadrature_energy = Vec::with_capacity(self.states.len());
        let mut deformation = Vec::with_capacity(self.states.len());
        let mut piola_pa = Vec::with_capacity(self.states.len());
        let (mut old_free, mut new_free, mut old_internal, mut new_internal) = (0.0, 0.0, 0.0, 0.0);
        let tolerance = controls.equilibrium_work_tolerance_j;
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let map = |source| FiniteThermoplasticSheetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let response = frozen
                    .response(e, q, &solved.displacement_m, &self.states)
                    .map_err(map)?;
                // Reconstruct the actual two-temperature account from ORIGINAL
                // accepted history, not a second return from the new state.
                let value = self
                    .family(e)
                    .checked_increment_for_step(
                        &self.deformation[index],
                        &self.states[index],
                        self.temperatures_k[e],
                        &response.deformation_gradient,
                        new_temperatures_k[e],
                        controls.yield_tolerance_pa,
                        dt_s,
                    )
                    .map_err(map)?;
                if response.material.state != solved.point_states[index]
                    || value.mechanical.state != solved.point_states[index]
                {
                    return Err(invalid("thermal sheet endpoint audit differs from its accepted frozen-history return").into());
                }
                let account = FiniteThermoplasticElementEnergy::from_increment(&value);
                let weight = point.volume_m3;
                energy_gate(
                    "thermal native point balance",
                    weight * account.balance_error_j_m3,
                    tolerance,
                )?;
                energy_gate(
                    "thermal native stored energy",
                    weight * (account.free_after_j_m3 - value.mechanical.stored_energy_j_m3),
                    tolerance,
                )?;
                if [
                    account.plastic_dissipation_j_m3,
                    account.creep_dissipation_j_m3,
                ]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
                {
                    return Err(invalid(format!(
                        "thermal sheet element {e} point {q} has invalid physical dissipation"
                    ))
                    .into());
                }
                old_free += weight * account.free_before_j_m3;
                new_free += weight * account.free_after_j_m3;
                old_internal += weight * account.internal_before_j_m3;
                new_internal += weight * account.internal_after_j_m3;
                energy.free_energy_change_j +=
                    weight * (account.free_after_j_m3 - account.free_before_j_m3);
                energy.internal_energy_change_j += weight * account.internal_energy_change_j_m3();
                energy.plastic_dissipation_j += weight * account.plastic_dissipation_j_m3;
                energy.creep_dissipation_j += weight * account.creep_dissipation_j_m3;
                energy.reversible_heat_j += weight * account.reversible_heat_j_m3;
                energy.caloric_heat_j += weight * account.caloric_heat_j_m3();
                energy.coefficient_storage_change_j +=
                    weight * account.coefficient_storage_change_j_m3;
                energy.thermoelastic_storage_change_j +=
                    weight * account.thermoelastic_storage_change_j_m3;
                energy.endpoint_work_defect_j += weight * account.endpoint_work_defect_j_m3;
                energy.internal_endpoint_work_j += weight * account.internal_endpoint_work_j_m3;
                energy.incremental_potential_j += weight * account.incremental_potential_j_m3;
                deformation.push(response.deformation_gradient);
                piola_pa.push(value.mechanical.first_piola_pa);
                quadrature_energy.push(account);
            }
        }
        let old_u = self.displacement_m.as_flattened();
        energy.dead_load_endpoint_work_j = external
            .iter()
            .enumerate()
            .map(|(dof, force)| force * (solved.displacement_m[dof] - old_u[dof]))
            .sum();
        energy.support_endpoint_work_j = prescribed_m
            .iter()
            .zip(&equilibrium.reactions_n)
            .map(|(&(dof, _), reaction)| reaction * (solved.displacement_m[dof] - old_u[dof]))
            .sum();
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        energy.balance_error_j = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j
            - energy.internal_energy_change_j
            - energy.caloric_heat_j
            - energy.endpoint_work_defect_j;
        for (quantity, error) in [
            (
                "thermal accepted free energy",
                old_free - self.free_energy_j,
            ),
            (
                "thermal accepted internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "thermal assembled incremental potential",
                assembled.potential - energy.incremental_potential_j,
            ),
            (
                "thermal free-energy increment",
                new_free - old_free - energy.free_energy_change_j,
            ),
            (
                "thermal internal-energy increment",
                new_internal - old_internal - energy.internal_energy_change_j,
            ),
            (
                "thermal caloric decomposition",
                energy.caloric_heat_j - energy.inelastic_dissipation_j() - energy.reversible_heat_j,
            ),
            (
                "thermal equilibrium endpoint work",
                energy.equilibrium_work_residual_j,
            ),
            ("thermal total endpoint balance", energy.balance_error_j),
        ] {
            energy_gate(quantity, error, tolerance)?;
        }
        energy_gate(
            "thermal endpoint work quadrature defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if [
            new_free,
            new_internal,
            energy.plastic_dissipation_j,
            energy.creep_dissipation_j,
            energy.reversible_heat_j,
            energy.caloric_heat_j,
            energy.coefficient_storage_change_j,
            energy.thermoelastic_storage_change_j,
            energy.incremental_potential_j,
            energy.internal_endpoint_work_j,
            energy.dead_load_endpoint_work_j,
            energy.support_endpoint_work_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "thermal sheet integrated physical accounts are not representable",
            )
            .into());
        }
        checkpoint(cx)?;
        Ok(FiniteThermoplasticSheetStep {
            sheet: Self {
                mesh: self.mesh,
                families: self.families,
                thickness_m: self.thickness_m,
                elements: Arc::clone(&self.elements),
                local: self.local,
                budget: self.budget,
                displacement_m: solved
                    .displacement_m
                    .chunks_exact(2)
                    .map(|u| [u[0], u[1]])
                    .collect(),
                states: solved.point_states,
                deformation,
                piola_pa,
                temperatures_k: new_temperatures_k.to_vec(),
                free_energy_j: new_free,
                internal_energy_j: new_internal,
                physical_time_s,
            },
            equilibrium,
            quadrature_energy,
            energy,
        })
    }
}
