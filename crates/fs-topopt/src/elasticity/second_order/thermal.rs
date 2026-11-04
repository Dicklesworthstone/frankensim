//! Prescribed-temperature topology design with the original elastic tensors.
//!
//! For dimensionless stiffness weights `m`, the same interpolation multiplies
//! stiffness and thermal stress: `K = sum m_e K_e`, `f = f_dead + sum m_e q_e`,
//! where `q_e = integral B^T C_e epsilon_thermal`. Temperature, integrated free
//! strain, material state, mesh and supports are held fixed during design.
//! This is not a derivative through a heat-conduction solve.
//!
//! The objective is equivalent-load compliance `c = f^T u`, not the physical
//! stored elastic energy `1/2 integral (epsilon-epsilon_thermal):C:(...)`.
//! Both are reported. Exact modulus derivatives are
//! `g_e = 2 q_e^T u - u^T K_e u` and
//! `(H v)_e = 2 q_e^T du - 2 du^T K_e u`, with
//! `K du = sum v_e (q_e-K_e u)`. One tangent solve per load suffices. Thermal
//! gradients can be positive, so the dead-load monotone OC update is unsuitable.

use super::{ComplianceError, ComplianceGoal, ComplianceSolveOptions, WeightedOp, finite, gather};
use crate::control::SolveControl;
use crate::elasticity::DensityElasticity;
use crate::pipeline::{DesignPipeline, LoadCase};
use fs_solid::linear3::{
    TetElasticMaterial, TetMaterialField, TetThermalStrainField, TetThermalStrainState,
};
use fs_solver::LinearOp;

/// One independently equilibrated thermal and mechanical load case.
#[derive(Debug, Clone, Copy)]
pub struct ThermalLoadCase<'a> {
    /// Dead nodal forces [N], including any forces applied at supports.
    pub force: &'a [f64],
    /// Nonnegative objective weight; loads are never summed across cases.
    pub weight: f64,
    /// Integrated, oriented free strain bound to each original elastic state.
    pub strains: TetThermalStrainField<'a>,
}

enum OwnedThermalField {
    Uniform(TetThermalStrainState),
    PerElement(Vec<TetThermalStrainState>),
}
impl OwnedThermalField {
    fn new(field: TetThermalStrainField<'_>, cells: usize) -> Result<Self, ComplianceError> {
        match field {
            TetThermalStrainField::Uniform(state) => Ok(Self::Uniform(*state)),
            TetThermalStrainField::PerElement(states) if states.len() == cells => {
                Ok(Self::PerElement(states.to_vec()))
            }
            TetThermalStrainField::PerElement(_) => Err(ComplianceError::InvalidInput(
                "one thermal state per elasticity cell is required",
            )),
        }
    }

    fn at(&self, cell: usize) -> &TetThermalStrainState {
        match self {
            Self::Uniform(state) => state,
            Self::PerElement(states) => &states[cell],
        }
    }

    fn borrowed(&self) -> TetThermalStrainField<'_> {
        match self {
            Self::Uniform(state) => TetThermalStrainField::Uniform(state),
            Self::PerElement(states) => TetThermalStrainField::PerElement(states),
        }
    }
}

/// Constant P1 cell response in global Mandel coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalElementResponse {
    pub total_strain_mandel: [f64; 6],
    pub elastic_strain_mandel: [f64; 6],
    /// Homogenized design stress `m_e C_e (epsilon-epsilon_thermal)` [Pa].
    /// This is not a solid-phase stress recovery in a perforated microstructure.
    pub effective_stress_mandel_pa: [f64; 6],
    pub elastic_energy_j: f64,
}

/// Physical response to one original load, with homogeneous fixed supports.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalCaseResponse {
    pub elements: Vec<ThermalElementResponse>,
    /// Force of each support on the body [N]; zero on free DOFs.
    /// Applied forces at fixed DOFs are included in this balance.
    pub support_reactions_n: Vec<f64>,
    pub elastic_energy_j: f64,
}

struct ThermalCase {
    field: OwnedThermalField,
    element_forces_n: Vec<[f64; 12]>,
    response: ThermalCaseResponse,
}

/// Immutable thermal design point, reusing the native compliance solver,
/// preconditioner, work budget and true-residual gate. Original material and
/// thermal states remain attached. Both total and elastic compatible strains
/// must satisfy each thermal state's declared small-strain domain.
pub struct ThermalComplianceGoal<'a> {
    base: ComplianceGoal<'a>,
    cases: Vec<ThermalCase>,
    gradient: Vec<f64>,
}

fn material_at(field: TetMaterialField<'_>, cell: usize) -> &TetElasticMaterial {
    match field {
        TetMaterialField::Uniform(material) => material,
        TetMaterialField::PerElement(materials) => &materials[cell],
    }
}

impl<'a> ThermalComplianceGoal<'a> {
    /// Solve all cases at positive dimensionless stiffness weights. Stiffness
    /// and thermal equivalent loads use the same weights, including any SIMP
    /// floor. The material-bearing constructor is required; an unbound legacy
    /// E/nu operator cannot silently acquire a thermal source identity.
    pub fn new(
        elasticity: &'a DensityElasticity,
        moduli: &[f64],
        loads: &[ThermalLoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("thermal-compliance-admission")?;
        let materials = elasticity
            .material_field()
            .ok_or(ComplianceError::InvalidInput(
                "thermal design requires an assigned native elastic material field",
            ))?;
        if moduli.len() != elasticity.cells()
            || moduli.is_empty()
            || moduli.iter().any(|m| !m.is_finite() || *m <= 0.0)
        {
            return Err(ComplianceError::InvalidInput(
                "one finite positive stiffness weight per cell is required",
            ));
        }
        if loads.is_empty() {
            return Err(ComplianceError::InvalidInput(
                "at least one thermal load case is required",
            ));
        }
        let mut fields = Vec::with_capacity(loads.len());
        let mut element_forces = Vec::with_capacity(loads.len());
        let mut total_forces = Vec::with_capacity(loads.len());
        for load in loads {
            if load.force.len() != elasticity.n || !load.weight.is_finite() || load.weight < 0.0 {
                return Err(ComplianceError::InvalidInput(
                    "load shape or weight is invalid",
                ));
            }
            let field = OwnedThermalField::new(load.strains, elasticity.cells())?;
            let mut total = load.force.to_vec();
            let mut local_forces = Vec::with_capacity(elasticity.cells());
            for (cell, tet) in elasticity.tets.iter().enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-compliance-load")?;
                }
                let state = field.at(cell);
                if state.elastic_material_state_identity()
                    != material_at(materials, cell).material_state_identity
                {
                    return Err(ComplianceError::InvalidInput(
                        "thermal strain and elastic material state identities differ",
                    ));
                }
                // CB maps engineering strain to physical stress. Its transpose
                // is work-conjugate to engineering eigenstrain, not Mandel strain.
                let strain = state.free_strain_mandel();
                let mut force = [0.0; 12];
                for (local, output) in force.iter_mut().enumerate() {
                    for (component, &epsilon) in strain.iter().enumerate() {
                        let engineering = epsilon
                            * if component < 3 {
                                1.0
                            } else {
                                core::f64::consts::SQRT_2
                            };
                        *output = elasticity.cb[cell][12 * component + local]
                            .mul_add(engineering, *output);
                    }
                    *output = finite(*output * elasticity.volumes[cell], "thermal element load")?;
                    let dof = 3 * tet[local / 3] as usize + local % 3;
                    total[dof] = finite(
                        moduli[cell].mul_add(*output, total[dof]),
                        "thermal nodal load",
                    )?;
                }
                local_forces.push(force);
            }
            fields.push(field);
            element_forces.push(local_forces);
            total_forces.push(total);
        }
        let effective_loads: Vec<_> = total_forces
            .iter()
            .zip(loads)
            .map(|(force, load)| LoadCase {
                force,
                weight: load.weight,
            })
            .collect();
        let base = ComplianceGoal::new(elasticity, moduli, &effective_loads, options, control)?;
        let mut cases = Vec::with_capacity(loads.len());
        let mut gradient = base.gradient.clone();
        for (index, (field, element_forces_n)) in fields.into_iter().zip(element_forces).enumerate()
        {
            let u = &base.displacements[index];
            let response = recover(
                elasticity,
                moduli,
                &field,
                &element_forces_n,
                u,
                loads[index].force,
                control,
            )?;
            for (cell, (tet, force)) in elasticity.tets.iter().zip(&element_forces_n).enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-compliance-gradient")?;
                }
                let local = gather(elasticity, tet, u);
                let work: f64 = force.iter().zip(local).map(|(f, u)| f * u).sum();
                gradient[cell] = finite(
                    gradient[cell] + 2.0 * loads[index].weight * work,
                    "thermal compliance gradient",
                )?;
            }
            cases.push(ThermalCase {
                field,
                element_forces_n,
                response,
            });
        }
        control.checkpoint("thermal-compliance-publish")?;
        Ok(Self {
            base,
            cases,
            gradient,
        })
    }

    /// Weighted equivalent-load compliance [J], with all thermal load derivatives.
    pub fn value(&self) -> f64 {
        self.base.value()
    }
    /// Derivative with respect to stiffness weights, not raw design densities.
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }
    pub fn displacements(&self) -> &[Vec<f64>] {
        self.base.displacements()
    }
    pub fn case_compliances(&self) -> &[f64] {
        self.base.case_compliances()
    }
    pub fn primal_relative_residuals(&self) -> &[f64] {
        self.base.primal_relative_residuals()
    }
    pub fn material_field(&self) -> TetMaterialField<'_> {
        self.base
            .elasticity
            .material_field()
            .expect("admitted thermal material field")
    }
    pub fn thermal_field(&self, index: usize) -> Option<TetThermalStrainField<'_>> {
        self.cases.get(index).map(|case| case.field.borrowed())
    }
    pub fn response(&self, index: usize) -> Option<&ThermalCaseResponse> {
        self.cases.get(index).map(|case| &case.response)
    }

    /// Exact unweighted derivative of one independently equilibrated case.
    pub fn load_case_gradient(
        &self,
        index: usize,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        let mut gradient = self.base.load_case_gradient(index, control)?;
        let u = &self.base.displacements[index];
        for (cell, tet) in self.base.elasticity.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-case-gradient")?;
            }
            let local = gather(self.base.elasticity, tet, u);
            let work: f64 = self.cases[index].element_forces_n[cell]
                .iter()
                .zip(local)
                .map(|(f, u)| f * u)
                .sum();
            gradient[cell] = finite(gradient[cell] + 2.0 * work, "thermal case gradient")?;
        }
        control.checkpoint("thermal-case-gradient-publish")?;
        Ok(gradient)
    }

    /// One native tangent solve per positive-weight case, no assembled Hessian.
    pub fn hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.admit_direction(direction, control)?;
        let mut product = vec![0.0; direction.len()];
        if direction.iter().all(|&v| v == 0.0) {
            return Ok(product);
        }
        for (index, &weight) in self.base.weights.iter().enumerate() {
            if weight == 0.0 {
                continue;
            }
            let action = self.case_hessian(index, direction, control)?;
            for (value, entry) in product.iter_mut().zip(action) {
                *value = finite(*value + weight * entry, "thermal compliance Hessian action")?;
            }
        }
        control.checkpoint("thermal-hessian-publish")?;
        Ok(product)
    }

    /// Unweighted Hessian action, including a case with zero aggregate weight.
    pub fn load_case_hessian_vector_product(
        &self,
        index: usize,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.admit_direction(direction, control)?;
        if index >= self.cases.len() {
            return Err(ComplianceError::InvalidInput(
                "load-case index out of range",
            ));
        }
        if direction.iter().all(|&v| v == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let product = self.case_hessian(index, direction, control)?;
        control.checkpoint("thermal-case-hessian-publish")?;
        Ok(product)
    }

    fn admit_direction(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), ComplianceError> {
        control.checkpoint("thermal-hessian-admission")?;
        if direction.len() != self.base.moduli.len() {
            return Err(ComplianceError::InvalidInput(
                "Hessian direction shape mismatch",
            ));
        }
        for (cell, &value) in direction.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-hessian-admission")?;
            }
            finite(value, "Hessian direction")?;
        }
        Ok(())
    }

    fn case_hessian(
        &self,
        index: usize,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        let el = self.base.elasticity;
        let u = &self.base.displacements[index];
        let forces = &self.cases[index].element_forces_n;
        let derivative = WeightedOp {
            elasticity: el,
            weights: direction,
            identity: false,
        };
        let mut rhs = vec![0.0; el.n];
        derivative.apply(u, &mut rhs);
        for value in &mut rhs {
            *value = -*value;
        }
        for (cell, tet) in el.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-tangent-rhs")?;
            }
            for local in 0..12 {
                let dof = 3 * tet[local / 3] as usize + local % 3;
                if el.free[dof] {
                    rhs[dof] = finite(
                        direction[cell].mul_add(forces[cell][local], rhs[dof]),
                        "thermal tangent load",
                    )?;
                }
            }
        }
        let (du, _) = self
            .base
            .solve(&rhs, control, "thermal-compliance-tangent")?;
        let mut result = self.base.bilinear(&du, u, control)?;
        for (cell, tet) in el.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-hessian-contraction")?;
            }
            let local = gather(el, tet, &du);
            let work: f64 = forces[cell].iter().zip(local).map(|(f, u)| f * u).sum();
            result[cell] = finite(2.0 * (work - result[cell]), "thermal case Hessian action")?;
        }
        Ok(result)
    }
}

fn recover(
    el: &DensityElasticity,
    moduli: &[f64],
    thermal: &OwnedThermalField,
    forces: &[[f64; 12]],
    u: &[f64],
    dead_force: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<ThermalCaseResponse, ComplianceError> {
    let materials = el
        .material_field()
        .expect("admitted thermal material field");
    let mut elements = Vec::with_capacity(el.cells());
    let mut support_reactions_n: Vec<f64> = dead_force.iter().map(|f| -*f).collect();
    let mut elastic_energy_j = 0.0;
    for (cell, tet) in el.tets.iter().enumerate() {
        if cell % 256 == 0 {
            control.checkpoint("thermal-response")?;
        }
        let local = gather(el, tet, u);
        let mut gradient = [[0.0; 3]; 3];
        for a in 0..4 {
            for i in 0..3 {
                for j in 0..3 {
                    gradient[i][j] =
                        local[3 * a + i].mul_add(el.gradients[cell][a][j], gradient[i][j]);
                }
            }
        }
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let total_strain_mandel = [
            gradient[0][0],
            gradient[1][1],
            gradient[2][2],
            (gradient[0][1] + gradient[1][0]) * q,
            (gradient[1][2] + gradient[2][1]) * q,
            (gradient[2][0] + gradient[0][2]) * q,
        ];
        let state = thermal.at(cell);
        let free = state.free_strain_mandel();
        let elastic_strain_mandel = core::array::from_fn(|i| total_strain_mandel[i] - free[i]);
        if total_strain_mandel
            .iter()
            .chain(&elastic_strain_mandel)
            .any(|v| !v.is_finite() || v.abs() > state.maximum_absolute_mandel_strain())
        {
            return Err(ComplianceError::InvalidInput(
                "total or elastic compatible strain exceeds the thermal state's small-strain domain",
            ));
        }
        let c = material_at(materials, cell).stiffness_mandel_pa();
        let effective_stress_mandel_pa: [f64; 6] = core::array::from_fn(|i| {
            moduli[cell]
                * c[i]
                    .iter()
                    .zip(elastic_strain_mandel)
                    .map(|(c, strain)| c * strain)
                    .sum::<f64>()
        });
        let energy = finite(
            0.5 * el.volumes[cell]
                * elastic_strain_mandel
                    .iter()
                    .zip(effective_stress_mandel_pa)
                    .map(|(e, s)| e * s)
                    .sum::<f64>(),
            "thermal elastic energy",
        )?;
        if energy < 0.0 || effective_stress_mandel_pa.iter().any(|s| !s.is_finite()) {
            return Err(ComplianceError::NonFinite(
                "thermal stress or nonnegative elastic energy",
            ));
        }
        elastic_energy_j = finite(elastic_energy_j + energy, "total thermal elastic energy")?;
        for (row, &free_force) in forces[cell].iter().enumerate() {
            let dof = 3 * tet[row / 3] as usize + row % 3;
            let internal: f64 = el.ke[cell][12 * row..12 * (row + 1)]
                .iter()
                .zip(local)
                .map(|(k, u)| k * u)
                .sum();
            support_reactions_n[dof] = finite(
                moduli[cell].mul_add(internal - free_force, support_reactions_n[dof]),
                "thermal reaction",
            )?;
        }
        elements.push(ThermalElementResponse {
            total_strain_mandel,
            elastic_strain_mandel,
            effective_stress_mandel_pa,
            elastic_energy_j: energy,
        });
    }
    for (reaction, &free) in support_reactions_n.iter_mut().zip(&el.free) {
        if free {
            *reaction = 0.0;
        }
    }
    Ok(ThermalCaseResponse {
        elements,
        support_reactions_n,
        elastic_energy_j,
    })
}

/// Complete raw-density pullback with the original physical linearization.
pub struct ThermalDesignCompliance<'a> {
    pub goal: ThermalComplianceGoal<'a>,
    /// Derivative through SIMP, projection and the transpose Helmholtz filter.
    pub gradient: Vec<f64>,
}

impl DesignPipeline {
    /// Evaluate thermal equivalent-load compliance through the original design
    /// pipeline. The filter, all independent equilibria and the final transpose
    /// share one work/cancellation budget. No public elasticity weights mutate.
    /// These sensitivities need a general constrained optimizer; the monotone
    /// dead-load OC driver must not be applied to a thermal objective.
    pub fn try_thermal_compliance_and_gradient<'a>(
        &self,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[ThermalLoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<ThermalDesignCompliance<'a>, ComplianceError> {
        control.checkpoint("thermal-design-admission")?;
        if rho.len() != elasticity.cells() || rho.len() != self.filter.cells() {
            return Err(ComplianceError::InvalidInput(
                "design, filter and elasticity cell counts must agree",
            ));
        }
        let (filtered, _, moduli) = self.try_forward(rho, control)?;
        let goal = ThermalComplianceGoal::new(elasticity, &moduli, loads, options, control)?;
        let gradient = self.try_pullback(&filtered, goal.gradient(), control)?;
        control.checkpoint("thermal-design-publish")?;
        Ok(ThermalDesignCompliance { goal, gradient })
    }
}
