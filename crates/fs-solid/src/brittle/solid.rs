//! Retained source-bound brittle history on the native P1 equilibrium solver.
//!
//! Each endpoint freezes the original accepted histories during Newton, then
//! independently checks the native constitutive work before publishing a new
//! solid. The fixed crack planes and widths remain bound to the original mesh.
//! Fracture dissipation is physical loss; the signed endpoint-work defect is
//! integration error and is never converted into damage or heat.
//!
//! Ordinary mechanical endpoints keep material conditions fixed. The explicit
//! [`temperature`] path transfers sourced coefficients and history before
//! re-equilibrating. The [`initiation`] path preserves elastic preload while
//! selecting fixed planes from a sourced tensile predictor. The [`dynamics`]
//! path retains reference mass and velocity through physical-time fracture.
//! The [`fracture`] path releases explicitly completed traction-free faces and
//! retains their history through native continuation. Fixed planes do not
//! rotate; closing released faces require an explicit contact response. The
//! [`continuation`] path follows proportional-force equilibrium branches through
//! tensile softening. Finite kinematics remain absent. The caller resolves
//! loading reversals and softening with physical endpoints.

/// Checked temperature/source transfer followed by native mechanical equilibration.
pub mod temperature;

/// Adaptive physical load intervals with bounded damage and endpoint-work error.
pub mod path;

/// Retained elastic prehistory and source-bound fixed-plane brittle handoff.
pub mod initiation;

/// Source-bound reference mass and retained physical-time brittle motion.
pub mod dynamics;

/// Joint sourced thermal eigenstrain and native mechanical equilibrium.
pub mod expansion;

/// Traction-free release of completed crack bands into retained physical fragments.
pub mod fracture;

/// History-preserving pseudo-arclength paths through tensile softening.
pub mod continuation;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::{BrittleEnergy, BrittleError, BrittleState};

use super::{TetCrackBandGeometry, TetCrackBandGeometryError};
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};

/// Independent, absolute whole-mesh work gates, all in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetControls {
    /// Native material work minus storage change and fracture loss.
    /// Also bounds incremental versus cumulative fracture-loss consistency.
    pub material_energy_tolerance_j: f64,
    /// External minus internal endpoint work at the converged equilibrium.
    pub equilibrium_work_tolerance_j: f64,
    /// Absolute internal endpoint work minus exact constitutive path work.
    /// A smaller physical step is needed when this integration gate refuses.
    pub maximum_endpoint_work_defect_j: f64,
}

/// Whole-mesh physical work and separately identified numerical defects.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BrittleTetEnergy {
    /// Recoverable elastic storage before this endpoint.
    pub stored_before_j: f64,
    /// Recoverable elastic storage after this endpoint.
    pub stored_after_j: f64,
    /// Exact native constitutive work along the straight strain increment.
    pub material_work_j: f64,
    /// Irreversible fracture loss during this endpoint, never healing.
    pub fracture_dissipation_j: f64,
    /// Material work minus storage change and fracture dissipation.
    pub material_balance_error_j: f64,
    /// Sum of reference volume times endpoint stress double-contracted with
    /// the actual strain change. Tensor shears carry their factor of two.
    pub internal_endpoint_work_j: f64,
    /// Endpoint applied nodal force dotted with the actual displacement change.
    pub dead_load_endpoint_work_j: f64,
    /// Endpoint support reaction dotted with actual prescribed motion.
    pub support_endpoint_work_j: f64,
    /// External endpoint work minus internal endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus exact constitutive path work. Signed
    /// quadrature error, not additional fracture dissipation or generated heat.
    pub endpoint_work_defect_j: f64,
}

impl BrittleTetEnergy {
    /// Applied/support endpoint work minus exact constitutive path work.
    #[must_use]
    pub fn external_endpoint_work_defect_j(&self) -> f64 {
        self.equilibrium_work_residual_j + self.endpoint_work_defect_j
    }
}

/// Original source owners, geometry and complete immutable accepted history.
///
/// Sources are explicit per element: even one material card may need distinct
/// prepared laws because the directional band widths differ. Construction
/// admits the virgin state; no independent displacement/history mutation or
/// source replacement is exposed.
#[derive(Debug)]
pub struct BrittleTetSolid<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    bands: Arc<[TetCrackBandGeometry]>,
    released_faces: Arc<[fracture::ReleasedBrittleFace]>,
    budget: TetAssemblyBudget,
    solution: SmallStrainTetSolution<BrittleState<'a>>,
    stored_energy_j: f64,
    fracture_dissipation_j: f64,
}

/// One accepted physical endpoint. Original equilibrium output is available
/// through `solid.solution()`; it is not duplicated or reconstructed.
#[derive(Debug)]
pub struct BrittleTetStep<'a> {
    pub solid: BrittleTetSolid<'a>,
    /// Original native work/storage/loss per reference volume, J/m3.
    pub element_energy: Vec<BrittleEnergy>,
    pub energy: BrittleTetEnergy,
}

/// Admission, original material/solver refusal, or a physical-work gate.
#[derive(Debug, Clone, PartialEq)]
pub enum BrittleTetError {
    Invalid(&'static str),
    MaterialCount {
        expected: usize,
        got: usize,
    },
    WidthMismatch {
        element: usize,
        required_m: f64,
        material_m: f64,
    },
    Geometry(TetCrackBandGeometryError),
    Material {
        element: usize,
        source: BrittleError,
    },
    Solve(SmallStrainTetError),
    EnergyGate {
        account: &'static str,
        defect_j: f64,
        tolerance_j: f64,
    },
    /// A released face needs a contact response before this endpoint can pass.
    ReleasedFaceClosure {
        face: u64,
        vertex: usize,
        opening_m: f64,
        tolerance_m: f64,
    },
    Cancelled,
}

impl core::fmt::Display for BrittleTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound brittle tetrahedron refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetError {}
impl From<SmallStrainTetError> for BrittleTetError {
    fn from(error: SmallStrainTetError) -> Self {
        match error {
            SmallStrainTetError::Cancelled => Self::Cancelled,
            error => Self::Solve(error),
        }
    }
}

impl<'a> BrittleTetSolid<'a> {
    /// Bind one original resolved law to each admitted reference tetrahedron.
    ///
    /// The law's fixed normal selects its footprint; its width must equal
    /// `V / A_projected` within 64 relative machine epsilons. The original
    /// law supplies each virgin response without assembling a matrix. Native
    /// mechanical geometry, positive-connectivity and free-DOF admission apply
    /// when advancing, under the retained assembly budget.
    pub fn try_new(
        geometry: &'a TetLoadGeometry3<'a>,
        materials: &'a [ResolvedBrittleCrackBand],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetError> {
        checkpoint(cx)?;
        if materials.len() != geometry.tetrahedra().len() {
            return Err(BrittleTetError::MaterialCount {
                expected: geometry.tetrahedra().len(),
                got: materials.len(),
            });
        }
        if geometry.nodes_m().len() > budget.maximum_nodes
            || materials.len() > budget.maximum_tetrahedra
        {
            return Err(BrittleTetError::Invalid(
                "reference mesh exceeds the assembly budget",
            ));
        }
        if geometry.nodes_m().len().checked_mul(3).is_none()
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
        {
            return Err(BrittleTetError::Invalid(
                "invalid reference geometry budget",
            ));
        }
        let mut bands = Vec::with_capacity(materials.len());
        let mut states = Vec::with_capacity(materials.len());
        let mut stresses = Vec::with_capacity(materials.len());
        for (element, source) in materials.iter().enumerate() {
            checkpoint(cx)?;
            let law = source.law();
            let band = TetCrackBandGeometry::from_element(geometry, element, law.normal())
                .map_err(BrittleTetError::Geometry)?;
            let required = band.band_width_m();
            let supplied = law.band_width_m();
            if (required - supplied).abs() > (64.0 * f64::EPSILON) * required.max(supplied) {
                return Err(BrittleTetError::WidthMismatch {
                    element,
                    required_m: required,
                    material_m: supplied,
                });
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
            released_faces: Arc::from([]),
            budget,
            solution,
            stored_energy_j: 0.0,
            fracture_dissipation_j: 0.0,
        })
    }

    /// Original immutable reference geometry, including exterior load faces.
    #[must_use]
    pub const fn geometry(&self) -> &'a TetLoadGeometry3<'a> {
        self.geometry
    }

    /// Complete original source bundles in connectivity order.
    #[must_use]
    pub const fn materials(&self) -> &'a [ResolvedBrittleCrackBand] {
        self.materials
    }

    /// Directional footprint and width bound to every accepted history.
    #[must_use]
    pub fn band_geometry(&self) -> &[TetCrackBandGeometry] {
        &self.bands
    }

    /// Accepted native displacement, stresses and original material histories.
    /// The increment record describes the latest endpoint only; virgin state
    /// has no physical increment record.
    #[must_use]
    pub const fn solution(&self) -> &SmallStrainTetSolution<BrittleState<'a>> {
        &self.solution
    }

    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }

    /// Cumulative original fracture dissipation over the accepted load path.
    #[must_use]
    pub const fn fracture_dissipation_j(&self) -> f64 {
        self.fracture_dissipation_j
    }

    /// Advance one physical endpoint through the existing Newton solve.
    ///
    /// Prescribed displacements and dead nodal forces are absolute target
    /// values, not increments. Force entries may repeat, but their count must
    /// not exceed the original DOF count. Native solver gates, immutable source
    /// history, irreversibility, work gates and cancellation all pass before a
    /// candidate is returned. Every refusal leaves `self` unchanged.
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetStep<'a>, BrittleTetError> {
        checkpoint(cx)?;
        if [
            controls.material_energy_tolerance_j,
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(BrittleTetError::Invalid(
                "finite nonnegative work tolerances required",
            ));
        }
        let old_u = self.solution.displacement_m.as_flattened();
        if prescribed_m.len() > old_u.len() || nodal_forces_n.len() > old_u.len() {
            return Err(BrittleTetError::Invalid(
                "endpoint load entries exceed reference DOF count",
            ));
        }
        let mut laws = Vec::with_capacity(self.materials.len());
        for source in self.materials {
            checkpoint(cx)?;
            laws.push(source.small_strain_law());
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
        let solution = problem.solve_from(old_u, &self.solution.element_states, cx)?;
        self.accept_equilibrium(solution, prescribed_m, nodal_forces_n, controls, cx)
    }

    /// Shared complete publication gate for native load and arclength solves.
    /// The caller must already have checked the actual free-force equation.
    fn accept_equilibrium(
        &self,
        solution: SmallStrainTetSolution<BrittleState<'a>>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetStep<'a>, BrittleTetError> {
        if [
            controls.material_energy_tolerance_j,
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(BrittleTetError::Invalid(
                "finite nonnegative work tolerances required",
            ));
        }
        let old_u = self.solution.displacement_m.as_flattened();
        self.require_released_faces(&solution.displacement_m, cx)?;
        let (mut energy, element_energy, cumulative_loss) = self.material_energy(&solution, cx)?;
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
        energy.material_balance_error_j = (energy.material_work_j
            - (energy.stored_after_j - energy.stored_before_j))
            - energy.fracture_dissipation_j;
        energy.equilibrium_work_residual_j = (energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j)
            - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - energy.material_work_j;
        energy_gate(
            "material work/storage/fracture",
            energy.material_balance_error_j,
            controls.material_energy_tolerance_j,
        )?;
        energy_gate(
            "cumulative fracture loss",
            (cumulative_loss - self.fracture_dissipation_j) - energy.fracture_dissipation_j,
            controls.material_energy_tolerance_j,
        )?;
        energy_gate(
            "support/dead-load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "signed endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if [
            energy.stored_after_j,
            energy.fracture_dissipation_j,
            cumulative_loss,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(BrittleTetError::Invalid(
                "finite nonnegative physical energy required",
            ));
        }
        checkpoint(cx)?;
        Ok(BrittleTetStep {
            solid: Self {
                geometry: self.geometry,
                materials: self.materials,
                bands: Arc::clone(&self.bands),
                released_faces: Arc::clone(&self.released_faces),
                budget: self.budget,
                solution,
                stored_energy_j: energy.stored_after_j,
                fracture_dissipation_j: cumulative_loss,
            },
            element_energy,
            energy,
        })
    }

    /// Recover original constitutive work once for static and inertial endpoints.
    /// Geometry, native state ownership and solver admission stay with callers.
    fn material_energy(
        &self,
        solution: &SmallStrainTetSolution<BrittleState<'a>>,
        cx: &Cx<'_>,
    ) -> Result<(BrittleTetEnergy, Vec<BrittleEnergy>, f64), BrittleTetError> {
        let mut energy = BrittleTetEnergy {
            stored_before_j: self.stored_energy_j,
            ..BrittleTetEnergy::default()
        };
        let mut cumulative_loss = 0.0;
        let mut element_energy = Vec::with_capacity(self.materials.len());
        for (element, source) in self.materials.iter().enumerate() {
            checkpoint(cx)?;
            let old = &self.solution.element_states[element];
            let strain = &solution.element_strain[element];
            let step = source
                .law()
                .trial(strain, old)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            if step.state != solution.element_states[element]
                || step.stress_pa != solution.element_stress[element]
                || step.state.damage() < old.damage()
                || step.state.dissipated_j_m3() < old.dissipated_j_m3()
                || step.energy.dissipated < 0.0
            {
                return Err(BrittleTetError::Invalid(
                    "native history/stress or irreversibility mismatch",
                ));
            }
            let volume = self.bands[element].reference_volume_m3();
            energy.stored_after_j =
                volume.mul_add(step.state.stored_energy_j_m3(), energy.stored_after_j);
            energy.material_work_j = volume.mul_add(step.energy.work, energy.material_work_j);
            energy.fracture_dissipation_j =
                volume.mul_add(step.energy.dissipated, energy.fracture_dissipation_j);
            cumulative_loss = volume.mul_add(step.state.dissipated_j_m3(), cumulative_loss);
            let mut endpoint_density = 0.0;
            for component in 0..6 {
                let weight = if component < 3 { 1.0 } else { 2.0 };
                let change = weight * (strain[component] - old.strain()[component]);
                endpoint_density = step.stress_pa[component].mul_add(change, endpoint_density);
            }
            energy.internal_endpoint_work_j =
                volume.mul_add(endpoint_density, energy.internal_endpoint_work_j);
            element_energy.push(step.energy);
        }
        Ok((energy, element_energy, cumulative_loss))
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), BrittleTetError> {
    cx.checkpoint().map_err(|_| BrittleTetError::Cancelled)
}

fn energy_gate(
    account: &'static str,
    defect_j: f64,
    tolerance_j: f64,
) -> Result<(), BrittleTetError> {
    if !defect_j.is_finite() || defect_j.abs() > tolerance_j {
        Err(BrittleTetError::EnergyGate {
            account,
            defect_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}
