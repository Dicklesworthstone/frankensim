//! Reciprocal heat and retained finite porous compaction on original owners.
//!
//! The native pressure-sensitive solid keeps its non-isochoric plastic metric,
//! total deformation, temperature and optional velocity. Every proposal starts
//! from the same accepted history. Its yield heat and signed reversible thermal
//! exchange enter the original caloric timeline exactly once; mechanical work
//! defects and backward-Euler kinetic losses remain separate accounts.
//!
//! This uses the shared enthalpy feedback driver, source admission and thermal
//! transport. Original cell volume and mass are reference inventory; current
//! solid density is rho0/det(F). Conductances belong to the supplied reference
//! thermal network, without an inferred moving-area transport approximation.

/// Independent sealed gas storage shares the original compaction/heat endpoint.
pub mod cavity;
pub mod material;
/// Current-area pressure drives the same retained thermal feedback.
pub mod pressure;
/// The original explicit cell/query protocol also binds porous families.
pub use super::FiniteThermoplasticCell as FiniteThermalPorousCell;
pub use cavity::{
    FiniteThermalPorousCavitiesCandidate, FiniteThermalPorousCavitiesEnergy,
    FiniteThermalPorousCavityEquilibrium,
};
pub use material::FiniteThermalPorousMaterialBinding;
pub use pressure::{FiniteThermalPorousPressureCandidate, FiniteThermalPorousPressureEnergy};

use std::collections::BTreeSet;

use crate::iqn_ils::driver::CouplingMethod;
use fs_exec::Cx;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::hyper3::plastic::porous::thermal::{
    FiniteThermalPorousTetSolid, FiniteThermoplasticControls,
};
use fs_solid::hyper3::{HyperTetSettings, HyperTetSolution};
use fs_solid::linear3::TetAssemblyBudget;
use fs_thermal::enthalpy::{
    DepositedCell, DepositionObserver, DepositionSimulation, StepConfig, ThermalLaw,
};

use super::feedback::{self, MechanicalProposal, ThermalBinding};
use super::{
    FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticReport, finite, poll,
};

/// Accepted source-bound finite porous solid and exact original thermal cells.
#[derive(Debug)]
pub struct FiniteThermalPorousSolid<'a> {
    solid: FiniteThermalPorousTetSolid<'a, FiniteThermalPorousMaterialBinding<'a>>,
    cells: Vec<DepositedCell>,
    time_s: f64,
}

impl<'a> FiniteThermalPorousSolid<'a> {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<FiniteThermalPorousCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        if sources.len() != tetrahedra.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "exactly one original source cell per porous tetrahedron required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        let mut materials = Vec::with_capacity(sources.len());
        let mut temperatures = Vec::with_capacity(sources.len());
        for source in sources {
            poll(cx)?;
            let cell =
                thermal
                    .cells()
                    .get(&source.cell)
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "porous structural cell is inactive",
                    ))?;
            if !used.insert(source.cell)
                || !thermal.uses_material_law(
                    cell.material_id,
                    source.query.caloric.curve() as &dyn ThermalLaw,
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "unique porous cell mapping and exact original caloric object required",
                ));
            }
            let material =
                FiniteThermalPorousMaterialBinding::try_new(source.query, cell.specific_h_j_kg)?;
            temperatures.push(
                material
                    .state_at_specific_enthalpy(cell.specific_h_j_kg)
                    .map_err(FiniteThermoplasticError::Caloric)?
                    .temperature_k(),
            );
            materials.push(material);
            cells.push(cell.clone());
        }
        let solid = FiniteThermalPorousTetSolid::try_new(
            nodes_m,
            tetrahedra,
            materials,
            &temperatures,
            budget,
            cx,
        )
        .map_err(FiniteThermoplasticError::Porous)?;
        for (i, cell) in cells.iter().enumerate() {
            poll(cx)?;
            let volume = solid
                .reference_volume_m3(i)
                .expect("retained tetrahedron index");
            if !super::super::same_measure(cell.geometry.volume_m3, volume)
                || !super::super::same_measure(
                    cell.mass_kg,
                    volume * solid.families()[i].reference_density_kg_m3(),
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "porous cell inventory must match original tet volume and source reference density",
                ));
            }
        }
        poll(cx)?;
        Ok(Self {
            solid,
            cells,
            time_s: thermal.time_s(),
        })
    }

    /// Use the original source density and retain velocity through heat proposals.
    pub fn with_inertia(
        mut self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        let densities: Vec<_> = self
            .solid
            .families()
            .iter()
            .map(FiniteThermalPorousMaterialBinding::reference_density_kg_m3)
            .collect();
        self.solid = self
            .solid
            .with_inertia(&densities, velocity_mps, cx)
            .map_err(FiniteThermoplasticError::Porous)?;
        Ok(self)
    }

    pub const fn mechanics(
        &self,
    ) -> &FiniteThermalPorousTetSolid<'a, FiniteThermalPorousMaterialBinding<'a>> {
        &self.solid
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }

    /// Current density from the original sourced inventory and actual total
    /// deformation. Caloric chart density metadata never replaces this mass.
    pub fn current_density_kg_m3(&self, cx: &Cx<'_>) -> Result<Vec<f64>, FiniteThermoplasticError> {
        let densities: Vec<_> = self
            .solid
            .families()
            .iter()
            .map(FiniteThermalPorousMaterialBinding::reference_density_kg_m3)
            .collect();
        self.solid
            .current_density_kg_m3(&densities, cx)
            .map_err(FiniteThermoplasticError::Porous)
    }

    fn binding(&self) -> ThermalBinding<'_, '_, FiniteThermalPorousMaterialBinding<'a>> {
        ThermalBinding {
            cells: &self.cells,
            materials: FinitePlasticLawField::PerElement(self.solid.families()),
        }
    }

    /// Advance mechanics and native transport together over the actual interval.
    /// Endpoint loads/supports are held; observer clones must retain independent
    /// mutable history as required by the original thermal owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousCandidate<'a, 't, O>, FiniteThermoplasticError<O::Error>> {
        self.advance_with_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            mechanical_controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Use the original bounded Picard or IQN-ILS enthalpy iteration. Numerical
    /// convergence does not establish constitutive calibration or time accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupling<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousCandidate<'a, 't, O>, FiniteThermoplasticError<O::Error>> {
        let next = feedback::advance_with_coupling(
            self.binding(),
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            method,
            cx,
            |temperatures, duration_s| {
                let step = self
                    .solid
                    .advance(
                        temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Porous)?;
                let heat_j = step
                    .element_energy
                    .iter()
                    .enumerate()
                    .map(|(i, energy)| {
                        finite(
                            self.solid
                                .reference_volume_m3(i)
                                .expect("retained tetrahedron index")
                                * energy.caloric_heat_j_m3(),
                        )
                    })
                    .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?;
                Ok(MechanicalProposal {
                    state: step.solid,
                    equilibrium: step.equilibrium,
                    energy: step.energy,
                    inertia: step.inertia,
                    heat_j,
                })
            },
        )?;
        Ok(FiniteThermalPorousCandidate {
            thermal: next.thermal,
            solid: Self {
                solid: next.state,
                cells: next.cells,
                time_s: target_s,
            },
            observer: next.observer,
            report: next.report,
        })
    }
}

/// Atomic native thermal/porous/observer proposal. A refusal mutates none of the
/// accepted inputs; consumers commit the complete candidate as one endpoint.
pub struct FiniteThermalPorousCandidate<
    'a,
    't,
    O,
    A = fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy,
    E = HyperTetSolution,
> {
    thermal: DepositionSimulation<'t>,
    solid: FiniteThermalPorousSolid<'a>,
    observer: O,
    report: FiniteThermoplasticReport<E, A>,
}
impl<'a, 't, O, A, E> FiniteThermalPorousCandidate<'a, 't, O, A, E> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &FiniteThermalPorousSolid<'a> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &FiniteThermoplasticReport<E, A> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteThermalPorousSolid<'a>,
        O,
        FiniteThermoplasticReport<E, A>,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}
