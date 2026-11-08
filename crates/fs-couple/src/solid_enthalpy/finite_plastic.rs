//! Finite-deformation plasticity and reciprocal heat on the original owners.
//!
//! Every coupled proposal repeats the native total-Lagrangian endpoint solve
//! from the SAME accepted deformation, plastic metric, accumulated flow,
//! temperature and optional velocity. The original heat timeline receives
//! physical plastic dissipation plus signed reversible exchange exactly once.
//! The accepted pair retains both histories, original source identities and
//! reference mass. A failed proposal changes neither original owner.
//!
//! The material free energy is per original reference volume. The explicitly
//! selected stress-free zero-pressure caloric chart supplies caloric internal
//! energy; mechanical entropy and pressure work remain in the finite law.
//! Signed endpoint work defects and backward-Euler kinetic losses are reported
//! separately and never deposited as heat.
//!
//! Thermal cells and links form a reference network with independently supplied
//! conductances. Their stored geometry is the reference geometry, not a moving
//! boundary reconstruction. Observers retain their original temperature-driven
//! sources, link/reservoir conductances and conservative exchanges. This owner
//! does not infer current-area convection or a deformed Fourier operator from
//! those reference admittances. Cell temperatures are uniform approximations;
//! load intervals and thermal substeps need independent refinement.
//! The explicit [`conduction`] entries instead request complete current normal
//! transport measures from the same proposed finite solid, retaining reference
//! storage and independently gating the final geometry discrepancy. With
//! `experimental-solid-contact-heat`, translating rigid tooling can participate
//! in that same solve; its mechanical work and numerical contact defect remain
//! separate from physical material heat.
//!
//! Mechanical limits are those of the original compressible displacement-only
//! P1 finite solid, including its volumetric-locking limitations. Optional
//! inertia uses the same fixed reference mass and endpoint solve.

/// Native source conduction rebuilt at each solved finite-deformation trial.
pub mod conduction;
/// Sustained rigid-tool contact through the original plastic/thermal owners.
#[cfg(feature = "experimental-solid-contact-heat")]
pub mod contact;
pub(super) mod feedback;
pub mod impact;
pub mod material;
/// Retained thermal porous compaction on the same conservative heat timeline.
pub mod porous;
/// Held current-area pressure through the same native thermal feedback.
pub mod pressure;
/// Reciprocal finite sheet forming and native caloric transport.
pub mod sheet;
#[cfg(feature = "experimental-solid-contact-heat")]
pub use contact::{FiniteThermoplasticContactCandidate, FiniteThermoplasticContactEnergy};
pub use impact::{
    FiniteThermoplasticEllipticImpactError, FiniteThermoplasticImpactError,
    FiniteThermoplasticSourceImpactError,
};
pub use material::{
    CaloricReference, FiniteRecoveringThermoplasticMaterial, FiniteThermoplasticMaterial,
    FiniteThermoplasticMaterialFamily, FiniteThermoplasticQuery, ReferenceCaloricMaterial,
    ReferenceCaloricQuery,
};
pub use pressure::{FiniteThermoplasticPressureCandidate, FiniteThermoplasticPressureEnergy};

use std::collections::BTreeSet;
use std::convert::Infallible;

use crate::iqn_ils::{IqnIlsError, driver::CouplingMethod};
use fs_exec::Cx;
use fs_matdb::ClaimSelection;
use fs_material::phase::PhaseStateError;
use fs_solid::hyper2d::plastic::thermal::FiniteThermoplasticSheetError;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetEnergy, FiniteThermoplasticTetError,
    FiniteThermoplasticTetSolid,
};
use fs_solid::hyper3::{HyperTetSettings, HyperTetSolution};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::dynamics::InertialEnergy;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    StepConfig, ThermalLaw,
};

/// Source, native mechanical/thermal, physical account or iteration refusal.
#[derive(Debug)]
pub enum FiniteThermoplasticError<E = Infallible> {
    Invalid(&'static str),
    Source(String),
    Caloric(PhaseStateError),
    Solid(FiniteThermoplasticTetError),
    #[cfg(feature = "experimental-solid-contact-heat")]
    Contact(fs_solid::hyper3::plastic::thermal::contact::FiniteThermoplasticPlaneContactError),
    DuctileSolid(fs_solid::hyper3::plastic::ductile::FiniteDuctileTetError),
    DuctileSheet(fs_solid::hyper2d::plastic::ductile::FiniteDuctileSheetError),
    Sheet(FiniteThermoplasticSheetError),
    Porous(fs_solid::hyper3::plastic::porous::thermal::FiniteThermalPorousTetError),
    Thermal(DepositionError),
    /// Original IQN-ILS policy or arithmetic refusal.
    Accelerator(IqnIlsError),
    Observer(E),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
    },
    /// The accepted-temperature geometry still disagrees with the geometry
    /// used for transport, even if scalar heat and temperature gates pass.
    TransportNotConverged {
        iterations: usize,
        geometry_relative_defect: f64,
        geometry_tolerance: f64,
    },
}
impl<E> From<DepositionError> for FiniteThermoplasticError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E: std::fmt::Debug> std::fmt::Display for FiniteThermoplasticError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "finite thermoplastic coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteThermoplasticError<E> {}

/// Exactly one original caloric cell and material family per tetrahedron.
#[derive(Clone, Debug)]
pub struct FiniteThermoplasticCell<'a> {
    pub cell: CellId,
    pub query: FiniteThermoplasticQuery<'a>,
}

/// Numerical coupling gates; neither material calibration nor time accuracy.
#[derive(Clone, Copy, Debug)]
pub struct FiniteThermoplasticConfig {
    /// Positive bounded number of within-interval proposals, at most 256.
    pub maximum_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    /// Picard or IQN-ILS startup relaxation of trial specific caloric energy in (0,1].
    pub enthalpy_relaxation: f64,
}

/// Accepted finite solid and its exact thermal cell identities.
#[derive(Debug)]
pub struct FiniteThermoplasticSolid<
    'a,
    M: FiniteThermoplasticMaterialFamily = FiniteThermoplasticMaterial<'a>,
> {
    solid: FiniteThermoplasticTetSolid<'a, M>,
    cells: Vec<DepositedCell>,
    time_s: f64,
}

impl<'a> FiniteThermoplasticSolid<'a> {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<FiniteThermoplasticCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        Self::prepare_sources(
            thermal,
            nodes_m,
            tetrahedra,
            sources,
            budget,
            cx,
            FiniteThermoplasticMaterial::try_new,
        )
    }

    /// Bind original thermal hardening surfaces and caloric cells to the same
    /// native equilibrium/reciprocal-heat owner. Every source must explicitly
    /// select a tabulated thermal model. Selection and preparation budgets apply
    /// independently to each original material card before mechanical stepping.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_tabulated_hardening(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<FiniteThermoplasticCell<'a>>,
        curve_selection: ClaimSelection,
        maximum_source_samples: usize,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        Self::prepare_sources(
            thermal,
            nodes_m,
            tetrahedra,
            sources,
            budget,
            cx,
            |query, enthalpy| {
                FiniteThermoplasticMaterial::try_new_with_tabulated_hardening(
                    query,
                    enthalpy,
                    curve_selection,
                    maximum_source_samples,
                )
            },
        )
    }
}

impl<'a> FiniteThermoplasticSolid<'a, FiniteRecoveringThermoplasticMaterial<'a>> {
    /// Advance permanent finite creep, yielding and independent resistance on
    /// the same original caloric timeline. The selected source schema retains
    /// recovery history through every mechanics/heat trial and accepted step.
    pub fn try_new_with_recovery(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<FiniteThermoplasticCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        Self::prepare_sources(
            thermal,
            nodes_m,
            tetrahedra,
            sources,
            budget,
            cx,
            FiniteRecoveringThermoplasticMaterial::try_new_recovering,
        )
    }
}

impl<'a, M: FiniteThermoplasticMaterialFamily> FiniteThermoplasticSolid<'a, M> {
    fn prepare_sources(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<FiniteThermoplasticCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
        mut prepare: impl FnMut(
            FiniteThermoplasticQuery<'a>,
            f64,
        ) -> Result<M, FiniteThermoplasticError>,
    ) -> Result<Self, FiniteThermoplasticError> {
        poll(cx)?;
        if sources.len() != tetrahedra.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "exactly one source cell per tetrahedron required",
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
                        "finite structural cell is inactive",
                    ))?;
            if !used.insert(source.cell)
                || !thermal.uses_material_law(
                    cell.material_id,
                    source.query.caloric.curve() as &dyn ThermalLaw,
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "unique cell mapping and exact original caloric object required",
                ));
            }
            let material = prepare(source.query, cell.specific_h_j_kg)?;
            temperatures.push(
                material
                    .caloric_source()
                    .curve()
                    .state_at_specific_enthalpy(cell.specific_h_j_kg)
                    .map_err(FiniteThermoplasticError::Caloric)?
                    .temperature_k(),
            );
            materials.push(material);
            cells.push(cell.clone());
        }
        let solid = FiniteThermoplasticTetSolid::try_new(
            nodes_m,
            tetrahedra,
            materials,
            &temperatures,
            budget,
            cx,
        )
        .map_err(FiniteThermoplasticError::Solid)?;
        for (i, cell) in cells.iter().enumerate() {
            let volume = solid
                .reference_volume_m3(i)
                .expect("retained tetrahedron index");
            if !super::same_measure(cell.geometry.volume_m3, volume)
                || !super::same_measure(
                    cell.mass_kg,
                    volume * solid.families()[i].reference_density_kg_m3(),
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "cell inventory must match original tet volume and source reference density",
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

    /// Bind original sourced reference density to the shared consistent mass.
    /// The returned owner retains velocity across all coupled heat iterations.
    pub fn with_inertia(
        mut self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticError> {
        let densities: Vec<_> = self
            .solid
            .families()
            .iter()
            .map(FiniteThermoplasticMaterialFamily::reference_density_kg_m3)
            .collect();
        self.solid = self
            .solid
            .with_inertia(&densities, velocity_mps, cx)
            .map_err(FiniteThermoplasticError::Solid)?;
        Ok(self)
    }

    pub const fn mechanics(&self) -> &FiniteThermoplasticTetSolid<'a, M> {
        &self.solid
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }

    fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticError<E>> {
        self.binding().enthalpies(thermal, accepted, cx)
    }

    fn binding(&self) -> feedback::ThermalBinding<'_, '_, M> {
        feedback::ThermalBinding {
            cells: &self.cells,
            materials: FinitePlasticLawField::PerElement(self.solid.families()),
        }
    }

    /// Advance actual mechanics and heat together. Loads and supports are held
    /// endpoint values; optional inertia uses the actual physical interval.
    /// The observer's mutable history must obey the original deep-clone contract.
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
    ) -> Result<
        FiniteThermoplasticCandidate<'a, 't, O, FiniteThermoplasticTetEnergy, M>,
        FiniteThermoplasticError<O::Error>,
    > {
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

    /// Advance the original owners with an explicit enthalpy coupling method.
    ///
    /// IQN-ILS changes only the next trial specific caloric energies [J/kg].
    /// Its bounded secant history belongs to this physical interval and never
    /// advances material, velocity, thermal or observer history. The original
    /// temperature, per-cell heat and both energy gates still decide acceptance.
    /// Source-domain and accelerator refusals propagate without clamping.
    /// Neither numerical convergence nor acceleration certifies physical accuracy.
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
    ) -> Result<
        FiniteThermoplasticCandidate<'a, 't, O, FiniteThermoplasticTetEnergy, M>,
        FiniteThermoplasticError<O::Error>,
    > {
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
                    .map_err(FiniteThermoplasticError::Solid)?;
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
                Ok(feedback::MechanicalProposal {
                    state: step.solid,
                    equilibrium: step.equilibrium,
                    energy: step.energy,
                    inertia: step.inertia,
                    heat_j,
                })
            },
        )?;
        Ok(FiniteThermoplasticCandidate {
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

/// Original mechanical and thermal reports plus independently checked coupling.
#[derive(Debug)]
pub struct FiniteThermoplasticReport<E = HyperTetSolution, A = FiniteThermoplasticTetEnergy> {
    pub iterations: usize,
    pub equilibrium: E,
    pub material_energy: A,
    pub inertia: Option<InertialEnergy>,
    pub thermal_audit: EnergyAudit,
    pub applied_heat_j: f64,
    /// Boundary/observer heat with this interval's plastic/reversible source removed.
    pub external_heat_j: f64,
    pub temperature_defect_k: f64,
    /// Sum of absolute cell heat discrepancies; opposing cells cannot cancel.
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    /// Caloric + mechanical internal + explicitly coupled gas storage + kinetic
    /// change + numerical defects, minus actual external load/pressure/tool/support
    /// work, boundary heat and incoming energy.
    pub combined_residual_j: f64,
}

/// Commit the complete triple together; borrowing it cannot mutate accepted owners.
pub struct FiniteThermoplasticCandidate<
    'm,
    't,
    O,
    A = FiniteThermoplasticTetEnergy,
    M: FiniteThermoplasticMaterialFamily = FiniteThermoplasticMaterial<'m>,
> {
    thermal: DepositionSimulation<'t>,
    solid: FiniteThermoplasticSolid<'m, M>,
    observer: O,
    report: FiniteThermoplasticReport<HyperTetSolution, A>,
}
impl<'m, 't, O, A, M: FiniteThermoplasticMaterialFamily>
    FiniteThermoplasticCandidate<'m, 't, O, A, M>
{
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &FiniteThermoplasticSolid<'m, M> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &FiniteThermoplasticReport<HyperTetSolution, A> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteThermoplasticSolid<'m, M>,
        O,
        FiniteThermoplasticReport<HyperTetSolution, A>,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteThermoplasticError<E>> {
    cx.checkpoint()
        .map_err(|_| FiniteThermoplasticError::Cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite<E>(value: f64) -> Result<f64, FiniteThermoplasticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FiniteThermoplasticError::Invalid(
            "unrepresentable finite thermoplastic heat/work arithmetic",
        ))
    }
}
