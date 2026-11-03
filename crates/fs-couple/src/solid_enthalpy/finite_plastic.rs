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
//!
//! Mechanical limits are those of the original compressible displacement-only
//! P1 finite solid, including its volumetric-locking limitations. Optional
//! inertia uses the same fixed reference mass and endpoint solve.

pub mod impact;
pub mod material;
pub use impact::FiniteThermoplasticImpactError;
pub use material::{CaloricReference, FiniteThermoplasticMaterial, FiniteThermoplasticQuery};

use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;

use fs_exec::Cx;
use fs_material::phase::PhaseStateError;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetEnergy, FiniteThermoplasticTetError,
    FiniteThermoplasticTetSolid,
};
use fs_solid::hyper3::{HyperTetSettings, HyperTetSolution};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::dynamics::InertialEnergy;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, HeatTransfer, StepConfig, ThermalLaw,
};

/// Source, native mechanical/thermal, physical account or iteration refusal.
#[derive(Debug)]
pub enum FiniteThermoplasticError<E = Infallible> {
    Invalid(&'static str),
    Source(String),
    Caloric(PhaseStateError),
    Solid(FiniteThermoplasticTetError),
    Thermal(DepositionError),
    Observer(E),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
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
    /// Relax the trial specific caloric energy in (0,1].
    pub enthalpy_relaxation: f64,
}

/// Accepted finite solid and its exact thermal cell identities.
#[derive(Debug)]
pub struct FiniteThermoplasticSolid<'a> {
    solid: FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
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
            let material =
                FiniteThermoplasticMaterial::try_new(source.query, cell.specific_h_j_kg)?;
            temperatures.push(
                material
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
            .map(FiniteThermoplasticMaterial::reference_density_kg_m3)
            .collect();
        self.solid = self
            .solid
            .with_inertia(&densities, velocity_mps, cx)
            .map_err(FiniteThermoplasticError::Solid)?;
        Ok(self)
    }

    pub const fn mechanics(
        &self,
    ) -> &FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>> {
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
        let mut result = Vec::with_capacity(self.cells.len());
        for (old, material) in self.cells.iter().zip(self.solid.families()) {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or(FiniteThermoplasticError::Invalid(
                    "finite structural cell became inactive",
                ))?;
            if cell.id != old.id
                || cell.segment_id != old.segment_id
                || cell.body_id != old.body_id
                || cell.material_id != old.material_id
                || cell.geometry != old.geometry
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || accepted && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                || !thermal.uses_material_law(
                    cell.material_id,
                    material.caloric_source().curve() as &dyn ThermalLaw,
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "finite structural identity, reference geometry, mass or accepted caloric state changed",
                ));
            }
            material
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FiniteThermoplasticError::Caloric)?;
            result.push(cell.specific_h_j_kg);
        }
        Ok(result)
    }

    fn temperatures<E>(&self, enthalpies: &[f64]) -> Result<Vec<f64>, FiniteThermoplasticError<E>> {
        self.solid
            .families()
            .iter()
            .zip(enthalpies)
            .map(|(material, h)| {
                material
                    .state_at_specific_enthalpy(*h)
                    .map(|s| s.temperature_k())
                    .map_err(FiniteThermoplasticError::Caloric)
            })
            .collect()
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
    ) -> Result<FiniteThermoplasticCandidate<'a, 't, O>, FiniteThermoplasticError<O::Error>> {
        poll(cx)?;
        let duration_s = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !positive(duration_s)
            || !target_s.is_finite()
            || config.maximum_iterations == 0
            || config.maximum_iterations > 256
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.enthalpy_relaxation)
            || config.enthalpy_relaxation > 1.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(FiniteThermoplasticError::Invalid(
                "matching clocks, positive interval and bounded independent thermal coupling controls required",
            ));
        }
        let mut enthalpies = self.enthalpies(thermal, true, cx)?;
        let mut defects = [f64::INFINITY; 4];
        // Every cached proposal still starts from this SAME accepted solid.
        // With unrelaxed thermal iteration, the preceding actual-temperature
        // endpoint is exactly the next trial; do not solve Newton twice for it.
        let mut next_trial = None;
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let trial_temperatures = self.temperatures(&enthalpies)?;
            let trial = if let Some(cached) = next_trial.take() {
                cached
            } else {
                self.solid
                    .advance(
                        &trial_temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Solid)?
            };
            let trial_heat: Vec<_> = trial
                .element_energy
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    finite(
                        self.solid
                            .reference_volume_m3(i)
                            .expect("retained tetrahedron index")
                            * e.caloric_heat_j_m3(),
                    )
                })
                .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?;
            let applied_heat_j = finite(trial_heat.iter().sum())?;
            let heat = FiniteHeat {
                solid: self,
                inner: observer.clone(),
                cx,
                watts: self
                    .cells
                    .iter()
                    .zip(&trial_heat)
                    .map(|(cell, q)| Ok((cell.id, finite(*q / duration_s)?)))
                    .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?,
            };
            let (next_thermal, heat, thermal_audit) =
                thermal.advance_observed(&heat, target_s, thermal_config, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })?;
            let final_h = self.enthalpies(&next_thermal, false, cx)?;
            let final_temperatures = self.temperatures(&final_h)?;
            let actual = if final_temperatures == trial_temperatures {
                trial
            } else {
                self.solid
                    .advance(
                        &final_temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Solid)?
            };
            let actual_heat: Vec<_> = actual
                .element_energy
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    finite(
                        self.solid
                            .reference_volume_m3(i)
                            .expect("retained tetrahedron index")
                            * e.caloric_heat_j_m3(),
                    )
                })
                .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?;
            let caloric_heat_j = finite(actual_heat.iter().sum())?;
            let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(
                thermal_audit.stored_energy_change_j
                    - thermal_audit.incoming_energy_j
                    - external_heat_j
                    - caloric_heat_j,
            )?;
            let inertia_account = actual.inertia.as_ref().map_or(0.0, |inertia| {
                inertia.kinetic_change_j + inertia.numerical_loss_j
            });
            let combined_residual_j = finite(
                actual.energy.internal_energy_change_j
                    + thermal_audit.stored_energy_change_j
                    + actual.energy.endpoint_work_defect_j
                    + inertia_account
                    - actual.energy.dead_load_endpoint_work_j
                    - actual.energy.support_endpoint_work_j
                    - external_heat_j
                    - thermal_audit.incoming_energy_j,
            )?;
            defects = [
                trial_temperatures
                    .iter()
                    .zip(&final_temperatures)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
                finite(
                    trial_heat
                        .iter()
                        .zip(&actual_heat)
                        .map(|(a, b)| (a - b).abs())
                        .sum(),
                )?,
                thermal_residual_j,
                combined_residual_j,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j
                && defects[3].abs() <= config.energy_tolerance_j
            {
                let cells = self
                    .cells
                    .iter()
                    .map(|cell| next_thermal.cells()[&cell.id].clone())
                    .collect();
                let solid = Self {
                    solid: actual.solid,
                    cells,
                    time_s: target_s,
                };
                let report = FiniteThermoplasticReport {
                    iterations: iteration,
                    equilibrium: actual.equilibrium,
                    material_energy: actual.energy,
                    inertia: actual.inertia,
                    thermal_audit,
                    applied_heat_j,
                    external_heat_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                    thermal_residual_j,
                    combined_residual_j,
                };
                poll(cx)?;
                return Ok(FiniteThermoplasticCandidate {
                    thermal: next_thermal,
                    solid,
                    observer: heat.inner,
                    report,
                });
            }
            if config.enthalpy_relaxation == 1.0 {
                enthalpies = final_h;
                next_trial = Some(actual);
            } else {
                for (h, new_h) in enthalpies.iter_mut().zip(final_h) {
                    *h = finite(
                        (1.0 - config.enthalpy_relaxation) * *h
                            + config.enthalpy_relaxation * new_h,
                    )?;
                }
            }
        }
        poll(cx)?;
        Err(FiniteThermoplasticError::NotConverged {
            iterations: config.maximum_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            thermal_residual_j: defects[2],
            combined_residual_j: defects[3],
        })
    }
}

/// Original mechanical and thermal reports plus independently checked coupling.
#[derive(Debug)]
pub struct FiniteThermoplasticReport {
    pub iterations: usize,
    pub equilibrium: HyperTetSolution,
    pub material_energy: FiniteThermoplasticTetEnergy,
    pub inertia: Option<InertialEnergy>,
    pub thermal_audit: EnergyAudit,
    pub applied_heat_j: f64,
    /// Boundary/observer heat with this interval's plastic/reversible source removed.
    pub external_heat_j: f64,
    pub temperature_defect_k: f64,
    /// Sum of absolute cell heat discrepancies; opposing cells cannot cancel.
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    /// Caloric + mechanical internal + kinetic change + numerical defects
    /// minus actual load/support work, boundary heat and incoming energy.
    pub combined_residual_j: f64,
}

/// Commit the complete triple together; borrowing it cannot mutate accepted owners.
pub struct FiniteThermoplasticCandidate<'m, 't, O> {
    thermal: DepositionSimulation<'t>,
    solid: FiniteThermoplasticSolid<'m>,
    observer: O,
    report: FiniteThermoplasticReport,
}
impl<'m, 't, O> FiniteThermoplasticCandidate<'m, 't, O> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &FiniteThermoplasticSolid<'m> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &FiniteThermoplasticReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        FiniteThermoplasticSolid<'m>,
        O,
        FiniteThermoplasticReport,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}

#[derive(Clone)]
struct FiniteHeat<'s, 'm, 'cx, O> {
    solid: &'s FiniteThermoplasticSolid<'m>,
    inner: O,
    cx: &'s Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
}
impl<O: DepositionObserver> DepositionObserver for FiniteHeat<'_, '_, '_, O> {
    type Error = FiniteThermoplasticError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner
            .begin()
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner
            .next_transition_time_s(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.solid.enthalpies(thermal, false, self.cx)?;
        let mut rates = self
            .inner
            .sources(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)?;
        for (&cell, &watts) in &self.watts {
            checkpoint()?;
            let rate = rates.entry(cell).or_insert(0.0);
            *rate = finite(*rate + watts)?;
        }
        Ok(rates)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .step_sources(thermal, start, dt, end, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }

    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .ambient_conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, dt, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(FiniteThermoplasticError::Observer)?;
        self.solid.enthalpies(thermal, false, self.cx).map(|_| ())
    }
    fn step(
        &mut self,
        start: f64,
        duration: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start, duration, old, thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)?;
        self.solid.enthalpies(thermal, false, self.cx).map(|_| ())
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
