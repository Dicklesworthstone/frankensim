//! Conservative finite-volume osmosis with pressure and dilution feedback.
//!
//! The original `fs_material::osmosis` constitutive law owns water activity,
//! selectivity and exponentially fitted solute transport. This module evolves
//! actual solvent volumes and solute amounts, including initially absent
//! species. Mechanical pressure is a separate, explicit volume response.
//!
//! Each membrane is a fully implicit two-compartment exchange. At a proposed
//! water transfer, the solute equations have an analytic conservative solution;
//! only one scalar water balance needs iteration. Ordered membrane sweeps are
//! first-order Lie splitting. One full sweep versus two half sweeps controls
//! local temporal error; the two-half-step endpoint is retained. No uniqueness
//! theorem is claimed for the nonlinear permeating-solute residual, nor is the
//! local error estimate a certified global error bound.
//!
//! Compartments are well mixed and isothermal, with negligible solute volume.
//! A common declared species ordering, solvent and temperature is mandatory.
//! Optional finite neutral solid stocks equilibrate with the original material
//! solubility owner. Only dissolved molecules cross a bare selective membrane
//! or contribute to osmotic pressure; each compartment retains its crystals.
//! Solute transfer can exhaust or replenish that stock during the same implicit
//! water step. Solid-bearing cells exclude concentration films and ion models.
//! Alternatively, a finite-rate owner retains actual dissolved and crystal
//! stocks. Original-law fixed-volume phase half-stages surround each membrane sweep;
//! crystals remain fixed during transport, and step doubling controls their
//! coupling to dilution and permeation. The implicit membrane sweep keeps the
//! overall method first order. Declared conductances, optional self-similar
//! particle geometry and growth policies remain with the original material
//! owner. Particle area follows actual local crystal stock and retains its
//! reference across transport. No equilibrium target replaces an actual
//! inventory; no phase heat, particle breakup or nucleation law is inferred.
//! [`electrolyte`] retains original ionic mineral kinetics on this same water
//! timeline. Fully reflecting membranes carry solvent only; fast aqueous
//! speciation follows every candidate volume. Its conserved amounts are
//! analytical cation/family/spectator pools, not changing protonation forms.
//! Ideal selective membranes can use quasi-steady concentration films, which
//! reconstruct actual membrane-face concentrations without hidden inventories.
//! Nonideal compartments use their shared free-energy pressure with fully
//! retained solutes and bare membranes; leaky and film transport need matching
//! chemical-potential laws. Separate owners provide finite-storage films in
//! [`transient`], ideal or nonideal fast-ion water exchange in [`donnan`], and
//! finite-rate charged swelling in [`electrodiffusion`]. These isothermal owners
//! supply neither solvent inertia nor mechanical inertia or heat evolution.
//! Ordinary water transport and constant-data pairwise diffusion allocate no
//! memory. Reactive ionic water trials allocate original speciation payloads.
//! Material-sourced solution diffusion allocates transient original queries and
//! receipts; its prepared-storage estimate is not a peak-memory bound. Failed,
//! cancelled or exhausted advances publish neither state nor time.
//! [`coupled`] combines retained-water and independent molecular pathways using
//! the original pressure-law primitive and mixing free energy in one account.
//! Molecular paths can select a floating zero-current ionic closure on a forest,
//! including a charge-constrained implicit binary-salt stage. General closed
//! graphs can instead solve globally compatible diffusion potentials and
//! circulating currents with native checked CG and simultaneous species updates.
//! Optional global backward Euler covers multiple salts and neutral co-solutes
//! through conserved ionic coordinates and the original Newton--Krylov owner.
//! Every original node/species equation remains an acceptance condition. The
//! global options allocate bounded transient numerical payloads; accepted-state
//! observations expose fresh potentials and native species fluxes. All retain
//! the original source-qualified forces; none supplies double layers or
//! electrochemical reactions.

/// Nonideal Maxwell–Stefan diffusion at fixed solvent volume on this timeline.
pub mod diffusion;

/// Joint nonideal water and molecular transport with mechanical work accounted.
pub mod coupled;

/// Finite-rate ionic membrane permeation with coupled water and fixed charge.
pub mod electrodiffusion;

/// Finite-storage concentration films coupled to membrane water transport.
pub mod transient;

/// Finite water exchange with conservative, fast-ion Donnan equilibrium.
pub mod donnan;

/// Original ionic mineral kinetics coupled to fully retaining water membranes.
pub mod electrolyte;

mod mechanics;
mod solubility;

#[cfg(test)]
mod nonideal_tests;

#[cfg(test)]
use fs_material::osmosis::IdealDiluteSolution;
use fs_material::osmosis::polarization::PolarizedMembraneTransport;
use fs_material::osmosis::solubility::kinetics::{
    KineticSolubilityError, KineticSolubilitySolution,
};
use fs_material::osmosis::solubility::{FiniteSolubilityError, FiniteSolubilitySolution};
use fs_material::osmosis::{AqueousSolution, MembraneTransport, OsmosisError};
use fs_qty::Time;

/// Maximum independently conserved species, including any retained solid stocks.
pub const MAX_OSMOTIC_SPECIES: usize = 64;
const MAX_COMPARTMENTS: usize = 65_536;
const MAX_MEMBRANES: usize = 262_144;

/// Shared constitutive refusal or finite-volume integration failure.
#[derive(Debug)]
pub enum OsmoticNetworkError {
    /// Malformed geometry, layout, limits or incompatible thermodynamic state.
    Invalid(&'static str),
    /// The original solution or membrane law refused the state.
    Material(OsmosisError),
    /// The original nonlinear mechanical material refused its volume branch.
    Mechanical(fs_material::MaterialError),
    /// The original finite solid/dissolved partition refused the state.
    Solubility(FiniteSolubilityError),
    /// Original ionic mineral/speciation source or analytical phase balance.
    ElectrolyteSolubility(
        fs_material::osmosis::solubility::electrolyte::ElectrolyteSolubilityError,
    ),
    /// Source-owned solution thermodynamics or original MS friction refused.
    SolutionDiffusion(fs_material::osmosis::diffusion::SolutionDiffusionError),
    /// A smaller physical diffusion step may satisfy withdrawal or passivity.
    DiffusionStep(diffusion::SolutionDiffusionStepRefusal),
    /// Actual closed-window species inventory failed its absolute molar budget.
    DiffusionBalance { species: usize, residual_mol: f64 },
    /// Closed fixed-volume free energy increased beyond its whole-window budget.
    /// Kinetic crystals use the original solid-relative chemical reference.
    DiffusionFreeEnergy { change_j: f64, maximum_j: f64 },
    /// Closed coupled transport failed its absolute solvent-volume budget.
    CoupledWaterBalance { residual_m3: f64, maximum_m3: f64 },
    /// A smaller finite water stage may satisfy the total-potential budget.
    CoupledStepFreeEnergy { change_j: f64, maximum_j: f64 },
    /// The complete coupled window increased mixing plus mechanical potential.
    CoupledFreeEnergy { change_j: f64, maximum_j: f64 },
    /// Prepared numerical workspace exceeds the caller's byte budget.
    Workspace {
        required_bytes: usize,
        budget_bytes: usize,
    },
    /// No admissible scalar water root was bracketed at this step size.
    Unbracketed,
    /// Scalar solve exhausted its iteration bound.
    NonlinearConvergence,
    /// An implicit molecular solve could not resolve BOTH original
    /// inventory equations within the declared amount budget.
    ImplicitDiffusionConvergence { residual_mol: f64, maximum_mol: f64 },
    /// The native multicomponent Newton--Krylov solve exhausted or stalled.
    ImplicitMulticomponentConvergence {
        scaled_residual: f64,
        diagnosis: Option<fs_solver::NewtonStallDiagnosis>,
    },
    /// The native graph potential solve failed its recomputed residual gate.
    ElectroneutralPotentialConvergence { relative_residual: f64 },
    /// Actual original ionic fluxes violate charge continuity at a compartment.
    ElectroneutralChargeRate {
        compartment: usize,
        residual_mol_s: f64,
        maximum_mol_s: f64,
    },
    /// Time integration exhausted its accepted-step or rejection budget.
    WorkBudget,
    /// Further halving cannot advance floating-point time.
    StepUnderflow,
    /// Cancellation left accepted state and clock unchanged.
    Cancelled,
}

impl core::fmt::Display for OsmoticNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-volume osmosis: {self:?}")
    }
}
impl std::error::Error for OsmoticNetworkError {}
impl From<OsmosisError> for OsmoticNetworkError {
    fn from(value: OsmosisError) -> Self {
        Self::Material(value)
    }
}
impl From<FiniteSolubilityError> for OsmoticNetworkError {
    fn from(value: FiniteSolubilityError) -> Self {
        Self::Solubility(value)
    }
}
impl From<KineticSolubilityError> for OsmoticNetworkError {
    fn from(value: KineticSolubilityError) -> Self {
        match value {
            KineticSolubilityError::Material(error) => Self::Solubility(error),
            KineticSolubilityError::Cancelled => Self::Cancelled,
        }
    }
}

/// Mechanical resistance to solvent-volume change, independent of osmosis.
///
/// [`Self::new`] supplies `p_ref + elastance * (V - V_ref)`.
/// A zero elastance is a prescribed pressure. An open vertical reservoir of
/// constant cross-section A has elastance rho*g/A; an elastic chamber can
/// provide its tangent pressure/volume response. [`Self::from_bulk_material`]
/// instead retains an original nonlinear convex material bulk energy, its
/// pressure and its actual volume-dependent tangent. This pressure excludes
/// osmotic pressure, which the membrane law accounts for exactly once.
/// Pressures may use any COMMON reference; this law does not admit cavitation.
#[derive(Clone, Copy, Debug)]
pub struct OsmoticPressureVolume {
    reference_volume_m3: f64,
    reference_pressure_pa: f64,
    elastance_pa_per_m3: f64,
    bulk: Option<mechanics::BulkPressure>,
}

impl OsmoticPressureVolume {
    /// Admit a finite reference pressure, positive volume and nonnegative stiffness.
    pub fn new(
        reference_volume_m3: f64,
        reference_pressure_pa: f64,
        elastance_pa_per_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        if !positive(reference_volume_m3)
            || !reference_pressure_pa.is_finite()
            || !nonnegative(elastance_pa_per_m3)
        {
            return Err(OsmoticNetworkError::Invalid(
                "invalid mechanical pressure/volume law",
            ));
        }
        Ok(Self {
            reference_volume_m3,
            reference_pressure_pa,
            elastance_pa_per_m3,
            bulk: None,
        })
    }

    /// Mechanical pressure at the supplied solvent volume [Pa].
    pub fn pressure_pa(self, volume_m3: f64) -> Result<f64, OsmoticNetworkError> {
        let pressure = self.reference_pressure_pa + self.response_pressure_pa(volume_m3)?;
        if !pressure.is_finite() {
            return Err(OsmoticNetworkError::Invalid(
                "nonfinite pressure/volume response",
            ));
        }
        Ok(pressure)
    }

    /// Difference of the SAME pressure laws, cancelling common pressure
    /// references before arithmetic can hide a small elastic response.
    pub fn pressure_difference_pa(
        self,
        volume_m3: f64,
        other: Self,
        other_volume_m3: f64,
    ) -> Result<f64, OsmoticNetworkError> {
        let left = self.response_pressure_pa(volume_m3)?;
        let right = other.response_pressure_pa(other_volume_m3)?;
        let difference =
            (self.reference_pressure_pa - other.reference_pressure_pa) + left - right;
        if !difference.is_finite()
            || !(self.reference_pressure_pa + left).is_finite()
            || !(other.reference_pressure_pa + right).is_finite()
        {
            return Err(OsmoticNetworkError::Invalid(
                "pressure difference is not representable",
            ));
        }
        Ok(difference)
    }

    /// Exact pressure primitive change, with a COMMON linear pressure
    /// reference removed. This is a closed-volume accounting primitive; the
    /// actual solvent-volume residual must be checked independently. It supplies
    /// no caloric energy. Prescribed loads contribute their pressure work here.
    pub fn reference_cancelled_potential_change_j(
        self,
        before_volume_m3: f64,
        after_volume_m3: f64,
        common_pressure_pa: f64,
    ) -> Result<f64, OsmoticNetworkError> {
        self.pressure_pa(before_volume_m3)?;
        self.pressure_pa(after_volume_m3)?;
        let delta = after_volume_m3 - before_volume_m3;
        let change = match self.bulk {
            Some(bulk) => {
                (self.reference_pressure_pa - common_pressure_pa) * delta
                    + bulk.energy_change_j(
                        before_volume_m3 - self.reference_volume_m3,
                        after_volume_m3 - self.reference_volume_m3,
                    )?
            }
            None => delta
                * ((self.reference_pressure_pa - common_pressure_pa)
                    + self.elastance_pa_per_m3
                        * ((before_volume_m3 - self.reference_volume_m3) + 0.5 * delta)),
        };
        if !common_pressure_pa.is_finite() || !change.is_finite() {
            return Err(OsmoticNetworkError::Invalid(
                "mechanical potential difference is not representable",
            ));
        }
        Ok(change)
    }
}

/// One finite, well-mixed solution and its mechanical volume response.
#[derive(Debug)]
pub struct OsmoticCompartment {
    solution: AqueousSolution,
    solubility: Option<FiniteSolubilitySolution>,
    kinetic_solubility: Option<KineticSolubilitySolution>,
    kinetic_electrolyte: Option<electrolyte::KineticElectrolyte>,
    pressure: OsmoticPressureVolume,
    minimum_volume_m3: f64,
}

impl Clone for OsmoticCompartment {
    fn clone(&self) -> Self {
        Self {
            solution: self.solution.clone(),
            solubility: self
                .solubility
                .as_ref()
                .map(FiniteSolubilitySolution::independent_clone),
            kinetic_solubility: self.kinetic_solubility.clone(),
            kinetic_electrolyte: self.kinetic_electrolyte.clone(),
            pressure: self.pressure,
            minimum_volume_m3: self.minimum_volume_m3,
        }
    }
}

impl OsmoticCompartment {
    /// The explicit positive volume floor is an applicability boundary, never
    /// a clipped state or invented residual water inventory.
    pub fn new(
        solution: impl Into<AqueousSolution>,
        pressure: OsmoticPressureVolume,
        minimum_volume_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        let solution = solution.into();
        if !positive(minimum_volume_m3) || minimum_volume_m3 >= solution.solvent_volume_m3() {
            return Err(OsmoticNetworkError::Invalid(
                "volume floor must be below initial volume",
            ));
        }
        pressure.pressure_pa(solution.solvent_volume_m3())?;
        Ok(Self {
            solution,
            solubility: None,
            kinetic_solubility: None,
            kinetic_electrolyte: None,
            pressure,
            minimum_volume_m3,
        })
    }
    /// Prepare finite neutral dissolved/solid inventories for a bare selective
    /// membrane. Equilibrium is instantaneous and isothermal; crystals stay in
    /// their compartment while dissolved molecules can permeate. The original
    /// pair laws remain the sole owners of saturation and phase partitioning.
    /// Storage is made independent here, even if the caller supplied a clone.
    pub fn with_solubility(
        solution: FiniteSolubilitySolution,
        pressure: OsmoticPressureVolume,
        minimum_volume_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        let mut cell = Self::new(solution.solution().clone(), pressure, minimum_volume_m3)?;
        cell.solubility = Some(solution.independent_clone());
        Ok(cell)
    }
    /// Prepare actual neutral dissolved/solid stocks with declared finite-rate
    /// exchange. Native adaptive membrane steps evolve this same material owner;
    /// dissolution and precipitation never bypass its conductance, selected
    /// particle geometry or growth policy. Only actual dissolved molecules
    /// permeate or exert osmotic pressure.
    /// This constructor selects kinetics instead of instantaneous equilibrium.
    pub fn with_kinetic_solubility(
        solution: KineticSolubilitySolution,
        pressure: OsmoticPressureVolume,
        minimum_volume_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        let mut cell = Self::new(solution.solution().clone(), pressure, minimum_volume_m3)?;
        cell.kinetic_solubility = Some(solution);
        Ok(cell)
    }
    /// Actual finite solvent and dissolved inventories.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// Allocate an independent checkpoint of complete conserved inventories
    /// and their dissolved/solid split. Cloning or retaining this checkpoint
    /// cannot share mutable storage with the prepared transport buffers.
    /// `solution()` provides the allocation-free dissolved-only view.
    pub fn solubility_snapshot(&self) -> Option<FiniteSolubilitySolution> {
        self.solubility
            .as_ref()
            .map(FiniteSolubilitySolution::independent_clone)
    }
    /// Independent complete kinetic checkpoint, including actual crystal stocks
    /// and original rate laws. Its last receipts describe the final phase
    /// half-stage, not the net transfer over an entire network advance.
    pub fn kinetic_solubility_snapshot(&self) -> Option<KineticSolubilitySolution> {
        self.kinetic_solubility.clone()
    }
    /// Current mechanical pressure [Pa].
    pub fn pressure_pa(&self) -> Result<f64, OsmoticNetworkError> {
        self.pressure.pressure_pa(self.solution.solvent_volume_m3())
    }

    fn total_moles(&self, component: usize) -> f64 {
        if let Some(phases) = &self.kinetic_solubility {
            return phases.last_steps()[component]
                .state()
                .total_solute()
                .value();
        }
        self.solubility.as_ref().map_or_else(
            || self.solution.solute_moles()[component],
            |phases| phases.phase_partitions()[component].total_solute.value(),
        )
    }

    fn has_phase_inventory(&self) -> bool {
        self.solubility.is_some()
            || self.kinetic_solubility.is_some()
            || self.kinetic_electrolyte.is_some()
    }

    /// Molecular/energy paths may opt into the original neutral kinetic owner
    /// only. Ionic reaction requires different analytical and Gibbs accounts.
    fn unsupported_diffusion_phases(&self, neutral_kinetics: bool) -> bool {
        self.has_phase_inventory() && !(neutral_kinetics && self.kinetic_solubility.is_some())
    }

    /// Equilibrium cells expose their exchangeable total to the active-set
    /// solve. Kinetic crystals are frozen during transport, so only the actual
    /// dissolved amount is its coordinate for a kinetic cell.
    fn transport_moles(&self, component: usize) -> f64 {
        if self.solubility.is_some() {
            self.total_moles(component)
        } else {
            self.solution.solute_moles()[component]
        }
    }

    fn capacity_at_volume(
        &self,
        component: usize,
        volume_m3: f64,
    ) -> Result<Option<f64>, OsmoticNetworkError> {
        self.solubility
            .as_ref()
            .map(|phases| {
                Ok(phases
                    .partition_at_inventory(component, volume_m3, self.total_moles(component))?
                    .capacity
                    .value())
            })
            .transpose()
    }

    fn dissolved_moles_at(
        &self,
        component: usize,
        volume_m3: f64,
        total_moles: f64,
    ) -> Result<f64, OsmoticNetworkError> {
        match &self.solubility {
            Some(phases) => Ok(phases
                .partition_at_inventory(component, volume_m3, total_moles)?
                .dissolved
                .value()),
            None => Ok(total_moles),
        }
    }

    fn osmotic_pressure_at_total(
        &self,
        volume_m3: f64,
        totals: &[f64],
    ) -> Result<f64, OsmoticNetworkError> {
        if self.solubility.is_none() {
            return Ok(self.solution.osmotic_pressure_at(volume_m3, totals)?);
        }
        let mut dissolved = [0.0; MAX_OSMOTIC_SPECIES];
        for (i, &total) in totals.iter().enumerate() {
            dissolved[i] = self.dissolved_moles_at(i, volume_m3, total)?;
        }
        Ok(self
            .solution
            .osmotic_pressure_at(volume_m3, &dissolved[..totals.len()])?)
    }

    fn replace_inventories(
        &mut self,
        volume_m3: f64,
        totals: &[f64],
    ) -> Result<(), OsmoticNetworkError> {
        if self.kinetic_electrolyte.is_some() {
            return Err(OsmoticNetworkError::Invalid(
                "reactive ionic inventories require the original analytical owner, not species replacement",
            ));
        }
        if let Some(phases) = &mut self.kinetic_solubility {
            phases.replace_dissolved_inventories(volume_m3, totals)?;
            self.solution
                .replace_inventories(volume_m3, phases.solution().solute_moles())?;
            return Ok(());
        }
        match &mut self.solubility {
            Some(phases) => {
                // The prepared owner validates the complete endpoint before
                // mutation. Its identical aqueous view is then copied in place.
                phases.replace_inventories(volume_m3, totals)?;
                self.solution
                    .replace_inventories(volume_m3, phases.solution().solute_moles())?;
            }
            None => self.solution.replace_inventories(volume_m3, totals)?,
        }
        Ok(())
    }

    fn heap_bytes(&self) -> usize {
        self.solution.heap_bytes()
            + self
                .solubility
                .as_ref()
                .map_or(0, FiniteSolubilitySolution::heap_bytes)
            + self
                .kinetic_solubility
                .as_ref()
                .map_or(0, KineticSolubilitySolution::heap_bytes)
            + self
                .kinetic_electrolyte
                .as_ref()
                .map_or(0, |phase| phase.material.workspace_heap_bytes())
    }
}

/// One oriented membrane with a fixed physical area and the original material law.
#[derive(Clone, Debug)]
pub struct OsmoticMembrane {
    left: usize,
    right: usize,
    area_m2: f64,
    material: MembraneBoundary,
}

#[derive(Clone, Debug)]
enum MembraneBoundary {
    WellMixed(MembraneTransport),
    Polarized(PolarizedMembraneTransport),
}

impl OsmoticMembrane {
    /// Positive transfer is left to right. Species follow the network ordering.
    pub fn new(
        left: usize,
        right: usize,
        area_m2: f64,
        material: MembraneTransport,
    ) -> Result<Self, OsmoticNetworkError> {
        Self::from_boundary(left, right, area_m2, MembraneBoundary::WellMixed(material))
    }

    /// Add species-specific concentration films on either membrane face.
    ///
    /// Films have negligible storage and prescribed mass-transfer coefficients.
    /// The original membrane owns selectivity and hydraulic resistance. Its
    /// osmotic force uses reconstructed face concentrations, while the same
    /// species flux debits and credits the finite bulk inventories. Film viscous
    /// losses and transient boundary-layer storage are outside this model.
    pub fn new_polarized(
        left: usize,
        right: usize,
        area_m2: f64,
        material: PolarizedMembraneTransport,
    ) -> Result<Self, OsmoticNetworkError> {
        Self::from_boundary(left, right, area_m2, MembraneBoundary::Polarized(material))
    }

    fn from_boundary(
        left: usize,
        right: usize,
        area_m2: f64,
        material: MembraneBoundary,
    ) -> Result<Self, OsmoticNetworkError> {
        if left == right || !positive(area_m2) {
            return Err(OsmoticNetworkError::Invalid(
                "membrane needs distinct cells and positive area",
            ));
        }
        Ok(Self {
            left,
            right,
            area_m2,
            material,
        })
    }

    fn membrane_law(&self) -> &MembraneTransport {
        match &self.material {
            MembraneBoundary::WellMixed(material) => material,
            MembraneBoundary::Polarized(material) => material.membrane(),
        }
    }

    fn heap_bytes(&self) -> usize {
        match &self.material {
            MembraneBoundary::WellMixed(material) => core::mem::size_of_val(material.solutes()),
            MembraneBoundary::Polarized(material) => material.heap_bytes(),
        }
    }
}

/// Explicit temporal accuracy and work limits; all quantities use SI units.
#[derive(Clone, Copy, Debug)]
pub struct OsmoticStepControls {
    /// Largest attempted substep [s].
    pub maximum_step_s: f64,
    /// Relative local step-doubling tolerance, in (0,1).
    pub relative_tolerance: f64,
    /// Absolute solvent-volume tolerance [m³].
    pub absolute_volume_tolerance_m3: f64,
    /// Absolute per-species amount tolerance [mol], applied to dissolved
    /// inventories, complete totals and independently to kinetic crystal stocks.
    pub absolute_amount_tolerance_mol: f64,
    /// Maximum accepted substeps across the requested interval.
    pub maximum_steps: usize,
    /// Maximum total rejected trials across the requested interval.
    pub maximum_rejections: usize,
    /// Maximum bracket-search and safeguarded scalar iterations per exchange.
    pub maximum_root_iterations: usize,
}

impl Default for OsmoticStepControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1.0,
            relative_tolerance: 1e-5,
            absolute_volume_tolerance_m3: 1e-18,
            absolute_amount_tolerance_mol: 1e-18,
            maximum_steps: 10_000,
            maximum_rejections: 100,
            maximum_root_iterations: 96,
        }
    }
}

/// Numerical work and actual closed-system inventory balances.
#[derive(Clone, Debug)]
pub struct OsmoticStepReport {
    /// Original and accepted time [s].
    pub interval_s: [f64; 2],
    /// Accepted two-half-step endpoints.
    pub accepted_steps: usize,
    /// Trials rejected for source, root or error-control reasons.
    pub rejected_steps: usize,
    /// Scalar water residuals or solution-diffusion face evaluations, including rejected trials.
    /// Material-source qualification queries are additional work: each sourced
    /// face evaluates three original claims per binary pair, and proposed/final
    /// sweep states repeat this qualification. They are not scalar residuals.
    /// Global electrical stages count one affine preparation and one native
    /// combined-force replay per active face; CG iterations are reported separately.
    pub residual_evaluations: usize,
    /// Positive-duration per-component kinetic phase advances, including
    /// discarded adaptive trials. Zero-duration transport target refreshes are
    /// additional work. Zero for networks without finite-rate phase owners.
    pub phase_evaluations: usize,
    /// Attempted original fast-aqueous volume re-speciations, including
    /// discarded water trials. Each uses its source's bounded pH controls.
    /// Native mineral advance calls are counted by `phase_evaluations`; their
    /// nested bounded aqueous solves are additional work, not included here.
    pub electrolyte_volume_speciations: usize,
    /// Largest accepted step-doubling error divided by its requested scale.
    pub maximum_scaled_error: f64,
    /// Final minus initial sum of compartment solvent volumes [m³].
    pub water_balance_residual_m3: f64,
    /// Final minus initial TOTAL amount for each species [mol], including
    /// retained solids; only the first `species_count` entries are active.
    /// Reactive ionic networks have no active entries: their conserved
    /// analytical inventories are reported by `electrolyte_balances` instead.
    pub solute_balance_residual_mol: [f64; MAX_OSMOTIC_SPECIES],
    /// Number of active entries in the amount balance.
    pub species_count: usize,
    /// Original analytical inventories checked per compartment for fully
    /// retaining reactive ionic water transport. `Some` distinguishes this
    /// conservation axis explicitly; individual H+/OH-/family-form amounts
    /// may change. Residuals use independent dissolved/solid endpoint changes
    /// and the outer absolute-plus-relative amount scale. No reaction Gibbs
    /// or heat balance is inferred from these amount balances.
    pub electrolyte_balances: Option<Vec<electrolyte::ElectrolyteCompartmentBalance>>,
    /// Reference-cancelled original-to-final free-energy change [J].
    /// With kinetic solids this uses their original solid-relative chemical
    /// potential and includes dissolution/precipitation as well as mixing.
    /// Some only for fixed-volume solution diffusion; membrane water work is
    /// not silently dropped from a mixed mechanical/chemical energy account.
    pub diffusion_free_energy_change_j: Option<f64>,
    /// Closed coupled-window source mixing and pressure-primitive changes.
    /// Present for ordinary and finite-rate crystal coupled transport;
    /// no heat balance. Crystal networks retain their common solid reference.
    pub coupled_potential_change: Option<coupled::CoupledPotentialChange>,
    /// Global electrical solves, including discarded adaptive probes. Present
    /// only for the deliberately selected electroneutral network closure.
    pub electroneutral_diffusion: Option<diffusion::ElectroneutralDiffusionReport>,
}

/// Prepared closed network. Public compartments are read-only checkpoints;
/// cloning the network retains complete physical and numerical state.
#[derive(Clone, Debug)]
pub struct OsmoticNetwork {
    species: Vec<String>,
    cells: Vec<OsmoticCompartment>,
    membranes: Vec<OsmoticMembrane>,
    diffusion: Vec<diffusion::SolutionDiffusionLink>,
    diffusion_controls: Option<diffusion::SolutionDiffusionControls>,
    implicit_diffusion: Option<diffusion::ImplicitBinaryDiffusion>,
    implicit_multicomponent: Option<diffusion::ImplicitMulticomponentDiffusion>,
    implicit_salt: Option<diffusion::ImplicitSaltDiffusion>,
    electroneutral_diffusion: Option<diffusion::ElectroneutralDiffusion>,
    implicit_electroneutral: Option<diffusion::ImplicitElectroneutralDiffusion>,
    coupled_controls: Option<coupled::CoupledSolutionControls>,
    work: Vec<OsmoticCompartment>,
    coarse: Vec<OsmoticCompartment>,
    fine: Vec<OsmoticCompartment>,
    workspace_bytes: usize,
    time_s: f64,
}

impl OsmoticNetwork {
    /// Prepare retained numerical storage, binding one explicit ordering of
    /// aqueous species. Finite solid-bearing cells must retain identical original pair
    /// laws and the declared solute identities; kinetic conductances and seeding
    /// policies may differ by compartment. Named ions in ordinary aqueous
    /// cells must already be separate species; this function never infers
    /// dissociation or a van't Hoff multiplier.
    /// Reactive ionic cells instead conserve original analytical pools; their
    /// report marks that different conservation axis explicitly. Their byte
    /// budget counts retained numerical state, excluding pre-existing shared
    /// immutable chemical sources and transient speciation/report allocations.
    pub fn new(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        membranes: Vec<OsmoticMembrane>,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if cells.is_empty()
            || cells.len() > MAX_COMPARTMENTS
            || membranes.len() > MAX_MEMBRANES
            || species.len() > MAX_OSMOTIC_SPECIES
            || species
                .iter()
                .enumerate()
                .any(|(i, name)| name.trim().is_empty() || species[..i].contains(name))
        {
            return Err(OsmoticNetworkError::Invalid(
                "invalid compartment, membrane or species layout",
            ));
        }
        let first = &cells[0].solution;
        if cells.iter().any(|cell| {
            cell.solution.solute_moles().len() != species.len()
                || cell.solution.temperature_k().to_bits() != first.temperature_k().to_bits()
                || cell.solution.solvent_molar_volume_m3_mol().to_bits()
                    != first.solvent_molar_volume_m3_mol().to_bits()
        }) || membranes.iter().any(|m| {
            m.left >= cells.len()
                || m.right >= cells.len()
                || m.membrane_law().solutes().len() != species.len()
        }) {
            return Err(OsmoticNetworkError::Invalid(
                "all cells need one solvent, temperature and species ordering",
            ));
        }
        electrolyte::admit(&species, &cells, &membranes)?;
        let reference_phases = cells.iter().find_map(|cell| cell.solubility.as_ref());
        for phases in cells.iter().filter_map(|cell| cell.solubility.as_ref()) {
            if phases
                .components()
                .iter()
                .zip(&species)
                .any(|(component, name)| component.law().solute().as_str() != name.as_str())
                || reference_phases.is_some_and(|reference| {
                    phases.solvent() != reference.solvent()
                        || phases
                            .components()
                            .iter()
                            .zip(reference.components())
                            .any(|(component, original)| component.law() != original.law())
                })
            {
                return Err(OsmoticNetworkError::Invalid(
                    "finite solid stocks require the declared species order and identical original solvent/solute pair laws",
                ));
            }
        }
        let reference_kinetics = cells
            .iter()
            .find_map(|cell| cell.kinetic_solubility.as_ref());
        for phases in cells
            .iter()
            .filter_map(|cell| cell.kinetic_solubility.as_ref())
        {
            if phases
                .models()
                .iter()
                .zip(&species)
                .any(|(model, name)| model.law().solute().as_str() != name.as_str())
                || reference_phases.is_some_and(|reference| {
                    phases.solvent() != reference.solvent()
                        || phases
                            .models()
                            .iter()
                            .zip(reference.components())
                            .any(|(model, component)| model.law() != component.law())
                })
                || reference_kinetics.is_some_and(|reference| {
                    phases.solvent() != reference.solvent()
                        || phases
                            .models()
                            .iter()
                            .zip(reference.models())
                            .any(|(model, original)| model.law() != original.law())
                })
            {
                return Err(OsmoticNetworkError::Invalid(
                    "kinetic solid stocks require the declared species order and identical original solvent/solute pair laws",
                ));
            }
        }
        for membrane in &membranes {
            if (cells[membrane.left].has_phase_inventory()
                || cells[membrane.right].has_phase_inventory())
                && !matches!(membrane.material, MembraneBoundary::WellMixed(_))
            {
                return Err(OsmoticNetworkError::Invalid(
                    "finite solid stocks require a bare membrane; surface precipitation in concentration films is not modeled",
                ));
            }
            let nonideal = !cells[membrane.left].solution.is_ideal()
                || !cells[membrane.right].solution.is_ideal();
            if nonideal
                && (!matches!(membrane.material, MembraneBoundary::WellMixed(_))
                    || membrane
                        .membrane_law()
                        .solutes()
                        .iter()
                        .any(|(reflection, permeance)| {
                            reflection.value() != 1.0 || permeance.value() != 0.0
                        }))
            {
                return Err(OsmoticNetworkError::Invalid(
                    "nonideal solutions require fully retained solutes and a bare membrane; leaking and film transport need matching chemical-potential laws",
                ));
            }
        }
        // Four complete state buffers; input excess capacities are included.
        // Ordinary material/species payloads are charged conservatively. Ionic
        // cells charge cloned numerical/publication payloads, excluding their
        // pre-existing immutable Arc sources and transient native speciation.
        // Scalar stack scratch and returned report allocations are additional;
        // this is retained numerical storage, not a peak-memory bound.
        let cell_bytes = core::mem::size_of::<OsmoticCompartment>();
        let required = cells
            .capacity()
            .checked_add(3 * cells.len())
            .and_then(|n| n.checked_mul(cell_bytes))
            .and_then(|n| {
                let heaps = cells
                    .iter()
                    .try_fold(0usize, |total, cell| total.checked_add(cell.heap_bytes()))?;
                n.checked_add(heaps.checked_mul(4)?)
            })
            .and_then(|n| {
                n.checked_add(membranes.capacity() * core::mem::size_of::<OsmoticMembrane>())
            })
            .and_then(|n| n.checked_add(membranes.iter().map(OsmoticMembrane::heap_bytes).sum()))
            .and_then(|n| n.checked_add(species.capacity() * core::mem::size_of::<String>()))
            .and_then(|n| n.checked_add(species.iter().map(String::capacity).sum()))
            .ok_or(OsmoticNetworkError::Invalid(
                "workspace byte count overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        let (water, amounts) = totals(&cells, species.len());
        if !water.is_finite() || amounts[..species.len()].iter().any(|n| !n.is_finite()) {
            return Err(OsmoticNetworkError::Invalid(
                "network inventory total is not representable",
            ));
        }
        Ok(Self {
            species,
            work: cells.clone(),
            coarse: cells.clone(),
            fine: cells.clone(),
            cells,
            membranes,
            diffusion: Vec::new(),
            diffusion_controls: None,
            implicit_diffusion: None,
            implicit_multicomponent: None,
            implicit_salt: None,
            electroneutral_diffusion: None,
            implicit_electroneutral: None,
            coupled_controls: None,
            workspace_bytes: required,
            time_s: 0.0,
        })
    }
    /// Declared ordered species names.
    pub fn species(&self) -> &[String] {
        &self.species
    }
    /// Current accepted solution compartments.
    pub fn compartments(&self) -> &[OsmoticCompartment] {
        &self.cells
    }
    /// Current accepted time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Admitted numerical and immutable-input storage estimate [bytes].
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Advance exactly the requested interval with conservative, positive
    /// exchange and bounded adaptive substeps. Membrane stages are implicit;
    /// molecular stages use conservative fluxes at their current volumes, with
    /// deliberately selected pairwise backward Euler or global electrical
    /// continuity when configured. Pass
    /// `|| cx.checkpoint().is_err()` to use the shared execution scope.
    /// Cancellation is polled per scalar iteration, membrane and publication.
    pub fn advance(
        &mut self,
        duration_s: f64,
        controls: OsmoticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OsmoticStepReport, OsmoticNetworkError> {
        admit_controls(duration_s, controls)?;
        let end_time = self.time_s + duration_s;
        if !end_time.is_finite() || end_time <= self.time_s {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        poll(&mut cancelled)?;
        copy_state(&self.cells, &mut self.work)?;
        if let Some(implicit) = &mut self.implicit_diffusion {
            implicit.begin_window();
        }
        if let Some(implicit) = &mut self.implicit_multicomponent {
            implicit.begin_window();
        }
        if let Some(implicit) = &mut self.implicit_salt {
            implicit.begin_window();
        }
        if let Some(electrical) = &mut self.electroneutral_diffusion {
            electrical.begin_window();
        }
        if let Some(implicit) = &mut self.implicit_electroneutral {
            implicit.begin_window();
        }
        let (water_before, amounts_before) = totals(&self.cells, self.species.len());
        let reactive_ionic = self
            .cells
            .iter()
            .any(|cell| cell.kinetic_electrolyte.is_some());
        let mut report = OsmoticStepReport {
            interval_s: [self.time_s, end_time],
            accepted_steps: 0,
            rejected_steps: 0,
            residual_evaluations: 0,
            phase_evaluations: 0,
            electrolyte_volume_speciations: 0,
            maximum_scaled_error: 0.0,
            water_balance_residual_m3: 0.0,
            solute_balance_residual_mol: [0.0; MAX_OSMOTIC_SPECIES],
            species_count: if reactive_ionic {
                0
            } else {
                self.species.len()
            },
            electrolyte_balances: None,
            diffusion_free_energy_change_j: None,
            coupled_potential_change: None,
            electroneutral_diffusion: None,
        };
        let mut elapsed = 0.0;
        let mut step = duration_s.min(controls.maximum_step_s);
        while elapsed < duration_s {
            poll(&mut cancelled)?;
            if report.accepted_steps >= controls.maximum_steps {
                return Err(OsmoticNetworkError::WorkBudget);
            }
            let remaining = duration_s - elapsed;
            step = step.min(remaining);
            if elapsed + step <= elapsed || 0.5 * step == 0.0 {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            copy_state(&self.work, &mut self.coarse)?;
            copy_state(&self.work, &mut self.fine)?;
            let trial = (|| {
                sweep_all(
                    &mut self.coarse,
                    &self.membranes,
                    &mut self.diffusion,
                    self.implicit_diffusion.as_mut(),
                    self.implicit_multicomponent.as_mut(),
                    self.implicit_salt.as_mut(),
                    self.electroneutral_diffusion.as_mut(),
                    self.implicit_electroneutral.as_mut(),
                    self.diffusion_controls,
                    self.coupled_controls,
                    step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut report.phase_evaluations,
                    &mut report.electrolyte_volume_speciations,
                    &mut cancelled,
                )?;
                sweep_all(
                    &mut self.fine,
                    &self.membranes,
                    &mut self.diffusion,
                    self.implicit_diffusion.as_mut(),
                    self.implicit_multicomponent.as_mut(),
                    self.implicit_salt.as_mut(),
                    self.electroneutral_diffusion.as_mut(),
                    self.implicit_electroneutral.as_mut(),
                    self.diffusion_controls,
                    self.coupled_controls,
                    0.5 * step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut report.phase_evaluations,
                    &mut report.electrolyte_volume_speciations,
                    &mut cancelled,
                )?;
                sweep_all(
                    &mut self.fine,
                    &self.membranes,
                    &mut self.diffusion,
                    self.implicit_diffusion.as_mut(),
                    self.implicit_multicomponent.as_mut(),
                    self.implicit_salt.as_mut(),
                    self.electroneutral_diffusion.as_mut(),
                    self.implicit_electroneutral.as_mut(),
                    self.diffusion_controls,
                    self.coupled_controls,
                    0.5 * step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut report.phase_evaluations,
                    &mut report.electrolyte_volume_speciations,
                    &mut cancelled,
                )
            })();
            let error = match trial {
                Ok(()) => scaled_error(&self.work, &self.coarse, &self.fine, controls),
                Err(OsmoticNetworkError::Cancelled) => return Err(OsmoticNetworkError::Cancelled),
                Err(OsmoticNetworkError::DiffusionStep(_)) => f64::INFINITY,
                Err(OsmoticNetworkError::CoupledStepFreeEnergy { .. }) => f64::INFINITY,
                Err(
                    OsmoticNetworkError::Unbracketed | OsmoticNetworkError::NonlinearConvergence,
                ) if self.coupled_controls.is_some() => f64::INFINITY,
                Err(error) if self.diffusion_controls.is_some() => return Err(error),
                Err(_) => f64::INFINITY,
            };
            if !error.is_finite() || error > 1.0 {
                if report.rejected_steps >= controls.maximum_rejections {
                    return Err(OsmoticNetworkError::WorkBudget);
                }
                report.rejected_steps += 1;
                step *= 0.5;
                continue;
            }
            core::mem::swap(&mut self.work, &mut self.fine);
            report.accepted_steps += 1;
            report.maximum_scaled_error = report.maximum_scaled_error.max(error);
            elapsed = if step == remaining {
                duration_s
            } else {
                elapsed + step
            };
            // Fixed binary growth preserves deterministic decision boundaries;
            // a first-order method's local difference scales as step².
            if error < 0.125 {
                step = (2.0 * step).min(controls.maximum_step_s);
            }
        }
        let (water_after, amounts_after) = totals(&self.work, self.species.len());
        if !water_after.is_finite()
            || amounts_after[..self.species.len()]
                .iter()
                .any(|n| !n.is_finite())
        {
            return Err(OsmoticNetworkError::Invalid(
                "network inventory total is not representable",
            ));
        }
        report.water_balance_residual_m3 = water_after - water_before;
        for i in 0..report.species_count {
            report.solute_balance_residual_mol[i] = amounts_after[i] - amounts_before[i];
        }
        if reactive_ionic {
            report.electrolyte_balances = Some(electrolyte::balances(
                &self.cells,
                &self.work,
                controls,
                &mut cancelled,
            )?);
        }
        let kinetic_phases = self
            .cells
            .iter()
            .any(|cell| cell.kinetic_solubility.is_some());
        if kinetic_phases {
            // Retain independent solid tails in the balance, even when adding
            // one to a large dissolved stock rounds to the same total. Sum
            // actual endpoint changes rather than subtracting large totals.
            report.solute_balance_residual_mol.fill(0.0);
            let mut correction = [0.0; MAX_OSMOTIC_SPECIES];
            for (before, after) in self.cells.iter().zip(&self.work) {
                poll(&mut cancelled)?;
                for i in 0..self.species.len() {
                    if let (Some(old), Some(next)) =
                        (&before.kinetic_solubility, &after.kinetic_solubility)
                    {
                        for delta in [
                            next.last_steps()[i].state().dissolved().value()
                                - old.last_steps()[i].state().dissolved().value(),
                            next.last_steps()[i].state().solid().value()
                                - old.last_steps()[i].state().solid().value(),
                        ] {
                            compensated_add(
                                &mut report.solute_balance_residual_mol[i],
                                &mut correction[i],
                                delta,
                            );
                        }
                    } else {
                        compensated_add(
                            &mut report.solute_balance_residual_mol[i],
                            &mut correction[i],
                            after.total_moles(i) - before.total_moles(i),
                        );
                    }
                }
            }
        }
        if let Some(physical) = self.diffusion_controls {
            // Sum actual endpoint changes, avoiding subtraction of two large
            // inventory totals. The same absolute budget covers the complete
            // window, independent of accepted/rejected stage counts.
            // Kinetic networks already accumulated BOTH phase changes above.
            // Replacing that account with dissolved-only changes would mistake
            // actual dissolution for creation of a transported component.
            if !kinetic_phases {
                report.solute_balance_residual_mol.fill(0.0);
                let mut correction = [0.0; MAX_OSMOTIC_SPECIES];
                for (old, next) in self.cells.iter().zip(&self.work) {
                    poll(&mut cancelled)?;
                    for i in 0..self.species.len() {
                        compensated_add(
                            &mut report.solute_balance_residual_mol[i],
                            &mut correction[i],
                            next.solution.solute_moles()[i] - old.solution.solute_moles()[i],
                        );
                    }
                }
            }
            for i in 0..self.species.len() {
                let residual = report.solute_balance_residual_mol[i];
                if !residual.is_finite() || residual.abs() > physical.amount_balance_tolerance_mol {
                    return Err(OsmoticNetworkError::DiffusionBalance {
                        species: i,
                        residual_mol: residual,
                    });
                }
            }
            if let Some(coupled) = self.coupled_controls {
                let mut water = 0.0;
                let mut correction = 0.0;
                for (before, after) in self.cells.iter().zip(&self.work) {
                    compensated_add(
                        &mut water,
                        &mut correction,
                        after.solution.solvent_volume_m3() - before.solution.solvent_volume_m3(),
                    );
                }
                coupled::admit_water_balance(water, coupled.water_balance_tolerance_m3)?;
                report.water_balance_residual_m3 = water;
                let change = coupled::potential_change(&self.cells, &self.work, &mut cancelled)?;
                if change.total_potential_change_j > physical.free_energy_tolerance_j {
                    return Err(OsmoticNetworkError::CoupledFreeEnergy {
                        change_j: change.total_potential_change_j,
                        maximum_j: physical.free_energy_tolerance_j,
                    });
                }
                report.coupled_potential_change = Some(change);
            } else {
                let change =
                    diffusion::free_energy_change(&self.cells, &self.work, &mut cancelled)?;
                if change > physical.free_energy_tolerance_j {
                    return Err(OsmoticNetworkError::DiffusionFreeEnergy {
                        change_j: change,
                        maximum_j: physical.free_energy_tolerance_j,
                    });
                }
                report.diffusion_free_energy_change_j = Some(change);
            }
        }
        poll(&mut cancelled)?;
        report.electroneutral_diffusion = self
            .electroneutral_diffusion
            .as_ref()
            .map(|electrical| electrical.report());
        core::mem::swap(&mut self.cells, &mut self.work);
        self.time_s = end_time;
        Ok(report)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), OsmoticNetworkError> {
    if cancelled() {
        Err(OsmoticNetworkError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_controls(duration: f64, c: OsmoticStepControls) -> Result<(), OsmoticNetworkError> {
    if !positive(duration)
        || !positive(c.maximum_step_s)
        || !positive(c.relative_tolerance)
        || c.relative_tolerance >= 1.0
        || !positive(c.absolute_volume_tolerance_m3)
        || !positive(c.absolute_amount_tolerance_mol)
        || c.maximum_steps == 0
        || c.maximum_root_iterations == 0
        || c.maximum_root_iterations > 256
    {
        return Err(OsmoticNetworkError::Invalid(
            "invalid temporal accuracy or work controls",
        ));
    }
    Ok(())
}
fn copy_state(
    from: &[OsmoticCompartment],
    to: &mut [OsmoticCompartment],
) -> Result<(), OsmoticNetworkError> {
    for (source, dest) in from.iter().zip(to) {
        if let Some(phases) = &source.kinetic_electrolyte {
            electrolyte::install(dest, phases.material.clone())?;
            continue;
        }
        if let Some(phases) = &source.kinetic_solubility {
            let target = dest
                .kinetic_solubility
                .as_mut()
                .ok_or(OsmoticNetworkError::Invalid("kinetic workspace layout"))?;
            target.copy_state_from(phases)?;
            dest.solution.replace_inventories(
                source.solution.solvent_volume_m3(),
                source.solution.solute_moles(),
            )?;
            continue;
        }
        let count = source.solution.solute_moles().len();
        let mut totals = [0.0; MAX_OSMOTIC_SPECIES];
        for (i, total) in totals[..count].iter_mut().enumerate() {
            *total = source.total_moles(i);
        }
        dest.replace_inventories(source.solution.solvent_volume_m3(), &totals[..count])?;
    }
    Ok(())
}
fn scaled_error(
    old: &[OsmoticCompartment],
    coarse: &[OsmoticCompartment],
    fine: &[OsmoticCompartment],
    controls: OsmoticStepControls,
) -> f64 {
    let mut error = 0.0_f64;
    for ((old, coarse), fine) in old.iter().zip(coarse).zip(fine) {
        let before = old.solution.solvent_volume_m3();
        let after = fine.solution.solvent_volume_m3();
        error = error.max(
            (after - coarse.solution.solvent_volume_m3()).abs()
                / (controls.absolute_volume_tolerance_m3
                    + controls.relative_tolerance * before.abs().max(after.abs())),
        );
        for ((before, coarse), after) in old
            .solution
            .solute_moles()
            .iter()
            .zip(coarse.solution.solute_moles())
            .zip(fine.solution.solute_moles())
        {
            error = error.max(
                (after - coarse).abs()
                    / (controls.absolute_amount_tolerance_mol
                        + controls.relative_tolerance * before.abs().max(after.abs())),
            );
        }
        if old.solubility.is_some() || old.kinetic_solubility.is_some() {
            // A saturated dissolved view can stay fixed while a finite solid
            // reservoir drains. Control the conserved stock independently.
            for i in 0..old.solution.solute_moles().len() {
                let before = old.total_moles(i);
                let after = fine.total_moles(i);
                error = error.max(
                    (after - coarse.total_moles(i)).abs()
                        / (controls.absolute_amount_tolerance_mol
                            + controls.relative_tolerance * before.abs().max(after.abs())),
                );
            }
        }
        if let (Some(before), Some(coarse), Some(after)) = (
            &old.kinetic_solubility,
            &coarse.kinetic_solubility,
            &fine.kinetic_solubility,
        ) {
            for ((old, coarse), fine) in before
                .last_steps()
                .iter()
                .zip(coarse.last_steps())
                .zip(after.last_steps())
            {
                let before = old.state().solid().value();
                let after = fine.state().solid().value();
                error = error.max(
                    (after - coarse.state().solid().value()).abs()
                        / (controls.absolute_amount_tolerance_mol
                            + controls.relative_tolerance * before.abs().max(after.abs())),
                );
            }
        }
        if let (Some(before), Some(coarse), Some(after)) = (
            &old.kinetic_electrolyte,
            &coarse.kinetic_electrolyte,
            &fine.kinetic_electrolyte,
        ) {
            let before = before.material.solid_amount().value();
            let after = after.material.solid_amount().value();
            error = error.max(
                (after - coarse.material.solid_amount().value()).abs()
                    / (controls.absolute_amount_tolerance_mol
                        + controls.relative_tolerance * before.abs().max(after.abs())),
            );
        }
    }
    error
}
fn totals(cells: &[OsmoticCompartment], species: usize) -> (f64, [f64; MAX_OSMOTIC_SPECIES]) {
    let mut volume = 0.0;
    let mut vc = 0.0;
    let mut amounts = [0.0; MAX_OSMOTIC_SPECIES];
    let mut corrections = [0.0; MAX_OSMOTIC_SPECIES];
    for cell in cells {
        compensated_add(&mut volume, &mut vc, cell.solution.solvent_volume_m3());
        // Reactive ionic inventories have a per-cell analytical axis. Their
        // changing H+/OH-/family forms are not conserved transport species.
        if cell.kinetic_electrolyte.is_some() {
            continue;
        }
        for i in 0..species {
            if let Some(phases) = &cell.kinetic_solubility {
                let state = phases.last_steps()[i].state();
                compensated_add(
                    &mut amounts[i],
                    &mut corrections[i],
                    state.dissolved().value(),
                );
                compensated_add(&mut amounts[i], &mut corrections[i], state.solid().value());
            } else {
                compensated_add(&mut amounts[i], &mut corrections[i], cell.total_moles(i));
            }
        }
    }
    (volume, amounts)
}
fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}

#[allow(clippy::too_many_arguments)]
fn sweep_all(
    cells: &mut [OsmoticCompartment],
    membranes: &[OsmoticMembrane],
    links: &mut [diffusion::SolutionDiffusionLink],
    implicit: Option<&mut diffusion::ImplicitBinaryDiffusion>,
    multicomponent: Option<&mut diffusion::ImplicitMulticomponentDiffusion>,
    salt: Option<&mut diffusion::ImplicitSaltDiffusion>,
    electrical: Option<&mut diffusion::ElectroneutralDiffusion>,
    implicit_electrical: Option<&mut diffusion::ImplicitElectroneutralDiffusion>,
    physical: Option<diffusion::SolutionDiffusionControls>,
    coupled: Option<coupled::CoupledSolutionControls>,
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    phase_evaluations: &mut usize,
    electrolyte_evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    if let Some(coupled) = coupled {
        let kinetic_phases = cells.iter().any(|cell| cell.kinetic_solubility.is_some());
        if kinetic_phases {
            advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)?;
        }
        coupled::sweep(
            cells,
            membranes,
            links,
            implicit,
            multicomponent,
            salt,
            electrical,
            implicit_electrical,
            coupled,
            dt,
            controls,
            evaluations,
            cancelled,
        )?;
        if kinetic_phases {
            advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)?;
            diffusion::qualify_sources(cells, links, cancelled)?;
        }
        return Ok(());
    }
    match physical {
        Some(physical) => {
            let kinetic_phases = cells.iter().any(|cell| cell.kinetic_solubility.is_some());
            if kinetic_phases {
                advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)?;
            }
            diffusion::sweep(
                cells,
                links,
                implicit,
                multicomponent,
                salt,
                electrical,
                implicit_electrical,
                dt,
                physical,
                evaluations,
                cancelled,
            )?;
            if kinetic_phases {
                advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)?;
                // Reaction changes the actual source query point after the
                // molecular stage; stale pair samples cannot qualify it.
                diffusion::qualify_sources(cells, links, cancelled)?;
            }
            Ok(())
        }
        None => {
            advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)?;
            sweep(
                cells,
                membranes,
                dt,
                controls,
                evaluations,
                electrolyte_evaluations,
                cancelled,
            )?;
            advance_phases(cells, 0.5 * dt, controls, phase_evaluations, cancelled)
        }
    }
}

/// Original-law phase evolution at each compartment's current solvent volume.
/// The native controller compares complete coupled trials; it never publishes
/// either this half-stage or a partially completed membrane sweep separately.
fn advance_phases(
    cells: &mut [OsmoticCompartment],
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    for cell in cells {
        if let Some(phases) = &mut cell.kinetic_solubility {
            poll(cancelled)?;
            *evaluations = evaluations
                .checked_add(phases.models().len())
                .ok_or(OsmoticNetworkError::WorkBudget)?;
            phases.advance(Time::new(dt), &mut *cancelled)?;
            cell.solution.replace_inventories(
                phases.solution().solvent_volume_m3(),
                phases.solution().solute_moles(),
            )?;
        }
        electrolyte::advance_phase(cell, dt, controls, evaluations, cancelled)?;
    }
    Ok(())
}

fn sweep(
    cells: &mut [OsmoticCompartment],
    membranes: &[OsmoticMembrane],
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    electrolyte_evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    for membrane in membranes {
        poll(cancelled)?;
        let endpoint = exchange_counted(
            &cells[membrane.left],
            &cells[membrane.right],
            membrane,
            dt,
            controls,
            evaluations,
            electrolyte_evaluations,
            cancelled,
        )?;
        if let Some([left, right]) = endpoint.electrolytes {
            electrolyte::install(&mut cells[membrane.left], left)?;
            electrolyte::install(&mut cells[membrane.right], right)?;
            continue;
        }
        let count = cells[membrane.left].solution.solute_moles().len();
        cells[membrane.left].replace_inventories(endpoint.volumes[0], &endpoint.left[..count])?;
        cells[membrane.right].replace_inventories(endpoint.volumes[1], &endpoint.right[..count])?;
    }
    Ok(())
}

struct PairEndpoint {
    volumes: [f64; 2],
    // TOTAL inventories for equilibrium cells; actual dissolved amounts for
    // kinetic and aqueous cells. Kinetic crystals are frozen in this substage.
    left: [f64; MAX_OSMOTIC_SPECIES],
    right: [f64; MAX_OSMOTIC_SPECIES],
    residual_m3: f64,
    // Original source-bound ionic endpoints, already re-speciated at the
    // candidate volumes. Carry these through acceptance without a second solve.
    electrolytes: Option<[
        fs_material::osmosis::solubility::electrolyte::kinetics::KineticElectrolyteSolubilitySolution;
        2
    ]>,
}

/// Eliminate implicit solute transfers analytically at a proposed water transfer.
fn pair_endpoint(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    area_time: f64,
    water: f64,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    let mut endpoint = PairEndpoint {
        volumes: [
            left.solution.solvent_volume_m3() - water,
            right.solution.solvent_volume_m3() + water,
        ],
        left: [0.0; MAX_OSMOTIC_SPECIES],
        right: [0.0; MAX_OSMOTIC_SPECIES],
        residual_m3: 0.0,
        electrolytes: None,
    };
    let [vl, vr] = endpoint.volumes;
    if !positive(vl) || !positive(vr) || vl < left.minimum_volume_m3 || vr < right.minimum_volume_m3
    {
        return Err(OsmoticNetworkError::Unbracketed);
    }
    let mut cl = [0.0; MAX_OSMOTIC_SPECIES];
    let mut cr = [0.0; MAX_OSMOTIC_SPECIES];
    let count = left.solution.solute_moles().len();
    let volume_flux_m_s = water / area_time;
    for i in 0..count {
        // Each film/membrane exponential fit is evaluated once for this trial.
        // The cached positive reconstruction is reused after the conservative
        // implicit species update. No nested water solve or film inventory.
        let polarization = match &membrane.material {
            MembraneBoundary::WellMixed(_) => None,
            MembraneBoundary::Polarized(material) => {
                Some(material.coefficients(i, volume_flux_m_s)?)
            }
        };
        let (kl, kr) = match &polarization {
            Some(coefficients) => coefficients.donor_coefficients(),
            None => membrane
                .membrane_law()
                .solute_stream_coefficients(i, volume_flux_m_s)?,
        };
        let a = area_time * kl / vl;
        let b = area_time * kr / vr;
        if !nonnegative(a) || !nonnegative(b) {
            return Err(OsmoticNetworkError::Invalid(
                "implicit solute coefficient overflow",
            ));
        }
        let nl = left.transport_moles(i);
        let nr = right.transport_moles(i);
        if left.solubility.is_some() || right.solubility.is_some() {
            // At fixed water transfer the original donor coefficients remain
            // constant. Each independent component therefore has one bounded
            // piecewise-affine inventory equation, not another transport solve.
            let capacities = [
                left.capacity_at_volume(i, vl)?,
                right.capacity_at_volume(i, vr)?,
            ];
            [endpoint.left[i], endpoint.right[i]] =
                solubility::finite_stock_exchange([nl, nr], capacities, [a, b])
                    .map_err(OsmoticNetworkError::Invalid)?;
        } else {
            let scale = 1.0_f64.max(a).max(b);
            let a = a / scale;
            let b = b / scale;
            let inverse_scale = 1.0 / scale;
            let denominator = inverse_scale + a + b;
            let transfer = (a * nl - b * nr) / denominator;
            // Use a shared signed transfer for resolved changes. When nearly all
            // of a donor transfers, subtraction would erase its small but finite
            // remainder (and possibly its significant concentration). Recover that
            // remainder from positive implicit weights, then debit the SAME amount
            // on the other side. No endpoint is independently clipped.
            if transfer > 0.5 * nl {
                endpoint.left[i] =
                    ((inverse_scale + b) / denominator) * nl + (b / denominator) * nr;
                endpoint.right[i] = nr + (nl - endpoint.left[i]);
            } else if -transfer > 0.5 * nr {
                endpoint.right[i] =
                    ((inverse_scale + a) / denominator) * nr + (a / denominator) * nl;
                endpoint.left[i] = nl + (nr - endpoint.right[i]);
            } else {
                endpoint.left[i] = nl - transfer;
                endpoint.right[i] = nr + transfer;
            }
        }
        let bulk_left = left.dissolved_moles_at(i, vl, endpoint.left[i])? / vl;
        let bulk_right = right.dissolved_moles_at(i, vr, endpoint.right[i])? / vr;
        (cl[i], cr[i]) = match &polarization {
            Some(coefficients) => {
                coefficients.membrane_face_concentrations(bulk_left, bulk_right)?
            }
            None => (bulk_left, bulk_right),
        };
    }
    let volume_flux = if left.solution.is_ideal() && right.solution.is_ideal() {
        membrane
            .membrane_law()
            .water_flux(
                left.solution.temperature_k(),
                left.pressure
                    .pressure_difference_pa(vl, right.pressure, vr)?,
                0.0,
                &cl[..count],
                &cr[..count],
            )?
            .volume_flux_m_s
    } else {
        // Admission guarantees complete retention and no dilute Fick films.
        // The same water chemical potential supplies pressure and activity;
        // independent ionic/solute leakage is not inferred from this force.
        let osmotic_drop = left.osmotic_pressure_at_total(vl, &endpoint.left[..count])?
            - right.osmotic_pressure_at_total(vr, &endpoint.right[..count])?;
        membrane.membrane_law().hydraulic_permeance().value()
            * (left
                .pressure
                .pressure_difference_pa(vl, right.pressure, vr)?
                - osmotic_drop)
    };
    endpoint.residual_m3 = water - area_time * volume_flux;
    if !endpoint.residual_m3.is_finite() {
        return Err(OsmoticNetworkError::Invalid(
            "nonfinite implicit water residual",
        ));
    }
    Ok(endpoint)
}

fn exchange(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    // Existing molecular/coupled callers refuse reactive ionic cells during
    // shared admission; retain their established entry point and arithmetic.
    exchange_counted(
        left,
        right,
        membrane,
        dt,
        controls,
        evaluations,
        &mut 0,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn exchange_counted(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    electrolyte_evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    let area_time = dt * membrane.area_m2;
    if !positive(area_time) {
        return Err(OsmoticNetworkError::StepUnderflow);
    }
    // The solve tolerance is tighter than temporal error control. Relative
    // residual scaling follows transfer, not compartment size: tiny exchange
    // must not disappear merely because the reservoirs are large.
    let atol = controls.absolute_volume_tolerance_m3 * 0.01;
    let rtol = controls.relative_tolerance.min(1e-8) * 0.01;
    let mut evaluate = |water| {
        poll(cancelled)?;
        *evaluations += 1;
        if left.kinetic_electrolyte.is_some() || right.kinetic_electrolyte.is_some() {
            electrolyte::pair_endpoint(
                left,
                right,
                membrane,
                area_time,
                water,
                electrolyte_evaluations,
                cancelled,
            )
        } else {
            pair_endpoint(left, right, membrane, area_time, water)
        }
    };
    let at_zero = evaluate(0.0)?;
    if at_zero.residual_m3.abs() <= atol {
        return Ok(at_zero);
    }
    // Continue outwards only along the initial residual's flow direction.
    // Starting at a distant volume floor can select a remote root or overflow
    // the constitutive law even when a nearby physical root is admissible.
    let direction = if at_zero.residual_m3 < 0.0 { 1.0 } else { -1.0 };
    let available = if direction > 0.0 {
        left.solution.solvent_volume_m3() - left.minimum_volume_m3
    } else {
        right.solution.solvent_volume_m3() - right.minimum_volume_m3
    }
    .next_down();
    let mut distance = at_zero.residual_m3.abs().min(0.125 * available);
    let mut previous = 0.0;
    let mut previous_residual = at_zero.residual_m3;
    let mut refused_bound = None;
    let mut bracket = None;
    for _ in 0..controls.maximum_root_iterations {
        if distance <= previous || distance > available {
            return Err(OsmoticNetworkError::Unbracketed);
        }
        match evaluate(direction * distance) {
            Ok(endpoint) => {
                let residual = endpoint.residual_m3;
                if residual.abs() <= atol + rtol * distance {
                    return Ok(endpoint);
                }
                if direction * residual >= 0.0 {
                    bracket = Some(if direction > 0.0 {
                        (previous, distance, previous_residual, residual)
                    } else {
                        (-distance, -previous, residual, previous_residual)
                    });
                    break;
                }
                previous = distance;
                previous_residual = residual;
                distance = refused_bound.map_or_else(
                    || (2.0 * distance).min(available),
                    |upper| distance + 0.5 * (upper - distance),
                );
            }
            Err(OsmoticNetworkError::Cancelled) => return Err(OsmoticNetworkError::Cancelled),
            Err(_) => {
                refused_bound = Some(distance);
                distance = previous + 0.5 * (distance - previous);
            }
        }
    }
    let (mut lo, mut hi, mut flo, mut fhi) = bracket.ok_or(OsmoticNetworkError::Unbracketed)?;
    let mut bisect_next = false;
    for _ in 0..controls.maximum_root_iterations {
        let width = hi - lo;
        let midpoint = lo + 0.5 * width;
        // A narrow sign-changing bracket controls the transferred volume even
        // when a steep polarized osmotic force amplifies floating-point noise
        // in the residual. This is a state-error stopping rule, not permission
        // to publish an arbitrary last iterate after work exhaustion.
        if width <= 2.0 * (atol + rtol * midpoint.abs()) {
            return evaluate(midpoint);
        }
        // Scaled secant weights avoid overflow in fhi-flo. A midpoint follows
        // any interpolation that failed to halve the bracket, so interpolation
        // cannot stagnate at one endpoint. Keep the near-origin bracket itself:
        // permeating-solute models do not have a global unique-root guarantee.
        let scale = (-flo).max(fhi);
        let low_weight = -flo / scale;
        let high_weight = fhi / scale;
        let fraction = low_weight / (low_weight + high_weight);
        let interpolate = !bisect_next && fraction > 0.01 && fraction < 0.99;
        let water = if interpolate {
            lo + fraction * width
        } else {
            midpoint
        };
        let endpoint = evaluate(water)?;
        let residual = endpoint.residual_m3;
        if residual.abs() <= atol + rtol * water.abs() {
            return Ok(endpoint);
        }
        if water == lo || water == hi {
            return Err(OsmoticNetworkError::NonlinearConvergence);
        }
        if residual < 0.0 {
            lo = water;
            flo = residual;
        } else {
            hi = water;
            fhi = residual;
        }
        bisect_next = interpolate && hi - lo > 0.5 * width;
        // Both endpoint signs stay explicit; no unbracketed Newton jump is used.
        debug_assert!(flo <= 0.0 && fhi >= 0.0);
    }
    Err(OsmoticNetworkError::NonlinearConvergence)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod polarization_tests;

#[cfg(test)]
mod kinetic_solubility_tests;
