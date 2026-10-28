//! Finite thermal feedback for the original retained partial-slip contact.
//!
//! This is an explicitly declared affine contact FREE-ENERGY model, not a
//! calibration inferred from a material name or from a bulk elastic modulus.
//! Each surface owns its own temperature and entropy contribution. For the
//! two original Jenkins regions, psi = sum(k(Ta,Tb) z^2/2), si = -dpsi/dTi,
//! and U_contact = psi + Ta*sa + Tb*sb. The caloric network receives physical
//! microslip heat and reversible exchange -Ci-d(Ti*si), never numerical loss.
//!
//! Fixed normal patch/load and an explicit prescribed or conduction-resolved
//! heat route. Coulomb coefficients are
//! fixed or re-queried through the original pinned material sliding binding.
//! no normal work, expansion, wear, adhesion, flash-temperature or resolved
//! traction-field claim. The existing caloric objects and transport sources
//! retain their own source domains. Every nonlinear trial starts from the same
//! accepted mechanical/thermal history. Endpoint coupling is first order;
//! nonlinear convergence is not a temporal or experimental accuracy claim.

use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_math::det;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, PartialSlipCheckpoint, PartialSlipError, PartialSlipKinematics,
    PartialSlipLaw, PartialSlipState, PartialSlipStep, TangentFrame,
};
use fs_tribo::{FrictionLaw, HeatPartition, InputAuthority};
use std::cell::{Cell, RefCell};
use std::convert::Infallible;

use super::{
    MaterialSlidingContact, SlidingError, SlidingHeatPartition, SlidingInterface,
    SlidingMaterialSample,
};

/// A common thermal interface partitions physical loss through the original
/// body network while retaining each body's own reversible entropy exchange.
pub mod resolved;

#[derive(Clone)]
enum CaloricHeatRouting {
    Prescribed(HeatPartition),
    Resolved(crate::friction_heat::FrictionHeatRoute),
}

/// Derivatives of the two complete-patch stiffnesses with respect to ONE
/// surface temperature, holding the other temperature and elastic strains fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StiffnessDerivative {
    pub tangential_n_per_m_k: f64,
    pub torsional_nm_per_rad_k: f64,
}

/// Immutable declared free-energy model and its closed rectangular T domain.
/// The reference law's source/model identity names this whole affine model.
/// Supplying an identity does not establish that its coefficients are measured.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineContactStiffness {
    reference: PartialSlipLaw,
    reference_temperature_k: [f64; 2],
    derivatives: [StiffnessDerivative; 2],
    temperature_domain_k: [[f64; 2]; 2],
}

impl AffineContactStiffness {
    pub fn new(
        reference: PartialSlipLaw,
        reference_temperature_k: [f64; 2],
        derivatives: [StiffnessDerivative; 2],
        temperature_domain_k: [[f64; 2]; 2],
    ) -> Result<Self, CaloricContactError> {
        for i in 0..2 {
            let [lo, hi] = temperature_domain_k[i];
            if !lo.is_finite()
                || lo <= 0.0
                || !hi.is_finite()
                || hi < lo
                || !reference_temperature_k[i].is_finite()
                || reference_temperature_k[i] < lo
                || reference_temperature_k[i] > hi
                || !derivatives[i].tangential_n_per_m_k.is_finite()
                || !derivatives[i].torsional_nm_per_rad_k.is_finite()
            {
                return Err(CaloricContactError::Invalid(
                    "invalid affine contact temperature domain",
                ));
            }
        }
        let model = Self {
            reference,
            reference_temperature_k,
            derivatives,
            temperature_domain_k,
        };
        // Affine stiffness is positive throughout the rectangle iff it is
        // positive at all four corners. Native law admission checks both modes.
        for a in temperature_domain_k[0] {
            for b in temperature_domain_k[1] {
                model.law_at([a, b])?;
            }
        }
        Ok(model)
    }

    pub const fn reference_law(&self) -> &PartialSlipLaw {
        &self.reference
    }
    pub const fn derivatives(&self) -> &[StiffnessDerivative; 2] {
        &self.derivatives
    }
    pub const fn temperature_domain_k(&self) -> &[[f64; 2]; 2] {
        &self.temperature_domain_k
    }

    pub fn law_at<E>(
        &self,
        temperature_k: [f64; 2],
    ) -> Result<PartialSlipLaw, CaloricContactError<E>> {
        let mut p = self.reference.parameters();
        for (i, &t) in temperature_k.iter().enumerate() {
            if !t.is_finite()
                || t < self.temperature_domain_k[i][0]
                || t > self.temperature_domain_k[i][1]
            {
                return Err(CaloricContactError::Invalid(
                    "contact temperature outside declared free-energy domain",
                ));
            }
            let dt = finite(t - self.reference_temperature_k[i])?;
            p.tangential_stiffness_n_per_m = finite(
                p.tangential_stiffness_n_per_m
                    + finite(self.derivatives[i].tangential_n_per_m_k * dt)?,
            )?;
            p.torsional_stiffness_nm_per_rad = finite(
                p.torsional_stiffness_nm_per_rad
                    + finite(self.derivatives[i].torsional_nm_per_rad_k * dt)?,
            )?;
        }
        PartialSlipLaw::new(self.reference.model_id(), self.reference.source_id(), p)
            .map_err(CaloricContactError::Contact)
    }

    fn thermodynamics<E>(
        &self,
        state: &PartialSlipState,
        temperature_k: [f64; 2],
    ) -> Result<ContactThermodynamics, CaloricContactError<E>> {
        let p = self.law_at::<E>(temperature_k)?.parameters();
        let [translation, twist] = quadratic_coordinates(state, p.partial_slip_hardening_fraction)?;
        let free_energy_j = finite(
            finite(p.tangential_stiffness_n_per_m * translation)?
                + finite(p.torsional_stiffness_nm_per_rad * twist)?,
        )?;
        let mut entropy_j_k = [0.0; 2];
        let mut internal_tangent = self.reference.parameters().tangential_stiffness_n_per_m;
        let mut internal_twist = self.reference.parameters().torsional_stiffness_nm_per_rad;
        for i in 0..2 {
            entropy_j_k[i] = finite(
                -finite(self.derivatives[i].tangential_n_per_m_k * translation)?
                    - finite(self.derivatives[i].torsional_nm_per_rad_k * twist)?,
            )?;
            // k - sum(Ti * dki/dTi) is temperature independent for this
            // affine model. Evaluate it from the reference, so a large psi and
            // opposing T*s are not subtracted at every nonlinear trial.
            internal_tangent = finite(
                internal_tangent
                    - finite(
                        self.reference_temperature_k[i] * self.derivatives[i].tangential_n_per_m_k,
                    )?,
            )?;
            internal_twist = finite(
                internal_twist
                    - finite(
                        self.reference_temperature_k[i]
                            * self.derivatives[i].torsional_nm_per_rad_k,
                    )?,
            )?;
        }
        Ok(ContactThermodynamics {
            free_energy_j,
            entropy_j_k,
            internal_energy_j: finite(
                finite(internal_tangent * translation)? + finite(internal_twist * twist)?,
            )?,
        })
    }
}

/// Mechanical contribution only; add the original finite caloric storage.
/// These entropies do not assert total entropy production or a heat-transfer law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactThermodynamics {
    pub free_energy_j: f64,
    pub entropy_j_k: [f64; 2],
    pub internal_energy_j: f64,
}

/// One admitted model plus COMPLETE accepted native contact history and clock.
#[derive(Clone)]
pub struct CaloricPartialSlip<'a> {
    model: AffineContactStiffness,
    checkpoint: PartialSlipCheckpoint,
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    cells: [CellId; 2],
    routing: CaloricHeatRouting,
    temperature_k: [f64; 2],
    thermodynamics: ContactThermodynamics,
    time_s: f64,
    // Reuse only the original immutable pinned binding/query owner. Native
    // partial-slip history above remains the sole mechanical history owner.
    material_friction: Option<MaterialSlidingContact<'a>>,
}

/// Explicit nonlinear work and independently enforced cell/global energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricContactSolve {
    pub maximum_cells: usize,
    pub maximum_contact_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Finite-difference perturbation of CELL energy [J], not temperature.
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum CaloricContactError<E = Infallible> {
    Invalid(&'static str),
    Contact(PartialSlipError),
    Material(SlidingError<E>),
    Phase(PhaseStateError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    Cancelled,
}
impl<E> From<DepositionError> for CaloricContactError<E> {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for CaloricContactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "partial-slip caloric coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CaloricContactError<E> {}

pub struct CaloricContactReport<S> {
    pub contact: PartialSlipStep,
    /// Exact pinned friction claims at the accepted two-temperature endpoint.
    /// Independent of the declared affine stiffness/free-energy source identity.
    pub material: Option<SlidingMaterialSample>,
    pub thermodynamics: ContactThermodynamics,
    pub contact_internal_energy_change_j: f64,
    /// Signed reversible heat into each surface; may heat or cool.
    pub reversible_heat_j: [f64; 2],
    /// Nonnegative physical microslip allocation, excluding numerical loss.
    pub physical_heat_j: [f64; 2],
    pub deposited_heat_j: [f64; 2],
    pub thermal: EnergyAudit,
    pub external_heat_j: f64,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub contact_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Delta(H_caloric)+Delta(U_contact)+D_numerical-W_wrench-Q_external.
    pub total_energy_residual_j: f64,
}

pub struct CaloricContactCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: CaloricPartialSlip<'a>,
    report: CaloricContactReport<S>,
}

struct MechanicalEndpoint {
    contact: PartialSlipStep,
    material: Option<SlidingMaterialSample>,
    thermodynamics: ContactThermodynamics,
    temperature_k: [f64; 2],
    internal_change: f64,
    reversible_heat_j: [f64; 2],
}
impl<'a, 't, S> CaloricContactCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &CaloricPartialSlip<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &CaloricContactReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CaloricPartialSlip<'a>,
        CaloricContactReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a> CaloricPartialSlip<'a> {
    /// Bind the ORIGINAL caloric objects and accepted native checkpoint.
    /// This does not reset loaded strains, kinetic branches or frame transport.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        model: AffineContactStiffness,
        checkpoint: PartialSlipCheckpoint,
        curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
        cells: [CellId; 2],
        partition: HeatPartition,
    ) -> Result<Self, CaloricContactError> {
        if cells[0] == cells[1] || partition.other() != 0.0 {
            return Err(CaloricContactError::Invalid(
                "two distinct finite caloric destinations required",
            ));
        }
        Self::from_checkpoint(
            thermal,
            model,
            checkpoint,
            curves,
            cells,
            CaloricHeatRouting::Prescribed(partition),
        )
    }

    fn from_checkpoint(
        thermal: &DepositionSimulation<'_>,
        model: AffineContactStiffness,
        checkpoint: PartialSlipCheckpoint,
        curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
        cells: [CellId; 2],
        routing: CaloricHeatRouting,
    ) -> Result<Self, CaloricContactError> {
        if cells[0] == cells[1] {
            return Err(CaloricContactError::Invalid(
                "two distinct finite caloric destinations required",
            ));
        }
        let temperature_k = temperatures(thermal, curves, cells)?;
        if model.law_at::<Infallible>(temperature_k)? != *checkpoint.law() {
            return Err(CaloricContactError::Invalid(
                "checkpoint differs from affine law at actual initial temperatures",
            ));
        }
        checkpoint
            .law()
            .restore_checkpoint(
                checkpoint.normal_patch(),
                checkpoint.interface(),
                &checkpoint,
            )
            .map_err(CaloricContactError::Contact)?;
        let thermodynamics = model.thermodynamics(checkpoint.state(), temperature_k)?;
        Ok(Self {
            model,
            checkpoint,
            curves,
            cells,
            routing,
            temperature_k,
            thermodynamics,
            time_s: thermal.time_s(),
            material_friction: None,
        })
    }

    /// Attach the ORIGINAL pinned Coulomb source resolver without resetting
    /// loaded strains, kinetic branches or the affine free-energy model.
    ///
    /// Ordered interface/history/source identities, caloric objects and cells
    /// must match. Prescribed routing also requires the original heat split;
    /// resolved routing derives its split from the actual interface conductance
    /// and does not use the binding's legacy partition metadata. The source's initial coefficients
    /// must match the accepted native law; changing them requires a real step.
    /// Rebinding an existing source is refused. The stiffness model's identity
    /// and the friction claim identities remain distinct.
    ///
    /// This ingress admits translational contact only: no retained twist,
    /// relative torsional motion, accumulated-distance law or effusivity model.
    /// The source speed is the actual norm of creepage times rolling speed,
    /// not rolling speed alone or an invented torsional equivalent. Pressure
    /// uses the held normal load and the binding's explicit nominal area.
    /// Reported travel is actual prescribed translational arc length since this
    /// attachment; it is diagnostic and cannot supply a prehistory-dependent
    /// source axis. Its counter advances only with an accepted candidate.
    pub fn with_material_friction(
        mut self,
        thermal: &DepositionSimulation<'_>,
        data: SlidingInterface<'a>,
        initial_speed_mps: f64,
        authority: InputAuthority,
    ) -> Result<Self, CaloricContactError> {
        if self.material_friction.is_some()
            || thermal.time_s() != self.time_s
            || temperatures(thermal, self.curves, self.cells)? != self.temperature_k
            || data.cells != self.cells
            || (0..2).any(|i| !std::ptr::eq(data.caloric[i], self.curves[i]))
            || data.accumulated_distance_axis.is_some()
            || !matches!(&data.partition, SlidingHeatPartition::Prescribed(p)
            if p.other() == 0.0 && match &self.routing {
                CaloricHeatRouting::Prescribed(original) => p == original,
                CaloricHeatRouting::Resolved(_) => true,
            })
            || self.checkpoint.state().region_elastic_twist_rad() != [0.0; 2]
        {
            return Err(CaloricContactError::Invalid(
                "material contact needs original ordered cells/curves, compatible heat ownership and translational history",
            ));
        }
        let source = MaterialSlidingContact::new(
            thermal,
            data,
            self.checkpoint.normal_patch().normal_load_n(),
            initial_speed_mps,
            authority,
        )
        .map_err(CaloricContactError::Material)?;
        let bound = source.interface().interface();
        let original = self.checkpoint.interface();
        if original.ordered_interface_id() != bound.ordered_system_id()
            || original.history_id() != bound.history_id()
            || original.source_id() != bound.provenance().source_id()
        {
            return Err(CaloricContactError::Invalid(
                "material friction differs from retained ordered contact identity",
            ));
        }
        let sample = source
            .sample(
                thermal,
                self.checkpoint.normal_patch().normal_load_n(),
                initial_speed_mps,
                0.0,
            )
            .map_err(CaloricContactError::Material)?;
        let [static_mu, kinetic_mu] = coulomb(&sample)?;
        let p = self.checkpoint.law().parameters();
        if p.static_mu != static_mu || p.kinetic_mu != kinetic_mu {
            return Err(CaloricContactError::Invalid(
                "material friction differs from accepted initial coefficients",
            ));
        }
        self.material_friction = Some(source);
        Ok(self)
    }

    /// Immutable original binding including exact pinned friction claims.
    pub fn material_friction_binding(&self) -> Option<&SlidingInterface<'a>> {
        self.material_friction
            .as_ref()
            .map(MaterialSlidingContact::binding)
    }

    pub const fn checkpoint(&self) -> &PartialSlipCheckpoint {
        &self.checkpoint
    }
    pub const fn model(&self) -> &AffineContactStiffness {
        &self.model
    }
    pub const fn thermodynamics(&self) -> &ContactThermodynamics {
        &self.thermodynamics
    }
    pub const fn temperature_k(&self) -> [f64; 2] {
        self.temperature_k
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    #[allow(clippy::too_many_arguments)]
    fn mechanical_endpoint<E>(
        &self,
        state: &DepositionSimulation<'_>,
        frame: TangentFrame,
        kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        controls: CaloricContactSolve,
    ) -> Result<MechanicalEndpoint, CaloricContactError<E>> {
        let dt = kinematics.dt_s;
        let temperature_k = temperatures(state, self.curves, self.cells)?;
        let mut law = self.model.law_at(temperature_k)?;
        let material = self
            .material_friction
            .as_ref()
            .map(|source| {
                let vx = finite(kinematics.creepage[0] * kinematics.rolling_speed_mps)?;
                let vy = finite(kinematics.creepage[1] * kinematics.rolling_speed_mps)?;
                let speed = finite(det::hypot(vx, vy))?;
                let distance = finite(speed * dt)?;
                let accumulated = finite(source.distance_m + distance)?;
                if distance > 0.0 && accumulated == source.distance_m {
                    return Err(CaloricContactError::Invalid(
                        "material contact travel is unrepresentable",
                    ));
                }
                source
                    .sample(
                        state,
                        self.checkpoint.normal_patch().normal_load_n(),
                        speed,
                        accumulated,
                    )
                    .map_err(CaloricContactError::Material)
            })
            .transpose()?;
        if let Some(sample) = &material {
            let [static_mu, kinetic_mu] = coulomb(sample)?;
            let mut parameters = law.parameters();
            parameters.static_mu = static_mu;
            parameters.kinetic_mu = kinetic_mu;
            law = PartialSlipLaw::new(law.model_id(), law.source_id(), parameters)
                .map_err(CaloricContactError::Contact)?;
        }
        let contact = law
            .advance_with_coefficients(
                self.checkpoint.normal_patch(),
                self.checkpoint.interface(),
                frame,
                kinematics,
                ownership,
                &self.checkpoint,
            )
            .map_err(CaloricContactError::Contact)?;
        let thermodynamics = self
            .model
            .thermodynamics(&contact.next_state, temperature_k)?;
        let internal_change =
            finite(thermodynamics.internal_energy_j - self.thermodynamics.internal_energy_j)?;
        if contact.dissipation.heat_j < 0.0
            || contact.generalized_work.numerical_dissipation_j < 0.0
            || contact.generalized_work.closure_residual_j.abs() > controls.total_energy_tolerance_j
            || finite(thermodynamics.free_energy_j - contact.stored_energy_j)?.abs()
                > controls.total_energy_tolerance_j
        {
            return Err(CaloricContactError::Invalid(
                "native contact free energy or work failed closure",
            ));
        }
        let mut reversible_heat_j = [0.0; 2];
        for i in 0..2 {
            // Retain each body's original -Ti*dsi exchange. It is neither
            // inter-body conduction nor irreversible microslip heating.
            reversible_heat_j[i] = finite(
                -temperature_k[i]
                    * finite(thermodynamics.entropy_j_k[i] - self.thermodynamics.entropy_j_k[i])?,
            )?;
        }
        let defect = finite(
            internal_change
                + contact.dissipation.heat_j
                + reversible_heat_j[0]
                + reversible_heat_j[1]
                + contact.generalized_work.numerical_dissipation_j
                - contact.generalized_work.work_into_interface_j,
        )?;
        if defect.abs() > controls.total_energy_tolerance_j {
            return Err(CaloricContactError::Invalid(
                "contact internal energy and reversible caloric exchange do not close",
            ));
        }
        Ok(MechanicalEndpoint {
            contact,
            material,
            thermodynamics,
            temperature_k,
            internal_change,
            reversible_heat_j,
        })
    }
    /// Solve prescribed relative motion and both actual surface temperatures in
    /// one original finite thermal network. The normal patch/load remains fixed.
    /// Callback supplies only OTHER sources and conductances, evaluated at this
    /// same proposed endpoint. Every trial is immutable; publish both returned
    /// owners together. Cut intervals at original source/topology events.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        frame: TangentFrame,
        mut kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        controls: CaloricContactSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CaloricContactCandidate<'a, 't, S>, CaloricContactError<E>> {
        let CaloricHeatRouting::Prescribed(partition) = &self.routing else {
            return Err(CaloricContactError::Invalid(
                "use the resolved-interface step for this contact heat owner",
            ));
        };
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        controls
            .newton
            .validate()
            .map_err(CaloricContactError::Newton)?;
        if thermal.time_s() != self.time_s
            || temperatures(thermal, self.curves, self.cells)? != self.temperature_k
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || guess_j.len() != thermal.cells().len()
            || controls.maximum_contact_trials == 0
            || controls.maximum_newton_iterations == 0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CaloricContactError::Invalid(
                "incompatible accepted state, thermal unknowns or controls",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, poll)?;
        let dt = step.duration_s();
        if kinematics.dt_s != dt {
            return Err(CaloricContactError::Invalid(
                "contact and caloric physical intervals must match",
            ));
        }
        if self.material_friction.is_some() && kinematics.torsional_spin_rad_per_s != 0.0 {
            return Err(CaloricContactError::Invalid(
                "material friction has no torsional-speed mapping",
            ));
        }
        // Keep the actual physical interval identical in all native calls.
        kinematics.dt_s = dt;
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |state| {
                    poll()?;
                    if trials.get() >= controls.maximum_contact_trials {
                        return Err(CaloricContactError::TrialBudget);
                    }
                    trials.set(trials.get() + 1);
                    let MechanicalEndpoint {
                        contact,
                        material,
                        thermodynamics,
                        temperature_k,
                        internal_change,
                        reversible_heat_j,
                    } = self.mechanical_endpoint(state, frame, kinematics, ownership, controls)?;
                    let physical_heat_j = [
                        finite(partition.surface_a() * contact.dissipation.heat_j)?,
                        finite(partition.surface_b() * contact.dissipation.heat_j)?,
                    ];
                    let deposited_heat_j = [
                        finite(physical_heat_j[0] + reversible_heat_j[0])?,
                        finite(physical_heat_j[1] + reversible_heat_j[1])?,
                    ];
                    let coupled_defect = finite(
                        internal_change
                            + deposited_heat_j[0]
                            + deposited_heat_j[1]
                            + contact.generalized_work.numerical_dissipation_j
                            - contact.generalized_work.work_into_interface_j,
                    )?;
                    if coupled_defect.abs() > controls.total_energy_tolerance_j {
                        return Err(CaloricContactError::Invalid(
                            "contact internal energy and caloric exchange do not close",
                        ));
                    }
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(state).map_err(CaloricContactError::Source)?;
                    let external_heat_j = finite(
                        dt * state.external_step_power_w(&sources_w, &conductances_w_k, poll)?,
                    )?;
                    for (id, heat) in self.cells.into_iter().zip(deposited_heat_j) {
                        let power = finite(heat / dt)?;
                        if heat != 0.0 && power == 0.0 {
                            return Err(CaloricContactError::Invalid(
                                "contact heat vanished in power conversion",
                            ));
                        }
                        let total = sources_w.entry(id).or_insert(0.0);
                        *total = finite(*total + power)?;
                    }
                    Ok(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: (
                            contact,
                            material,
                            thermodynamics,
                            temperature_k,
                            internal_change,
                            reversible_heat_j,
                            physical_heat_j,
                            deposited_heat_j,
                            external_heat_j,
                            payload,
                        ),
                    })
                },
                poll,
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec()),
        )
        .map_err(difference_error)?;
        let construction = NewtonKrylovState::new(&problem, guess_j.to_vec(), controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        let mut solver = construction.map_err(CaloricContactError::Newton)?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(CaloricContactError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(CaloricContactError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            next,
            (
                contact,
                material,
                thermodynamics,
                temperature_k,
                contact_internal_energy_change_j,
                reversible_heat_j,
                physical_heat_j,
                deposited_heat_j,
                external_heat_j,
                source,
            ),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            poll,
        )?;
        let total_energy_residual_j = finite(
            thermal_audit.stored_energy_change_j
                + contact_internal_energy_change_j
                + contact.generalized_work.numerical_dissipation_j
                - contact.generalized_work.work_into_interface_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(CaloricContactError::Invalid(
                "finite caloric and mechanical internal energy do not balance",
            ));
        }
        let mut next_contact = self.clone();
        next_contact.checkpoint = contact.checkpoint.clone();
        next_contact.temperature_k = temperature_k;
        next_contact.thermodynamics = thermodynamics;
        next_contact.time_s = target_s;
        if let (Some(binding), Some(sample)) = (&mut next_contact.material_friction, &material) {
            binding.distance_m = sample.accumulated_distance_m;
            binding.time_s = target_s;
        }
        poll()?;
        Ok(CaloricContactCandidate {
            thermal: next,
            contact: next_contact,
            report: CaloricContactReport {
                contact,
                material,
                thermodynamics,
                contact_internal_energy_change_j,
                reversible_heat_j,
                physical_heat_j,
                deposited_heat_j,
                thermal: thermal_audit,
                external_heat_j,
                source,
                nonlinear: result.report,
                contact_trials: trials.get(),
                maximum_cell_residual_j,
                total_energy_residual_j,
            },
        })
    }
}

fn temperatures<E>(
    thermal: &DepositionSimulation<'_>,
    curves: [&EquilibriumEnthalpyPhaseCurve; 2],
    cells: [CellId; 2],
) -> Result<[f64; 2], CaloricContactError<E>> {
    let mut temperatures = [0.0; 2];
    for i in 0..2 {
        let cell = thermal
            .cells()
            .get(&cells[i])
            .ok_or(CaloricContactError::Invalid(
                "inactive contact caloric cell",
            ))?;
        if !thermal.uses_material_law(cell.material_id, curves[i]) {
            return Err(CaloricContactError::Invalid(
                "contact must retain original caloric objects",
            ));
        }
        let phase = curves[i]
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(CaloricContactError::Phase)?;
        if phase.liquid_mass_fraction() != 0.0 {
            return Err(CaloricContactError::Invalid(
                "dry loaded contact refuses liquid caloric states",
            ));
        }
        temperatures[i] = phase.temperature_k();
    }
    Ok(temperatures)
}

fn coulomb<E>(sample: &SlidingMaterialSample) -> Result<[f64; 2], CaloricContactError<E>> {
    match sample.law {
        FrictionLaw::Coulomb {
            static_mu,
            kinetic_mu,
        } => Ok([static_mu, kinetic_mu]),
        _ => Err(CaloricContactError::Invalid(
            "source must retain the original resolved Coulomb law",
        )),
    }
}

fn quadratic_coordinates<E>(
    state: &PartialSlipState,
    hard: f64,
) -> Result<[f64; 2], CaloricContactError<E>> {
    let shear = state.region_elastic_displacement_m();
    let twist = state.region_elastic_twist_rad();
    let mut translation = 0.0;
    let mut torsion = 0.0;
    for (i, weight) in [1.0 - hard, hard].into_iter().enumerate() {
        translation = finite(
            translation
                + finite(
                    0.5 * weight
                        * finite(
                            finite(shear[i][0] * shear[i][0])? + finite(shear[i][1] * shear[i][1])?,
                        )?,
                )?,
        )?;
        torsion = finite(torsion + finite(0.5 * weight * finite(twist[i] * twist[i])?)?)?;
    }
    Ok([translation, torsion])
}
fn finite<E>(value: f64) -> Result<f64, CaloricContactError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CaloricContactError::Invalid(
            "unrepresentable contact energy or coefficient",
        ))
    }
}
fn difference_error<E>(e: FiniteDifferenceError<CaloricContactError<E>>) -> CaloricContactError<E> {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => CaloricContactError::Invalid(s),
    }
}
