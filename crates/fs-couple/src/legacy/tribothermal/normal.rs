//! Reciprocal finite-caloric Hertz/Hunt--Crossley normal contact.
//!
//! The original point-contact owner supplies geometry, force, pressure, potential
//! and applicability. The reduced-modulus model supplies its original
//! two-temperature free-energy derivatives at fixed geometry. The included
//! affine model is explicitly declared, not a measured material calibration.
//!
//! Indentation follows one prescribed linear segment. Its mean force and
//! Hunt--Crossley work integrate the native Hertz potential exactly at the
//! proposed endpoint temperatures. The finite thermal network receives physical
//! damping heat plus reversible -Ci-d(Ti*si), where si=-dpsi/dTi,
//! U=psi+sum(Ti*si), and Ci is each surface's coefficient free-energy change.
//! No coefficient/free-energy change is silently converted to heat. Thermal
//! composition remains first order and needs timestep refinement.
//!
//! Only nonnegative indentation, sphere/elliptic half-space geometry and constant
//! declared Hunt--Crossley damping are admitted. No flight/impact event location,
//! changing curvature, thermal expansion, adhesion, flash temperature or
//! simultaneous tangential/body solve is claimed. Every trial restarts from the
//! accepted state; publish the returned mechanical and thermal owners together.

use fs_contact::normal_patch::{
    NormalPatchError, NormalPatchLaw, NormalPatchReceipt, NormalPatchRequest,
    PointNormalPatchReceipt,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::HeatPartition;
use std::cell::{Cell, RefCell};
use std::convert::Infallible;

/// Source or declared thermodynamic reduced-modulus response at one ordered
/// pair of surface temperatures. Derivatives hold geometry fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalModulusState {
    pub modulus_pa: f64,
    pub gradient_pa_k: [f64; 2],
    pub hessian_pa_k2: [[f64; 2]; 2],
}

/// Immutable constitutive model. Implementations must return derivatives of
/// the SAME free-energy modulus and retain their original source domains.
/// The solver also queries the two mixed-temperature corners of an interval
/// for a symmetric coordinate-path split of reversible thermal exchange.
pub trait NormalModulus: Clone {
    fn reference_request(&self) -> &NormalPatchRequest;
    fn evaluate(&self, temperature_k: [f64; 2]) -> Result<NormalModulusState, NormalThermalError>;
    /// A source-backed implementation can require its exact original caloric
    /// objects. Bare declared models retain the objects supplied at binding.
    fn expected_caloric_curves(&self) -> Option<[&EquilibriumEnthalpyPhaseCurve; 2]> {
        None
    }
}

/// Immutable declared reduced-modulus law, retaining the original normal
/// request's interface, geometry, source identity and applicability envelope.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineNormalModulus {
    reference: NormalPatchRequest,
    reference_temperature_k: [f64; 2],
    derivatives_pa_k: [f64; 2],
    temperature_domain_k: [[f64; 2]; 2],
    internal_modulus_pa: f64,
}

impl AffineNormalModulus {
    pub fn new(
        reference: NormalPatchRequest,
        reference_temperature_k: [f64; 2],
        derivatives_pa_k: [f64; 2],
        temperature_domain_k: [[f64; 2]; 2],
    ) -> Result<Self, NormalThermalError> {
        let (reference_modulus, _) = modulus_and_damping(&reference.law)?;
        let mut internal_modulus_pa = reference_modulus;
        for i in 0..2 {
            let [lo, hi] = temperature_domain_k[i];
            if !lo.is_finite()
                || lo <= 0.0
                || !hi.is_finite()
                || hi < lo
                || !reference_temperature_k[i].is_finite()
                || reference_temperature_k[i] < lo
                || reference_temperature_k[i] > hi
                || !derivatives_pa_k[i].is_finite()
                || lo < reference.limits.min_temperature_k
                || hi > reference.limits.max_temperature_k
            {
                return Err(NormalThermalError::Invalid(
                    "invalid normal-contact temperature domain",
                ));
            }
            internal_modulus_pa = finite(
                internal_modulus_pa - finite(reference_temperature_k[i] * derivatives_pa_k[i])?,
            )?;
        }
        // This bounded rung admits positive elastic internal storage as well as
        // positive tangent stiffness; caloric heat capacities remain unchanged.
        if internal_modulus_pa <= 0.0 {
            return Err(NormalThermalError::Invalid(
                "nonpositive contact internal-energy modulus",
            ));
        }
        let model = Self {
            reference,
            reference_temperature_k,
            derivatives_pa_k,
            temperature_domain_k,
            internal_modulus_pa,
        };
        for a in temperature_domain_k[0] {
            for b in temperature_domain_k[1] {
                model.modulus_at::<Infallible>([a, b])?;
            }
        }
        sample::<_, Infallible>(
            &model,
            model.reference.indentation_m,
            0.0,
            reference_temperature_k,
        )?;
        Ok(model)
    }

    pub const fn reference_request(&self) -> &NormalPatchRequest {
        &self.reference
    }
    pub const fn derivatives_pa_k(&self) -> [f64; 2] {
        self.derivatives_pa_k
    }
    pub const fn temperature_domain_k(&self) -> [[f64; 2]; 2] {
        self.temperature_domain_k
    }

    pub fn modulus_at<E>(&self, temperature_k: [f64; 2]) -> Result<f64, NormalThermalError<E>> {
        let (mut modulus, _) = modulus_and_damping(&self.reference.law)?;
        for i in 0..2 {
            let t = temperature_k[i];
            if !t.is_finite()
                || t < self.temperature_domain_k[i][0]
                || t > self.temperature_domain_k[i][1]
            {
                return Err(NormalThermalError::Invalid(
                    "normal temperature outside declared domain",
                ));
            }
            modulus = finite(
                modulus
                    + finite(
                        self.derivatives_pa_k[i] * finite(t - self.reference_temperature_k[i])?,
                    )?,
            )?;
        }
        if modulus <= 0.0 {
            return Err(NormalThermalError::Invalid("nonpositive normal modulus"));
        }
        Ok(modulus)
    }
}
impl NormalModulus for AffineNormalModulus {
    fn reference_request(&self) -> &NormalPatchRequest {
        &self.reference
    }
    fn evaluate(&self, temperature_k: [f64; 2]) -> Result<NormalModulusState, NormalThermalError> {
        Ok(NormalModulusState {
            modulus_pa: self.modulus_at(temperature_k)?,
            gradient_pa_k: self.derivatives_pa_k,
            hessian_pa_k2: [[0.0; 2]; 2],
        })
    }
}

fn modulus_state<E>(
    model: &impl NormalModulus,
    temperature_k: [f64; 2],
) -> Result<NormalModulusState, NormalThermalError<E>> {
    let state = model.evaluate(temperature_k).map_err(model_error)?;
    let limits = model.reference_request().limits;
    if temperature_k.iter().any(|t| {
        !t.is_finite()
            || *t <= 0.0
            || *t < limits.min_temperature_k
            || *t > limits.max_temperature_k
    }) || !state.modulus_pa.is_finite()
        || state.modulus_pa <= 0.0
        || state
            .gradient_pa_k
            .iter()
            .chain(state.hessian_pa_k2.iter().flatten())
            .any(|v| !v.is_finite())
        || state.hessian_pa_k2[0][1] != state.hessian_pa_k2[1][0]
    {
        return Err(NormalThermalError::Invalid(
            "invalid normal free-energy response",
        ));
    }
    Ok(state)
}

fn sample<M: NormalModulus, E>(
    model: &M,
    indentation_m: f64,
    rate_mps: f64,
    temperature_k: [f64; 2],
) -> Result<(PointNormalPatchReceipt, NormalModulusState), NormalThermalError<E>> {
    let modulus = modulus_state(model, temperature_k)?;
    let mut request = model.reference_request().clone();
    match &mut request.law {
        NormalPatchLaw::HertzSpherePlane {
            reduced_modulus_pa, ..
        }
        | NormalPatchLaw::HuntCrossleySphere {
            reduced_modulus_pa, ..
        }
        | NormalPatchLaw::HertzEllipticParaboloid {
            reduced_modulus_pa, ..
        }
        | NormalPatchLaw::HuntCrossleyEllipticParaboloid {
            reduced_modulus_pa, ..
        } => *reduced_modulus_pa = modulus.modulus_pa,
        _ => {
            return Err(NormalThermalError::Invalid(
                "normal law is not a Hertz point potential",
            ));
        }
    }
    request.indentation_m = indentation_m;
    request.indentation_rate_m_per_s = rate_mps;
    // Instantaneous native sample, not endpoint-power quadrature.
    request.step_s = 0.0;
    // Both temperatures were admitted independently above. The mean only
    // labels this constant-damping request; it is not a material source query.
    request.applicability.temperature_k = finite(0.5 * temperature_k[0] + 0.5 * temperature_k[1])?;
    match request.evaluate().map_err(NormalThermalError::Normal)? {
        NormalPatchReceipt::Point(point) => Ok((point, modulus)),
        _ => Err(NormalThermalError::Invalid("point normal receipt required")),
    }
}

fn potentials<E>(
    point: &PointNormalPatchReceipt,
    modulus: NormalModulusState,
    temperature_k: [f64; 2],
) -> Result<NormalThermodynamics, NormalThermalError<E>> {
    let basis = finite(point.reversible_energy_j / modulus.modulus_pa)?;
    let entropy_j_k = [
        finite(-modulus.gradient_pa_k[0] * basis)?,
        finite(-modulus.gradient_pa_k[1] * basis)?,
    ];
    Ok(NormalThermodynamics {
        free_energy_j: point.reversible_energy_j,
        entropy_j_k,
        internal_energy_j: finite(
            basis
                * finite(
                    modulus.modulus_pa
                        - finite(temperature_k[0] * modulus.gradient_pa_k[0])?
                        - finite(temperature_k[1] * modulus.gradient_pa_k[1])?,
                )?,
        )?,
    })
}

/// Sufficient local thermodynamic stability: the full two-temperature
/// Helmholtz Hessian must be negative definite. Use the original caloric
/// chart's one-sided derivative at knots; no large latent Cp is fabricated.
fn capacity<E>(
    thermal: &DepositionSimulation<'_>,
    curves: [&EquilibriumEnthalpyPhaseCurve; 2],
    cells: [CellId; 2],
    point: &PointNormalPatchReceipt,
    modulus: NormalModulusState,
    temperature_k: [f64; 2],
) -> Result<(), NormalThermalError<E>> {
    let basis = finite(point.reversible_energy_j / modulus.modulus_pa)?;
    let mut diagonal = [0.0; 2];
    for i in 0..2 {
        let cell = &thermal.cells()[&cells[i]];
        let slope = curves[i]
            .temperature_derivative_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(NormalThermalError::Phase)?;
        if !slope.is_finite() || slope <= 0.0 {
            return Err(NormalThermalError::Invalid(
                "normal contact needs positive finite caloric capacity",
            ));
        }
        diagonal[i] = finite(
            finite(cell.mass_kg / slope)? / temperature_k[i]
                - finite(basis * modulus.hessian_pa_k2[i][i])?,
        )?;
        if diagonal[i] <= 0.0 {
            return Err(NormalThermalError::Invalid(
                "nonpositive total normal-contact thermal capacity",
            ));
        }
    }
    let off = finite(basis * modulus.hessian_pa_k2[0][1])?;
    let correlation = finite(finite(off.abs() / diagonal[0].sqrt())? / diagonal[1].sqrt())?;
    if correlation >= 1.0 {
        return Err(NormalThermalError::Invalid(
            "normal-contact thermal capacity is not positive definite",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalThermodynamics {
    pub free_energy_j: f64,
    pub entropy_j_k: [f64; 2],
    pub internal_energy_j: f64,
}

/// Exact line-segment work for the admitted Hertz potential at endpoint
/// temperatures. Mean force, not endpoint force, is conjugate to displacement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalContactWork {
    pub indentation_change_m: f64,
    pub indentation_rate_mps: f64,
    pub mean_elastic_force_n: f64,
    pub mean_normal_force_n: f64,
    pub elastic_work_j: f64,
    pub physical_dissipation_j: f64,
    pub work_into_contact_j: f64,
    /// Change of psi due only to temperature at the OLD indentation.
    /// This is neither external mechanical work nor dissipated heat.
    pub coefficient_free_energy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub closure_residual_j: f64,
}

#[derive(Clone)]
pub struct CaloricNormalContact<'a, M: NormalModulus = AffineNormalModulus> {
    model: M,
    modulus: NormalModulusState,
    mass_kg: [f64; 2],
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    cells: [CellId; 2],
    partition: HeatPartition,
    endpoint: PointNormalPatchReceipt,
    temperature_k: [f64; 2],
    thermodynamics: NormalThermodynamics,
    time_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalThermalSolve {
    pub maximum_cells: usize,
    pub maximum_contact_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Perturbation of cell ENERGY [J], not temperature.
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum NormalThermalError<E = Infallible> {
    Invalid(&'static str),
    Normal(NormalPatchError),
    Phase(PhaseStateError),
    Thermal(DepositionError),
    Source(E),
    /// Original structured source error, retained for downcast/diagnosis.
    MaterialModel(Box<dyn std::error::Error + Send + Sync>),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    Cancelled,
}
impl<E> From<DepositionError> for NormalThermalError<E> {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for NormalThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-caloric normal contact: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for NormalThermalError<E> {}

pub struct NormalThermalReport<S> {
    pub modulus: NormalModulusState,
    /// Native endpoint sample with step_s=0: irreversible_work_j is zero.
    /// physical_dissipation_j in work is the actual segment-integrated heat.
    pub endpoint: PointNormalPatchReceipt,
    pub thermodynamics: NormalThermodynamics,
    pub work: NormalContactWork,
    pub reversible_heat_j: [f64; 2],
    pub physical_heat_j: [f64; 2],
    pub deposited_heat_j: [f64; 2],
    pub thermal: EnergyAudit,
    pub external_heat_j: f64,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub contact_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Delta(H_caloric)+Delta(U_normal)-W_normal-Q_external.
    pub total_energy_residual_j: f64,
}

pub struct NormalThermalCandidate<'a, 't, S, M: NormalModulus = AffineNormalModulus> {
    thermal: DepositionSimulation<'t>,
    contact: CaloricNormalContact<'a, M>,
    report: NormalThermalReport<S>,
}
impl<'a, 't, S, M: NormalModulus> NormalThermalCandidate<'a, 't, S, M> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &CaloricNormalContact<'a, M> {
        &self.contact
    }
    pub const fn report(&self) -> &NormalThermalReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CaloricNormalContact<'a, M>,
        NormalThermalReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a, M: NormalModulus> CaloricNormalContact<'a, M> {
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        model: M,
        curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
        cells: [CellId; 2],
        partition: HeatPartition,
    ) -> Result<Self, NormalThermalError> {
        if cells[0] == cells[1] || partition.other() != 0.0 {
            return Err(NormalThermalError::Invalid(
                "two distinct finite caloric destinations required",
            ));
        }
        if model
            .expected_caloric_curves()
            .is_some_and(|expected| (0..2).any(|i| !std::ptr::eq(expected[i], curves[i])))
        {
            return Err(NormalThermalError::Invalid(
                "normal model requires its original caloric objects",
            ));
        }
        let temperature_k = temperatures(thermal, curves, cells)?;
        let (endpoint, modulus) = sample(
            &model,
            model.reference_request().indentation_m,
            0.0,
            temperature_k,
        )?;
        let thermodynamics = potentials(&endpoint, modulus, temperature_k)?;
        capacity(thermal, curves, cells, &endpoint, modulus, temperature_k)?;
        let mass_kg = cells.map(|id| thermal.cells()[&id].mass_kg);
        Ok(Self {
            model,
            modulus,
            mass_kg,
            curves,
            cells,
            partition,
            endpoint,
            temperature_k,
            thermodynamics,
            time_s: thermal.time_s(),
        })
    }

    pub const fn model(&self) -> &M {
        &self.model
    }
    pub const fn modulus(&self) -> NormalModulusState {
        self.modulus
    }
    pub const fn endpoint(&self) -> &PointNormalPatchReceipt {
        &self.endpoint
    }
    pub const fn thermodynamics(&self) -> &NormalThermodynamics {
        &self.thermodynamics
    }
    pub const fn temperature_k(&self) -> [f64; 2] {
        self.temperature_k
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Solve one monotone linear indentation segment and the original finite
    /// thermal network together. Resolve supplies only OTHER heat/conductances
    /// at this same proposed endpoint. Cut intervals at original thermal events.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        indentation_m: f64,
        controls: NormalThermalSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NormalThermalCandidate<'a, 't, S, M>, NormalThermalError<E>> {
        self.step_candidate_with_patch(
            thermal,
            target_s,
            indentation_m,
            controls,
            guess_j,
            |state, _, _| resolve(state),
            cancelled,
        )
    }

    /// Compose another physical channel inside the SAME thermal endpoint.
    /// The original normal owner supplies its already evaluated patch and work.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_candidate_with_patch<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        indentation_m: f64,
        controls: NormalThermalSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
            &PointNormalPatchReceipt,
            &NormalContactWork,
        ) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NormalThermalCandidate<'a, 't, S, M>, NormalThermalError<E>> {
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
            .map_err(NormalThermalError::Newton)?;
        if thermal.time_s() != self.time_s
            || temperatures(thermal, self.curves, self.cells)? != self.temperature_k
            || self.cells.map(|id| thermal.cells()[&id].mass_kg) != self.mass_kg
            || !indentation_m.is_finite()
            || indentation_m < 0.0
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
            return Err(NormalThermalError::Invalid(
                "incompatible normal/caloric state or controls",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, poll)?;
        let dt = step.duration_s();
        let delta = finite(indentation_m - self.endpoint.approach_m)?;
        let rate = finite(delta / dt)?;
        if delta != 0.0 && rate == 0.0 {
            return Err(NormalThermalError::Invalid(
                "normal displacement rate is unrepresentable",
            ));
        }
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |state| {
                    poll()?;
                    if trials.get() >= controls.maximum_contact_trials {
                        return Err(NormalThermalError::TrialBudget);
                    }
                    trials.set(trials.get() + 1);
                    let temperature_k = temperatures(state, self.curves, self.cells)?;
                    let (start, _) =
                        sample(&self.model, self.endpoint.approach_m, rate, temperature_k)?;
                    let (endpoint, modulus) =
                        sample(&self.model, indentation_m, rate, temperature_k)?;
                    let thermodynamics = potentials(&endpoint, modulus, temperature_k)?;
                    capacity(
                        state,
                        self.curves,
                        self.cells,
                        &endpoint,
                        modulus,
                        temperature_k,
                    )?;
                    let (_, damping) = modulus_and_damping(&self.model.reference_request().law)?;
                    // psi=C*E*d^(5/2). Factor the divided difference using
                    // r=sqrt(min(d0,d1)/max(d0,d1)); never subtract nearby powers.
                    let largest = if indentation_m >= start.approach_m {
                        &endpoint
                    } else {
                        &start
                    };
                    let mean_elastic_force_n = if largest.approach_m == 0.0 {
                        0.0
                    } else {
                        let smaller = indentation_m.min(start.approach_m);
                        let r = finite((smaller / largest.approach_m).sqrt())?;
                        let polynomial = 1.0 + r * (1.0 + r * (1.0 + r * (1.0 + r)));
                        finite(
                            finite(largest.reversible_energy_j / largest.approach_m)?
                                * finite(polynomial / (1.0 + r))?,
                        )?
                    };
                    let elastic_work_j = finite(mean_elastic_force_n * delta)?;
                    let basis_change = finite(elastic_work_j / modulus.modulus_pa)?;
                    let physical_dissipation_j = finite(
                        finite(damping * rate.abs())? * finite(mean_elastic_force_n * delta.abs())?,
                    )?;
                    let work_into_contact_j = finite(elastic_work_j + physical_dissipation_j)?;
                    let mean_normal_force_n =
                        finite(mean_elastic_force_n * finite(1.0 + damping * rate)?)?;
                    // Average both temperature coordinate paths. Their coefficient
                    // increments sum to the exact modulus change, so unequal-body
                    // exchange has no arbitrary surface ordering. This is still a
                    // first-order temperature/mechanics composition.
                    let a = modulus_state(&self.model, [temperature_k[0], self.temperature_k[1]])?;
                    let b = modulus_state(&self.model, [self.temperature_k[0], temperature_k[1]])?;
                    let old_basis =
                        finite(self.endpoint.reversible_energy_j / self.modulus.modulus_pa)?;
                    let coefficient_work = [
                        finite(
                            old_basis
                                * finite(
                                    0.5 * finite(a.modulus_pa - self.modulus.modulus_pa)?
                                        + 0.5 * finite(modulus.modulus_pa - b.modulus_pa)?,
                                )?,
                        )?,
                        finite(
                            old_basis
                                * finite(
                                    0.5 * finite(b.modulus_pa - self.modulus.modulus_pa)?
                                        + 0.5 * finite(modulus.modulus_pa - a.modulus_pa)?,
                                )?,
                        )?,
                    ];
                    let coefficient_free_energy_change_j =
                        finite(coefficient_work[0] + coefficient_work[1])?;
                    let mut entropy_energy_change = [0.0; 2];
                    let mut reversible_heat_j = [0.0; 2];
                    let mut physical_heat_j = [0.0; 2];
                    let mut deposited_heat_j = [0.0; 2];
                    for (i, fraction) in [self.partition.surface_a(), self.partition.surface_b()]
                        .into_iter()
                        .enumerate()
                    {
                        // Factor d(T*s) into a temperature increment and an
                        // entropy increment; the latter splits geometric and
                        // coefficient changes instead of subtracting large stores.
                        let entropy_change = finite(
                            -finite(modulus.gradient_pa_k[i] * basis_change)?
                                - finite(
                                    old_basis
                                        * finite(
                                            modulus.gradient_pa_k[i]
                                                - self.modulus.gradient_pa_k[i],
                                        )?,
                                )?,
                        )?;
                        entropy_energy_change[i] = finite(
                            finite(temperature_k[i] * entropy_change)?
                                + finite(
                                    finite(temperature_k[i] - self.temperature_k[i])?
                                        * self.thermodynamics.entropy_j_k[i],
                                )?,
                        )?;
                        reversible_heat_j[i] =
                            finite(-coefficient_work[i] - entropy_energy_change[i])?;
                        physical_heat_j[i] = finite(fraction * physical_dissipation_j)?;
                        deposited_heat_j[i] = finite(reversible_heat_j[i] + physical_heat_j[i])?;
                    }
                    let internal_energy_change_j = finite(
                        elastic_work_j
                            + coefficient_free_energy_change_j
                            + entropy_energy_change[0]
                            + entropy_energy_change[1],
                    )?;
                    let closure_residual_j = finite(
                        internal_energy_change_j + deposited_heat_j[0] + deposited_heat_j[1]
                            - work_into_contact_j,
                    )?;
                    if physical_dissipation_j < 0.0
                        || closure_residual_j.abs() > controls.total_energy_tolerance_j
                        || finite(
                            thermodynamics.internal_energy_j
                                - self.thermodynamics.internal_energy_j
                                - internal_energy_change_j,
                        )?
                        .abs()
                            > controls.total_energy_tolerance_j
                    {
                        return Err(NormalThermalError::Invalid(
                            "normal potential/work failed closure",
                        ));
                    }
                    let work = NormalContactWork {
                        indentation_change_m: delta,
                        indentation_rate_mps: rate,
                        mean_elastic_force_n,
                        mean_normal_force_n,
                        elastic_work_j,
                        physical_dissipation_j,
                        work_into_contact_j,
                        coefficient_free_energy_change_j,
                        internal_energy_change_j,
                        closure_residual_j,
                    };
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(state, &endpoint, &work).map_err(NormalThermalError::Source)?;
                    let external_heat_j = finite(
                        dt * state.external_step_power_w(&sources_w, &conductances_w_k, poll)?,
                    )?;
                    for (id, heat) in self.cells.into_iter().zip(deposited_heat_j) {
                        let power = finite(heat / dt)?;
                        if heat != 0.0 && power == 0.0 {
                            return Err(NormalThermalError::Invalid(
                                "normal heat vanished in power conversion",
                            ));
                        }
                        let total = sources_w.entry(id).or_insert(0.0);
                        *total = finite(*total + power)?;
                    }
                    Ok(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: (
                            endpoint,
                            modulus,
                            thermodynamics,
                            temperature_k,
                            work,
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
        let mut solver = construction.map_err(NormalThermalError::Newton)?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(NormalThermalError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(NormalThermalError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            next,
            (
                endpoint,
                modulus,
                thermodynamics,
                temperature_k,
                work,
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
            thermal_audit.stored_energy_change_j + work.internal_energy_change_j
                - work.work_into_contact_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(NormalThermalError::Invalid(
                "normal and caloric internal energy do not balance",
            ));
        }
        let mut contact = self.clone();
        contact.endpoint = endpoint.clone();
        contact.modulus = modulus;
        contact.temperature_k = temperature_k;
        contact.thermodynamics = thermodynamics;
        contact.time_s = target_s;
        poll()?;
        Ok(NormalThermalCandidate {
            thermal: next,
            contact,
            report: NormalThermalReport {
                endpoint,
                modulus,
                thermodynamics,
                work,
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

fn modulus_and_damping<E>(law: &NormalPatchLaw) -> Result<(f64, f64), NormalThermalError<E>> {
    match *law {
        NormalPatchLaw::HertzSpherePlane {
            reduced_modulus_pa, ..
        }
        | NormalPatchLaw::HertzEllipticParaboloid {
            reduced_modulus_pa, ..
        } => Ok((reduced_modulus_pa, 0.0)),
        NormalPatchLaw::HuntCrossleySphere {
            reduced_modulus_pa,
            dissipation_s_per_m,
            ..
        }
        | NormalPatchLaw::HuntCrossleyEllipticParaboloid {
            reduced_modulus_pa,
            dissipation_s_per_m,
            ..
        } => Ok((reduced_modulus_pa, dissipation_s_per_m)),
        _ => Err(NormalThermalError::Invalid(
            "normal law is not a Hertz point potential",
        )),
    }
}
fn temperatures<E>(
    thermal: &DepositionSimulation<'_>,
    curves: [&EquilibriumEnthalpyPhaseCurve; 2],
    cells: [CellId; 2],
) -> Result<[f64; 2], NormalThermalError<E>> {
    let mut temperatures = [0.0; 2];
    for i in 0..2 {
        let cell = thermal
            .cells()
            .get(&cells[i])
            .ok_or(NormalThermalError::Invalid(
                "inactive normal-contact caloric cell",
            ))?;
        if !thermal.uses_material_law(cell.material_id, curves[i]) {
            return Err(NormalThermalError::Invalid(
                "normal contact must retain original caloric objects",
            ));
        }
        let phase = curves[i]
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(NormalThermalError::Phase)?;
        if phase.liquid_mass_fraction() != 0.0 {
            return Err(NormalThermalError::Invalid(
                "solid normal contact refuses liquid caloric states",
            ));
        }
        temperatures[i] = phase.temperature_k();
    }
    Ok(temperatures)
}
fn finite<E>(value: f64) -> Result<f64, NormalThermalError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NormalThermalError::Invalid(
            "unrepresentable normal contact energy or coefficient",
        ))
    }
}
fn difference_error<E>(e: FiniteDifferenceError<NormalThermalError<E>>) -> NormalThermalError<E> {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => NormalThermalError::Invalid(s),
    }
}

fn model_error<E>(error: NormalThermalError) -> NormalThermalError<E> {
    match error {
        NormalThermalError::Invalid(s) => NormalThermalError::Invalid(s),
        NormalThermalError::Normal(e) => NormalThermalError::Normal(e),
        NormalThermalError::Phase(e) => NormalThermalError::Phase(e),
        NormalThermalError::Thermal(e) => NormalThermalError::Thermal(e),
        NormalThermalError::Source(e) => match e {},
        NormalThermalError::MaterialModel(e) => NormalThermalError::MaterialModel(e),
        NormalThermalError::Newton(e) => NormalThermalError::Newton(e),
        NormalThermalError::NotConverged(e) => NormalThermalError::NotConverged(e),
        NormalThermalError::Interrupted(e) => NormalThermalError::Interrupted(e),
        NormalThermalError::TrialBudget => NormalThermalError::TrialBudget,
        NormalThermalError::Cancelled => NormalThermalError::Cancelled,
    }
}

/// Original material-card compliance and finite caloric-source binding.
pub mod material;

/// Two-body inertial normal motion with the same reciprocal finite thermal solve.
pub mod dynamics;

/// Evolving Hertz load/area, source directional sliding and one finite heat solve.
pub mod directional;
