//! Dissipation/heat feedback using the existing Maxwell update and thermal owner.
//!
//! One accepted mechanical history and one accepted caloric state are held fixed
//! through all within-step trials. Physical dashpot dissipation and signed
//! reversible thermoelastic exchange enter the SAME transport owner.
//! Backward-Euler numerical dissipation remains a separate, nonthermal loss;
//! the opt-in exact strain-ramp law has no such constitutive loss.
//!
//! This rung admits globally affine common spring-modulus and free-strain laws,
//! a strictly increasing single-solid-phase caloric chart, and constant density.
//! With reference spring potential Xi, psi=S(T)*Xi,
//! s=f'*sigma-S'*Xi and U=psi+T*s. Every original knot participates in affine
//! admission; derivative jumps refuse. The actual caloric chart must exceed a
//! conservative bound on the expansion contribution to fixed-strain heat capacity.
//! Relaxation times and caloric heat capacity may depend on temperature.
//! At fixed reference pressure/density, chart h changes are caloric internal
//! energy changes. Phase changes and finite-deformation volume work remain
//! unsupported. Endpoint integration can create numerical entropy; it is not
//! physical dashpot heating or a proof of exact adiabatic entropy conservation.

use core::{convert::Infallible, fmt};

use super::{Energy, HistoryIntegration, Point};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput};

/// Hard work ceiling; a failed fixed point is not silently accepted.
pub const MAX_FEEDBACK_ITERATIONS: usize = 256;

/// Independent numerical gates, not a physical accuracy certificate.
#[derive(Clone, Copy, Debug)]
pub struct FeedbackConfig {
    /// Thermal proposals allowed in this step, in 1..=256.
    pub max_iterations: usize,
    /// Absolute endpoint-temperature consistency tolerance [K].
    pub temperature_tolerance_k: f64,
    /// Absolute allowance for heat-source and combined energy defects [J].
    pub energy_tolerance_j: f64,
    /// Fixed-point temperature relaxation in (0,1].
    pub relaxation: f64,
}

/// Work-conjugate uniaxial endpoint control, independent of the material law.
/// Stress is work per reference volume per total engineering strain, matching
/// the existing small-strain Maxwell owner; it is not a finite-strain load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MechanicalControl {
    /// Imposed total engineering strain. The material supplies the reaction.
    TotalStrain(f64),
    /// Imposed endpoint stress [Pa], positive in tension. Strain is solved from
    /// the existing affine constitutive tangent at EACH thermal trial. A changing
    /// load uses endpoint work, as in the existing backward-Euler material law.
    Stress {
        stress_pa: f64,
        /// Independent absolute endpoint stress residual gate [Pa].
        tolerance_pa: f64,
    },
    /// Endpoint load line `stress = stress_at_reference - stiffness *
    /// (strain - reference_strain)`. A nonnegative stiffness describes a
    /// passive spring or the condensed inertia of an implicit structural step.
    /// The structural owner supplies its geometry, inertia and work accounts;
    /// this material control neither invents mass nor deposits algorithmic loss.
    AffineStress {
        reference_strain: f64,
        stress_at_reference_pa: f64,
        stiffness_pa: f64,
        tolerance_pa: f64,
    },
}

impl MechanicalControl {
    fn validate<E>(self) -> Result<(), FeedbackError<E>> {
        match self {
            Self::TotalStrain(strain) if strain.is_finite() => Ok(()),
            Self::Stress { stress_pa, tolerance_pa }
                if stress_pa.is_finite() && positive(tolerance_pa) => Ok(()),
            Self::AffineStress { reference_strain, stress_at_reference_pa, stiffness_pa, tolerance_pa }
                if reference_strain.is_finite() && stress_at_reference_pa.is_finite()
                    && stiffness_pa.is_finite() && stiffness_pa >= 0.0 && positive(tolerance_pa) => Ok(()),
            _ => Err(FeedbackError::Invalid("control needs finite strain/stress, nonnegative affine stiffness and positive stress tolerance")),
        }
    }

    fn target_stress<E>(self, strain: f64) -> Result<Option<(f64, f64)>, FeedbackError<E>> {
        match self {
            Self::TotalStrain(_) => Ok(None),
            Self::Stress { stress_pa, tolerance_pa } => Ok(Some((stress_pa, tolerance_pa))),
            Self::AffineStress { reference_strain, stress_at_reference_pa, stiffness_pa, tolerance_pa } =>
                Ok(Some((finite(stress_at_reference_pa - finite(stiffness_pa * finite(strain-reference_strain)?)?)?, tolerance_pa))),
        }
    }

    fn propose<'a, E>(self, old: &Point<'a>, temperature: f64, duration_s: f64,
        energy_tolerance: f64) -> Result<super::Step<'a>, FeedbackError<E>> {
        let strain = match self {
            Self::TotalStrain(strain) => strain,
            Self::Stress { stress_pa, .. } => {
                // This inversion changes the boundary condition, not the
                // constitutive update or its branch-history integration.
                let affine = old.linearize(temperature, duration_s)
                    .map_err(FeedbackError::Mechanical)?;
                finite(finite(stress_pa - affine.offset_pa)? / affine.tangent_pa)?
            }
            Self::AffineStress { reference_strain, stress_at_reference_pa, stiffness_pa, .. } => {
                let affine = old.linearize(temperature, duration_s)
                    .map_err(FeedbackError::Mechanical)?;
                // Center the load at its reference strain. Condensed inertia
                // can be very stiff, so do not subtract two large k*strain
                // offsets to obtain the small displacement increment.
                let at_reference = finite(finite(affine.tangent_pa * reference_strain)? + affine.offset_pa)?;
                finite(reference_strain + finite(finite(stress_at_reference_pa-at_reference)?
                    / finite(affine.tangent_pa+stiffness_pa)?)?)?
            }
        };
        let step = old.advance(temperature, strain, duration_s, energy_tolerance)
            .map_err(FeedbackError::Mechanical)?;
        if let Some((stress_pa, tolerance_pa)) = self.target_stress(step.point.total_strain())? {
            let residual = finite(step.stress_pa - stress_pa)?;
            if residual.abs() > tolerance_pa {
                return Err(FeedbackError::StressResidual { target_pa: stress_pa,
                    actual_pa: step.stress_pa, tolerance_pa });
            }
        }
        Ok(step)
    }
}

/// The original material, caloric or transport failure remains distinguishable.
#[derive(Debug)]
pub enum FeedbackError<E = Infallible> {
    /// Original Maxwell law refusal.
    Mechanical(super::Error),
    /// Original caloric-chart refusal.
    Caloric(PhaseStateError),
    /// Original transport-owner refusal, retaining its concrete type.
    Transport(E),
    /// Inconsistent state, unsupported physics or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The old mechanical and caloric histories have not been changed.
    Cancelled,
    /// Actual constitutive stress missed the independently prescribed load gate.
    StressResidual { target_pa: f64, actual_pa: f64, tolerance_pa: f64 },
    /// Work budget exhausted without a consistent physical endpoint.
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
    },
}
impl<E: fmt::Display> fmt::Display for FeedbackError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Mechanical(e) => write!(f, "thermovisco heat feedback: {e}"),
            Self::Caloric(e) => write!(f, "thermovisco heat feedback: {e}"),
            Self::Transport(e) => write!(f, "thermovisco heat transport: {e}"),
            Self::Invalid(message) => write!(f, "thermovisco heat feedback: {message}"),
            Self::Cancelled => f.write_str("thermovisco heat feedback: cancelled"),
            Self::StressResidual { target_pa, actual_pa, tolerance_pa } => write!(f,
                "thermovisco load control: target={target_pa} Pa, actual={actual_pa} Pa, allowance={tolerance_pa} Pa"),
            Self::NotConverged { iterations, temperature_defect_k, heat_defect_j,
                thermal_residual_j, combined_residual_j } => write!(f,
                "thermovisco heat feedback did not converge in {iterations} proposals: dT={temperature_defect_k} K, dQ={heat_defect_j} J, thermal={thermal_residual_j} J, combined={combined_residual_j} J"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FeedbackError<E> {}

/// One uniformly strained/isothermal specimen with fixed mass and volume.
///
/// The mechanical law and caloric chart are independently declared models;
/// binding them does not invent shared material-card provenance or qualify a
/// specimen. Transport owns the isothermal/Biot admission for its boundary law.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricPoint<'a> {
    point: Point<'a>,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    time_s: f64,
    modulus_slope_per_k: f64,
    free_strain_slope_per_k: f64,
    mechanical_entropy_j_k: f64,
    mechanical_internal_energy_j: f64,
}

/// A complete proposed endpoint and the actual accepted transport-owner result.
/// No external reservoir may be mutated until this proposal is accepted.
#[derive(Debug)]
pub struct FeedbackStep<'a, R> {
    pub next: CaloricPoint<'a>,
    /// The imposed boundary condition, not a fitted material parameter.
    pub control: MechanicalControl,
    /// Actual minus prescribed endpoint/load-line stress [Pa]; None for strain control.
    pub stress_residual_pa: Option<f64>,
    /// Boundary work [J]. For stress control this uses the prescribed
    /// load, not a possibly imperfect constitutive reaction. Under strain
    /// control it is the original law's mechanical work times reference volume,
    /// including interval-mean stress when exact strain ramps are selected.
    pub external_mechanical_work_j: f64,
    /// Original Maxwell endpoint audit, in J/m3; not recomputed by this adapter.
    /// Its spring storage is FREE energy. The terms below complete the caloric account.
    pub mechanical_energy: Energy,
    /// Change in mechanical internal energy psi+T*s, including all retained branches [J].
    pub mechanical_internal_energy_change_j: f64,
    /// Change in the mechanical entropy contribution T*s to internal energy [J].
    pub mechanical_entropy_energy_change_j: f64,
    /// Signed reversible exchange with caloric storage, distinct from dashpot loss [J].
    pub reversible_heat_j: f64,
    /// Final thermal proposal, including its concrete report and boundary heat.
    pub transport: UniformEnthalpyStep<R>,
    pub iterations: usize,
    /// Physical dashpot heat PLUS signed reversible exchange supplied to the
    /// final thermal proposal [J]. May be negative. Numerical loss is excluded.
    pub applied_heat_j: f64,
    /// Mismatch against total material heat reevaluated at the final temperature [J].
    pub heat_defect_j: f64,
    pub temperature_defect_k: f64,
    /// m delta_h - dashpot heat - reversible exchange - external heat [J].
    pub thermal_residual_j: f64,
    /// delta_internal + m delta_h + numerical_loss - imposed work - external heat [J].
    pub combined_residual_j: f64,
}

impl<'a> CaloricPoint<'a> {
    /// Bind an existing mechanical history without resetting its branch strains.
    /// Its temperature must exactly match the supplied caloric coordinate.
    ///
    /// # Errors
    /// Refuses non-solid/variable-density caloric charts, non-affine spring or
    /// free-strain laws, inadequate total heat capacity, inconsistent mass/volume,
    /// and invalid state. Expansion is a reference eigenstrain, not changing mass
    /// or finite-volume kinematics.
    pub fn try_new(point: Point<'a>, curve: &'a EquilibriumEnthalpyPhaseCurve,
        specific_enthalpy_j_kg: f64, mass_kg: f64, volume_m3: f64,
        surface_area_m2: f64) -> Result<Self, FeedbackError> {
        for value in [mass_kg, volume_m3, surface_area_m2] {
            if !positive(value) { return Err(FeedbackError::Invalid("mass, volume and area must be finite and positive")); }
        }
        if curve.knots().len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(FeedbackError::Invalid("caloric chart exceeds the bounded knot count"));
        }
        let knots = &point.law.knots;
        let reference = knots[0];
        let modulus_knots: Vec<_> = knots.iter()
            .map(|k| (k.temperature_k, k.modulus_scale)).collect();
        let modulus_slope_per_k = crate::state_point::affine_temperature_slope(
            &modulus_knots, reference.temperature_k, reference.modulus_scale)
            .map_err(FeedbackError::Invalid)?;
        let free_strain_knots: Vec<_> = knots.iter()
            .map(|k| (k.temperature_k, k.free_strain)).collect();
        let free_strain_slope_per_k = crate::state_point::affine_temperature_slope(
            &free_strain_knots, reference.temperature_k, reference.free_strain)
            .map_err(|_| FeedbackError::Invalid("free thermal strain must be globally affine over every original knot"))?;
        let density = curve.knots()[0].bulk_density_kg_m3;
        let represented_mass = product(volume_m3, density)?;
        if (represented_mass / mass_kg - 1.0).abs() > 32.0 * f64::EPSILON
            || curve.knots().iter().any(|k| k.liquid_mass_fraction != 0.0
                || k.bulk_density_kg_m3 != density)
            || curve.knots().windows(2).any(|pair| pair[0].temperature_k >= pair[1].temperature_k) {
            return Err(FeedbackError::Invalid("heat feedback needs a fixed-density, strictly increasing solid chart matching reference mass"));
        }
        let state = curve.state_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map_err(FeedbackError::Caloric)?;
        if state.temperature_k() != point.temperature_k()
            || !point.stored_energy_j_m3().is_finite() {
            return Err(FeedbackError::Invalid("mechanical history and caloric temperature must match"));
        }
        admit_expansion_capacity(&point, curve, modulus_slope_per_k, free_strain_slope_per_k)?;
        let (mechanical_entropy_j_k, mechanical_internal_energy_j) =
            thermodynamics(&point, modulus_slope_per_k, free_strain_slope_per_k, volume_m3)?;
        Ok(Self { point, curve, state, mass_kg, volume_m3, surface_area_m2, time_s: 0.0,
            modulus_slope_per_k, free_strain_slope_per_k,
            mechanical_entropy_j_k, mechanical_internal_energy_j })
    }

    #[must_use]
    pub const fn mechanical(&self) -> &Point<'a> { &self.point }
    #[must_use]
    pub const fn caloric_state(&self) -> &EquilibriumPhaseState { &self.state }
    #[must_use]
    pub const fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.curve }
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }
    /// Mechanical entropy from the admitted spring and expansion free energy, excluding the
    /// independent caloric chart [J/K]. This is a constitutive contribution,
    /// not the total specimen entropy or an entropy-production receipt.
    #[must_use]
    pub const fn mechanical_entropy_j_k(&self) -> f64 { self.mechanical_entropy_j_k }
    /// Mechanical internal energy psi+T*s [J], including retained branch state.
    /// Add the caloric energy with its declared reference to form total storage.
    #[must_use]
    pub const fn mechanical_internal_energy_j(&self) -> f64 { self.mechanical_internal_energy_j }
    /// Original declared material strain-rate ceiling [1/s], also applicable
    /// to the initial velocity of a structural owner using this specimen.
    #[must_use]
    pub const fn max_strain_rate_s(&self) -> f64 { self.point.law.limits.max_strain_rate_s }

    /// Propose one prescribed-strain step with implicit damping/heat feedback.
    ///
    /// Pass an existing thermal owner, e.g. `|input| environment.advance(cx,input)`.
    /// The callback must be side-effect-free with respect to physical history:
    /// EVERY trial uses the same old enthalpy, old Maxwell history and duration.
    /// Only its temperature guess changes. The final mechanical update uses the
    /// actual returned thermal temperature, never a stale fixed-point iterate.
    ///
    /// This is bounded fixed-point iteration, not unconditional convergence or
    /// stability. Temperature and physical-heat defects are independently gated.
    /// Numerical dissipation is never deposited as heat. Failure or cancellation
    /// returns no candidate; the original specimen and histories are immutable.
    pub fn advance<R, E>(&self, total_strain: f64, duration_s: f64,
        config: FeedbackConfig,
        transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        cancelled: impl FnMut() -> bool) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        self.advance_controlled(MechanicalControl::TotalStrain(total_strain), duration_s,
            config, transport, cancelled)
    }

    /// Solve strain or stress control and heat feedback in the SAME iteration.
    /// Stress control enables creep, unloading and signed load reversal without
    /// resetting viscoelastic history. Every trial and the returned endpoint
    /// solve the affine stress law from the same old point. The original affine
    /// tangent is a fixed-temperature constitutive tangent, not a derivative of
    /// this coupled thermal solve. This method supplies no structural dynamics.
    ///
    /// Prescribed endpoint work enters the combined energy gate independently
    /// of the stress gate; an allowed force residual is not relabeled as heat.
    /// Affine stress control is solved at every thermal trial as well; its
    /// supplied load-line parameters stay fixed for the whole proposed step.
    /// The original transport transaction, source-domain, work and cancellation
    /// limits apply to all control modes. A rejected proposal changes neither
    /// strain, branch memory, enthalpy nor the accepted clock.
    /// Exact strain-ramp integration admits prescribed strain here. Stress and
    /// affine-stress controls retain endpoint-load work and therefore refuse
    /// that integration mode; their complete loading path is not supplied by
    /// an endpoint stress alone.
    pub fn advance_controlled<R, E>(&self, control: MechanicalControl, duration_s: f64,
        config: FeedbackConfig,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        poll(&mut cancelled)?;
        control.validate()?;
        if self.point.law.history_integration == HistoryIntegration::ExactStrainRamp
            && !matches!(control, MechanicalControl::TotalStrain(_)) {
            return Err(FeedbackError::Invalid(
                "exact strain-ramp history requires prescribed strain; endpoint load work is a different boundary history"));
        }
        let final_time = self.time_s + duration_s;
        let material_tolerance = config.energy_tolerance_j / self.volume_m3;
        if config.max_iterations == 0 || config.max_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k) || !positive(config.energy_tolerance_j)
            || !positive(config.relaxation) || config.relaxation > 1.0
            || !positive(duration_s) || !positive(material_tolerance)
            || !final_time.is_finite() || final_time <= self.time_s {
            return Err(FeedbackError::Invalid("invalid feedback controls, strain or physical interval"));
        }
        let mut temperature = self.point.temperature_k();
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.max_iterations {
            poll(&mut cancelled)?;
            let trial = control.propose(&self.point, temperature, duration_s, material_tolerance)?;
            let applied_heat_j = self.transfer(&trial)?.heat_j;
            let proposal = transport(UniformEnthalpyStepInput { curve: self.curve, initial: self.state,
                mass_kg: self.mass_kg, volume_m3: self.volume_m3, surface_area_m2: self.surface_area_m2,
                internal_heat_j: applied_heat_j, duration_s }).map_err(FeedbackError::Transport)?;
            poll(&mut cancelled)?;
            let checked = self.curve.state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg())
                .map_err(FeedbackError::Caloric)?;
            if checked != proposal.state || !proposal.external_heat_j.is_finite()
                || !proposal.energy_residual_tolerance_j.is_finite()
                || proposal.energy_residual_tolerance_j < 0.0
                || proposal.energy_residual_tolerance_j > config.energy_tolerance_j {
                return Err(FeedbackError::Invalid("transport returned a foreign state or an inadmissible energy allowance"));
            }
            let thermal_change = product(self.mass_kg,
                proposal.state.specific_enthalpy_j_kg() - self.state.specific_enthalpy_j_kg())?;
            let transport_residual = finite(thermal_change - applied_heat_j - proposal.external_heat_j)?;
            if transport_residual.abs() > proposal.energy_residual_tolerance_j {
                return Err(FeedbackError::Invalid("transport proposal violates its own physical heat balance"));
            }
            let final_temperature = proposal.state.temperature_k();
            let mechanical = control.propose(&self.point, final_temperature, duration_s, material_tolerance)?;
            let exchange = self.transfer(&mechanical)?;
            let actual_heat = exchange.heat_j;
            let thermal_residual = finite(thermal_change - actual_heat - proposal.external_heat_j)?;
            let target = control.target_stress(mechanical.point.total_strain())?;
            let stress_residual_pa = target.map(|(stress, _)| finite(mechanical.stress_pa-stress)).transpose()?;
            let work_density = match target {
                None => mechanical.energy.mechanical_work,
                Some((stress_pa, _)) => product(stress_pa,
                    mechanical.point.total_strain() - self.point.total_strain())?,
            };
            let external_mechanical_work_j = product(self.volume_m3, work_density)?;
            let combined_residual = finite(exchange.internal_change_j
                + thermal_change + product(self.volume_m3, mechanical.energy.numerical_dissipation)?
                - external_mechanical_work_j - proposal.external_heat_j)?;
            defects = [(final_temperature - temperature).abs(), (actual_heat - applied_heat_j).abs(),
                thermal_residual, combined_residual];
            if defects[0] <= config.temperature_tolerance_k && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j && defects[3].abs() <= config.energy_tolerance_j {
                poll(&mut cancelled)?;
                let next = Self { point: mechanical.point, curve: self.curve, state: proposal.state,
                    mass_kg: self.mass_kg, volume_m3: self.volume_m3,
                    surface_area_m2: self.surface_area_m2, time_s: final_time,
                    modulus_slope_per_k: self.modulus_slope_per_k,
                    free_strain_slope_per_k: self.free_strain_slope_per_k,
                    mechanical_entropy_j_k: exchange.entropy_j_k,
                    mechanical_internal_energy_j: exchange.internal_energy_j };
                return Ok(FeedbackStep { next, control, stress_residual_pa, external_mechanical_work_j,
                    mechanical_energy: mechanical.energy, transport: proposal,
                    mechanical_internal_energy_change_j: exchange.internal_change_j,
                    mechanical_entropy_energy_change_j: exchange.entropy_energy_change_j,
                    reversible_heat_j: exchange.reversible_heat_j,
                    iterations: iteration, applied_heat_j, heat_defect_j: defects[1],
                    temperature_defect_k: defects[0], thermal_residual_j: defects[2], combined_residual_j: defects[3] });
            }
            temperature = finite((1.0 - config.relaxation) * temperature + config.relaxation * final_temperature)?;
        }
        poll(&mut cancelled)?;
        Err(FeedbackError::NotConverged { iterations: config.max_iterations,
            temperature_defect_k: defects[0], heat_defect_j: defects[1],
            thermal_residual_j: defects[2], combined_residual_j: defects[3] })
    }

    fn transfer<E>(&self, step: &super::Step<'_>) -> Result<ThermalExchange, FeedbackError<E>> {
        let (entropy_j_k, internal_energy_j) = thermodynamics(&step.point,
            self.modulus_slope_per_k, self.free_strain_slope_per_k, self.volume_m3)?;
        let entropy_energy_change_j = if self.modulus_slope_per_k == 0.0
            && self.free_strain_slope_per_k == 0.0 { 0.0 } else {
            finite(product(step.point.temperature_k(), entropy_j_k)?
                - product(self.point.temperature_k(), self.mechanical_entropy_j_k)?)?
        };
        // The original native balance is delta(psi)+D+N=W+C_psi+W_thermal,
        // with W_thermal=-sigma_work*delta(f) from the original eigenstrain:
        // endpoint stress for BE, interval-mean stress for exact strain ramps.
        // Adding delta(T*s) gives internal energy; caloric exchange cancels
        // C_psi+W_thermal+delta(T*s), with N never supplied to transport.
        let reversible_heat_j = finite(-product(self.volume_m3, step.energy.modulus_work)?
            - product(self.volume_m3, step.energy.thermal_strain_work)?
            - entropy_energy_change_j)?;
        let dissipated_j = product(self.volume_m3, step.energy.dissipated)?;
        if dissipated_j < 0.0 {
            return Err(FeedbackError::Invalid("native Maxwell dissipation must be nonnegative"));
        }
        Ok(ThermalExchange { entropy_j_k, internal_energy_j,
            internal_change_j: finite(internal_energy_j - self.mechanical_internal_energy_j)?,
            entropy_energy_change_j, reversible_heat_j,
            heat_j: finite(dissipated_j + reversible_heat_j)? })
    }
}

struct ThermalExchange {
    entropy_j_k: f64,
    internal_energy_j: f64,
    internal_change_j: f64,
    entropy_energy_change_j: f64,
    reversible_heat_j: f64,
    heat_j: f64,
}

fn thermodynamics<E>(point: &Point<'_>, slope_per_k: f64, expansion_per_k: f64, volume_m3: f64)
    -> Result<(f64, f64), FeedbackError<E>> {
    let modulus_entropy_j_k = if slope_per_k == 0.0 { 0.0 } else {
        product(volume_m3, product(-slope_per_k, point.unscaled_storage())?)?
    };
    let entropy_j_k = if expansion_per_k == 0.0 { modulus_entropy_j_k } else {
        // With v_j=e-q_j, BE gives v_new=(1-a)*e_new+a*v_old.
        // Exact ramps give a*v_old+(b-a)*e_old+(1-b)*e_new, where
        // a=exp(-r) <= b=(1-exp(-r))/r <= 1. Both are convex combinations;
        // born v=0 and admitted |e|<=limit therefore imply |q|<=2*limit.
        // Enforce that envelope on the actual history as well: the conservative
        // heat-capacity admission must never silently cover a larger state.
        let bound = product(2.0, point.law.limits.max_strain)?;
        if point.q.iter().any(|q| !q.is_finite() || q.abs() > bound) {
            return Err(FeedbackError::Invalid("thermal expansion exceeds the admitted branch-strain capacity envelope"));
        }
        finite(modulus_entropy_j_k
            + product(volume_m3, product(expansion_per_k, point.stress_pa())?)?)?
    };
    let internal_energy_j = finite(product(volume_m3, point.stored_energy_j_m3())?
        + product(point.temperature_k(), entropy_j_k)?)?;
    Ok((entropy_j_k, internal_energy_j))
}

fn admit_expansion_capacity<E>(point: &Point<'_>, curve: &EquilibriumEnthalpyPhaseCurve,
    modulus_slope: f64, expansion: f64) -> Result<(), FeedbackError<E>> {
    if expansion == 0.0 { return Ok(()); }
    // At fixed total strain and viscous strains, psi_TT=alpha²*S*K0-2*alpha*b*A0,
    // K0=E_inf+sum(E_j), A0=E_inf*e+sum(E_j*q_j). Total volumetric capacity is
    // rho*dh/dT-T*psi_TT. Bound every admitted strain/history over the complete
    // shared temperature range, not merely the currently stress-free endpoint.
    let law = point.law;
    let branch_modulus = law.branches.iter().try_fold(0.0,
        |sum, &(e, _)| finite::<E>(sum+e))?;
    let instantaneous = finite(law.equilibrium_pa+branch_modulus)?;
    let stress_bound = product(law.limits.max_strain,
        finite(law.equilibrium_pa+product(2.0, branch_modulus)?)?)?;
    let scale = law.knots.iter().map(|k| k.modulus_scale).fold(0.0, f64::max);
    let temperature = law.knots.last().expect("bounded law").temperature_k
        .min(curve.knots().last().expect("bounded caloric chart").temperature_k);
    let expansion_term = product(product(expansion, expansion)?, product(scale, instantaneous)?)?;
    let mixed_term = product(product(2.0, product(expansion, modulus_slope)?.abs())?, stress_bound)?;
    let capacity_bound = product(temperature, finite(expansion_term+mixed_term)?)?;
    let mut caloric_capacity = f64::INFINITY;
    for pair in curve.knots().windows(2) {
        let dh = finite(pair[1].specific_enthalpy_j_kg-pair[0].specific_enthalpy_j_kg)?;
        let dt = finite(pair[1].temperature_k-pair[0].temperature_k)?;
        let cp = finite(dh/dt)?;
        if !positive(cp) { return Err(FeedbackError::Invalid("positive caloric segment capacity required")); }
        caloric_capacity = caloric_capacity.min(product(pair[0].bulk_density_kg_m3, cp)?);
    }
    if !caloric_capacity.is_finite() || caloric_capacity <= capacity_bound {
        return Err(FeedbackError::Invalid("caloric capacity does not exceed the full thermal-expansion coupling bound"));
    }
    Ok(())
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn finite<E>(value: f64) -> Result<f64, FeedbackError<E>> {
    if value.is_finite() { Ok(value) } else { Err(FeedbackError::Invalid("unrepresentable feedback arithmetic")) }
}
fn product<E>(a: f64, b: f64) -> Result<f64, FeedbackError<E>> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(FeedbackError::Invalid("nonzero energy or mass underflowed"));
    }
    Ok(value)
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FeedbackError<E>> {
    if cancelled() { Err(FeedbackError::Cancelled) } else { Ok(()) }
}
