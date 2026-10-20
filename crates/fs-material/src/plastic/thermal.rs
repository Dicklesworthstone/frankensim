//! Source-bound uniform thermoplastic heat feedback using native J2 and transport.
//!
//! The declared reduction has energetic linear isotropic hardening, fixed
//! reference geometry/density, temperature-independent elastic/hardening coefficients,
//! optional constant sourced isotropic expansion, and a strictly increasing single-solid caloric
//! chart. Yield stress and Cp may depend on temperature. Under these fixed
//! reference-density assumptions, the stress-free caloric enthalpy is the caloric
//! internal-energy chart. Expansion includes reversible stress/temperature energy.
//! Variable elastic free energy, phase/volume work and finite
//! deformation require additional terms and are outside this reduction.
//!
//! Physical yield dissipation and signed reversible thermoelastic exchange are
//! passed separately to the existing thermal owner.
//! Hardening storage and backward-Euler numerical loss are never heated. Every
//! trial retains the same old strain, plastic history, enthalpy and clock. A
//! candidate is published only after source, heat and combined work gates pass.

mod expansion;
pub use expansion::{ThermalJ2, ThermoplasticIncrement};

use super::J2State;
use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use super::source::{
    HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity, j2_requirements,
    resolve_j2_plasticity,
};
use crate::phase::{
    EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use crate::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use crate::{SmallStrainLaw, Voigt};
use core::{convert::Infallible, fmt};
use fs_matdb::{ClaimSelection, EvaluationDecision, MaterialCard, PropertyValue, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Frequency, QuantitySpec, Temperature};

/// Hard ceiling on coupled within-step proposals; caller controls may tighten it.
pub const MAX_FEEDBACK_ITERATIONS: usize = 256;

/// Fixed external source context. The initial temperature must agree with h0.
/// Selection applies to the four mechanical roles; Cp and density retain the
/// exact claims already admitted by the supplied caloric owner. Optional
/// expansion has its own explicit source selection and retains that claim.
#[derive(Clone, Debug)]
pub struct SourceQuery<'a> {
    /// One immutable card owning the mechanical and caloric claims.
    pub card: &'a MaterialCard,
    /// Original source-derived Cp chart and density receipts.
    pub caloric: &'a ResolvedHeatCapacityCurve,
    /// Complete reference query, with fixed nonthermal coordinates.
    pub context: QueryPoint,
    /// Explicit absolute-temperature coordinate.
    pub temperature_axis: String,
    /// Initial mechanical claim selection; subsequent queries pin those IDs.
    pub selection: MaterialPropertySelection,
    /// Independent H or uniaxial Et source convention.
    pub hardening_basis: HardeningPropertyBasis,
    /// None explicitly selects the zero-expansion reduction. Some selects the
    /// canonical instantaneous linear expansion coefficient from this same card.
    /// It must be a constant over the complete caloric span; no missing-data fallback.
    pub thermal_expansion: Option<MaterialPropertySelection>,
    /// Positive tensor strain-norm ceiling, at most 0.05.
    pub strain_limit: f64,
}

/// Independent numerical gates, not a material-accuracy certificate.
#[derive(Clone, Copy, Debug)]
pub struct FeedbackConfig {
    /// Coupled proposal limit, bounded by MAX_FEEDBACK_ITERATIONS.
    pub max_iterations: usize,
    /// Positive absolute temperature fixed-point allowance [K].
    pub temperature_tolerance_k: f64,
    /// Positive independent allowance for each accepted heat/work gate [J].
    pub energy_tolerance_j: f64,
    /// Positive native endpoint yield allowance [Pa].
    pub yield_tolerance_pa: f64,
    /// Temperature fixed-point relaxation in (0, 1].
    pub relaxation: f64,
}

/// Source, transport, physical-account or bounded-iteration refusal.
#[derive(Debug)]
pub enum FeedbackError<E = Infallible> {
    Source(J2SourceError),
    Support {
        property: String,
        reason: String,
    },
    Mechanical(EnergyError),
    Caloric(PhaseStateError),
    Transport(E),
    Invalid(&'static str),
    Cancelled,
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
            Self::Source(e) => write!(f, "thermoplastic source: {e}"),
            Self::Support { property, reason } => {
                write!(f, "thermoplastic source {property}: {reason}")
            }
            Self::Mechanical(e) => write!(f, "thermoplastic mechanical update: {e}"),
            Self::Caloric(e) => write!(f, "thermoplastic caloric update: {e}"),
            Self::Transport(e) => write!(f, "thermoplastic transport: {e}"),
            Self::Invalid(s) => write!(f, "thermoplastic feedback: {s}"),
            Self::Cancelled => f.write_str("thermoplastic feedback: cancelled"),
            Self::NotConverged {
                iterations,
                temperature_defect_k,
                heat_defect_j,
                thermal_residual_j,
                combined_residual_j,
            } => write!(
                f,
                "thermoplastic feedback did not converge in {iterations} proposals: dT={temperature_defect_k} K, dQ={heat_defect_j} J, thermal={thermal_residual_j} J, combined={combined_residual_j} J"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FeedbackError<E> {}

/// Source admission shared by material points and spatial solid elements.
/// This owns no mechanical history, thermal inventory, geometry or integrator.
/// The complete caloric interval has the same constant elastic/hardening energy
/// and density assumptions as [`CaloricPoint`]. Only yield stress and Cp vary;
/// optional constant expansion contributes reversible thermoelastic energy.
#[derive(Clone, Debug)]
pub struct ThermoplasticSource<'a> {
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    context: QueryPoint,
    temperature_axis: String,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
    source: ResolvedJ2Plasticity,
    density_kg_m3: f64,
    expansion_per_k: f64,
    reference_temperature_k: f64,
}

/// One initially uniform specimen with immutable accepted native history.
/// Transport owns the lumped/Biot admission; uniform internal heating and fixed
/// reference geometry remain assumptions of this small-strain specimen model.
#[derive(Clone, Debug)]
pub struct CaloricPoint<'a> {
    binding: ThermoplasticSource<'a>,
    source: ResolvedJ2Plasticity,
    strain: Voigt,
    stress_pa: Voigt,
    history: J2State,
    state: EquilibriumPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    time_s: f64,
}

/// Complete candidate and the actual transport result. The caller accepts the
/// returned specimen and any finite-reservoir report together.
#[derive(Debug)]
pub struct FeedbackStep<'a, R> {
    /// Immutable specimen candidate, including the accepted native history.
    pub next: CaloricPoint<'a>,
    /// Native mechanical proposal evaluated at the actual thermal endpoint.
    /// Its strain argument is total strain minus thermal strain.
    pub mechanical: EnergyIncrement,
    /// Original thermal proposal and its concrete report.
    pub transport: UniformEnthalpyStep<R>,
    /// Number of within-step transport proposals used.
    pub iterations: usize,
    /// Native endpoint stress work over the reference specimen volume [J].
    pub mechanical_work_j: f64,
    /// Final-temperature physical dissipation [J].
    pub physical_heat_j: f64,
    /// Signed reversible thermoelastic exchange into the caloric chart [J].
    pub reversible_heat_j: f64,
    /// Change in T alpha tr(sigma) over the reference volume [J].
    pub thermoelastic_storage_change_j: f64,
    /// Plastic heat plus reversible exchange passed to the thermal proposal [J].
    pub applied_heat_j: f64,
    /// Absolute difference between applied and final net caloric heat [J].
    pub heat_defect_j: f64,
    /// Absolute temperature fixed-point defect [K].
    pub temperature_defect_k: f64,
    /// m delta_h - physical dissipation - reversible exchange - boundary heat [J].
    pub thermal_residual_j: f64,
    /// delta_elastic + delta_hardening + delta(T alpha tr(sigma)) + m delta_h + numerical_loss
    /// - endpoint mechanical work - boundary heat [J].
    pub combined_residual_j: f64,
}

impl<'a> ThermoplasticSource<'a> {
    /// Bind one source at its initial enthalpy, without creating a material point.
    /// The six source roles and optional expansion share one card, context and temperature. Every
    /// selected claim must continuously cover the complete caloric chart.
    /// E, nu, H (or Et), density and expansion must be validity-wide scalar constants;
    /// equal endpoint samples of a varying energetic coefficient do not suffice.
    pub fn try_new(
        query: SourceQuery<'a>,
        specific_enthalpy_j_kg: f64,
    ) -> Result<Self, FeedbackError> {
        let curve = query.caloric.curve();
        let knots = curve.knots();
        if curve.material_card_identity() != query.card.content_hash()
            || knots.len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
            || knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
            || knots
                .windows(2)
                .any(|p| p[0].temperature_k >= p[1].temperature_k)
        {
            return Err(FeedbackError::Invalid(
                "the shared caloric chart must be strictly increasing and fully solid",
            ));
        }
        let state = curve
            .state_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let declared_temperature = query
            .context
            .axis_quantities()
            .get(&query.temperature_axis)
            .copied();
        let supported_temperature = declared_temperature == Some(absolute)
            || query.temperature_axis == "T"
                && (declared_temperature.is_none()
                    || declared_temperature == Some(QuantitySpec::dimensional(Temperature::DIMS)));
        if query.context.axes().get(&query.temperature_axis).copied() != Some(state.temperature_k())
            || !supported_temperature
        {
            return Err(FeedbackError::Invalid(
                "initial source temperature must match h0; renamed axes require explicit absolute-temperature semantics",
            ));
        }
        let first = resolve_j2_plasticity(
            query.card,
            &query.context,
            query.selection,
            query.hardening_basis,
            query.strain_limit,
        )
        .map_err(FeedbackError::Source)?;
        let receipts = query.caloric.receipts();
        if receipts.len() < 4 || receipts.len() % 2 != 0 {
            return Err(FeedbackError::Invalid(
                "caloric owner lacks paired Cp and density source receipts",
            ));
        }
        let mut requirements = j2_requirements(query.hardening_basis)
            .map_err(|e| FeedbackError::Source(e.into()))?
            .to_vec();
        let mut pins: Vec<_> = first
            .material()
            .properties()
            .iter()
            .map(|p| {
                (
                    p.requirement().name().to_owned(),
                    p.answer().receipt.selected,
                )
            })
            .collect();
        for (index, dims) in [Dims([2, 0, -2, -1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])]
            .into_iter()
            .enumerate()
        {
            let receipt = &receipts[index];
            let claim =
                query
                    .card
                    .claims()
                    .claim(receipt.selected)
                    .ok_or(FeedbackError::Invalid(
                        "caloric source claim is absent from the shared card",
                    ))?;
            if claim.key.quantity() != QuantitySpec::dimensional(dims)
                || requirements.iter().any(|r| r.name() == claim.key.name())
            {
                return Err(FeedbackError::Invalid(
                    "six distinct roles with exact SI quantities are required",
                ));
            }
            requirements.push(
                ScalarPropertyRequirement::try_with_key(
                    &claim.key,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .map_err(|e| FeedbackError::Source(e.into()))?,
            );
            pins.push((claim.key.name().to_owned(), receipt.selected));
        }
        if let Some(selection) = query.thermal_expansion {
            use crate::state_point::{
                INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            };
            let requirement = ScalarPropertyRequirement::try_new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite,
            )
            .map_err(|e| FeedbackError::Source(e.into()))?;
            let expansion = resolve_material_state_point(
                query.card,
                &query.context,
                std::slice::from_ref(&requirement),
                selection,
            )
            .map_err(|e| FeedbackError::Source(e.into()))?;
            pins.push((
                requirement.name().to_owned(),
                expansion.properties()[0].answer().receipt.selected,
            ));
            requirements.push(requirement);
        }
        let selection = MaterialPropertySelection::PinnedByProperty(pins);
        let bundle = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            selection.clone(),
        )
        .map_err(|e| FeedbackError::Source(e.into()))?;
        if receipts.iter().any(|r| {
            r.axis_quantities != *bundle.axis_quantities()
                || !same_context(
                    &r.query_point,
                    bundle.query_point(),
                    &query.temperature_axis,
                )
        }) || receipts
            .chunks_exact(2)
            .any(|p| p[0].selected != receipts[0].selected || p[1].selected != receipts[1].selected)
        {
            return Err(FeedbackError::Invalid(
                "caloric and mechanical sources must retain one fixed context and the same pinned Cp/density claims",
            ));
        }
        let density_name = &receipts[1].property;
        let density = bundle
            .property(density_name)
            .ok_or(FeedbackError::Invalid(
                "density source is missing from the shared snapshot",
            ))?
            .value_si();
        if knots.iter().any(|k| k.bulk_density_kg_m3 != density) {
            return Err(FeedbackError::Invalid(
                "density must be constant over the entire solid caloric chart",
            ));
        }
        let lower = point_at(
            &query.context,
            &query.temperature_axis,
            knots[0].temperature_k,
        )?;
        let upper = point_at(
            &query.context,
            &query.temperature_axis,
            knots[knots.len() - 1].temperature_k,
        )?;
        for property in bundle.properties() {
            let receipt = &property.answer().receipt;
            let name = property.requirement().name();
            let claim =
                query
                    .card
                    .claims()
                    .claim(receipt.selected)
                    .ok_or(FeedbackError::Invalid(
                        "selected coefficient claim is absent",
                    ))?;
            let fixed = name == YOUNG_MODULUS_PROPERTY
                || name == POISSON_RATIO_PROPERTY
                || name == query.hardening_basis.property_name()
                || name == density_name
                || name == crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY;
            match &claim.value {
                PropertyValue::Scalar { .. }
                    if receipt.decision == EvaluationDecision::ConstantWithinValidity =>
                {
                    ()
                }
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    knots: samples,
                    ..
                } if !fixed
                    && abscissa == &query.temperature_axis
                    && *abscissa_dims == Temperature::DIMS =>
                {
                    if samples.iter().any(|&(t, value)| {
                        t >= knots[0].temperature_k
                            && t <= knots[knots.len() - 1].temperature_k
                            && !positive(value)
                    }) {
                        return Err(FeedbackError::Invalid(
                            "yield and Cp curves must stay positive throughout the admitted thermal span",
                        ));
                    }
                }
                _ => {
                    return Err(FeedbackError::Invalid(
                        "elasticity, hardening, density and expansion require scalar constants; only yield and Cp may be temperature curves",
                    ));
                }
            }
            // Mechanics, caloric data and linked transport use the same
            // rate-independent source admission.
            admit_rate_independent_property(property)?;
            let support = query
                .card
                .claims()
                .query_envelope_typed(
                    &claim.key,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(receipt.selected),
                )
                .map_err(|e| FeedbackError::Support {
                    property: name.to_owned(),
                    reason: e.to_string(),
                })?;
            if !fixed
                && (!positive(support.lower.evidence.value.value)
                    || !positive(support.upper.evidence.value.value))
            {
                return Err(FeedbackError::Invalid(
                    "yield and Cp must be positive at both thermal support boundaries",
                ));
            }
        }
        let source = ResolvedJ2Plasticity::try_from_material(
            bundle,
            query.hardening_basis,
            query.strain_limit,
        )
        .map_err(FeedbackError::Source)?;
        let expansion_per_k = source
            .material()
            .property(crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .map_or(0.0, |p| p.value_si());
        let result = Self {
            card: query.card,
            caloric: query.caloric,
            context: query.context,
            temperature_axis: query.temperature_axis,
            requirements,
            selection,
            source,
            density_kg_m3: density,
            expansion_per_k,
            reference_temperature_k: state.temperature_k(),
        };
        result.admit_thermoelastic_capacity()?;
        Ok(result)
    }

    /// Source bundle at the admitted reference temperature.
    #[must_use]
    pub const fn reference(&self) -> &ResolvedJ2Plasticity {
        &self.source
    }

    /// Original caloric owner, which must also back the coupled inventory.
    #[must_use]
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric
    }

    /// Constant delivered reference density [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Transport binding with the same immutable card and complete reference
    /// context used by mechanics and the caloric chart. Only temperature evolves
    /// in this fully solid, fixed-density reduction.
    #[must_use]
    pub fn thermal_material_point(&self) -> crate::thermal::source::ThermalMaterialPoint<'_> {
        crate::thermal::source::ThermalMaterialPoint {
            material: self.card,
            context: &self.context,
            temperature_axis: &self.temperature_axis,
            liquid_fraction_axis: None,
        }
    }

    /// Select one nonnegative static transport coefficient at this source's
    /// actual caloric state, retaining its claim identity for later heat steps.
    /// Temperature dependence is allowed. Mechanical-coordinate or frequency
    /// bounds need an evolving constitutive response and are refused here.
    /// The consuming transport owner still checks the exact property quantity.
    pub fn pin_thermal_transport(
        &self,
        specific_enthalpy_j_kg: f64,
        key: &fs_matdb::PropertyKey,
        selection: MaterialPropertySelection,
    ) -> Result<MaterialPropertySelection, FeedbackError> {
        let (state, _) = self.resolve_enthalpy(specific_enthalpy_j_kg)?;
        let thermal = self.thermal_material_point();
        let support_error =
            |e: crate::thermal::source::ThermalMaterialError| FeedbackError::Support {
                property: key.name().to_owned(),
                reason: e.to_string(),
            };
        let point = thermal.at(state).map_err(support_error)?;
        let requirement =
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                .map_err(|e| FeedbackError::Source(e.into()))?;
        let material = resolve_material_state_point(self.card, &point, &[requirement], selection)
            .map_err(|e| FeedbackError::Source(e.into()))?;
        let property = &material.properties()[0];
        thermal.admit_static(property).map_err(support_error)?;
        admit_rate_independent_property(property)?;
        Ok(MaterialPropertySelection::PinnedByProperty(vec![(
            key.name().to_owned(),
            property.answer().receipt.selected,
        )]))
    }

    /// Refresh all pinned roles at the actual caloric state. This does not
    /// advance plastic history or supply a second independently authored T.
    pub fn resolve_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<(EquilibriumPhaseState, ResolvedJ2Plasticity), FeedbackError> {
        let state = self
            .caloric
            .curve()
            .state_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map_err(FeedbackError::Caloric)?;
        Ok((state, self.resolve_temperature(state.temperature_k())?))
    }

    fn resolve_temperature<E>(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedJ2Plasticity, FeedbackError<E>> {
        let point = point_at(&self.context, &self.temperature_axis, temperature_k)?;
        let material = resolve_material_state_point(
            self.card,
            &point,
            &self.requirements,
            self.selection.clone(),
        )
        .map_err(|e| FeedbackError::Source(e.into()))?;
        let source = ResolvedJ2Plasticity::try_from_material(
            material,
            self.source.hardening_basis(),
            self.source.law().elastic.strain_limit,
        )
        .map_err(FeedbackError::Source)?;
        if source.law().elastic != self.source.law().elastic
            || source.law().hardening != self.source.law().hardening
        {
            return Err(FeedbackError::Invalid(
                "elastic and hardening energies changed without reversible caloric terms",
            ));
        }
        Ok(source)
    }
}

impl<'a> CaloricPoint<'a> {
    /// Birth at zero total strain with native virgin plastic history. A loaded
    /// history is reached only through accepted native increments. Source
    /// admission is shared with spatial solid coupling through ThermoplasticSource.
    pub fn try_new(
        query: SourceQuery<'a>,
        specific_enthalpy_j_kg: f64,
        volume_m3: f64,
        surface_area_m2: f64,
    ) -> Result<Self, FeedbackError> {
        if !positive(volume_m3) || !positive(surface_area_m2) {
            return Err(FeedbackError::Invalid(
                "positive finite reference volume and surface area are required",
            ));
        }
        let binding = ThermoplasticSource::try_new(query, specific_enthalpy_j_kg)?;
        let (state, source) = binding.resolve_enthalpy(specific_enthalpy_j_kg)?;
        let history = source.law().initial_state();
        let mass_kg = product(volume_m3, binding.density_kg_m3)?;
        Ok(Self {
            binding,
            source,
            strain: [0.0; 6],
            stress_pa: [0.0; 6],
            history,
            state,
            mass_kg,
            volume_m3,
            surface_area_m2,
            time_s: 0.0,
        })
    }

    /// Original complete source snapshot at the accepted temperature.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedJ2Plasticity {
        &self.source
    }
    /// Accepted total strain in the native tensor-Voigt convention.
    #[must_use]
    pub const fn strain(&self) -> &Voigt {
        &self.strain
    }
    /// Accepted native stress [Pa].
    #[must_use]
    pub const fn stress_pa(&self) -> &Voigt {
        &self.stress_pa
    }
    /// Accepted native plastic strain and accumulated history.
    #[must_use]
    pub const fn history(&self) -> &J2State {
        &self.history
    }
    /// Accepted state from the original caloric curve.
    #[must_use]
    pub const fn caloric_state(&self) -> &EquilibriumPhaseState {
        &self.state
    }
    /// Original source-derived Cp/density owner.
    #[must_use]
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.binding.caloric
    }
    /// Fixed reference mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Fixed reference volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Accepted elapsed physical time [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Prescribed tensor-strain increment with implicit physical-heat feedback.
    /// The callback must propose transport from the supplied accepted state;
    /// it must not mutate a physical reservoir during within-step trials.
    /// Its reported boundary heat is a separate signed transfer, not a residual
    /// inferred from internal dissipation. The transport allowance is capped by
    /// the caller's independent coupled energy gate.
    ///
    /// This is bounded fixed-point coupling of existing owners, not a new heat
    /// integrator or an unconditional convergence/stability claim. Every trial
    /// starts from the same old history. The final mechanics uses the actual
    /// returned temperature; rejection/cancellation leaves this point unchanged.
    pub fn advance<R, E>(
        &self,
        strain: Voigt,
        duration_s: f64,
        config: FeedbackConfig,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        poll(&mut cancelled)?;
        let final_time = self.time_s + duration_s;
        let material_tolerance = config.energy_tolerance_j / self.volume_m3;
        if config.max_iterations == 0
            || config.max_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.yield_tolerance_pa)
            || !positive(material_tolerance)
            || !positive(config.relaxation)
            || config.relaxation > 1.0
            || !positive(duration_s)
            || !final_time.is_finite()
            || final_time <= self.time_s
        {
            return Err(FeedbackError::Invalid(
                "invalid feedback controls or physical interval",
            ));
        }
        self.source
            .law()
            .validate_committed_history(
                &self
                    .binding
                    .mechanical_strain_at(self.state.temperature_k(), &self.strain)?,
                &self.history,
                config.yield_tolerance_pa,
            )
            .map_err(FeedbackError::Mechanical)?;
        let tolerances = EnergyTolerances {
            energy_j_m3: material_tolerance,
            yield_pa: config.yield_tolerance_pa,
        };
        let mut temperature = self.state.temperature_k();
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.max_iterations {
            poll(&mut cancelled)?;
            let (_, trial) = self.propose_at(temperature, &strain, tolerances)?;
            let applied_heat_j = product(self.volume_m3, trial.caloric_heat_j_m3())?;
            let proposal = transport(UniformEnthalpyStepInput {
                curve: self.binding.caloric.curve(),
                initial: self.state,
                mass_kg: self.mass_kg,
                volume_m3: self.volume_m3,
                surface_area_m2: self.surface_area_m2,
                internal_heat_j: applied_heat_j,
                duration_s,
            })
            .map_err(FeedbackError::Transport)?;
            poll(&mut cancelled)?;
            let checked = self
                .binding
                .caloric
                .curve()
                .state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg())
                .map_err(FeedbackError::Caloric)?;
            if checked != proposal.state
                || !proposal.external_heat_j.is_finite()
                || !proposal.energy_residual_tolerance_j.is_finite()
                || proposal.energy_residual_tolerance_j < 0.0
                || proposal.energy_residual_tolerance_j > config.energy_tolerance_j
            {
                return Err(FeedbackError::Invalid(
                    "transport returned a foreign state or an inadmissible energy allowance",
                ));
            }
            let caloric_change = product(
                self.mass_kg,
                proposal.state.specific_enthalpy_j_kg() - self.state.specific_enthalpy_j_kg(),
            )?;
            let transport_residual =
                finite(caloric_change - applied_heat_j - proposal.external_heat_j)?;
            if transport_residual.abs() > proposal.energy_residual_tolerance_j {
                return Err(FeedbackError::Invalid(
                    "transport proposal violates its own physical heat balance",
                ));
            }
            let final_temperature = proposal.state.temperature_k();
            let (source, increment) = self.propose_at(final_temperature, &strain, tolerances)?;
            let reversible_heat_j = product(self.volume_m3, increment.reversible_heat_j_m3)?;
            let thermoelastic_storage_change_j =
                product(self.volume_m3, increment.thermoelastic_storage_change_j_m3)?;
            let mechanical_work_j = product(self.volume_m3, increment.mechanical_work_j_m3)?;
            let mechanical = increment.mechanical;
            let physical_heat_j =
                product(self.volume_m3, mechanical.energy.physical_dissipation_j_m3)?;
            let caloric_heat_j = finite(physical_heat_j + reversible_heat_j)?;
            let thermal_residual =
                finite(caloric_change - caloric_heat_j - proposal.external_heat_j)?;
            let combined_residual = finite(
                product(self.volume_m3, mechanical.energy.stored_change_j_m3())?
                    + thermoelastic_storage_change_j
                    + caloric_change
                    + product(self.volume_m3, mechanical.energy.numerical_dissipation_j_m3)?
                    - mechanical_work_j
                    - proposal.external_heat_j,
            )?;
            defects = [
                (final_temperature - temperature).abs(),
                (caloric_heat_j - applied_heat_j).abs(),
                thermal_residual,
                combined_residual,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j
                && defects[3].abs() <= config.energy_tolerance_j
            {
                poll(&mut cancelled)?;
                let mut next = self.clone();
                next.source = source;
                next.strain = strain;
                next.stress_pa = mechanical.stress_pa;
                next.history = mechanical.state.clone();
                next.state = proposal.state;
                next.time_s = final_time;
                return Ok(FeedbackStep {
                    next,
                    mechanical,
                    transport: proposal,
                    iterations: iteration,
                    mechanical_work_j,
                    physical_heat_j,
                    reversible_heat_j,
                    thermoelastic_storage_change_j,
                    applied_heat_j,
                    heat_defect_j: defects[1],
                    temperature_defect_k: defects[0],
                    thermal_residual_j: defects[2],
                    combined_residual_j: defects[3],
                });
            }
            temperature = finite(
                (1.0 - config.relaxation) * temperature + config.relaxation * final_temperature,
            )?;
        }
        poll(&mut cancelled)?;
        Err(FeedbackError::NotConverged {
            iterations: config.max_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            thermal_residual_j: defects[2],
            combined_residual_j: defects[3],
        })
    }

    fn propose_at<E>(
        &self,
        temperature_k: f64,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(ResolvedJ2Plasticity, ThermoplasticIncrement), FeedbackError<E>> {
        self.binding.increment_at_temperature(
            self.state.temperature_k(),
            temperature_k,
            &self.strain,
            &self.history,
            strain,
            tolerances,
        )
    }
}

fn admit_rate_independent_property<E>(
    property: &crate::state_point::ResolvedScalarProperty,
) -> Result<(), FeedbackError<E>> {
    let receipt = &property.answer().receipt;
    for (axis, _) in &receipt.query_point {
        if property
            .answer()
            .evidence
            .model
            .validity
            .bound(axis)
            .is_none()
        {
            continue;
        }
        let quantity = receipt.axis_quantities.get(axis).copied();
        let kind = quantity
            .and_then(QuantitySpec::semantic_type)
            .map(|ty| ty.kind());
        let dynamic_kind = matches!(
            kind,
            Some(
                QuantityKind::Frequency(_)
                    | QuantityKind::AngularVelocity(_)
                    | QuantityKind::Strain { .. }
                    | QuantityKind::Stress
            )
        );
        let rate_dimension = quantity.is_some_and(|q| q.dims() == Frequency::DIMS);
        let legacy_dynamic_axis = kind.is_none()
            && [
                "omega",
                "frequency",
                "strain",
                "epsilon",
                "eps",
                "shear_strain",
                "strain_rate",
                "shear_rate",
                "plastic_strain",
                "equivalent_plastic_strain",
                "stress",
            ]
            .contains(&axis.as_str());
        if dynamic_kind || rate_dimension || legacy_dynamic_axis {
            return Err(FeedbackError::Invalid(
                "mechanically bounded coefficients require their evolving strain, stress or rate coordinate",
            ));
        }
    }
    if property
        .angular_frequency_band_rad_s()
        .map_err(|e| FeedbackError::Source(e.into()))?
        .is_some()
    {
        return Err(FeedbackError::Invalid(
            "frequency-bounded data cannot define rate-independent thermoplastic coefficients",
        ));
    }
    Ok(())
}

fn same_context(a: &[(String, f64)], b: &[(String, f64)], axis: &str) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|((ka, va), (kb, vb))| ka == kb && (ka == axis || va.to_bits() == vb.to_bits()))
}
fn point_at<E>(
    context: &QueryPoint,
    axis: &str,
    temperature_k: f64,
) -> Result<QueryPoint, FeedbackError<E>> {
    if !positive(temperature_k) {
        return Err(FeedbackError::Invalid(
            "source temperature must be positive finite kelvin",
        ));
    }
    context
        .clone()
        .with(axis, temperature_k)
        .map_err(|e| FeedbackError::Support {
            property: axis.to_owned(),
            reason: e.to_string(),
        })
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite<E>(value: f64) -> Result<f64, FeedbackError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FeedbackError::Invalid(
            "unrepresentable thermoplastic feedback arithmetic",
        ))
    }
}
fn product<E>(a: f64, b: f64) -> Result<f64, FeedbackError<E>> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(FeedbackError::Invalid("nonzero mass or energy underflowed"));
    }
    Ok(value)
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FeedbackError<E>> {
    if cancelled() {
        Err(FeedbackError::Cancelled)
    } else {
        Ok(())
    }
}
