//! Source-bound uniform thermoplastic heat feedback using native J2 and transport.
//!
//! The declared reduction has energetic linear isotropic hardening and optional
//! Prager kinematic hardening, fixed
//! reference geometry/density and Poisson ratio, globally affine E(T) and native
//! energetic H_iso(T)/H_k(T), optional sourced isotropic expansion, and a strictly
//! increasing single-solid caloric chart. Yield stress, Cp and expansion may be
//! piecewise linear.
//! Explicit Voce hardening adds affine Q(T) and a fixed positive b to this same
//! retained history. Its recoverable saturation energy and entropy enter the
//! reversible account; saturation storage is never deposited as plastic heat.
//! Under these fixed-reference-density assumptions, stress-free caloric enthalpy
//! is the caloric internal-energy chart. Coefficient work and temperature derivatives
//! of elastic/hardening free energy enter the same reversible heat/storage account.
//! Kinked moduli require strain/history-dependent entropy jumps and are refused.
//! Optional linear or power-law J2 overstress uses sourced kinetics and the actual
//! physical interval in the same radial return, tangent and work account.
//! Without that explicit source selection the original rate-independent law
//! is unchanged. A duration-free response refuses when overstress is selected.
//! Phase/volume work, general nonlinear coefficient laws and finite deformation
//! remain outside this reduction.
//!
//! Physical yield/overstress dissipation and signed reversible thermoelastic exchange are
//! passed separately to the existing thermal owner.
//! Hardening storage and backward-Euler numerical loss are never heated. Every
//! trial retains the same old strain, plastic history, enthalpy and clock. A
//! candidate is published only after source, heat and combined work gates pass.

mod expansion;
pub use expansion::{ThermalJ2, ThermoplasticIncrement};
mod voce;
use voce::VoceSelection;

use super::J2State;
use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use super::source::{HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity,
    KINEMATIC_HARDENING_MODULUS_PROPERTY, j2_requirements};
use crate::phase::{
    EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use crate::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, affine_temperature_slope,
    resolve_material_state_point,
};
use crate::{SmallStrainLaw, Voigt};
use core::{convert::Infallible, fmt};
use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

/// Hard ceiling on coupled within-step proposals; caller controls may tighten it.
pub const MAX_FEEDBACK_ITERATIONS: usize = 256;

/// Fixed external source context. The initial temperature must agree with h0.
/// Selection applies to the four or five mechanical roles; Cp and density retain the
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
    /// Independent H, uniaxial Et, or separately sourced H_iso/H_k convention.
    pub hardening_basis: HardeningPropertyBasis,
    /// None explicitly selects the zero-expansion reduction. Some selects the
    /// canonical instantaneous linear expansion coefficient from this same card.
    /// Its complete temperature path is integrated exactly; scalar, piecewise-linear
    /// and fixed-context surface claims are admitted. No missing-data fallback.
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
/// The complete caloric interval admits globally affine Young's modulus and native
/// energetic hardening, with fixed Poisson ratio and density. Yield stress, Cp and
/// optional expansion may be piecewise linear. Modulus derivatives and integrated
/// expansion contribute both reversible heat and their internal-energy storage.
/// An explicitly selected viscosity is refreshed from the same complete bundle
/// at each trial and accepted temperature; it changes kinetics, not stored energy.
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
    // Empty for a constant/absent coefficient, preserving its original arithmetic.
    // Otherwise the source owner's complete fixed-context line, including cell cuts.
    expansion_knots: Vec<(f64, f64)>,
    minimum_source_heat_capacity_j_kg_k: f64,
    youngs_temperature_slope_pa_k: f64,
    hardening_temperature_slope_pa_k: f64,
    kinematic_hardening_temperature_slope_pa_k: f64,
    maximum_youngs_modulus_pa: f64,
    reference_temperature_k: f64,
    overstress: Option<OverstressSelection>,
    voce: Option<VoceSelection>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OverstressSelection { Linear, PowerLaw }

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
    /// Frozen-history free-energy change caused by sourced E(T)/H_iso(T)/H_k(T) [J].
    /// Add this to mechanical.energy.stored_change_j_m3() times reference volume.
    pub coefficient_storage_change_j: f64,
    /// Change in T times mechanical entropy, separate from free energy [J].
    pub thermoelastic_storage_change_j: f64,
    /// Plastic heat plus reversible exchange passed to the thermal proposal [J].
    pub applied_heat_j: f64,
    /// Absolute difference between applied and final net caloric heat [J].
    pub heat_defect_j: f64,
    /// Absolute temperature fixed-point defect [K].
    pub temperature_defect_k: f64,
    /// m delta_h - physical dissipation - reversible exchange - boundary heat [J].
    pub thermal_residual_j: f64,
    /// Actual delta_free_energy + delta(T s_mechanical) + m delta_h + numerical_loss
    /// - endpoint mechanical work - boundary heat [J].
    pub combined_residual_j: f64,
}

impl<'a> ThermoplasticSource<'a> {
    /// Bind one source at its initial enthalpy, without creating a material point.
    /// Mechanical/caloric source roles and optional expansion share one card, context and temperature. Every
    /// selected claim must continuously cover the complete caloric chart.
    /// E and both native energetic hardening moduli may vary affinely over the
    /// complete temperature span; nu and reference density remain fixed.
    /// A nonzero uniaxial Et reduction remains
    /// constant: its generally rational H(T) conversion is not an affine law.
    /// Expansion is integrated along its complete pinned path. Every crossed
    /// source knot/surface edge is checked; equal endpoints alone do not suffice.
    /// Modulus derivative jumps refuse because this caloric chart has no
    /// strain/history-dependent entropy discontinuity.
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
        let mut requirements = j2_requirements(query.hardening_basis)
            .map_err(|e| FeedbackError::Source(e.into()))?;
        let first_material = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            query.selection,
        )
        .map_err(|e| FeedbackError::Source(e.into()))?;
        for property in first_material.properties() {
            admit_rate_independent_property(property, &query.temperature_axis)?;
        }
        let first = ResolvedJ2Plasticity::try_from_material(
            first_material,
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
                    "distinct mechanical and caloric roles with exact SI quantities are required",
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
        let mut expansion_knots = Vec::new();
        let mut youngs_knots = Vec::new();
        let mut hardening_knots = Vec::new();
        let mut kinematic_hardening_knots = Vec::new();
        let mut minimum_source_heat_capacity_j_kg_k = f64::INFINITY;
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
            let expansion =
                name == crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY;
            let youngs = name == YOUNG_MODULUS_PROPERTY;
            let hardening = name == query.hardening_basis.property_name();
            let kinematic = name == KINEMATIC_HARDENING_MODULUS_PROPERTY;
            let fixed = name == POISSON_RATIO_PROPERTY || name == density_name;
            // Mechanics, caloric data and linked transport use the same
            // admission of instantaneous coefficients at fixed external context.
            admit_rate_independent_property(property, &query.temperature_axis)?;
            let support = query
                .card
                .claims()
                .query_line_typed(
                    &claim.key,
                    &query.temperature_axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(receipt.selected),
                )
                .map_err(|e| FeedbackError::Support {
                    property: name.to_owned(),
                    reason: e.to_string(),
                })?;
            for (temperature, sample) in support.knots {
                let value = sample.evidence.value.value;
                if expansion {
                    if !value.is_finite() {
                        return Err(FeedbackError::Invalid(
                            "expansion must be finite over the complete source path",
                        ));
                    }
                    expansion_knots.push((temperature, value));
                }
                if youngs {
                    youngs_knots.push((temperature, value));
                }
                if hardening || kinematic {
                    if !value.is_finite() || value < 0.0 {
                        return Err(FeedbackError::Invalid(
                            "hardening must be finite and nonnegative over the complete source path",
                        ));
                    }
                    if hardening {
                        hardening_knots.push((temperature, value));
                    } else {
                        kinematic_hardening_knots.push((temperature, value));
                    }
                }
                if name == receipts[0].property {
                    minimum_source_heat_capacity_j_kg_k =
                        minimum_source_heat_capacity_j_kg_k.min(value);
                }
                if fixed && value != property.value_si() {
                    return Err(FeedbackError::Invalid(
                        "Poisson ratio and reference density must be constant on every source-line segment",
                    ));
                }
                if !fixed && !expansion && !hardening && !kinematic && !positive(value) {
                    return Err(FeedbackError::Invalid(
                        "Young's modulus, yield and Cp must be positive at every source-line breakpoint",
                    ));
                }
            }
        }
        let source = ResolvedJ2Plasticity::try_from_material(
            bundle,
            query.hardening_basis,
            query.strain_limit,
        )
        .map_err(FeedbackError::Source)?;
        let youngs_temperature_slope_pa_k = affine_temperature_slope(
            &youngs_knots,
            state.temperature_k(),
            source.law().elastic.youngs,
        ).map_err(FeedbackError::Invalid)?;
        let hardening_source_slope = affine_temperature_slope(
            &hardening_knots,
            state.temperature_k(),
            source.hardening_source().value_si(),
        ).map_err(FeedbackError::Invalid)?;
        let hardening_temperature_slope_pa_k = match query.hardening_basis {
            HardeningPropertyBasis::PlasticStrainModulus
                | HardeningPropertyBasis::MixedPlasticStrainModuli => hardening_source_slope,
            HardeningPropertyBasis::UniaxialTangentModulus => {
                if (youngs_temperature_slope_pa_k != 0.0 || hardening_source_slope != 0.0)
                    && hardening_knots.iter().any(|(_, value)| *value != 0.0)
                {
                    return Err(FeedbackError::Invalid(
                        "temperature-dependent energetic hardening requires the native plastic-strain modulus; E*Et/(E-Et) is not this affine model",
                    ));
                }
                0.0
            }
        };
        let maximum_youngs_modulus_pa = youngs_knots
            .iter()
            .map(|(_, value)| *value)
            .fold(0.0_f64, f64::max);
        let kinematic_hardening_temperature_slope_pa_k =
            if query.hardening_basis == HardeningPropertyBasis::MixedPlasticStrainModuli {
                affine_temperature_slope(
                    &kinematic_hardening_knots,
                    state.temperature_k(),
                    source.law().kinematic_hardening,
                )
                .map_err(FeedbackError::Invalid)?
            } else {
                0.0
            };
        let expansion_per_k = source
            .material()
            .property(crate::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .map_or(0.0, |p| p.value_si());
        if expansion_knots
            .iter()
            .all(|(_, alpha)| *alpha == expansion_per_k)
        {
            expansion_knots.clear();
        }
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
            expansion_knots,
            minimum_source_heat_capacity_j_kg_k,
            youngs_temperature_slope_pa_k,
            hardening_temperature_slope_pa_k,
            kinematic_hardening_temperature_slope_pa_k,
            maximum_youngs_modulus_pa,
            reference_temperature_k: state.temperature_k(),
            overstress: None,
            voce: None,
        };
        // Affine extrema occur at the span ends; the native owner checks the
        // elastic tensor, hardening conversion and return-map denominator there.
        result.resolve_temperature(knots[0].temperature_k)?;
        result.resolve_temperature(knots[knots.len() - 1].temperature_k)?;
        result.admit_thermoelastic_capacity()?;
        Ok(result)
    }

    /// Select linear J2 overstress viscosity from this same immutable material.
    /// Its original claim is pinned and must continuously support the caloric
    /// interval with positive finite values. Temperature may change viscosity;
    /// additional rate/stress/history source axes are not silently frozen.
    /// The actual physical duration is supplied by each later step.
    pub fn with_linear_overstress(self, selection: MaterialPropertySelection)
        -> Result<Self, FeedbackError> {
        self.with_overstress(selection, OverstressSelection::Linear)
    }

    /// Select the associative power law alpha_dot=(s0/eta)*max(f/s0,0)^n.
    /// Viscosity [Pa s], stress scale [Pa] and exponent n>=1 are selected
    /// together from this card, pinned over its complete caloric interval, and
    /// refreshed at every trial temperature. No rate axis is silently frozen.
    pub fn with_power_law_overstress(self, selection: MaterialPropertySelection)
        -> Result<Self, FeedbackError> {
        self.with_overstress(selection, OverstressSelection::PowerLaw)
    }

    fn with_overstress(mut self, selection: MaterialPropertySelection, model: OverstressSelection)
        -> Result<Self, FeedbackError> {
        use super::rate::{OVERSTRESS_VISCOSITY_DIMS, OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_EXPONENT_PROPERTY};
        if self.overstress.is_some() {
            return Err(FeedbackError::Invalid("an overstress law is already selected"));
        }
        let mut roles = vec![(OVERSTRESS_VISCOSITY_PROPERTY, OVERSTRESS_VISCOSITY_DIMS)];
        if model == OverstressSelection::PowerLaw {
            roles.push((OVERSTRESS_STRESS_SCALE_PROPERTY, Dims([-1, 1, -2, 0, 0, 0])));
            roles.push((OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6])));
        }
        let requirements = roles.iter().map(|(name, dims)|
            ScalarPropertyRequirement::try_new(*name, *dims, ScalarAdmissibility::StrictlyPositive))
            .collect::<Result<Vec<_>, _>>().map_err(|e| FeedbackError::Source(e.into()))?;
        let resolved = resolve_material_state_point(self.card, &self.context,
            &requirements, selection).map_err(|e| FeedbackError::Source(e.into()))?;
        let knots = self.caloric.curve().knots();
        let lower = point_at(&self.context, &self.temperature_axis, knots[0].temperature_k)?;
        let upper = point_at(&self.context, &self.temperature_axis, knots[knots.len()-1].temperature_k)?;
        let mut additional_pins = Vec::with_capacity(roles.len());
        for property in resolved.properties() {
            admit_rate_independent_property(property, &self.temperature_axis)?;
            let name = property.requirement().name();
            let selected = property.answer().receipt.selected;
            let claim = self.card.claims().claim(selected)
                .ok_or(FeedbackError::Invalid("overstress claim left its original card"))?;
            let support = self.card.claims().query_line_typed(&claim.key, &self.temperature_axis,
                &lower, &upper, ClaimSelection::Pinned(selected))
                .map_err(|e| FeedbackError::Support { property: name.into(), reason: e.to_string() })?;
            if support.knots.iter().any(|(_, sample)| {
                let value = sample.evidence.value.value;
                !positive(value) || name == OVERSTRESS_EXPONENT_PROPERTY && value < 1.0
            }) {
                return Err(FeedbackError::Invalid(
                    "overstress coefficients must stay positive and the power exponent >= 1 throughout the source interval"));
            }
            additional_pins.push((name.to_owned(), selected));
        }
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut self.selection else {
            return Err(FeedbackError::Invalid("thermoplastic coefficients lost original claim pins"));
        };
        pins.extend(additional_pins);
        self.requirements.extend(requirements);
        self.overstress = Some(model);
        self.source = self.resolve_temperature(self.reference_temperature_k)?;
        Ok(self)
    }

    /// Whether the original n=1 rate law was explicitly selected.
    pub const fn has_linear_overstress(&self) -> bool {
        matches!(self.overstress, Some(OverstressSelection::Linear))
    }
    /// Whether the source-selected power law was explicitly selected.
    pub const fn has_power_law_overstress(&self) -> bool {
        matches!(self.overstress, Some(OverstressSelection::PowerLaw))
    }
    /// Whether constitutive updates require the actual physical interval.
    pub const fn has_overstress(&self) -> bool { self.overstress.is_some() }

    /// Admit accepted strain/history using the selected kinetic interpretation.
    /// Viscous material may carry real positive overstress into the next step.
    pub fn validate_committed_history(&self, source: &ResolvedJ2Plasticity,
        mechanical_strain: &Voigt, history: &J2State, tolerance_pa: f64)
        -> Result<(), EnergyError> {
        if self.overstress.is_some() {
            source.law().validate_overstress_history(mechanical_strain, history)
        } else if let Some(law) = self.voce_law::<Infallible>(source)
            .map_err(|_| EnergyError::Invalid("invalid retained Voce source"))? {
            law.validate_committed_history(mechanical_strain, history, tolerance_pa)
        } else {
            source.law().validate_committed_history(mechanical_strain, history, tolerance_pa)
        }
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
        admit_rate_independent_property(property, &self.temperature_axis)?;
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
        if source.law().elastic.poisson != self.source.law().elastic.poisson
            || self.kinematic_hardening_temperature_slope_pa_k == 0.0
                && source.law().kinematic_hardening != self.source.law().kinematic_hardening
            || self.youngs_temperature_slope_pa_k == 0.0
                && source.law().elastic.youngs != self.source.law().elastic.youngs
            || self.hardening_temperature_slope_pa_k == 0.0
                && source.law().hardening != self.source.law().hardening
        {
            return Err(FeedbackError::Invalid(
                "pinned coefficients left the admitted affine thermal free-energy model",
            ));
        }
        self.voce_law::<E>(&source)?;
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

    /// Attach source-defined overstress before any loading or time advance.
    /// The same candidate/transport owner supplies physical dt and finite heat.
    pub fn with_linear_overstress(mut self, selection: MaterialPropertySelection)
        -> Result<Self, FeedbackError> {
        if self.time_s != 0.0 || self.strain != [0.0; 6] || self.history.alpha != 0.0 {
            return Err(FeedbackError::Invalid("select the overstress law before advancing the specimen"));
        }
        self.binding = self.binding.with_linear_overstress(selection)?;
        self.source = self.binding.resolve_temperature(self.state.temperature_k())?;
        Ok(self)
    }

    /// Attach nonlinear power-law flow before loading. The existing uniform
    /// caloric/transport owner supplies its actual dt and receives physical heat.
    pub fn with_power_law_overstress(mut self, selection: MaterialPropertySelection)
        -> Result<Self, FeedbackError> {
        if self.time_s != 0.0 || self.strain != [0.0; 6] || self.history.alpha != 0.0 {
            return Err(FeedbackError::Invalid("select the overstress law before advancing the specimen"));
        }
        self.binding = self.binding.with_power_law_overstress(selection)?;
        self.source = self.binding.resolve_temperature(self.state.temperature_k())?;
        Ok(self)
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
        self.binding
            .validate_committed_history(
                &self.source,
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
            let (_, trial) = self.propose_at(temperature, &strain, tolerances, duration_s)?;
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
            let (source, increment) = self.propose_at(final_temperature, &strain, tolerances, duration_s)?;
            let reversible_heat_j = product(self.volume_m3, increment.reversible_heat_j_m3)?;
            let thermoelastic_storage_change_j =
                product(self.volume_m3, increment.thermoelastic_storage_change_j_m3)?;
            let mechanical_work_j = product(self.volume_m3, increment.mechanical_work_j_m3)?;
            let coefficient_storage_change_j =
                product(self.volume_m3, increment.coefficient_storage_change_j_m3)?;
            let stored_change_j = product(self.volume_m3, increment.stored_change_j_m3())?;
            let mechanical = increment.mechanical;
            let physical_heat_j =
                product(self.volume_m3, mechanical.energy.physical_dissipation_j_m3)?;
            let caloric_heat_j = finite(physical_heat_j + reversible_heat_j)?;
            let thermal_residual =
                finite(caloric_change - caloric_heat_j - proposal.external_heat_j)?;
            let combined_residual = finite(
                stored_change_j
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
                    coefficient_storage_change_j,
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
        duration_s: f64,
    ) -> Result<(ResolvedJ2Plasticity, ThermoplasticIncrement), FeedbackError<E>> {
        self.binding.increment_at_temperature(
            self.state.temperature_k(),
            temperature_k,
            &self.strain,
            &self.history,
            strain,
            tolerances,
            Some(duration_s),
        )
    }
}

fn admit_rate_independent_property<E>(
    property: &crate::state_point::ResolvedScalarProperty,
    temperature_axis: &str,
) -> Result<(), FeedbackError<E>> {
    if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
        return Err(FeedbackError::Invalid(
            "coefficients require the owned temperature and fixed external conditions; mechanical, rate and additional temperature axes need an evolving owner",
        ));
    }
    if property
        .angular_frequency_band_rad_s()
        .map_err(|e| FeedbackError::Source(e.into()))?
        .is_some()
    {
        return Err(FeedbackError::Invalid(
            "frequency-bounded data cannot define instantaneous thermoplastic coefficients",
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
