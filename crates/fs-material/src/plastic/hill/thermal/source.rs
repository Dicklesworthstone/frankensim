//! Complete pinned source paths for the original thermal Hill law.

use super::*;
use crate::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY, HILL_KINEMATIC_MODULUS_PROPERTY,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY, HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
    HILL_VOCE_SATURATION_RATE_PROPERTY, hill48_requirements, voce_requirements,
};
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, YOUNG_MODULUS_PROPERTY, affine_temperature_slope,
    pinned_temperature_line, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

/// Fixed source context and independently selected Hill directional data.
/// The initial query temperature must agree with the supplied initial enthalpy.
#[derive(Clone, Debug)]
pub struct ThermalHillSourceQuery<'a> {
    pub card: &'a MaterialCard,
    pub caloric: &'a ResolvedHeatCapacityCurve,
    pub context: QueryPoint,
    pub temperature_axis: String,
    /// E, nu, Hill reference radius, energetic H and all six Hill coefficients.
    /// Their original claim IDs remain pinned for every later temperature.
    pub selection: MaterialPropertySelection,
    /// None explicitly selects zero expansion. Some selects the canonical
    /// constant isotropic coefficient; absent directional data is never inferred.
    pub thermal_expansion: Option<MaterialPropertySelection>,
    /// Tensor strain norm ceiling in (0,0.05].
    pub strain_limit: f64,
    /// Columns are the immutable material axes in spatial coordinates.
    pub material_to_spatial: [[f64; 3]; 3],
}

impl<'a> ThermalHillSource<'a> {
    /// Admit every crossed source-line knot, preserving original Cp/density
    /// receipts. No extrapolation, changing material metric, hidden rate axis,
    /// or derivative jump in stored-energy coefficients is admitted.
    pub fn try_new(
        query: ThermalHillSourceQuery<'a>,
        initial_h: f64,
    ) -> Result<Self, ThermalHillError> {
        let invalid = ThermalHillError::Invalid;
        let curve = query.caloric.curve();
        let knots = curve.knots();
        if curve.material_card_identity() != query.card.content_hash()
            || knots.len() < 2
            || knots.len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
            || knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
            || knots
                .windows(2)
                .any(|p| p[0].temperature_k >= p[1].temperature_k)
        {
            return Err(invalid(
                "one strictly increasing, fully solid caloric chart required",
            ));
        }
        let phase = curve.state_at_specific_enthalpy(initial_h)?;
        let reference_temperature_k = phase.temperature_k();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let declared = query
            .context
            .axis_quantities()
            .get(&query.temperature_axis)
            .copied();
        let supported = declared == Some(absolute)
            || query.temperature_axis == "T"
                && (declared.is_none()
                    || declared == Some(QuantitySpec::dimensional(Temperature::DIMS)));
        if !supported
            || query.context.axes().get(&query.temperature_axis).copied()
                != Some(reference_temperature_k)
        {
            return Err(invalid(
                "initial temperature must match enthalpy; renamed axes require absolute-temperature semantics",
            ));
        }
        let thermal = ThermalMaterialPoint {
            material: query.card,
            context: &query.context,
            temperature_axis: &query.temperature_axis,
            liquid_fraction_axis: None,
        };
        thermal.at(phase)?;
        let mut requirements = hill48_requirements()?;
        let initial = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            query.selection,
        )?;
        let mut pins: Vec<_> = initial
            .properties()
            .iter()
            .map(|p| {
                (
                    p.requirement().name().to_owned(),
                    p.answer().receipt.selected,
                )
            })
            .collect();
        let receipts = query.caloric.receipts();
        if receipts.len() < 4 || receipts.len() % 2 != 0 {
            return Err(invalid(
                "paired original Cp and reference-density receipts required",
            ));
        }
        for (index, dims) in [Dims([2, 0, -2, -1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])]
            .into_iter()
            .enumerate()
        {
            let receipt = &receipts[index];
            let claim = query.card.claims().claim(receipt.selected).ok_or(invalid(
                "caloric claim is absent from the shared material card",
            ))?;
            if claim.key.quantity() != QuantitySpec::dimensional(dims)
                || requirements.iter().any(|r| r.name() == claim.key.name())
            {
                return Err(invalid(
                    "distinct mechanical and caloric roles with exact SI units required",
                ));
            }
            requirements.push(ScalarPropertyRequirement::try_with_key(
                &claim.key,
                ScalarAdmissibility::StrictlyPositive,
            )?);
            pins.push((claim.key.name().to_owned(), receipt.selected));
        }
        if let Some(selection) = query.thermal_expansion {
            let requirement = ScalarPropertyRequirement::try_new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite,
            )?;
            let selected = resolve_material_state_point(
                query.card,
                &query.context,
                std::slice::from_ref(&requirement),
                selection,
            )?;
            pins.push((
                requirement.name().to_owned(),
                selected.properties()[0].answer().receipt.selected,
            ));
            requirements.push(requirement);
        }
        let selection = MaterialPropertySelection::PinnedByProperty(pins);
        let bundle = resolve_material_state_point(
            query.card,
            &query.context,
            &requirements,
            selection.clone(),
        )?;
        if receipts.iter().any(|receipt| {
            receipt.axis_quantities != *bundle.axis_quantities()
                || !same_context(
                    &receipt.query_point,
                    bundle.query_point(),
                    &query.temperature_axis,
                )
        }) || receipts.chunks_exact(2).any(|pair| {
            pair[0].selected != receipts[0].selected || pair[1].selected != receipts[1].selected
        }) {
            return Err(invalid(
                "mechanical and caloric receipts must retain the same fixed context and Cp/density claims",
            ));
        }
        let density_name = &receipts[1].property;
        let density = bundle
            .property(density_name)
            .ok_or(invalid("original density missing from source bundle"))?
            .value_si();
        if knots.iter().any(|k| k.bulk_density_kg_m3 != density) {
            return Err(invalid(
                "reference density must be constant throughout the caloric chart",
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
        let mut youngs_slope = 0.0;
        let mut hardening_slope = 0.0;
        let mut maximum_youngs = 0.0_f64;
        let mut minimum_cp = f64::INFINITY;
        for property in bundle.properties() {
            thermal.admit_static(property)?;
            if !property.has_only_fixed_external_conditions(Some(&query.temperature_axis), false)
                || property.angular_frequency_band_rad_s()?.is_some()
            {
                return Err(invalid(
                    "only the owned temperature and fixed external conditions may drive instantaneous Hill coefficients",
                ));
            }
            let name = property.requirement().name();
            let claim = query
                .card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(invalid("selected Hill coefficient left its source card"))?;
            let path = pinned_temperature_line(
                query.card,
                &claim.key,
                &query.temperature_axis,
                &lower,
                &upper,
                property.answer().receipt.selected,
            )?;
            let fixed = name == POISSON_RATIO_PROPERTY
                || name == density_name
                || HILL_COEFFICIENT_PROPERTIES.contains(&name)
                || name == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY;
            for &(_, value) in &path {
                if !value.is_finite() || fixed && value != property.value_si() {
                    return Err(invalid(
                        "Poisson ratio, Hill metric, density and selected expansion must remain fixed over every source-line cell",
                    ));
                }
                if name == HILL_HARDENING_MODULUS_PROPERTY && value < 0.0
                    || (name == YOUNG_MODULUS_PROPERTY
                        || name == HILL_REFERENCE_YIELD_STRESS_PROPERTY
                        || name == receipts[0].property)
                        && value <= 0.0
                {
                    return Err(invalid(
                        "E, reference yield and Cp must remain positive; energetic hardening must remain nonnegative",
                    ));
                }
                if name == YOUNG_MODULUS_PROPERTY {
                    maximum_youngs = maximum_youngs.max(value);
                }
                if name == receipts[0].property {
                    minimum_cp = minimum_cp.min(value);
                }
            }
            if name == YOUNG_MODULUS_PROPERTY || name == HILL_HARDENING_MODULUS_PROPERTY {
                let slope =
                    affine_temperature_slope(&path, reference_temperature_k, property.value_si())
                        .map_err(invalid)?;
                if name == YOUNG_MODULUS_PROPERTY {
                    youngs_slope = slope;
                } else {
                    hardening_slope = slope;
                }
            }
        }
        let reference = ResolvedHill48::try_from_material(
            bundle,
            query.strain_limit,
            query.material_to_spatial,
        )?;
        let expansion_per_k = reference
            .material()
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .map_or(0.0, |p| p.value_si());
        // Affine stored moduli have no second-T derivative at fixed history.
        // Bound the expansion/elastic contribution for ALL admitted elastic
        // tensors, including plastic unloads, using |tr(e)| <= sqrt(3)||e||.
        let denominator = 3.0 * (1.0 - 2.0 * reference.local_law().elastic().poisson);
        let bulk_max = maximum_youngs / denominator;
        let bulk_slope = youngs_slope / denominator;
        let correction_per_k = finite(
            product(9.0 * bulk_max, product(expansion_per_k, expansion_per_k)?)?
                + product(
                    6.0 * fs_math::det::sqrt(3.0),
                    product(
                        product(expansion_per_k, bulk_slope)?.abs(),
                        query.strain_limit,
                    )?,
                )?,
        )?;
        let correction = product(knots[knots.len() - 1].temperature_k, correction_per_k)?;
        for index in 0..knots.len() - 1 {
            minimum_cp = minimum_cp.min(curve.sensible_heat_capacity_bounds_j_kg_k(index)?[0]);
        }
        if !minimum_cp.is_finite()
            || minimum_cp <= 0.0
            || finite(product(density, minimum_cp)? - correction)? <= 0.0
        {
            return Err(invalid(
                "combined caloric and mechanical heat capacity must stay positive over the complete source span",
            ));
        }
        let result = Self {
            card: query.card,
            caloric: query.caloric,
            context: query.context,
            temperature_axis: query.temperature_axis,
            requirements,
            selection,
            reference,
            reference_temperature_k,
            density_kg_m3: density,
            expansion_per_k,
            youngs_slope_pa_k: youngs_slope,
            hardening_slope_pa_k: hardening_slope,
            kinematic_slope_pa_k: None,
            voce_slope_pa_k: None,
        };
        // The native owner checks elasticity and the return-map denominator at
        // both affine extrema. Convex combinations keep the fixed Hill metric.
        result.resolve_temperature(knots[0].temperature_k)?;
        result.resolve_temperature(knots[knots.len() - 1].temperature_k)?;
        Ok(result)
    }

    /// Select independently sourced energetic Prager hardening before retaining
    /// any material histories. The source is C [Pa] in X=C*eps_p and
    /// psi_k=C*eps_p:eps_p/2, not the differently normalized J2 H_k modulus.
    /// C must be affine and nonnegative throughout the original caloric span;
    /// zero is admitted without dividing by C or changing the strain-like state.
    /// The fixed-history second temperature derivative is zero, so this adds
    /// no hidden hardening term to the admitted effective heat capacity.
    pub fn with_linear_kinematic_hardening(
        mut self,
        selection: MaterialPropertySelection,
    ) -> Result<Self, ThermalHillError> {
        let invalid = ThermalHillError::Invalid;
        if self.kinematic_slope_pa_k.is_some()
            || self
                .requirements
                .iter()
                .any(|r| r.name() == HILL_KINEMATIC_MODULUS_PROPERTY)
        {
            return Err(invalid(
                "Hill kinematic hardening must be selected once as a distinct material role",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            HILL_KINEMATIC_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let resolved = resolve_material_state_point(
            self.card,
            &self.context,
            std::slice::from_ref(&requirement),
            selection,
        )?;
        let property = &resolved.properties()[0];
        self.thermal_material_point().admit_static(property)?;
        if !property.has_only_fixed_external_conditions(Some(&self.temperature_axis), false)
            || property.angular_frequency_band_rad_s()?.is_some()
        {
            return Err(invalid(
                "kinematic hardening requires only the owned temperature and fixed external conditions, without rate or frequency dependence",
            ));
        }
        let selected = property.answer().receipt.selected;
        let claim = self.card.claims().claim(selected).ok_or(invalid(
            "selected Hill kinematic source is absent from the original material card",
        ))?;
        let knots = self.caloric.curve().knots();
        let lower = point_at(
            &self.context,
            &self.temperature_axis,
            knots[0].temperature_k,
        )?;
        let upper = point_at(
            &self.context,
            &self.temperature_axis,
            knots[knots.len() - 1].temperature_k,
        )?;
        let path = pinned_temperature_line(
            self.card,
            &claim.key,
            &self.temperature_axis,
            &lower,
            &upper,
            selected,
        )?;
        if path
            .iter()
            .any(|(_, value)| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "Hill kinematic modulus must remain finite and nonnegative on every source-line cell",
            ));
        }
        let slope =
            affine_temperature_slope(&path, self.reference_temperature_k, property.value_si())
                .map_err(invalid)?;
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut self.selection else {
            return Err(invalid("thermal Hill source lost its original claim pins"));
        };
        pins.push((HILL_KINEMATIC_MODULUS_PROPERTY.to_owned(), selected));
        self.requirements.push(requirement);
        let bundle = resolve_material_state_point(
            self.card,
            &self.context,
            &self.requirements,
            self.selection.clone(),
        )?;
        self.kinematic_slope_pa_k = Some(slope);
        self.reference = self.apply_selected_hardening(ResolvedHill48::try_from_material(
            bundle,
            self.reference.local_law().elastic().strain_limit,
            *self.reference.law().material_to_spatial(),
        )?)?;
        // Endpoint checks use the original scalar return denominator. Affine
        // C and E cannot hide an interior maximum or minimum between them.
        self.resolve_temperature(knots[0].temperature_k)?;
        self.resolve_temperature(knots[knots.len() - 1].temperature_k)?;
        Ok(self)
    }

    /// Select independent Hill Q(T) and b before retaining specimen history.
    /// Q is nonnegative and globally affine over the complete caloric span;
    /// b is positive and fixed per Hill work-conjugate alpha. Every source-line
    /// knot and surface crossing is admitted, with the original claims pinned.
    /// The extra free energy Q(T)*g(alpha) contributes entropy -Q'*g(alpha),
    /// but no frozen-history heat-capacity correction for affine Q and fixed b.
    pub fn with_voce_isotropic_hardening(
        mut self,
        selection: MaterialPropertySelection,
    ) -> Result<Self, ThermalHillError> {
        let invalid = ThermalHillError::Invalid;
        if self.voce_slope_pa_k.is_some()
            || self.requirements.iter().any(|r| {
                r.name() == HILL_VOCE_SATURATION_INCREMENT_PROPERTY
                    || r.name() == HILL_VOCE_SATURATION_RATE_PROPERTY
            })
        {
            return Err(invalid(
                "Hill Voce hardening must be selected once as independent material roles",
            ));
        }
        let requirements = voce_requirements()?;
        let resolved =
            resolve_material_state_point(self.card, &self.context, &requirements, selection)?;
        let saturation = resolved
            .property(HILL_VOCE_SATURATION_INCREMENT_PROPERTY)
            .expect("resolved Hill Voce requirements retain Q")
            .value_si();
        let rate = resolved
            .property(HILL_VOCE_SATURATION_RATE_PROPERTY)
            .expect("resolved Hill Voce requirements retain b")
            .value_si();
        let knots = self.caloric.curve().knots();
        let lower = point_at(
            &self.context,
            &self.temperature_axis,
            knots[0].temperature_k,
        )?;
        let upper = point_at(
            &self.context,
            &self.temperature_axis,
            knots[knots.len() - 1].temperature_k,
        )?;
        let mut saturation_path = Vec::new();
        let mut additional_pins = Vec::with_capacity(2);
        for property in resolved.properties() {
            self.thermal_material_point().admit_static(property)?;
            if !property.has_only_fixed_external_conditions(Some(&self.temperature_axis), false)
                || property.angular_frequency_band_rad_s()?.is_some()
            {
                return Err(invalid(
                    "Hill Voce coefficients require the owned temperature and fixed external conditions, without rate or frequency dependence",
                ));
            }
            let name = property.requirement().name();
            let selected = property.answer().receipt.selected;
            let claim = self.card.claims().claim(selected).ok_or(invalid(
                "selected Hill Voce source is absent from its original material card",
            ))?;
            let path = pinned_temperature_line(
                self.card,
                &claim.key,
                &self.temperature_axis,
                &lower,
                &upper,
                selected,
            )?;
            if name == HILL_VOCE_SATURATION_INCREMENT_PROPERTY {
                if path
                    .iter()
                    .any(|(_, value)| !value.is_finite() || *value < 0.0)
                {
                    return Err(invalid(
                        "Hill Voce Q must remain finite and nonnegative on every source-line cell",
                    ));
                }
                saturation_path = path;
            } else if path
                .iter()
                .any(|(_, value)| !value.is_finite() || *value <= 0.0 || *value != rate)
            {
                return Err(invalid(
                    "Hill Voce b must remain positive and constant on every source-line cell",
                ));
            }
            additional_pins.push((name.to_owned(), selected));
        }
        let slope =
            affine_temperature_slope(&saturation_path, self.reference_temperature_k, saturation)
                .map_err(invalid)?;
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut self.selection else {
            return Err(invalid("thermal Hill source lost its original claim pins"));
        };
        pins.extend(additional_pins);
        self.requirements.extend(requirements);
        self.voce_slope_pa_k = Some(slope);
        let bundle = resolve_material_state_point(
            self.card,
            &self.context,
            &self.requirements,
            self.selection.clone(),
        )?;
        self.reference = self.apply_selected_hardening(ResolvedHill48::try_from_material(
            bundle,
            self.reference.local_law().elastic().strain_limit,
            *self.reference.law().material_to_spatial(),
        )?)?;
        // E, H, C and Q are affine; the native owner checks representable
        // combined return moduli at both extrema, including exactly Q=0.
        self.resolve_temperature(knots[0].temperature_k)?;
        self.resolve_temperature(knots[knots.len() - 1].temperature_k)?;
        Ok(self)
    }

    /// Complete initial source bundle, including original caloric receipts.
    pub const fn reference(&self) -> &ResolvedHill48 {
        &self.reference
    }

    pub(super) fn resolve_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedHill48, ThermalHillError> {
        let knots = self.caloric.curve().knots();
        if !temperature_k.is_finite()
            || temperature_k < knots[0].temperature_k
            || temperature_k > knots[knots.len() - 1].temperature_k
        {
            return Err(ThermalHillError::Invalid(
                "temperature is outside the original caloric source span",
            ));
        }
        if temperature_k == self.reference_temperature_k {
            return Ok(self.reference.clone());
        }
        let point = point_at(&self.context, &self.temperature_axis, temperature_k)?;
        let bundle = resolve_material_state_point(
            self.card,
            &point,
            &self.requirements,
            self.selection.clone(),
        )?;
        let resolved = ResolvedHill48::try_from_material(
            bundle,
            self.reference.local_law().elastic().strain_limit,
            *self.reference.law().material_to_spatial(),
        )?;
        self.apply_selected_hardening(resolved)
    }

    // Rebuild both explicitly selected owners after every source refresh.
    // Selection order cannot drop previously admitted Q/b or Prager C.
    fn apply_selected_hardening(
        &self,
        mut resolved: ResolvedHill48,
    ) -> Result<ResolvedHill48, ThermalHillError> {
        if self.kinematic_slope_pa_k.is_some() {
            resolved = resolved.with_linear_kinematic_hardening()?;
        }
        if self.voce_slope_pa_k.is_some() {
            resolved = resolved.with_voce_isotropic_hardening()?;
        }
        Ok(resolved)
    }
}

fn point_at(
    context: &QueryPoint,
    axis: &str,
    temperature_k: f64,
) -> Result<QueryPoint, ThermalHillError> {
    context
        .clone()
        .with(axis, temperature_k)
        .map_err(|e| ThermalHillError::Transport(ThermalMaterialError::Coordinate(e.to_string())))
}
fn same_context(a: &[(String, f64)], b: &[(String, f64)], axis: &str) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|((ka, va), (kb, vb))| ka == kb && (ka == axis || va.to_bits() == vb.to_bits()))
}
