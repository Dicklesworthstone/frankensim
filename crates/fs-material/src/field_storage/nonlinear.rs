//! Reversible isotropic H(|B|) and its energy, from one immutable source path.
//!
//! This is a single-valued, zero-remanence law, not a hysteresis loop or an
//! inference from a material name. The canonical matdb line query establishes
//! continuous support and pins the same claim over the complete field range.
//! Its positive segment slopes make both eigenvalues of dH_vector/dB_vector
//! positive. Energy is integral_0^B H(b) db; B*H/2 is wrong after saturation.

use super::{
    FieldStorageError,
    permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
    resolve_material_state_point,
};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyLineSlice, QueryPoint};
use fs_qty::QuantitySpec;

mod convexity;
/// Explicit source-scaled H(B,T) and its reversible thermal energy.
pub mod thermal;

/// Magnitude H [A/m] of an explicitly reversible, isotropic, zero-remanence law.
/// A major hysteresis loop, coercivity or secant permeability is not this key.
pub const REVERSIBLE_FIELD_STRENGTH_PROPERTY: &str = "reversible_magnetic_field_strength";

/// One pinned physical context and continuously supported H(|B|) source path.
/// Temperature, pressure, composition and all other external state stay fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct ReversibleMagneticMaterial {
    card: MaterialCard,
    material: ResolvedMaterialStatePoint,
    path: PropertyLineSlice,
    energy_at_knots_j_m3: Vec<f64>,
    coenergy_at_knots_j_m3: Vec<f64>,
    slopes_a_m_t: Vec<f64>,
    // Only the explicit separable thermal constructor creates this binding.
    // The original path/receipts remain H_ref, never synthesized H(B,T) data.
    temperature_scale: Option<thermal::TemperatureScale>,
}

impl ReversibleMagneticMaterial {
    /// Admit an electrical/material refresh without changing this magnetic law.
    /// The same card, selected H(B) source and coordinate schema remain pinned.
    /// Only the owned absolute-temperature coordinate may change, and an H(B)
    /// source that depends on that temperature refuses even coincident endpoints.
    /// Additional electrical properties in the refreshed bundle are retained by
    /// the caller. The complete original field range must still be supported.
    pub fn admit_fixed_material_refresh(
        &self,
        refreshed: &ResolvedMaterialStatePoint,
    ) -> Result<(), FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        if let Some(scale) = &self.temperature_scale {
            let original = self
                .material
                .property(&scale.property)
                .ok_or(invalid("prepared magnetic scale source is missing"))?;
            let actual = refreshed.property(&scale.property).ok_or(invalid(
                "refreshed conductor lost the magnetic scale source",
            ))?;
            if original.requirement() != actual.requirement()
                || original.answer().receipt.selected != actual.answer().receipt.selected
                || original.source_axes() != actual.source_axes()
                || actual.value_si().to_bits() != scale.value.to_bits()
            {
                return Err(invalid(
                    "fixed magnetic refresh cannot change its source scale without parameter work and entropy",
                ));
            }
        }
        let original = self
            .material
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .expect("prepared reversible field source");
        let property = refreshed
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .ok_or(invalid(
                "refreshed conductor must retain the reversible magnetic source",
            ))?;
        if self.material.card_identity() != refreshed.card_identity()
            || self.material.material() != refreshed.material()
            || original.requirement() != property.requirement()
            || original.source_axes() != property.source_axes()
            || original.answer().receipt.selected != property.answer().receipt.selected
        {
            return Err(invalid(
                "nonlinear magnetic refresh changed the original card or source",
            ));
        }
        // Rebuild just the original H(B) bundle so common field continuity
        // admission does not discard or miscompare extra conductor properties.
        let query = refreshed.query();
        let checked = self.material.at_query_point(&self.card, &query)?;
        super::admit_source_temperature_refresh(&self.material, &checked)?;
        if property.value_si().to_bits()
            != checked
                .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
                .expect("rechecked reversible source")
                .value_si()
                .to_bits()
            || self
                .material
                .query_point()
                .iter()
                .zip(checked.query_point())
                .any(|((axis, a), (_, b))| {
                    a.to_bits() != b.to_bits()
                        && original.source_axes().iter().any(|(name, _)| name == axis)
                })
        {
            return Err(invalid(
                "temperature-dependent H(B) requires magnetic parameter work and entropy",
            ));
        }
        let upper = query
            .clone()
            .with(&self.path.axis, self.maximum_flux_density_t())
            .map_err(|_| invalid("refreshed magnetic field range is unrepresentable"))?;
        let selected = original.answer().receipt.selected;
        let claim = self
            .card
            .claims()
            .claim(selected)
            .expect("prepared source claim");
        let path = self
            .card
            .claims()
            .query_line_typed(
                &claim.key,
                &self.path.axis,
                &query,
                &upper,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: REVERSIBLE_FIELD_STRENGTH_PROPERTY.to_owned(),
                source: Box::new(source),
            })?;
        if path.knots.len() != self.path.knots.len()
            || path
                .knots
                .iter()
                .zip(&self.path.knots)
                .any(|((b, h), (old_b, old_h))| {
                    b.to_bits() != old_b.to_bits()
                        || h.evidence.value.value.to_bits() != old_h.evidence.value.value.to_bits()
                        || h.receipt.selected != selected
                })
        {
            return Err(invalid(
                "refreshed magnetic range changes the prepared constitutive law",
            ));
        }
        Ok(())
    }

    /// Admit an explicit separable thermal refresh of this original curve.
    /// Reference path, scale claim, affine law, physical temperature range and
    /// every other source role remain fixed. Native owners must report the
    /// resulting parameter work; a finite thermal owner must close field entropy.
    pub fn admit_temperature_refresh(&self, after: &Self) -> Result<(), FieldStorageError> {
        if self == after {
            return Ok(());
        }
        let (Some(a), Some(b)) = (&self.temperature_scale, &after.temperature_scale) else {
            return Err(FieldStorageError::Invalid(
                "nonlinear magnetic temperature refresh requires an explicit separable source law",
            ));
        };
        if a.property != b.property
            || a.slope_per_k.to_bits() != b.slope_per_k.to_bits()
            || a.reference_temperature_k.to_bits() != b.reference_temperature_k.to_bits()
            || a.temperature_range_k != b.temperature_range_k
            || self.path != after.path
            || self.card.content_hash() != after.card.content_hash()
        {
            return Err(FieldStorageError::Invalid(
                "separable magnetic refresh changed its original reference or scale law",
            ));
        }
        super::admit_source_temperature_refresh(&self.material, &after.material)
    }

    /// Nonnegative backward-Euler magnetic loss density [J/m3],
    /// H(after).(after-before) - (psi(after)-psi(before)).
    /// This is numerical dissipation, NEVER constitutive Joule heat.
    /// Integrating original positive slopes and angular convexity avoids
    /// cancellation when a small step changes a large stored field energy.
    pub fn bregman_density_j_m3(
        &self,
        before_flux_density_t: [f64; 2],
        after_flux_density_t: [f64; 2],
    ) -> Result<f64, FieldStorageError> {
        self.bregman_density(before_flux_density_t, after_flux_density_t)
    }

    /// The same source-exact backward-Euler convexity loss for a spatial field.
    /// Includes changes of both magnitude and direction, including rotations
    /// crossing the original saturation segments. This is numerical loss, not
    /// hysteretic or conductor heating; it must remain a separate work account.
    pub fn bregman_density_3d_j_m3(
        &self,
        before_flux_density_t: [f64; 3],
        after_flux_density_t: [f64; 3],
    ) -> Result<f64, FieldStorageError> {
        self.bregman_density(before_flux_density_t, after_flux_density_t)
    }

    fn bregman_density<const D: usize>(
        &self,
        before_flux_density_t: [f64; D],
        after_flux_density_t: [f64; D],
    ) -> Result<f64, FieldStorageError> {
        let before_magnitude = before_flux_density_t
            .iter()
            .fold(0.0_f64, |m, x| m.hypot(*x));
        let after_magnitude = after_flux_density_t
            .iter()
            .fold(0.0_f64, |m, x| m.hypot(*x));
        self.evaluate(before_magnitude)?;
        let after = self.evaluate(after_magnitude)?;
        let first = self.segment_at_flux_density(before_magnitude.min(after_magnitude));
        let last = self.segment_at_flux_density(before_magnitude.max(after_magnitude));
        convexity::bregman_density(
            before_flux_density_t,
            after_flux_density_t,
            after.field_strength_a_m,
            self.path.knots[first..=last + 1]
                .windows(2)
                .zip(&self.slopes_a_m_t[first..=last])
                .map(|(pair, &slope)| (pair[0].0, pair[1].0, slope * self.field_scale())),
        )
        .ok_or(FieldStorageError::Invalid(
            "magnetic numerical dissipation is unrepresentable",
        ))
    }
}

/// Response of the SAME source law at one admitted flux-density magnitude.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReversibleMagneticResponse {
    /// Field-strength magnitude [A/m].
    pub field_strength_a_m: f64,
    /// H/B [m/H], using the first segment's limit at zero.
    pub secant_reluctivity_m_h: f64,
    /// dH/dB [m/H]. Interior knots use their right-hand segment.
    pub differential_reluctivity_m_h: f64,
    /// Magnetic energy integral_0^B H(b) db [J/m³], relative to B=0.
    pub energy_density_j_m3: f64,
    /// Legendre dual integral_0^H B(h) dh [J/m³], relative to H=0.
    /// Evaluated directly to avoid cancellation in the equivalent B*H-energy.
    pub coenergy_density_j_m3: f64,
}

impl ReversibleMagneticMaterial {
    /// Compile one source-selected curve or fixed-context surface slice.
    ///
    /// `flux_density_axis` must already have explicit tesla dimensions in the
    /// query. Only it varies, from zero to `maximum_flux_density_t`; the entire
    /// span must have continuous canonical source support. H(0)=0 and strictly
    /// increasing H are required. A changed material/context creates a new law;
    /// this operation does not transfer a live field or pay parameter work.
    pub fn resolve(
        card: &MaterialCard,
        context: &QueryPoint,
        flux_density_axis: &str,
        maximum_flux_density_t: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        if context.axis_quantities().get(flux_density_axis)
            != Some(&QuantitySpec::dimensional(FLUX_DENSITY_DIMS))
            || !context.axes().contains_key(flux_density_axis)
            || !maximum_flux_density_t.is_finite()
            || maximum_flux_density_t <= 0.0
        {
            return Err(invalid(
                "reversible magnetic law needs an existing tesla coordinate and positive finite field range",
            ));
        }
        let point = |b| {
            context
                .clone()
                .with(flux_density_axis, b)
                .map_err(|_| invalid("magnetic field query is unrepresentable"))
        };
        let lower = point(0.0)?;
        let upper = point(maximum_flux_density_t)?;
        let requirement = ScalarPropertyRequirement::try_new(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let material = resolve_material_state_point(card, &lower, &[requirement], selection)?;
        let property = material
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .ok_or(invalid("reversible field source is missing"))?;
        let selected = property.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or(invalid("selected reversible field claim is absent"))?;
        if !property
            .source_axes()
            .iter()
            .any(|(axis, dims)| axis == flux_density_axis && *dims == FLUX_DENSITY_DIMS)
        {
            return Err(invalid(
                "reversible H(B) source must depend on the owned flux-density axis",
            ));
        }
        // Field, frequency, strain and rate conditions cannot masquerade as
        // fixed thermodynamic context, including applicability-only bounds.
        let external = |axis: &str| {
            is_fixed_external_coordinate(
                axis,
                lower.axis_quantities().get(axis).copied(),
                lower.axes().get(axis).copied(),
            )
        };
        if property.source_axes().iter().any(|(axis, dims)| {
            if axis == flux_density_axis {
                *dims != FLUX_DENSITY_DIMS
            } else {
                !external(axis)
                    || lower
                        .axis_quantities()
                        .get(axis)
                        .is_none_or(|quantity| quantity.dims() != *dims)
            }
        }) || claim
            .validity
            .bounds()
            .keys()
            .any(|axis| axis != flux_density_axis && !external(axis))
        {
            return Err(invalid(
                "reversible magnetic law requires one owned field axis and fixed external conditions",
            ));
        }
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                flux_density_axis,
                &lower,
                &upper,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: REVERSIBLE_FIELD_STRENGTH_PROPERTY.to_owned(),
                source: Box::new(source),
            })?;
        if path.knots.len() < 2
            || path.knots[0].0 != 0.0
            || path.knots[0].1.evidence.value.value != 0.0
            || path.knots.last().map(|k| k.0) != Some(maximum_flux_density_t)
        {
            return Err(invalid(
                "reversible magnetic path must include H(0)=0 and the complete declared range",
            ));
        }
        let mut energy_at_knots_j_m3 = vec![0.0];
        let mut coenergy_at_knots_j_m3 = vec![0.0];
        let mut slopes_a_m_t = Vec::with_capacity(path.knots.len() - 1);
        for pair in path.knots.windows(2) {
            let width = pair[1].0 - pair[0].0;
            let a = pair[0].1.evidence.value.value;
            let b = pair[1].1.evidence.value.value;
            let slope = (b - a) / width;
            let increment = (0.5 * a + 0.5 * b) * width;
            let energy = energy_at_knots_j_m3.last().copied().unwrap_or(0.0) + increment;
            // Integrate B dH directly. Subtracting energy from B*H loses the
            // entire coenergy of nearly flat, high-susceptibility H(B) spans.
            let coenergy_increment = (0.5 * pair[0].0 + 0.5 * pair[1].0) * (b - a);
            let coenergy =
                coenergy_at_knots_j_m3.last().copied().unwrap_or(0.0) + coenergy_increment;
            if !width.is_finite()
                || width <= 0.0
                || !slope.is_finite()
                || slope <= 0.0
                || !increment.is_finite()
                || increment <= 0.0
                || !energy.is_finite()
                || !coenergy_increment.is_finite()
                || coenergy_increment <= 0.0
                || !coenergy.is_finite()
            {
                return Err(invalid(
                    "reversible magnetic source must have representably positive differential reluctivity and finite energy",
                ));
            }
            slopes_a_m_t.push(slope);
            energy_at_knots_j_m3.push(energy);
            coenergy_at_knots_j_m3.push(coenergy);
        }
        Ok(Self {
            card: card.clone(),
            material,
            path,
            energy_at_knots_j_m3,
            coenergy_at_knots_j_m3,
            slopes_a_m_t,
            temperature_scale: None,
        })
    }

    /// Original zero-field material snapshot and source usage receipt.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Canonical reference-curve answers at every source breakpoint. An
    /// explicitly thermal-scaled law retains these original H_ref receipts;
    /// `evaluate` applies its separately sourced `field_scale` without inventing
    /// new measured H(B,T) answers.
    pub fn source_path(&self) -> &PropertyLineSlice {
        &self.path
    }
    /// Maximum admitted flux-density magnitude [T]. Never extrapolated.
    pub fn maximum_flux_density_t(&self) -> f64 {
        self.path.knots.last().unwrap().0
    }

    /// Maximum admitted field-strength magnitude [A/m] on the same source path.
    pub fn maximum_field_strength_a_m(&self) -> f64 {
        self.field_scale() * self.path.knots.last().unwrap().1.evidence.value.value
    }

    /// Positive source scale of an explicitly prepared separable thermal law;
    /// ordinary reversible curves retain exactly one.
    pub fn field_scale(&self) -> f64 {
        self.temperature_scale
            .as_ref()
            .map_or(1.0, |scale| scale.value)
    }

    /// Derivative of an explicitly declared separable scale, if present.
    /// A nonzero value carries field entropy even at fixed temperature when B
    /// changes, so finite thermal consumers require its caloric receiver.
    pub fn field_scale_slope_per_k(&self) -> Option<f64> {
        self.temperature_scale
            .as_ref()
            .map(|scale| scale.slope_per_k)
    }

    /// Invert the SAME strictly increasing source curve, returning |B| [T] for
    /// |H| [A/m]. This is a segment solve, not an independently fitted B(H) law.
    /// Neither end of the original field/source range is extrapolated.
    pub fn flux_density_at_field_strength(
        &self,
        field_strength_a_m: f64,
    ) -> Result<f64, FieldStorageError> {
        if !field_strength_a_m.is_finite()
            || field_strength_a_m < 0.0
            || field_strength_a_m > self.maximum_field_strength_a_m()
        {
            return Err(FieldStorageError::Invalid(
                "magnetic field strength leaves the source-backed reversible range",
            ));
        }
        if field_strength_a_m == self.maximum_field_strength_a_m() {
            return Ok(self.maximum_flux_density_t());
        }
        let reference_field = field_strength_a_m / self.field_scale();
        if !reference_field.is_finite() || field_strength_a_m > 0.0 && reference_field == 0.0 {
            return Err(FieldStorageError::Invalid(
                "scaled inverse magnetic response is unrepresentable",
            ));
        }
        let field_strength_a_m = reference_field;
        let segment = self
            .path
            .knots
            .partition_point(|(_, answer)| answer.evidence.value.value <= field_strength_a_m)
            .saturating_sub(1)
            .min(self.slopes_a_m_t.len() - 1);
        let (b0, ref answer) = self.path.knots[segment];
        let b =
            (field_strength_a_m - answer.evidence.value.value) / self.slopes_a_m_t[segment] + b0;
        if !b.is_finite() || (field_strength_a_m > 0.0 && b <= 0.0) {
            return Err(FieldStorageError::Invalid(
                "inverse reversible magnetic response is unrepresentable",
            ));
        }
        // The admitted H is inside this exact segment. Limit only rounding of
        // the affine inverse at its endpoints, never an out-of-domain input.
        Ok(b.clamp(b0, self.path.knots[segment + 1].0))
    }

    fn segment_at_flux_density(&self, flux_density_t: f64) -> usize {
        self.path
            .knots
            .partition_point(|(b, _)| *b <= flux_density_t)
            .saturating_sub(1)
            .min(self.slopes_a_m_t.len() - 1)
    }

    /// Evaluate the compiled pinned source path and its piecewise-linear
    /// derivative/integrals without querying or allocating in the field loop.
    /// The canonical line admission already proves continuous support and
    /// retains every source breakpoint, including fixed-context surface cuts.
    /// Binary search costs O(log knots); no independent permeability fit or
    /// approximation to the admitted piecewise-affine law is introduced.
    pub fn evaluate(
        &self,
        flux_density_t: f64,
    ) -> Result<ReversibleMagneticResponse, FieldStorageError> {
        let mut response = self.reference_response(flux_density_t)?;
        let scale = self.field_scale();
        if scale == 1.0 {
            return Ok(response);
        }
        response.field_strength_a_m *= scale;
        response.secant_reluctivity_m_h *= scale;
        response.differential_reluctivity_m_h *= scale;
        response.energy_density_j_m3 *= scale;
        response.coenergy_density_j_m3 *= scale;
        if [
            response.field_strength_a_m,
            response.secant_reluctivity_m_h,
            response.differential_reluctivity_m_h,
            response.energy_density_j_m3,
            response.coenergy_density_j_m3,
        ]
        .iter()
        .any(|value| !value.is_finite())
            || response.secant_reluctivity_m_h <= 0.0
            || response.differential_reluctivity_m_h <= 0.0
            || flux_density_t > 0.0
                && (response.field_strength_a_m <= 0.0
                    || response.energy_density_j_m3 <= 0.0
                    || response.coenergy_density_j_m3 <= 0.0)
        {
            return Err(FieldStorageError::Invalid(
                "scaled magnetic response or energy is unrepresentable",
            ));
        }
        Ok(response)
    }

    fn reference_response(
        &self,
        flux_density_t: f64,
    ) -> Result<ReversibleMagneticResponse, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        if !flux_density_t.is_finite()
            || flux_density_t < 0.0
            || flux_density_t > self.maximum_flux_density_t()
        {
            return Err(invalid(
                "magnetic flux density leaves the source-backed reversible range",
            ));
        }
        let segment = self.segment_at_flux_density(flux_density_t);
        let (b0, ref answer) = self.path.knots[segment];
        let h0 = answer.evidence.value.value;
        let slope = self.slopes_a_m_t[segment];
        let width = flux_density_t - b0;
        let (h, energy, coenergy) = if flux_density_t == self.maximum_flux_density_t() {
            // Preserve canonical endpoint values and the prepared integrals.
            (
                self.path.knots.last().unwrap().1.evidence.value.value,
                *self.energy_at_knots_j_m3.last().unwrap(),
                *self.coenergy_at_knots_j_m3.last().unwrap(),
            )
        } else {
            let h = slope.mul_add(width, h0);
            let energy = (0.5 * h0 + 0.5 * h).mul_add(width, self.energy_at_knots_j_m3[segment]);
            // Use slope*width rather than H-H0: a small positive dH can round
            // away in the total H while its work against an existing B matters.
            let coenergy =
                (b0 + 0.5 * width).mul_add(slope * width, self.coenergy_at_knots_j_m3[segment]);
            (h, energy, coenergy)
        };
        let secant = if flux_density_t == 0.0 {
            slope
        } else {
            h / flux_density_t
        };
        if !h.is_finite()
            || !secant.is_finite()
            || secant <= 0.0
            || !energy.is_finite()
            || !coenergy.is_finite()
            || (flux_density_t > 0.0 && (h <= 0.0 || energy <= 0.0 || coenergy <= 0.0))
        {
            return Err(invalid(
                "reversible magnetic response or energy is unrepresentable",
            ));
        }
        Ok(ReversibleMagneticResponse {
            field_strength_a_m: h,
            secant_reluctivity_m_h: secant,
            differential_reluctivity_m_h: slope,
            energy_density_j_m3: energy,
            coenergy_density_j_m3: coenergy,
        })
    }
}
