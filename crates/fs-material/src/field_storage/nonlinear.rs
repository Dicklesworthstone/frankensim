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
    slopes_a_m_t: Vec<f64>,
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
    /// Legendre dual B*H-energy [J/m³], relative to H=0.
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
        let mut slopes_a_m_t = Vec::with_capacity(path.knots.len() - 1);
        for pair in path.knots.windows(2) {
            let width = pair[1].0 - pair[0].0;
            let a = pair[0].1.evidence.value.value;
            let b = pair[1].1.evidence.value.value;
            let slope = (b - a) / width;
            let increment = (0.5 * a + 0.5 * b) * width;
            let energy = energy_at_knots_j_m3.last().copied().unwrap_or(0.0) + increment;
            if !width.is_finite()
                || width <= 0.0
                || !slope.is_finite()
                || slope <= 0.0
                || !increment.is_finite()
                || increment <= 0.0
                || !energy.is_finite()
            {
                return Err(invalid(
                    "reversible magnetic source must have representably positive differential reluctivity and finite energy",
                ));
            }
            slopes_a_m_t.push(slope);
            energy_at_knots_j_m3.push(energy);
        }
        Ok(Self {
            card: card.clone(),
            material,
            path,
            energy_at_knots_j_m3,
            slopes_a_m_t,
        })
    }

    /// Original zero-field material snapshot and source usage receipt.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Canonical point answers at every source breakpoint, retaining provenance.
    pub fn source_path(&self) -> &PropertyLineSlice {
        &self.path
    }
    /// Maximum admitted flux-density magnitude [T]. Never extrapolated.
    pub fn maximum_flux_density_t(&self) -> f64 {
        self.path.knots.last().unwrap().0
    }

    /// Evaluate the pinned canonical property and its exact piecewise-linear
    /// derivative/integral. No independent permeability fit is introduced.
    pub fn evaluate(
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
        let query = self
            .material
            .query()
            .with(&self.path.axis, flux_density_t)
            .map_err(|_| invalid("magnetic field query is unrepresentable"))?;
        let material = self.material.at_query_point(&self.card, &query)?;
        let h = material
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .ok_or(invalid("reversible field source is missing"))?
            .value_si();
        let segment = self
            .path
            .knots
            .partition_point(|(b, _)| *b <= flux_density_t)
            .saturating_sub(1)
            .min(self.slopes_a_m_t.len() - 1);
        let (b0, ref answer) = self.path.knots[segment];
        let h0 = answer.evidence.value.value;
        let slope = self.slopes_a_m_t[segment];
        let secant = if flux_density_t == 0.0 {
            slope
        } else {
            h / flux_density_t
        };
        let energy =
            (0.5 * h0 + 0.5 * h).mul_add(flux_density_t - b0, self.energy_at_knots_j_m3[segment]);
        let coenergy = flux_density_t.mul_add(h, -energy);
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
