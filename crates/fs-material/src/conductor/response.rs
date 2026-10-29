//! Original scalar conductor sources evaluated at owned physical field magnitudes.
//!
//! This pins the complete material bundle; it does not supply a second current
//! law. A spatial owner must solve J=sigma(|E|,|J|,context) E and independently
//! check that these queried coordinates agree with its actual recovered fields.

use super::{ConductorError, ResolvedConductorState};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};
use std::sync::Arc;

/// Rotationally invariant coordinate of an isotropic memoryless source.
/// Signed components, Hall response, anisotropy and hysteresis are not inferred.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConductorMagnitude {
    ElectricField,
    CurrentDensity,
}
impl ConductorMagnitude {
    pub const fn quantity(self) -> QuantitySpec {
        QuantitySpec::dimensional(match self {
            Self::ElectricField => Dims([1, 1, -3, 0, -1, 0]),
            Self::CurrentDensity => Dims([-2, 0, 0, 0, 1, 0]),
        })
    }
}

/// One original scalar law and its explicitly named E/J operating coordinates.
/// The immutable card, independent conductivity/resistivity basis, selected
/// claims, other property requirements and all coordinate quantities survive.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarConductorResponse {
    card: Arc<MaterialCard>,
    original: ResolvedConductorState,
    axes: Vec<(String, ConductorMagnitude)>,
}

fn invalid(quantity: &'static str) -> ConductorError {
    ConductorError::InvalidInput { quantity }
}

impl ScalarConductorResponse {
    /// Bind a liquid MASS fraction supplied by a physical enthalpy owner while
    /// this response owner supplies the actual spatial E/J magnitudes. The
    /// source's existing dimensional phase coordinate must lie in [0,1].
    /// This retains the original source roles, pins and support; it neither
    /// invents a mixing law nor adds a higher-dimensional source table.
    pub fn new_with_liquid_mass_fraction_axis(
        card: &MaterialCard,
        original: ResolvedConductorState,
        axes: &[(&str, ConductorMagnitude)],
        phase_axis: impl Into<String>,
    ) -> Result<Self, ConductorError> {
        if axes.len() > 2 {
            return Err(invalid("at most two conductor operating coordinates"));
        }
        let owned: Vec<_> = axes
            .iter()
            .map(|(axis, magnitude)| (*axis, magnitude.quantity()))
            .collect();
        let original = original.with_phase_axis_and_coordinates(phase_axis.into(), &owned)?;
        Self::new(card, original, axes)
    }

    /// At most one electric-field and one current-density magnitude coordinate.
    /// Empty axes retain an ordinary source in a mixed nonlinear specimen. Each
    /// field coordinate already present in the complete query must be owned;
    /// frequency, strain, stress and unowned response restrictions still refuse.
    pub fn new(
        card: &MaterialCard,
        original: ResolvedConductorState,
        axes: &[(&str, ConductorMagnitude)],
    ) -> Result<Self, ConductorError> {
        if axes.len() > 2 || card.content_hash() != original.material().card_identity() {
            return Err(invalid(
                "original conductor card and at most two operating coordinates",
            ));
        }
        let law = Self {
            card: Arc::new(card.clone()),
            original,
            axes: axes
                .iter()
                .map(|(name, coordinate)| ((*name).to_owned(), *coordinate))
                .collect(),
        };
        for (index, (axis, coordinate)) in law.axes.iter().enumerate() {
            if axis.is_empty()
                || law.axes[..index]
                    .iter()
                    .any(|(other, kind)| other == axis || kind == coordinate)
            {
                return Err(invalid(
                    "unique conductor magnitude names and physical coordinates",
                ));
            }
        }
        law.admit(&law.original)?;
        // Replaying the original pins also refuses a foreign card or a lost
        // selected claim; honest pin-policy receipts need not be identical.
        let replay = law.at_query_point(&law.original.material().query())?;
        if replay
            .material()
            .properties()
            .iter()
            .zip(law.original.material().properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "original conductor values must replay their pinned source",
            ));
        }
        Ok(law)
    }

    pub fn original(&self) -> &ResolvedConductorState {
        &self.original
    }
    pub fn card(&self) -> &MaterialCard {
        &self.card
    }
    pub fn axes(&self) -> &[(String, ConductorMagnitude)] {
        &self.axes
    }

    /// Admit a context refresh with unchanged original source roles and pins.
    /// Temperature/phase owners still own their physical histories and must
    /// account for them. This check never reinterprets a field as temperature.
    pub fn admit(&self, conductor: &ResolvedConductorState) -> Result<(), ConductorError> {
        let original = self.original.material();
        let material = conductor.material();
        if conductor.source_basis() != self.original.source_basis()
            || conductor.source_key() != self.original.source_key()
            || conductor.liquid_mass_fraction_axis() != self.original.liquid_mass_fraction_axis()
            || material.card_identity() != self.card.content_hash()
            || material.axis_quantities() != original.axis_quantities()
            || material.query_point().len() != original.query_point().len()
            || material
                .query_point()
                .iter()
                .zip(original.query_point())
                .any(|((a, _), (b, _))| a != b)
            || material.properties().len() != original.properties().len()
            || material
                .properties()
                .iter()
                .zip(original.properties())
                .any(|(a, b)| {
                    a.requirement() != b.requirement()
                        || a.answer().receipt.selected != b.answer().receipt.selected
                })
        {
            return Err(invalid(
                "conductor response must retain the complete original source bundle",
            ));
        }
        let mut owned: Vec<_> = self
            .axes
            .iter()
            .map(|(name, coordinate)| (name.as_str(), coordinate.quantity()))
            .collect();
        for (axis, coordinate) in &self.axes {
            if material.axis_quantities().get(axis) != Some(&coordinate.quantity())
                || !material
                    .query_point()
                    .iter()
                    .find(|(name, _)| name == axis)
                    .is_some_and(|(_, value)| value.is_finite() && *value >= 0.0)
            {
                return Err(invalid(
                    "existing typed nonnegative conductor magnitude coordinate",
                ));
            }
        }
        if let Some(axis) = conductor.liquid_mass_fraction_axis() {
            owned.push((axis, QuantitySpec::dimensional(Dims([0; 6]))));
        }
        if !conductor
            .source_property()
            .has_only_owned_coordinates_and_fixed_external_conditions(&owned)
        {
            return Err(invalid(
                "conductor law has an unowned response or applicability coordinate",
            ));
        }
        let electrical = |dims| {
            [
                ConductorMagnitude::ElectricField.quantity().dims(),
                ConductorMagnitude::CurrentDensity.quantity().dims(),
            ]
            .contains(&dims)
        };
        if material.axis_quantities().iter().any(|(axis, quantity)| {
            electrical(quantity.dims()) && !self.axes.iter().any(|(name, _)| name == axis)
        }) || material.properties().iter().any(|property| {
            property.source_axes().iter().any(|(axis, dims)| {
                electrical(*dims) && !self.axes.iter().any(|(name, _)| name == axis)
            })
        }) {
            return Err(invalid(
                "every original electrical coordinate must be solved, never frozen",
            ));
        }
        Ok(())
    }

    /// Refresh all original requirements at an explicit context, with no source
    /// fallback. The caller owns every changed non-electrical coordinate.
    pub fn at_query_point(
        &self,
        point: &QueryPoint,
    ) -> Result<ResolvedConductorState, ConductorError> {
        let material = self.original.material().at_query_point(&self.card, point)?;
        let owned: Vec<_> = self
            .axes
            .iter()
            .map(|(axis, magnitude)| (axis.as_str(), magnitude.quantity()))
            .collect();
        let conductor = self
            .original
            .with_material_and_coordinates(material, &owned)?;
        self.admit(&conductor)?;
        Ok(conductor)
    }

    /// Evaluate the SAME complete source at actual spatial E and J. Other
    /// coordinates come from the explicitly admitted seed, e.g. one thermal
    /// Newton endpoint. No scalar fit, constitutive clipping or extrapolation.
    pub fn at_fields(
        &self,
        seed: &ResolvedConductorState,
        electric_field_v_m: [f64; 2],
        current_density_a_m2: [f64; 2],
    ) -> Result<ResolvedConductorState, ConductorError> {
        self.admit(seed)?;
        let e = electric_field_v_m[0].hypot(electric_field_v_m[1]);
        let j = current_density_a_m2[0].hypot(current_density_a_m2[1]);
        if !e.is_finite() || !j.is_finite() {
            return Err(invalid(
                "finite physical electric field and current density",
            ));
        }
        let mut query = seed.material().query();
        for (axis, coordinate) in &self.axes {
            let value = match coordinate {
                ConductorMagnitude::ElectricField => e,
                ConductorMagnitude::CurrentDensity => j,
            };
            query = query
                .with(axis, value)
                .map_err(|_| invalid("representable conductor operating coordinate"))?;
        }
        self.at_query_point(&query)
    }
}
