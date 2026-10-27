//! Finite isotropic thermal reference from a source-declared mean coefficient.
//!
//! The source convention is engineering length: s = L(T)/L(T_ref) =
//! 1 + mean_alpha * (T - T_ref). The existing mean-expansion owner admits the
//! source-declared endpoint interval. An instantaneous small-strain alpha integral is
//! not silently reinterpreted as this finite law.
//!
//! F = Fe Fth with Fth = s I. The native law is per thermally relaxed volume;
//! returning to the original reference gives W0 = s^3 W(Fe), P0 = s^2 Pe,
//! and A0 = s Ae. Density belongs to the original reference at T_ref.
//! Temperature is held per element during mechanics. This preparation supplies
//! no heat equation, entropy evolution, temperature-change work or heat claim.

use super::{FiniteSolidError, FiniteSolidMaterial};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::MaterialError;
use fs_material::hyper::cards::ResolvedHyperelastic;
use fs_material::state_point::{
    MaterialPropertySelection, MeanThermalExpansionStrain,
    is_fixed_external_coordinate, resolve_mean_thermal_expansion_strain,
};

/// A same-card, explicitly sourced isotropic stress-free length ratio.
///
/// Reference and endpoint temperature roles are explicit in the mean source.
/// The remaining physical coordinates and their exact quantity descriptors
/// are copied unchanged to both mechanical contexts.
#[derive(Debug, Clone, PartialEq)]
pub struct MeanThermalReference {
    source: MeanThermalExpansionStrain,
    reference_context: QueryPoint,
    current_context: QueryPoint,
    stretch: f64,
    area_ratio: f64,
    volume_ratio: f64,
}

impl MeanThermalReference {
    /// Resolve a property documented as a MEAN engineering-length coefficient.
    ///
    /// The original source owner requires both typed, positive temperatures
    /// and both named roles in source validity. Non-temperature source axes
    /// must be supported fixed external conditions. There is no material-name
    /// inference, instantaneous-alpha conversion or extrapolation.
    pub fn resolve(
        card: &MaterialCard,
        coefficient: &PropertyKey,
        point: &QueryPoint,
        reference_axis: &str,
        endpoint_axis: &str,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FiniteSolidError> {
        if point.axes().contains_key("T") && reference_axis != "T" && endpoint_axis != "T" {
            return Err(FiniteSolidError::Invalid(
                "thermal source has an ambiguous extra canonical temperature coordinate",
            ));
        }
        let source = resolve_mean_thermal_expansion_strain(
            card, coefficient, point, reference_axis, endpoint_axis, selection,
        ).map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        let selected = source.resolved().properties()[0].answer().receipt.selected;
        let claim = card.claims().claim(selected).ok_or(FiniteSolidError::Invalid(
            "mean expansion lost its immutable selected source",
        ))?;
        // Curves may declare a value-dependence axis without a validity bound.
        // Admit the union of applicability and actual source coordinates.
        for axis in claim.validity.bounds().keys().chain(
            source.resolved().properties()[0].source_axes().iter().map(|(axis, _)| axis),
        ) {
            if axis != reference_axis && axis != endpoint_axis
                && !is_fixed_external_coordinate(
                    axis, point.axis_quantities().get(axis).copied(),
                    point.axes().get(axis).copied(),
                )
            {
                return Err(FiniteSolidError::Source(format!(
                    "mean expansion axis {axis} requires an evolving physical owner",
                )));
            }
        }
        let context = |temperature_axis: &str| -> Result<QueryPoint, FiniteSolidError> {
            let mut context = QueryPoint::new();
            for (axis, value) in point.axes() {
                if axis == reference_axis || axis == endpoint_axis { continue; }
                context = match point.axis_quantities().get(axis) {
                    Some(quantity) => context.with_quantity(axis, *quantity, *value),
                    None => context.with(axis, *value),
                }.map_err(|error| FiniteSolidError::Source(error.to_string()))?;
            }
            let quantity = point.axis_quantities().get(temperature_axis).copied()
                .ok_or(FiniteSolidError::Invalid("mean source temperature is not typed"))?;
            let temperature = point.axes().get(temperature_axis).copied()
                .ok_or(FiniteSolidError::Invalid("mean source temperature is missing"))?;
            context.with_quantity("T", quantity, temperature)
                .map_err(|error| FiniteSolidError::Source(error.to_string()))
        };
        let reference_context = context(reference_axis)?;
        let current_context = context(endpoint_axis)?;
        let stretch = 1.0 + source.engineering_strain();
        let area_ratio = stretch * stretch;
        let volume_ratio = area_ratio * stretch;
        if [stretch, area_ratio, volume_ratio].iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(FiniteSolidError::Invalid(
                "thermal reference length/area/volume ratio is not representable",
            ));
        }
        Ok(Self {
            source, reference_context, current_context, stretch, area_ratio, volume_ratio,
        })
    }

    /// Original mean coefficient, source interval, receipts and identity.
    pub const fn source(&self) -> &MeanThermalExpansionStrain { &self.source }
    /// Exact physical context for the original reference density.
    pub const fn reference_context(&self) -> &QueryPoint { &self.reference_context }
    /// Exact physical context for the held-temperature constitutive law.
    pub const fn current_context(&self) -> &QueryPoint { &self.current_context }
    /// Isotropic stress-free length ratio relative to the original mesh.
    pub const fn stretch(&self) -> f64 { self.stretch }
    /// Area Jacobian used in the Piola pullback.
    pub const fn area_ratio(&self) -> f64 { self.area_ratio }
    /// Relaxed-volume/original-volume Jacobian used in stored energy.
    pub const fn volume_ratio(&self) -> f64 { self.volume_ratio }

    pub(super) fn elastic_deformation(&self, deformation: &[f64; 9]) -> [f64; 9] {
        deformation.map(|value| value / self.stretch)
    }

    pub(super) fn check_response<'a>(
        &self,
        mut values: impl Iterator<Item = &'a f64>,
    ) -> Result<(), MaterialError> {
        if values.any(|value| !value.is_finite()) {
            Err(MaterialError::State {
                what: "nonfinite thermal-reference stress or tangent".into(),
            })
        } else {
            Ok(())
        }
    }
}

impl FiniteSolidMaterial {
    /// Bind density at the original reference temperature and the constitutive
    /// law at the source's endpoint temperature, from the SAME material card.
    ///
    /// Both endpoint model contexts are admitted. Mass is subsequently fixed
    /// on the original mesh. The optional stress-free thermal reference enters
    /// every native energy, Piola and tangent evaluation; no geometry scaling
    /// or thermal force vector substitutes for the finite constitutive response.
    pub fn resolve_mean_thermal(
        card: &MaterialCard,
        model_index: usize,
        thermal: MeanThermalReference,
        density_selection: MaterialPropertySelection,
    ) -> Result<Self, FiniteSolidError> {
        if thermal.source().resolved().card_identity() != card.content_hash()
            || thermal.source().resolved().material() != card.id()
        {
            return Err(FiniteSolidError::Source(
                "thermal reference and mechanics require the same immutable material card".into(),
            ));
        }
        let mut material = Self::resolve(
            card, model_index, thermal.reference_context(), density_selection,
        )?;
        material.law = ResolvedHyperelastic::try_from_model_card_at(
            &card.models()[model_index], thermal.current_context(),
        ).map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        material.thermal = Some(thermal);
        Ok(material)
    }

    /// Present only when this material has a source-declared thermal reference.
    pub const fn thermal_reference(&self) -> Option<&MeanThermalReference> {
        self.thermal.as_ref()
    }
}
