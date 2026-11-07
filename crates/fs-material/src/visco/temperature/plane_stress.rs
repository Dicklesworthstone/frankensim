//! Calibrated fixed-temperature kinetics on the original reduced source tensor.

use super::{ShiftedViscoelastic, TemperatureShift, invalid};
use crate::plane_stress_viscosity::PlaneStressViscosity;
use crate::visco::ViscoError;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

/// Existing source-bound plane-stress viscosity with its exact temperature
/// declaration. The original resolved material record remains the reference;
/// the prepared matrix and frequency band describe the selected temperature.
pub type ShiftedPlaneStressViscosity = ShiftedViscoelastic<PlaneStressViscosity>;

impl TemperatureShift {
    /// Prepare the original reduced Kelvin viscosity at a fixed temperature.
    /// Every coefficient and cached passive channel uses the existing positive
    /// kinetic shift: `V(T) = a_T V_ref`, with support `[lo/a_T, hi/a_T]`.
    /// Elasticity, geometry, source receipts and all other held environmental
    /// coordinates are unchanged. No live thermal feedback is supplied.
    ///
    /// The exact resolved reference temperature must match this declaration.
    /// Its axis must carry AbsoluteTemperature/Static semantics, or be legacy
    /// `T` with absent or dimensional kelvin metadata. Selected coefficients
    /// may depend on other admitted fixed external conditions, but must not
    /// already depend on this temperature axis; doing both would double-shift
    /// the same source kinetics. A second temperature owner also refuses.
    ///
    /// Call this on the original reference law every time, as for the tensor
    /// shift. The caller qualifies this calibrated superposition against that
    /// source; this method does not infer calibration from material names or
    /// reinterpret the source's original reference-temperature receipt.
    pub fn shift_plane_stress_viscosity(
        &self,
        reference: &PlaneStressViscosity,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<ShiftedPlaneStressViscosity, ViscoError> {
        let source = reference.source();
        let quantity = source.axis_quantities().get(temperature_axis).copied();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let temperature_owner = quantity == Some(absolute)
            || temperature_axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if !temperature_owner
            || source.state_coordinate(temperature_axis) != Some(self.reference_temperature_k())
        {
            return Err(invalid(
                "viscosity shift needs its exact original reference-temperature coordinate",
            ));
        }
        for name in reference.properties() {
            let property = source.property(name).ok_or_else(|| {
                invalid("viscosity shift is missing an original selected coefficient")
            })?;
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), true)
                || property
                    .source_axes()
                    .iter()
                    .any(|(axis, _)| axis == temperature_axis)
            {
                return Err(invalid(
                    "viscosity shift needs held external context without an existing temperature-dependent coefficient",
                ));
            }
        }
        let state = self.evaluate(temperature_k)?;
        let law = reference.with_kinetic_scale(state.factor())?;
        Ok(ShiftedViscoelastic { law, state })
    }
}
