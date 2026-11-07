//! Original material emissivity on the existing finite-body gray surface seam.

use super::{
    ConductionError, DerivativeMode, EquilibriumEnthalpyPhaseCurve, MaterialEmissivitySource,
    MaterialPropertySelection, MaterialRadiationNodeState, SurfaceEmissivity, ThermalMaterialPoint,
    invalid,
};
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// A pinned original emissivity source for native gray boundaries and links.
///
/// The existing radiation source owner evaluates static scalars, thermal/phase
/// curves and joint temperature/liquid-fraction surfaces. Its analytic chain
/// rule retains emissivity evolution on latent plateaus. The SAME caloric
/// object must be used by the native cell's conductivity law and thermal store.
/// No material-name inference, frozen representative sample, new interpolation
/// or source switching is introduced.
///
/// This explicitly declares the native gray, diffuse, bulk-state closure.
/// A total hemispherical measurement alone does not establish gray/directional
/// independence. The source owner's original (0,1] emissivity restriction
/// remains; convection is zero and belongs to the separately declared gap or
/// boundary. The existing native face owner alone supplies heat and entropy.
#[derive(Debug)]
pub struct MaterialGraySurfaceLaw<'a> {
    source: MaterialEmissivitySource<'a>,
    bounds_j_kg: [f64; 2],
}

impl<'a> MaterialGraySurfaceLaw<'a> {
    /// Bind the original source and finite energy interval without extrapolation.
    /// Both endpoint samples must be supported. Bounds are outer intervals;
    /// every later query still checks its complete material/source conditions.
    pub fn new(
        name: &str,
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        material: ThermalMaterialPoint<'a>,
        selection: MaterialPropertySelection,
        bounds_j_kg: [f64; 2],
    ) -> Result<Self, ConductionError> {
        let [lo, hi] = bounds_j_kg;
        if !lo.is_finite() || !hi.is_finite() || lo >= hi {
            return Err(invalid(
                "native gray source needs finite increasing enthalpy bounds",
            ));
        }
        for receipt in curve.source_receipts() {
            if receipt.query_point.iter().any(|(axis, value)| {
                axis != material.temperature_axis
                    && Some(axis.as_str()) != material.liquid_fraction_axis
                    && (material.context.axes().get(axis) != Some(value)
                        || material.context.axis_quantities().get(axis)
                            != receipt.axis_quantities.get(axis))
            }) {
                return Err(invalid(
                    "native emissivity and caloric sources need the same fixed conditions",
                ));
            }
        }
        let source = MaterialEmissivitySource::prepare(name, curve, material, selection)?;
        let law = Self {
            source,
            bounds_j_kg,
        };
        law.evaluate(lo)?;
        law.evaluate(hi)?;
        Ok(law)
    }

    /// Initial original material/claim snapshot; never a frozen runtime value.
    #[must_use]
    pub fn reference_emissivity(&self) -> &SurfaceEmissivity {
        self.source.snapshot()
    }

    /// Resolve emissivity, its full d(epsilon)/dh and the original source receipt.
    /// The returned single sample has index zero; it invents no mesh vertex.
    /// Native callbacks use this same query and publish no physical state.
    pub fn evaluate(
        &self,
        h: f64,
    ) -> Result<(ImplicitGraySurfaceState, MaterialRadiationNodeState), ConductionError> {
        if !h.is_finite() || !(self.bounds_j_kg[0]..=self.bounds_j_kg[1]).contains(&h) {
            return Err(invalid(
                "enthalpy outside original native gray source bounds",
            ));
        }
        let response = self.source.resolve(0, h, None, DerivativeMode::Primal)?;
        let state = ImplicitGraySurfaceState {
            emissivity: response.state.emissivity.value(),
            emissivity_enthalpy_derivative_kg_j: response.depsilon_dh,
            convection_w_m2_k: 0.0,
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
        };
        Ok((state, response.state))
    }
}

impl ImplicitGraySurfaceLaw for MaterialGraySurfaceLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.source.curve()
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.bounds_j_kg
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.evaluate(h)
            .map(|(state, _)| state)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
