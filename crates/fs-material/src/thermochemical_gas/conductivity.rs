//! Original fixed-composition Wilke/WMS conductivity over a supported T interval.
//!
//! Preparation retains the selected pure viscosity/conductivity lines and their
//! receipts. Queries interpolate ONLY these original affine source laws, then
//! reevaluate the existing Wilke interactions and WMS sum. Mixture conductivity
//! is never refitted, frozen, or sampled at an average temperature as an integral.
use super::{
    MaterialError, ThermochemicalGasModel, coordinates, failure, positive, wms_conductivity,
};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point,
};
use fs_matdb::{
    ClaimSelection, PropertyLineSlice, PropertySupportError, QueryPoint, SelectionPolicy,
};
use fs_qty::{DynViscosity, MolarMass, Temperature};
use fs_thermochem::transport::{WilkeGasViscosityComponentV1, WilkeGasViscosityMixtureV1};
use fs_thermochem::{Composition, SpeciesId};

mod quadrature;
pub use quadrature::{GasConductivityMean, GasConductivityQuadrature};

#[derive(Debug)]
pub enum GasConductivityError {
    Invalid(&'static str),
    Material(MaterialError),
    Source(PropertySupportError),
    Cancelled,
    WorkLimit,
    Quadrature(&'static str),
}
impl core::fmt::Display for GasConductivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) | Self::Quadrature(reason) => {
                write!(f, "gas conductivity: {reason}")
            }
            Self::Material(error) => write!(f, "gas conductivity material: {error}"),
            Self::Source(error) => write!(f, "gas conductivity source: {error}"),
            Self::Cancelled => write!(f, "gas conductivity cancelled"),
            Self::WorkLimit => write!(f, "gas conductivity work allowance exhausted"),
        }
    }
}
impl std::error::Error for GasConductivityError {}
impl From<MaterialError> for GasConductivityError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}

/// Original active species in canonical order. Only immutable views escape
/// the prepared law; full-bank dormant species remain owned by the caller.
#[derive(Clone, Debug)]
pub struct GasConductivityComponentSources {
    pub species: SpeciesId,
    pub molar_mass: MolarMass,
    pub material_identity: fs_blake3::ContentHash,
    pub viscosity: PropertyLineSlice,
    pub conductivity: PropertyLineSlice,
}

/// Owned original transport paths: no borrowed model and no model/card clone
/// per query. Pressure and every auxiliary source coordinate remain fixed.
/// This is a transport law; the caller still owns NASA enthalpy and gas phase.
#[derive(Clone, Debug)]
pub struct PreparedThermochemicalConductivity {
    components: Vec<GasConductivityComponentSources>,
    composition: Composition,
    context: QueryPoint,
    temperatures: Vec<f64>,
    source_packs: Vec<super::packs::ThermochemicalPackReceipt>,
}

fn reserve<T>(count: usize) -> Result<Vec<T>, GasConductivityError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| GasConductivityError::Invalid("source preparation allocation refused"))?;
    Ok(values)
}

impl PreparedThermochemicalConductivity {
    /// Prepare each active original source once on a caller-declared branch.
    ///
    /// Each source line has the matdb 4096-knot cap and proves continuous
    /// support with unchanged selection, including sampled-surface interiors.
    /// Source queries are synchronous: checkpoints bracket each operation and
    /// occur every 256 retained knots. Use a cancellation-only checkpoint here;
    /// preparation checkpoints do not represent WMS point evaluations.
    pub fn new(
        model: &ThermochemicalGasModel,
        context: &QueryPoint,
        temperature_range_k: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), GasConductivityError>,
    ) -> Result<Self, GasConductivityError> {
        checkpoint()?;
        let [lower_k, upper_k] = temperature_range_k;
        if !positive(lower_k) || !positive(upper_k) || lower_k >= upper_k {
            return Err(GasConductivityError::Invalid(
                "a positive finite temperature interval is required",
            ));
        }
        coordinates(context)?;
        let lower = context
            .clone()
            .with("temperature", lower_k)
            .map_err(|error| failure("conductivity lower coordinate", error))?;
        let upper = context
            .clone()
            .with("temperature", upper_k)
            .map_err(|error| failure("conductivity upper coordinate", error))?;
        let requirements = [
            ScalarPropertyRequirement::try_new(
                "dynamic_viscosity",
                DynViscosity::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(|error| failure("conductivity viscosity requirement", error))?,
            ScalarPropertyRequirement::try_new(
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(|error| failure("conductivity requirement", error))?,
        ];
        let mut components = reserve(model.transport_materials().len())?;
        let mut temperatures = reserve(2)?;
        temperatures.extend([lower_k, upper_k]);
        for (data, species) in model
            .transport_materials()
            .iter()
            .zip(model.thermodynamics().components())
        {
            checkpoint()?;
            let bundle = resolve_material_state_point(
                &data.material,
                &lower,
                &requirements,
                data.selection.clone(),
            )
            .map_err(|error| failure("conductivity source admission", error))?;
            let mut lines = reserve(2)?;
            for property in bundle.properties() {
                if !property.has_only_fixed_external_conditions(Some("temperature"), false) {
                    return Err(GasConductivityError::Invalid(
                        "conductivity paths require static temperature and fixed external conditions",
                    ));
                }
                let selected = property.answer().receipt.selected;
                let claim =
                    data.material
                        .claims()
                        .claim(selected)
                        .ok_or(GasConductivityError::Invalid(
                            "selected transport claim is absent",
                        ))?;
                let selection = match &data.selection {
                    MaterialPropertySelection::SingleClaimOnly => {
                        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)
                    }
                    MaterialPropertySelection::PreferObservationBacked => {
                        ClaimSelection::Policy(SelectionPolicy::PreferObservationBacked)
                    }
                    MaterialPropertySelection::PinnedByProperty(_) => {
                        ClaimSelection::Pinned(selected)
                    }
                };
                checkpoint()?;
                let line = data
                    .material
                    .claims()
                    .query_line_typed(&claim.key, "temperature", &lower, &upper, selection)
                    .map_err(GasConductivityError::Source)?;
                checkpoint()?;
                temperatures
                    .try_reserve(line.knots.len())
                    .map_err(|_| GasConductivityError::Invalid("source knot allocation refused"))?;
                for (index, &(t, ref answer)) in line.knots.iter().enumerate() {
                    if index.is_multiple_of(256) {
                        checkpoint()?;
                    }
                    if !positive(answer.evidence.value.value) || answer.receipt.selected != selected
                    {
                        return Err(GasConductivityError::Invalid(
                            "original transport line must stay positive with unchanged source selection",
                        ));
                    }
                    temperatures.push(t);
                }
                lines.push(line);
            }
            if lines.len() != 2 {
                return Err(GasConductivityError::Invalid(
                    "both original transport paths are required",
                ));
            }
            let conductivity = lines.pop().expect("two admitted original paths");
            let viscosity = lines.pop().expect("two admitted original paths");
            components.push(GasConductivityComponentSources {
                species: species.species().clone(),
                molar_mass: species.molar_mass(),
                material_identity: data.material.content_hash(),
                viscosity,
                conductivity,
            });
        }
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        checkpoint()?;
        Ok(Self {
            components,
            composition: model.thermodynamics().mole_composition().clone(),
            context: context.clone(),
            temperatures,
            source_packs: model.source_packs().to_vec(),
        })
    }

    pub fn temperature_range_k(&self) -> [f64; 2] {
        [
            self.temperatures[0],
            self.temperatures[self.temperatures.len() - 1],
        ]
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures
    }
    pub fn components(&self) -> &[GasConductivityComponentSources] {
        &self.components
    }
    pub fn mole_composition(&self) -> &Composition {
        &self.composition
    }
    pub fn context(&self) -> &QueryPoint {
        &self.context
    }
    pub fn source_packs(&self) -> &[super::packs::ThermochemicalPackReceipt] {
        &self.source_packs
    }

    fn admit(&self, temperature: f64) -> Result<(), GasConductivityError> {
        let [lo, hi] = self.temperature_range_k();
        if !temperature.is_finite() || !(lo..=hi).contains(&temperature) {
            Err(GasConductivityError::Invalid(
                "temperature leaves the original source interval",
            ))
        } else {
            Ok(())
        }
    }

    /// One original WMS point. The checkpoint runs exactly once immediately
    /// before source interpolation and the bounded original Wilke operation.
    /// No NASA inverse, material query, receipt selection or source fit occurs.
    pub fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), GasConductivityError>,
    ) -> Result<f64, GasConductivityError> {
        self.admit(temperature_k)?;
        checkpoint()?;
        let mut viscosities = reserve(self.components.len())?;
        for component in &self.components {
            viscosities.push(
                WilkeGasViscosityComponentV1::new(
                    component.species.clone(),
                    component.molar_mass,
                    DynViscosity::new(source_value(&component.viscosity, temperature_k)),
                )
                .map_err(|error| failure("conductivity species viscosity", error))?,
            );
        }
        let viscosity = WilkeGasViscosityMixtureV1::new(
            Temperature::new(temperature_k),
            viscosities,
            self.composition.clone(),
        )
        .map_err(|error| failure("conductivity Wilke mixture", error))?
        .evaluate()
        .map_err(|error| failure("conductivity Wilke evaluation", error))?;
        wms_conductivity(
            self.components
                .iter()
                .map(|component| Ok(source_value(&component.conductivity, temperature_k))),
            self.composition.fractions(),
            viscosity.receipt().denominator_bits(),
        )
        .map_err(Into::into)
    }

    /// Integrate the actual original WMS function across every crossed μ/k
    /// breakpoint. Reversing endpoints is bit-identical; equality is one point.
    /// Work/depth exhaustion refuses instead of returning an unconverged mean.
    /// The embedded error is an estimate in W/(m K), not an energy defect.
    pub fn mean_conductivity(
        &self,
        from_k: f64,
        to_k: f64,
        controls: GasConductivityQuadrature,
        mut checkpoint: impl FnMut() -> Result<(), GasConductivityError>,
    ) -> Result<GasConductivityMean, GasConductivityError> {
        self.admit(from_k)?;
        self.admit(to_k)?;
        quadrature::mean(&self.temperatures, from_k, to_k, controls, |temperature| {
            self.conductivity_at_temperature(temperature, &mut checkpoint)
        })
        .map_err(|error| match error {
            quadrature::Error::Invalid => {
                GasConductivityError::Invalid("invalid conductivity quadrature controls")
            }
            quadrature::Error::Arithmetic => {
                GasConductivityError::Quadrature("conductivity quadrature is unrepresentable")
            }
            quadrature::Error::WorkLimit => GasConductivityError::WorkLimit,
            quadrature::Error::Resolution => GasConductivityError::Quadrature(
                "conductivity quadrature did not resolve its requested tolerance",
            ),
            quadrature::Error::Integrand(error) => error,
        })
    }
}

// Exactly the admitted source's affine interpolation, with original knots
// returned verbatim. Nonnegative endpoint interpolation remains positive.
fn source_value(line: &PropertyLineSlice, temperature: f64) -> f64 {
    let upper = line.knots.partition_point(|(t, _)| *t < temperature);
    let (right, value) = &line.knots[upper];
    if *right == temperature {
        return value.evidence.value.value;
    }
    let (left, lower) = &line.knots[upper - 1];
    let low = lower.evidence.value.value;
    low + (temperature - left) / (right - left) * (value.evidence.value.value - low)
}

#[cfg(test)]
mod tests;
