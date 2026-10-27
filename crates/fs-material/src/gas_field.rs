//! A prescribed, continuous primitive gas field at common static pressure.
//!
//! Temperature and canonical mole fractions are linearly interpolated between
//! supplied stations. Density, sound speed, heat capacities and transport are
//! always reevaluated by `ThermochemicalGasModel`; interpolating those derived
//! coefficients independently would violate the common EOS and transport law.
//! This is a frozen field, not an advection, heat-transfer or condensation
//! solver. The caller retains the measurements and uncertainty of the supplied
//! stations. Material and thermochemistry source receipts survive every sample.

/// Shared explicit species-column admission for authored primitive fields.
pub mod input;

/// Conservative stationary heat/species generation using the shared thermal owner.
pub mod steady;

use fs_matdb::QueryPoint;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::{Composition, CompositionBasis};

use crate::MaterialError;
use crate::gas::GasState;
use crate::gas_species::input::ThermochemicalPhase;
use crate::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

/// Bounds construction and a single sampling request; no silent decimation.
pub const MAX_FIELD_STATIONS: usize = 4096;

/// Primitive coordinates in the field model's canonical species order.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimitiveGasState {
    /// Absolute temperature [K].
    pub temperature_k: f64,
    /// Explicit mole fractions in the complete source model's canonical order.
    /// Zero means absent locally; species may not be added or reordered outside
    /// the supplied source set.
    pub composition: Composition,
}

/// One reusable axial field, independent of any particular duct geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialGasField {
    stations_m: Vec<f64>,
    states: Vec<PrimitiveGasState>,
    model: ThermochemicalGasModel,
    source_conditions: QueryPoint,
    phase: ThermochemicalPhase,
    uniform: Option<ResolvedThermochemicalGasState>,
}

fn invalid(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn resolve_primitive(
    source: &ThermochemicalGasModel, source_conditions: &QueryPoint,
    phase: ThermochemicalPhase, state: &PrimitiveGasState,
) -> Result<ResolvedThermochemicalGasState, MaterialError> {
    if state.composition.basis() != CompositionBasis::MoleFraction {
        return Err(invalid("gas field stations require explicit mole fractions"));
    }
    let model = source.with_composition(state.composition.clone())?;
    let point = source_conditions.clone()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), state.temperature_k)
        .map_err(|error| invalid(format!("gas field temperature: {error}")))?;
    match phase {
        ThermochemicalPhase::FrozenGas => model.evaluate(&point),
        ThermochemicalPhase::LiquidWaterLimit => {
            if source.water_species().is_none() {
                return Err(invalid("liquid-water field admission needs an explicit source water species"));
            }
            if model.water_species().is_some() {
                model.evaluate_with_liquid_water_limit(&point)
            } else {
                // Explicit zero water preserves the source identity and is unsaturated.
                model.evaluate(&point)
            }
        }
    }
}

impl AxialGasField {
    /// Admit 2..=4096 strictly increasing stations starting at x=0 metres.
    /// The final station is the field length; extrapolation is refused.
    /// `source_conditions` supplies additional common material coordinates.
    /// Temperature and pressure must be absent there because this constructor
    /// binds them explicitly. The phase policy is shared with gas-file input.
    ///
    /// Every supplied station is evaluated before the field is returned. A
    /// nonconvex source-validity domain or intermediate supersaturation may
    /// still refuse at sampling; no untested phase interpolation is asserted.
    pub fn new(
        stations_m: Vec<f64>,
        static_pressure_pa: f64,
        states: Vec<PrimitiveGasState>,
        model: ThermochemicalGasModel,
        source_conditions: QueryPoint,
        phase: ThermochemicalPhase,
    ) -> Result<Self, MaterialError> {
        if !(2..=MAX_FIELD_STATIONS).contains(&stations_m.len())
            || stations_m.len() != states.len()
            || stations_m[0] != 0.0
            || stations_m.iter().any(|x| !x.is_finite())
            || stations_m.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(invalid(
                "gas field needs 2..=4096 increasing matched stations starting at zero",
            ));
        }
        if source_conditions.axes().contains_key("temperature")
            || source_conditions.axes().contains_key("pressure")
        {
            return Err(invalid(
                "gas field owns temperature and common pressure; duplicate source coordinates refuse",
            ));
        }
        let source_conditions = source_conditions
            .with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                static_pressure_pa,
            )
            .map_err(|error| invalid(format!("gas field pressure: {error}")))?;
        let mut field = Self {
            stations_m,
            states,
            model,
            source_conditions,
            phase,
            uniform: None,
        };
        let first = field.evaluate_primitive(&field.states[0])?;
        let mut uniform = true;
        for pair in field.states.windows(2) {
            if pair[0] != pair[1] {
                field.evaluate_primitive(&pair[1])?;
                uniform = false;
            }
        }
        if uniform {
            field.uniform = Some(first);
        }
        Ok(field)
    }

    /// Original supplied station positions [m], retained without resampling.
    #[must_use]
    pub fn stations_m(&self) -> &[f64] {
        &self.stations_m
    }

    /// Original supplied temperatures and canonical mole compositions.
    #[must_use]
    pub fn primitive_states(&self) -> &[PrimitiveGasState] {
        &self.states
    }

    /// Shared source model; original material cards and pack identities remain
    /// here. Its initial composition does not override the station inputs.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalGasModel {
        &self.model
    }

    /// Original frozen-gas or liquid-water admission policy. Runtime consumers
    /// must retain this policy when rebuilding a field from a new gas state.
    #[must_use]
    pub const fn phase(&self) -> ThermochemicalPhase {
        self.phase
    }

    /// Length of the supplied field [m].
    #[must_use]
    pub fn length_m(&self) -> f64 {
        *self.stations_m.last().expect("validated stations")
    }

    /// Exact common static pressure [Pa].
    #[must_use]
    pub fn static_pressure_pa(&self) -> f64 {
        self.source_conditions.axes()["pressure"]
    }

    /// A uniform field can preserve the existing uniform consumer's operation
    /// order instead of unnecessarily introducing discretization interfaces.
    #[must_use]
    pub fn uniform_state(&self) -> Option<&ResolvedThermochemicalGasState> {
        self.uniform.as_ref()
    }

    fn evaluate_primitive(
        &self, state: &PrimitiveGasState,
    ) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        resolve_primitive(&self.model, &self.source_conditions, self.phase, state)
    }

    /// Reevaluate the physical medium at one position, preserving exact supplied
    /// knots and uniform spans. Only primitive T and mole fractions interpolate.
    /// A bounded roundoff residual is added to the largest interpolated fraction
    /// so the existing thermochemistry owner receives a unit simplex. Authored
    /// station fractions are never normalized or repaired.
    pub fn state_at(&self, x_m: f64) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        if !x_m.is_finite() || x_m < 0.0 || x_m > self.length_m() {
            return Err(invalid(
                "gas field sample lies outside its supplied axial interval",
            ));
        }
        if let Some(state) = &self.uniform {
            return Ok(state.clone());
        }
        let right = self.stations_m.partition_point(|x| *x < x_m);
        if self.stations_m[right] == x_m {
            return self.evaluate_primitive(&self.states[right]);
        }
        let left = right - 1;
        let a = &self.states[left];
        let b = &self.states[right];
        if a == b {
            return self.evaluate_primitive(a);
        }
        let w = (x_m - self.stations_m[left]) / (self.stations_m[right] - self.stations_m[left]);
        let temperature_k = a.temperature_k + w * (b.temperature_k - a.temperature_k);
        let composition = if a.composition == b.composition {
            a.composition.clone()
        } else {
            let mut fractions = a
                .composition
                .fractions()
                .iter()
                .zip(b.composition.fractions())
                .map(|(&a, &b)| a + w * (b - a))
                .collect::<Vec<_>>();
            let largest = fractions
                .iter()
                .enumerate()
                .max_by(|(ia, a), (ib, b)| a.total_cmp(b).then_with(|| ib.cmp(ia)))
                .map(|(i, _)| i)
                .expect("validated active composition");
            for _ in 0..4 {
                let sum: f64 = fractions.iter().sum();
                if sum == 1.0 {
                    break;
                }
                let residual = 1.0 - sum;
                if !sum.is_finite() || residual.abs() > 32.0 * f64::EPSILON * fractions.len() as f64
                {
                    return Err(invalid(
                        "gas field interpolation left the composition simplex",
                    ));
                }
                fractions[largest] += residual;
            }
            Composition::new(CompositionBasis::MoleFraction, fractions)
                .map_err(|error| invalid(format!("gas field interpolated composition: {error}")))?
        };
        self.evaluate_primitive(&PrimitiveGasState {
            temperature_k,
            composition,
        })
    }

    /// Sample the shared acoustic coefficient view at caller-supplied positions.
    /// The field retains the source models; `state_at` exposes full per-sample
    /// receipts when required. A refused sample publishes no partial vector.
    pub fn sample_states(&self, positions_m: &[f64]) -> Result<Vec<GasState>, MaterialError> {
        if positions_m.len() > MAX_FIELD_STATIONS {
            return Err(invalid("gas field sample count exceeds MAX_FIELD_STATIONS"));
        }
        positions_m
            .iter()
            .map(|&x| self.state_at(x).map(|state| *state.state()))
            .collect()
    }
}

#[cfg(test)]
mod tests;
