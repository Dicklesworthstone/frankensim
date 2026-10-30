//! Temperature and buoyancy expansion of the SAME evaluated fluid state.

use super::*;
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Immutable thermal reservoir view at one physical state point.
///
/// Original material inputs retain their source receipts. `try_new` is also an
/// explicit authored-state path; numerical values alone are not material data
/// evidence. Unknown expansion is distinct from a fluid with exactly zero
/// expansion, so consumers cannot silently suppress natural convection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFluidTransportState {
    transport: FluidTransportState,
    temperature_k: f64,
    volumetric_expansion_per_k: Option<f64>,
}

impl ThermalFluidTransportState {
    /// Declare a coherent transport tuple, absolute temperature and optional
    /// isobaric volumetric expansion coefficient [1/K]. Negative expansion is
    /// supported; no liquid value is inferred from temperature or its name.
    pub fn try_new(
        transport: FluidTransportState,
        temperature_k: f64,
        volumetric_expansion_per_k: Option<f64>,
    ) -> Result<Self, MaterialError> {
        if !positive(temperature_k)
            || volumetric_expansion_per_k.is_some_and(|beta| !beta.is_finite())
        {
            return Err(invalid(
                "thermal fluid needs positive finite temperature and finite known expansion",
            ));
        }
        Ok(Self {
            transport,
            temperature_k,
            volumetric_expansion_per_k,
        })
    }

    /// Coefficients at this same temperature/composition/pressure state.
    #[must_use]
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// Absolute reservoir temperature [K].
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Isobaric volumetric expansion [1/K], or unknown (never presumed zero).
    #[must_use]
    pub const fn volumetric_expansion_per_k(&self) -> Option<f64> {
        self.volumetric_expansion_per_k
    }
}

/// Thermal extension of an existing fluid source, without a second registry.
/// Temperature, transport and known expansion must describe the SAME source
/// point. Consumers freezing multiple views must verify their tuples agree.
pub trait ThermalFluidTransport: FluidTransport {
    /// Resolve one coherent snapshot. Missing source temperature refuses.
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError>;
}

impl FluidTransport for ThermalFluidTransportState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}
impl ThermalFluidTransport for ThermalFluidTransportState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        Ok(*self)
    }
}
impl ThermalFluidTransport for GasState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        // Admission checks the original EOS/caloric/acoustic state too. These
        // original gas owners use a frozen-composition ideal-gas EOS, hence
        // beta = -(1/rho) (d rho/dT)_p = 1/T, not a generic-fluid assumption.
        let transport = self.coherent_transport()?;
        let beta = 1.0 / self.temperature;
        if !positive(beta) {
            return Err(invalid("ideal-gas expansion is unrepresentable"));
        }
        ThermalFluidTransportState::try_new(transport, self.temperature, Some(beta))
    }
}
impl ThermalFluidTransport for ResolvedGasState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        self.state().thermal_transport_state()
    }
}
impl ThermalFluidTransport for ResolvedGasMixtureState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        self.state().thermal_transport_state()
    }
}
impl ThermalFluidTransport for ResolvedMoistAirState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        self.state().thermal_transport_state()
    }
}
impl ThermalFluidTransport for ResolvedThermochemicalGasState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        let gas = self.state().thermal_transport_state()?;
        if gas.transport() != self.transport_state()? {
            return Err(invalid(
                "thermochemical gas thermal and momentum transport disagree",
            ));
        }
        Ok(gas)
    }
}
impl ThermalFluidTransport for ResolvedLiquidState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        let query = self.material().query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let mut temperatures = query
            .axis_quantities()
            .iter()
            .filter(|(_, meaning)| **meaning == absolute)
            .filter_map(|(axis, _)| query.axes().get(axis).copied());
        let temperature = temperatures.next().ok_or_else(|| {
            invalid(
                "liquid thermal reservoir needs its source's explicitly typed absolute temperature",
            )
        })?;
        if temperatures.next().is_some() {
            return Err(invalid(
                "liquid source has ambiguous absolute-temperature coordinates",
            ));
        }
        // Only the optional, explicitly sourced isobaric expansion of this
        // SAME bundle may supply beta. Four transport claims alone leave it unknown.
        ThermalFluidTransportState::try_new(
            self.transport_state()?,
            temperature,
            self.volumetric_expansion_per_k(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gas::GasSpec;
    use crate::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, resolve_liquid_state};
    use crate::state_point::{
        DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
        SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_qty::Density;

    fn liquid(point: &QueryPoint) -> ResolvedLiquidState {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            (DENSITY_PROPERTY, Density::DIMS, 1000.0),
            (
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
                4000.0,
            ),
            (
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
                0.6,
            ),
            (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.001),
        ] {
            let key = if name == THERMAL_CONDUCTIVITY_PROPERTY {
                PropertyKey::with_quantity(
                    name,
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::ThermalConductivity,
                        ValueForm::Static,
                    )),
                )
            } else {
                PropertyKey::new(name, dims)
            };
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic thermal-fluid fixture".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "liquid".into(),
                process: "fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        resolve_liquid_state(&card, point, MaterialPropertySelection::SingleClaimOnly).unwrap()
    }

    #[test]
    fn gas_thermal_view_is_the_same_physical_state() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, 200_000.0).unwrap();
        let state = gas.thermal_transport_state().unwrap();
        assert_eq!(state.transport(), gas.transport_state().unwrap());
        assert_eq!(state.temperature_k(), gas.temperature);
        assert_eq!(
            state.volumetric_expansion_per_k(),
            Some(1.0 / gas.temperature)
        );
        let mut stale = gas;
        stale.temperature *= 2.0;
        assert!(stale.thermal_transport_state().is_err());
    }

    #[test]
    fn authored_unknown_and_zero_expansion_remain_distinct() {
        let transport = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
        let unknown = ThermalFluidTransportState::try_new(transport, 300.0, None).unwrap();
        let zero = ThermalFluidTransportState::try_new(transport, 300.0, Some(0.0)).unwrap();
        assert_ne!(unknown, zero);
        assert_eq!(unknown.transport_state().unwrap(), transport);
        assert!(ThermalFluidTransportState::try_new(transport, 0.0, None).is_err());
        assert!(ThermalFluidTransportState::try_new(transport, 300.0, Some(f64::NAN)).is_err());
    }

    #[test]
    fn liquid_temperature_is_explicitly_typed_unique_and_retains_unknown_expansion() {
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = QueryPoint::new()
            .with_quantity("bulk_temperature", absolute.clone(), 321.0)
            .unwrap();
        let source = liquid(&point);
        let thermal = source.thermal_transport_state().unwrap();
        assert_eq!(thermal.temperature_k(), 321.0);
        assert_eq!(thermal.transport(), source.transport_state().unwrap());
        assert_eq!(thermal.volumetric_expansion_per_k(), None);
        assert!(
            liquid(&QueryPoint::new())
                .thermal_transport_state()
                .is_err()
        );
        assert!(
            liquid(&QueryPoint::new().with("T", 321.0).unwrap())
                .thermal_transport_state()
                .is_err()
        );
        let ambiguous = point
            .with_quantity("other_temperature", absolute, 321.0)
            .unwrap();
        assert!(liquid(&ambiguous).thermal_transport_state().is_err());
    }
}
