//! Material-bound density similarity for the existing fan/network solvers.
//!
//! At fixed geometry and speed, fan static pressure scales with inlet density,
//! while the volume-flow axis does not. The same density factor enters
//! `FanBank`'s pressure, slope, composition and operating-point paths. See the
//! original fan-law source in [`crate::fan_law_source`]. This is an explicit
//! incompressible similarity model, not a new measured curve. Its viscosity
//! applicability is expressed by rotor Reynolds ratio `N/N_ref * nu_ref/nu`.
//! Changes in Mach number, efficiency, stall or chemical compatibility are not
//! inferred from gas identity. Equipment speed and similarity domains remain
//! caller/source requirements.

use core::ops::Deref;

use fs_material::MaterialError;
use fs_material::gas::GasState;

use super::FrozenFlowMedium;
use crate::FanBank;

/// Explicit applicability relative to the original curve's reference medium
/// and reference speed. These are model limits, not physical error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanGasDomain {
    /// Actual inlet density / reference density; inclusive positive bounds.
    pub density_ratio: [f64; 2],
    /// Rotor Reynolds ratio `speed_ratio * nu_reference / nu_inlet` at fixed
    /// geometry. Inclusive positive bounds must be supplied by the caller.
    pub reynolds_ratio: [f64; 2],
    /// Maximum static fan rise / absolute inlet pressure, in `(0, 0.1]`.
    /// Admission checks the whole retained curve including its pressure
    /// tolerance, so every subsequent operating-point trial shares this limit.
    pub maximum_pressure_fraction: f64,
}

/// Source-material, similarity or numerical refusal. No input is modified.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialFanError {
    /// Invalid applicability, non-gas input, repeated binding or arithmetic.
    Invalid(&'static str),
    /// The original material owner refused a stale or incoherent snapshot.
    Material(MaterialError),
    /// A data substitution lies outside the declared equipment similarity.
    OutsideDomain {
        /// Physical quantity that exceeded its support.
        field: &'static str,
        /// Actual dimensionless ratio.
        value: f64,
        /// Inclusive lower support.
        lower: f64,
        /// Inclusive upper support.
        upper: f64,
    },
}
impl core::fmt::Display for MaterialFanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material fan: {self:?}")
    }
}
impl std::error::Error for MaterialFanError {}
impl From<MaterialError> for MaterialFanError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn bounds(range: [f64; 2]) -> bool {
    positive(range[0]) && positive(range[1]) && range[0] <= range[1]
}
fn admit_ratio(field: &'static str, value: f64, range: [f64; 2]) -> Result<f64, MaterialFanError> {
    if !positive(value) || !(range[0]..=range[1]).contains(&value) {
        return Err(MaterialFanError::OutsideDomain {
            field,
            value,
            lower: range[0],
            upper: range[1],
        });
    }
    Ok(value)
}

/// Shared density-affinity admission for frozen material fans and evolving
/// compartment fans. Source admission remains with their original owners.
pub(crate) fn density_affinity_ratio(
    reference: f64,
    inlet: f64,
    domain: [f64; 2],
) -> Result<f64, MaterialFanError> {
    if !positive(reference) || !positive(inlet) || !bounds(domain) {
        return Err(MaterialFanError::Invalid(
            "positive finite fan densities and ordered domain required",
        ));
    }
    admit_ratio("density ratio", inlet / reference, domain)
}

fn gas(medium: FrozenFlowMedium<'_>) -> Result<GasState, MaterialFanError> {
    let state = match medium {
        FrozenFlowMedium::Gas(value) => value.state(),
        FrozenFlowMedium::GasMixture(value) => value.state(),
        FrozenFlowMedium::MoistAir(value) => value.state(),
        FrozenFlowMedium::ThermochemicalGas(value) => value.state(),
        FrozenFlowMedium::DeclaredGas(value) => value,
        FrozenFlowMedium::Liquid(_) => {
            return Err(MaterialFanError::Invalid(
                "gas fan similarity cannot reclassify a liquid as an admitted fan medium",
            ));
        }
    };
    state.validate_coherence()?;
    Ok(*state)
}

/// Immutable physical inputs behind a density-scaled `FanBank`. This numerical
/// binding is carried by graph results and derivatives. Original source owners
/// remain on [`MaterialFanBank`], as for other material-bound flow adapters.
#[derive(Clone, Debug, PartialEq)]
pub struct FanGasBinding {
    reference: GasState,
    inlet: GasState,
    domain: FanGasDomain,
    density_ratio: f64,
    reynolds_ratio: f64,
}

// A composed bank has its factors folded into synthetic sample pressures.
// Keeping that distinction prevents a later material binding from applying
// density a second time after composition discarded the direct binding.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum FanGasScaling {
    Unbound,
    Medium(Box<FanGasBinding>),
    Composed {
        inlet: Box<GasState>,
        maximum_pressure_fraction: f64,
    },
}
impl FanGasBinding {
    /// Actual shared inlet state used by pressure, mass flow and heat transport.
    pub const fn inlet(&self) -> &GasState {
        &self.inlet
    }
    /// Original pressure curve's measurement/reference state, at speed ratio 1.
    pub const fn reference(&self) -> &GasState {
        &self.reference
    }
    /// Actual immutable static-pressure multiplier due to material density.
    pub const fn density_ratio(&self) -> f64 {
        self.density_ratio
    }
    /// Rotor Reynolds similarity at this bank's declared speed.
    pub const fn reynolds_ratio(&self) -> f64 {
        self.reynolds_ratio
    }
    /// Explicit physical applicability; no new source accuracy is asserted.
    pub const fn domain(&self) -> FanGasDomain {
        self.domain
    }

    pub(crate) fn identity_bytes(&self) -> Vec<u8> {
        let mut bytes = b"material-fan-gas-v1\0".to_vec();
        for state in [&self.reference, &self.inlet] {
            for value in [
                state.temperature,
                state.pressure,
                state.density,
                state.sound_speed,
                state.dynamic_viscosity,
                state.thermal_conductivity,
                state.gamma,
                state.specific_gas_constant,
                state.specific_heat_cp,
                state.prandtl,
                state.characteristic_impedance,
                state.water_mole_fraction,
            ] {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
        }
        for value in [
            self.domain.density_ratio[0],
            self.domain.density_ratio[1],
            self.domain.reynolds_ratio[0],
            self.domain.reynolds_ratio[1],
            self.domain.maximum_pressure_fraction,
            self.density_ratio,
            self.reynolds_ratio,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        bytes
    }
}

/// Original reference and inlet material owners beside the ordinary fan bank.
/// Dereferencing supplies the SAME `FanBank` to existing enclosure and graph
/// solvers, fan compositions and sensitivities. No alternate root solver or
/// detached pressure-curve rewrite is involved.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFanBank<'a> {
    bank: FanBank,
    reference_medium: FrozenFlowMedium<'a>,
    inlet_medium: FrozenFlowMedium<'a>,
}
impl<'a> MaterialFanBank<'a> {
    /// Original source medium at which the pressure/volume-flow data were rated.
    pub const fn reference_medium(&self) -> FrozenFlowMedium<'a> {
        self.reference_medium
    }
    /// Same source medium supplied to losses, convection and advected heat.
    pub const fn inlet_medium(&self) -> FrozenFlowMedium<'a> {
        self.inlet_medium
    }
    /// Original curve, arrangement and speed, with the admitted density factor.
    pub const fn bank(&self) -> &FanBank {
        &self.bank
    }
}
impl Deref for MaterialFanBank<'_> {
    type Target = FanBank;
    fn deref(&self) -> &Self::Target {
        &self.bank
    }
}

impl<'a> FrozenFlowMedium<'a> {
    /// Bind this inlet to a fan rated at `reference_medium` and reference speed.
    ///
    /// The original `reference_bank` is unchanged. Its source curve, stall
    /// boundary, speed law and tolerance remain intact; only its shared static
    /// pressure factor gains `rho_inlet/rho_reference`. Build the loss network
    /// and heat path from this SAME medium, then pass the result directly to
    /// `solve_operating_point` or `LossGraph::solve_with_fan`.
    ///
    /// The whole retained static-pressure envelope must remain in the declared
    /// low-pressure range. Reference and inlet must be coherent gas snapshots;
    /// liquid substitutions, already-bound banks, unsupported Reynolds/density
    /// similarity and unrepresentable factors refuse without fallback. This
    /// does not supply fan total pressure, efficiency or new measurement data.
    pub fn bind_fan_bank(
        self,
        reference_bank: &FanBank,
        reference_medium: FrozenFlowMedium<'a>,
        domain: FanGasDomain,
    ) -> Result<MaterialFanBank<'a>, MaterialFanError> {
        if reference_bank.has_material_density() {
            return Err(MaterialFanError::Invalid(
                "bind the original reference bank, not an already material-bound bank",
            ));
        }
        if !bounds(domain.density_ratio)
            || !bounds(domain.reynolds_ratio)
            || !positive(domain.maximum_pressure_fraction)
            || domain.maximum_pressure_fraction > 0.1
        {
            return Err(MaterialFanError::Invalid(
                "invalid gas fan similarity domain",
            ));
        }
        let reference = gas(reference_medium)?;
        let inlet = gas(self)?;
        let density_ratio =
            density_affinity_ratio(reference.density, inlet.density, domain.density_ratio)?;
        let reynolds_ratio = admit_ratio(
            "rotor Reynolds ratio",
            reference_bank.speed_ratio() * (reference.dynamic_viscosity / reference.density)
                / (inlet.dynamic_viscosity / inlet.density),
            domain.reynolds_ratio,
        )?;
        let binding = FanGasBinding {
            reference,
            inlet,
            domain,
            density_ratio,
            reynolds_ratio,
        };
        let mut bank = reference_bank.clone();
        bank.gas_scaling = FanGasScaling::Medium(Box::new(binding));
        let pressure_factor = bank.pressure_factor();
        let flow_factor = bank.flow_factor();
        if !positive(pressure_factor) || !positive(flow_factor) {
            return Err(MaterialFanError::Invalid(
                "fan factors are not representable",
            ));
        }
        let source_peak = bank.curve().points()[0].pressure.value();
        let maximum_rise =
            source_peak * pressure_factor * (1.0 + bank.curve().pressure_tolerance_rel());
        let fraction = maximum_rise / inlet.pressure;
        if !fraction.is_finite()
            || fraction > domain.maximum_pressure_fraction
            || (source_peak > 0.0 && maximum_rise <= 0.0)
            || (maximum_rise > 0.0 && fraction <= 0.0)
        {
            return Err(MaterialFanError::OutsideDomain {
                field: "static pressure rise / inlet pressure",
                value: fraction,
                lower: 0.0,
                upper: domain.maximum_pressure_fraction,
            });
        }
        Ok(MaterialFanBank {
            bank,
            reference_medium,
            inlet_medium: self,
        })
    }
}

impl FanBank {
    // Existing speed and pressure-sample derivatives need a two-sided
    // neighborhood inside the material similarity limits as well as the
    // original equipment curve's numerical domain.
    pub(crate) fn material_controls_are_interior(&self) -> bool {
        self.gas_density_binding().is_none_or(|gas| {
            let domain = gas.domain();
            let pressure_fraction = self.curve().points()[0].pressure.value()
                * self.pressure_factor()
                * (1.0 + self.curve().pressure_tolerance_rel())
                / gas.inlet().pressure;
            gas.reynolds_ratio() > domain.reynolds_ratio[0]
                && gas.reynolds_ratio() < domain.reynolds_ratio[1]
                && pressure_fraction < domain.maximum_pressure_fraction
        })
    }
}
