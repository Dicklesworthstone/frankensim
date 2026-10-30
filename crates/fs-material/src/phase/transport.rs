//! Phase-first source snapshots for a fixed-density stationary melting model.
//!
//! The original h chart decides solid/mushy/liquid topology before a liquid
//! property is requested. Every state retains sourced density, conductivity
//! and explicit additional properties (for example electrical resistivity).
//! Exactly solid material has no fluid transport tuple: viscosity and liquid
//! Cp are absent, not extrapolated, zeroed or replaced by arbitrary constants.
//!
//! The existing ThermalMaterialPoint updates T and an optional liquid MASS
//! fraction coordinate from that same chart. The explicit shear entry also
//! consumes the native physical rate; other context remains fixed.
//! This is a constant equal-phase-density reduction, not a volume-change law,
//! equilibrium solver, moving solid skeleton or free-surface model.

use crate::MaterialError;
use crate::fluid_transport::{FluidTransport, FluidTransportState};
use crate::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY,
};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::state_point::{
    DENSITY_PROPERTY, MAX_MATERIAL_STATE_PROPERTIES, MaterialPropertySelection,
    ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use crate::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Frequency, QuantitySpec};

use super::EquilibriumPhaseState;
use super::mushy::{MushyPhaseLaw, PhaseMobility};

/// A retained source object with an exact caloric phase and conductivity even
/// when its fluid transport is undefined. Consumers classify the phase first:
/// never call FluidTransport::transport_state on an immobile solid.
pub trait PhaseTransport: FluidTransport {
    fn phase_state(&self) -> &EquilibriumPhaseState;
    /// Density proposed for transported/inertial mass, available in exact solid
    /// without asking for liquid properties. An explicit bounded projection
    /// reports its retained reference here and keeps raw source rho separately.
    fn transport_density_kg_m3(&self) -> f64;
    fn thermal_conductivity_w_m_k(&self) -> Result<f64, MaterialError>;
}

/// One source binding shared by thermal, electrical and momentum responses.
/// The ordinary entry evolves only temperature/fraction; the explicit shear
/// entry also owns its declared rate. Source pressure and composition must
/// describe the same fixed context as the enthalpy chart.
#[derive(Clone, Debug)]
pub struct PhaseTransportSourceQuery<'a> {
    pub source: ThermalMaterialPoint<'a>,
    pub viscosity_basis: LiquidViscosityBasis,
    /// Additional static scalar responses, retained in the SAME material
    /// snapshot. Canonical density/k/Cp/viscosity requirements are owned here;
    /// callers may not override or duplicate them through this list.
    pub extra_requirements: &'a [ScalarPropertyRequirement],
    /// Applies to exactly the requirements of the classified branch. A pinned
    /// selection therefore names no liquid-only property for an immobile solid.
    pub selection: MaterialPropertySelection,
}

/// Original equilibrium state and all source receipts from one atomic query.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPhaseTransportState {
    mobility: PhaseMobility,
    material: ResolvedMaterialStatePoint,
    thermal_conductivity_w_m_k: f64,
    mobile_transport: Option<FluidTransportState>,
}
impl ResolvedPhaseTransportState {
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    #[must_use]
    pub const fn mobility(&self) -> &PhaseMobility {
        &self.mobility
    }
    #[must_use]
    pub const fn phase_state(&self) -> &EquilibriumPhaseState {
        match &self.mobility {
            PhaseMobility::ImmobileSolid { phase, .. } => phase,
            PhaseMobility::Mobile(state) => &state.phase,
        }
    }
    #[must_use]
    pub const fn mobile_transport(&self) -> Option<&FluidTransportState> {
        self.mobile_transport.as_ref()
    }
    /// Unprojected source density, including when no liquid transport exists.
    /// A fixed-mass approximation must retain the chart density separately.
    #[must_use]
    pub fn source_density_kg_m3(&self) -> f64 {
        self.material
            .property(DENSITY_PROPERTY)
            .expect("phase transport always resolves density")
            .value_si()
    }
}
impl FluidTransport for ResolvedPhaseTransportState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.mobile_transport.ok_or_else(|| invalid(
            "exactly solid material has no liquid transport; preserve its stationary phase branch",
        ))
    }
}
impl PhaseTransport for ResolvedPhaseTransportState {
    fn phase_state(&self) -> &EquilibriumPhaseState {
        ResolvedPhaseTransportState::phase_state(self)
    }
    fn transport_density_kg_m3(&self) -> f64 {
        self.source_density_kg_m3()
    }
    fn thermal_conductivity_w_m_k(&self) -> Result<f64, MaterialError> {
        Ok(self.thermal_conductivity_w_m_k)
    }
}

/// Resolve source properties only after classification on the ORIGINAL h chart.
///
/// A fixed-volume consumer requires equal, constant constituent and chart
/// densities. The common source density must match that exact density at each
/// state. Source k and optional electrical responses remain available for a
/// solid, including an all-solid domain. Liquid Cp is an auxiliary transport
/// coefficient only; it never replaces dT/dh or adds latent heat a second time.
///
/// The supplied ThermalMaterialPoint binds the SAME material card as the
/// chart, writes exact T(h) and (when declared) liquid mass fraction, and checks
/// every selected property as a static law on its owned/fixed coordinates.
/// A missing mobile-only source, unsupported coordinate or extrapolation
/// refuses atomically; no partial state or substitute fluid tuple is returned.
pub fn resolve_phase_transport_state(
    law: &MushyPhaseLaw<'_>,
    specific_enthalpy_j_kg: f64,
    query: PhaseTransportSourceQuery<'_>,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    resolve_phase_transport_state_inner(law, specific_enthalpy_j_kg, query, None, None)
}

/// Resolve bounded source-density departures for an explicit reference-mass model.
///
/// The original chart still owns phase, latent energy and its constant reference
/// density. The returned source and mobile tuple retain ACTUAL rho(T,phi), mu
/// and k; this function neither changes mass nor projects momentum/storage.
/// A consumer must explicitly supply that projection (for example Boussinesq
/// buoyancy). The ordinary fixed-density phase-flow entry rejects this tuple
/// when its density differs from the registered caloric mass density, including
/// exactly solid states; a raw sample cannot silently choose a projection.
///
/// The bound must lie in [0,1) and be justified as small for the intended model.
/// Exactly solid material is bounded too, without requesting liquid Cp or mu.
/// Unequal constituent/chart densities remain unsupported: this is not a
/// shrinkage law or a phase-volume source.
pub fn resolve_phase_transport_state_with_density_departure(
    law: &MushyPhaseLaw<'_>,
    specific_enthalpy_j_kg: f64,
    query: PhaseTransportSourceQuery<'_>,
    maximum_relative_density_departure: f64,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    resolve_phase_transport_state_inner(
        law,
        specific_enthalpy_j_kg,
        query,
        None,
        Some(maximum_relative_density_departure),
    )
}

/// Resolve the SAME phase/material bundle at an owned physical shear rate.
///
/// The native momentum owner must supply `sqrt(2 D:D)` in 1/s, not a lattice
/// relaxation parameter, cell speed or cyclic frequency. Only this explicitly
/// typed existing coordinate and the original h-derived T/fraction may evolve;
/// the supplied selection policy, card, pressure and composition stay intact.
/// Source viscosity and conductivity can depend on that rate, so the same
/// snapshot drives momentum, molecular heat, Darcy response and conduction.
///
/// Exact solids still resolve before liquid properties are requested and must
/// receive zero rate. Rate-dependent density or Cp requires a different
/// thermodynamic storage model and refuses, even if its current value happens
/// to match the fixed caloric chart. This adds no constitutive solver, elastic
/// memory, yield history or independent temperature/phase state.
pub fn resolve_phase_transport_state_with_shear_rate(
    law: &MushyPhaseLaw<'_>,
    specific_enthalpy_j_kg: f64,
    query: PhaseTransportSourceQuery<'_>,
    shear_rate_axis: &str,
    shear_rate_per_s: f64,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    resolve_phase_transport_state_inner(
        law,
        specific_enthalpy_j_kg,
        query,
        Some((shear_rate_axis, shear_rate_per_s)),
        None,
    )
}

fn resolve_phase_transport_state_inner(
    law: &MushyPhaseLaw<'_>,
    specific_enthalpy_j_kg: f64,
    query: PhaseTransportSourceQuery<'_>,
    shear: Option<(&str, f64)>,
    maximum_relative_density_departure: Option<f64>,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    let density = law.fixed_density_kg_m3().ok_or_else(|| {
        invalid("phase transport requires exactly equal constant phase and chart densities")
    })?;
    if maximum_relative_density_departure
        .is_some_and(|limit| !limit.is_finite() || !(0.0..1.0).contains(&limit))
    {
        return Err(invalid(
            "phase source density-departure bound must lie in [0,1)",
        ));
    }
    // This must precede even the liquid-property requirement list.
    let mobility = law.classify_at_specific_enthalpy(specific_enthalpy_j_kg)?;
    let mut point = query.source.at(mobility.phase()).map_err(source_error)?;
    let mobile = mobility.mobile().is_some();
    if let Some((axis, rate)) = shear {
        if axis == query.source.temperature_axis
            || Some(axis) == query.source.liquid_fraction_axis
            || query.source.context.axis_quantities().get(axis)
                != Some(&QuantitySpec::dimensional(Frequency::DIMS))
            || !query.source.context.axes().contains_key(axis)
            || !rate.is_finite()
            || rate < 0.0
            || !mobile && rate != 0.0
        {
            return Err(invalid(
                "phase rheology needs an existing typed physical rate, finite nonnegative in liquid and zero in exact solid",
            ));
        }
        point = point.with(axis, rate).map_err(source_error)?;
    }
    let required_count = if mobile { 4 } else { 2 };
    if query.extra_requirements.len() > MAX_MATERIAL_STATE_PROPERTIES - required_count {
        return Err(invalid(
            "phase transport source exceeds the material property budget",
        ));
    }
    let (viscosity_name, viscosity_dims) = match query.viscosity_basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS),
    };
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(source_error)?,
        ScalarPropertyRequirement::try_with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(source_error)?,
    ];
    if mobile {
        requirements.push(
            ScalarPropertyRequirement::try_new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(source_error)?,
        );
        requirements.push(
            ScalarPropertyRequirement::try_new(
                viscosity_name,
                viscosity_dims,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(source_error)?,
        );
    }
    for requirement in query.extra_requirements {
        if [
            DENSITY_PROPERTY,
            THERMAL_CONDUCTIVITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            DYNAMIC_VISCOSITY_PROPERTY,
            KINEMATIC_VISCOSITY_PROPERTY,
        ]
        .contains(&requirement.name())
        {
            return Err(invalid(
                "extra phase properties cannot replace canonical transport requirements",
            ));
        }
        requirements.push(requirement.clone());
    }
    let material = resolve_material_state_point(
        query.source.material,
        &point,
        &requirements,
        query.selection,
    )
    .map_err(source_error)?;
    for property in material.properties() {
        if let Some((axis, _)) = shear {
            if [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY]
                .contains(&property.requirement().name())
                && property.source_axes().iter().any(|(name, _)| name == axis)
            {
                return Err(invalid(
                    "shear-dependent density or heat capacity needs another thermodynamic storage owner",
                ));
            }
            query
                .source
                .admit_static_with_axes(
                    property,
                    &[(axis, QuantitySpec::dimensional(Frequency::DIMS))],
                )
                .map_err(source_error)?;
        } else {
            query.source.admit_static(property).map_err(source_error)?;
        }
    }
    let value = |name| {
        material
            .property(name)
            .expect("required phase transport property was resolved")
            .value_si()
    };
    let source_density = value(DENSITY_PROPERTY);
    let departure = ((source_density - density) / density).abs();
    let density_admitted = match maximum_relative_density_departure {
        Some(limit) => departure.is_finite() && departure <= limit,
        None => source_density.to_bits() == density.to_bits(),
    };
    if !density_admitted || mobility.phase().bulk_density_kg_m3().to_bits() != density.to_bits() {
        return Err(invalid(
            "phase transport source density differs from the original fixed thermal mass density",
        ));
    }
    let thermal_conductivity_w_m_k = value(THERMAL_CONDUCTIVITY_PROPERTY);
    let mobile_transport = if mobile {
        let viscosity = match query.viscosity_basis {
            LiquidViscosityBasis::Dynamic => value(viscosity_name),
            LiquidViscosityBasis::Kinematic => source_density * value(viscosity_name),
        };
        Some(FluidTransportState::try_new(
            source_density,
            value(SPECIFIC_HEAT_CAPACITY_PROPERTY),
            thermal_conductivity_w_m_k,
            viscosity,
        )?)
    } else {
        None
    };
    Ok(ResolvedPhaseTransportState {
        mobility,
        material,
        thermal_conductivity_w_m_k,
        mobile_transport,
    })
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}
fn source_error(error: impl core::fmt::Display) -> MaterialError {
    MaterialError::State {
        what: format!("phase transport source: {error}"),
    }
}

#[cfg(test)]
mod tests;
