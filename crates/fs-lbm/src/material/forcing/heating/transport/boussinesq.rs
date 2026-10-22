//! Material density drives buoyancy without changing conserved reference mass.
//!
//! The Boussinesq momentum equation uses `(rho(T) - rho0) * gravity`, after
//! removing the hydrostatic reference load. Inertia and sensible storage use
//! the retained `rho0` and `cp0`; viscosity and conductivity use the SAME local
//! material evaluation as buoyancy. Explicit bounds admit the neglected
//! density and heat-capacity departures. This is not variable-mass transport.

use super::{
    FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError, MaterialSlab,
    MaterialTransportError, MaterialTransportStep, Slab, ThermalLbm, ThermalLbm3,
    candidate_with_prepared_input, poll, sample,
};
use fs_exec::CancelGate;

/// Explicit Boussinesq approximation and physical gravity for one step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoussinesqConfig<const D: usize> {
    /// Signed gravitational acceleration [m/s²], in the solver's axes.
    /// For ordinary downward gravity the final component is negative.
    pub gravity_m_s2: [f64; D],
    /// Maximum `abs(rho(T) - rho0) / rho0`, checked at both query epochs.
    /// Must be in `[0,1)`; the caller must choose a physically small allowance.
    pub max_relative_density_departure: f64,
    /// Maximum `abs(cp(T) - cp0) / cp0` neglected by fixed sensible storage.
    /// Must be in `[0,1)`. Zero requires exact reference heat capacity.
    pub max_relative_specific_heat_departure: f64,
}

impl<const D: usize> BoussinesqConfig<D> {
    /// Check the declared approximation bounds without evaluating a source.
    pub fn validate(self) -> Result<(), MaterialTransportError> {
        if !self.gravity_m_s2.iter().all(|g| g.is_finite())
            || [
                self.max_relative_density_departure,
                self.max_relative_specific_heat_departure,
            ]
            .iter()
            .any(|v| !v.is_finite() || !(0.0..1.0).contains(v))
        {
            return Err(MaterialTransportError::Invalid(
                "Boussinesq gravity must be finite and both relative departure bounds must lie in [0,1)",
            ));
        }
        Ok(())
    }
}

/// Original source response and its explicit Boussinesq transport projection.
///
/// The source density is used only in buoyancy. The transport view retains
/// reference density and heat capacity so changing a local source evaluation
/// cannot reinterpret mass or stored sensible energy. This projection carries
/// no source authority beyond the original material's own documented model.
#[derive(Debug, Clone, PartialEq)]
pub struct BoussinesqMaterial<F> {
    material: F,
    resolved: FluidTransportState,
    transport: FluidTransportState,
}

impl<F> BoussinesqMaterial<F> {
    /// Original concrete source owner, including its retained query receipts.
    pub const fn material(&self) -> &F {
        &self.material
    }
    /// Unmodified local density, heat capacity, viscosity and conductivity.
    pub const fn resolved_transport(&self) -> &FluidTransportState {
        &self.resolved
    }
    /// Actual numerical transport: reference rho/cp, local viscosity/conductivity.
    pub const fn reference_storage_transport(&self) -> &FluidTransportState {
        &self.transport
    }
}

impl<F: FluidTransport> BoussinesqMaterial<F> {
    /// Resolve one source once, admit its density/Cp departures and retain
    /// reference storage with this SAME sample's molecular viscosity and k.
    /// This is shared by the lattice and finite-volume enthalpy owners.
    pub fn try_new<const D: usize>(
        material: F,
        reference: &FluidTransportState,
        config: BoussinesqConfig<D>,
    ) -> Result<Self, MaterialError> {
        config.validate().map_err(|error| MaterialError::State {
            what: error.to_string(),
        })?;
        resolve_bounded(material, reference, config)
    }
}

impl<F> FluidTransport for BoussinesqMaterial<F> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}

fn resolve_bounded<F: FluidTransport, const D: usize>(
    material: F,
    reference: &FluidTransportState,
    config: BoussinesqConfig<D>,
) -> Result<BoussinesqMaterial<F>, MaterialError> {
    // Query the original owner once. Further sampling reads this frozen view,
    // never an interior-mutable material owner a second time.
    let resolved = material.transport_state()?;
    for (name, value, retained, limit) in [
        (
            "density",
            resolved.density_kg_m3(),
            reference.density_kg_m3(),
            config.max_relative_density_departure,
        ),
        (
            "specific heat",
            resolved.specific_heat_j_kg_k(),
            reference.specific_heat_j_kg_k(),
            config.max_relative_specific_heat_departure,
        ),
    ] {
        let departure = ((value - retained) / retained).abs();
        if !departure.is_finite() || departure > limit {
            return Err(MaterialError::State {
                what: format!(
                    "Boussinesq {name} relative departure {departure} exceeds admitted bound {limit}",
                ),
            });
        }
    }
    let transport = FluidTransportState::try_new(
        reference.density_kg_m3(),
        reference.specific_heat_j_kg_k(),
        resolved.thermal_conductivity_w_m_k(),
        resolved.dynamic_viscosity_pa_s(),
    )?;
    Ok(BoussinesqMaterial {
        material,
        resolved,
        transport,
    })
}

pub(in crate::material::forcing) trait BackgroundBuoyancy {
    fn coefficient(&self) -> f64;
    fn set_coefficient(&mut self, value: f64);
}
impl BackgroundBuoyancy for ThermalLbm {
    fn coefficient(&self) -> f64 {
        self.gbeta
    }
    fn set_coefficient(&mut self, value: f64) {
        self.gbeta = value;
    }
}
impl BackgroundBuoyancy for ThermalLbm3 {
    fn coefficient(&self) -> f64 {
        self.gbeta
    }
    fn set_coefficient(&mut self, value: f64) {
        self.gbeta = value;
    }
}

#[allow(clippy::too_many_arguments)]
fn advance<
    S: Slab<D, Q> + BackgroundBuoyancy,
    F: FluidTransport,
    const D: usize,
    const Q: usize,
>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    config: BoussinesqConfig<D>,
    gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<MaterialTransportStep<BoussinesqMaterial<F>>, MaterialTransportError> {
    poll(gate)?;
    config.validate()?;
    let (cells, plane) = slab.lattice.layout()?;
    let input = slab.heat_inputs(power, forces, cells, depth_m, limits)?;
    if slab.lattice.walls() != (1.0, 0.0) || !slab.lattice.flow_valid() {
        return Err(MaterialTransportError::Invalid(
            "invalid material slab or independently changed normalized walls",
        ));
    }
    if power[..plane]
        .iter()
        .chain(&power[cells - plane..])
        .any(|p| *p != 0.0)
    {
        return Err(MaterialTransportError::Invalid("wall power must be zero"));
    }
    let mut bounded = |cell, temperature, reference: &FluidTransportState| {
        resolve_bounded(resolve(cell, temperature, reference)?, reference, config)
    };
    let applied = sample::<S, BoussinesqMaterial<F>, D, Q>(
        slab,
        &slab.lattice,
        cells,
        plane,
        limits.temperature_range_k,
        gate,
        &mut bounded,
    )?;
    let (trial, heat) = candidate_with_buoyancy::<S, BoussinesqMaterial<F>, D, Q>(
        slab, input, forces, &applied, config, gate,
        |material| material.resolved.density_kg_m3(),
    )?;
    let endpoint = sample::<S, BoussinesqMaterial<F>, D, Q>(
        slab,
        &trial,
        cells,
        plane,
        limits.temperature_range_k,
        gate,
        &mut bounded,
    )?;
    poll(gate)?;
    slab.lattice = trial;
    Ok(MaterialTransportStep {
        heat,
        applied,
        endpoint,
    })
}

// Joule/Lorentz and prescribed-source callers share the SAME material force
// projection and tentative transport step. The supplied density is from the
// already admitted BoussinesqMaterial, never a second constitutive query.
#[allow(clippy::too_many_arguments)]
pub(in crate::material::forcing) fn candidate_with_buoyancy<
    S: Slab<D, Q> + BackgroundBuoyancy, F, const D: usize, const Q: usize,
>(
    slab: &MaterialSlab<S>, mut input: super::super::Inputs<D>, forces: &[[f64; D]],
    applied: &[super::TransportSample<F>], config: BoussinesqConfig<D>,
    gate: &CancelGate, density: impl Fn(&F) -> f64,
) -> Result<(S, super::HeatBalance), MaterialTransportError> {
    let units = slab.force_units;
    for value in applied {
        poll(gate)?;
        let departure = density(&value.material) - units.transport.density_kg_m3();
        for axis in 0..D {
            let force = forces[value.cell][axis] + departure * config.gravity_m_s2[axis];
            input.force[value.cell][axis] = units.component(force)
                .map_err(|source| MaterialTransportError::State {
                    cell: value.cell, temperature_k: value.temperature_k, source,
                })?;
        }
    }
    // Replace, rather than superpose, the legacy beta*(T-Tref) contribution.
    let (mut trial, heat) = candidate_with_prepared_input::<S, F, D, Q>(
        slab, input, applied, gate, |candidate| candidate.set_coefficient(0.0),
    )?;
    trial.set_coefficient(slab.lattice.coefficient());
    Ok((trial, heat))
}

impl MaterialSlab<ThermalLbm> {
    /// Advance 2-D heat/flow with density-resolved material buoyancy.
    ///
    /// Each initial local temperature resolves ONE material evaluation. Its
    /// rho supplies `(rho-rho0)*gravity`; its mu and k supply relaxation with
    /// reference inertia/storage: `nu=mu/rho0`, `alpha=k/(rho0*cp0)`. Original
    /// source responses and this projection are both retained in the report.
    /// External force [N/m³] and prescribed heat [W/m³] retain their existing
    /// meanings. The slab's independent beta load is replaced for this step.
    ///
    /// `rho0`, `cp0`, physical units and accepted populations are never rebound.
    /// Both initial and proposed endpoint material departures must satisfy the
    /// explicit bounds. Neither clamping nor extrapolation is performed.
    /// Failure or cancellation preserves the entire accepted slab. Source
    /// callbacks must not advance external physical history during queries.
    ///
    /// Coefficients/forces are held at the initial state for one explicit step;
    /// endpoint sampling checks support, not convergence. Call this path each
    /// step for ongoing feedback. The neglected density/cp changes must remain
    /// physically small; their bounds are not an error or validity certificate.
    /// This does not add compressible work, species/mass advection, phase change
    /// or gravitational/mechanical energy to the sensible-heat account.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        config: BoussinesqConfig<2>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<BoussinesqMaterial<F>>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            config,
            gate,
            resolve,
        )
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// 3-D counterpart with actual dx³ volumes and D3Q7 thermal relaxation.
    /// Shares all sampling, approximation bounds, force ownership, thermal
    /// accounting and rollback with the 2-D Boussinesq material step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: BoussinesqConfig<3>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<BoussinesqMaterial<F>>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            config,
            gate,
            resolve,
        )
    }
}
