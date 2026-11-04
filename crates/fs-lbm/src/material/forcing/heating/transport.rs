//! Local state-resolved viscosity and conductivity with retained thermal storage.
//!
//! The existing material interface supplies the laws. The existing physical
//! heating adapter supplies SI conversions. The existing thermal/flow kernels
//! evolve the populations. No constitutive interpolation is invented here.

/// Density-resolved buoyancy with explicit, fixed reference storage.
pub mod boussinesq;
/// Source-resolved shear-dependent viscosity on the same thermal/flow path.
pub mod rheology;

use super::{HeatBalance, HeatSourceLimits, MaterialSlab, ThermalLbm, ThermalLbm3, ThermalSourceError};
use crate::material::{MaterialLbmError, plan_scaling_from_material};
use crate::thermal::sources::{Slab, source_step_transport_mapped, transport::RelaxationField};
use fs_exec::CancelGate;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};

/// One complete material evaluation, including the original concrete source
/// owner (and its receipts when source-backed), not merely an anonymous scalar.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportSample<F> {
    /// Original linear fluid-cell index, x fastest. Wall cells are omitted.
    pub cell: usize,
    /// Actual evaluation temperature [K].
    pub temperature_k: f64,
    /// The value returned by the constitutive owner, retained without erasure.
    pub material: F,
    /// Exact single transport snapshot extracted from that value.
    pub transport: FluidTransportState,
    /// Derived momentum relaxation used for this evaluation.
    pub momentum_relaxation: f64,
    /// Derived thermal relaxation, with the actual stencil's conversion.
    pub thermal_relaxation: f64,
}

/// Accepted sensible-heat balance and both source-evaluation epochs.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialTransportStep<F> {
    /// Independent source/storage/wall accounting in joules.
    pub heat: HeatBalance,
    /// Initial-state material evaluations actually used during transport.
    pub applied: Vec<TransportSample<F>>,
    /// Evaluations at the accepted endpoint, checked BEFORE publishing it.
    /// These are domain checks/diagnostics, not a second transport step.
    pub endpoint: Vec<TransportSample<F>>,
}

/// No refusal changes the accepted slab. Externally captured callback receipts
/// remain caller-owned, including records of rejected proposals.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialTransportError {
    /// Invalid layout, wall convention, or unsupported storage change.
    Invalid(&'static str),
    /// Material law or its lattice conversion refused a specific query.
    State {
        /// Original linear cell index.
        cell: usize,
        /// Query temperature [K].
        temperature_k: f64,
        /// Original material/conversion error.
        source: MaterialLbmError,
    },
    /// Local shear rate and material viscosity did not become self-consistent.
    RheologyConvergence {
        /// Original linear fluid-cell index.
        cell: usize,
        /// Material evaluations attempted at this cell and time level.
        iterations: usize,
        /// Final absolute relaxation-time fixed-point residual.
        relaxation_defect: f64,
    },
    /// Existing source, physical-unit or heat-balance admission failed.
    Thermal(ThermalSourceError),
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for MaterialTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material transport: {s}"),
            Self::State { cell, temperature_k, source } => write!(f, "material transport at cell {cell}, {temperature_k} K: {source}"),
            Self::RheologyConvergence { cell, iterations, relaxation_defect } => write!(f,
                "material rheology at cell {cell} did not converge in {iterations} iterations (relaxation defect {relaxation_defect})"),
            Self::Thermal(e) => write!(f, "material transport: {e}"),
            Self::Cancelled => write!(f, "material transport cancelled"),
        }
    }
}
impl std::error::Error for MaterialTransportError {}
impl From<ThermalSourceError> for MaterialTransportError {
    fn from(e: ThermalSourceError) -> Self { Self::Thermal(e) }
}
fn poll(gate: &CancelGate) -> Result<(), MaterialTransportError> {
    if gate.is_requested() { Err(MaterialTransportError::Cancelled) } else { Ok(()) }
}

#[allow(clippy::too_many_arguments)]
pub(in crate::material::forcing) fn sample<S: Slab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    mapping: &MaterialSlab<S>, state: &S, cells: usize, plane: usize,
    range: [f64; 2], gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<Vec<TransportSample<F>>, MaterialTransportError> {
    let temperatures = temperatures::<S, D, Q>(mapping, state, range, gate)?;
    sample_temperatures::<S, F, D, Q>(mapping, &temperatures, cells, plane, range, gate, resolve)
}

fn temperatures<S: Slab<D, Q>, const D: usize, const Q: usize>(
    mapping: &MaterialSlab<S>, state: &S, range: [f64; 2], gate: &CancelGate,
) -> Result<Vec<f64>, MaterialTransportError> {
    // Admit the whole field before calling a material owner; no late bad
    // temperature causes an otherwise avoidable partially evaluated field.
    let mut temperatures = Vec::with_capacity(state.populations().len());
    for (cell, pop) in state.populations().iter().enumerate() {
        if cell % 256 == 0 { poll(gate)?; }
        if !pop.iter().all(|x| x.is_finite()) {
            return Err(MaterialTransportError::Invalid("nonfinite temperature population"));
        }
        let t = mapping.temperature_si(pop.iter().sum()).map_err(ThermalSourceError::from)?;
        if t < range[0] || t > range[1] {
            return Err(MaterialTransportError::Invalid("temperature outside the admitted interval"));
        }
        temperatures.push(t);
    }
    Ok(temperatures)
}

// The same sampler accepts physical midpoint coordinates without fabricating
// temperature populations or resetting the accepted nonequilibrium state.
#[allow(clippy::too_many_arguments)]
pub(in crate::material::forcing) fn sample_temperatures<
    S: Slab<D, Q>, F: FluidTransport, const D: usize, const Q: usize,
>(
    mapping: &MaterialSlab<S>, temperatures: &[f64], cells: usize, plane: usize,
    range: [f64; 2], gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<Vec<TransportSample<F>>, MaterialTransportError> {
    sample_temperatures_with::<S, F, D, Q>(mapping, temperatures, cells, plane, range, gate,
        &mut |cell, temperature_k, retained| resolve(cell, temperature_k, retained)
            .map_err(|source| MaterialTransportError::State {
                cell, temperature_k, source: source.into(),
            }))
}

#[allow(clippy::too_many_arguments)]
fn sample_temperatures_with<S: Slab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    mapping: &MaterialSlab<S>, temperatures: &[f64], cells: usize, plane: usize,
    range: [f64; 2], gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialTransportError>,
) -> Result<Vec<TransportSample<F>>, MaterialTransportError> {
    poll(gate)?;
    if temperatures.len() != cells || plane == 0 || plane > cells / 2
        || !range.iter().all(|t| t.is_finite() && *t > 0.0) || range[0] > range[1]
    {
        return Err(MaterialTransportError::Invalid("invalid material temperature coordinates"));
    }
    for (cell, &t) in temperatures.iter().enumerate() {
        if cell % 256 == 0 { poll(gate)?; }
        if !t.is_finite() || t < range[0] || t > range[1] {
            return Err(MaterialTransportError::Invalid("temperature outside the admitted interval"));
        }
    }
    let units = mapping.force_units;
    let mut samples = Vec::with_capacity(cells-2*plane);
    for cell in plane..cells-plane {
        poll(gate)?;
        let temperature_k = temperatures[cell];
        let error = |source| MaterialTransportError::State { cell, temperature_k, source };
        let material = resolve(cell, temperature_k, &units.transport)?;
        let transport = material.transport_state().map_err(|e| error(e.into()))?;
        if transport.density_kg_m3().to_bits() != units.transport.density_kg_m3().to_bits()
            || transport.specific_heat_j_kg_k().to_bits() != units.transport.specific_heat_j_kg_k().to_bits()
        {
            return Err(MaterialTransportError::Invalid(
                "local density or heat-capacity changes require conservative mass/enthalpy evolution",
            ));
        }
        // Reuse the physical conversion owner on the exact sampled value;
        // do not query the possibly stateful source object a second time.
        let scale = plan_scaling_from_material(&transport, units.cell_size_m, units.timestep_s, 0.0)
            .map_err(error)?;
        let thermal_relaxation = if D == 3 { scale.thermal_relaxation_time3().map_err(error)? }
            else { scale.thermal_relaxation_time() };
        samples.push(TransportSample {
            cell, temperature_k, material, transport,
            momentum_relaxation: scale.flow().tau, thermal_relaxation,
        });
    }
    poll(gate)?;
    Ok(samples)
}

#[allow(clippy::too_many_arguments)]
fn advance<S: Slab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, power: &[f64], forces: &[[f64; D]],
    depth_m: f64, limits: HeatSourceLimits, gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<MaterialTransportStep<F>, MaterialTransportError> {
    advance_with_sampler::<S, F, D, Q>(slab, power, forces, depth_m, limits, gate,
        |mapping, state, cells, plane| sample::<S, F, D, Q>(mapping, state, cells, plane,
            limits.temperature_range_k, gate, &mut resolve))
}

// One preflight, numerical step, endpoint admission and commit for every
// constitutive coordinate set. A rate-aware material is not a second solver.
#[allow(clippy::too_many_arguments)]
fn advance_with_sampler<S: Slab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    gate: &CancelGate,
    sample: impl FnMut(
        &MaterialSlab<S>,
        &S,
        usize,
        usize,
    ) -> Result<Vec<TransportSample<F>>, MaterialTransportError>,
) -> Result<MaterialTransportStep<F>, MaterialTransportError> {
    advance_with_sampler_and_sources::<S, F, (), D, Q>(
        slab,
        power,
        forces,
        depth_m,
        limits,
        gate,
        sample,
        |_, input, _| Ok((input, ())),
    )
    .map(|(step, ())| step)
}

// A physical loss owner may add heat from the SAME already-resolved samples.
// Source construction and its accounting must pass before the candidate is
// stepped; source, endpoint and cancellation refusals all precede publication.
#[allow(clippy::too_many_arguments)]
fn advance_with_sampler_and_sources<
    S: Slab<D, Q>,
    F: FluidTransport,
    A,
    const D: usize,
    const Q: usize,
>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    gate: &CancelGate,
    sample: impl FnMut(
        &MaterialSlab<S>,
        &S,
        usize,
        usize,
    ) -> Result<Vec<TransportSample<F>>, MaterialTransportError>,
    sources: impl FnOnce(
        &MaterialSlab<S>,
        super::Inputs<D>,
        &[TransportSample<F>],
    ) -> Result<(super::Inputs<D>, A), MaterialTransportError>,
) -> Result<(MaterialTransportStep<F>, A), MaterialTransportError> {
    advance_with_sampler_sources_and_candidate::<S, F, A, D, Q>(
        slab,
        power,
        forces,
        depth_m,
        limits,
        gate,
        sample,
        sources,
        |mapping, input, applied| candidate_with_input::<S, F, D, Q>(mapping, input, applied, gate),
    )
}

// Composition seam for a force owner AND a physical heat owner. Both consume
// the same resolved samples, and neither may bypass endpoint admission or
// publish a partial step. The ordinary transport path uses the same driver.
#[allow(clippy::too_many_arguments)]
fn advance_with_sampler_sources_and_candidate<
    S: Slab<D, Q>,
    F: FluidTransport,
    A,
    const D: usize,
    const Q: usize,
>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    gate: &CancelGate,
    mut sample: impl FnMut(
        &MaterialSlab<S>,
        &S,
        usize,
        usize,
    ) -> Result<Vec<TransportSample<F>>, MaterialTransportError>,
    sources: impl FnOnce(
        &MaterialSlab<S>,
        super::Inputs<D>,
        &[TransportSample<F>],
    ) -> Result<(super::Inputs<D>, A), MaterialTransportError>,
    candidate: impl FnOnce(
        &MaterialSlab<S>,
        super::Inputs<D>,
        &[TransportSample<F>],
    ) -> Result<(S, HeatBalance), MaterialTransportError>,
) -> Result<(MaterialTransportStep<F>, A), MaterialTransportError> {
    poll(gate)?;
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
    let applied = sample(slab, &slab.lattice, cells, plane)?;
    let (input, audit) = sources(slab, input, &applied)?;
    let (trial, heat) = candidate(slab, input, &applied)?;
    let endpoint = sample(slab, &trial, cells, plane)?;
    poll(gate)?;
    slab.lattice = trial;
    Ok((
        MaterialTransportStep {
            heat,
            applied,
            endpoint,
        },
        audit,
    ))
}

// Shared by prescribed and electrical heat. The returned solver is tentative;
// material endpoint, circuit accounting and cancellation owners publish later.
#[allow(clippy::too_many_arguments)]
pub(in crate::material::forcing) fn candidate_with_samples<
    S: Slab<D, Q>, F, const D: usize, const Q: usize,
>(
    slab: &MaterialSlab<S>, power: &[f64], forces: &[[f64; D]],
    depth_m: f64, limits: HeatSourceLimits, applied: &[TransportSample<F>],
    gate: &CancelGate,
) -> Result<(S, HeatBalance), MaterialTransportError> {
    poll(gate)?;
    let (cells, _) = slab.lattice.layout()?;
    let input = slab.heat_inputs(power, forces, cells, depth_m, limits)?;
    candidate_with_input::<S, F, D, Q>(slab, input, applied, gate)
}

fn candidate_with_input<S: Slab<D, Q>, F, const D: usize, const Q: usize>(
    slab: &MaterialSlab<S>, input: super::Inputs<D>, applied: &[TransportSample<F>],
    gate: &CancelGate,
) -> Result<(S, HeatBalance), MaterialTransportError> {
    candidate_with_prepared_input::<S, F, D, Q>(slab, input, applied, gate, |_| {})
}

// A material force owner may replace background buoyancy on this tentative
// state. It still uses the same transport kernel, SI balance and rollback.
fn candidate_with_prepared_input<S: Slab<D, Q>, F, const D: usize, const Q: usize>(
    slab: &MaterialSlab<S>, input: super::Inputs<D>, applied: &[TransportSample<F>],
    gate: &CancelGate, prepare: impl FnOnce(&mut S),
) -> Result<(S, HeatBalance), MaterialTransportError> {
    let (cells, plane) = slab.lattice.layout()?;
    if applied.len() != cells - 2 * plane
        || applied.iter().enumerate().any(|(i, s)| s.cell != plane + i)
        || slab.lattice.walls() != (1.0, 0.0) || !slab.lattice.flow_valid()
    {
        return Err(MaterialTransportError::Invalid("invalid material samples or slab"));
    }
    // Wall entries are admitted but unused by the fluid kernels. Do not query
    // a fluid material law for the separate fixed-temperature reservoir.
    let mut momentum = vec![1.0;cells];
    let mut thermal = vec![1.0;cells];
    for value in applied {
        momentum[value.cell] = value.momentum_relaxation;
        thermal[value.cell] = value.thermal_relaxation;
    }
    poll(gate)?;
    // Endpoint source errors/cancellation remain atomic without erasing their
    // original type. As in electrical coupling, original, trial and thermal
    // tentative grids coexist; max_cells is a cell budget, not a byte budget.
    let (mut trial, _) = slab.lattice.snapshot();
    prepare(&mut trial);
    let heat = source_step_transport_mapped::<S, D, Q, _>(
        &mut trial, &input.heat, &input.force, input.limits,
        Some(RelaxationField { momentum: &momentum, thermal: &thermal }),
        |_, raw| input.balance(raw),
    )?;
    Ok((trial, heat))
}

impl MaterialSlab<ThermalLbm> {
    /// Advance 2D heat/flow with local material transport at actual accepted K.
    /// `resolve(cell, kelvin, retained_transport)` returns any existing
    /// `FluidTransport` owner, retained in the result with its concrete receipts.
    /// Pressure, composition, specimen assignments and source policy may be
    /// captured by the callback; no material name or fallback is inferred.
    ///
    /// Power [W/m3] and force [N/m3] are prescribed. Depth [m] sets represented
    /// energy, not 2D dynamics. Viscosity/conductivity vary locally and between
    /// calls; density and Cp must equal the retained snapshot bit-for-bit.
    /// Buoyancy/expansion retain the existing Boussinesq model. Initial-state
    /// coefficients are held for one step; endpoint evaluation checks support,
    /// not implicit convergence or continuous support between query points.
    ///
    /// Failed inputs, source queries, heat checks or cancellation preserve the
    /// slab. Cancellation is checked between queries and around whole kernel
    /// steps, not within every lattice update. Sources must themselves return.
    /// Existing solver-state/stability preconditions remain caller obligations.
    /// To continue local feedback call this method each step; legacy methods
    /// deliberately continue using their unchanged background transport setup.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_transport<F: FluidTransport>(
        &mut self, power_w_m3: &[f64], forces_n_m3: &[[f64;2]], depth_m: f64,
        limits: HeatSourceLimits, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<F>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(self, power_w_m3, forces_n_m3, depth_m, limits, gate, resolve)
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart using retained dx^3 volumes and D3Q7 thermal scaling.
    /// All material evaluation, source/force conversion, transport, endpoint
    /// checks and rollback are shared with the 2D method.
    pub fn step_with_material_transport<F: FluidTransport>(
        &mut self, power_w_m3: &[f64], forces_n_m3: &[[f64;3]],
        limits: HeatSourceLimits, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<F>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(self, power_w_m3, forces_n_m3, depth_m, limits, gate, resolve)
    }
}
