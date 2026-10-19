//! Local material feedback on the original conservative enthalpy timeline.
//!
//! Each actual thermal substep resolves material transport at its cells' own
//! enthalpy-derived temperatures, solves the original steady MAC momentum
//! system with that viscosity field, and advects h with those exact face fluxes.
//! This is a quasisteady, constant-density liquid approximation: momentum
//! equilibrates between thermal substeps. It does not advance stored momentum.

use std::collections::BTreeMap;

use fs_conduction::transient::deposition::advection::{Advected, AdvectionAudit};
use fs_conduction::transient::deposition::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer,
};
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::{FlowEnthalpyInput, admit_map, invalid, mass_graph};
use crate::conjugate::{
    ChtError, FluidProperties, FvFlow, SimpleConfig, Turbulence, VoxelDomain,
    simple_flow_with_viscosity,
};

/// A coupled refusal retains the originating material, flow or history error.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialFlowError<E> {
    Thermal(DepositionError),
    Material { cell: CellId, source: MaterialError },
    Momentum(ChtError),
    Inner(E),
}

impl<E> From<DepositionError> for MaterialFlowError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<ChtError> for MaterialFlowError<E> {
    fn from(error: ChtError) -> Self {
        if error == ChtError::Cancelled {
            Self::Thermal(DepositionError::Cancelled)
        } else {
            Self::Momentum(error)
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for MaterialFlowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "material flow thermal owner: {error}"),
            Self::Material { cell, source } => write!(f, "material flow cell {cell:?}: {source}"),
            Self::Momentum(error) => write!(f, "material flow momentum owner: {error}"),
            Self::Inner(error) => write!(f, "material flow coupled history: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MaterialFlowError<E> {}

/// One concrete source result and the exact local query/transport snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowMaterialSample<F> {
    pub voxel: usize,
    pub cell: CellId,
    pub specific_enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub material: F,
    pub transport: FluidTransportState,
}

/// The source epoch that produced the LAST applied face-flux field.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialFlowEpoch<F> {
    pub time_s: f64,
    pub samples: Vec<FlowMaterialSample<F>>,
    pub flow: FvFlow,
}

/// Live local viscosity/momentum/enthalpy coupling around any existing observer.
///
/// `resolve(voxel, actual_cell, actual_kelvin)` returns an existing concrete
/// `FluidTransport` owner, including its source receipts. This receives h as
/// well as T, so a caller can explicitly admit a phase-dependent closure rather
/// than caching properties at a latent plateau's constant temperature. The
/// callback owns its material selection, fixed pressure/composition and phase
/// mobility admission. It must not mutate shared external physical history.
///
/// Density must remain bitwise equal to the initial snapshot. T(h), mass and
/// conductivity/contact evolution stay with the actual registered caloric law
/// and wrapped observer; the source Cp is never substituted for that law.
/// In particular this does not turn a solid/liquid fraction into a porous drag
/// law, or assume that a partially solid material can move like a liquid.
///
/// Each source epoch performs one laminar, quasisteady full-stress momentum
/// solve. The shared `Advected` owner contributes its actual face mass rates,
/// boundary h flux, internal transfers and combined conductive/advective CFL.
/// Endpoint source checks run before publication. The last flow's time is its
/// explicit starting epoch, not a claimed momentum solve at the returned T.
/// Refine the thermal step and establish rapid viscous equilibration relative
/// to thermal evolution before using this approximation. There is no kinetic
/// energy storage, pressure-work heat, viscous heating or free-surface motion.
#[derive(Clone)]
pub struct MaterialFlowFeedback<'a, F, R, O> {
    domain: &'a VoxelDomain,
    input: FlowEnthalpyInput<'a>,
    config: &'a SimpleConfig,
    gate: &'a CancelGate,
    resolve: R,
    reference: FluidTransportState,
    advected: Option<Advected<O>>,
    applied: Option<MaterialFlowEpoch<F>>,
    endpoint: Vec<FlowMaterialSample<F>>,
    audit: AdvectionAudit,
    solves: usize,
    max_solves: usize,
}

impl<'a, F: FluidTransport + Clone, R, O: DepositionObserver> MaterialFlowFeedback<'a, F, R, O>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Bind one existing timeline/domain and bound momentum solves per advance.
    /// The map and reservoirs use the same admission as frozen solved flow.
    /// Each thermal substep costs one steady solve; `max_solves` is separate
    /// from the original SIMPLEC iterations and thermal substep budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        input: FlowEnthalpyInput<'a>,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        mut resolve: R,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(gate, &mut || Ok(()))?;
        if max_solves == 0 || config.turbulence != Turbulence::Laminar {
            return Err(invalid(
                "live material flow needs laminar momentum and a positive solve budget",
            )
            .into());
        }
        let endpoint = sample(thermal, domain, input, gate, &mut resolve, &mut || Ok(()))?;
        let reference = endpoint
            .first()
            .ok_or(invalid("live material flow needs an active fluid cell"))?
            .transport;
        admit_density(&endpoint, reference.density_kg_m3())?;
        admit_map(
            thermal,
            domain,
            reference.density_kg_m3(),
            input,
            &mut || Ok(()),
        )?;
        let advected = Advected::new(inner, Vec::new(), input.max_connections)?;
        Ok(Self {
            domain,
            input,
            config,
            gate,
            resolve,
            reference,
            advected: Some(advected),
            applied: None,
            endpoint,
            audit: AdvectionAudit::default(),
            solves: 0,
            max_solves,
        })
    }

    /// Source-backed inputs and actual MAC result for the last applied step.
    pub const fn applied(&self) -> Option<&MaterialFlowEpoch<F>> {
        self.applied.as_ref()
    }
    /// Material evaluations at the accepted thermal endpoint, including h.
    pub fn endpoint(&self) -> &[FlowMaterialSample<F>] {
        &self.endpoint
    }
    /// Summed original advection accounting over the last successful advance.
    /// Its boundary energy is already in the thermal audit; never add it twice.
    pub const fn audit(&self) -> AdvectionAudit {
        self.audit
    }
    /// Number of actual momentum solves in the last successful advance.
    pub const fn momentum_solves(&self) -> usize {
        self.solves
    }
    /// The original accepted source, contact and observation history.
    pub fn inner(&self) -> &O {
        self.advected.as_ref().expect("accepted observer").inner()
    }
    pub fn into_inner(self) -> O {
        self.advected.expect("accepted observer").into_inner()
    }

    fn check_map(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        check_gate(self.gate, check)?;
        admit_map(
            thermal,
            self.domain,
            self.reference.density_kg_m3(),
            self.input,
            check,
        )?;
        Ok(())
    }
}

impl<F: FluidTransport + Clone, R, O: DepositionObserver> DepositionObserver
    for MaterialFlowFeedback<'_, F, R, O>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    type Error = MaterialFlowError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("unfinished material-flow source epoch"))?
            .begin()
            .map_err(MaterialFlowError::Inner)?;
        self.audit = AdvectionAudit::default();
        self.solves = 0;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_map(thermal, check)?;
        if self.solves >= self.max_solves {
            return Err(DepositionError::Budget("live material momentum solves").into());
        }
        let samples = sample(
            thermal,
            self.domain,
            self.input,
            self.gate,
            &mut self.resolve,
            check,
        )?;
        admit_density(&samples, self.reference.density_kg_m3())?;
        let mut viscosities =
            vec![self.reference.dynamic_viscosity_pa_s(); self.domain.cell_count()];
        for value in &samples {
            viscosities[value.voxel] = value.transport.dynamic_viscosity_pa_s();
        }
        // All caloric evolution remains in the thermal law. These other
        // fields satisfy the original flow owner's parameter admission only.
        let fluid = FluidProperties {
            density_kg_m3: self.reference.density_kg_m3(),
            specific_heat_j_kg_k: self.reference.specific_heat_j_kg_k(),
            conductivity_w_m_k: self.reference.thermal_conductivity_w_m_k(),
            kinematic_viscosity_m2_s: self.reference.kinematic_viscosity_m2_s(),
        };
        let flow =
            simple_flow_with_viscosity(self.domain, &fluid, &viscosities, self.config, self.gate)?;
        check_gate(self.gate, check)?;
        let (flows, boundaries) = mass_graph(
            self.domain,
            &flow.field,
            self.reference.density_kg_m3(),
            self.input,
            check,
        )?;
        let old = self
            .advected
            .take()
            .ok_or(invalid("unfinished material-flow source epoch"))?;
        // Rebuilding the graph keeps the accepted wrapped physical history.
        // It does not call begin() again or reset another owner's work counter.
        let mut advected = Advected::new(old.into_inner(), flows, self.input.max_connections)?
            .with_boundaries(boundaries, self.input.max_connections)?;
        let sources = advected
            .sources(thermal, check)
            .map_err(MaterialFlowError::Inner)?;
        self.advected = Some(advected);
        self.applied = Some(MaterialFlowEpoch {
            time_s: thermal.time_s(),
            samples,
            flow,
        });
        self.solves += 1;
        Ok(sources)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .conductances(thermal, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .heat_transfers(thermal, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.check_map(thermal, &mut || Ok(()))?;
        self.endpoint = sample(
            thermal,
            self.domain,
            self.input,
            self.gate,
            &mut self.resolve,
            &mut || Ok(()),
        )?;
        admit_density(&self.endpoint, self.reference.density_kg_m3())?;
        self.advected
            .as_mut()
            .ok_or(invalid("unfinished material-flow source epoch"))?
            .event(event, thermal)
            .map_err(MaterialFlowError::Inner)
    }

    fn step(
        &mut self,
        start_s: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_map(thermal, check)?;
        self.endpoint = sample(
            thermal,
            self.domain,
            self.input,
            self.gate,
            &mut self.resolve,
            check,
        )?;
        admit_density(&self.endpoint, self.reference.density_kg_m3())?;
        let advected = self
            .advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?;
        advected
            .step(start_s, dt, old, thermal, check)
            .map_err(MaterialFlowError::Inner)?;
        add_audit(&mut self.audit, advected.audit())?;
        check_gate(self.gate, check)?;
        Ok(())
    }
}

fn check_gate(
    gate: &CancelGate,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    check()?;
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

fn sample<F: FluidTransport, R, E>(
    thermal: &DepositionSimulation<'_>,
    domain: &VoxelDomain,
    input: FlowEnthalpyInput<'_>,
    gate: &CancelGate,
    resolve: &mut R,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<FlowMaterialSample<F>>, MaterialFlowError<E>>
where
    R: FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    check_gate(gate, check)?;
    if domain.cell_count() > input.max_cells {
        return Err(DepositionError::Budget("live material flow voxels").into());
    }
    if input.cell_ids.len() != domain.cell_count() {
        return Err(invalid("live material flow cell map differs from its domain").into());
    }
    let mut samples = Vec::with_capacity(domain.fluid_count());
    for (voxel, id) in input.cell_ids.iter().enumerate() {
        if !domain.is_fluid(voxel) {
            continue;
        }
        check_gate(gate, check)?;
        let cell = id
            .and_then(|id| thermal.cells().get(&id))
            .ok_or(invalid("every flowing voxel needs an active enthalpy cell"))?;
        let temperature_k = thermal.temperature_k(cell.id)?;
        let material =
            resolve(voxel, cell, temperature_k).map_err(|source| MaterialFlowError::Material {
                cell: cell.id,
                source,
            })?;
        let transport =
            material
                .transport_state()
                .map_err(|source| MaterialFlowError::Material {
                    cell: cell.id,
                    source,
                })?;
        samples.push(FlowMaterialSample {
            voxel,
            cell: cell.id,
            specific_enthalpy_j_kg: cell.specific_h_j_kg,
            temperature_k,
            material,
            transport,
        });
    }
    check_gate(gate, check)?;
    Ok(samples)
}

fn admit_density<F>(samples: &[FlowMaterialSample<F>], rho: f64) -> Result<(), DepositionError> {
    if samples
        .iter()
        .any(|value| value.transport.density_kg_m3().to_bits() != rho.to_bits())
    {
        Err(invalid(
            "live flow density changes require compatible mass/volume and momentum evolution",
        ))
    } else {
        Ok(())
    }
}

fn add_audit(total: &mut AdvectionAudit, step: AdvectionAudit) -> Result<(), DepositionError> {
    for (sum, value) in [
        (&mut total.transported_mass_kg, step.transported_mass_kg),
        (
            &mut total.absolute_enthalpy_transfer_j,
            step.absolute_enthalpy_transfer_j,
        ),
        (&mut total.inlet_mass_kg, step.inlet_mass_kg),
        (&mut total.outlet_mass_kg, step.outlet_mass_kg),
        (&mut total.external_enthalpy_j, step.external_enthalpy_j),
    ] {
        *sum += value;
        if !sum.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
    }
    total.maximum_courant = total.maximum_courant.max(step.maximum_courant);
    total.substeps = total
        .substeps
        .checked_add(step.substeps)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}
