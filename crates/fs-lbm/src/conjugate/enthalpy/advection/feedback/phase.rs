//! Phase-first material sampling and exact stationary solid constraints.
//!
//! The existing equilibrium h chart owns phase and latent energy. This module
//! changes only momentum occupancy; caloric cells, source snapshots, electrical
//! paths and conductive links survive both freezing and melting.

use fs_material::phase::transport::PhaseTransport;

use super::*;
use crate::conjugate::{FvBoundary, Voxel};

/// Bounded source-density convection through exact solid, mushy and liquid states.
pub mod buoyancy;

/// Transport information consumed by the common material-flow owner.
///
/// The original liquid API uses `FluidTransportState` directly. The phase API
/// retains conductivity for every cell, with no liquid tuple for exact solids.
pub trait FlowTransportSnapshot: Clone {
    fn mobile_transport(&self) -> Option<FluidTransportState>;
    fn thermal_conductivity_w_m_k(&self) -> f64;
    fn phase_mobility(&self) -> Option<&PhaseMobility> {
        None
    }
}

impl FlowTransportSnapshot for FluidTransportState {
    fn mobile_transport(&self) -> Option<FluidTransportState> {
        Some(*self)
    }
    fn thermal_conductivity_w_m_k(&self) -> f64 {
        FluidTransportState::thermal_conductivity_w_m_k(self)
    }
}

/// One phase/source epoch, without fabricated solid viscosity or heat capacity.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseFlowTransport {
    mobility: PhaseMobility,
    mobile: Option<FluidTransportState>,
    conductivity_w_m_k: f64,
}

impl PhaseFlowTransport {
    pub const fn mobility(&self) -> &PhaseMobility {
        &self.mobility
    }
    pub const fn mobile_transport(&self) -> Option<FluidTransportState> {
        self.mobile
    }
    pub const fn thermal_conductivity_w_m_k(&self) -> f64 {
        self.conductivity_w_m_k
    }
}

impl FlowTransportSnapshot for PhaseFlowTransport {
    fn mobile_transport(&self) -> Option<FluidTransportState> {
        self.mobile
    }
    fn thermal_conductivity_w_m_k(&self) -> f64 {
        self.conductivity_w_m_k
    }
    fn phase_mobility(&self) -> Option<&PhaseMobility> {
        Some(&self.mobility)
    }
}

pub(super) fn fluid_snapshot<F: FluidTransport>(
    source: &F,
    _: Option<PhaseMobility>,
) -> Result<FluidTransportState, MaterialError> {
    source.transport_state()
}

fn phase_snapshot<F: PhaseTransport>(
    source: &F,
    mobility: Option<PhaseMobility>,
) -> Result<PhaseFlowTransport, MaterialError> {
    let invalid = |what: &str| MaterialError::State { what: what.into() };
    let mobility = mobility.ok_or_else(|| invalid("phase flow lacks a classified h state"))?;
    let phase = match &mobility {
        PhaseMobility::ImmobileSolid { phase, .. } => phase,
        PhaseMobility::Mobile(state) => &state.phase,
    };
    if source.phase_state() != phase {
        return Err(invalid(
            "phase transport source differs from the original caloric state",
        ));
    }
    if source.transport_density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits() {
        return Err(invalid(
            "phase source density needs an explicit reference-mass projection or volume-change owner",
        ));
    }
    let conductivity_w_m_k = source.thermal_conductivity_w_m_k()?;
    if !conductivity_w_m_k.is_finite() || conductivity_w_m_k <= 0.0 {
        return Err(invalid(
            "phase source conductivity must be positive and finite",
        ));
    }
    // This branch is deliberately before transport_state(): liquid-only
    // source properties need not exist at a cold, completely solid state.
    let mobile = match mobility {
        PhaseMobility::ImmobileSolid { .. } => None,
        PhaseMobility::Mobile(_) => {
            let transport = source.transport_state()?;
            if transport.thermal_conductivity_w_m_k().to_bits() != conductivity_w_m_k.to_bits() {
                return Err(invalid(
                    "mobile phase conductivity differs from its transport source",
                ));
            }
            Some(transport)
        }
    };
    Ok(PhaseFlowTransport {
        mobility,
        mobile,
        conductivity_w_m_k,
    })
}

impl<'a, F: PhaseTransport + Clone, R, O: DepositionObserver>
    MaterialFlowFeedback<'a, F, R, O, PhaseFlowTransport>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Continue fixed-volume material through exact solid, mushy and liquid states.
    ///
    /// All original `domain` fluid voxels remain mapped caloric/material cells.
    /// At each source and endpoint, the SAME registered equilibrium h chart
    /// classifies mobility before the resolver runs. Exact solid cells become
    /// native stationary no-slip voxels; every positive admitted liquid fraction
    /// retains the declared Carman–Kozeny drag with its actual liquid viscosity.
    /// No epsilon permeability, threshold locking or fictitious liquid tuple is
    /// introduced. Freezing and remelting rebuild the same native face topology.
    ///
    /// The original source snapshots supply conductivity in every phase and
    /// remain available to the existing current adapter. Source k on full voxel
    /// faces is attached with `Self::with_fluid_conduction`; its name refers to
    /// the original potential-flow volume, not to the current mobility mask.
    /// The existing enthalpy owner alone stores mass*h, including latent energy.
    ///
    /// Equal constant phase densities, a fixed stationary skeleton and closed
    /// walls/symmetry are mandatory. Unsupported liquid fractions, unequal
    /// phase densities, openings, fans and authored extra resistance refuse.
    /// A moving wall touching frozen material or a body force on a constrained
    /// voxel also refuses: no omitted solid friction/work is silently removed.
    /// An all-solid zero-drive epoch has exactly zero flow and heat without
    /// a momentum solve or invented viscosity. Its report has zero iterations.
    ///
    /// This retains the existing homogeneous single-velocity enthalpy/porosity
    /// approximation. It does not supply moving solid mechanics, residual stress,
    /// shrinkage, separate phase velocities or free surfaces. Retained native
    /// momentum and its physical/numerical energy account may be selected with
    /// `with_unsteady_momentum` and the unsteady heating attachments. Newly
    /// frozen native faces must already have exactly zero speed; otherwise an
    /// explicit momentum/kinetic-energy handoff is still required. Source
    /// support and all enclosing histories keep the original
    /// transactional `advance_observed` refusal/publication boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_phase(
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        input: FlowEnthalpyInput<'a>,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        law: &'a MushyPhaseLaw<'a>,
        mut resolve: R,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(gate, &mut || Ok(()))?;
        admit_phase_config(config)?;
        if max_solves == 0 || domain.fluid_count() == 0 {
            return Err(invalid(
                "phase flow needs a positive solve budget and mapped material volume",
            )
            .into());
        }
        let density_kg_m3 = law.fixed_density_kg_m3().ok_or(invalid(
            "phase flow needs equal constant phase density without a volume source",
        ))?;
        if !law.admits_immobile_solid()
            || !core::ptr::addr_eq(
                input.caloric as *const dyn fs_conduction::transient::deposition::ThermalLaw,
                law.curve() as *const fs_material::phase::EquilibriumEnthalpyPhaseCurve,
            )
        {
            return Err(invalid(
                "phase flow needs the original caloric object and exact-solid mobility law",
            )
            .into());
        }
        admit_map(thermal, domain, density_kg_m3, input, &mut || Ok(()))?;
        let endpoint = sample(
            thermal,
            domain,
            input,
            gate,
            &mut resolve,
            Some(law),
            phase_snapshot::<F>,
            &mut || Ok(()),
        )?;
        admit_density(&endpoint, density_kg_m3)?;
        let coupled = Self {
            domain,
            input,
            config,
            gate,
            resolve,
            reference: None,
            density_kg_m3,
            snapshot: phase_snapshot::<F>,
            exact_phase: true,
            support_capture_endpoint_s: None,
            buoyancy: None,
            external_force_n_m3: None,
            mushy: Some(law),
            phase_conduction: None,
            viscous_heating: None,
            unsteady: None,
            rheology: None,
            bound_time_s: thermal.time_s(),
            fluid_conduction: BTreeMap::new(),
            advected: Some(Advected::new(inner, Vec::new(), input.max_connections)?),
            applied: None,
            endpoint,
            audit: AdvectionAudit::default(),
            solves: 0,
            max_solves,
        };
        coupled.mushy_drag(&coupled.endpoint, &mut || Ok(()))?;
        let active = coupled.momentum_domain(&coupled.endpoint, &mut || Ok(()))?;
        coupled.check_phase_drives(active.as_ref().unwrap_or(domain), &mut || Ok(()))?;
        Ok(coupled)
    }
}

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    pub(super) fn momentum_domain(
        &self,
        samples: &[FlowMaterialSample<F, T>],
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<VoxelDomain>, MaterialFlowError<O::Error>> {
        if !self.exact_phase {
            return Ok(None);
        }
        let mut active = self.domain.clone();
        for sample in samples {
            check_gate(self.gate, check)?;
            if sample.transport.mobile_transport().is_none() {
                let [x, y, z] = active.coords(sample.voxel);
                // The tag never enters caloric/material selection: this clone
                // belongs only to native momentum and its face-flux graph.
                active.set(x, y, z, Voxel::Solid(0));
            }
        }
        Ok(Some(active))
    }

    pub(super) fn check_phase_drives(
        &self,
        active: &VoxelDomain,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        if !self.exact_phase {
            return Ok(());
        }
        admit_phase_config(self.config)?;
        for (side, face) in self.config.faces.iter().enumerate() {
            if let FvBoundary::Wall { velocity } = face {
                if velocity.iter().all(|&v| v == 0.0) {
                    continue;
                }
                if active.fluid_count() == 0 {
                    return Err(invalid(
                        "moving-wall work on an all-solid phase volume is unsupported",
                    )
                    .into());
                }
                for c in 0..active.cell_count() {
                    check_gate(self.gate, check)?;
                    if self.domain.is_fluid(c)
                        && !active.is_fluid(c)
                        && active.neighbor(c, side).is_none()
                    {
                        return Err(invalid(
                            "a moving wall contacts frozen material without a solid contact owner",
                        )
                        .into());
                    }
                }
            }
        }
        if let Some(forces) = &self.external_force_n_m3 {
            for (c, force) in forces.iter().enumerate() {
                check_gate(self.gate, check)?;
                if !active.is_fluid(c) && force.iter().any(|&v| v != 0.0) {
                    return Err(
                        invalid("prescribed body force acts on a constrained solid voxel").into(),
                    );
                }
            }
        }
        Ok(())
    }
}

fn admit_phase_config(config: &SimpleConfig) -> Result<(), DepositionError> {
    if config.turbulence != Turbulence::Laminar
        || config.fan.is_some()
        || !config.internal_fans.is_empty()
        || !config.resistances.is_empty()
        || config
            .faces
            .iter()
            .any(|face| !matches!(face, FvBoundary::Wall { .. } | FvBoundary::Symmetry))
        || !(config.velocity_relaxation > 0.0 && config.velocity_relaxation < 1.0)
        || [
            config.tolerance,
            config.momentum_tolerance,
            config.pressure_tolerance,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(invalid(
            "phase continuation needs a closed laminar wall/symmetry domain without fans or configured resistance",
        ));
    }
    for (side, face) in config.faces.iter().enumerate() {
        if let FvBoundary::Wall { velocity } = face {
            if velocity.iter().any(|v| !v.is_finite()) || velocity[side / 2] != 0.0 {
                return Err(invalid(
                    "phase-domain walls must have finite tangential velocity",
                ));
            }
        }
    }
    Ok(())
}

/// Exact empty constrained system, only after phase/drive admission.
pub(super) fn immobile_flow(domain: &VoxelDomain) -> FvFlow {
    FvFlow::stationary(domain)
}

pub(super) fn immobile_power(domain: &VoxelDomain, darcy: bool, body: bool) -> ViscousPower {
    ViscousPower::stationary(domain, darcy, body)
}
