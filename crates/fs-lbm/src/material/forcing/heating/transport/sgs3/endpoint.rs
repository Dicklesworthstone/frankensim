//! Source-density endpoint quadrature composed with the existing SGS owners.

use super::*;
use super::super::boussinesq::BoussinesqEndpointStep;
use crate::thermal::sources::sgs3::admit_endpoint;
use crate::thermal::sources::transport::source_step_resolved_mapped;

/// Three original material query epochs, accepted heat, and the applied SGS.
/// Molecular source records remain distinct from the modeled transport report.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialSgsEndpointStep<F, R> {
    pub material: BoussinesqEndpointStep<F>,
    pub sgs: R,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<F: FluidTransport, R: Default>(
    slab: &mut MaterialSlab<ThermalLbm3>,
    power: &[f64], forces: &[[f64; 3]], limits: HeatSourceLimits,
    buoyancy: BoussinesqConfig<3>, gate: &CancelGate,
    resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    step: impl FnOnce(
        &mut ThermalLbm3, &Inputs<3>, RelaxationField<'_>,
        &mut dyn FnMut(&ThermalLbm3, Vec<f64>) -> Result<Vec<[f64; 3]>, MaterialTransportError>,
    ) -> Result<(HeatBalance, R), MaterialTransportError>,
) -> Result<MaterialSgsEndpointStep<F, R>, MaterialTransportError> {
    let depth = slab.force_units.cell_size_m;
    let mut sgs = R::default();
    let material = boussinesq::advance_with_candidate::<ThermalLbm3, F, 3, 7>(
        slab, power, forces, depth, limits, buoyancy, gate, resolve,
        |mapping, input, applied, mut bounded| {
            boussinesq::endpoint::candidate_with_step::<ThermalLbm3, F, 3, 7>(
                mapping, input, forces, applied, limits.temperature_range_k,
                buoyancy, gate, &mut bounded,
                |candidate, input, molecular, force| {
                    let (heat, report) = step(candidate, input, molecular, force)?;
                    sgs = report;
                    Ok(heat)
                },
            )
        },
    )?;
    Ok(MaterialSgsEndpointStep { material, sgs })
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance local SGS heat/momentum with endpoint-averaged source-density
    /// buoyancy. The original local SGS owner prepares total rates once, from
    /// accepted Guo-centered velocity and freshly sampled molecular mu/k.
    /// Those rates stay fixed for the one actual heat/flow step; scratch storage
    /// is reused, and no augmented rate becomes next step's molecular input.
    ///
    /// After heat transport under those total rates, the original bounded
    /// material owner resolves density at the scalar endpoint predictor. The
    /// existing endpoint projector averages initial/predicted density departures
    /// and adds external force once. Its distinct force epoch and the separately
    /// sampled actual accepted endpoint are retained in `material`.
    /// All three material queries, source/wall heat, endpoint Mach and cancellation
    /// checks precede publication. On refusal only workspace scratch may change.
    ///
    /// Reference rho0/Cp0, resolved walls, explicit SGS usage limits and sensible
    /// storage are unchanged. Neither SGS dissipation nor gravitational work is
    /// added as heat. This is force quadrature with frozen molecular/SGS rates,
    /// not an implicit solve or a global second-order coupled-solver claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_sgs_endpoint_averaged<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64], forces_n_m3: &[[f64; 3]], limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>, config: ThermalSgsConfig,
        workspace: &mut ThermalSgsWorkspace3, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsEndpointStep<F, ThermalSgsReport>, MaterialTransportError> {
        advance(
            self, power_w_m3, forces_n_m3, limits, buoyancy, gate, resolve,
            |candidate, input, molecular, force| {
                // The accepted fext is still bound here. Preparing after a
                // force replacement would corrupt the resolved velocity fit.
                let sgs = candidate.prepare_thermal_sgs(config, Some(molecular), workspace, gate)?;
                let heat = source_step_resolved_mapped::<ThermalLbm3, 3, 7, _, MaterialTransportError>(
                    candidate, &input.heat, &input.force, input.limits,
                    workspace.relaxation(), force,
                    |endpoint, raw| {
                        admit_endpoint(endpoint, config.maximum_mach, gate)?;
                        input.balance(raw).map_err(MaterialTransportError::from)
                    },
                )?;
                Ok((heat, sgs))
            },
        )
    }

    /// Independently fitted plane-channel momentum and heat with the same
    /// endpoint-density quadrature, three source epochs and atomic publication.
    /// Both existing dynamic fits use the accepted state exactly once before
    /// thermal transport. No turbulent Prandtl number is imposed, and a zero
    /// momentum fit does not disable a positive heat fit. Existing scalar wall
    /// laws, separate ratio bounds, workspace byte budget and usage limits apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_channel_sgs_endpoint_averaged<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64], forces_n_m3: &[[f64; 3]], limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>, config: ThermalChannelSgsConfig,
        workspace: &mut ThermalChannelSgsWorkspace3, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsEndpointStep<F, ThermalChannelSgsReport>, MaterialTransportError> {
        advance(
            self, power_w_m3, forces_n_m3, limits, buoyancy, gate, resolve,
            |candidate, input, molecular, force| {
                let sgs = candidate.prepare_thermal_channel_sgs(
                    config, Some(molecular), workspace, gate,
                )?;
                let heat = source_step_resolved_mapped::<ThermalLbm3, 3, 7, _, MaterialTransportError>(
                    candidate, &input.heat, &input.force, input.limits,
                    workspace.relaxation(), force,
                    |endpoint, raw| {
                        admit_endpoint(endpoint, config.maximum_mach, gate)?;
                        input.balance(raw).map_err(MaterialTransportError::from)
                    },
                )?;
                Ok((heat, sgs))
            },
        )
    }
}
