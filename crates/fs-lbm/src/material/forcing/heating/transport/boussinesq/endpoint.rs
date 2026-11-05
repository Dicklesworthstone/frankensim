//! Density quadrature from one transported scalar endpoint prediction.

use super::super::{RelaxationField, TransportSample, candidate_with_transport, sample_temperatures};
use super::{
    BackgroundBuoyancy, BoussinesqConfig, BoussinesqMaterial, FluidTransport, FluidTransportState,
    HeatSourceLimits, MaterialError, MaterialSlab, MaterialTransportError, MaterialTransportStep,
    Slab, ThermalLbm, ThermalLbm3, advance, poll, project_force,
};
use crate::thermal::sources::transport::source_step_resolved_mapped;
use crate::thermal::sources::{HeatBalance, ThermalSourceError};
use fs_exec::CancelGate;

/// One accepted material step and the distinct source epoch used by buoyancy.
#[derive(Debug, Clone, PartialEq)]
pub struct BoussinesqEndpointStep<F> {
    /// Original initial/actual-endpoint source samples and sensible-heat audit.
    /// `endpoint` is queried from the accepted population sums, not relabelled
    /// from the pre-flow prediction, even when their temperatures differ by ulps.
    pub transport: MaterialTransportStep<BoussinesqMaterial<F>>,
    /// Original bounded material queries at `T_transport + heat_per_step/2`,
    /// converted through the retained temperature coordinate. This predicts
    /// the complete scalar endpoint without executing the final source kick
    /// early. Its densities, averaged with the initial densities, drive flow.
    /// Its viscosity/conductivity are checked but do not replace initial rates.
    pub force_endpoint: Vec<TransportSample<BoussinesqMaterial<F>>>,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn candidate<
    S: Slab<D, Q> + BackgroundBuoyancy,
    F: FluidTransport,
    const D: usize,
    const Q: usize,
>(
    slab: &MaterialSlab<S>,
    input: super::super::super::Inputs<D>,
    forces: &[[f64; D]],
    applied: &[TransportSample<BoussinesqMaterial<F>>],
    range_k: [f64; 2],
    config: BoussinesqConfig<D>,
    gate: &CancelGate,
    resolve: &mut impl FnMut(
        usize,
        f64,
        &FluidTransportState,
    ) -> Result<BoussinesqMaterial<F>, MaterialError>,
) -> Result<(S, HeatBalance, Vec<TransportSample<BoussinesqMaterial<F>>>), MaterialTransportError> {
    candidate_with_step::<S, F, D, Q>(
        slab, input, forces, applied, range_k, config, gate, resolve,
        |trial, input, molecular, force| {
            source_step_resolved_mapped::<S, D, Q, _, MaterialTransportError>(
                trial, &input.heat, &input.force, input.limits, molecular, force,
                |_, raw| input.balance(raw).map_err(MaterialTransportError::from),
            )
        },
    )
}

// SGS augments fresh molecular rates, then consumes this exact density
// projector through the original one-heat/one-flow source transaction.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(in crate::material::forcing::heating::transport) fn candidate_with_step<
    S: Slab<D, Q> + BackgroundBuoyancy, F: FluidTransport, const D: usize, const Q: usize,
>(
    slab: &MaterialSlab<S>, mut input: super::super::super::Inputs<D>, forces: &[[f64; D]],
    applied: &[TransportSample<BoussinesqMaterial<F>>], range_k: [f64; 2],
    config: BoussinesqConfig<D>, gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<BoussinesqMaterial<F>, MaterialError>,
    step: impl FnOnce(
        &mut S, &super::super::super::Inputs<D>, RelaxationField<'_>,
        &mut dyn FnMut(&S, Vec<f64>) -> Result<Vec<[f64; D]>, MaterialTransportError>,
    ) -> Result<HeatBalance, MaterialTransportError>,
) -> Result<(S, HeatBalance, Vec<TransportSample<BoussinesqMaterial<F>>>), MaterialTransportError> {
    project_force::<S, BoussinesqMaterial<F>, D, Q>(
        slab,
        &mut input,
        forces,
        applied,
        config,
        gate,
        |material| material.resolved.density_kg_m3(),
    )?;
    let units = slab.force_units;
    let mut force_endpoint = Vec::new();
    let (mut trial, heat) = candidate_with_transport::<S, BoussinesqMaterial<F>, D, Q>(
        slab,
        input,
        applied,
        gate,
        |candidate| candidate.set_coefficient(0.0),
        |trial, input, molecular| {
            step(
                trial,
                input,
                molecular,
                &mut |transported, mut predicted| {
                    let (cells, plane) = transported.layout()?;
                    for value in &mut predicted[plane..cells - plane] {
                        *value = slab
                            .temperature_si(*value)
                            .map_err(ThermalSourceError::from)?;
                    }
                    let endpoint = sample_temperatures::<S, BoussinesqMaterial<F>, D, Q>(
                        slab, &predicted, cells, plane, range_k, gate, resolve,
                    )?;
                    let mut resolved_forces = input.force.clone();
                    for (initial, next) in applied.iter().zip(&endpoint) {
                        poll(gate)?;
                        let rho0 = units.transport.density_kg_m3();
                        let departure = 0.5 * (initial.material.resolved.density_kg_m3() - rho0)
                            + 0.5 * (next.material.resolved.density_kg_m3() - rho0);
                        for axis in 0..D {
                            let force =
                                forces[next.cell][axis] + departure * config.gravity_m_s2[axis];
                            resolved_forces[next.cell][axis] =
                                units.component(force).map_err(|source| {
                                    MaterialTransportError::State {
                                        cell: next.cell,
                                        temperature_k: next.temperature_k,
                                        source,
                                    }
                                })?;
                        }
                    }
                    poll(gate)?;
                    force_endpoint = endpoint;
                    Ok(resolved_forces)
                },
            )
        },
    )?;
    trial.set_coefficient(slab.lattice.coefficient());
    Ok((trial, heat, force_endpoint))
}

impl MaterialSlab<ThermalLbm> {
    /// Advance one native heat/flow step with endpoint-averaged source-density
    /// buoyancy. Initial source mu/k set local relaxation; reference rho0/cp0
    /// continue to own inertia and storage. No beta is inferred or added.
    ///
    /// After the only heat transport, query the original bounded source at
    /// `T_transport + heat_per_step/2`. Average its density departure with the
    /// original initial departure, add the prescribed external force once and
    /// advance flow once in the existing Guo frame. The last equilibrium heat
    /// half-kick retains the actual new velocity and advective-flux history.
    /// The extra `force_endpoint` epoch is explicit; the ordinary report still
    /// samples the actual accepted endpoint separately. All three source epochs,
    /// wall/storage checks and cancellation gates precede publication.
    ///
    /// This is a density-force quadrature on a scalar endpoint predictor, not
    /// an implicit coupled solve or a global second-order claim. Mu/k remain
    /// initial-state coefficients, and no gravitational/mechanical work enters
    /// the sensible-heat balance. Original Boussinesq bounds remain required.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_endpoint_averaged<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        config: BoussinesqConfig<2>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqEndpointStep<F>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            config,
            true,
            gate,
            resolve,
        )
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart retaining dx³ storage and the original
    /// D3Q19 collision selection. Uses the same one-pass heat/flow sequencing,
    /// distinct force/accepted endpoint epochs, bounds and atomic publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_endpoint_averaged<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: BoussinesqConfig<3>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqEndpointStep<F>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            config,
            true,
            gate,
            resolve,
        )
    }
}
