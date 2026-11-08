//! Original molecular source epochs composed with pending turbulence history.

use super::*;
use crate::thermal::sources::sgs3::admit_endpoint;
use crate::thermal::sources::sgs3::lagrangian::{
    ThermalLagrangianSgsConfig, ThermalLagrangianSgsReport, ThermalLagrangianSgsWorkspace3,
    step_lagrangian_sgs_mapped,
};
use crate::thermal::sources::transport::source_step_resolved_mapped;

impl MaterialSlab<ThermalLbm3> {
    /// Advance source-resolved molecular viscosity/conduction and independent
    /// Lagrangian momentum/heat fits in one material transaction.
    ///
    /// The fit reads the accepted Guo-centered state before newly authored
    /// forces are applied. Molecular mu(T), k(T), density and Cp retain their
    /// original provider, SI queries and source admission. Total relaxation is
    /// step-local; it never becomes next step's molecular source. Both source
    /// epochs, original wall heat, sensible balance, Mach and cancellation must
    /// accept before either pathline history publishes. A late endpoint source
    /// failure therefore permits retry at the same sample_time_lu.
    ///
    /// The caller supplies actual sampling time in lattice timesteps, using the
    /// same physical trajectory and time scale as the configured memory. This
    /// retains the native low-Mach, fixed-reference sensible-storage scope;
    /// modeled kinetic transfer is not deposited as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_transport_and_lagrangian_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        sample_time_lu: f64,
        config: ThermalLagrangianSgsConfig,
        workspace: &mut ThermalLagrangianSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsTransportStep<F, ThermalLagrangianSgsReport>, MaterialTransportError>
    {
        let accepted = advance(
            self,
            power_w_m3,
            forces_n_m3,
            limits,
            gate,
            resolve,
            false,
            |_, _, _| Ok(()),
            move |candidate, input, molecular| {
                step_lagrangian_sgs_mapped(
                    candidate,
                    &input.heat,
                    &input.force,
                    sample_time_lu,
                    config,
                    Some(molecular),
                    workspace,
                    input.limits,
                    gate,
                    |raw| input.balance(raw),
                )
            },
        )?;
        let sgs = accepted.sgs.report();
        accepted.sgs.commit();
        Ok(MaterialSgsTransportStep {
            material: accepted.material,
            sgs,
        })
    }

    /// Source-density Boussinesq forcing and joint pathline turbulence share
    /// the same freshly sampled molecular material and accepted velocity.
    /// Original rho0/Cp0 storage and authored density/Cp bounds are retained;
    /// source-density departure replaces the background beta force once.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_lagrangian_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>,
        sample_time_lu: f64,
        config: ThermalLagrangianSgsConfig,
        workspace: &mut ThermalLagrangianSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        MaterialSgsTransportStep<BoussinesqMaterial<F>, ThermalLagrangianSgsReport>,
        MaterialTransportError,
    > {
        let accepted = advance_boussinesq(
            self,
            power_w_m3,
            forces_n_m3,
            limits,
            buoyancy,
            gate,
            resolve,
            move |candidate, input, molecular| {
                step_lagrangian_sgs_mapped(
                    candidate,
                    &input.heat,
                    &input.force,
                    sample_time_lu,
                    config,
                    Some(molecular),
                    workspace,
                    input.limits,
                    gate,
                    |raw| input.balance(raw),
                )
            },
        )?;
        let sgs = accepted.sgs.report();
        accepted.sgs.commit();
        Ok(MaterialSgsTransportStep {
            material: accepted.material,
            sgs,
        })
    }

    /// Endpoint-averaged source-density buoyancy with the same joint history.
    ///
    /// All three original material epochs remain: initial transport, thermal
    /// endpoint force quadrature, and actual accepted endpoint. The history
    /// token survives every query and native/outer check. No intermediate
    /// velocity or density probe advances a turbulence clock, and the fit is
    /// never recomputed in the newly authored force frame. This is the native
    /// force quadrature, not an implicit or globally second-order thermal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_lagrangian_sgs_endpoint_averaged<
        F: FluidTransport,
    >(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>,
        sample_time_lu: f64,
        config: ThermalLagrangianSgsConfig,
        workspace: &mut ThermalLagrangianSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsEndpointStep<F, ThermalLagrangianSgsReport>, MaterialTransportError>
    {
        let accepted = super::endpoint::advance(
            self,
            power_w_m3,
            forces_n_m3,
            limits,
            buoyancy,
            gate,
            resolve,
            move |candidate, input, molecular, force| {
                let pending = candidate.prepare_thermal_lagrangian_sgs(
                    config,
                    sample_time_lu,
                    Some(molecular),
                    workspace,
                    gate,
                )?;
                let heat =
                    source_step_resolved_mapped::<ThermalLbm3, 3, 7, _, MaterialTransportError>(
                        candidate,
                        &input.heat,
                        &input.force,
                        input.limits,
                        pending.relaxation(),
                        force,
                        |endpoint, raw| {
                            admit_endpoint(endpoint, config.maximum_mach, gate)?;
                            input.balance(raw).map_err(MaterialTransportError::from)
                        },
                    )?;
                Ok((heat, pending))
            },
        )?;
        let sgs = accepted.sgs.report();
        accepted.sgs.commit();
        Ok(MaterialSgsEndpointStep {
            material: accepted.material,
            sgs,
        })
    }
}
