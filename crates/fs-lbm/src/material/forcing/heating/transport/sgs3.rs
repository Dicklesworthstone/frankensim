//! Source-resolved molecular transport with native thermal SGS.
//!
//! Initial and endpoint material queries retain their concrete owners. Only
//! step-local relaxation is augmented; reference density, sensible storage and
//! the accepted Guo force frame are never rebound. Both closures use the
//! original source/wall heat transaction and material publication checkpoint.

use super::super::Inputs;
use super::boussinesq::{self, BoussinesqConfig, BoussinesqMaterial};
use super::{
    FluidTransport, FluidTransportState, HeatBalance, HeatSourceLimits, MaterialError,
    MaterialSlab, MaterialTransportError, MaterialTransportStep, RelaxationField,
    ThermalLbm3, TransportSample, advance_with_sampler_sources_and_candidate,
    candidate_with_transport, sample,
};
use crate::thermal::sources::sgs3::channel::{
    ThermalChannelSgsConfig, ThermalChannelSgsReport, ThermalChannelSgsWorkspace3,
    step_channel_sgs_mapped,
};
use crate::thermal::sources::sgs3::{
    ThermalSgsConfig, ThermalSgsError, ThermalSgsReport, ThermalSgsWorkspace3, step_sgs_mapped,
};
use fs_exec::CancelGate;

/// Original molecular material records and SI heat balance, together with the
/// closure actually applied to this step. Material records never contain SGS.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialSgsTransportStep<F, R> {
    pub material: MaterialTransportStep<F>,
    pub sgs: R,
}

#[allow(clippy::too_many_arguments)]
fn advance<F: FluidTransport, R: Default>(
    slab: &mut MaterialSlab<ThermalLbm3>,
    power: &[f64],
    forces: &[[f64; 3]],
    limits: HeatSourceLimits,
    gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    replace_buoyancy: bool,
    force: impl FnOnce(
        &MaterialSlab<ThermalLbm3>, &mut Inputs<3>, &[TransportSample<F>],
    ) -> Result<(), MaterialTransportError>,
    step: impl FnOnce(
        &mut ThermalLbm3, &Inputs<3>, RelaxationField<'_>,
    ) -> Result<(HeatBalance, R), ThermalSgsError>,
) -> Result<MaterialSgsTransportStep<F, R>, MaterialTransportError> {
    let depth = slab.force_units.cell_size_m;
    let background = slab.lattice.gbeta;
    let mut report = R::default();
    let (material, ()) = advance_with_sampler_sources_and_candidate::<ThermalLbm3, F, (), 3, 7>(
        slab, power, forces, depth, limits, gate,
        |mapping, state, cells, plane| {
            sample::<ThermalLbm3, F, 3, 7>(
                mapping, state, cells, plane, limits.temperature_range_k, gate, &mut resolve,
            )
        },
        |mapping, mut input, applied| {
            force(mapping, &mut input, applied)?;
            Ok((input, ()))
        },
        |mapping, input, applied| {
            let (mut trial, heat) = candidate_with_transport::<ThermalLbm3, F, 3, 7>(
                mapping, input, applied, gate,
                |candidate| {
                    if replace_buoyancy {
                        candidate.gbeta = 0.0;
                    }
                },
                |candidate, input, molecular| {
                    // The snapshot still carries the accepted fext. SGS reads
                    // that force-centered velocity before the source kernel
                    // applies this step's newly prepared external/buoyant load.
                    let (heat, applied_report) = step(candidate, input, molecular)?;
                    report = applied_report;
                    Ok(heat)
                },
            )?;
            trial.gbeta = background;
            Ok((trial, heat))
        },
    )?;
    Ok(MaterialSgsTransportStep { material, sgs: report })
}

#[allow(clippy::too_many_arguments)]
fn advance_boussinesq<F: FluidTransport, R: Default>(
    slab: &mut MaterialSlab<ThermalLbm3>,
    power: &[f64],
    forces: &[[f64; 3]],
    limits: HeatSourceLimits,
    buoyancy: BoussinesqConfig<3>,
    gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    step: impl FnOnce(
        &mut ThermalLbm3, &Inputs<3>, RelaxationField<'_>,
    ) -> Result<(HeatBalance, R), ThermalSgsError>,
) -> Result<MaterialSgsTransportStep<BoussinesqMaterial<F>, R>, MaterialTransportError> {
    super::poll(gate)?;
    buoyancy.validate()?;
    advance(
        slab, power, forces, limits, gate,
        |cell, temperature, retained| {
            BoussinesqMaterial::try_new(
                resolve(cell, temperature, retained)?, retained, buoyancy,
            )
        },
        true,
        |mapping, input, applied| {
            boussinesq::project_force::<ThermalLbm3, BoussinesqMaterial<F>, 3, 7>(
                mapping, input, forces, applied, buoyancy, gate,
                |material| material.resolved_transport().density_kg_m3(),
            )
        },
        step,
    )
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance source-resolved molecular mu(T), k(T) and explicit local SGS.
    ///
    /// The original sampler resolves each fluid cell once at the accepted
    /// temperature and once at the candidate endpoint, retaining both concrete
    /// source owners. Density and Cp must equal the retained storage snapshot.
    /// Fresh molecular rates supply the local SGS ratio bound; augmented rates
    /// are never accumulated into the retained material or next step.
    ///
    /// Power [W/m3], force [N/m3], source splitting, stationary halfway plates,
    /// independently measured wall heat and SI balance keep their existing
    /// owners. SGS uses accepted Guo-centered velocity, checks endpoint Mach,
    /// and remains diagnostic rather than adding dissipation heat.
    /// Material/closure refusal or cancellation publishes no candidate.
    /// Workspaces are reusable scratch and may change on refusal.
    ///
    /// Coefficients are frozen for one explicit step; endpoint queries check
    /// source support, not implicit convergence. This is low-Mach sensible heat,
    /// without variable storage, viscous heating or total-energy qualification.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_transport_and_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: ThermalSgsConfig,
        workspace: &mut ThermalSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsTransportStep<F, ThermalSgsReport>, MaterialTransportError> {
        advance(
            self, power_w_m3, forces_n_m3, limits, gate, resolve, false,
            |_, _, _| Ok(()),
            |candidate, input, molecular| {
                step_sgs_mapped(
                    candidate, &input.heat, &input.force, config, Some(molecular),
                    workspace, input.limits, gate, |raw| input.balance(raw),
                )
            },
        )
    }

    /// Independently fit plane-channel momentum and heat using the SAME local
    /// molecular source epochs and transaction as the local SGS method.
    /// Molecular mu/rho0 and k/(rho0*Cp0) provide their respective ratio bounds.
    /// A zero momentum fit does not disable modeled heat. The existing channel
    /// topology, scalar wall conditions, byte budget and usage limits apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_transport_and_channel_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        config: ThermalChannelSgsConfig,
        workspace: &mut ThermalChannelSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialSgsTransportStep<F, ThermalChannelSgsReport>, MaterialTransportError> {
        advance(
            self, power_w_m3, forces_n_m3, limits, gate, resolve, false,
            |_, _, _| Ok(()),
            |candidate, input, molecular| {
                step_channel_sgs_mapped(
                    candidate, &input.heat, &input.force, config, Some(molecular),
                    workspace, input.limits, gate, |raw| input.balance(raw),
                )
            },
        )
    }

    /// Source-density buoyancy and local SGS from each SAME material sample.
    /// The original bounded projection uses (rho(T)-rho0)*gravity, replaces
    /// the slab's background beta force for this step, and retains rho0/Cp0
    /// for inertia and sensible storage. Both material epochs must satisfy
    /// the authored density/Cp departure limits. Original material responses
    /// remain available through BoussinesqMaterial; no extra source query is made.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>,
        config: ThermalSgsConfig,
        workspace: &mut ThermalSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        MaterialSgsTransportStep<BoussinesqMaterial<F>, ThermalSgsReport>,
        MaterialTransportError,
    > {
        advance_boussinesq(
            self, power_w_m3, forces_n_m3, limits, buoyancy, gate, resolve,
            |candidate, input, molecular| {
                step_sgs_mapped(
                    candidate, &input.heat, &input.force, config, Some(molecular),
                    workspace, input.limits, gate, |raw| input.balance(raw),
                )
            },
        )
    }

    /// Bounded source-density buoyancy with independently fitted channel stress
    /// and heat. Retains the material, force frame, storage and publication
    /// rules of the local Boussinesq SGS method, including endpoint admission.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_transport_and_channel_sgs<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        buoyancy: BoussinesqConfig<3>,
        config: ThermalChannelSgsConfig,
        workspace: &mut ThermalChannelSgsWorkspace3,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        MaterialSgsTransportStep<BoussinesqMaterial<F>, ThermalChannelSgsReport>,
        MaterialTransportError,
    > {
        advance_boussinesq(
            self, power_w_m3, forces_n_m3, limits, buoyancy, gate, resolve,
            |candidate, input, molecular| {
                step_channel_sgs_mapped(
                    candidate, &input.heat, &input.force, config, Some(molecular),
                    workspace, input.limits, gate, |raw| input.balance(raw),
                )
            },
        )
    }
}

#[cfg(test)]
mod tests;
