//! Generalized-Newtonian stress power on the existing thermal/flow step.
//!
//! The same initial material evaluation supplies collision viscosity and
//! nonnegative stress power. In 3-D, the source's shear-rate query excludes
//! dilation, while heat retains `2 mu D_dev:D_dev + zeta (div u)^2`.
//! The ordinary collision family ties `zeta = 2 mu / 3`; independent bulk
//! relaxation requires the explicitly sourced physical bulk entry point.
//! The 2-D incompressible law retains `2 mu D:D = mu gamma_dot^2`.
//! No numerical kinetic-energy residual, pressure work, turbulence closure or
//! elastic energy is relabelled as thermal power.

use super::super::super::Inputs;
use super::super::boussinesq::{
    BackgroundBuoyancy, BoussinesqConfig, BoussinesqMaterial, candidate_with_buoyancy,
};
use super::super::{advance_with_sampler_and_sources, advance_with_sampler_sources_and_candidate,
    candidate_with_prepared_input};
use super::{
    DilationRate, FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError, MaterialSlab,
    MaterialTransportError, MaterialTransportStep, RheologyMaterial, RheologyPoint, RheologySolve,
    StressSlab, ThermalLbm, ThermalLbm3, TransportSample, poll, sample, sample_with_strain,
};
use fs_exec::CancelGate;

/// Coherent local shear/bulk material laws on the same thermal transaction.
pub mod bulk;
pub use bulk::BulkViscousMaterial;

fn physical_bulk_tau(mapping: &MaterialSlab<ThermalLbm3>, bulk_viscosity_pa_s: f64) -> Result<f64, &'static str> {
    let units = mapping.force_units;
    crate::material::physical_bulk_tau(
        units.transport.density_kg_m3(), units.velocity_scale, units.cell_size_m,
        bulk_viscosity_pa_s,
    )
}

/// Explicit source partition for an accepted physical viscous-heating step.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousHeating {
    /// Constitutive stress power [W/m3], with zero wall entries. The 3-D APIs
    /// separate deviatoric shear and volumetric power without counting
    /// dilation twice. The ordinary collision's bulk/shear tie is retained;
    /// an independent bulk coefficient requires its explicit physical source.
    /// These are the applied initial-state values, held over one time step.
    pub power_w_m3: Vec<f64>,
    /// Integrated constitutive stress work converted into sensible heat [J].
    /// This is not an independently measured decrement of lattice kinetic energy.
    pub viscous_heat_j: f64,
    /// Separately authored signed volumetric heat integrated over the step [J].
    /// It must exclude the molecular viscous heating owned by this method.
    pub authored_heat_j: f64,
    /// Combined thermal source minus authored heat minus viscous heat [J].
    /// Checked against the same caller-specified thermal balance tolerance.
    pub source_partition_defect_j: f64,
}

/// Original material/thermal receipts plus the explicitly owned viscous heat.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousRheologyStep<F> {
    /// The same transport step and both source epochs as the unheated API.
    /// Its heat source is the sum of authored and physical viscous heat.
    pub transport: MaterialTransportStep<RheologyMaterial<F>>,
    /// Separate physical loss and external-source accounting.
    pub viscous: ViscousHeating,
}

#[allow(clippy::too_many_arguments)]
fn source_inputs<S, F, const D: usize>(
    mapping: &MaterialSlab<S>,
    authored: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    original: Inputs<D>,
    applied: &[TransportSample<RheologyMaterial<F>>],
    bulk_viscosity_pa_s: impl Fn(&RheologyMaterial<F>) -> Option<f64>,
    gate: &CancelGate,
) -> Result<(Inputs<D>, ViscousHeating), MaterialTransportError> {
    let units = mapping.force_units;
    let volume = (units.cell_size_m * units.cell_size_m) * depth_m;
    let mut power_w_m3 = vec![0.0; authored.len()];
    let mut combined = authored.to_vec();
    let mut viscous_heat_j = 0.0;
    for value in applied {
        poll(gate)?;
        let rate = value.material.point().shear_rate_per_s;
        let viscosity = value.transport.dynamic_viscosity_pa_s();
        // This multiplication order avoids overflowing rate^2 before a small
        // viscosity can scale it. Both values come from the cached query.
        let shear_power = (viscosity * rate) * rate;
        let (dilation, bulk_power) = match (
            bulk_viscosity_pa_s(&value.material), value.material.dilation_rate_per_s,
        ) {
            (None, None) => (0.0, 0.0),
            (None, Some(DilationRate::Tied(dilation))) => {
                (dilation, (viscosity * dilation) * dilation)
            }
            (Some(bulk), Some(DilationRate::Independent(dilation))) => {
                (dilation, ((bulk * dilation) * dilation) * 1.5)
            }
            _ => return Err(MaterialTransportError::Invalid(
                "independent bulk relaxation requires explicitly sourced bulk heating",
            )),
        };
        if (rate != 0.0 && shear_power == 0.0) || (dilation != 0.0 && bulk_power == 0.0) {
            return Err(MaterialTransportError::Invalid("physical viscous heating is unrepresentable"));
        }
        let power = shear_power + bulk_power;
        let energy = (power * volume) * units.timestep_s;
        if !power.is_finite()
            || power < 0.0
            || (rate != 0.0 && power == 0.0)
            || !energy.is_finite()
            || (power != 0.0 && energy == 0.0)
        {
            return Err(MaterialTransportError::Invalid(
                "physical viscous heating is unrepresentable",
            ));
        }
        power_w_m3[value.cell] = power;
        // Preserve the existing zero-source path, including signed zero.
        if power != 0.0 {
            combined[value.cell] += power;
        }
        viscous_heat_j += energy;
    }
    if !viscous_heat_j.is_finite() {
        return Err(MaterialTransportError::Invalid(
            "summed physical viscous heat is unrepresentable",
        ));
    }
    let authored_heat_j = original.prescribed_j;
    let input = mapping.heat_inputs(&combined, forces, authored.len(), depth_m, limits)?;
    let source_partition_defect_j = (input.prescribed_j - authored_heat_j) - viscous_heat_j;
    if !source_partition_defect_j.is_finite()
        || source_partition_defect_j.abs() > limits.balance_tolerance_j
    {
        return Err(MaterialTransportError::Invalid(
            "authored and viscous heat partition exceeds the requested joule tolerance",
        ));
    }
    Ok((
        input,
        ViscousHeating {
            power_w_m3,
            viscous_heat_j,
            authored_heat_j,
            source_partition_defect_j,
        },
    ))
}

#[allow(clippy::too_many_arguments)]
fn advance<S: StressSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    solve: RheologySolve,
    gate: &CancelGate,
    mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
    solve.validate()?;
    let (transport, viscous) =
        advance_with_sampler_and_sources::<S, RheologyMaterial<F>, ViscousHeating, D, Q>(
            slab,
            power,
            forces,
            depth_m,
            limits,
            gate,
            |mapping, state, cells, plane| {
                sample::<S, F, D, Q>(
                    mapping,
                    state,
                    cells,
                    plane,
                    limits.temperature_range_k,
                    solve,
                    gate,
                    &mut resolve,
                )
            },
            |mapping, original, applied| {
                source_inputs(
                    mapping, power, forces, depth_m, limits, original, applied, |_| None, gate,
                )
            },
        )?;
    Ok(ViscousRheologyStep { transport, viscous })
}

#[allow(clippy::too_many_arguments)]
fn advance_boussinesq<
    S: StressSlab<D, Q> + BackgroundBuoyancy,
    F: FluidTransport,
    const D: usize,
    const Q: usize,
>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    solve: RheologySolve,
    config: BoussinesqConfig<D>,
    gate: &CancelGate,
    mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<ViscousRheologyStep<BoussinesqMaterial<F>>, MaterialTransportError> {
    poll(gate)?;
    solve.validate()?;
    config.validate()?;
    let mut bounded = |point, reference: &FluidTransportState| {
        BoussinesqMaterial::try_new(resolve(point, reference)?, reference, config)
    };
    let (transport, viscous) = advance_with_sampler_sources_and_candidate::<
        S,
        RheologyMaterial<BoussinesqMaterial<F>>,
        ViscousHeating,
        D,
        Q,
    >(
        slab,
        power,
        forces,
        depth_m,
        limits,
        gate,
        |mapping, state, cells, plane| {
            sample::<S, BoussinesqMaterial<F>, D, Q>(
                mapping,
                state,
                cells,
                plane,
                limits.temperature_range_k,
                solve,
                gate,
                &mut bounded,
            )
        },
        |mapping, original, applied| {
            source_inputs(
                mapping, power, forces, depth_m, limits, original, applied, |_| None, gate,
            )
        },
        |mapping, input, applied| {
            candidate_with_buoyancy::<S, RheologyMaterial<BoussinesqMaterial<F>>, D, Q>(
                mapping,
                input,
                forces,
                applied,
                config,
                gate,
                |material| material.source().resolved_transport().density_kg_m3(),
            )
        },
    )?;
    Ok(ViscousRheologyStep { transport, viscous })
}

impl MaterialSlab<ThermalLbm> {
    /// Advance with material viscosity AND its physical viscous heat.
    ///
    /// `power_w_m3` supplies external heat only: this method adds the molecular
    /// `mu * gamma_dot^2` loss once, using the exact initial viscosity/shear
    /// samples that drive momentum. The original two source half-kicks advance
    /// sensible heat, and endpoint material queries see the resulting actual
    /// temperature. Repeated calls close explicit temperature/viscosity feedback.
    /// Source, conversion, temperature, heat-balance, endpoint and cancellation
    /// refusals all preserve the accepted slab.
    /// Independent central-moment bulk relaxation is refused: this shear-only
    /// source has no separately authored coefficient for volumetric heat.
    ///
    /// The incompressible generalized-Newtonian, fixed-density/Cp, low-Mach
    /// assumptions of `step_with_material_rheology` still apply. Loss power is
    /// frozen at the initial state: this is not an implicit thermal solve or a
    /// second-order accuracy claim. The thermal audit does not prove transient
    /// total-energy conservation; numerical momentum dissipation and pressure
    /// work are neither measured nor converted into heat here.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }

    /// Couple material buoyancy, shear-dependent viscosity and viscous heating.
    ///
    /// ONE source evaluation at each actual `(temperature, shear rate)` query
    /// supplies density, viscosity and conductivity. The original bounded
    /// constitutive solve resolves shear and viscosity together. Its converged
    /// sample supplies both `mu * rate^2` heat and `(rho-rho0) * gravity` force,
    /// replacing the slab's independent beta load for this step. Authored
    /// external forces remain additive; authored power must exclude this
    /// molecular viscous loss. Original source receipts, unprojected transport
    /// and actual constitutive coordinates survive in the returned samples.
    ///
    /// Reference density/Cp retain mass and sensible storage. Both initial and
    /// endpoint queries must satisfy the explicit Boussinesq departure bounds
    /// and source domains. All owners use the same tentative step and publish
    /// together after endpoint admission; any refusal or cancellation preserves
    /// the accepted slab. Call this method every step to retain feedback.
    ///
    /// This is explicit, low-Mach, generalized-Newtonian Boussinesq evolution;
    /// local constitutive convergence is not implicit time integration. The
    /// sensible-heat account does not include gravitational potential, resolved
    /// kinetic storage, compressible work or numerical momentum dissipation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        config: BoussinesqConfig<2>,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<BoussinesqMaterial<F>>, MaterialTransportError> {
        advance_boussinesq::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            config,
            gate,
            resolve,
        )
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Couple local temperature/shear-dependent material transport to explicit
    /// physical shear AND bulk viscous heat on the original D3Q19/D3Q7 slab.
    ///
    /// `bulk_viscosity_pa_s` is a declared physical material coefficient,
    /// spatially uniform and held for this step. It supplies both
    /// `tau_bulk = 1/2 + (9/2)(zeta/rho_ref) dt/dx^2` and the volumetric
    /// power `zeta (div u)^2`. The material callback still supplies the local
    /// `mu(T, sqrt(2 dev(D):dev(D)))`, conductivity and retained source at both epochs.
    /// Deviatoric power is `2 mu D_dev:D_dev`, so dilation is counted once.
    /// Authored power must exclude both constitutive losses.
    ///
    /// The independent central-moment selector, source half-kicks, original
    /// material receipts and endpoint temperature/energy checks publish
    /// together. Every refusal leaves the prior populations, force frame and
    /// collision selector intact. Accepted steps retain the selected bulk and
    /// higher-order rates; subsequent calls can explicitly change physical
    /// zeta. The legacy heated API refuses that selector without this binding.
    ///
    /// Zero or unrepresentably small zeta is refused: this collision family
    /// requires `0 < bulk_rate < 2`, with no artificial positive viscosity floor.
    /// `higher_order_rate` also lies in `(0, 2)` and never supplies heat.
    /// This is explicit low-Mach sensible-heat transport, not a compressible
    /// total-energy solver, variable-density storage or implicit heating solve.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rheology_and_bulk_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        bulk_viscosity_pa_s: f64,
        higher_order_rate: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
        poll(gate)?;
        solve.validate()?;
        let units = self.force_units;
        let bulk_tau = physical_bulk_tau(self, bulk_viscosity_pa_s).map_err(MaterialTransportError::Invalid)?;
        let bulk_rate = 1.0 / bulk_tau;
        let model = crate::d3q19::CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0, bulk_rate, higher_order_rate,
        };
        if model.validate().is_err() {
            return Err(MaterialTransportError::Invalid("physical bulk viscosity or relaxation is unrepresentable"));
        }
        let depth_m = units.cell_size_m;
        let (transport, viscous) = advance_with_sampler_sources_and_candidate::<
            ThermalLbm3, RheologyMaterial<F>, ViscousHeating, 3, 7,
        >(
            self, power_w_m3, forces_n_m3, depth_m, limits, gate,
            |mapping, state, cells, plane| {
                sample_with_strain::<ThermalLbm3, F, 3, 7>(
                    mapping, state, cells, plane, limits.temperature_range_k,
                    solve, gate, &mut resolve, true,
                    |cell| {
                        let (rho, velocity) = state.grid.moments3(cell);
                        crate::d3q19::rheology_stress_components3(
                            &state.grid.f[cell], rho, velocity, state.grid.fext[cell], model,
                        )
                    },
                )
            },
            |mapping, original, applied| source_inputs(
                mapping, power_w_m3, forces_n_m3, depth_m, limits, original,
                applied, |_| Some(bulk_viscosity_pa_s), gate,
            ),
            |mapping, input, applied| candidate_with_prepared_input::<ThermalLbm3, _, 3, 7>(
                mapping, input, applied, gate,
                |trial| {
                    trial.grid.trt_magic = None;
                    trial.grid.central_moment_bulk = Some(crate::d3q19::BulkRelaxation3::Uniform([bulk_rate, higher_order_rate]));
                },
            ),
        )?;
        Ok(ViscousRheologyStep { transport, viscous })
    }

    /// 3D counterpart with a deviatoric material shear query and the collision's
    /// existing `zeta = 2 mu / 3` volume heating. Independent bulk relaxation
    /// requires the explicitly sourced bulk method. Cell volume is the retained
    /// `dx^3`; the same source feedback, atomicity and low-Mach limits apply.
    pub fn step_with_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<F>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }

    /// Three-dimensional counterpart of the combined material-buoyancy,
    /// rheology and viscous-heat step, with the retained `dx^3` cell volume.
    /// The same source samples, balances and atomic publication are used.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boussinesq_material_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        solve: RheologySolve,
        config: BoussinesqConfig<3>,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<BoussinesqMaterial<F>>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance_boussinesq::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            config,
            gate,
            resolve,
        )
    }
}
