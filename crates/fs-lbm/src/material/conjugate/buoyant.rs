//! Steady and transient natural/mixed convection from one local material response.
//!
//! The existing Boussinesq projection owns reference mass and heat capacity;
//! the existing variable-viscosity SIMPLEC and conductivity-field energy solvers
//! own the discrete equations. This driver only closes their temperature
//! feedback. It never replaces density response with an ideal-gas beta.

mod unsteady;
pub use unsteady::MaterialBuoyancyMarch;

use super::{
    EnergyConfig, EnergySolution, FvFlow, MaterialChtFluid, SimpleConfig, SolidMaterial,
    ThermalSetup, VoxelDomain,
};
use crate::conjugate::{
    ChtError, ThermalFace, Turbulence, Voxel, simple_flow_with_viscosity_and_force,
};
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
    state_point::MaterialStatePointError,
};
use super::solid::{SolidConductivityLaw, SolidConductivitySample};

/// Source-caloric h(T) transport on the same steady material fixed point.
pub mod enthalpy;
use enthalpy::{CaloricCoupling, MaterialEnthalpyReport, MaterialEnthalpyState};
use crate::conjugate::AdvectiveEnthalpyField;

/// Physical approximation limits and bounded steady fixed-point iteration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialBuoyancyConfig {
    /// Shared gravity and maximum departures from reference rho and Cp.
    pub boussinesq: BoussinesqConfig<3>,
    /// Maximum complete momentum/energy solves.
    pub maximum_passes: usize,
    /// Absolute unrelaxed fixed-point temperature defect [K].
    pub temperature_tolerance_k: f64,
    /// Relative tolerance on fluid viscosity and all fluid/solid conductivities.
    pub transport_relative_tolerance: f64,
    /// Absolute tolerance on each buoyancy-force component [N/m3].
    pub force_tolerance_n_m3: f64,
    /// Relaxation of the next trial temperature, in (0, 1]. The unrelaxed
    /// defect is checked, so a small value cannot manufacture convergence.
    pub temperature_relaxation: f64,
}

/// One original source evaluation at one explicitly retained temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialBuoyancySample<F> {
    /// Original domain cell, not a fluid-only index.
    pub cell: usize,
    /// Absolute query temperature [K].
    pub temperature_k: f64,
    /// Original source and the SAME admitted reference-storage projection.
    pub material: BoussinesqMaterial<F>,
    /// Actual `(rho(T)-rho0) * gravity` [N/m3].
    pub force_n_m3: [f64; 3],
}

/// Measured nonlinear defects, independent of the inner solver residuals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialBuoyancyReport {
    /// Complete, converged momentum/energy pairs, excluding an optional opening seed.
    pub passes: usize,
    /// Maximum unrelaxed output-minus-trial temperature magnitude [K].
    pub temperature_defect_k: f64,
    /// Maximum relative change of fluid viscosity or any fluid/solid conductivity.
    pub transport_relative_defect: f64,
    /// Maximum componentwise output-minus-trial buoyancy magnitude [N/m3].
    pub force_defect_n_m3: f64,
}

/// Existing flow and energy results with the material epochs that produced them.
pub struct MaterialBuoyancySolution<F> {
    /// Fixed inertial and sensible-heat storage, shared by both solvers.
    pub reference: FluidTransportState,
    /// Original momentum owner, conservative fluxes and convergence report.
    pub flow: FvFlow,
    /// Original thermal owner on exactly `flow.field`.
    pub energy: EnergySolution,
    /// Actual coefficients used in both returned discrete solves.
    pub applied: Vec<MaterialBuoyancySample<F>>,
    /// Fresh source support at the returned temperatures; differences from
    /// `applied` passed the explicit nonlinear tolerances in `report`.
    pub endpoint: Vec<MaterialBuoyancySample<F>>,
    /// Original solid sources actually used in the heat solve. Empty for the
    /// explicitly fixed-coefficient entry point.
    pub solid_applied: Vec<SolidConductivitySample>,
    /// Those same solid sources at the returned temperature field. Their
    /// conductivity changes pass the same transport convergence gate.
    pub solid_endpoint: Vec<SolidConductivitySample>,
    /// Actual defects admitted at the returned fixed point.
    pub report: MaterialBuoyancyReport,
}

/// Source identity/domain failures remain distinguishable from solver failure.
#[derive(Debug)]
pub enum MaterialBuoyancyError {
    /// Conflicting physical inputs or inadmissible coupling controls.
    Invalid(&'static str),
    /// Source query or Boussinesq bound failed at this exact cell/temperature.
    Material {
        cell: usize,
        temperature_k: f64,
        source: MaterialError,
    },
    /// Original solid conductivity source refused at this exact cell/state.
    SolidMaterial {
        cell: usize,
        temperature_k: f64,
        source: MaterialStatePointError,
    },
    /// Original single-phase caloric law refused this trial or endpoint.
    Caloric {
        temperature_k: f64,
        source: MaterialError,
    },
    /// Original momentum, thermal or cancellation refusal.
    Solver(ChtError),
    /// All available passes completed without a consistent material endpoint.
    NotConverged(MaterialBuoyancyReport),
    /// Local nonlinear enthalpy transport did not close within the same budget.
    EnthalpyNotConverged {
        material: MaterialBuoyancyReport,
        enthalpy: MaterialEnthalpyReport,
    },
}
impl core::fmt::Display for MaterialBuoyancyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material natural convection: {self:?}")
    }
}
impl std::error::Error for MaterialBuoyancyError {}
impl From<ChtError> for MaterialBuoyancyError {
    fn from(error: ChtError) -> Self {
        Self::Solver(error)
    }
}

fn poll(gate: &CancelGate) -> Result<(), MaterialBuoyancyError> {
    if gate.is_requested() {
        Err(ChtError::Cancelled.into())
    } else {
        Ok(())
    }
}

impl MaterialBuoyancyConfig {
    fn validate(self) -> Result<(), MaterialBuoyancyError> {
        self.boussinesq.validate().map_err(|_| {
            MaterialBuoyancyError::Invalid(
                "gravity and reference density/heat-capacity departure bounds must be admitted",
            )
        })?;
        if self.maximum_passes == 0
            || [
                self.temperature_tolerance_k,
                self.transport_relative_tolerance,
                self.force_tolerance_n_m3,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || self.transport_relative_tolerance >= 1.0
            || !self.temperature_relaxation.is_finite()
            || self.temperature_relaxation <= 0.0
            || self.temperature_relaxation > 1.0
        {
            return Err(MaterialBuoyancyError::Invalid(
                "positive pass budget, finite nonnegative tolerances, relative tolerance below one and relaxation in (0,1] required",
            ));
        }
        Ok(())
    }
}

fn samples<F: FluidTransport>(
    domain: &VoxelDomain,
    reference: &FluidTransportState,
    temperature: &[f64],
    config: BoussinesqConfig<3>,
    gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<Vec<MaterialBuoyancySample<F>>, MaterialBuoyancyError> {
    let mut output = Vec::with_capacity(domain.fluid_count());
    for (cell, &temperature_k) in temperature.iter().enumerate() {
        poll(gate)?;
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(MaterialBuoyancyError::Invalid(
                "every trial temperature must be positive and finite",
            ));
        }
        if !domain.is_fluid(cell) {
            continue;
        }
        let error = |source| MaterialBuoyancyError::Material {
            cell,
            temperature_k,
            source,
        };
        let source = resolve(cell, temperature_k, reference).map_err(error)?;
        let material = BoussinesqMaterial::try_new(source, reference, config).map_err(error)?;
        let departure = material.resolved_transport().density_kg_m3() - reference.density_kg_m3();
        let force_n_m3 = config.gravity_m_s2.map(|g| departure * g);
        if force_n_m3.iter().any(|v| !v.is_finite()) {
            return Err(MaterialBuoyancyError::Invalid(
                "material buoyancy force is unrepresentable",
            ));
        }
        output.push(MaterialBuoyancySample {
            cell,
            temperature_k,
            material,
            force_n_m3,
        });
    }
    Ok(output)
}

fn solid_samples(
    domain: &VoxelDomain, temperature: &[f64],
    sources: Option<&[SolidConductivityLaw<'_>]>, gate: &CancelGate,
) -> Result<Vec<SolidConductivitySample>, MaterialBuoyancyError> {
    let Some(sources) = sources else { return Ok(Vec::new()); };
    let mut samples = Vec::new();
    for (cell, &temperature_k) in temperature.iter().enumerate() {
        poll(gate)?;
        let Voxel::Solid(index) = domain.voxel_at(cell) else { continue; };
        let source = sources.get(usize::from(index)).ok_or(MaterialBuoyancyError::Invalid(
            "solid voxel refers to an absent material source"))?;
        samples.push(source.at_temperature(cell, index, temperature_k)
            .map_err(|source| MaterialBuoyancyError::SolidMaterial { cell, temperature_k, source })?);
    }
    Ok(samples)
}

fn admit_open_boundary_materials<F: FluidTransport>(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    reference: &FluidTransportState,
    config: BoussinesqConfig<3>,
    gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<(), MaterialBuoyancyError> {
    for (face, rule) in setup.faces.iter().enumerate() {
        let temperature_k = match *rule {
            ThermalFace::Inflow { temperature }
                | ThermalFace::Outflow { backflow_temperature: temperature } => temperature,
            _ => continue,
        };
        for cell in 0..domain.cell_count() {
            poll(gate)?;
            if !domain.is_fluid(cell) || domain.neighbor(cell, face).is_some() { continue; }
            // The incident cell identifies the SAME caller-owned fluid/source
            // context. A supported downstream temperature cannot admit an
            // unsupported inlet, including potential pressure-opening backflow.
            let error = |source| MaterialBuoyancyError::Material { cell, temperature_k, source };
            let source = resolve(cell, temperature_k, reference).map_err(error)?;
            BoussinesqMaterial::try_new(source, reference, config).map_err(error)?;
        }
    }
    Ok(())
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// Close material feedback through the original steady flow/heat solvers.
    ///
    /// `resolve` evaluates the SAME chosen fluid at each absolute temperature;
    /// pressure, composition and other conditions remain the caller's explicit
    /// source context. Each evaluation supplies rho for buoyancy and mu/k for
    /// both operators. Negative expansion and nonmonotone density need no special
    /// solver code. The inherited reference rho/Cp stay in inertia, continuity
    /// and heat-capacity advection. Their neglected departures are bounded by
    /// `config.boussinesq` at every trial and returned endpoint.
    ///
    /// `initial_temperature_k` provides one positive value per original cell
    /// and supplies the first momentum and thermal solves. Only if that flow is
    /// exactly zero and the thermal owner refuses an unanchored temperature
    /// problem do pressure openings receive one conductive seed with their
    /// backflow temperatures held. Useful warm starts and pressure-driven flows
    /// therefore retain their actual transport and source-domain admission.
    /// All complete passes use the original open thermal boundaries.
    /// Each pass solves momentum to `flow_config` convergence, then heat on
    /// exactly those conservative fluxes. Acceptance also checks the UNRELAXED
    /// temperature defect, local mu/k change, and absolute force change after
    /// re-evaluating the proposed endpoint. A discontinuous source cannot hide
    /// behind a small temperature update. All returned source epochs are retained.
    ///
    /// Solids use their explicit fixed grid-aligned conductivity. Stateful solid
    /// response, turbulent closure and independently authored eddy conductivity
    /// refuse instead of competing with this material owner. Ordinary thermal
    /// contacts, sources, fixed boundaries and sinks retain their existing owner.
    /// No transient momentum, phase change, variable-density enthalpy advection,
    /// gravitational storage or viscous heating is added. This is a bounded
    /// numerical steady estimate, not experimental validation of a fluid model.
    ///
    /// Source callbacks must be repeatable queries, never commit physical history.
    /// Any material, convergence or cancellation failure returns no partial
    /// solution and leaves every caller-owned input untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_natural_convection<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: MaterialBuoyancyConfig,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialBuoyancySolution<S>, MaterialBuoyancyError> {
        self.solve_natural_convection_inner(domain, solids, None, flow_config, thermal,
            energy_config, initial_temperature_k, config, None, gate, resolve)
            .map(|(solution, _)| solution)
    }

    /// Close BOTH solid and fluid material response on the same temperature field.
    ///
    /// Every original solid assignment names one source-bound conductivity law.
    /// Its complete source bundle is refreshed at the actual trial, seed and
    /// returned endpoint temperatures, with original claims pinned. Scalar or
    /// explicitly grid-aligned orthotropic conductivity flows through the same
    /// conservative half-cell/interface heat assembly. The original flow solver
    /// then sees the resulting fluid temperature, density and viscosity.
    ///
    /// Acceptance includes every solid principal conductivity in the existing
    /// transport-defect gate; small temperature changes cannot conceal a large
    /// material change. Source refusals, cancellation and nonlinear exhaustion
    /// return no partial solution and mutate neither inputs nor material laws.
    /// Steady fixed geometry and the original Boussinesq/laminar approximation
    /// remain; this supplies no solid storage, deformation, melting or off-axis
    /// conductivity discretization. No second thermal solver is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_natural_convection_with_material_solids<S: FluidTransport>(
        &self, domain: &VoxelDomain, sources: &[SolidConductivityLaw<'_>],
        flow_config: &SimpleConfig, thermal: &ThermalSetup, energy_config: &EnergyConfig,
        initial_temperature_k: &[f64], config: MaterialBuoyancyConfig, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialBuoyancySolution<S>, MaterialBuoyancyError> {
        poll(gate)?;
        let solids: Vec<_> = sources.iter().map(|source| {
            let k = source.reference_conductivity_w_m_k();
            SolidMaterial::new("source-bound solid", k[0]).with_orthotropic(k)
        }).collect();
        self.solve_natural_convection_inner(domain, &solids, Some(sources), flow_config,
            thermal, energy_config, initial_temperature_k, config, None, gate, resolve)
            .map(|(solution, _)| solution)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn solve_natural_convection_inner<S: FluidTransport>(
        &self, domain: &VoxelDomain, solids: &[SolidMaterial],
        solid_sources: Option<&[SolidConductivityLaw<'_>]>,
        flow_config: &SimpleConfig, thermal: &ThermalSetup, energy_config: &EnergyConfig,
        initial_temperature_k: &[f64], config: MaterialBuoyancyConfig,
        caloric: Option<CaloricCoupling<'_>>, gate: &CancelGate,
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<(MaterialBuoyancySolution<S>, Option<MaterialEnthalpyState>), MaterialBuoyancyError> {
        poll(gate)?;
        config.validate()?;
        if let Some(caloric) = caloric { caloric.validate(energy_config)?; }
        if initial_temperature_k.len() != domain.cell_count() || domain.fluid_count() == 0 {
            return Err(MaterialBuoyancyError::Invalid(
                "one initial temperature per cell and at least one fluid cell required",
            ));
        }
        if flow_config.turbulence != Turbulence::Laminar
            || !thermal.eddy_conductivity_w_m_k.is_empty()
            || solids.iter().any(|s| !s.conductivity_table.is_empty())
        {
            return Err(MaterialBuoyancyError::Invalid(
                "material steady feedback requires laminar flow and no independently authored eddy or solid conductivity tables",
            ));
        }
        let mut conductivity =
            vec![[self.transport.thermal_conductivity_w_m_k(); 3]; domain.cell_count()];
        for (cell, values) in conductivity.iter_mut().enumerate() {
            poll(gate)?;
            if let Voxel::Solid(index) = domain.voxel_at(cell) {
                *values = solids
                    .get(usize::from(index))
                    .ok_or(MaterialBuoyancyError::Invalid(
                        "solid voxel refers to an absent material",
                    ))?
                    .axis_conductivity();
            }
        }
        let mut temperature = initial_temperature_k.to_vec();
        let mut applied = samples(
            domain,
            &self.transport,
            &temperature,
            config.boussinesq,
            gate,
            &mut resolve,
        )?;
        let mut solid_applied = solid_samples(domain, &temperature, solid_sources, gate)?;
        if caloric.is_some() {
            admit_open_boundary_materials(domain, thermal, &self.transport,
                config.boussinesq, gate, &mut resolve)?;
        }
        // Boundary data belong to the same immutable caloric source and common
        // datum. They stay at their actual inlet/backflow temperatures.
        let boundary_h = caloric.map(|law| law.boundaries(thermal, &self.transport, config, gate))
            .transpose()?.unwrap_or([None; 6]);
        let mut caloric_applied = caloric.map(|law|
            law.field(domain, &temperature, &self.transport, config, gate)).transpose()?;
        let mut viscosity = vec![self.transport.dynamic_viscosity_pa_s(); domain.cell_count()];
        let mut force = vec![[0.0; 3]; domain.cell_count()];
        let mut seeded = false;
        for pass in 1..=config.maximum_passes {
            poll(gate)?;
            let (flow, energy) = loop {
                for sample in &applied {
                    let transport = sample.material.reference_storage_transport();
                    viscosity[sample.cell] = transport.dynamic_viscosity_pa_s();
                    conductivity[sample.cell] = [transport.thermal_conductivity_w_m_k(); 3];
                    force[sample.cell] = sample.force_n_m3;
                }
                for sample in &solid_applied {
                    conductivity[sample.cell] = sample.conductivity_w_m_k;
                }
                let flow = simple_flow_with_viscosity_and_force(
                    domain,
                    &self.properties,
                    &viscosity,
                    &force,
                    flow_config,
                    gate,
                )?;
                let solved = if let Some(cells) = &caloric_applied {
                    energy_config.solve_with_enthalpy(
                        domain, &self.properties, solids, &flow.field, thermal, &conductivity,
                        AdvectiveEnthalpyField { cells, boundary_j_kg: boundary_h }, gate,
                    )
                } else {
                    energy_config.solve_with_conductivity(
                        domain, &self.properties, solids, &flow.field, thermal, &conductivity, gate,
                    )
                };
                let error = match solved {
                    Ok(energy) => break (flow, energy),
                    Err(error) => error,
                };
                // Keep anchor admission with the thermal owner. Do not replace a
                // useful initial field, or discard fan/pressure-driven advection
                // merely because its material buoyancy force happens to be zero.
                let no_anchor = matches!(&error,
                    ChtError::InvalidInput { field: "thermal.faces", reason }
                        if reason.starts_with("no temperature anchor:")
                            || (reason.starts_with("cell ")
                                && reason.contains(" has no temperature anchor (")));
                if pass != 1
                    || seeded
                    || flow.field.max_flux() != 0.0
                    || !no_anchor
                    || !thermal.faces.iter().any(|face| matches!(face, ThermalFace::Outflow { .. }))
                {
                    return Err(error.into());
                }
                let mut seed = thermal.clone();
                for face in &mut seed.faces {
                    if let ThermalFace::Outflow { backflow_temperature } = *face {
                        *face = ThermalFace::Temperature(backflow_temperature);
                    }
                }
                temperature = energy_config
                    .solve_with_conductivity(
                        domain,
                        &self.properties,
                        solids,
                        &flow.field,
                        &seed,
                        &conductivity,
                        gate,
                    )?
                    .temperature;
                applied = samples(
                    domain,
                    &self.transport,
                    &temperature,
                    config.boussinesq,
                    gate,
                    &mut resolve,
                )?;
                solid_applied = solid_samples(domain, &temperature, solid_sources, gate)?;
                caloric_applied = caloric.map(|law|
                    law.field(domain, &temperature, &self.transport, config, gate)).transpose()?;
                seeded = true;
            };
            let endpoint = samples(
                domain,
                &self.transport,
                &energy.temperature,
                config.boussinesq,
                gate,
                &mut resolve,
            )?;
            let solid_endpoint = solid_samples(domain, &energy.temperature, solid_sources, gate)?;
            let caloric_endpoint = caloric.map(|law|
                law.field(domain, &energy.temperature, &self.transport, config, gate)).transpose()?;
            let caloric_report = if let (Some(old), Some(next)) = (&caloric_applied, &caloric_endpoint) {
                Some(enthalpy::defects(domain, &flow.field, self.transport.density_kg_m3(),
                    &energy.temperature, old, next, gate)?)
            } else { None };
            let mut report = MaterialBuoyancyReport {
                passes: pass,
                temperature_defect_k: 0.0,
                transport_relative_defect: 0.0,
                force_defect_n_m3: 0.0,
            };
            for (old, next) in temperature.iter().zip(&energy.temperature) {
                report.temperature_defect_k = report.temperature_defect_k.max((next - old).abs());
            }
            for (old, next) in applied.iter().zip(&endpoint) {
                poll(gate)?;
                let a = old.material.reference_storage_transport();
                let b = next.material.reference_storage_transport();
                for (x, y) in [
                    (a.dynamic_viscosity_pa_s(), b.dynamic_viscosity_pa_s()),
                    (
                        a.thermal_conductivity_w_m_k(),
                        b.thermal_conductivity_w_m_k(),
                    ),
                ] {
                    report.transport_relative_defect = report
                        .transport_relative_defect
                        .max((x - y).abs() / x.max(y));
                }
                for (a, b) in old.force_n_m3.iter().zip(next.force_n_m3) {
                    report.force_defect_n_m3 = report.force_defect_n_m3.max((b - a).abs());
                }
            }
            for (old, next) in solid_applied.iter().zip(&solid_endpoint) {
                poll(gate)?;
                for (x, y) in old.conductivity_w_m_k.iter().zip(next.conductivity_w_m_k) {
                    report.transport_relative_defect = report.transport_relative_defect
                        .max((x-y).abs()/x.max(y));
                }
            }
            poll(gate)?;
            if report.temperature_defect_k <= config.temperature_tolerance_k
                && report.transport_relative_defect <= config.transport_relative_tolerance
                && report.force_defect_n_m3 <= config.force_tolerance_n_m3
                && caloric_report.zip(caloric).is_none_or(|(report, law)| {
                    report.maximum_face_flux_defect_w <= law.config.maximum_flux_defect_w
                        && report.maximum_cell_flux_defect_w <= law.config.maximum_flux_defect_w
                })
            {
                let enthalpy = caloric_report.map(|report| MaterialEnthalpyState {
                    applied: caloric_applied.expect("caloric report has an applied field"),
                    endpoint: caloric_endpoint.expect("caloric report has an endpoint field"),
                    boundary_j_kg: boundary_h,
                    report,
                });
                return Ok((MaterialBuoyancySolution {
                    reference: self.transport,
                    flow,
                    energy,
                    applied,
                    endpoint,
                    solid_applied,
                    solid_endpoint,
                    report,
                }, enthalpy));
            }
            if pass == config.maximum_passes {
                return Err(match caloric_report {
                    Some(enthalpy) => MaterialBuoyancyError::EnthalpyNotConverged {
                        material: report, enthalpy,
                    },
                    None => MaterialBuoyancyError::NotConverged(report),
                });
            }
            if config.temperature_relaxation == 1.0 {
                temperature = energy.temperature;
                applied = endpoint;
                solid_applied = solid_endpoint;
                caloric_applied = caloric_endpoint;
            } else {
                for (old, next) in temperature.iter_mut().zip(energy.temperature) {
                    *old = (1.0 - config.temperature_relaxation) * *old
                        + config.temperature_relaxation * next;
                }
                applied = samples(
                    domain,
                    &self.transport,
                    &temperature,
                    config.boussinesq,
                    gate,
                    &mut resolve,
                )?;
                solid_applied = solid_samples(domain, &temperature, solid_sources, gate)?;
                caloric_applied = caloric.map(|law|
                    law.field(domain, &temperature, &self.transport, config, gate)).transpose()?;
            }
        }
        unreachable!("positive bounded pass count either accepts or refuses its final pass")
    }
}
