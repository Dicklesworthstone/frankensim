//! Steady Boussinesq natural convection with the finite-volume SIMPLEC flow
//! ([`super::simple_flow`]) and the conservative conjugate energy equation
//! ([`super::solve_energy`]), including OPEN boundaries: `FvBoundary::Outlet`
//! faces are pressure openings that pass flow in either direction (pair them
//! with `ThermalFace::Outflow { backflow_temperature }` for the ambient).
//!
//! # Coupling
//!
//! The momentum equations carry the Boussinesq body force
//! `f = -rho beta (T - T_ref) g` on fluid cells, with pressure measured from
//! the reference hydrostatic state, so an opening at pressure zero is ambient
//! at `T_ref`. The coupling starts from conduction through still fluid with
//! the openings held at their ambient temperature. Each coupling runs `sweeps_per_coupling` SIMPLEC iterations
//! under the current force, then re-solves the energy equation on the
//! current fluxes. Convergence requires, at the same coupling, the SIMPLEC
//! mass and momentum residuals below `flow.tolerance` and the largest
//! temperature change of the energy re-solve below `temperature_tolerance`
//! times the temperature span. The returned flow is projected tightly and the
//! returned energy solution is computed on exactly those fluxes.
//!
//! # No-claim boundaries
//!
//! Steady laminar Boussinesq flow (`beta |T - T_ref|` small); a configuration
//! above the transition to unsteady convection refuses as `FlowNotSteady`
//! rather than returning a time average. Radiation only as the optional
//! emission to the surroundings; turbulence only through the optional LVEL
//! closure. The SIMPLEC no-claims apply (staircase walls, power-law
//! momentum convection).

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, VoxelDomain};
use super::energy::{EnergyConfig, EnergySolution, ThermalFace, ThermalSetup, solve_energy};
use super::flow::FlowField;
use super::radiation::{RadiationConfig, Radiators};
use super::simple::{FvFlow, SimpleConfig, Solver, Turbulence, admit};
use super::turbulence::TURBULENT_PRANDTL;
use super::{ChtError, finite, finite_positive, poll};

/// Boussinesq coupling controls.
#[derive(Debug, Clone, PartialEq)]
pub struct FvBuoyancyConfig {
    /// Gravitational acceleration vector, m/s^2.
    pub gravity_m_s2: [f64; 3],
    /// Volumetric expansion coefficient `beta`, 1/K.
    pub expansion_per_k: f64,
    /// Temperature of the reference (ambient) state, K.
    pub reference_temperature_k: f64,
    /// Flow rules and SIMPLEC controls; `tolerance` is the flow criterion.
    pub flow: SimpleConfig,
    /// Energy solver controls.
    pub energy: EnergyConfig,
    /// SIMPLEC iterations between energy re-solves.
    pub sweeps_per_coupling: usize,
    /// Coupling budget.
    pub max_couplings: usize,
    /// Largest temperature change of one energy re-solve over the
    /// temperature span, at convergence.
    pub temperature_tolerance: f64,
    /// Optional surface radiation to the surroundings (each coupling's
    /// energy solve carries the sinks linearized about the current
    /// temperature, so the coupling criterion also converges radiation).
    pub radiation: Option<RadiationConfig>,
}

impl FvBuoyancyConfig {
    /// Defaults for gravity `gravity_m_s2`, expansion `beta` about
    /// `reference_temperature_k`, with the given flow rules.
    #[must_use]
    pub fn new(
        gravity_m_s2: [f64; 3],
        expansion_per_k: f64,
        reference_temperature_k: f64,
        flow: SimpleConfig,
    ) -> Self {
        Self {
            gravity_m_s2,
            expansion_per_k,
            reference_temperature_k,
            flow,
            energy: EnergyConfig::default(),
            sweeps_per_coupling: 5,
            max_couplings: 4000,
            temperature_tolerance: 1e-7,
            radiation: None,
        }
    }
}

/// Run evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct FvNaturalConvectionReport {
    /// Energy re-solves.
    pub couplings: usize,
    /// SIMPLEC iterations in total.
    pub sweeps: usize,
    /// Mass residual of the last sweep.
    pub mass_residual: f64,
    /// Momentum residual of the last sweep.
    pub momentum_residual: f64,
    /// Largest temperature change of the last energy re-solve over the
    /// temperature span.
    pub temperature_change: f64,
    /// Net heat radiated to the surroundings at the returned temperatures,
    /// W, when radiation is declared.
    pub radiated_w: Option<f64>,
}

/// Steady natural convection: flow, temperature and evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct FvNaturalConvection {
    /// The converged flow (its report covers the final projection).
    pub flow: FvFlow,
    /// Temperature on exactly `flow.field`.
    pub energy: EnergySolution,
    /// Coupling evidence.
    pub report: FvNaturalConvectionReport,
}

fn body_force(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    config: &FvBuoyancyConfig,
    temperature: &[f64],
) -> Vec<[f64; 3]> {
    (0..domain.cell_count())
        .map(|c| {
            if !domain.is_fluid(c) {
                return [0.0; 3];
            }
            let scale = -fluid.density_kg_m3
                * config.expansion_per_k
                * (temperature[c] - config.reference_temperature_k);
            config.gravity_m_s2.map(|g| scale * g)
        })
        .collect()
}

/// Solve steady Boussinesq natural (or mixed) convection by alternating
/// SIMPLEC sweeps and conjugate energy solves.
///
/// # Errors
/// Input refusals, [`ChtError::FlowNotSteady`] when the coupling budget ends
/// first, [`ChtError::FlowDiverged`], energy-solver refusals (for example
/// flow through a closed thermal face), or [`ChtError::Cancelled`].
pub fn fv_natural_convection(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &FvBuoyancyConfig,
    gate: &CancelGate,
) -> Result<FvNaturalConvection, ChtError> {
    admit(domain, fluid, &config.flow)?;
    for g in config.gravity_m_s2 {
        finite("buoyancy.gravity_m_s2", g)?;
    }
    finite("buoyancy.expansion_per_k", config.expansion_per_k)?;
    if config.expansion_per_k < 0.0 {
        return Err(ChtError::InvalidInput {
            field: "buoyancy.expansion_per_k",
            reason: "must be non-negative".into(),
        });
    }
    finite(
        "buoyancy.reference_temperature_k",
        config.reference_temperature_k,
    )?;
    finite_positive(
        "buoyancy.temperature_tolerance",
        config.temperature_tolerance,
    )?;
    if config.sweeps_per_coupling == 0 {
        return Err(ChtError::InvalidInput {
            field: "buoyancy.sweeps_per_coupling",
            reason: "must be at least one".into(),
        });
    }
    let exposed = match &config.radiation {
        Some(radiation) => Some(Radiators::build(domain, solids, radiation, gate)?),
        None => None,
    };
    // Each energy solve carries the current radiative sinks and, with a
    // turbulence closure, the current eddy conductivity.
    let energy_setup = |temperature: &[f64], eddy_viscosity: &[f64]| {
        let mut step = setup.clone();
        if let Some(radiators) = &exposed {
            step.cell_sinks.extend(radiators.sinks(domain, temperature));
        }
        if config.flow.turbulence != Turbulence::Laminar {
            step.eddy_conductivity_w_m_k = eddy_viscosity
                .iter()
                .map(|nu_t| fluid.volumetric_heat_capacity() * nu_t / TURBULENT_PRANDTL)
                .collect();
        }
        step
    };
    let mut solver = Solver::new(domain, fluid, &config.flow);
    // Seed the coupling with conduction through still fluid, the openings
    // held at their ambient temperature: a uniform reference field would
    // carry no buoyancy, hence no flow, and openings without flow anchor no
    // temperature.
    let mut seed = setup.clone();
    for face in &mut seed.faces {
        if let ThermalFace::Outflow {
            backflow_temperature,
        } = *face
        {
            *face = ThermalFace::Temperature(backflow_temperature);
        }
    }
    let mut temperature = solve_energy(
        domain,
        fluid,
        solids,
        &FlowField::quiescent(domain),
        &seed,
        &config.energy,
        gate,
    )?
    .temperature;
    let mut residuals = (f64::INFINITY, f64::INFINITY);
    let mut temperature_change = f64::INFINITY;
    let (mut couplings, mut sweeps) = (0usize, 0usize);
    let tolerance = config.flow.tolerance;
    while couplings < config.max_couplings {
        poll(gate)?;
        couplings += 1;
        solver.set_force(body_force(domain, fluid, config, &temperature));
        for _ in 0..config.sweeps_per_coupling {
            sweeps += 1;
            residuals = solver.sweep(gate)?;
            if !(residuals.0.is_finite() && residuals.1.is_finite()) {
                return Err(ChtError::FlowDiverged { step: sweeps });
            }
        }
        let energy = solve_energy(
            domain,
            fluid,
            solids,
            &solver.field(),
            &energy_setup(&temperature, &solver.heat_eddy_viscosity()),
            &config.energy,
            gate,
        )?;
        let (lo, hi) = energy
            .temperature
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(*t), hi.max(*t))
            });
        let change = energy
            .temperature
            .iter()
            .zip(&temperature)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        temperature_change = change / (hi - lo).max(f64::MIN_POSITIVE);
        temperature = energy.temperature;
        if residuals.0 <= tolerance
            && residuals.1 <= tolerance
            && temperature_change <= config.temperature_tolerance
        {
            break;
        }
    }
    if residuals.0 > tolerance
        || residuals.1 > tolerance
        || temperature_change > config.temperature_tolerance
    {
        return Err(ChtError::FlowNotSteady {
            steps: sweeps,
            last_change: residuals.0.max(residuals.1).max(temperature_change),
            tolerance,
        });
    }
    let flow = solver.finish(fluid, sweeps, residuals, gate)?;
    let energy = solve_energy(
        domain,
        fluid,
        solids,
        &flow.field,
        &energy_setup(&temperature, &flow.heat_eddy_viscosity_m2_s),
        &config.energy,
        gate,
    )?;
    let radiated_w = exposed
        .as_ref()
        .map(|radiators| radiators.radiated(domain, &energy.temperature));
    Ok(FvNaturalConvection {
        flow,
        energy,
        report: FvNaturalConvectionReport {
            couplings,
            sweeps,
            mass_residual: residuals.0,
            momentum_residual: residuals.1,
            temperature_change,
            radiated_w,
        },
    })
}
