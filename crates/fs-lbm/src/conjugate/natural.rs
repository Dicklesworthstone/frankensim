//! Two-way Boussinesq coupling: buoyancy-driven D3Q19 flow in a closed
//! enclosure and the conservative conjugate energy equation.
//!
//! The coupled problem is marched in pseudo time to its steady fixed point:
//!
//! 1. the Guo body force of every fluid cell is the Boussinesq buoyancy
//!    `-beta (T - T_ref) g`, converted to lattice units;
//! 2. the lattice advances `coupling_interval` steps;
//! 3. the cell MASS fluxes are projected onto divergence-free face fluxes;
//! 4. one implicit (backward-Euler) energy step of length
//!    the actual lattice batch updates the temperature. Solid cells
//!    use the fluid's volumetric heat capacity as a PSEUDO capacity: it sets
//!    the relaxation path, never the steady answer.
//!
//! A per-coupling change understates the distance to the fixed point when
//! each coupling moves only a small fraction of the way, so convergence is
//! judged on the geometric extrapolation of the remaining change,
//! `delta_k r / (1 - r)` with the observed contraction `r = delta_k /
//! delta_(k-1)` (no estimate while the sequence is not contracting): both
//! the relative mass-flux estimate and the temperature estimate (relative to
//! the temperature range) must fall below the steady tolerance. A candidate
//! STEADY energy solve on those same projected fluxes must also agree with
//! both the temperature that drove the lattice batch and its pseudo-time
//! endpoint, to that tolerance. A failed endpoint check continues the coupled
//! march within the original step budget; it cannot publish a new thermal
//! field that has never driven the accepted flow. The accepted field, energy
//! solve and projection are returned together without another solve.
//!
//! No-claims: steady laminar natural convection only (no claim of physical
//! steadiness above the transition Rayleigh number; a run that does not
//! settle refuses as `FlowNotSteady`); Boussinesq with constant properties;
//! closed enclosures (walls or periodic axes); BGK collision, so the
//! relaxation time must sit comfortably above one half.

use fs_exec::{CancelGate, TilePool};

use super::domain::{FluidProperties, SolidMaterial, VoxelDomain};
use super::energy::{
    EnergyConfig, EnergySolution, PseudoStep, ThermalFace, ThermalSetup, solve_energy_inner,
};
use super::flow::{FlowFace, FlowField, ProjectionReport};
use super::{ChtError, finite, finite_positive, poll};
use crate::d3q19::{
    BoundaryGrid3, BoundarySpec3, BoundaryStepError3, CollisionModel3, E3, FaceBoundary3, TILE,
};

/// Natural-convection run configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuoyancyConfig {
    /// Gravitational acceleration vector, m/s^2 (for example `[0, 0, -9.81]`).
    pub gravity_m_s2: [f64; 3],
    /// Volumetric thermal expansion coefficient `beta`, 1/K (`1/T` for an
    /// ideal gas).
    pub expansion_per_k: f64,
    /// Boussinesq reference temperature, K.
    pub reference_temperature_k: f64,
    /// BGK relaxation time; fixes the lattice viscosity and so the time step.
    pub tau: f64,
    /// Axes the LATTICE treats as periodic. The energy equation treats the
    /// same faces as adiabatic, which is exact only for solutions invariant
    /// along that axis (for example a 2-D cavity extruded along it); the
    /// corresponding thermal faces must be `Adiabatic`.
    pub periodic: [bool; 3],
    /// Lattice steps per energy update.
    pub coupling_interval: usize,
    /// Lattice step budget.
    pub max_steps: usize,
    /// Steady criterion on the EXTRAPOLATED remaining change (relative
    /// mass flux, and temperature over the temperature range), and on the
    /// steady thermal endpoint's difference from the driving and pseudo-time
    /// temperatures, normalized by the returned field's temperature range.
    pub steady_tolerance: f64,
    /// Largest admitted lattice speed (compressibility guard).
    pub max_lattice_speed: f64,
    /// Pool workers; `0` means every host core.
    pub workers: usize,
    /// Projection tolerance.
    pub projection_tolerance: f64,
    /// Energy solver settings.
    pub energy: EnergyConfig,
}

impl Default for BuoyancyConfig {
    fn default() -> Self {
        Self {
            gravity_m_s2: [0.0, 0.0, -9.81],
            expansion_per_k: 1.0 / 300.0,
            reference_temperature_k: 300.0,
            tau: 0.8,
            periodic: [false; 3],
            coupling_interval: 20,
            max_steps: 400_000,
            steady_tolerance: 1e-7,
            max_lattice_speed: 0.15,
            workers: 0,
            projection_tolerance: 1e-12,
            energy: EnergyConfig::default(),
        }
    }
}

/// Run evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct NaturalConvectionReport {
    /// Lattice steps executed.
    pub steps: usize,
    /// Energy updates executed.
    pub couplings: usize,
    /// Relative mass-flux change at the last coupling.
    pub last_flux_change: f64,
    /// Temperature change at the last coupling over the temperature range.
    pub last_temperature_change: f64,
    /// Extrapolated remaining relative mass-flux change at the stop.
    pub remaining_flux_change: f64,
    /// Extrapolated remaining temperature change over the range at the stop.
    pub remaining_temperature_change: f64,
    /// Largest difference of the returned steady thermal field from either
    /// the temperature that drove the final lattice batch or its pseudo-time
    /// endpoint, K. Dividing by the returned temperature range satisfies
    /// `BuoyancyConfig::steady_tolerance` (exact equality for a uniform field).
    pub coupling_residual_k: f64,
    /// Lattice kinematic viscosity `(tau - 1/2)/3`.
    pub lattice_viscosity: f64,
    /// Metres per second per lattice velocity unit.
    pub velocity_scale_m_s: f64,
    /// Seconds per lattice step.
    pub time_step_s: f64,
    /// Largest lattice speed in the converged field.
    pub max_lattice_speed: f64,
    /// Projection of the final field.
    pub projection: ProjectionReport,
}

/// Steady natural convection with conjugate conduction.
#[derive(Debug, Clone, PartialEq)]
pub struct NaturalConvection {
    /// Projected face fluxes of the final flow.
    pub field: FlowField,
    /// Cell mass flux over reference density, m/s (zero in solids).
    pub velocity_m_s: Vec<[f64; 3]>,
    /// Steady energy solution on the final fluxes.
    pub energy: EnergySolution,
    /// Run evidence.
    pub report: NaturalConvectionReport,
}

/// Solve steady Boussinesq natural convection with conjugate conduction in a
/// closed enclosure.
///
/// Requirements: every dimension a positive multiple of four; no `Inflow` /
/// `Outflow` thermal faces; thermal faces on periodic lattice axes must be
/// `Adiabatic`.
///
/// # Errors
/// Input refusals, [`ChtError::LatticeResolution`] when the buoyant flow
/// would exceed `max_lattice_speed`, [`ChtError::FlowDiverged`],
/// [`ChtError::FlowNotSteady`], solver refusals, or [`ChtError::Cancelled`].
#[allow(clippy::too_many_lines)] // admission -> conduction start -> coupled march -> endpoint acceptance
pub fn natural_convection(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &BuoyancyConfig,
    gate: &CancelGate,
) -> Result<NaturalConvection, ChtError> {
    domain.require_uniform("the lattice-Boltzmann natural convection")?;
    fluid.validate()?;
    for g in config.gravity_m_s2 {
        finite("buoyancy.gravity_m_s2", g)?;
    }
    finite("buoyancy.expansion_per_k", config.expansion_per_k)?;
    finite(
        "buoyancy.reference_temperature_k",
        config.reference_temperature_k,
    )?;
    finite_positive("buoyancy.steady_tolerance", config.steady_tolerance)?;
    finite_positive("buoyancy.max_lattice_speed", config.max_lattice_speed)?;
    if !(config.tau.is_finite() && config.tau >= 0.51) {
        return Err(ChtError::InvalidInput {
            field: "buoyancy.tau",
            reason: format!("BGK buoyancy needs tau >= 0.51, got {}", config.tau),
        });
    }
    if config.coupling_interval == 0 {
        return Err(ChtError::InvalidInput {
            field: "buoyancy.coupling_interval",
            reason: "must be positive".into(),
        });
    }
    let [nx, ny, nz] = domain.dims();
    if [nx, ny, nz].iter().any(|n| !n.is_multiple_of(TILE)) {
        return Err(ChtError::InvalidDomain {
            reason: format!("LBM dimensions must be multiples of {TILE}, got {nx}x{ny}x{nz}"),
        });
    }
    if domain.fluid_count() == 0 {
        return Err(ChtError::InvalidDomain {
            reason: "no fluid cell".into(),
        });
    }
    for (f, rule) in setup.faces.iter().enumerate() {
        if matches!(
            rule,
            ThermalFace::Inflow { .. } | ThermalFace::Outflow { .. }
        ) {
            return Err(ChtError::InvalidInput {
                field: "thermal.faces",
                reason: "natural convection runs in a closed enclosure; open faces refuse".into(),
            });
        }
        if config.periodic[f / 2] && *rule != ThermalFace::Adiabatic {
            return Err(ChtError::InvalidInput {
                field: "thermal.faces",
                reason: format!("face {f} lies on a periodic lattice axis and must be Adiabatic"),
            });
        }
    }

    let lattice_viscosity = (config.tau - 0.5) / 3.0;
    let dx = domain.dx();
    let velocity_scale_m_s = fluid.kinematic_viscosity_m2_s / (lattice_viscosity * dx);
    let time_step_s = dx / velocity_scale_m_s;
    // Lattice acceleration per kelvin of excess temperature.
    let accel_per_k = config
        .gravity_m_s2
        .map(|g| -config.expansion_per_k * g * dx / (velocity_scale_m_s * velocity_scale_m_s));

    let wall = FaceBoundary3::stationary_wall();
    let face = |axis: usize| {
        if config.periodic[axis] {
            FaceBoundary3::Periodic
        } else {
            wall
        }
    };
    let spec = BoundarySpec3::new([face(0), face(0), face(1), face(1), face(2), face(2)]);
    let mut grid = BoundaryGrid3::with_collision_model(
        nx,
        ny,
        nz,
        CollisionModel3::Bgk { tau: config.tau },
        [0.0; 3],
        spec,
    );
    grid.voxelize_sdf(|p| {
        let (x, y, z) = (p[0] as usize, p[1] as usize, p[2] as usize);
        if domain.is_fluid(domain.index(x, y, z)) {
            1.0
        } else {
            -1.0
        }
    });
    let fluid_cells: Vec<usize> = (0..domain.cell_count())
        .filter(|&c| domain.is_fluid(c))
        .collect();
    let flow_faces = [FlowFace::Wall; 6];

    // Start from steady conduction.
    let start = solve_energy_inner(
        domain,
        fluid,
        solids,
        &FlowField::quiescent(domain),
        setup,
        &config.energy,
        None,
        gate,
    )?;
    let mut temperature = start.temperature;
    let cell_capacity_j_k = fluid.volumetric_heat_capacity() * dx * dx * dx;
    let mut storage = vec![
        cell_capacity_j_k / (config.coupling_interval as f64 * time_step_s);
        domain.cell_count()
    ];

    let workers = if config.workers == 0 {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    } else {
        config.workers
    };
    let pool = TilePool::for_host(workers, 0);
    let sample = |grid: &BoundaryGrid3, step: usize| -> Result<(Vec<[f64; 3]>, f64), ChtError> {
        let mut velocities = vec![[0.0f64; 3]; domain.cell_count()];
        let mut max_speed_squared = 0.0f64;
        for &c in &fluid_cells {
            let [x, y, z] = domain.coords(c);
            let mut m = [0.0f64; 4];
            for (q, fq) in grid.populations(x, y, z).into_iter().enumerate() {
                m[0] += f64::from(E3[q].0) * fq;
                m[1] += f64::from(E3[q].1) * fq;
                m[2] += f64::from(E3[q].2) * fq;
                m[3] += fq;
            }
            if !(m.iter().all(|v| v.is_finite()) && m[3] > 0.0) {
                return Err(ChtError::FlowDiverged { step });
            }
            // Projection consumes mass flux over reference density (rho0=1),
            // while the low-Mach admission applies to physical velocity.
            // Both must include the Guo half-source correction.
            let force = grid.force_density(x, y, z);
            let momentum = core::array::from_fn(|axis| m[axis] + 0.5 * force[axis]);
            let velocity = momentum.map(|value| value / m[3]);
            if !momentum
                .iter()
                .chain(&velocity)
                .all(|value| value.is_finite())
            {
                return Err(ChtError::FlowDiverged { step });
            }
            max_speed_squared = max_speed_squared.max(
                velocity
                    .iter()
                    .fold(0.0, |sum, value| value.mul_add(*value, sum)),
            );
            velocities[c] = momentum;
        }
        Ok((velocities, fs_math::det::sqrt(max_speed_squared)))
    };
    let mut steps = 0usize;
    let mut couplings = 0usize;
    let mut last_flux_change = f64::INFINITY;
    let mut last_temperature_change = f64::INFINITY;
    let mut lattice_momentum = vec![[0.0f64; 3]; domain.cell_count()];
    let mut remaining_flux_change = f64::INFINITY;
    let mut remaining_temperature_change = f64::INFINITY;
    // Geometric tail of a contracting sequence of per-coupling changes.
    let remaining = |current: f64, previous: f64| {
        if !previous.is_finite() {
            return f64::INFINITY;
        }
        if current == 0.0 {
            return 0.0;
        }
        let ratio = current / previous;
        if ratio.is_finite() && ratio < 1.0 {
            current * ratio / (1.0 - ratio)
        } else {
            f64::INFINITY
        }
    };
    pool.with_parked_crew_local(|parked| -> Result<NaturalConvection, ChtError> {
        while steps < config.max_steps {
            poll(gate)?;
            let force = |cell: [usize; 3]| {
                let excess = temperature[domain.index(cell[0], cell[1], cell[2])]
                    - config.reference_temperature_k;
                accel_per_k.map(|a| a * excess)
            };
            grid.rebind_force_field(force);
            let batch = config.coupling_interval.min(config.max_steps - steps);
            if batch != config.coupling_interval {
                // The last batch can be shortened by the original step
                // budget. Its heat step must use that same physical duration.
                storage.fill(cell_capacity_j_k / (batch as f64 * time_step_s));
            }
            for _ in 0..batch {
                grid.step_pooled(parked, gate)
                    .map_err(|error| match error {
                        BoundaryStepError3::Cancelled => ChtError::Cancelled,
                        BoundaryStepError3::Collision { .. }
                        | BoundaryStepError3::Unphysical { .. }
                        | BoundaryStepError3::WallMomentum { .. }
                        | BoundaryStepError3::FiniteMassWall { .. }
                        | BoundaryStepError3::EndpointForce { .. } => {
                            ChtError::FlowDiverged { step: steps }
                        }
                        BoundaryStepError3::Pool(detail) => ChtError::Executor { detail },
                    })?;
                steps += 1;
            }
            let (momentum, speed) = sample(&grid, steps)?;
            let (mut diff, mut norm) = (0.0f64, 0.0f64);
            for (a, b) in momentum.iter().zip(&lattice_momentum) {
                for k in 0..3 {
                    diff = (a[k] - b[k]).mul_add(a[k] - b[k], diff);
                    norm = a[k].mul_add(a[k], norm);
                }
            }
            if speed > config.max_lattice_speed {
                return Err(ChtError::LatticeResolution {
                    tau: config.tau,
                    cell_reynolds: speed * dx * velocity_scale_m_s / fluid.kinematic_viscosity_m2_s,
                    remedy: format!(
                        "buoyant lattice speed {speed:.3} exceeds {}; lower tau or refine dx",
                        config.max_lattice_speed
                    ),
                });
            }
            let flux_change = fs_math::det::sqrt(diff / norm.max(f64::MIN_POSITIVE));
            // The geometric ratio presumes equal elapsed times. A shortened
            // final batch must not manufacture a small tail just by moving
            // less time; an exactly unchanged field retains its zero limit.
            remaining_flux_change = if batch == config.coupling_interval || flux_change == 0.0 {
                remaining(flux_change, last_flux_change)
            } else {
                f64::INFINITY
            };
            last_flux_change = flux_change;
            lattice_momentum = momentum;
            let physical: Vec<[f64; 3]> = lattice_momentum
                .iter()
                .map(|m| m.map(|v| v * velocity_scale_m_s))
                .collect();
            let (field, projection) = FlowField::from_cell_velocities(
                domain,
                &physical,
                flow_faces,
                config.projection_tolerance,
                gate,
            )?;
            let next = solve_energy_inner(
                domain,
                fluid,
                solids,
                &field,
                setup,
                &config.energy,
                Some(&PseudoStep {
                    coefficient_w_k: &storage,
                    previous: &temperature,
                    junction_coefficient_w_k: &[],
                }),
                gate,
            )?;
            let scale = temperature_range(&next.temperature);
            let temperature_change = next
                .temperature
                .iter()
                .zip(&temperature)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
                / scale;
            remaining_temperature_change =
                if batch == config.coupling_interval || temperature_change == 0.0 {
                    remaining(temperature_change, last_temperature_change)
                } else {
                    f64::INFINITY
                };
            last_temperature_change = temperature_change;
            couplings += 1;
            if remaining_flux_change <= config.steady_tolerance
                && remaining_temperature_change <= config.steady_tolerance
            {
                let energy = solve_energy_inner(
                    domain,
                    fluid,
                    solids,
                    &field,
                    setup,
                    &config.energy,
                    None,
                    gate,
                )?;
                // The lattice still carries the force from `temperature`,
                // not from `next`. Checking only the last pseudo-time state
                // would miss the force change of this very coupling.
                let coupling_residual_k =
                    steady_temperature_defect(&energy.temperature, &temperature, &next.temperature);
                let coupling_residual =
                    coupling_residual_k / temperature_range(&energy.temperature);
                if coupling_residual.is_finite() && coupling_residual <= config.steady_tolerance {
                    poll(gate)?;
                    return Ok(NaturalConvection {
                        field,
                        velocity_m_s: physical,
                        energy,
                        report: NaturalConvectionReport {
                            steps,
                            couplings,
                            last_flux_change,
                            last_temperature_change,
                            remaining_flux_change,
                            remaining_temperature_change,
                            coupling_residual_k,
                            lattice_viscosity,
                            velocity_scale_m_s,
                            time_step_s,
                            max_lattice_speed: speed,
                            projection,
                        },
                    });
                }
                // A geometric tail estimate is not an endpoint residual.
                // Keep the stable pseudo-time trajectory and restart its
                // contraction history instead of jumping to an unaccepted
                // steady temperature or extending the caller's step budget.
                remaining_temperature_change = remaining_temperature_change.max(coupling_residual);
                last_flux_change = f64::INFINITY;
                last_temperature_change = f64::INFINITY;
            }
            temperature = next.temperature;
        }
        Err(ChtError::FlowNotSteady {
            steps,
            last_change: remaining_flux_change.max(remaining_temperature_change),
            tolerance: config.steady_tolerance,
        })
    })
}

fn temperature_range(temperature: &[f64]) -> f64 {
    let (lo, hi) = temperature
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &value| {
            (lo.min(value), hi.max(value))
        });
    (hi - lo).max(f64::MIN_POSITIVE)
}

// All three fields come from the same admitted domain and have equal lengths.
fn steady_temperature_defect(steady: &[f64], driving: &[f64], pseudo: &[f64]) -> f64 {
    steady
        .iter()
        .zip(driving)
        .zip(pseudo)
        .map(|((&endpoint, &before), &after)| {
            (endpoint - before).abs().max((endpoint - after).abs())
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::{steady_temperature_defect, temperature_range};

    #[test]
    fn g0_steady_endpoint_includes_the_temperature_that_drove_buoyancy() {
        let steady = [301.0, 303.0, 305.0];
        let previous = [300.0, 303.0, 306.0];
        // The heat equation has already equilibrated at the projected flow,
        // but the force carried by that flow is still one kelvin out of date.
        assert_eq!(steady_temperature_defect(&steady, &previous, &steady), 1.0);
        // Conversely, a still-relaxing pseudo-time state cannot be omitted.
        assert_eq!(steady_temperature_defect(&steady, &steady, &previous), 1.0);
        assert_eq!(temperature_range(&steady), 4.0);
    }

    #[test]
    fn g3_thermal_endpoint_normalization_uses_the_span_not_an_absolute_floor() {
        let endpoint = [0.25, 0.5, 0.75];
        let before = [0.0, 0.5, 1.0];
        for (scale, offset) in [(1.0, 0.0), (1.0, 300.0), (1e-12, 0.0)] {
            let steady = endpoint.map(|value| offset + scale * value);
            let driving = before.map(|value| offset + scale * value);
            let relative =
                steady_temperature_defect(&steady, &driving, &steady) / temperature_range(&steady);
            assert!(
                (relative - 0.5).abs() < 1e-15,
                "{scale} {offset}: {relative}"
            );
        }
        let uniform = [300.0; 3];
        assert_eq!(
            steady_temperature_defect(&uniform, &uniform, &uniform) / temperature_range(&uniform),
            0.0,
        );
    }
}
