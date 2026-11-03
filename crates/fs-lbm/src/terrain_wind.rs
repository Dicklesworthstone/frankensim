//! Native D3Q19 wind over voxelized terrain (E10.3).
//!
//! The logarithmic profile is prescribed at the inlet and used as a diagnostic
//! reference; every reported velocity comes from actual BGK evolution. This is
//! a finite-time, laminar, isothermal numerical experiment, not atmospheric LES,
//! a converged steady solution, or a validated Tier-D/fs-atmo comparison.

use fs_blake3::hash_domain;
use fs_exec::{CancelGate, KernelRunner, TilePool};

use crate::d3q19::{BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3, InletFlux3, TILE};

/// Schema version for the native terrain-flow diagnostic receipt.
pub const V04C_RECEIPT_SCHEMA_V2: &str = "org.frankensim.wf.v04c.discrepancy.v2";

/// Configuration for D3Q19 wind-over-terrain simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct TerrainWindConfig {
    /// Physical domain extent [Lx, Ly, Lz] [m].
    pub domain_size_m: [f64; 3],
    /// Inflow mean speed at 10m reference height [m/s].
    pub inflow_speed_mps: f64,
    /// Surface aerodynamic roughness length $z_0$ [m].
    pub roughness_z0_m: f64,
    /// Grid resolution [Nx, Ny, Nz].
    pub resolution: [usize; 3],
    /// Simulation timestep count.
    pub num_steps: usize,
    /// Lattice Boltzmann relaxation time $\tau$ ($\tau > 0.5$).
    pub tau: f64,
    /// Physical duration of one native lattice step [s]. Together with cubic
    /// spacing and tau this fixes the simulated kinematic viscosity.
    pub timestep_s: f64,
    /// Explicit allocation ceiling in lattice cells, including terrain solids.
    pub max_cells: usize,
    /// Declared acceptable discrepancy band [m/s].
    pub declared_band_mps: f64,
}

impl Default for TerrainWindConfig {
    fn default() -> Self {
        Self {
            domain_size_m: [200.0, 200.0, 50.0],
            inflow_speed_mps: 8.0,
            roughness_z0_m: 0.005,
            // Standard D3Q19 needs cubic cells: all three spacings are 25/6 m.
            resolution: [48, 48, 12],
            num_steps: 100,
            tau: 0.58,
            timestep_s: 0.025,
            max_cells: 262_144,
            declared_band_mps: 0.8,
        }
    }
}

/// Finite-time native LBM discrepancy against the declared inlet profile.
/// This receipt alone does not establish a matched fs-atmo comparison.
#[derive(Clone, Debug, PartialEq)]
pub struct V04cDiscrepancyReceipt {
    /// Schema version string.
    pub schema_version: &'static str,
    /// Content digest of the input run configuration.
    pub run_manifest_digest: String,
    /// RMS vector discrepancy against the prescribed flat-ground logarithmic
    /// profile, over interior fluid cells at the final time [m/s].
    pub mean_wind_rmse_mps: f64,
    /// Largest local velocity-vector discrepancy [m/s].
    pub max_discrepancy_mps: f64,
    /// Negative streamwise velocity observed in the downstream half-domain.
    /// This is a local reverse-flow diagnostic, not a separation certificate.
    pub recirculation_detected: bool,
    /// Declared tolerance band [m/s].
    pub declared_band_mps: f64,
    /// Whether both RMS and maximum discrepancies meet the declared band.
    pub within_declared_band: bool,
    /// Physical cubic-cell spacing [m].
    pub cell_size_m: f64,
    /// Actual physical lattice timestep [s].
    pub timestep_s: f64,
    /// Duration actually simulated from rest [s].
    pub elapsed_time_s: f64,
    /// Laminar kinematic viscosity implied by tau, dx and dt [m²/s].
    pub kinematic_viscosity_m2_s: f64,
    /// Number of interior fluid samples in the diagnostic.
    pub sampled_cells: usize,
    /// Arithmetic mean of the actual final physical velocity samples [m/s].
    pub mean_velocity_mps: [f64; 3],
    /// Largest sampled physical speed [m/s].
    pub maximum_speed_mps: f64,
    /// Smallest sampled streamwise velocity [m/s].
    pub minimum_streamwise_velocity_mps: f64,
    /// Cryptographic digest of this receipt.
    pub receipt_digest: String,
}

/// Run the native solver over a centered Gaussian hill of height `0.15 Lz`
/// and standard deviations `Lx/8, Ly/8`. The default is a reproducible synthetic
/// terrain, not a measured DEM. Use [`run_v04c_terrain_wind_with_heightfield`]
/// for authored terrain and an existing runner/cancellation scope.
///
/// # Errors
/// Refuses invalid physical scaling, grid/terrain admission or native evolution.
pub fn run_v04c_terrain_wind(config: &TerrainWindConfig) -> Result<V04cDiscrepancyReceipt, String> {
    let dx = admit_config(config)?;
    let [nx, ny, _] = config.resolution;
    let [lx, ly, lz] = config.domain_size_m;
    let mut heights = Vec::with_capacity(nx * ny);
    for y in 0..ny {
        for x in 0..nx {
            let sx = ((x as f64 + 0.5) * dx - 0.5 * lx) / (lx / 8.0);
            let sy = ((y as f64 + 0.5) * dx - 0.5 * ly) / (ly / 8.0);
            heights.push(0.15 * lz * (-0.5 * (sx * sx + sy * sy)).exp());
        }
    }
    run_v04c_terrain_wind_with_heightfield(
        config,
        &heights,
        &TilePool::for_host(1, 0),
        &CancelGate::new(),
    )
}

/// Evolve from rest above an authored row-major `(y * nx + x)` heightfield [m].
/// Cell centers below the terrain are stationary bounce-back solids. The x-min
/// face prescribes a roughness-shifted logarithmic profile measured above the
/// local inlet terrain; x-max prescribes reference density one. The y direction
/// is periodic. The bottom is stationary and the top moves at the flat-ground
/// profile's top speed. Thus the top is an explicit moving lid, not an open sky.
///
/// The existing regularized boundary retains its partial reconstruction at
/// wall/open seams. The diagnostic excludes the inlet/outlet and their adjacent
/// cell columns. No wall function, turbulence closure, equilibrium spin-up or
/// fs-atmo statistical calibration is inferred from the roughness length.
///
/// `tau` supplies numerical viscosity: `nu = (tau - 1/2) dx²/(3 dt)`.
/// Choosing a physically appropriate Reynolds number remains the caller's job.
/// Native populations, topology and reusable work buffers are created once.
///
/// # Errors
/// Refuses nonfinite/unsupported input, noncubic spacing, an exceeded cell
/// budget, obstructed open-face neighbors, low-Mach departure or cancellation.
pub fn run_v04c_terrain_wind_with_heightfield<P: KernelRunner>(
    config: &TerrainWindConfig,
    heights_m: &[f64],
    runner: &P,
    gate: &CancelGate,
) -> Result<V04cDiscrepancyReceipt, String> {
    let dx = admit_config(config)?;
    let [nx, ny, nz] = config.resolution;
    if gate.is_requested() {
        return Err("terrain wind cancelled before initialization".into());
    }
    if heights_m.len() != nx * ny
        || heights_m
            .iter()
            .any(|h| !h.is_finite() || *h < 0.0 || *h > config.domain_size_m[2] - 2.0 * dx)
    {
        return Err(
            "heightfield needs nx*ny finite nonnegative heights with two clear top cells".into(),
        );
    }
    // The native open boundary samples its first interior neighbor. Admit the
    // exact voxel relation before its initialization-only topology API is called.
    for y in 0..ny {
        for z in 0..nz {
            let height = (z as f64 + 0.5) * dx;
            for (edge, neighbor) in [(0, 1), (nx - 1, nx - 2)] {
                if height > heights_m[y * nx + edge] && height <= heights_m[y * nx + neighbor] {
                    return Err(
                        "terrain obstructs an open face's first interior fluid neighbor".into(),
                    );
                }
            }
        }
    }
    let velocity_unit = dx / config.timestep_s;
    let top_speed = profile(config, config.domain_size_m[2]) / velocity_unit;
    let boundaries = BoundarySpec3::new([
        FaceBoundary3::Velocity { velocity: [0.0; 3] },
        FaceBoundary3::Pressure { density: 1.0 },
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::Wall {
            velocity: [top_speed, 0.0, 0.0],
        },
    ]);
    let mut grid = BoundaryGrid3::new(nx, ny, nz, config.tau, [0.0; 3], boundaries);
    grid.voxelize_sdf(|point| {
        point[2] * dx - heights_m[point[1] as usize * nx + point[0] as usize]
    });
    grid.set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, gate, |cell| {
        let height = (cell[2] as f64 + 0.5) * dx - heights_m[cell[1] * nx];
        Ok([profile(config, height) / velocity_unit, 0.0, 0.0])
    })
    .map_err(|error| error.to_string())?;
    for step in 0..config.num_steps {
        grid.step_pooled(runner, gate)
            .map_err(|error| format!("terrain wind step {step}: {error}"))?;
        // A necessary validity gate on the state actually evolved, including
        // boundary cells. A finite but supersonic lattice cannot mint a receipt.
        for z in 0..nz {
            if gate.is_requested() {
                return Err("terrain wind cancelled during velocity admission".into());
            }
            for y in 0..ny {
                for x in 0..nx {
                    if grid.is_solid(x, y, z) {
                        continue;
                    }
                    let rho = grid.density(x, y, z);
                    if !rho.is_finite() || rho <= 0.0 {
                        return Err(format!(
                            "terrain wind lost positive density at [{x},{y},{z}]"
                        ));
                    }
                    let u = grid.velocity(x, y, z);
                    let speed_squared: f64 = u.iter().map(|v| v * v).sum();
                    if !speed_squared.is_finite() || speed_squared >= 0.03 {
                        return Err(format!(
                            "terrain wind exceeded the low-Mach bound at [{x},{y},{z}]"
                        ));
                    }
                }
            }
        }
    }

    // Bind exact inputs, including tau/dt and the supplied terrain samples.
    // Decimal formatting formerly lost distinctions and omitted live controls.
    let mut manifest = Vec::new();
    for value in config
        .domain_size_m
        .into_iter()
        .chain([
            config.inflow_speed_mps,
            config.roughness_z0_m,
            config.tau,
            config.timestep_s,
            config.declared_band_mps,
        ])
        .chain(heights_m.iter().copied())
    {
        manifest.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    for value in config
        .resolution
        .into_iter()
        .chain([config.num_steps, config.max_cells])
    {
        manifest.extend_from_slice(&(value as u64).to_le_bytes());
    }
    let run_manifest_digest = hash_domain("org.frankensim.wf.v04c.manifest.v2", &manifest)
        .to_hex()
        .to_string();
    let mut sum_sq_err = 0.0;
    let mut max_err: f64 = 0.0;
    let mut mean_velocity = [0.0; 3];
    let mut maximum_speed: f64 = 0.0;
    let mut minimum_streamwise = f64::INFINITY;
    let mut sampled_cells = 0;
    let mut recirculation = false;
    for z in 0..nz {
        if gate.is_requested() {
            return Err("terrain wind cancelled during final velocity observation".into());
        }
        let target = profile(config, (z as f64 + 0.5) * dx);
        for y in 0..ny {
            for x in 2..nx - 2 {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let velocity = grid.velocity(x, y, z).map(|u| u * velocity_unit);
                let error_squared = (velocity[0] - target).powi(2)
                    + velocity[1] * velocity[1]
                    + velocity[2] * velocity[2];
                sum_sq_err += error_squared;
                max_err = max_err.max(error_squared.sqrt());
                for axis in 0..3 {
                    mean_velocity[axis] += velocity[axis];
                }
                maximum_speed =
                    maximum_speed.max(velocity.iter().map(|u| u * u).sum::<f64>().sqrt());
                minimum_streamwise = minimum_streamwise.min(velocity[0]);
                recirculation |=
                    x >= nx / 2 && velocity[0] < -64.0 * f64::EPSILON * config.inflow_speed_mps;
                sampled_cells += 1;
            }
        }
    }
    let mean_wind_rmse_mps = (sum_sq_err / sampled_cells as f64).sqrt();
    let mean_velocity_mps = mean_velocity.map(|value| value / sampled_cells as f64);
    if ![
        mean_wind_rmse_mps,
        max_err,
        maximum_speed,
        minimum_streamwise,
    ]
    .into_iter()
    .chain(mean_velocity_mps)
    .all(f64::is_finite)
    {
        return Err("terrain velocity diagnostics are not representable in physical units".into());
    }
    let within_declared_band =
        mean_wind_rmse_mps <= config.declared_band_mps && max_err <= config.declared_band_mps;
    let mut receipt_bytes = run_manifest_digest.as_bytes().to_vec();
    for value in [
        mean_wind_rmse_mps,
        max_err,
        maximum_speed,
        minimum_streamwise,
    ]
    .into_iter()
    .chain(mean_velocity_mps)
    {
        receipt_bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    receipt_bytes.extend_from_slice(&(sampled_cells as u64).to_le_bytes());
    receipt_bytes.extend_from_slice(&[u8::from(recirculation), u8::from(within_declared_band)]);
    let receipt_digest = hash_domain("org.frankensim.wf.v04c.receipt.v2", &receipt_bytes)
        .to_hex()
        .to_string();
    Ok(V04cDiscrepancyReceipt {
        schema_version: V04C_RECEIPT_SCHEMA_V2,
        run_manifest_digest,
        mean_wind_rmse_mps,
        max_discrepancy_mps: max_err,
        recirculation_detected: recirculation,
        declared_band_mps: config.declared_band_mps,
        within_declared_band,
        cell_size_m: dx,
        timestep_s: config.timestep_s,
        elapsed_time_s: config.num_steps as f64 * config.timestep_s,
        kinematic_viscosity_m2_s: (config.tau - 0.5) / 3.0 * dx * velocity_unit,
        sampled_cells,
        mean_velocity_mps,
        maximum_speed_mps: maximum_speed,
        minimum_streamwise_velocity_mps: minimum_streamwise,
        receipt_digest,
    })
}

fn profile(config: &TerrainWindConfig, height_m: f64) -> f64 {
    config.inflow_speed_mps * (height_m / config.roughness_z0_m).ln_1p()
        / (10.0 / config.roughness_z0_m).ln_1p()
}

fn admit_config(config: &TerrainWindConfig) -> Result<f64, String> {
    if !config.tau.is_finite() || config.tau <= 0.5 {
        return Err("tau must be strictly greater than 0.5 for LBM stability".into());
    }
    if !config.inflow_speed_mps.is_finite() || config.inflow_speed_mps <= 0.0 {
        return Err("inflow_speed_mps must be strictly positive".into());
    }
    if config
        .resolution
        .iter()
        .any(|n| *n < 2 * TILE || !n.is_multiple_of(TILE))
    {
        return Err(
            "resolution must contain multiples of four with at least eight cells per axis".into(),
        );
    }
    let cells = config
        .resolution
        .iter()
        .try_fold(1usize, |count, n| count.checked_mul(*n))
        .ok_or("terrain grid cell count overflows")?;
    if cells >= u32::MAX as usize {
        return Err("terrain grid exceeds native D3Q19 pull-table addressing".into());
    }
    if cells > config.max_cells {
        return Err("terrain grid exceeds the declared max_cells budget".into());
    }
    if config
        .domain_size_m
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !config.timestep_s.is_finite()
        || config.timestep_s <= 0.0
        || !config.roughness_z0_m.is_finite()
        || config.roughness_z0_m <= 0.0
        || !config.declared_band_mps.is_finite()
        || config.declared_band_mps < 0.0
    {
        return Err("terrain dimensions, timestep and roughness must be finite and positive; band nonnegative".into());
    }
    let dx = config.domain_size_m[0] / config.resolution[0] as f64;
    for axis in 0..3 {
        let spacing = config.domain_size_m[axis] / config.resolution[axis] as f64;
        if !spacing.is_finite()
            || spacing <= 0.0
            || (spacing / dx - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(
                "native D3Q19 terrain requires cubic cells with equal physical spacing".into(),
            );
        }
    }
    let velocity_unit = dx / config.timestep_s;
    let top_speed = profile(config, config.domain_size_m[2]) / velocity_unit;
    let viscosity = (config.tau - 0.5) / 3.0 * dx * velocity_unit;
    if !velocity_unit.is_finite()
        || velocity_unit <= 0.0
        || !top_speed.is_finite()
        || top_speed <= 0.0
        || top_speed * top_speed >= 0.03
        || !viscosity.is_finite()
        || viscosity <= 0.0
        || !(config.num_steps as f64 * config.timestep_s).is_finite()
    {
        return Err(
            "terrain physical scaling is unrepresentable or exceeds the low-Mach bound".into(),
        );
    }
    Ok(dx)
}
