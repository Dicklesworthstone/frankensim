//! Gray-diffuse radiation from exposed solid surfaces to the surroundings,
//! coupled to the conjugate energy equation.
//!
//! # Model
//!
//! Every solid voxel face that borders a fluid voxel is an emitting patch of
//! area `dx^2` with the emissivity of its material. Its escape factor `F` is
//! the fraction of its diffuse (cosine-weighted) hemisphere that leaves the
//! domain through a face declared as surroundings (an opening or inlet,
//! typically) without first hitting a solid voxel or a non-surroundings
//! domain face; it is estimated by `rays_per_face` deterministic rays
//! marched voxel by voxel (Amanatides–Woo traversal) with counter-based
//! random streams keyed by `(seed, cell, face, ray)`, so the estimate is
//! bit-reproducible and its standard error is `sqrt(F (1 - F) / rays)`.
//! The patch then loses
//!
//! ```text
//! q = eps sigma F A (T^4 - T_amb^4)
//! ```
//!
//! to the surroundings temperature (the mean over its escaping rays). The
//! energy equation carries `q` as a cell sink Newton-linearized about the
//! previous iterate (`G = 4 eps sigma F A T^3`), iterated to convergence.
//!
//! # No-claim boundaries
//!
//! Radiation exchange BETWEEN surfaces (fin to fin, solid to a warm wall)
//! and re-radiation from walls are not modelled: a ray that hits a solid or
//! a wall is simply not escaping, which is exact when the obstructing
//! surfaces are at the emitter's temperature and overstates the loss
//! otherwise only through the omitted exchange. Air is transparent; gray,
//! diffuse, opaque surfaces; the escape factor is a Monte Carlo estimate.

use fs_exec::CancelGate;

use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::energy::{CellSink, EnergyConfig, EnergySolution, ThermalSetup, solve_energy};
use super::flow::FlowField;
use super::{ChtError, finite, poll};

/// Stefan–Boltzmann constant, W/(m^2 K^4) (CODATA 2018, exact).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Radiation controls.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationConfig {
    /// Emissivity of each solid material in `solids` order (`0` disables).
    pub emissivity: Vec<f64>,
    /// Surroundings temperature seen through each domain face, K, in
    /// `Face3::ALL` order; `None` is an opaque, non-exchanging face.
    pub surroundings_k: [Option<f64>; 6],
    /// Rays per exposed face.
    pub rays_per_face: usize,
    /// Seed of the counter-based ray streams.
    pub seed: u64,
    /// Picard budget of [`solve_energy_radiating`].
    pub max_iterations: usize,
    /// Largest temperature change between Picard iterates at convergence,
    /// relative to the temperature span.
    pub tolerance: f64,
}

impl RadiationConfig {
    /// Defaults for the given emissivities and surroundings.
    #[must_use]
    pub fn new(emissivity: Vec<f64>, surroundings_k: [Option<f64>; 6]) -> Self {
        Self {
            emissivity,
            surroundings_k,
            rays_per_face: 256,
            seed: 0x5EED_0FA1,
            max_iterations: 50,
            tolerance: 1e-9,
        }
    }
}

/// One emitting solid face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExposedFace {
    /// Solid cell.
    pub cell: usize,
    /// Local face of the cell (`Face3` order) bordering fluid.
    pub face: usize,
    /// Emissivity of the cell's material.
    pub emissivity: f64,
    /// Escape factor to the surroundings.
    pub escape: f64,
    /// Mean surroundings temperature of the escaping rays, K.
    pub surroundings_k: f64,
}

/// Radiation evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationReport {
    /// Emitting faces with a non-zero escape factor.
    pub exposed_faces: usize,
    /// Rays traced from the faces with a non-zero escape factor.
    pub rays: usize,
    /// Picard iterations.
    pub iterations: usize,
    /// Net heat radiated to the surroundings at the final temperatures
    /// (nonlinear law), W.
    pub radiated_w: f64,
    /// Largest temperature change of the last Picard iterate over the span.
    pub temperature_change: f64,
}

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

fn unit(x: u64) -> f64 {
    // 53 random bits in [0, 1).
    (x >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// March one ray from `origin` (inside fluid cell `start`) along `dir`;
/// returns the surroundings temperature if it leaves through a
/// surroundings face before hitting a solid or an opaque face.
fn march(
    domain: &VoxelDomain,
    surroundings: &[Option<f64>; 6],
    start: [usize; 3],
    origin: [f64; 3],
    dir: [f64; 3],
) -> Option<f64> {
    let dx = domain.dx();
    let dims = domain.dims();
    let mut cell = start.map(|v| v as i64);
    let mut t_max = [f64::INFINITY; 3];
    let mut t_delta = [f64::INFINITY; 3];
    let mut step = [0i64; 3];
    for a in 0..3 {
        if dir[a] > 0.0 {
            step[a] = 1;
            t_max[a] = ((cell[a] + 1) as f64 * dx - origin[a]) / dir[a];
            t_delta[a] = dx / dir[a];
        } else if dir[a] < 0.0 {
            step[a] = -1;
            t_max[a] = (cell[a] as f64 * dx - origin[a]) / dir[a];
            t_delta[a] = -dx / dir[a];
        }
    }
    let budget = 3 * (dims[0] + dims[1] + dims[2]) + 8;
    for _ in 0..budget {
        let axis = if t_max[0] <= t_max[1] && t_max[0] <= t_max[2] {
            0
        } else if t_max[1] <= t_max[2] {
            1
        } else {
            2
        };
        cell[axis] += step[axis];
        t_max[axis] += t_delta[axis];
        if cell[axis] < 0 || cell[axis] >= dims[axis] as i64 {
            let side = 2 * axis + usize::from(step[axis] > 0);
            return surroundings[side];
        }
        let index = domain.index(cell[0] as usize, cell[1] as usize, cell[2] as usize);
        if !domain.is_fluid(index) {
            return None;
        }
    }
    None
}

/// Escape factors of every emitting solid face (non-zero emissivity,
/// bordering fluid).
///
/// # Errors
/// Input refusals or [`ChtError::Cancelled`].
pub fn escape_factors(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    config: &RadiationConfig,
    gate: &CancelGate,
) -> Result<Vec<ExposedFace>, ChtError> {
    if config.emissivity.len() != solids.len() {
        return Err(ChtError::InvalidInput {
            field: "radiation.emissivity",
            reason: format!(
                "expected {} emissivities (one per solid), got {}",
                solids.len(),
                config.emissivity.len()
            ),
        });
    }
    for &e in &config.emissivity {
        finite("radiation.emissivity", e)?;
        if !(0.0..=1.0).contains(&e) {
            return Err(ChtError::InvalidInput {
                field: "radiation.emissivity",
                reason: format!("{e} is outside [0, 1]"),
            });
        }
    }
    for t in config.surroundings_k.iter().flatten() {
        finite("radiation.surroundings_k", *t)?;
        if *t <= 0.0 {
            return Err(ChtError::InvalidInput {
                field: "radiation.surroundings_k",
                reason: "absolute temperatures must be positive".into(),
            });
        }
    }
    if config.rays_per_face == 0 {
        return Err(ChtError::InvalidInput {
            field: "radiation.rays_per_face",
            reason: "must be at least one".into(),
        });
    }
    let dx = domain.dx();
    let mut out = Vec::new();
    for c in 0..domain.cell_count() {
        if c.is_multiple_of(1024) {
            poll(gate)?;
        }
        let Voxel::Solid(material) = domain.voxel_at(c) else {
            continue;
        };
        let emissivity = config.emissivity[usize::from(material)];
        if emissivity == 0.0 {
            continue;
        }
        let at = domain.coords(c);
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            if !domain.is_fluid(n) {
                continue;
            }
            let axis = face / 2;
            let sign = if face % 2 == 1 { 1.0 } else { -1.0 };
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let plane = if face % 2 == 1 {
                (at[axis] + 1) as f64 * dx
            } else {
                at[axis] as f64 * dx
            };
            let start = domain.coords(n);
            let (mut escaped, mut ambient_sum) = (0usize, 0.0f64);
            let stream = config
                .seed
                .wrapping_add((c as u64).wrapping_mul(6).wrapping_add(face as u64) << 20);
            for ray in 0..config.rays_per_face {
                let base = splitmix(stream.wrapping_add(ray as u64));
                let (r1, r2, r3, r4) = (
                    unit(splitmix(base ^ 1)),
                    unit(splitmix(base ^ 2)),
                    unit(splitmix(base ^ 3)),
                    unit(splitmix(base ^ 4)),
                );
                // Cosine-weighted direction about the outward normal.
                let radius = r1.sqrt();
                let phi = std::f64::consts::TAU * r2;
                let mut dir = [0.0; 3];
                dir[axis] = sign * (1.0 - r1).max(0.0).sqrt();
                dir[u] = radius * phi.cos();
                dir[v] = radius * phi.sin();
                // Uniform origin on the face, nudged into the fluid cell.
                let mut origin = [0.0; 3];
                origin[axis] = sign.mul_add(1e-9 * dx, plane);
                origin[u] = (at[u] as f64 + r3) * dx;
                origin[v] = (at[v] as f64 + r4) * dx;
                if let Some(t) = march(domain, &config.surroundings_k, start, origin, dir) {
                    escaped += 1;
                    ambient_sum += t;
                }
            }
            if escaped > 0 {
                out.push(ExposedFace {
                    cell: c,
                    face,
                    emissivity,
                    escape: escaped as f64 / config.rays_per_face as f64,
                    surroundings_k: ambient_sum / escaped as f64,
                });
            }
        }
    }
    Ok(out)
}

/// Radiative sinks of `faces` linearized about `temperature` by Newton:
/// `q(T) ~ q(T_k) + q'(T_k) (T - T_k)` with `q' = 4 eps sigma F A T_k^3`,
/// written as the sink `G (T - T_s)` with `G = q'` and
/// `T_s = T_k - q(T_k) / q'` (exact at `T = T_k`, quadratic convergence).
#[must_use]
pub fn radiative_sinks(
    domain: &VoxelDomain,
    faces: &[ExposedFace],
    temperature: &[f64],
) -> Vec<CellSink> {
    let area = domain.dx() * domain.dx();
    faces
        .iter()
        .map(|f| {
            let (t, ta) = (temperature[f.cell].max(1.0), f.surroundings_k);
            let scale = f.emissivity * STEFAN_BOLTZMANN * f.escape * area;
            let slope = 4.0 * scale * t * t * t;
            let q = scale * (t.powi(4) - ta.powi(4));
            CellSink {
                cell: f.cell,
                conductance_w_k: slope,
                temperature_k: t - q / slope.max(f64::MIN_POSITIVE),
            }
        })
        .collect()
}

/// Net radiated power at `temperature` by the nonlinear law, W.
#[must_use]
pub fn radiated_power(domain: &VoxelDomain, faces: &[ExposedFace], temperature: &[f64]) -> f64 {
    let area = domain.dx() * domain.dx();
    faces
        .iter()
        .map(|f| {
            let (t, ta) = (temperature[f.cell], f.surroundings_k);
            f.emissivity * STEFAN_BOLTZMANN * f.escape * area * (t.powi(4) - ta.powi(4))
        })
        .sum()
}

/// Steady conjugate energy with surface radiation to the surroundings on a
/// given flow, by Picard iteration of the linearized radiative sinks.
///
/// # Errors
/// Input and solver refusals, [`ChtError::SolverNotConverged`] (system
/// `"radiation"`) when the Picard budget ends first, or
/// [`ChtError::Cancelled`].
#[allow(clippy::too_many_arguments)] // physics inputs + both configurations
pub fn solve_energy_radiating(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    energy: &EnergyConfig,
    radiation: &RadiationConfig,
    gate: &CancelGate,
) -> Result<(EnergySolution, RadiationReport), ChtError> {
    let faces = escape_factors(domain, solids, radiation, gate)?;
    // Start from the sinks linearized about the surroundings temperature:
    // without radiation the first iterate can lack the dominant heat path.
    let surroundings = if faces.is_empty() {
        0.0
    } else {
        faces.iter().map(|f| f.surroundings_k).sum::<f64>() / faces.len() as f64
    };
    let mut first = setup.clone();
    first.cell_sinks.extend(radiative_sinks(
        domain,
        &faces,
        &vec![surroundings; domain.cell_count()],
    ));
    let mut solution = solve_energy(domain, fluid, solids, flow, &first, energy, gate)?;
    let mut change = f64::INFINITY;
    let mut iterations = 0usize;
    while iterations < radiation.max_iterations {
        poll(gate)?;
        iterations += 1;
        let mut step = setup.clone();
        step.cell_sinks
            .extend(radiative_sinks(domain, &faces, &solution.temperature));
        let next = solve_energy(domain, fluid, solids, flow, &step, energy, gate)?;
        let (lo, hi) = next
            .temperature
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(*t), hi.max(*t))
            });
        let largest = next
            .temperature
            .iter()
            .zip(&solution.temperature)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        change = largest / (hi - lo).max(f64::MIN_POSITIVE);
        solution = next;
        if change <= radiation.tolerance {
            break;
        }
    }
    if change > radiation.tolerance {
        return Err(ChtError::SolverNotConverged {
            system: "radiation",
            iterations,
            relative_residual: change,
            tolerance: radiation.tolerance,
        });
    }
    let report = RadiationReport {
        exposed_faces: faces.len(),
        rays: faces.len() * radiation.rays_per_face,
        iterations,
        radiated_w: radiated_power(domain, &faces, &solution.temperature),
        temperature_change: change,
    };
    Ok((solution, report))
}
