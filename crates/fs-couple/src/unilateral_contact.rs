//! Distributed unilateral contact as a nameless obstacle.
//!
//! `fs-dcontact` owns the power-law potential. A fretboard, a reed
//! lay, a snare, and a cable against a stay are fillings.

use fs_dcontact::{ContactStorage, DContactError, Obstacle, string_collocation};
use fs_math::det;
use fs_phs::Storage;
use fs_scenario::{PrestressedString, UnilateralObstacle};

/// Build a taut-span obstacle for sine-mode collocation.
///
/// # Errors
/// Station/gap mismatch or dcontact admission.
pub fn span_obstacle(
    string: PrestressedString,
    spec: &UnilateralObstacle,
) -> Result<Obstacle, DContactError> {
    if spec.stations.len() != spec.gaps_m.len() || spec.stations.is_empty() {
        return Err(DContactError::Shape {
            what: "obstacle stations vs gaps",
        });
    }
    let points: Vec<f64> = spec
        .stations
        .iter()
        .map(|&s| s.clamp(1.0e-6, 1.0 - 1.0e-6) * string.length_m)
        .collect();
    let phi = string_collocation(
        string.length_m,
        string.lin_density_kg_m,
        &points,
        string.n_modes,
    )?;
    let n = spec.stations.len();
    Obstacle::new(
        phi,
        n,
        string.n_modes,
        spec.gaps_m.clone(),
        vec![1.0 / n as f64; n],
        spec.stiffness,
        spec.alpha,
        spec.provenance.clone(),
    )?
    .with_internal_loss(spec.internal_loss)
}

/// Reed/valve lay: one collocation point on a 1-DOF opening.
///
/// # Errors
/// dcontact admission.
pub fn slit_lay(stiffness: f64, alpha: f64) -> Result<Obstacle, DContactError> {
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        stiffness,
        alpha,
        "slit-lay".to_string(),
    )
}

/// Contact force on a 1-DOF opening (`q` is the opening coordinate).
///
/// # Errors
/// Storage shape.
pub fn slit_contact_force(obstacle: &Obstacle, opening_m: f64) -> Result<f64, DContactError> {
    let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![obstacle.clone()])?;
    let x = [opening_m, 0.0];
    let mut g = [0.0, 0.0];
    storage.gradient(&x, &mut g);
    Ok(-g[0])
}

/// Modal contact forces `f_k = −∂V/∂q_k` for an interleaved `[q, p]` state.
///
/// # Errors
/// Storage shape.
pub fn modal_contact_forces(
    string: PrestressedString,
    obstacles: &[UnilateralObstacle],
    x: &[f64],
) -> Result<Vec<f64>, DContactError> {
    if obstacles.is_empty() {
        return Ok(vec![0.0; string.n_modes]);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    let obs = obs?;
    let storage = ContactStorage::new(Box::new(ZeroStorage), string.n_modes, obs.clone())?;
    let mut g = vec![0.0; 2 * string.n_modes];
    storage.gradient(x, &mut g);
    let mut forces: Vec<f64> = (0..string.n_modes).map(|k| -g[2 * k]).collect();
    let velocities: Vec<f64> = (0..string.n_modes).map(|k| x[2 * k + 1]).collect();
    for extra in obs
        .iter()
        .map(|o| o.dissipative_modal_forces(string.n_modes, x, &velocities))
    {
        for (f, e) in forces.iter_mut().zip(extra) {
            *f += e;
        }
    }
    let mass_scale = (string.lin_density_kg_m * string.length_m / 2.0).sqrt();
    let pi = core::f64::consts::PI;
    for spec in obstacles {
        if !(spec.mu_kinetic > 0.0) {
            continue;
        }
        let law = fs_tribo::FrictionLaw::Coulomb {
            static_mu: spec.mu_kinetic,
            kinetic_mu: spec.mu_kinetic,
        };
        for (i, &s) in spec.stations.iter().enumerate() {
            let mut y = 0.0;
            let mut v = 0.0;
            for k in 0..string.n_modes {
                let ph = det::sin((k + 1) as f64 * pi * s) / mass_scale;
                y += x[2 * k] * ph;
                v += x[2 * k + 1] * ph;
            }
            let gap = spec.gaps_m[i];
            let pen = (y - gap).max(0.0);
            if pen <= 0.0 {
                continue;
            }
            let normal = spec.stiffness * det::pow(pen, spec.alpha) / spec.stations.len() as f64;
            let ft = law
                .regularized_traction_1d(-v, normal, 1.0e-3)
                .unwrap_or(0.0);
            #[allow(clippy::needless_range_loop)] // the modal index spans forces and shapes
            for k in 0..string.n_modes {
                let ph = det::sin((k + 1) as f64 * pi * s) / mass_scale;
                forces[k] += ft * ph;
            }
        }
    }
    Ok(forces)
}

/// Wrap modal storage with the conservative contact potential.
///
/// Friction is not a gradient of `H` and stays a port force
/// ([`modal_friction_forces`]). An empty obstacle list returns the
/// inner storage unchanged.
///
/// # Errors
/// Obstacle shape or dcontact admission.
pub fn wrap_modal_contact(
    inner: Box<dyn Storage>,
    string: PrestressedString,
    obstacles: &[UnilateralObstacle],
) -> Result<Box<dyn Storage>, DContactError> {
    if obstacles.is_empty() {
        return Ok(inner);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    Ok(Box::new(ContactStorage::new(inner, string.n_modes, obs?)?))
}

/// Tangential Coulomb traction at contacting stations, as modal forces.
///
/// Conservative contact lives in [`wrap_modal_contact`]. This is only
/// the non-gradient `fs-tribo` remainder.
///
/// # Errors
/// None today; `Result` matches [`modal_contact_forces`].
pub fn modal_friction_forces(
    string: PrestressedString,
    obstacles: &[UnilateralObstacle],
    x: &[f64],
) -> Result<Vec<f64>, DContactError> {
    let mut forces = vec![0.0; string.n_modes];
    let mass_scale = (string.lin_density_kg_m * string.length_m / 2.0).sqrt();
    let pi = core::f64::consts::PI;
    for spec in obstacles {
        if !(spec.mu_kinetic > 0.0) {
            continue;
        }
        let law = fs_tribo::FrictionLaw::Coulomb {
            static_mu: spec.mu_kinetic,
            kinetic_mu: spec.mu_kinetic,
        };
        for (i, &s) in spec.stations.iter().enumerate() {
            let mut y = 0.0;
            let mut v = 0.0;
            for k in 0..string.n_modes {
                let ph = det::sin((k + 1) as f64 * pi * s) / mass_scale;
                y += x[2 * k] * ph;
                v += x[2 * k + 1] * ph;
            }
            let gap = spec.gaps_m[i];
            let pen = (y - gap).max(0.0);
            if pen <= 0.0 {
                continue;
            }
            let normal = spec.stiffness * det::pow(pen, spec.alpha) / spec.stations.len() as f64;
            let ft = law
                .regularized_traction_1d(-v, normal, 1.0e-3)
                .unwrap_or(0.0);
            #[allow(clippy::needless_range_loop)] // the modal index spans forces and shapes
            for k in 0..string.n_modes {
                let ph = det::sin((k + 1) as f64 * pi * s) / mass_scale;
                forces[k] += ft * ph;
            }
        }
    }
    Ok(forces)
}

/// Hunt–Crossley port forces only. Conservative contact lives in
/// [`wrap_modal_contact`].
///
/// # Errors
/// Obstacle shape or dcontact admission.
pub fn modal_hunt_crossley_forces(
    string: PrestressedString,
    obstacles: &[UnilateralObstacle],
    x: &[f64],
) -> Result<Vec<f64>, DContactError> {
    let mut forces = vec![0.0; string.n_modes];
    if obstacles.iter().all(|o| !(o.internal_loss > 0.0)) {
        return Ok(forces);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    let velocities: Vec<f64> = (0..string.n_modes).map(|k| x[2 * k + 1]).collect();
    for extra in obs?
        .iter()
        .map(|o| o.dissipative_modal_forces(string.n_modes, x, &velocities))
    {
        for (f, e) in forces.iter_mut().zip(extra) {
            *f += e;
        }
    }
    Ok(forces)
}

struct ZeroStorage;

impl Storage for ZeroStorage {
    fn hamiltonian(&self, _x: &[f64]) -> f64 {
        0.0
    }

    fn gradient(&self, _x: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slit_lay_pushes_a_penetrating_opening_back() {
        let lay = slit_lay(1.0e6, 2.0).expect("lay");
        let closed = slit_contact_force(&lay, 0.0).expect("force");
        let into = slit_contact_force(&lay, -1.0e-4).expect("force");
        assert!(closed.abs() < 1.0e-12);
        assert!(into > 0.0);
    }
}
