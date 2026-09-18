//! Distributed unilateral contact as a nameless obstacle.
//!
//! `fs-dcontact` owns the power-law potential. A fretboard, a reed
//! lay, a snare, and a cable against a stay are fillings.

use fs_dcontact::{ContactStorage, DContactError, Obstacle, string_collocation};
use fs_math::det;
use fs_phs::Storage;
use fs_scenario::{PrestressedString, UnilateralObstacle};

fn validate_span_stations(spec: &UnilateralObstacle) -> Result<(), DContactError> {
    if spec.stations.len() != spec.gaps_m.len() || spec.stations.is_empty() {
        return Err(DContactError::Shape {
            what: "obstacle stations vs gaps",
        });
    }
    if spec
        .stations
        .iter()
        .any(|&s| !s.is_finite() || s <= 0.0 || s >= 1.0)
    {
        return Err(DContactError::Parameter {
            what: "obstacle stations must be finite and strictly inside (0, 1)",
        });
    }
    Ok(())
}

/// Build a taut-span obstacle for sine-mode collocation.
/// Stations must lie strictly inside the span; they are never relocated.
///
/// # Errors
/// Station/gap mismatch or dcontact admission.
pub fn span_obstacle(
    string: &PrestressedString,
    spec: &UnilateralObstacle,
) -> Result<Obstacle, DContactError> {
    validate_span_stations(spec)?;
    let points: Vec<f64> = spec.stations.iter().map(|&s| s * string.length_m).collect();
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

/// Frozen-opening single-slit response: total force is
/// `max(elastic - damping * opening_velocity, 0)`. These are the
/// existing obstacle's elastic reaction and Hunt–Crossley coefficient,
/// for the unit opening coordinate used by `slit_lay`.
pub(crate) fn slit_contact_coefficients(
    obstacle: &Obstacle,
    opening_m: f64,
) -> Result<(f64, f64), DContactError> {
    if obstacle.n_points() != 1 || obstacle.collocation() != [-1.0] {
        return Err(DContactError::Shape {
            what: "slit response requires one unit opening coordinate",
        });
    }
    let elastic = slit_contact_force(obstacle, opening_m)?;
    let damping = elastic * obstacle.internal_loss();
    if !elastic.is_finite() || elastic < 0.0 || !damping.is_finite() || damping < 0.0 {
        return Err(DContactError::Parameter {
            what: "slit response coefficients must be finite and nonnegative",
        });
    }
    Ok((elastic, damping))
}

/// Modal contact forces `f_k = −∂V/∂q_k` for an interleaved `[q, p]` state.
///
/// # Errors
/// Storage shape.
pub fn modal_contact_forces(
    string: &PrestressedString,
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
    add_modal_friction_forces(string.n_modes, obstacles, &obs, x, &mut forces);
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
    string: &PrestressedString,
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
/// Obstacle shape or dcontact admission.
pub fn modal_friction_forces(
    string: &PrestressedString,
    obstacles: &[UnilateralObstacle],
    x: &[f64],
) -> Result<Vec<f64>, DContactError> {
    let mut forces = vec![0.0; string.n_modes];
    if obstacles.iter().all(|o| o.mu_kinetic == 0.0) {
        for obstacle in obstacles {
            validate_span_stations(obstacle)?;
        }
        return Ok(forces);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    add_modal_friction_forces(string.n_modes, obstacles, &obs?, x, &mut forces);
    Ok(forces)
}

fn add_modal_friction_forces(
    n_modes: usize,
    specs: &[UnilateralObstacle],
    obstacles: &[Obstacle],
    x: &[f64],
    forces: &mut [f64],
) {
    for (spec, obstacle) in specs.iter().zip(obstacles) {
        if !(spec.mu_kinetic > 0.0) {
            continue;
        }
        let law = fs_tribo::FrictionLaw::Coulomb {
            static_mu: spec.mu_kinetic,
            kinetic_mu: spec.mu_kinetic,
        };
        for i in 0..obstacle.n_points() {
            let row = &obstacle.collocation()[i * n_modes..(i + 1) * n_modes];
            let mut y = 0.0;
            let mut v = 0.0;
            for (k, &ph) in row.iter().enumerate() {
                y += x[2 * k] * ph;
                v += x[2 * k + 1] * ph;
            }
            let gap = obstacle.gaps()[i];
            let pen = (y - gap).max(0.0);
            if pen <= 0.0 {
                continue;
            }
            let normal =
                obstacle.stiffness() * det::pow(pen, obstacle.alpha()) * obstacle.weights()[i];
            // The obstacle is stationary; the string's relative velocity is v.
            // fs-tribo already returns the force opposing its velocity argument.
            let ft = law
                .regularized_traction_1d(v, normal, 1.0e-3)
                .unwrap_or(0.0);
            for (force, &ph) in forces.iter_mut().zip(row) {
                *force += ft * ph;
            }
        }
    }
}

/// Hunt–Crossley port forces only. Conservative contact lives in
/// [`wrap_modal_contact`].
///
/// # Errors
/// Obstacle shape or dcontact admission.
pub fn modal_hunt_crossley_forces(
    string: &PrestressedString,
    obstacles: &[UnilateralObstacle],
    x: &[f64],
) -> Result<Vec<f64>, DContactError> {
    let mut forces = vec![0.0; string.n_modes];
    if obstacles.iter().all(|o| o.internal_loss == 0.0) {
        for obstacle in obstacles {
            validate_span_stations(obstacle)?;
        }
        return Ok(forces);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    let obs = obs?;
    let velocities: Vec<f64> = (0..string.n_modes).map(|k| x[2 * k + 1]).collect();
    for extra in obs
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

    fn string() -> PrestressedString {
        PrestressedString {
            length_m: 0.8,
            tension_n: 20.0,
            lin_density_kg_m: 0.005,
            axial_stiffness_n: 0.0,
            width_m: 0.001,
            n_modes: 1,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
        }
    }

    fn obstacle(stations: Vec<f64>) -> UnilateralObstacle {
        UnilateralObstacle {
            gaps_m: vec![0.0; stations.len()],
            stations,
            stiffness: 1.0e6,
            alpha: 1.0,
            mu_kinetic: 0.0,
            internal_loss: 0.0,
            provenance: "analytical string contact fixture".into(),
        }
    }

    #[test]
    fn g1_span_contact_preserves_near_endpoint_geometry() {
        let string = string();
        let q = 0.01;
        for station in [1.0e-8, 0.37, 1.0 - 1.0e-8] {
            let spec = obstacle(vec![station]);
            let contact = span_obstacle(&string, &spec).expect("valid station");
            // Independent one-mode spring reference: V = K (phi q)^2 / 2.
            // Standard sine is the oracle here, not fs-dcontact collocation.
            let phi = (core::f64::consts::PI * station).sin()
                / (string.lin_density_kg_m * string.length_m / 2.0).sqrt();
            let expected = -spec.stiffness * phi * phi * q;
            let actual = modal_contact_forces(&string, &[spec], &[q, 0.0]).expect("force")[0];
            // Near pi, argument rounding is amplified by the small sine.
            assert!((contact.collocation()[0] - phi).abs() <= 5.0e-8 * phi.abs());
            assert!(
                (actual - expected).abs() <= 1.0e-7 * expected.abs(),
                "station {station}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn g3_span_contact_split_uses_identical_geometry() {
        let mut string = string();
        string.n_modes = 2;
        let mut spec = obstacle(vec![1.0e-8, 0.37, 1.0 - 1.0e-8]);
        spec.gaps_m.fill(-0.001);
        spec.mu_kinetic = 0.4;
        spec.internal_loss = 0.2;
        let specs = [spec];
        let x = [0.0001, 0.002, -0.00002, 0.003];
        let total = modal_contact_forces(&string, &specs, &x).expect("combined");
        let friction = modal_friction_forces(&string, &specs, &x).expect("friction");
        let loss = modal_hunt_crossley_forces(&string, &specs, &x).expect("internal loss");
        let storage = wrap_modal_contact(Box::new(ZeroStorage), &string, &specs)
            .expect("conservative storage");
        let mut gradient = [0.0; 4];
        storage.gradient(&x, &mut gradient);
        for k in 0..string.n_modes {
            let expected = -gradient[2 * k] + loss[k] + friction[k];
            assert!((total[k] - expected).abs() <= 1.0e-12 * expected.abs());
        }
        assert!(friction[0] * x[1] + friction[1] * x[3] < 0.0);
    }

    #[test]
    fn g0_span_contact_rejects_invalid_stations_on_all_paths() {
        let string = string();
        for station in [
            f64::NAN,
            f64::NEG_INFINITY,
            f64::INFINITY,
            -0.01,
            0.0,
            1.0,
            1.01,
        ] {
            let specs = [obstacle(vec![station])];
            assert!(span_obstacle(&string, &specs[0]).is_err());
            assert!(modal_contact_forces(&string, &specs, &[0.0; 2]).is_err());
            // Zero loss/friction does not bypass geometric admission.
            assert!(modal_friction_forces(&string, &specs, &[0.0; 2]).is_err());
            assert!(modal_hunt_crossley_forces(&string, &specs, &[0.0; 2]).is_err());
            assert!(wrap_modal_contact(Box::new(ZeroStorage), &string, &specs).is_err());
        }
    }

    #[test]
    fn slit_lay_pushes_a_penetrating_opening_back() {
        let lay = slit_lay(1.0e6, 2.0).expect("lay");
        let closed = slit_contact_force(&lay, 0.0).expect("force");
        let into = slit_contact_force(&lay, -1.0e-4).expect("force");
        assert!(closed.abs() < 1.0e-12);
        assert!(into > 0.0);
    }
}
