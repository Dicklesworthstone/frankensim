//! Distributed unilateral contact as a nameless obstacle.
//!
//! `fs-dcontact` owns the power-law potential. A fretboard, a reed
//! lay, a snare, and a cable against a stay are fillings.

/// Existing card and finite-gap owners lowered to work-conjugate contact storage.
pub mod binding;
pub(crate) mod distributed;

use fs_dcontact::{ContactStorage, DContactError, Obstacle, string_collocation};
use fs_math::det;
use fs_phs::Storage;
use fs_scenario::{PrestressedString, UnilateralObstacle};

fn validate_span_stations(spec: &UnilateralObstacle) -> Result<(), DContactError> {
    if spec.stations.len() != spec.gaps_m.len()
        || spec.stations.is_empty()
        || spec.stations.len() != spec.measures_m.len()
    {
        return Err(DContactError::Shape {
            what: "obstacle stations vs gaps and physical line measures",
        });
    }
    if spec.measures_m.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(DContactError::Parameter {
            what: "line quadrature measures must be positive finite metres",
        });
    }
    friction_law(spec)?;
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
        spec.measures_m.clone(),
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

/// Reed lay that starts to touch the reed before the slit is closed: the same
/// one-point power law as [`slit_lay`], but the penetration is `onset_m - q`, so
/// contact begins at the positive opening `onset_m` (the curved lay of a real
/// mouthpiece) and `onset_m = 0` is [`slit_lay`].
///
/// # Errors
/// A negative or non-finite onset, or dcontact admission.
pub fn slit_lay_with_onset(
    stiffness: f64,
    alpha: f64,
    onset_m: f64,
) -> Result<Obstacle, DContactError> {
    if !onset_m.is_finite() || onset_m < 0.0 {
        return Err(DContactError::Parameter {
            what: "lay onset opening must be finite and nonnegative",
        });
    }
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![-onset_m],
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

/// Finite-step slit reaction from the existing contact potential.
///
/// The conservative coefficient is the negative discrete gradient of the
/// `fs-dcontact` storage, so its work accounts for the change in contact
/// potential, including entry into and exit from contact. The second
/// coefficient is the same reaction times the obstacle's Hunt–Crossley loss.
/// At midpoint velocity `v`, use `max(elastic - damping * v, 0)`; the
/// unloading clamp discards released contact energy rather than attracting
/// the opening. This is a numerical contact adapter, not a new material law.
///
/// # Errors
/// Nonfinite openings, incompatible obstacle coordinates, or nonfinite or
/// negative reaction coefficients. No invalid trial is published as a force.
pub fn slit_contact_discrete_coefficients(
    obstacle: &Obstacle,
    opening_before_m: f64,
    opening_after_m: f64,
) -> Result<(f64, f64), DContactError> {
    SlitContactStep::new(obstacle, opening_before_m)?.coefficients(opening_after_m)
}

/// Retain one contact storage across all trials of an implicit step.
pub(crate) struct SlitContactStep {
    storage: ContactStorage,
    /// Opening at the start of the step, measured from the lay's onset.
    opening_before_m: f64,
    /// Opening at which the lay starts to touch the reed (minus its gap); zero for the original slit lay.
    onset_m: f64,
    internal_loss: f64,
}

impl SlitContactStep {
    pub(crate) fn new(obstacle: &Obstacle, opening_before_m: f64) -> Result<Self, DContactError> {
        if !opening_before_m.is_finite() {
            return Err(DContactError::Parameter {
                what: "discrete slit contact requires finite openings",
            });
        }
        slit_contact_coefficients(obstacle, opening_before_m)?;
        // A reed held in contact (a lay that starts before the slit closes, or a reed pressed shut on
        // the zero-gap lay) takes tiny steps dx at its turning points and while it rests, where the
        // discrete gradient drops its correction below dx = 1e-14 |q|. Two things then break the
        // momentum balance of the implicit reed step: the default endpoint-energy subtraction has
        // force noise eps V / dx (1e-3 of the force at the cited stiffness near that threshold), and a
        // penetration formed as (onset - q) carries a rounding mismatch of one ulp of the onset
        // against dx, which the remainder's affine term turns into 1e-2 of the force. Both are removed
        // by evaluating the same power law about the onset (penetration = -(q - onset), a negation
        // with no rounding, one shared increment) with an inner storage that states its zero
        // remainder, so fs-dcontact's cancellation-stable remainder is used for every lay. A zero
        // gap needs no shift (frankensim-0ll17: the modal slit pressed shut by 2 kPa refused).
        let onset_m = obstacle.gaps().first().map_or(0.0, |gap| -*gap);
        let shifted = if onset_m == 0.0 {
            obstacle.clone()
        } else {
            Obstacle::new(
                vec![-1.0],
                1,
                1,
                vec![0.0],
                obstacle.weights().to_vec(),
                obstacle.stiffness(),
                obstacle.alpha(),
                "slit-lay-about-onset".to_string(),
            )?
            .with_internal_loss(obstacle.internal_loss())?
        };
        Ok(Self {
            storage: ContactStorage::new(Box::new(CompleteZeroStorage), 1, vec![shifted])?,
            opening_before_m: opening_before_m - onset_m,
            onset_m,
            internal_loss: obstacle.internal_loss(),
        })
    }

    pub(crate) fn coefficients(&self, opening_after_m: f64) -> Result<(f64, f64), DContactError> {
        if !opening_after_m.is_finite() {
            return Err(DContactError::Parameter {
                what: "discrete slit contact requires finite openings",
            });
        }
        // Hold the dummy momentum at zero. The one-coordinate Gonzalez
        // gradient is the potential secant (or the derivative at coincidence).
        // Reuse the PHS owner's pinned discretization, not a private power law.
        let gradient = fs_phs::discrete_gradient(
            &self.storage,
            &[self.opening_before_m, 0.0],
            &[opening_after_m - self.onset_m, 0.0],
        );
        let elastic = -gradient[0];
        let damping = elastic * self.internal_loss;
        if !elastic.is_finite() || elastic < 0.0 || !damping.is_finite() || damping < 0.0 {
            return Err(DContactError::Parameter {
                what: "discrete slit response coefficients must be finite and nonnegative",
            });
        }
        Ok((elastic, damping))
    }
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
    validate_modal_state(string.n_modes, x)?;
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
        for (f, e) in forces.iter_mut().zip(extra?) {
            *f += e;
        }
    }
    add_modal_friction_forces(string.n_modes, obstacles, &obs, x, &mut forces)?;
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
    validate_modal_state(string.n_modes, x)?;
    let mut forces = vec![0.0; string.n_modes];
    if obstacles.iter().all(|o| o.mu_kinetic == 0.0) {
        for obstacle in obstacles {
            validate_span_stations(obstacle)?;
        }
        return Ok(forces);
    }
    let obs: Result<Vec<_>, _> = obstacles.iter().map(|o| span_obstacle(string, o)).collect();
    add_modal_friction_forces(string.n_modes, obstacles, &obs?, x, &mut forces)?;
    Ok(forces)
}

fn add_modal_friction_forces(
    n_modes: usize,
    specs: &[UnilateralObstacle],
    obstacles: &[Obstacle],
    x: &[f64],
    forces: &mut [f64],
) -> Result<(), DContactError> {
    if specs.len() != obstacles.len()
        || forces.len() != n_modes
        || n_modes.checked_mul(2).is_none_or(|n| x.len() < n)
        || x.iter().any(|v| !v.is_finite())
    {
        return Err(DContactError::Shape {
            what: "friction state, obstacle and force dimensions",
        });
    }
    for (spec, obstacle) in specs.iter().zip(obstacles) {
        let Some((law, regularization)) = friction_law(spec)? else {
            continue;
        };
        for i in 0..obstacle.n_points() {
            let row = &obstacle.collocation()[i * n_modes..(i + 1) * n_modes];
            let mut y = 0.0;
            let mut v = 0.0;
            for (k, &ph) in row.iter().enumerate() {
                y += x[2 * k] * ph;
                v += x[2 * k + 1] * ph;
            }
            if !y.is_finite() || !v.is_finite() {
                return Err(DContactError::Parameter {
                    what: "friction kinematic projection overflow",
                });
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
                .regularized_traction_1d(v, normal, regularization)
                .map_err(friction_error)?;
            for (force, &ph) in forces.iter_mut().zip(row) {
                *force += ft * ph;
            }
        }
    }
    if forces.iter().any(|f| !f.is_finite()) {
        return Err(DContactError::Parameter {
            what: "friction force projection overflow",
        });
    }
    Ok(())
}

fn friction_error(error: fs_tribo::TriboError) -> DContactError {
    let field = match error {
        fs_tribo::TriboError::InvalidInput { field } => field,
        _ => "tangential constitutive state",
    };
    DContactError::Constitutive {
        owner: "fs-tribo::FrictionLaw",
        field,
    }
}

fn friction_law(
    spec: &UnilateralObstacle,
) -> Result<Option<(fs_tribo::FrictionLaw, f64)>, DContactError> {
    let law = fs_tribo::FrictionLaw::Coulomb {
        static_mu: spec.mu_static,
        kinetic_mu: spec.mu_kinetic,
    };
    law.validate().map_err(friction_error)?;
    if spec.mu_kinetic == 0.0 {
        return Ok(None);
    }
    let regularization = spec
        .friction_regularization_m_s
        .ok_or(DContactError::Constitutive {
            owner: "fs-tribo::FrictionLaw",
            field: "regularization_m_s",
        })?;
    law.regularized_traction_1d(0.0, 0.0, regularization)
        .map_err(friction_error)?;
    Ok(Some((law, regularization)))
}

fn validate_modal_state(n_modes: usize, x: &[f64]) -> Result<(), DContactError> {
    if n_modes.checked_mul(2).is_none_or(|n| x.len() < n) || x.iter().any(|v| !v.is_finite()) {
        return Err(DContactError::Shape {
            what: "finite interleaved modal state prefix",
        });
    }
    Ok(())
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
    validate_modal_state(string.n_modes, x)?;
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
        for (f, e) in forces.iter_mut().zip(extra?) {
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

/// [`ZeroStorage`] that states its midpoint remainder: it stores no energy, so the remainder is
/// exactly zero and every stored-energy contribution is the contact's own.
struct CompleteZeroStorage;

impl Storage for CompleteZeroStorage {
    fn hamiltonian(&self, _x: &[f64]) -> f64 {
        0.0
    }

    fn gradient(&self, _x: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }

    fn midpoint_energy_remainder(&self, _a: &[f64], _b: &[f64], _midpoint: &[f64]) -> Option<f64> {
        Some(0.0)
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
            tension_loss_factor: None,
            bending_loss_factor: None,
        }
    }

    fn obstacle(stations: Vec<f64>) -> UnilateralObstacle {
        UnilateralObstacle {
            gaps_m: vec![0.0; stations.len()],
            measures_m: vec![string().length_m / stations.len() as f64; stations.len()],
            stations,
            stiffness: 1.0e6,
            alpha: 1.0,
            mu_kinetic: 0.0,
            mu_static: 0.4,
            friction_regularization_m_s: Some(1.0e-3),
            internal_loss: 0.0,
            provenance: "analytical string contact fixture".into(),
        }
    }

    // A lay with an onset opening pushes back below that opening with the power
    // law of the penetration `onset - q`, and not above it; zero onset is the
    // original slit lay.
    #[test]
    fn slit_lay_onset_starts_contact_before_the_slit_closes() {
        let (k, alpha, onset) = (7.9e6, 2.0, 2.0e-4);
        let lay = slit_lay_with_onset(k, alpha, onset).unwrap();
        let force = |q: f64| slit_contact_force(&lay, q).unwrap();
        assert_eq!(force(onset), 0.0);
        assert_eq!(force(1.5 * onset), 0.0);
        for q in [1.5e-4, 1.0e-4, 0.0, -1.0e-4] {
            let expect = k * (onset - q).powf(alpha);
            assert!(
                (force(q) - expect).abs() <= 1e-9 * expect,
                "opening {q}: {} against {expect}",
                force(q)
            );
        }
        let plain = slit_lay(k, alpha).unwrap();
        let zero = slit_lay_with_onset(k, alpha, 0.0).unwrap();
        for q in [1.0e-4, 0.0, -1.0e-4, -3.0e-4] {
            assert_eq!(
                slit_contact_force(&plain, q).unwrap(),
                slit_contact_force(&zero, q).unwrap()
            );
        }
        assert!(slit_lay_with_onset(k, alpha, -1.0e-6).is_err());
        assert!(slit_lay_with_onset(k, alpha, f64::NAN).is_err());
    }

    // The reed sits in deep contact (0.048 mm below an onset of 0.2 mm, K 7.9e6) when its velocity
    // passes through zero. The secant must stay the force of the midpoint penetration to O(z^2) for
    // every step, in particular across dx = 1e-14 |q| where the discrete gradient drops its correction;
    // the endpoint-energy subtraction made a 1e-3 jump there and left the implicit reed step without a root.
    #[test]
    fn curved_lay_secant_has_no_jump_at_tiny_steps() {
        let (k, alpha, onset, y) = (7.9e6, 2.0, 2.0e-4, 1.5194044950898812e-4);
        let lay = slit_lay_with_onset(k, alpha, onset).unwrap();
        let scale = k * (onset - y).powf(alpha);
        for dx in [
            1.0e-19, 1.4e-18, 1.5e-18, 1.52e-18, 1.54e-18, 1.6e-18, 1.0e-16, 1.0e-12, 1.0e-8,
            -1.5e-18, -1.54e-18, -1.0e-12,
        ] {
            let (elastic, damping) = slit_contact_discrete_coefficients(&lay, y, y + dx).unwrap();
            let midpoint = k * (onset - (y + 0.5 * dx)).powf(alpha);
            assert!(
                (elastic - midpoint).abs() <= 1e-8 * scale,
                "dx {dx:e}: secant {elastic:e} against midpoint force {midpoint:e}"
            );
            assert_eq!(damping, 0.0, "an elastic lay has no loss term");
        }
    }

    #[test]
    fn g0_obstacle_friction_refusal_propagates() {
        let mut invalid = obstacle(vec![0.5]);
        invalid.mu_static = 0.2;
        invalid.mu_kinetic = 0.4;
        assert!(matches!(
            modal_friction_forces(&string(), &[invalid], &[0.001, 0.1]),
            Err(DContactError::Constitutive {
                owner: "fs-tribo::FrictionLaw",
                field: "kinetic_mu > static_mu"
            })
        ));
        let mut missing = obstacle(vec![0.5]);
        missing.mu_kinetic = 0.4;
        missing.friction_regularization_m_s = None;
        assert!(matches!(
            modal_friction_forces(&string(), &[missing], &[0.001, 0.1]),
            Err(DContactError::Constitutive {
                field: "regularization_m_s",
                ..
            })
        ));
    }

    #[test]
    fn g0_span_obstacle_retains_physical_measures_without_normalization() {
        let mut spec = obstacle(vec![0.2, 0.5, 0.8]);
        spec.measures_m = vec![0.08, 0.24, 0.48];
        let law = span_obstacle(&string(), &spec).unwrap();
        assert_eq!(law.weights(), spec.measures_m);
        let mut bad = spec;
        bad.measures_m.pop();
        assert!(span_obstacle(&string(), &bad).is_err());
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
            let expected = -spec.stiffness * spec.measures_m[0] * phi * phi * q;
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

    #[test]
    fn g1_discrete_slit_reaction_accounts_for_contact_entry_and_release() {
        for alpha in [1.0, 1.5, 2.0, 2.3] {
            let lay = slit_lay(1e8, alpha).unwrap();
            for before in [-1e-4, 0.0, 2e-5] {
                let step = SlitContactStep::new(&lay, before).unwrap();
                for after in [-2e-4, -5e-5, 0.0, 1e-4, before] {
                    let (force, damping) = step.coefficients(after).unwrap();
                    let initial = step.storage.hamiltonian(&[before, 0.0]);
                    let final_energy = step.storage.hamiltonian(&[after, 0.0]);
                    let work = force * (after - before);
                    let scale =
                        (initial.abs() + final_energy.abs() + work.abs()).max(f64::MIN_POSITIVE);
                    assert!((work + final_energy - initial).abs() <= 1e-12 * scale);
                    assert!(force >= 0.0);
                    assert_eq!(damping, 0.0);
                    if before == after {
                        assert_eq!(force, slit_contact_force(&lay, before).unwrap());
                    }
                }
            }
        }
    }

    #[test]
    fn g1_discrete_slit_unloading_and_hunt_crossley_loss_are_passive() {
        let lay = slit_lay(1e8, 2.0).unwrap().with_internal_loss(5.0).unwrap();
        let before = -1e-4;
        let step = SlitContactStep::new(&lay, before).unwrap();
        for after in [-2e-4, -9e-5, 0.0, 1e-4] {
            let velocity = (after - before) / 1e-4;
            let (elastic, damping) = step.coefficients(after).unwrap();
            assert_eq!(damping, 5.0 * elastic);
            let force = (elastic - damping * velocity).max(0.0);
            let released =
                step.storage.hamiltonian(&[before, 0.0]) - step.storage.hamiltonian(&[after, 0.0]);
            let loss = released - force * (after - before);
            assert!(loss >= -1e-12 * released.abs().max(f64::MIN_POSITIVE));
            if velocity > 0.2 {
                assert_eq!(force, 0.0, "unloading must not attract the opening");
            }
        }
    }

    #[test]
    fn g0_discrete_slit_rejects_nonfinite_trials_and_nonunit_coordinates() {
        let lay = slit_lay(1e8, 2.0).unwrap();
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(slit_contact_discrete_coefficients(&lay, invalid, 0.0).is_err());
            assert!(slit_contact_discrete_coefficients(&lay, 0.0, invalid).is_err());
        }
        let other = Obstacle::new(
            vec![-2.0],
            1,
            1,
            vec![0.0],
            vec![1.0],
            1e8,
            2.0,
            "nonunit".into(),
        )
        .unwrap();
        assert!(slit_contact_discrete_coefficients(&other, -1e-4, 0.0).is_err());
    }
}
