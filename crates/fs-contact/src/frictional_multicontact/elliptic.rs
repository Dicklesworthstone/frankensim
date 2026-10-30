//! Simultaneous directional impulses by exact work-conjugate coordinates.
//!
//! Each physical tangent ellipse is mapped to the original scalar Coulomb
//! disk. With p=B*q and v'=B^T*v, the COMPLETE mobility becomes B^T*W*B.
//! Normal/tangent and intercontact blocks stay coupled in the same original
//! non-associated solver; no normal dilation or sequential impulse clipping
//! is introduced. Original-coordinate law, kinematics and work are admitted
//! independently before returning a candidate.
//!
//! The native determinacy proof is retained explicitly in transformed
//! coordinates. An invertible B preserves uniqueness of impulses or velocities,
//! but its contraction weights and PSD witness vectors are not physical-basis
//! quantities. No timestep, contact pressure or impact-heat allocation is inferred.

use super::*;
use fs_math::det;

/// Caller-declared instantaneous ellipse in the original contact tangent frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticFriction {
    pub coefficients_mu: [f64; 2],
    /// First principal direction; finite nonzero input is normalized.
    pub principal_axis: [f64; 2],
}

/// Original geometry and normal law, with one explicit directional impulse law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticFrictionalContact {
    pub id: u64,
    pub gap: f64,
    pub activation_tolerance: f64,
    pub relative_velocity: [f64; 3],
    pub normal_law: NormalVelocityLaw,
    pub friction: EllipticFriction,
}
impl From<FrictionalContact> for EllipticFrictionalContact {
    fn from(c: FrictionalContact) -> Self {
        Self {
            id: c.id,
            gap: c.gap,
            activation_tolerance: c.activation_tolerance,
            relative_velocity: c.relative_velocity,
            normal_law: c.normal_law,
            friction: EllipticFriction {
                coefficients_mu: [c.friction_coefficient; 2],
                principal_axis: [1.0, 0.0],
            },
        }
    }
}

/// Explicit basis/scaling for the unchanged native solver result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticContactCoordinates {
    /// A in physical p_t=A*q_t and native v'_t=A^T*v_t.
    /// Normal coordinates are unchanged. Equal coefficients use exact identity.
    pub tangent_impulse_map: [[f64; 2]; 2],
    /// Original scalar disk radius is this coefficient times normal impulse.
    pub scalar_friction_coefficient: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EllipticContactResiduals {
    pub normal_impulse_violation: f64,
    pub normal_velocity_violation: f64,
    pub normal_complementarity: f64,
    /// max(mu) times excess of the ORIGINAL principal impulse gauge over p_n.
    pub ellipse_impulse_violation: f64,
    pub sliding_impulse_residual: f64,
    pub tangential_duality_gap: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub coordinate_work_residual_j: f64,
}

/// Physical quantities are in original caller contact order and coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticFrictionalContactResult {
    pub impulses: Vec<[f64; 3]>,
    pub relative_velocities_after: Vec<[f64; 3]>,
    pub normal_targets: Vec<f64>,
    pub active_indices: Vec<usize>,
    pub relative_work: f64,
    pub dissipated_energy: f64,
    pub passivity_residual: f64,
    pub residuals: EllipticContactResiduals,
    /// Unchanged native result, INCLUDING its proof vectors, in q/v' coordinates.
    pub transformed: FrictionalContactResult,
    pub coordinates: Vec<EllipticContactCoordinates>,
    /// Records the conservative velocity tolerance used in native coordinates.
    pub transformed_config: FrictionalContactConfig,
}

#[derive(Debug)]
pub enum EllipticFrictionalContactError {
    Configuration(FrictionalContactError),
    Invalid(&'static str),
    UnsupportedDegenerateEllipse {
        index: usize,
    },
    /// Candidate impulses or determinacy diagnostics inside error use q/v';
    /// the supplied coordinates identify their physical interpretation.
    Transformed {
        error: FrictionalContactError,
        coordinates: Vec<EllipticContactCoordinates>,
    },
    NonFinite,
    ResidualNotSatisfied(EllipticContactResiduals),
    EnergyIncrease {
        work: f64,
        allowed: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EllipticFrictionalContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous elliptic contact: {self:?}")
    }
}
impl std::error::Error for EllipticFrictionalContactError {}

#[derive(Clone, Copy)]
struct Prepared {
    coordinates: EllipticContactCoordinates,
    axis: [f64; 2],
    ratios: [f64; 2],
    inverse_velocity: [[f64; 2]; 2],
    minimum_scale: f64,
}
impl Prepared {
    fn new(
        friction: EllipticFriction,
        index: usize,
    ) -> Result<Self, EllipticFrictionalContactError> {
        use EllipticFrictionalContactError as E;
        let mu = friction.coefficients_mu;
        let axis = friction.principal_axis;
        if mu.into_iter().any(|v| !v.is_finite() || v < 0.0)
            || axis.into_iter().any(|v| !v.is_finite())
        {
            return Err(E::Invalid(
                "finite nonnegative coefficients and texture required",
            ));
        }
        if (mu[0] == 0.0) != (mu[1] == 0.0) {
            return Err(E::UnsupportedDegenerateEllipse { index });
        }
        let scale = axis[0].abs().max(axis[1].abs());
        if scale == 0.0 {
            return Err(E::Invalid("nonzero texture direction required"));
        }
        let normalized = axis.map(|v| v / scale);
        let length = det::hypot(normalized[0], normalized[1]);
        let axis = normalized.map(|v| v / length);
        let reference = mu[0].max(mu[1]);
        let identity = [[1.0, 0.0], [0.0, 1.0]];
        if mu[0] == mu[1] {
            return Ok(Self {
                coordinates: EllipticContactCoordinates {
                    tangent_impulse_map: identity,
                    scalar_friction_coefficient: reference,
                },
                axis,
                ratios: [1.0; 2],
                inverse_velocity: identity,
                minimum_scale: 1.0,
            });
        }
        let ratios = mu.map(|v| v / reference);
        if ratios.into_iter().any(|v| !v.is_finite() || v <= 0.0) {
            return Err(E::Invalid(
                "representable positive directional ratios required",
            ));
        }
        let [c, s] = axis;
        let map = [
            [c * ratios[0], -s * ratios[1]],
            [s * ratios[0], c * ratios[1]],
        ];
        let length2 = c * c + s * s;
        let inverse = [
            [c / (ratios[0] * length2), -s / (ratios[1] * length2)],
            [s / (ratios[0] * length2), c / (ratios[1] * length2)],
        ];
        if inverse.iter().flatten().any(|v| !v.is_finite()) {
            return Err(E::Invalid(
                "representable inverse directional coordinates required",
            ));
        }
        Ok(Self {
            coordinates: EllipticContactCoordinates {
                tangent_impulse_map: map,
                scalar_friction_coefficient: reference,
            },
            axis,
            ratios,
            inverse_velocity: inverse,
            minimum_scale: ratios[0].min(ratios[1]),
        })
    }
    fn block(self) -> [[f64; 3]; 3] {
        let a = self.coordinates.tangent_impulse_map;
        [
            [1.0, 0.0, 0.0],
            [0.0, a[0][0], a[0][1]],
            [0.0, a[1][0], a[1][1]],
        ]
    }
    fn impulse(self, q: [f64; 3]) -> [f64; 3] {
        let a = self.coordinates.tangent_impulse_map;
        [
            q[0],
            a[0][0] * q[1] + a[0][1] * q[2],
            a[1][0] * q[1] + a[1][1] * q[2],
        ]
    }
    fn velocity(self, v: [f64; 3]) -> [f64; 3] {
        let a = self.coordinates.tangent_impulse_map;
        [
            v[0],
            a[0][0] * v[1] + a[1][0] * v[2],
            a[0][1] * v[1] + a[1][1] * v[2],
        ]
    }
}

/// Reuse the original coupled contact solve on one complete work-preserving map.
///
/// Positive principal coefficients and the exact both-zero frictionless limit
/// are admitted. Singly-zero ellipses refuse explicitly; no epsilon traction is
/// introduced. Equal coefficients use the original scalar coordinates exactly.
/// The complete row-major 3N-by-3N mobility is transformed without dropping
/// normal/shear or cross-contact blocks. Original work/evaluation/contact
/// budgets remain bounded; coordinate transforms poll at each contact row.
///
/// Native velocity tolerance is tightened by the smallest tangent scale, so a
/// transformed near-zero slip cannot conceal resolved original physical slip.
/// Impulse maps have operator norm at most one; work is invariant. Independent
/// original-coordinate residuals still decide final physical acceptance.
pub fn solve_elliptic_frictional_contacts(
    contacts: &[EllipticFrictionalContact],
    delassus: &[f64],
    config: FrictionalContactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EllipticFrictionalContactResult, EllipticFrictionalContactError> {
    use EllipticFrictionalContactError as E;
    poll(&mut cancelled)?;
    super::validate_config(config).map_err(E::Configuration)?;
    let n = contacts.len();
    if n > config.maximum_contacts {
        return Err(E::Invalid("bounded contact count required"));
    }
    let d = 3 * n;
    if delassus.len() != d * d || delassus.iter().any(|v| !v.is_finite()) {
        return Err(E::Invalid(
            "bounded contact count and finite complete mobility required",
        ));
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
    {
        return Err(E::Invalid("duplicate contact ID"));
    }
    let mut prepared = Vec::with_capacity(n);
    let mut native_contacts = Vec::with_capacity(n);
    let mut minimum_scale = 1.0_f64;
    for (i, c) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        let p = Prepared::new(c.friction, i)?;
        minimum_scale = minimum_scale.min(p.minimum_scale);
        native_contacts.push(FrictionalContact {
            id: c.id,
            gap: c.gap,
            activation_tolerance: c.activation_tolerance,
            relative_velocity: p.velocity(c.relative_velocity),
            normal_law: c.normal_law,
            friction_coefficient: p.coordinates.scalar_friction_coefficient,
        });
        prepared.push(p);
    }
    let coordinates: Vec<_> = prepared.iter().map(|p| p.coordinates).collect();
    let mut transformed_config = config;
    transformed_config.velocity_tolerance *= minimum_scale;
    if !transformed_config.velocity_tolerance.is_finite()
        || transformed_config.velocity_tolerance <= 0.0
    {
        return Err(E::Invalid(
            "representable transformed velocity tolerance required",
        ));
    }
    let identity = [[1.0, 0.0], [0.0, 1.0]];
    let matrix = if coordinates
        .iter()
        .all(|c| c.tangent_impulse_map == identity)
    {
        delassus.to_vec()
    } else {
        congruence(delassus, &prepared, &order, &mut cancelled)?
    };
    let transformed = match solve_frictional_contacts(
        &native_contacts,
        &matrix,
        transformed_config,
        &mut cancelled,
    ) {
        Ok(result) => result,
        Err(FrictionalContactError::Cancelled) => return Err(E::Cancelled),
        Err(error) => return Err(E::Transformed { error, coordinates }),
    };
    poll(&mut cancelled)?;
    let impulses: Vec<_> = prepared
        .iter()
        .zip(&transformed.impulses)
        .map(|(p, &q)| p.impulse(q))
        .collect();
    let mut velocities: Vec<_> = contacts.iter().map(|c| c.relative_velocity).collect();
    for &i in &order {
        poll(&mut cancelled)?;
        for a in 0..3 {
            for &j in &transformed.active_indices {
                for b in 0..3 {
                    velocities[i][a] += delassus[(3 * i + a) * d + 3 * j + b] * impulses[j][b];
                }
            }
        }
    }
    if impulses
        .iter()
        .chain(&velocities)
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err(E::NonFinite);
    }
    let mut residuals = EllipticContactResiduals::default();
    for &i in &order {
        poll(&mut cancelled)?;
        let qv = transformed.relative_velocities_after[i];
        let inv = prepared[i].inverse_velocity;
        let predicted = [
            qv[0],
            inv[0][0] * qv[1] + inv[0][1] * qv[2],
            inv[1][0] * qv[1] + inv[1][1] * qv[2],
        ];
        let diff = [0, 1, 2].map(|a| velocities[i][a] - predicted[a]);
        residuals.maximum_velocity_residual_m_s =
            residuals.maximum_velocity_residual_m_s.max(norm3(diff)?);
    }
    for &i in &transformed.active_indices {
        poll(&mut cancelled)?;
        physical_law(
            contacts[i],
            prepared[i],
            impulses[i],
            velocities[i],
            transformed.normal_targets[i],
            config,
            &mut residuals,
        )?;
    }
    let mut work = 0.0;
    let mut work_scale = 0.0;
    for &i in &order {
        for a in 0..3 {
            let term =
                impulses[i][a] * (0.5 * contacts[i].relative_velocity[a] + 0.5 * velocities[i][a]);
            work += term;
            work_scale += term.abs();
        }
    }
    checked(work)?;
    checked(work_scale)?;
    residuals.coordinate_work_residual_j = checked(work - transformed.relative_work)?;
    let work_allowance = checked(config.energy_tolerance + 128.0 * f64::EPSILON * work_scale)?;
    if residuals.normal_impulse_violation > config.impulse_tolerance
        || residuals.normal_velocity_violation > config.velocity_tolerance
        || residuals.normal_complementarity > config.complementarity_tolerance
        || residuals.ellipse_impulse_violation > config.impulse_tolerance
        || residuals.sliding_impulse_residual > config.impulse_tolerance
        || residuals.tangential_duality_gap > config.complementarity_tolerance
        || residuals.maximum_velocity_residual_m_s > config.velocity_tolerance
        || residuals.coordinate_work_residual_j.abs() > work_allowance
    {
        return Err(E::ResidualNotSatisfied(residuals));
    }
    let allowed = checked(config.energy_tolerance + 64.0 * f64::EPSILON * work_scale)?;
    if work > allowed {
        return Err(E::EnergyIncrease { work, allowed });
    }
    poll(&mut cancelled)?;
    Ok(EllipticFrictionalContactResult {
        impulses,
        relative_velocities_after: velocities,
        normal_targets: transformed.normal_targets.clone(),
        active_indices: transformed.active_indices.clone(),
        relative_work: work,
        dissipated_energy: (-work).max(0.0),
        passivity_residual: work.max(0.0),
        residuals,
        transformed,
        coordinates,
        transformed_config,
    })
}

fn congruence(
    original: &[f64],
    prepared: &[Prepared],
    order: &[usize],
    cancel: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, EllipticFrictionalContactError> {
    let d = 3 * prepared.len();
    let symmetric = (0..d).all(|i| (0..i).all(|j| original[i * d + j] == original[j * d + i]));
    let mut result = vec![0.0; d * d];
    // Canonical contact order makes the symmetric rounding choice invariant
    // under caller permutations. Nonsymmetric input is never symmetrized.
    for (rank_i, &i) in order.iter().enumerate() {
        poll(cancel)?;
        let left = prepared[i].block();
        for a in 0..3 {
            for (rank_j, &j) in order.iter().enumerate() {
                let right = prepared[j].block();
                for b in 0..3 {
                    let row = 3 * i + a;
                    let col = 3 * j + b;
                    if symmetric && 3 * rank_j + b < 3 * rank_i + a {
                        result[row * d + col] = result[col * d + row];
                        continue;
                    }
                    let mut sum = 0.0;
                    for x in 0..3 {
                        for y in 0..3 {
                            if left[x][a] != 0.0 && right[y][b] != 0.0 {
                                sum += left[x][a]
                                    * original[(3 * i + x) * d + 3 * j + y]
                                    * right[y][b];
                            }
                        }
                    }
                    result[row * d + col] = checked(sum)?;
                }
            }
        }
    }
    Ok(result)
}

fn physical_law(
    contact: EllipticFrictionalContact,
    prepared: Prepared,
    impulse: [f64; 3],
    velocity: [f64; 3],
    target: f64,
    config: FrictionalContactConfig,
    r: &mut EllipticContactResiduals,
) -> Result<(), EllipticFrictionalContactError> {
    let pn = impulse[0];
    let g = checked(velocity[0] - target)?;
    let [c, s] = prepared.axis;
    let p = [
        c * impulse[1] + s * impulse[2],
        -s * impulse[1] + c * impulse[2],
    ];
    let v = [
        c * velocity[1] + s * velocity[2],
        -s * velocity[1] + c * velocity[2],
    ];
    let mu = contact.friction.coefficients_mu;
    let reference = prepared.coordinates.scalar_friction_coefficient;
    let weighted = [prepared.ratios[0] * v[0], prepared.ratios[1] * v[1]];
    let weighted_speed = norm2(weighted)?;
    let speed = norm2([velocity[1], velocity[2]])?;
    let (violation, support_work, expected) = if reference == 0.0 {
        (norm2(p)?, 0.0, [0.0; 2])
    } else {
        let gauge = norm2([p[0] / mu[0], p[1] / mu[1]])?;
        let violation = checked(reference * (gauge - pn).max(0.0))?;
        let support_work = checked((reference * pn) * weighted_speed)?;
        let expected_material = if weighted_speed == 0.0 {
            [0.0; 2]
        } else {
            [
                -mu[0] * pn * (weighted[0] / weighted_speed),
                -mu[1] * pn * (weighted[1] / weighted_speed),
            ]
        };
        let expected = [
            c * expected_material[0] - s * expected_material[1],
            s * expected_material[0] + c * expected_material[1],
        ];
        (violation, support_work, expected)
    };
    r.normal_impulse_violation = r.normal_impulse_violation.max((-pn).max(0.0));
    r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
    r.normal_complementarity = r.normal_complementarity.max(checked((pn * g).abs())?);
    r.ellipse_impulse_violation = r.ellipse_impulse_violation.max(violation);
    let duality =
        checked((support_work + impulse[1] * velocity[1] + impulse[2] * velocity[2]).abs())?;
    r.tangential_duality_gap = r.tangential_duality_gap.max(duality);
    if speed > config.velocity_tolerance {
        r.sliding_impulse_residual = r
            .sliding_impulse_residual
            .max(norm2([impulse[1] - expected[0], impulse[2] - expected[1]])?);
    }
    Ok(())
}
fn checked(value: f64) -> Result<f64, EllipticFrictionalContactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EllipticFrictionalContactError::NonFinite)
    }
}
fn norm2(v: [f64; 2]) -> Result<f64, EllipticFrictionalContactError> {
    checked(det::hypot(v[0], v[1]))
}
fn norm3(v: [f64; 3]) -> Result<f64, EllipticFrictionalContactError> {
    checked(det::hypot(v[0], det::hypot(v[1], v[2])))
}
fn poll(cancel: &mut impl FnMut() -> bool) -> Result<(), EllipticFrictionalContactError> {
    if cancel() {
        Err(EllipticFrictionalContactError::Cancelled)
    } else {
        Ok(())
    }
}
