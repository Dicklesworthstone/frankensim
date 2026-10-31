//! Mutual Newtonian gravity on the original rigid-body states and mass properties.
//!
//! Each unordered pair is evaluated once, with equal and opposite COM forces and
//! one potential contribution. All bodies receive simultaneous endpoint loads;
//! advancing a source before evaluating its partner would break this symmetry.
//! This is an O(N²) direct monopole model, with O(N) scratch, not an FMM, softened
//! particle model, finite-body multipole expansion, or contact solver. Native
//! full-inertia rotation is retained but receives no monopole gravity torque.

use super::{GravityError, finite, poll, same_frame};
use crate::{
    DynamicsDiagnostics, DynamicsError, Gravity, MassProperties, RigidBodyIntegrator,
    RigidBodyState, Vec3, Wrench,
};

/// Explicit inertial frame, Newtonian constant, separation domain and work cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MutualGravity {
    frame: u64,
    gravitational_constant: f64,
    minimum_separation_m: f64,
    maximum_pairs: usize,
}

/// One simultaneous evaluation, in exactly the supplied body order.
#[derive(Clone, Debug, PartialEq)]
pub struct MutualGravityLoad {
    /// World COM forces; the reaction on every source is included.
    pub forces_world: Vec<Vec3>,
    /// Sum over i < j of -G m_i m_j / r_ij [J], zero at infinite separation.
    pub potential_energy_j: f64,
    /// Number of unordered interactions evaluated once in this sample.
    pub pair_count: usize,
}

/// A complete closed-system step; no candidate escapes on refusal/cancellation.
#[derive(Clone, Debug, PartialEq)]
pub struct MutualGravityStepReceipt {
    /// Positive fixed duration [s].
    pub duration_seconds: f64,
    /// Original native states in the supplied body order after the complete step.
    pub states_after: Vec<RigidBodyState>,
    /// Whole-system energies and momenta before the step. Potential is counted once.
    pub diagnostics_before: DynamicsDiagnostics,
    /// Whole-system energies and momenta after the step, about the same world origin.
    pub diagnostics_after: DynamicsDiagnostics,
    /// Interactions per pass: two force samples and one drift-domain check.
    pub pair_count: usize,
}

impl MutualGravity {
    /// Construct an explicit SI model; G has units m³ kg⁻¹ s⁻².
    ///
    /// G must be positive and finite. Every pair must stay strictly farther than
    /// the finite nonnegative separation cutoff. Zero selects the point-source
    /// singularity only; a positive cutoff can declare a resolved near-field limit.
    /// `maximum_pairs` bounds N(N-1)/2 before any allocation or pair work. Empty
    /// and one-body systems have no interactions and admit a zero pair budget.
    pub fn new(
        frame: u64,
        gravitational_constant: f64,
        minimum_separation_m: f64,
        maximum_pairs: usize,
    ) -> Result<Self, GravityError> {
        if !gravitational_constant.is_finite()
            || gravitational_constant <= 0.0
            || !minimum_separation_m.is_finite()
            || minimum_separation_m < 0.0
        {
            return Err(GravityError::Invalid(
                "mutual gravity needs positive G and nonnegative separation cutoff",
            ));
        }
        Ok(Self {
            frame,
            gravitational_constant,
            minimum_separation_m,
            maximum_pairs,
        })
    }

    /// Evaluate conservative pair forces and potential with no inferred source.
    ///
    /// States and properties must have equal lengths and matching indices. Input
    /// order fixes the reduction order; compensated sums reduce cancellation, but
    /// permuting bodies is not a claim of bit-identical floating-point results.
    /// Cancellation is checked per body and per pair, and before publication.
    pub fn load(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MutualGravityLoad, GravityError> {
        poll(&mut cancelled)?;
        same_frame(frame, self.frame)?;
        if states.len() != properties.len() {
            return Err(GravityError::Invalid(
                "mutual gravity body/property count mismatch",
            ));
        }
        let count = states.len();
        // Divide an even factor first, avoiding a spurious intermediate overflow.
        let previous = count.saturating_sub(1);
        let pair_count = if count % 2 == 0 {
            (count / 2).checked_mul(previous)
        } else {
            count.checked_mul(previous / 2)
        }
        .filter(|&pairs| pairs <= self.maximum_pairs)
        .ok_or(GravityError::Invalid("mutual gravity pair budget exceeded"))?;
        let mut force_sums = storage(count)?;
        for (&state, &mass) in states.iter().zip(properties) {
            poll(&mut cancelled)?;
            state.validate()?;
            mass.validate()?;
            force_sums.push(VectorSum::default());
        }
        let mut potential = Sum::default();
        for i in 0..count {
            for j in i + 1..count {
                poll(&mut cancelled)?;
                let displacement = states[j]
                    .pose()
                    .position_world()
                    .sub(states[i].pose().position_world());
                let radius = self.admit_separation(displacement)?;
                let (magnitude, force_magnitude) = pair_magnitudes(
                    self.gravitational_constant,
                    properties[i].mass(),
                    properties[j].mass(),
                    radius,
                )?;
                let force = Vec3::new(
                    force_component(force_magnitude, displacement.x, radius)?,
                    force_component(force_magnitude, displacement.y, radius)?,
                    force_component(force_magnitude, displacement.z, radius)?,
                );
                force_sums[i].add(force)?;
                force_sums[j].add(force.scale(-1.0))?;
                potential.add(-magnitude)?;
            }
        }
        let mut forces_world = storage(count)?;
        for sum in force_sums {
            poll(&mut cancelled)?;
            forces_world.push(sum.value()?);
        }
        let potential_energy_j = potential.value()?;
        poll(&mut cancelled)?;
        Ok(MutualGravityLoad {
            forces_world,
            potential_energy_j,
            pair_count,
        })
    }

    /// Advance all closed-system bodies by the same second-order fixed step.
    ///
    /// Half-kicks at the common initial configuration surround the original
    /// native force-free translation/rotation map and a second simultaneous
    /// pair sample. The final half-kicks use those endpoint forces. No source is
    /// held fixed, and no energy or momentum is projected after integration.
    ///
    /// The conservative smooth map is symplectic and time reversible up to
    /// roundoff. Energy has a second-order timestep defect, not exact conservation;
    /// bounded long-time behavior requires adequately resolved bounded motion.
    /// Linear and orbital-plus-spin angular momenta are conserved to roundoff.
    ///
    /// Both endpoints AND the complete straight relative drift segments must
    /// satisfy the separation cutoff. This prevents tunneling through a source
    /// during the split drift; it is not an enclosure of the true curved orbit.
    /// Close encounters still require caller-owned timestep refinement. Inputs
    /// are immutable; failures, exhausted budgets and late cancellation return
    /// no partly advanced system. External loads, collisions, variable masses,
    /// adaptivity and finite-body gravity are outside this closed-system map.
    pub fn step(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        duration_seconds: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MutualGravityStepReceipt, GravityError> {
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration.into());
        }
        let before = self.load(frame, states, properties, &mut cancelled)?;
        let diagnostics_before = diagnostics(
            states,
            properties,
            before.potential_energy_j,
            &mut cancelled,
        )?;
        let native = RigidBodyIntegrator::new(Gravity::ZERO);
        let half = 0.5 * duration_seconds;
        let mut candidates = storage(states.len())?;
        for ((&state, &mass), &force) in states.iter().zip(properties).zip(&before.forces_world) {
            poll(&mut cancelled)?;
            let kicked = kick(state, force, half)?;
            candidates.push(
                native
                    .step_symplectic(kicked, mass, Wrench::ZERO, duration_seconds)?
                    .state_after,
            );
        }
        for i in 0..states.len() {
            for j in i + 1..states.len() {
                poll(&mut cancelled)?;
                let start = states[j]
                    .pose()
                    .position_world()
                    .sub(states[i].pose().position_world());
                let end = candidates[j]
                    .pose()
                    .position_world()
                    .sub(candidates[i].pose().position_world());
                let change = end.sub(start);
                let length = change.stable_norm("mutual gravity drift displacement")?;
                if length > 0.0 {
                    let direction =
                        Vec3::new(change.x / length, change.y / length, change.z / length);
                    let fraction = (-finite(start.dot(direction))? / length).clamp(0.0, 1.0);
                    if fraction > 0.0 && fraction < 1.0 {
                        // Interpolating two opposite endpoints can leave a
                        // rounding residual at an exact source crossing. Use
                        // |start × end| / |end-start| for an interior minimum.
                        self.admit_radius(segment_line_distance(start, end, length)?)?;
                    }
                }
            }
        }
        let after = self.load(frame, &candidates, properties, &mut cancelled)?;
        for (state, &force) in candidates.iter_mut().zip(&after.forces_world) {
            poll(&mut cancelled)?;
            *state = kick(*state, force, half)?;
        }
        let diagnostics_after = diagnostics(
            &candidates,
            properties,
            after.potential_energy_j,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MutualGravityStepReceipt {
            duration_seconds,
            states_after: candidates,
            diagnostics_before,
            diagnostics_after,
            pair_count: before.pair_count,
        })
    }

    fn admit_separation(self, displacement: Vec3) -> Result<f64, GravityError> {
        let radius = displacement.stable_norm("mutual gravity separation")?;
        self.admit_radius(radius)
    }

    fn admit_radius(self, radius: f64) -> Result<f64, GravityError> {
        if radius <= self.minimum_separation_m {
            Err(GravityError::Invalid(
                "mutual gravity separation cutoff reached",
            ))
        } else {
            Ok(radius)
        }
    }
}

fn kick(state: RigidBodyState, force: Vec3, duration: f64) -> Result<RigidBodyState, GravityError> {
    Ok(RigidBodyState::new(
        state.pose(),
        state.linear_momentum_world().add(force.scale(duration)),
        state.angular_momentum_body(),
    )?)
}

fn diagnostics(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    potential: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DynamicsDiagnostics, GravityError> {
    let native = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut translation = Sum::default();
    let mut rotation = Sum::default();
    let mut momentum = VectorSum::default();
    let mut angular = VectorSum::default();
    for (&state, &mass) in states.iter().zip(properties) {
        poll(cancelled)?;
        let values = native.diagnostics(state, mass)?;
        translation.add(values.translational_kinetic_energy)?;
        rotation.add(values.rotational_kinetic_energy)?;
        momentum.add(values.linear_momentum_world)?;
        angular.add(values.angular_momentum_world)?;
    }
    let translational_kinetic_energy = translation.value()?;
    let rotational_kinetic_energy = rotation.value()?;
    Ok(DynamicsDiagnostics {
        translational_kinetic_energy,
        rotational_kinetic_energy,
        gravitational_potential_energy: potential,
        mechanical_energy: finite(
            translational_kinetic_energy + rotational_kinetic_energy + potential,
        )?,
        linear_momentum_world: momentum.value()?,
        angular_momentum_world: angular.value()?,
    })
}

fn positive(value: f64) -> Result<f64, GravityError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(GravityError::Invalid(
            "unrepresentable mutual gravity pair arithmetic",
        ))
    }
}

// Evaluate G*m1*m2/r and G*m1*m2/r² with bounded mantissas and integer binary
// exponents. Forming G*m1 first can underflow/overflow despite representable
// answers, and even change admission when the two body indices are exchanged.
// Sorting the masses fixes multiplication order for either pair enumeration.
fn pair_magnitudes(g: f64, m1: f64, m2: f64, radius: f64) -> Result<(f64, f64), GravityError> {
    let (g, eg) = binary_parts(g);
    let (small, es) = binary_parts(m1.min(m2));
    let (large, el) = binary_parts(m1.max(m2));
    let (r, er) = binary_parts(radius);
    let potential = ((g * small) * large) / r;
    let exponent = eg + es + el - er;
    Ok((
        scaled_positive(potential, exponent)?,
        scaled_positive(potential / r, exponent - er)?,
    ))
}

// Divide only after rescaling the product. A tiny direction component may
// underflow although its product with a large force is still representable.
fn force_component(magnitude: f64, displacement: f64, radius: f64) -> Result<f64, GravityError> {
    let (product, exponent) = product_parts(magnitude, displacement);
    let (r, er) = binary_parts(radius);
    Ok(scaled_nonnegative(product.abs() / r, exponent - er)?.copysign(displacement))
}

fn segment_line_distance(start: Vec3, end: Vec3, length: f64) -> Result<f64, GravityError> {
    let a = [start.x, start.y, start.z];
    let b = [end.x, end.y, end.z];
    let (denominator, denominator_exponent) = binary_parts(length);
    let mut perpendicular = [0.0; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let (first, e1) = product_parts(a[j], b[k]);
        let (second, e2) = product_parts(a[k], b[j]);
        let exponent = if first == 0.0 {
            e2
        } else if second == 0.0 {
            e1
        } else {
            e1.max(e2)
        };
        let first = scaled_nonnegative(first.abs(), e1 - exponent)?.copysign(first);
        let second = scaled_nonnegative(second.abs(), e2 - exponent)?.copysign(second);
        // Ordinary symmetric subtraction makes equal products cancel exactly;
        // a one-sided FMA would retain only one product's rounding error.
        perpendicular[i] = scaled_nonnegative(
            (first - second).abs() / denominator,
            exponent - denominator_exponent,
        )?;
    }
    Ok(
        Vec3::new(perpendicular[0], perpendicular[1], perpendicular[2])
            .stable_norm("mutual gravity drift closest separation")?,
    )
}

fn product_parts(first: f64, second: f64) -> (f64, i32) {
    if first == 0.0 || second == 0.0 {
        return (0.0, 0);
    }
    let (a, ea) = binary_parts(first.abs());
    let (b, eb) = binary_parts(second.abs());
    let sign = if first.is_sign_negative() == second.is_sign_negative() {
        1.0
    } else {
        -1.0
    };
    ((a * b).copysign(sign), ea + eb)
}

// Positive finite values only; exact powers of two also normalize subnormals.
fn binary_parts(value: f64) -> (f64, i32) {
    let (normal, shift) = if value < f64::MIN_POSITIVE {
        (value * 4_503_599_627_370_496.0, -52)
    } else {
        (value, 0)
    };
    let bits = normal.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + shift;
    (
        f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52)),
        exponent,
    )
}

fn scaled_positive(mantissa: f64, exponent: i32) -> Result<f64, GravityError> {
    positive(scaled_nonnegative(mantissa, exponent)?)
}

fn scaled_nonnegative(mantissa: f64, exponent: i32) -> Result<f64, GravityError> {
    if mantissa == 0.0 {
        return Ok(0.0);
    }
    let (unit, adjustment) = binary_parts(mantissa);
    let exponent = exponent + adjustment;
    if exponent > 1023 {
        return positive(0.0);
    }
    if exponent < -1075 {
        return Ok(0.0);
    }
    if exponent >= -1022 {
        let scale = f64::from_bits(((exponent + 1023) as u64) << 52);
        finite(unit * scale)
    } else {
        // Round just once into the subnormal range, including values above half
        // the smallest subnormal. Scalar pair magnitudes refuse a rounded zero;
        // individual transverse components may be genuinely unrepresentable.
        let scale = f64::from_bits(((exponent + 1074 + 1023) as u64) << 52);
        finite((unit * scale) * f64::from_bits(1))
    }
}

fn storage<T>(count: usize) -> Result<Vec<T>, GravityError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| GravityError::Invalid("mutual gravity allocation refused"))?;
    Ok(values)
}

// Neumaier reductions in fixed body/pair order retain small contributions when
// large opposing forces or momenta nearly cancel. They do not enforce invariants.
#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), GravityError> {
        let next = finite(self.sum + value)?;
        self.correction = finite(
            self.correction
                + if self.sum.abs() >= value.abs() {
                    (self.sum - next) + value
                } else {
                    (value - next) + self.sum
                },
        )?;
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64, GravityError> {
        finite(self.sum + self.correction)
    }
}
#[derive(Clone, Copy, Default)]
struct VectorSum([Sum; 3]);
impl VectorSum {
    fn add(&mut self, value: Vec3) -> Result<(), GravityError> {
        for (sum, component) in self.0.iter_mut().zip([value.x, value.y, value.z]) {
            sum.add(component)?;
        }
        Ok(())
    }
    fn value(self) -> Result<Vec3, GravityError> {
        Ok(Vec3::new(
            self.0[0].value()?,
            self.0[1].value()?,
            self.0[2].value()?,
        ))
    }
}
