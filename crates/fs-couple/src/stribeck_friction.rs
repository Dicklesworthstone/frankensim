//! Coupling-layer 1-D ODE embedding of the fs-tribo Stribeck rung.
//!
//! `fs_tribo::FrictionLaw::evaluate` returns zero traction at rest
//! (it never invents a stick reaction). An explicit modal stepper
//! needs a continuous force, so this type is the same regularized
//! ramp `FrictionLaw::regularized_traction_1d` exposes, with the
//! driven-body sign used by a coupling port (`+` when the driver is
//! faster). A bow, a brake, and a fault are the same law.
//!
//! Evaluation delegates to that shared owner, including its deterministic
//! exponential and input/overflow checks. Coefficients remain caller supplied;
//! this adapter does not admit a sourced interface or invent a stick reaction.

/// Regularized friction coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StribeckFriction {
    /// Static coefficient.
    pub mu_static: f64,
    /// Dynamic coefficient (`<= mu_static`).
    pub mu_dynamic: f64,
    /// Physical Stribeck decay speed [m/s], independent of the ODE ramp.
    pub stiction_m_s: f64,
    /// Explicit ODE ramp speed [m/s]; never a static capture threshold.
    pub regularization_m_s: f64,
}

impl StribeckFriction {
    /// Validate friction coefficients.
    ///
    /// # Errors
    /// Refuses non-finite or negative coefficients and a non-positive
    /// decay or regularization scale (fresh-eyes review bead 9svup: the former
    /// silent-zero fallbacks in [`StribeckFriction::traction`] masked
    /// misconfiguration as physics).
    pub fn try_new(
        mu_static: f64,
        mu_dynamic: f64,
        stiction_m_s: f64,
        regularization_m_s: f64,
    ) -> Result<Self, &'static str> {
        if !mu_static.is_finite() || mu_static < 0.0 {
            return Err("mu_static must be finite and non-negative");
        }
        if !mu_dynamic.is_finite() || mu_dynamic < 0.0 || mu_dynamic > mu_static {
            return Err("mu_dynamic must be finite, non-negative, and <= mu_static");
        }
        if !stiction_m_s.is_finite() || stiction_m_s <= 0.0 {
            return Err("Stribeck decay scale must be finite and strictly positive");
        }
        if !regularization_m_s.is_finite() || regularization_m_s <= 0.0 {
            return Err("regularization scale must be finite and strictly positive");
        }
        Ok(Self {
            mu_static,
            mu_dynamic,
            stiction_m_s,
            regularization_m_s,
        })
    }

    /// Tangential force [N] on the driven body for relative velocity
    /// `v_rel = v_driver − v_driven` and normal `n`.
    ///
    /// # Errors
    /// Refuses invalid coefficients (including direct struct construction),
    /// negative/non-finite normal loads, non-finite velocities, and
    /// unrepresentable results through the shared friction owner.
    pub fn traction(self, v_rel: f64, normal_n: f64) -> Result<f64, &'static str> {
        fs_tribo::FrictionLaw::Stribeck {
            static_mu: self.mu_static,
            kinetic_mu: self.mu_dynamic,
            characteristic_speed: self.stiction_m_s,
            viscous_per_speed: 0.0,
        }
        // fs-tribo takes body-minus-driver velocity. Its opposing traction
        // then acts in the driver's direction, as this coupling port requires.
        .regularized_traction_1d(-v_rel, normal_n, self.regularization_m_s)
        .map_err(|error| match error {
            fs_tribo::TriboError::InvalidInput { field } => field,
            _ => "shared friction evaluation refused",
        })
    }
}

/// A force/velocity pair satisfying a discrete collocated contact equation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitFrictionStep {
    /// Traction on the driven body [N].
    pub traction_n: f64,
    /// Driver-minus-body velocity after the held-force step [m/s].
    pub relative_velocity_m_s: f64,
    /// Signed velocity residual of the constitutive/kinematic solve [m/s].
    pub residual_m_s: f64,
}

/// A contact refusal never substitutes an explicit force or a different law.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplicitFrictionError<E> {
    /// Invalid input or an unrepresentable intermediate value.
    InvalidInput(&'static str),
    /// The caller's constitutive owner refused an evaluation.
    Constitutive(E),
    /// The supplied sliding law would add relative mechanical energy.
    NonDissipative,
    /// A discontinuity or unresolved root prevented a consistent step.
    NoConvergence,
}

/// Solve `v = v_free - mobility * F(v)` for an odd, dissipative sliding law.
///
/// `mobility` is the *discrete* collocated endpoint velocity response to one
/// newton held over this step [m/(N s)], not continuous-time inverse mass.
/// Both quantities must come from the same transition that the caller commits.
/// A string, a brake and a sliding mechanical port use the same scalar solve;
/// the constitutive callback retains ownership of materials and state validity.
///
/// The bracket does not allow traction to reverse the slip direction. Static
/// reaction/capture remains the caller's separate constraint solve. A
/// velocity-weakening law can have multiple roots; this deterministic bracket
/// selects one and does not certify uniqueness or continuous-time passivity.
/// In particular, an interval shrinking onto a constitutive jump is NOT
/// convergence: the velocity residual must also pass. There is no allocation.
///
/// # Errors
/// Invalid response data, constitutive refusal, active traction, or failure to
/// resolve a root to floating-point tolerance within 96 bisections.
pub fn solve_implicit_friction<E>(
    free_relative_velocity_m_s: f64,
    mobility_m_per_n_s: f64,
    mut traction: impl FnMut(f64) -> Result<f64, E>,
) -> Result<ImplicitFrictionStep, ImplicitFrictionError<E>> {
    if !free_relative_velocity_m_s.is_finite() {
        return Err(ImplicitFrictionError::InvalidInput("free_relative_velocity"));
    }
    if !mobility_m_per_n_s.is_finite() || mobility_m_per_n_s < 0.0 {
        return Err(ImplicitFrictionError::InvalidInput("discrete_mobility"));
    }
    let sign = if free_relative_velocity_m_s < 0.0 { -1.0 } else { 1.0 };
    let speed = free_relative_velocity_m_s.abs();
    let mut evaluate = |relative_speed: f64| {
        let force = traction(sign * relative_speed)
            .map_err(ImplicitFrictionError::Constitutive)?;
        if !force.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("traction"));
        }
        if force * sign < 0.0 || (relative_speed == 0.0 && force != 0.0) {
            return Err(ImplicitFrictionError::NonDissipative);
        }
        let velocity_change = mobility_m_per_n_s * force;
        if !velocity_change.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("contact_velocity_change"));
        }
        let endpoint = free_relative_velocity_m_s - velocity_change;
        let residual = sign * relative_speed - endpoint;
        if !endpoint.is_finite() || !residual.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("contact_residual"));
        }
        Ok(ImplicitFrictionStep {
            traction_n: force,
            relative_velocity_m_s: endpoint,
            residual_m_s: residual,
        })
    };
    // Admit the zero-slip endpoint too: this helper is not a static-reaction
    // solver, and cannot bracket a law that invents a force at rest.
    let at_zero = evaluate(0.0)?;
    if speed == 0.0 {
        return Ok(at_zero);
    }
    let at_free = evaluate(speed)?;
    if mobility_m_per_n_s == 0.0 || at_free.traction_n == 0.0 {
        return Ok(at_free);
    }
    let tolerance = (64.0 * f64::EPSILON * speed).max(f64::MIN_POSITIVE);
    let (mut lo, mut hi) = (0.0, speed);
    for _ in 0..96 {
        let mid = lo + 0.5 * (hi - lo);
        let candidate = evaluate(mid)?;
        if candidate.residual_m_s.abs() <= tolerance
            && candidate.relative_velocity_m_s * sign >= 0.0
        {
            return Ok(candidate);
        }
        if mid == lo || mid == hi {
            break;
        }
        if candidate.residual_m_s * sign > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Err(ImplicitFrictionError::NoConvergence)
}

/// Contact history that selects among the roots of one set-valued step.
///
/// A velocity-weakening law can meet the discrete load line at a stick reaction
/// and at two stable sliding roots. A contact stays on the branch it occupies
/// and leaves it only where that branch ends. There is no capture window: every
/// transition is fixed by the law, the normal load and the discrete mobility.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ContactBranch {
    /// No retained history: the first step after the load is applied.
    #[default]
    Fresh,
    /// Endpoint no-slip inside the static cone.
    Stuck,
    /// Sliding on the low-slip stable branch, below the residual's local maximum.
    SlidingLow {
        /// Slip direction; `true` when the driver is faster.
        positive: bool,
    },
    /// Sliding on the high-slip stable branch, above the residual's local minimum.
    SlidingHigh {
        /// Slip direction; `true` when the driver is faster.
        positive: bool,
    },
}

/// One committed set-valued contact step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SetValuedFrictionStep {
    /// Traction on the driven body [N].
    pub traction_n: f64,
    /// Driver-minus-body velocity after the held-force step [m/s].
    pub relative_velocity_m_s: f64,
    /// Branch the contact occupies after this step.
    pub branch: ContactBranch,
    /// Intersections of the discrete load line with the friction graph:
    /// the static cone (when it admits the held no-slip force) plus every
    /// sliding root. A static/kinetic junction-strength drop permits up to four.
    pub roots: u8,
    /// Zero-based index of the committed intersection, ordered by slip speed
    /// with the static cone first.
    pub selected: u8,
}

/// Discrete load line of the fs-tribo Stribeck law
/// `mu(v) = mu_k + (mu_s - mu_k) exp(-(v/c)^2) + sigma v` with its static cone
/// `|F| <= mu_s N`, at one normal load and one discrete mobility.
///
/// In slip speed `v` along the free direction, a held-force step solves
/// `g(v) = v - |v_free| + m N mu(v) = 0` on `(0, |v_free|]`. Since
/// `g'(v) = 1 + m N sigma - m N (mu_s - mu_k) (2v/c^2) exp(-(v/c)^2)` falls and
/// then rises, `g` has at most one local maximum `v_a` and one local minimum
/// `v_b`. There is at most one stable root below `v_a` (the low branch) and one
/// above `v_b` (the high branch). A root between them is unstable and is never
/// selected. Sticking is admissible exactly when the held no-slip force lies
/// inside the cone.
///
/// Root selection (McIntyre–Schumacher–Woodhouse hysteresis):
/// - From a fresh touch, from rest, or after the free slip reverses direction,
///   the contact sticks when admissible. Otherwise it takes the low root, or
///   the high root when the low root does not exist.
/// - A low-branch slide stays on its root while that root exists. At a finite
///   slip fold it jumps to the high root; unloading to zero slip permits arrest.
/// - A high-branch slide stays on its root while that root exists. When the root
///   folds away, the contact falls to the low root, or else sticks.
///
/// The retained memory is the occupied sliding branch and its direction.
///
/// The same line admits the fs-tribo hyperbolic law
/// `mu(v) = mu_k + (mu_s - mu_k)/(1 + v/v0)`. Its residual is convex, with
/// `g'(v) = 1 - m N (mu_s - mu_k) v0/(v + v0)^2`. When this derivative starts
/// negative, the critical slips are `0` and `sqrt(m N (mu_s - mu_k) v0) - v0`;
/// otherwise the residual is monotone and has only its low sliding branch.
///
/// [`Self::for_adhesive`] adds the actual junction forces `tau_static A` and
/// `tau_kinetic A`. The kinetic offset shifts `g` without moving its critical
/// points. When the two strengths differ, a stable low sliding root can coexist
/// with the static cone; its retained sliding branch continues until that root
/// ends. No signed net normal load or fictitious friction coefficient is used.
#[derive(Debug, Clone, PartialEq)]
pub struct StribeckLoadLine {
    law: fs_tribo::FrictionLaw,
    normal_n: f64,
    mobility_m_per_n_s: f64,
    static_capacity_n: f64,
    kinetic_junction_n: f64,
    zero_speed_capacity_n: f64,
    local_max_m_s: f64,
    local_min_m_s: f64,
}

type SetValuedError = ImplicitFrictionError<fs_tribo::TriboError>;

impl StribeckLoadLine {
    /// Admit a Stribeck or hyperbolic law at positive load and discrete mobility,
    /// and locate the residual's critical points once.
    ///
    /// # Errors
    /// An unsupported or invalid law, a non-positive or non-finite load or
    /// mobility, or an unrepresentable static capacity.
    pub fn new(
        law: &fs_tribo::FrictionLaw,
        normal_n: f64,
        mobility_m_per_n_s: f64,
    ) -> Result<Self, SetValuedError> {
        if !(normal_n.is_finite() && normal_n > 0.0) {
            return Err(ImplicitFrictionError::InvalidInput("normal_force"));
        }
        Self::with_junction(law, normal_n, mobility_m_per_n_s, 0.0, 0.0)
    }

    /// Continue a frozen adhesive junction's complete Stribeck/hyperbolic law.
    /// The real area and nonnegative repulsive load stay with the normal owner.
    /// Zero repulsive load is admitted: the independent junction still carries
    /// shear. This adds no normal evolution, shear peeling or elastic storage.
    ///
    /// # Errors
    /// An unsupported pressure law, non-positive or non-finite mobility, or
    /// an unrepresentable force capacity or load-line scale.
    pub fn for_adhesive(
        contact: &fs_tribo::adhesive::AdhesiveTangentialContact,
        mobility_m_per_n_s: f64,
    ) -> Result<Self, SetValuedError> {
        let law = contact.law();
        let area = contact.contact().area_m2();
        Self::with_junction(
            law.pressure_law(),
            contact.contact().repulsive_normal_force_n(),
            mobility_m_per_n_s,
            law.static_shear_strength_pa() * area,
            law.kinetic_shear_strength_pa() * area,
        )
    }

    fn with_junction(
        law: &fs_tribo::FrictionLaw,
        normal_n: f64,
        mobility_m_per_n_s: f64,
        static_junction_n: f64,
        kinetic_junction_n: f64,
    ) -> Result<Self, SetValuedError> {
        let (static_mu, kinetic_mu, characteristic_speed, viscous_per_speed, hyperbolic) = match *law {
            fs_tribo::FrictionLaw::Stribeck { static_mu, kinetic_mu, characteristic_speed, viscous_per_speed } =>
                (static_mu, kinetic_mu, characteristic_speed, viscous_per_speed, false),
            fs_tribo::FrictionLaw::Hyperbolic { static_mu, kinetic_mu, characteristic_speed } =>
                (static_mu, kinetic_mu, characteristic_speed, 0.0, true),
            _ => return Err(ImplicitFrictionError::InvalidInput("set-valued contact requires a Stribeck or hyperbolic law")),
        };
        // Validates the law through its owner before any root work.
        law.kinetic_coefficient(characteristic_speed)
            .map_err(ImplicitFrictionError::Constitutive)?;
        if !(normal_n.is_finite() && normal_n >= 0.0) {
            return Err(ImplicitFrictionError::InvalidInput("normal_force"));
        }
        if !(mobility_m_per_n_s.is_finite() && mobility_m_per_n_s > 0.0) {
            return Err(ImplicitFrictionError::InvalidInput("discrete_mobility"));
        }
        let static_capacity_n = static_junction_n + static_mu * normal_n;
        let zero_speed_capacity_n = kinetic_junction_n + static_mu * normal_n;
        if !static_capacity_n.is_finite() || !zero_speed_capacity_n.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("static_capacity"));
        }
        // With x = v/c, g'(v) = 0 exactly where x exp(-x^2) = b c / (2a); the
        // left side rises to its peak at x = 1/sqrt(2), then decays.
        let a = mobility_m_per_n_s * normal_n * (static_mu - kinetic_mu);
        let b = 1.0 + mobility_m_per_n_s * normal_n * viscous_per_speed;
        let shape = |x: f64| x * fs_math::det::exp(-(x * x));
        let peak = core::f64::consts::FRAC_1_SQRT_2;
        let (local_max_m_s, local_min_m_s) = if hyperbolic {
            // (v + v0)^2 = m N (mu_s - mu_k) v0 at the one minimum.
            let bottom = fs_math::det::sqrt(a * characteristic_speed) - characteristic_speed;
            if !bottom.is_finite() {
                return Err(ImplicitFrictionError::InvalidInput("load_line_scale"));
            }
            if bottom > 0.0 { (0.0, bottom) } else { (f64::INFINITY, f64::INFINITY) }
        } else if a > 0.0 && a.is_finite() && b.is_finite() {
            let target = b * characteristic_speed / (2.0 * a);
            if !(target > 0.0 && target.is_finite()) {
                return Err(ImplicitFrictionError::InvalidInput("load_line_scale"));
            }
            if target < shape(peak) {
                let rising = bisect_crossing(0.0, peak, |x| shape(x) < target);
                let mut far = 2.0 * peak;
                while shape(far) >= target {
                    far *= 2.0;
                }
                let falling = bisect_crossing(peak, far, |x| shape(x) > target);
                (characteristic_speed * rising, characteristic_speed * falling)
            } else {
                (f64::INFINITY, f64::INFINITY)
            }
        } else {
            (f64::INFINITY, f64::INFINITY)
        };
        Ok(Self {
            law: law.clone(),
            normal_n,
            mobility_m_per_n_s,
            static_capacity_n,
            kinetic_junction_n,
            zero_speed_capacity_n,
            local_max_m_s,
            local_min_m_s,
        })
    }

    /// Whether this line was built for exactly this load and mobility.
    #[must_use]
    pub fn matches(&self, normal_n: f64, mobility_m_per_n_s: f64) -> bool {
        self.normal_n.to_bits() == normal_n.to_bits()
            && self.mobility_m_per_n_s.to_bits() == mobility_m_per_n_s.to_bits()
    }

    /// Slip speeds of the residual's local maximum and minimum [m/s]. Both are
    /// infinite when the residual is monotone (no sliding hysteresis).
    #[must_use]
    pub const fn critical_slips_m_s(&self) -> (f64, f64) {
        (self.local_max_m_s, self.local_min_m_s)
    }

    /// Commit one held-force step from `previous`.
    ///
    /// # Errors
    /// A non-finite free slip, a constitutive refusal, or a root that does not
    /// resolve to floating-point tolerance.
    pub fn step(
        &self,
        free_relative_velocity_m_s: f64,
        previous: ContactBranch,
    ) -> Result<SetValuedFrictionStep, SetValuedError> {
        self.step_within(free_relative_velocity_m_s, previous, usize::MAX).map(|(step, _)| step)
    }

    /// The same step under a caller's bisection budget. Returns the step and
    /// the bisections its roots took; every root is resolved to adjacent
    /// floats, so a budget too small for that refuses.
    ///
    /// # Errors
    /// As [`Self::step`], and [`ImplicitFrictionError::NoConvergence`] when the
    /// roots need more than `max_bisections`.
    pub fn step_within(
        &self,
        free_relative_velocity_m_s: f64,
        previous: ContactBranch,
        max_bisections: usize,
    ) -> Result<(SetValuedFrictionStep, usize), SetValuedError> {
        let bisections = core::cell::Cell::new(0_usize);
        if !free_relative_velocity_m_s.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("free_relative_velocity"));
        }
        let mobility = self.mobility_m_per_n_s;
        let pin = free_relative_velocity_m_s / mobility;
        if !pin.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("sticking_force"));
        }
        let positive = free_relative_velocity_m_s > 0.0;
        let speed = free_relative_velocity_m_s.abs();
        let stick_admissible = pin.abs() <= self.static_capacity_n;
        // The static cone and the one-sided kinetic limit differ for a
        // junction-strength drop. A low sliding root can exist inside the cone.
        let starts_below = pin.abs() > self.zero_speed_capacity_n;
        let force = |v: f64| -> Result<f64, SetValuedError> {
            let mu = self.law.kinetic_coefficient(v).map_err(ImplicitFrictionError::Constitutive)?;
            Ok(self.kinetic_junction_n + mu * self.normal_n)
        };
        let residual = |v: f64| -> Result<f64, SetValuedError> {
            let r = v - speed + mobility * force(v)?;
            if r.is_finite() { Ok(r) } else { Err(ImplicitFrictionError::InvalidInput("contact_residual")) }
        };
        let root_on = |lo: f64, hi: f64| -> Result<f64, SetValuedError> {
            let mut failure = None;
            let root = bisect_crossing(lo, hi, |v| {
                bisections.set(bisections.get() + 1);
                match residual(v) {
                    Ok(r) => r < 0.0,
                    Err(e) => {
                        failure.get_or_insert(e);
                        false
                    }
                }
            });
            failure.map_or(Ok(root), Err)
        };
        let low = || -> Result<Option<f64>, SetValuedError> {
            let top = self.local_max_m_s.min(speed);
            if !starts_below || top <= 0.0 || residual(top)? < 0.0 {
                return Ok(None);
            }
            root_on(0.0, top).map(Some)
        };
        let high = || -> Result<Option<f64>, SetValuedError> {
            let bottom = self.local_min_m_s;
            if bottom >= speed || residual(bottom)? > 0.0 {
                return Ok(None);
            }
            root_on(bottom, speed).map(Some)
        };
        let reversed = matches!(previous,
            ContactBranch::SlidingLow { positive: p } | ContactBranch::SlidingHigh { positive: p }
                if p != positive);
        let stays_high = matches!(previous, ContactBranch::SlidingHigh { .. }) && !reversed;
        let stays_low = matches!(previous, ContactBranch::SlidingLow { .. }) && !reversed;
        let choice = if stays_high {
            match high()? {
                Some(v) => Some((v, true)),
                None => low()?.map(|v| (v, false)),
            }
        } else if stays_low {
            match low()? {
                Some(v) => Some((v, false)),
                // Loss at a finite-slip fold goes to the high branch, even
                // inside the larger static cone. Only reaching v=0 on
                // unloading permits arrest without reversing the drive.
                None if stick_admissible && !starts_below => None,
                None => high()?.map(|v| (v, true)),
            }
        } else if stick_admissible {
            None
        } else {
            match low()? {
                Some(v) => Some((v, false)),
                None => high()?.map(|v| (v, true)),
            }
        };
        // Outside the cone the residual is negative at rest and positive at the
        // free slip, so a root exists. A branch test lost to rounding at a fold
        // falls back to the first crossing instead of refusing the step.
        let choice = match choice {
            None if !stick_admissible => Some((root_on(0.0, speed)?, false)),
            other => other,
        };
        // Disclosure. The residual rises on (0, a], falls on (a, b] and rises
        // on (b, speed], with a and b its critical slips clipped to the free
        // slip. The low root exists when g(0+) < 0 and g(a) > 0, the
        // middle one when g(a) > 0 > g(b), and the high one when g(b) <= 0
        // (the same existence tests the branch selection uses).
        let a = self.local_max_m_s.min(speed);
        let b = self.local_min_m_s.min(speed);
        // An empty rising piece uses the KINETIC one-sided capacity; the
        // static junction strength is a separate set-valued reaction bound.
        let g_a = if a > 0.0 { residual(a)? } else { mobility * self.zero_speed_capacity_n - speed };
        let g_b = if b < speed { residual(b)? } else { 1.0 };
        let has_low = starts_below && g_a >= 0.0;
        let has_middle = b < speed && g_a > 0.0 && g_b < 0.0;
        let has_high = b < speed && g_b <= 0.0;
        let counted = u8::from(stick_admissible)
            + u8::from(has_low)
            + u8::from(has_middle)
            + u8::from(has_high);
        let used = bisections.get();
        if used > max_bisections {
            return Err(ImplicitFrictionError::NoConvergence);
        }
        let Some((slip, on_high)) = choice else {
            return Ok((SetValuedFrictionStep {
                traction_n: pin,
                relative_velocity_m_s: 0.0,
                branch: ContactBranch::Stuck,
                roots: counted,
                selected: 0,
            }, used));
        };
        // The high root is the last intersection. Otherwise this is the first
        // sliding root: the low one, or the fold fallback above.
        let selected = if on_high {
            counted.saturating_sub(1)
        } else {
            u8::from(stick_admissible)
        };
        let roots = counted.max(selected + 1);
        let sign = if positive { 1.0 } else { -1.0 };
        let traction = sign * force(slip)?;
        let endpoint = free_relative_velocity_m_s - mobility * traction;
        let tolerance = 64.0 * f64::EPSILON * speed.max(mobility * traction.abs());
        if !endpoint.is_finite() || endpoint * sign < 0.0 || (endpoint - sign * slip).abs() > tolerance {
            return Err(ImplicitFrictionError::NoConvergence);
        }
        Ok((SetValuedFrictionStep {
            traction_n: traction,
            relative_velocity_m_s: endpoint,
            branch: if on_high {
                ContactBranch::SlidingHigh { positive }
            } else {
                ContactBranch::SlidingLow { positive }
            },
            roots,
            selected,
        }, used))
    }
}

/// Bisect `[lo, hi]` to adjacent floating-point values, keeping `below(lo)`
/// true and `below(hi)` false; returns `hi`. Deterministic and allocation-free.
/// Terminates because every pass halves the interval until no midpoint exists.
fn bisect_crossing(mut lo: f64, mut hi: f64, mut below: impl FnMut(f64) -> bool) -> f64 {
    loop {
        let mid = lo + 0.5 * (hi - lo);
        if mid <= lo || mid >= hi {
            return hi;
        }
        if below(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
}

#[cfg(test)]
mod set_valued_contact_tests {
    use super::*;

    fn law(static_mu: f64, kinetic_mu: f64, speed: f64, viscous: f64) -> fs_tribo::FrictionLaw {
        fs_tribo::FrictionLaw::Stribeck {
            static_mu,
            kinetic_mu,
            characteristic_speed: speed,
            viscous_per_speed: viscous,
        }
    }

    fn adhesive_line(normal: f64, mobility: f64, static_force: f64, kinetic_force: f64)
        -> StribeckLoadLine
    {
        use fs_tribo::adhesive::{AdhesiveContactState, AdhesiveFrictionLaw, AdhesiveTangentialContact};
        let area = 1e-6;
        let junction = AdhesiveTangentialContact::new(
            AdhesiveFrictionLaw::new(law(1.2, 0.35, 0.05, 0.0),
                static_force / area, kinetic_force / area).unwrap(),
            AdhesiveContactState::new(area, normal).unwrap(),
        ).unwrap();
        StribeckLoadLine::for_adhesive(&junction, mobility).unwrap()
    }

    /// G1: a real junction retains kinetic sliding inside its larger static
    /// cone. An algebraic no-slip test on every step would falsely arrest it.
    #[test]
    fn adhesive_strength_drop_preserves_sliding_at_zero_repulsive_load() {
        let line = adhesive_line(0.0, 2.0, 0.8, 0.4);
        let stuck = line.step(1.2, ContactBranch::Fresh).unwrap();
        assert_eq!(stuck.branch, ContactBranch::Stuck);
        assert_eq!(stuck.roots, 2);
        assert_eq!(stuck.traction_n, 0.6);
        let sliding = line.step(1.2, ContactBranch::SlidingLow { positive: true }).unwrap();
        assert_eq!(sliding.branch, ContactBranch::SlidingLow { positive: true });
        assert_eq!((sliding.roots, sliding.selected), (2, 1));
        assert!((sliding.traction_n - 0.4).abs() < 1e-15);
        assert!((sliding.relative_velocity_m_s - 0.4).abs() < 1e-15);
        assert_eq!(line.step(0.7, sliding.branch).unwrap().branch, ContactBranch::Stuck);
        let reversed = line.step(-1.2, sliding.branch).unwrap();
        assert_eq!(reversed.branch, ContactBranch::Stuck);
        assert_eq!(reversed.traction_n, -0.6);

        let frictionless_slide = adhesive_line(0.0, 2.0, 0.8, 0.0)
            .step(1.2, sliding.branch).unwrap();
        assert_eq!(frictionless_slide.traction_n, 0.0);
        assert_eq!(frictionless_slide.relative_velocity_m_s, 1.2);
        assert_eq!((frictionless_slide.roots, frictionless_slide.selected), (2, 1));
    }

    /// G1: independently scan the physical equation including the junction
    /// offset. Here the cone coexists with THREE sliding roots; the middle
    /// root is unstable. Source forces are never divided by the normal load.
    #[test]
    fn adhesive_four_intersections_keep_both_stable_sliding_branches() {
        let (normal, mobility, speed) = (1.0, 0.1, 0.128);
        let line = adhesive_line(normal, mobility, 0.1, 0.05);
        let pressure = law(1.2, 0.35, 0.05, 0.0);
        let g = |v: f64| {
            let mu = if v == 0.0 { 1.2 } else { pressure.kinetic_coefficient(v).unwrap() };
            v - speed + mobility * (0.05 + normal * mu)
        };
        let mut roots = Vec::new();
        let mut left = 0.0;
        let mut previous = g(left);
        for i in 1..=20_000_u32 {
            let right = speed * f64::from(i) / 20_000.0;
            let current = g(right);
            if previous.is_sign_negative() != current.is_sign_negative() {
                let (mut lo, mut hi) = (left, right);
                for _ in 0..60 {
                    let mid = f64::midpoint(lo, hi);
                    if g(mid).is_sign_negative() == previous.is_sign_negative() { lo = mid; } else { hi = mid; }
                }
                roots.push(f64::midpoint(lo, hi));
            }
            left = right;
            previous = current;
        }
        assert_eq!(roots.len(), 3);
        for (branch, selected, slip) in [
            (ContactBranch::Fresh, 0, 0.0),
            (ContactBranch::SlidingLow { positive: true }, 1, roots[0]),
            (ContactBranch::SlidingHigh { positive: true }, 3, roots[2]),
        ] {
            let step = line.step(speed, branch).unwrap();
            assert_eq!((step.roots, step.selected), (4, selected));
            assert!((step.relative_velocity_m_s - slip).abs() < 1e-14);
            assert!((speed - mobility * step.traction_n - slip).abs() < 1e-14);
            assert!(step.traction_n * step.relative_velocity_m_s >= 0.0);
        }
    }

    /// G1: increasing drive past the low sliding fold cannot arrest a moving
    /// junction merely because its larger static cone still contains the force.
    #[test]
    fn adhesive_low_branch_jumps_at_fold_before_static_threshold() {
        let line = adhesive_line(1.0, 0.1, 0.2, 0.05);
        let low = line.step(0.128, ContactBranch::SlidingLow { positive: true }).unwrap();
        assert_eq!(low.branch, ContactBranch::SlidingLow { positive: true });
        let high = line.step(0.133, low.branch).unwrap();
        assert_eq!(high.branch, ContactBranch::SlidingHigh { positive: true });
        assert_eq!(high.selected, high.roots - 1);
        assert!(high.relative_velocity_m_s > low.relative_velocity_m_s);
        let remains_high = line.step(0.128, high.branch).unwrap();
        assert_eq!(remains_high.branch, high.branch);
        assert_eq!(line.step(0.119, remains_high.branch).unwrap().branch, ContactBranch::Stuck);
        assert_eq!(line.step(0.133, ContactBranch::Stuck).unwrap().branch, ContactBranch::Stuck);
    }

    /// G3: zero junction strength reproduces the original pressure graph and
    /// its branch selection exactly, including reversal and fold crossings.
    #[test]
    fn zero_junction_strength_recovers_dry_continuation() {
        let pressure = law(1.2, 0.35, 0.05, 0.0);
        let dry = StribeckLoadLine::new(&pressure, 1.0, 0.1).unwrap();
        let adhesive = adhesive_line(1.0, 0.1, 0.0, 0.0);
        for previous in [ContactBranch::Fresh, ContactBranch::Stuck,
            ContactBranch::SlidingLow { positive: true }, ContactBranch::SlidingHigh { positive: true }] {
            for free in [-0.3, -0.128, 0.0, 0.1, 0.12, 0.128, 0.3] {
                assert_eq!(dry.step(free, previous).unwrap(), adhesive.step(free, previous).unwrap());
            }
        }
    }

    // Independent oracle: every sign change of g on a uniform grid of (0, S].
    fn upward_crossings(law: &fs_tribo::FrictionLaw, n: f64, m: f64, s: f64) -> Vec<f64> {
        let cells = 200_000;
        let g = |v: f64| {
            let mu = if v == 0.0 {
                match *law {
                    fs_tribo::FrictionLaw::Stribeck { static_mu, .. } => static_mu,
                    _ => unreachable!(),
                }
            } else {
                law.kinetic_coefficient(v).unwrap()
            };
            v - s + m * mu * n
        };
        let mut out = Vec::new();
        let mut prev = g(0.0);
        for i in 1..=cells {
            let v = s * f64::from(i) / f64::from(cells);
            let cur = g(v);
            if prev < 0.0 && cur >= 0.0 {
                out.push(v);
            }
            prev = cur;
        }
        out
    }

    fn lcg(state: &mut u64) -> f64 {
        *state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        f64::from(u32::try_from(*state >> 40).unwrap()) / f64::from(1_u32 << 24)
    }

    /// G1: from rest the step takes the first upward crossing; a high-branch
    /// slide takes the last one; sticking happens exactly inside the cone.
    #[test]
    fn set_valued_roots_match_a_dense_scan_oracle() {
        let mut seed = 0x5eed_u64;
        let mut hysteretic = 0;
        for _ in 0..400 {
            let kinetic = 0.1 + 0.5 * lcg(&mut seed);
            let static_mu = kinetic + 0.6 * lcg(&mut seed);
            let c = 0.005 + 0.1 * lcg(&mut seed);
            let viscous = if lcg(&mut seed) < 0.3 { 2.0 * lcg(&mut seed) } else { 0.0 };
            let f = law(static_mu, kinetic, c, viscous);
            let n = 0.05 + 5.0 * lcg(&mut seed);
            let m = 0.01 + 0.5 * lcg(&mut seed);
            let line = StribeckLoadLine::new(&f, n, m).unwrap();
            let s = 1.5 * lcg(&mut seed) * (m * static_mu * n + 0.2);
            let roots = upward_crossings(&f, n, m, s);
            let cell = s / 200_000.0;
            let fresh = line.step(s, ContactBranch::Fresh).unwrap();
            let high = line.step(s, ContactBranch::SlidingHigh { positive: true }).unwrap();
            if s / m <= static_mu * n {
                assert_eq!(fresh.branch, ContactBranch::Stuck);
                assert_eq!(fresh.relative_velocity_m_s, 0.0);
                assert_eq!(fresh.traction_n, s / m);
                if !roots.is_empty() {
                    hysteretic += 1;
                }
            } else {
                assert!((fresh.relative_velocity_m_s - roots[0]).abs() <= cell, "{fresh:?} vs {roots:?}");
            }
            match roots.last() {
                Some(&last) => assert!((high.relative_velocity_m_s - last).abs() <= cell, "{high:?} vs {roots:?}"),
                None => assert_eq!(high.branch, ContactBranch::Stuck),
            }
            // Mirror symmetry is exact.
            let mirror = line.step(-s, ContactBranch::Fresh).unwrap();
            assert_eq!(mirror.traction_n.to_bits(), (-fresh.traction_n).to_bits());
        }
        assert!(hysteretic >= 10, "the sample must exercise the multi-root region: {hysteretic}");
    }

    /// G1: sweeping the free slip up then down releases at the cone and
    /// recaptures where the high root folds away, below the cone for a weakening law.
    #[test]
    fn hysteresis_loop_releases_at_the_cone_and_recaptures_at_the_fold() {
        let f = law(0.8, 0.4, 0.04, 0.0);
        let (n, m) = (3.9, 0.2);
        let line = StribeckLoadLine::new(&f, n, m).unwrap();
        let cone = m * 0.8 * n;
        let mut branch = ContactBranch::Fresh;
        let ladder: Vec<f64> = (0..=4000).map(|i| f64::from(i) * 1.0e-3).collect();
        let mut release = None;
        for &s in &ladder {
            branch = line.step(s, branch).unwrap().branch;
            if release.is_none() && branch != ContactBranch::Stuck {
                release = Some(s);
            }
        }
        let mut recapture = None;
        for &s in ladder.iter().rev() {
            branch = line.step(s, branch).unwrap().branch;
            if recapture.is_none() && branch == ContactBranch::Stuck {
                recapture = Some(s);
            }
        }
        let (release, recapture) = (release.unwrap(), recapture.unwrap());
        assert!(release >= cone - 1.0e-12 && release - cone <= 1.0e-3 + 1.0e-12, "{release} vs {cone}");
        // Fold: the high root exists while S >= v + m N mu(v) for some slip v.
        // Its smallest value comes from an independent dense scan over (0, 1] m/s.
        let fold = (1..=1_000_000)
            .map(|i| f64::from(i) * 1.0e-6)
            .map(|v| v + m * n * f.kinetic_coefficient(v).unwrap())
            .fold(f64::INFINITY, f64::min);
        assert!(recapture < fold && fold - recapture <= 1.0e-3 + 1.0e-9, "{recapture} vs {fold}");
        assert!(recapture < 0.8 * release, "weakening law must show hysteresis");
    }

    /// G1: equal coefficients (Coulomb) have no hysteresis and a closed form.
    #[test]
    fn coulomb_limit_is_unique_and_closed_form() {
        let f = law(0.5, 0.5, 0.04, 0.0);
        let (n, m) = (2.0, 0.3);
        let line = StribeckLoadLine::new(&f, n, m).unwrap();
        assert_eq!(line.critical_slips_m_s(), (f64::INFINITY, f64::INFINITY));
        for s in [0.0, 0.1, 0.29, 0.3, 0.31, 1.0, 7.0] {
            for previous in [ContactBranch::Fresh, ContactBranch::SlidingHigh { positive: true }] {
                let step = line.step(s, previous).unwrap();
                if s <= m * 0.5 * n {
                    assert_eq!(step.branch, ContactBranch::Stuck);
                } else {
                    assert_eq!(step.traction_n, 0.5 * n);
                    assert!((step.relative_velocity_m_s - (s - m * 0.5 * n)).abs() < 1.0e-15);
                }
            }
        }
    }

    /// G0: refusals are typed; nothing is defaulted.
    #[test]
    fn set_valued_line_refuses_bad_inputs() {
        let f = law(0.8, 0.4, 0.04, 0.0);
        let coulomb = fs_tribo::FrictionLaw::Coulomb { static_mu: 0.5, kinetic_mu: 0.4 };
        assert!(matches!(StribeckLoadLine::new(&coulomb, 1.0, 0.1), Err(ImplicitFrictionError::InvalidInput(_))));
        for (n, m) in [(0.0, 0.1), (-1.0, 0.1), (f64::NAN, 0.1), (1.0, 0.0), (1.0, -0.1), (1.0, f64::INFINITY)] {
            assert!(matches!(StribeckLoadLine::new(&f, n, m), Err(ImplicitFrictionError::InvalidInput(_))));
        }
        assert!(matches!(StribeckLoadLine::new(&law(0.4, 0.8, 0.04, 0.0), 1.0, 0.1),
            Err(ImplicitFrictionError::Constitutive(_))));
        let line = StribeckLoadLine::new(&f, 1.0, 0.1).unwrap();
        assert!(matches!(line.step(f64::NAN, ContactBranch::Fresh), Err(ImplicitFrictionError::InvalidInput(_))));
        let next = f64::from_bits(0.1_f64.to_bits() + 1);
        assert!(line.matches(1.0, 0.1) && !line.matches(1.0, next) && !line.matches(1.5, 0.1));
    }
}

#[cfg(test)]
mod implicit_contact_tests {
    use super::*;

    #[test]
    fn implicit_linear_contact_matches_closed_form_across_mechanical_scales() {
        for free in [-3.0, -0.01, 0.0, 0.01, 3.0] {
            for mobility in [0.0, 1.0e-6, 0.01, 1.0, 100.0] {
                let damping = 4.0;
                let step = solve_implicit_friction(free, mobility,
                    |v| Ok::<_, ()>(damping * v)).unwrap();
                let exact = free / (1.0 + mobility * damping);
                let tol = 2.0e-12 * free.abs().max(1.0e-20);
                assert!((step.relative_velocity_m_s - exact).abs() < tol);
                assert!((step.traction_n - damping * exact).abs() < 4.0 * tol);
                assert!(step.traction_n * step.relative_velocity_m_s >= 0.0);
            }
        }
    }

    #[test]
    fn implicit_material_law_retains_constitutive_ownership_and_odd_symmetry() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.05).unwrap();
        for free in [0.01, 0.2, 2.0] {
            for mobility in [1.0e-5, 0.01, 0.2] {
                let solve = |sign| solve_implicit_friction(sign * free, mobility,
                    |v| law.traction(v, 0.7)).unwrap();
                let a = solve(1.0);
                let b = solve(-1.0);
                assert_eq!(a.traction_n.to_bits(), (-b.traction_n).to_bits());
                let owned = law.traction(a.relative_velocity_m_s, 0.7).unwrap();
                assert!((owned - a.traction_n).abs() < 1.0e-11);
                assert!((free - mobility * a.traction_n
                    - a.relative_velocity_m_s).abs() < 1.0e-14);
            }
        }
    }

    #[test]
    fn implicit_contact_prevents_explicit_velocity_reversal() {
        // An explicit viscous force gives v_next = 1 - 10 * 2 = -19 m/s.
        let step = solve_implicit_friction(1.0, 10.0,
            |v| Ok::<_, ()>(2.0 * v)).unwrap();
        assert!((step.relative_velocity_m_s - 1.0 / 21.0).abs() < 1.0e-12);
        assert!(step.traction_n * step.relative_velocity_m_s > 0.0);
    }

    #[test]
    fn implicit_contact_refuses_active_laws_bad_response_and_constitutive_jumps() {
        assert!(matches!(solve_implicit_friction(1.0, -1.0,
            |v| Ok::<_, ()>(v)), Err(ImplicitFrictionError::InvalidInput(_))));
        assert!(matches!(solve_implicit_friction(f64::NAN, 1.0,
            |v| Ok::<_, ()>(v)), Err(ImplicitFrictionError::InvalidInput(_))));
        assert!(matches!(solve_implicit_friction(1.0, 1.0,
            |v| Ok::<_, ()>(-v)), Err(ImplicitFrictionError::NonDissipative)));
        assert_eq!(solve_implicit_friction(1.0, 1.0,
            |_| Err::<f64, _>("outside material domain")),
            Err(ImplicitFrictionError::Constitutive("outside material domain")));
        // There is no root: left limit is -0.5 and right limit is +0.5.
        assert_eq!(solve_implicit_friction(1.0, 1.0,
            |v| Ok::<_, ()>(if v < 0.5 { 0.0 } else { 1.0 })),
            Err(ImplicitFrictionError::NoConvergence));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_in_velocity_and_falls_to_kinetic() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.05).unwrap();
        let f_p = law.traction(0.01, 2.0).unwrap();
        let f_m = law.traction(-0.01, 2.0).unwrap();
        assert!((f_p + f_m).abs() < 1.0e-15);
        let slope = law.traction(0.01, 1.0).unwrap() / 0.01;
        let endpoint_mu = 0.3 + 0.5 * fs_math::det::exp(-1.0);
        assert!((slope - endpoint_mu / 0.05).abs() < 1.0e-12);
        assert!((law.traction(10.0, 1.0).unwrap() - 0.3).abs() < 0.02);
    }

    #[test]
    fn decay_speed_and_ode_ramp_are_independent() {
        let narrow = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.001).unwrap();
        let wide = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.01).unwrap();
        // Both are beyond their ramp: the same physical decay law applies.
        assert_eq!(
            narrow.traction(0.1, 2.0).unwrap().to_bits(),
            wide.traction(0.1, 2.0).unwrap().to_bits()
        );
        assert_ne!(
            narrow.traction(0.0005, 2.0).unwrap().to_bits(),
            wide.traction(0.0005, 2.0).unwrap().to_bits()
        );
        for law in [narrow, wide] {
            let scale = law.regularization_m_s;
            let left = law.traction(scale * (1.0 - 1e-12), 2.0).unwrap();
            let right = law.traction(scale * (1.0 + 1e-12), 2.0).unwrap();
            assert!((left - right).abs() <= 1e-12 * 0.8 * 2.0);
        }
    }

    /// Regression (bead 9svup): non-finite inputs REFUSE instead of
    /// publishing zero traction; misconfigured coefficients refuse at
    /// construction.
    #[test]
    fn nonfinite_inputs_refuse_as_typed_errors() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.05).unwrap();
        assert!(law.traction(f64::NAN, 1.0).is_err());
        assert!(law.traction(f64::INFINITY, 1.0).is_err());
        assert!(law.traction(0.01, f64::NAN).is_err());
        assert_eq!(law.traction(0.01, 0.0).unwrap().to_bits(), 0.0f64.to_bits());
        assert!(StribeckFriction::try_new(0.3, 0.8, 0.05, 0.05).is_err());
        assert!(StribeckFriction::try_new(f64::NAN, 0.3, 0.05, 0.05).is_err());
        assert!(StribeckFriction::try_new(0.8, 0.3, 0.0, 0.05).is_err());
        assert!(StribeckFriction::try_new(0.8, 0.3, 0.05, 0.0).is_err());
    }

    /// G1/G3: equal-and-opposite contact forces remove relative kinetic
    /// energy; the independent formula also checks the port's sign and SI load.
    #[test]
    fn shared_friction_matches_reference_and_dissipates_relative_work() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.05).unwrap();
        let owner = fs_tribo::FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.0,
        };
        for velocity in [-10.0_f64, -0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 10.0] {
            for load in [0.0, 0.1, 2.0, 100.0] {
                let force = law.traction(velocity, load).unwrap();
                assert_eq!(
                    force.to_bits(),
                    owner
                        .regularized_traction_1d(-velocity, load, 0.05)
                        .unwrap()
                        .to_bits()
                );
                let ratio = velocity.abs() / 0.05;
                let mu = if ratio < 1.0 {
                    // DEF-06: continuous kinetic endpoint rather than μ_static.
                    (0.3 + 0.5 * (-1.0_f64).exp()) * ratio
                } else {
                    0.3 + 0.5 * (-ratio * ratio).exp()
                };
                let reference = mu * load * velocity.signum();
                assert!((force - reference).abs() <= 1.0e-12 * load.max(1.0));
                let body_velocity = 0.7;
                let driver_velocity = body_velocity + velocity;
                let pair_power = force * body_velocity - force * driver_velocity;
                assert!(pair_power <= 1.0e-12);
                assert!((pair_power + force * velocity).abs() <= 1.0e-12);
            }
        }
    }

    #[test]
    fn shared_friction_refuses_forged_coefficients_negative_load_and_overflow() {
        let valid = StribeckFriction::try_new(0.8, 0.3, 0.05, 0.05).unwrap();
        for invalid in [
            StribeckFriction {
                mu_static: -0.8,
                ..valid
            },
            StribeckFriction {
                mu_dynamic: 0.9,
                ..valid
            },
            StribeckFriction {
                mu_dynamic: f64::NAN,
                ..valid
            },
            StribeckFriction {
                stiction_m_s: 0.0,
                ..valid
            },
            StribeckFriction {
                regularization_m_s: 0.0,
                ..valid
            },
        ] {
            assert!(invalid.traction(0.01, 2.0).is_err());
        }
        assert_eq!(valid.traction(0.01, -1.0), Err("normal_force"));
        let large = StribeckFriction::try_new(4.0, 4.0, 0.05, 0.05).unwrap();
        assert_eq!(large.traction(0.1, f64::MAX), Err("regularized_traction"));
    }
}
