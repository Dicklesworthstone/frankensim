//! Forced-response sensitivities of the original mechanical coefficients.
//!
//! At fixed physical frequency and fixed load, A q = F gives
//! A q_theta = F_theta - K_theta q + w^2 M_theta q - C_theta v.
//! The existing factor solves this equation; no eigenmode differentiation,
//! finite perturbation, second factorization or derivative of numerical scaling
//! is involved. The full signed velocity direction includes gyroscopic terms.

use super::{
    HarmonicResponse, HarmonicResponseError, PreparedHarmonicResponse, finite, invalid, poll,
    scale, scale_complex, symmetric_entry,
};
use crate::quadratic_modes::QuadraticParameterDirection;
use crate::{quadratic_fraction, quadratic_power_scale};
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};

mod load;
pub use load::HarmonicPortParameterDirection;

/// Additional resources for one physical parameter direction and its response.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicParameterControls {
    /// Conservative added numerical payload: 512*n + 384*p bytes, excluding
    /// the prepared factor, borrowed direction/source arrays and allocator metadata.
    /// A changing port map/impedance or loaded frequency adds 192*p bytes
    /// for its factored action.
    pub maximum_additional_payload_bytes: usize,
    /// Work score: 256*(n*n+n*p+p*p)+256*(n+p+k), plus 64*(n*n*p+k*k)
    /// for a compact loaded factor of dimension k. Includes both solves,
    /// audits and binary-normalized physical power contractions.
    /// A changing port map/impedance or loaded frequency adds
    /// 256*(n*p+p*p)+256*p; changing frequency also adds 192*n*n+256*n.
    pub maximum_work: u64,
}

/// Derivatives of the original physical load variables per unit parameter.
#[derive(Clone, Debug)]
pub struct HarmonicPortParameterResponse {
    /// B^T v_theta + B_theta^T v; the second term is zero for a fixed attachment.
    pub flow_derivative: Vec<C64>,
    /// Z flow_theta + Z_theta flow; the second term is zero for a fixed impedance.
    pub effort_derivative: Vec<C64>,
    /// d Re(flow^H effort)/2, signed boundary work per unit parameter.
    pub average_absorbed_power_derivative: f64,
}

/// Forced motion and physical work sensitivities in the original coordinates.
#[derive(Clone, Debug)]
pub struct HarmonicParameterResponse {
    /// Original response, with the same physical meaning as `solve`.
    pub response: HarmonicResponse,
    /// dq/dtheta, without changing or normalizing the source coordinates.
    pub displacement_derivative: Vec<C64>,
    /// -i(omega dq/dtheta + omega_theta q), retaining the DC kinematic term.
    pub velocity_derivative: Vec<C64>,
    /// Load flow/effort derivatives when a physical load is attached.
    pub port_load: Option<HarmonicPortParameterResponse>,
    /// d Re(F^H v)/2, including the authored force derivative.
    pub average_input_power_derivative: f64,
    /// d Re(v^H D v)/2, including D_theta=(C_theta+C_theta^T)/2.
    pub average_viscous_power_derivative: f64,
    /// Input minus independent viscous and boundary-work derivatives.
    pub power_balance_residual_derivative: f64,
    /// Maximum tangent-solve and published-coordinate equation backward error.
    /// This is not a forward-error bound near an ill-conditioned resonance.
    pub relative_derivative_residual: f64,
}

// Keep products in binary-normalized form until their COMMON RHS amplitude is
// known. Physical unit changes must not overflow w*w or dM*q before cancellation.
// DD accumulation retains weak coefficient directions beside a large preload;
// it is not a relative-accuracy certificate at an arbitrary exact cancellation.
#[derive(Clone, Copy)]
struct Term {
    fraction: Dd,
    exponent: i32,
}

impl Term {
    fn product(factors: [f64; 4], shift: i32) -> Self {
        if factors.contains(&0.0) {
            return Self {
                fraction: Dd::ZERO,
                exponent: 0,
            };
        }
        let mut fraction = Dd::ONE;
        let mut exponent = shift;
        for factor in factors {
            let (f, e) = quadratic_fraction(factor);
            fraction = fraction * Dd::from_f64(f.copysign(factor));
            exponent += e;
        }
        Self { fraction, exponent }
    }

    fn scaled(self, exponent: i32) -> Result<Dd, HarmonicResponseError> {
        let shift = self.exponent - exponent;
        Ok(Dd::from_pair(
            scale(self.fraction.hi, shift)?,
            quadratic_power_scale(self.fraction.lo, shift),
        ))
    }
}

// Physical power terms may individually exceed f64 even when their signed
// sum is zero or finite (source-phase rotations and damping nullspaces).
// Preserve the product exponent until AFTER cancellation, without another
// traversal or temporary term array. The bounded number of normalized terms
// fits in DD's finite fraction range; this is not an exact summation theorem.
#[derive(Default)]
struct PowerSum {
    fraction: Dd,
    exponent: i32,
}

impl PowerSum {
    fn add(&mut self, term: Term) -> Result<(), HarmonicResponseError> {
        if term.fraction == Dd::ZERO {
            return Ok(());
        }
        if self.fraction == Dd::ZERO {
            self.fraction = term.fraction;
            self.exponent = term.exponent;
            return Ok(());
        }
        if term.exponent > self.exponent {
            self.fraction = Term {
                fraction: self.fraction,
                exponent: self.exponent,
            }
            .scaled(term.exponent)?;
            self.exponent = term.exponent;
        }
        self.fraction = self.fraction + term.scaled(self.exponent)?;
        Ok(())
    }

    // c * Re(a^H b)/2. Even c/2 stays inside the factored product so a
    // subnormal physical damping coefficient is not prematurely rounded away.
    fn half_inner(&mut self, a: C64, b: C64, c: f64) -> Result<(), HarmonicResponseError> {
        self.add(Term::product([0.5, c, a.re, b.re], 0))?;
        self.add(Term::product([0.5, c, a.im, b.im], 0))
    }

    fn value(self) -> Result<f64, HarmonicResponseError> {
        scale(self.fraction.to_f64(), self.exponent)
    }
}

/// Differentiate signed cycle-mean work of original power-conjugate phasors.
///
/// Returns `Re(flow_derivative^H*effort + flow^H*effort_derivative)/2`, summed
/// over every original port. The derivative can be with respect to frequency
/// or any authored physical parameter; all four arrays use the same ordering.
/// Work exchanged through a port is not necessarily local dissipative heat.
///
/// Reuses the structural harmonic owner's binary-normalized double-double
/// products. Individual products may exceed f64 when their signed sum is
/// finite, including purely reactive changes and cancelling mutual exchanges.
/// Rounding occurs after the complete contraction. This is not an exact-sum
/// or forward-error certificate at an arbitrary cancellation.
///
/// Equal finite slices are required; an empty port set returns zero. Work is
/// linear in the borrowed port count with constant scratch and no allocation.
/// Cancellation is checked before, during and after contraction. A nonzero
/// final value outside f64's representable range refuses, without a partial sum.
pub fn harmonic_mean_power_derivative_with_cancel(
    flow: &[C64],
    effort: &[C64],
    flow_derivative: &[C64],
    effort_derivative: &[C64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<f64, HarmonicResponseError> {
    if cancelled() {
        return Err(HarmonicResponseError::Cancelled);
    }
    let n = flow.len();
    if effort.len() != n || flow_derivative.len() != n || effort_derivative.len() != n {
        return Err(invalid("equal original port and derivative lengths required"));
    }
    let mut sum = PowerSum::default();
    for i in 0..n {
        if cancelled() {
            return Err(HarmonicResponseError::Cancelled);
        }
        if [flow[i], effort[i], flow_derivative[i], effort_derivative[i]]
            .iter()
            .any(|&value| !finite(value))
        {
            return Err(invalid("finite original port phasors and derivatives required"));
        }
        sum.half_inner(flow_derivative[i], effort[i], 1.0)?;
        sum.half_inner(flow[i], effort_derivative[i], 1.0)?;
    }
    let value = sum.value()?;
    if cancelled() {
        return Err(HarmonicResponseError::Cancelled);
    }
    Ok(value)
}

impl PreparedHarmonicResponse {
    /// Differentiate the forced response along original M/C/K and F directions.
    ///
    /// Direction matrices are complete row-major n*n arrays per caller-declared
    /// unit of theta, using the same coordinates as the prepared source. M_theta
    /// and K_theta must be exactly symmetric, but need not be positive: a
    /// derivative is signed. C_theta may contain nonproportional viscous and
    /// skew gyroscopic changes. None for F_theta explicitly holds force fixed.
    ///
    /// Frequency, physical coordinate basis, B and Z are held fixed. This is
    /// therefore a mechanical-pencil derivative even for a loaded response,
    /// not a total derivative of changing attachment geometry or material Z.
    /// The caller supplies constitutive/assembly derivatives; they are never
    /// inferred from material names, eigenvalue slopes or neighboring samples.
    ///
    /// Reuses this owner's LU/scalar/compact factor. A direction may introduce
    /// coupling into an exactly diagonal base. No second dense operator is
    /// allocated. Independent RHS normalization retains original units and
    /// high-Q sensitivity. Input, viscous and boundary work are differentiated
    /// separately, including the direct change of viscous coefficients.
    /// Budgets precede allocation; refusal/cancellation leaves the owner intact.
    pub fn solve_with_parameter_derivative(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_derivative: Option<&[C64]>,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicParameterResponse, HarmonicResponseError> {
        self.solve_with_total_parameter_derivative(
            force,
            direction,
            force_derivative,
            None,
            0.0,
            controls,
            cx,
        )
    }

    /// Differentiate motion and power when physical attachments and loads change.
    ///
    /// The mechanical and source directions have the same meaning as
    /// `solve_with_parameter_derivative`. This loaded-only extension also takes
    /// complete authored derivatives of the real n*p attachment map B and the
    /// complex p*p impedance Z, at the SAME physical frequency and port units.
    /// Neither derivative is required to be passive or reciprocal: it is a
    /// signed tangent of the admitted load, not a second load to factor.
    ///
    /// With y=B^T v and e=Z y, the same-factor mechanical RHS adds
    /// `-B_theta e - B (Z_theta y + Z B_theta^T v)`. Both occurrences of B
    /// remain. Published port tangents are `y_theta=B^T v_theta+B_theta^T v`
    /// and `e_theta=Z y_theta+Z_theta y`, with independently differentiated
    /// signed boundary work. Full mutual and nonreciprocal entries are retained.
    ///
    /// Factored normalized port actions take O(n*p+p*p) work and O(p) added
    /// memory; no n*n load derivative or second factor is constructed. Their
    /// uncancelled envelopes join the original-coordinate equation audit.
    /// The physical coordinate basis, frequency and port count remain fixed;
    /// this does not infer geometry/material derivatives or a moving-basis map.
    pub fn solve_with_parameter_and_port_derivative(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_derivative: Option<&[C64]>,
        port_direction: HarmonicPortParameterDirection<'_>,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicParameterResponse, HarmonicResponseError> {
        self.solve_with_total_parameter_derivative(
            force,
            direction,
            force_derivative,
            Some(port_direction),
            0.0,
            controls,
            cx,
        )
    }

    /// Total physical parameter tangent, including a changing drive frequency.
    ///
    /// Supply the TOTAL authored M/C/K/F and optional B/Z derivatives along
    /// theta, including their constitutive dispersion when frequency changes.
    /// `angular_frequency_derivative_rad_s` is d(omega)/d(theta), finite and
    /// signed. None for the port direction fixes B and Z, but does not remove
    /// the frequency derivative of their velocity-load action.
    ///
    /// The mechanical RHS also includes
    /// `omega_theta*(2*omega*M*q + i*C*q + i*B*Z*B^T*q)`; published velocity
    /// is `-i*(omega*q_theta+omega_theta*q)`, including at DC. These are the
    /// explicit kinematic frequency terms, not substitutes for authored
    /// derivatives of dispersive stiffness, damping, mass or impedance.
    /// The same cached factor, normalized products and physical power audit
    /// serve fixed-frequency, load-changing and dispersive total tangents.
    #[allow(clippy::too_many_arguments)] // one physical direction with independent source/load/frequency data
    pub fn solve_with_total_parameter_derivative(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_derivative: Option<&[C64]>,
        port_direction: Option<HarmonicPortParameterDirection<'_>>,
        angular_frequency_derivative_rad_s: f64,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicParameterResponse, HarmonicResponseError> {
        poll(cx)?;
        self.admit_parameter(
            direction,
            force_derivative,
            port_direction,
            angular_frequency_derivative_rad_s,
            controls,
            cx,
        )?;
        let response = self.solve(force, cx)?;
        let port_action = load::direct_action(
            self,
            &response,
            port_direction,
            angular_frequency_derivative_rad_s,
            cx,
        )?;
        let mut rhs_exponent = i32::MIN;
        self.parameter_terms(
            &response,
            direction,
            force_derivative,
            angular_frequency_derivative_rad_s,
            cx,
            |_, _, term| {
                if term.fraction.hi != 0.0 {
                    rhs_exponent = rhs_exponent.max(term.exponent);
                }
                Ok(())
            },
        )?;
        if let Some(action) = &port_action {
            action.visit_rhs(self, &response, cx, |_, _, term, envelope| {
                for value in [term, envelope] {
                    if value.fraction.hi != 0.0 {
                        rhs_exponent = rhs_exponent.max(value.exponent);
                    }
                }
                Ok(())
            })?;
        }
        // Four binary fractions per term, six terms per matrix column, and
        // both quadratures: leave room for their uncancelled row envelope too.
        let columns = self
            .n
            .checked_add(port_action.as_ref().map_or(0, |action| action.ports()))
            .ok_or_else(|| invalid("parameter source extent overflow"))?;
        rhs_exponent = if rhs_exponent == i32::MIN {
            0
        } else {
            rhs_exponent + columns.next_power_of_two().ilog2() as i32 + 8
        };
        let mut sums = vec![[Dd::ZERO; 2]; self.n];
        let mut envelopes = vec![Dd::ZERO; self.n];
        self.parameter_terms(
            &response,
            direction,
            force_derivative,
            angular_frequency_derivative_rad_s,
            cx,
            |row, axis, term| {
                let value = term.scaled(rhs_exponent)?;
                sums[row][axis] = sums[row][axis] + value;
                envelopes[row] = envelopes[row] + value.abs();
                Ok(())
            },
        )?;
        if let Some(action) = &port_action {
            action.visit_rhs(self, &response, cx, |row, axis, term, envelope| {
                sums[row][axis] = sums[row][axis] + term.scaled(rhs_exponent)?;
                envelopes[row] = envelopes[row] + envelope.scaled(rhs_exponent)?;
                Ok(())
            })?;
        }
        let rhs: Vec<_> = sums
            .iter()
            .map(|v| C64::new(v[0].to_f64(), v[1].to_f64()))
            .collect();
        let (tangent, solve_residual) = self.solve_normalized(&rhs, cx)?;
        let mut displacement_derivative = Vec::with_capacity(self.n);
        let mut velocity_derivative = Vec::with_capacity(self.n);
        for (i, &value) in tangent.iter().enumerate() {
            poll(cx)?;
            let exponent = rhs_exponent - self.coordinate_exponents[i];
            let dq = scale_complex(value, exponent - 2 * self.time_exponent)?;
            displacement_derivative.push(dq);
            let dv = if angular_frequency_derivative_rad_s == 0.0 {
                scale_complex(
                    C64::new(value.im, -value.re).scale(self.scaled_frequency),
                    exponent - self.time_exponent,
                )?
            } else {
                let q = response.displacement[i];
                let mut real = PowerSum::default();
                let mut imaginary = PowerSum::default();
                for (rate, state) in [
                    (self.angular_frequency_rad_s, dq),
                    (angular_frequency_derivative_rad_s, q),
                ] {
                    real.add(Term::product([rate, state.im, 1.0, 1.0], 0))?;
                    imaginary.add(Term::product([-rate, state.re, 1.0, 1.0], 0))?;
                }
                C64::new(real.value()?, imaginary.value()?)
            };
            velocity_derivative.push(dv);
        }
        let residual =
            self.parameter_residual(&displacement_derivative, &rhs, &envelopes, rhs_exponent, cx)?;
        let mut input = PowerSum::default();
        let mut viscous = PowerSum::default();
        for i in 0..self.n {
            poll(cx)?;
            let df = force_derivative.map_or(C64::ZERO, |v| v[i]);
            input.half_inner(df, response.velocity[i], 1.0)?;
            input.half_inner(force[i], velocity_derivative[i], 1.0)?;
            for j in 0..self.n {
                let c = if self.layout.columns(i, self.n).contains(&j) {
                    self.viscous_entry(i, j)
                } else {
                    0.0
                };
                let dc = symmetric_entry(direction.damping, self.n, i, j);
                // Symmetric parts are formed BEFORE contraction: a changed
                // lossless gyroscope affects motion without manufacturing heat.
                viscous.half_inner(velocity_derivative[i], response.velocity[j], c)?;
                viscous.half_inner(response.velocity[i], velocity_derivative[j], c)?;
                viscous.half_inner(response.velocity[i], response.velocity[j], dc)?;
            }
        }
        let port_load = load::observe(self, &response, &velocity_derivative, port_direction, cx)?;
        let input = input.value()?;
        let viscous = viscous.value()?;
        let load = port_load
            .as_ref()
            .map_or(0.0, |v| v.average_absorbed_power_derivative);
        let mut balance = PowerSum::default();
        for value in [input, -viscous, -load] {
            balance.add(Term::product([value, 1.0, 1.0, 1.0], 0))?;
        }
        let balance = balance.value()?;
        if ![input, viscous, balance].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "parameter derivative of physical work is not representable",
            ));
        }
        poll(cx)?;
        Ok(HarmonicParameterResponse {
            response,
            displacement_derivative,
            velocity_derivative,
            port_load,
            average_input_power_derivative: input,
            average_viscous_power_derivative: viscous,
            power_balance_residual_derivative: balance,
            relative_derivative_residual: solve_residual.max(residual),
        })
    }

    fn admit_parameter(
        &self,
        direction: QuadraticParameterDirection<'_>,
        force: Option<&[C64]>,
        port_direction: Option<HarmonicPortParameterDirection<'_>>,
        frequency_derivative: f64,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        let square = self
            .n
            .checked_mul(self.n)
            .ok_or(invalid("parameter matrix extent overflow"))?;
        if !frequency_derivative.is_finite()
            || [direction.mass, direction.damping, direction.stiffness]
                .iter()
                .any(|v| v.len() != square)
            || force.is_some_and(|v| v.len() != self.n || v.iter().any(|&z| !finite(z)))
        {
            return Err(invalid(
                "complete parameter matrices and finite force/frequency directions required",
            ));
        }
        let n = self.n as u128;
        let p = self.port_load.as_ref().map_or(0, |v| v.ports) as u128;
        let k = self
            .compact_port_factor
            .as_ref()
            .map_or(0, |v| v.dimension()) as u128;
        let port_action =
            port_direction.is_some() || frequency_derivative != 0.0 && self.port_load.is_some();
        let additional_port_bytes = if port_action { 192 * p } else { 0 };
        let bytes = usize::try_from(512 * n + 384 * p + additional_port_bytes)
            .map_err(|_| invalid("parameter payload overflow"))?;
        if bytes > controls.maximum_additional_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "parameter tangent additional payload bytes",
                requested: bytes,
                cap: controls.maximum_additional_payload_bytes,
            });
        }
        let work = 256 * (n * n + n * p + p * p)
            + 256 * (n + p + k)
            + if k == 0 { 0 } else { 64 * (n * n * p + k * k) }
            + if port_action {
                256 * (n * p + p * p) + 256 * p
            } else {
                0
            }
            + if frequency_derivative != 0.0 {
                192 * n * n + 256 * n
            } else {
                0
            };
        let requested = u64::try_from(work).map_err(|_| invalid("parameter work overflow"))?;
        if requested > controls.maximum_work {
            return Err(HarmonicResponseError::WorkBudget {
                requested,
                cap: controls.maximum_work,
            });
        }
        load::admit_direction(self, port_direction, cx)?;
        for i in 0..self.n {
            poll(cx)?;
            for j in 0..self.n {
                let a = i * self.n + j;
                let b = j * self.n + i;
                if ![
                    direction.mass[a],
                    direction.damping[a],
                    direction.stiffness[a],
                ]
                .into_iter()
                .all(f64::is_finite)
                    || direction.mass[a] != direction.mass[b]
                    || direction.stiffness[a] != direction.stiffness[b]
                {
                    return Err(invalid(
                        "finite parameter directions with symmetric mass and stiffness required",
                    ));
                }
            }
        }
        Ok(())
    }

    fn parameter_terms(
        &self,
        response: &HarmonicResponse,
        direction: QuadraticParameterDirection<'_>,
        force: Option<&[C64]>,
        frequency_derivative: f64,
        cx: &Cx<'_>,
        mut visit: impl FnMut(usize, usize, Term) -> Result<(), HarmonicResponseError>,
    ) -> Result<(), HarmonicResponseError> {
        let w = self.angular_frequency_rad_s;
        for i in 0..self.n {
            poll(cx)?;
            let shift = -self.coordinate_exponents[i];
            let f = force.map_or(C64::ZERO, |v| v[i]);
            for (axis, x) in [f.re, f.im].into_iter().enumerate() {
                visit(i, axis, Term::product([x, 1.0, 1.0, 1.0], shift))?;
            }
            for j in 0..self.n {
                let at = i * self.n + j;
                let q = response.displacement[j];
                let v = response.velocity[j];
                for (axis, (q, v)) in [(q.re, v.re), (q.im, v.im)].into_iter().enumerate() {
                    visit(
                        i,
                        axis,
                        Term::product([-direction.stiffness[at], q, 1.0, 1.0], shift),
                    )?;
                    visit(i, axis, Term::product([direction.mass[at], q, w, w], shift))?;
                    visit(
                        i,
                        axis,
                        Term::product([-direction.damping[at], v, 1.0, 1.0], shift),
                    )?;
                }
            }
            if frequency_derivative != 0.0 {
                for j in self.layout.columns(i, self.n) {
                    let at = self.layout.index(i, j, self.n);
                    let q = response.displacement[j];
                    // M_hat=M/2^(coordinate_i+coordinate_j). Recover the
                    // original 2*omega*M*q action without materializing M
                    // or overflowing its products in the original units.
                    for (axis, value) in [q.re, q.im].into_iter().enumerate() {
                        visit(
                            i,
                            axis,
                            Term::product(
                                [self.scaled_mass[at], value, w, frequency_derivative],
                                self.coordinate_exponents[j] + 1,
                            ),
                        )?;
                    }
                    visit(
                        i,
                        0,
                        Term::product([-frequency_derivative, self.damping[at], q.im, 1.0], shift),
                    )?;
                    visit(
                        i,
                        1,
                        Term::product([frequency_derivative, self.damping[at], q.re, 1.0], shift),
                    )?;
                }
            }
        }
        Ok(())
    }

    fn parameter_residual(
        &self,
        displacement: &[C64],
        rhs: &[C64],
        envelopes: &[Dd],
        rhs_exponent: i32,
        cx: &Cx<'_>,
    ) -> Result<f64, HarmonicResponseError> {
        let mut state = Vec::with_capacity(self.n);
        for (i, &q) in displacement.iter().enumerate() {
            state.push(scale_complex(
                q,
                self.coordinate_exponents[i] + 2 * self.time_exponent - rhs_exponent,
            )?);
        }
        let exponent = state
            .iter()
            .flat_map(|v| [v.re, v.im])
            .filter(|v| *v != 0.0)
            .map(|v| quadratic_fraction(v).1)
            .max()
            .unwrap_or(0)
            .max(0);
        for q in &mut state {
            *q = scale_complex(*q, -exponent)?;
        }
        let mut maximum = 0.0_f64;
        for i in 0..self.n {
            poll(cx)?;
            let mut residual = -scale_complex(rhs[i], -exponent)?;
            let mut envelope = scale(envelopes[i].to_f64(), -exponent)?;
            for j in self.operator_columns(i) {
                let a = self.operator_entry(i, j)?;
                residual = residual + a * state[j];
                envelope += a.abs() * state[j].abs();
            }
            let relative = if envelope == 0.0 {
                0.0
            } else {
                residual.abs() / envelope
            };
            if !finite(residual) || !envelope.is_finite() || !relative.is_finite() {
                return Err(invalid(
                    "physical parameter tangent audit is not representable",
                ));
            }
            maximum = maximum.max(relative);
        }
        if maximum > self.maximum_relative_residual {
            return Err(HarmonicResponseError::Residual {
                measured: maximum,
                limit: self.maximum_relative_residual,
            });
        }
        Ok(maximum)
    }
}
