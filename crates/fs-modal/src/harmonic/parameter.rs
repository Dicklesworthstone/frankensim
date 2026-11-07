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

/// Additional resources for one physical parameter direction and its response.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicParameterControls {
    /// Conservative added numerical payload: 512*n + 384*p bytes, excluding
    /// the prepared factor, borrowed direction/source arrays and allocator metadata.
    pub maximum_additional_payload_bytes: usize,
    /// Work score: 256*(n*n+n*p+p*p)+256*(n+p+k), plus 64*(n*n*p+k*k)
    /// for a compact loaded factor of dimension k. Includes both solves,
    /// audits and binary-normalized physical power contractions.
    pub maximum_work: u64,
}

/// Derivatives of the original physical load variables per unit parameter.
#[derive(Clone, Debug)]
pub struct HarmonicPortParameterResponse {
    /// B^T v_theta. The physical attachment map is held fixed.
    pub flow_derivative: Vec<C64>,
    /// Z flow_theta. The supplied load impedance is held fixed.
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
    /// -i omega dq/dtheta, since the physical forcing frequency is fixed.
    pub velocity_derivative: Vec<C64>,
    /// Load flow/effort derivatives when a fixed load is attached.
    pub port_load: Option<HarmonicPortParameterResponse>,
    /// d Re(F^H v)/2, including the authored force derivative.
    pub average_input_power_derivative: f64,
    /// d Re(v^H D v)/2, including D_theta=(C_theta+C_theta^T)/2.
    pub average_viscous_power_derivative: f64,
    /// Input minus independent viscous and fixed-load work derivatives.
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
        poll(cx)?;
        self.admit_parameter(direction, force_derivative, controls, cx)?;
        let response = self.solve(force, cx)?;
        let mut rhs_exponent = i32::MIN;
        self.parameter_terms(&response, direction, force_derivative, cx, |_, _, term| {
            if term.fraction.hi != 0.0 {
                rhs_exponent = rhs_exponent.max(term.exponent);
            }
            Ok(())
        })?;
        // Four binary fractions per term, six terms per matrix column, and
        // both quadratures: leave room for their uncancelled row envelope too.
        rhs_exponent = if rhs_exponent == i32::MIN {
            0
        } else {
            rhs_exponent + self.n.next_power_of_two().ilog2() as i32 + 8
        };
        let mut sums = vec![[Dd::ZERO; 2]; self.n];
        let mut envelopes = vec![Dd::ZERO; self.n];
        self.parameter_terms(
            &response,
            direction,
            force_derivative,
            cx,
            |row, axis, term| {
                let value = term.scaled(rhs_exponent)?;
                sums[row][axis] = sums[row][axis] + value;
                envelopes[row] = envelopes[row] + value.abs();
                Ok(())
            },
        )?;
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
            displacement_derivative.push(scale_complex(value, exponent - 2 * self.time_exponent)?);
            velocity_derivative.push(scale_complex(
                C64::new(value.im, -value.re).scale(self.scaled_frequency),
                exponent - self.time_exponent,
            )?);
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
        let port_load = self.parameter_port_response(&response, &velocity_derivative, cx)?;
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
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        let square = self
            .n
            .checked_mul(self.n)
            .ok_or(invalid("parameter matrix extent overflow"))?;
        if [direction.mass, direction.damping, direction.stiffness]
            .iter()
            .any(|v| v.len() != square)
            || force.is_some_and(|v| v.len() != self.n || v.iter().any(|&z| !finite(z)))
        {
            return Err(invalid(
                "complete parameter matrices and finite force direction required",
            ));
        }
        let n = self.n as u128;
        let p = self.port_load.as_ref().map_or(0, |v| v.ports) as u128;
        let k = self
            .compact_port_factor
            .as_ref()
            .map_or(0, |v| v.dimension()) as u128;
        let bytes = usize::try_from(512 * n + 384 * p)
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
            + if k == 0 { 0 } else { 64 * (n * n * p + k * k) };
        let requested = u64::try_from(work).map_err(|_| invalid("parameter work overflow"))?;
        if requested > controls.maximum_work {
            return Err(HarmonicResponseError::WorkBudget {
                requested,
                cap: controls.maximum_work,
            });
        }
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

    fn parameter_port_response(
        &self,
        response: &HarmonicResponse,
        velocity: &[C64],
        cx: &Cx<'_>,
    ) -> Result<Option<HarmonicPortParameterResponse>, HarmonicResponseError> {
        let (Some(load), Some(primal)) = (&self.port_load, &response.port_load) else {
            return Ok(None);
        };
        let p = load.ports;
        let mut flow = vec![C64::ZERO; p];
        let mut effort = vec![C64::ZERO; p];
        for a in 0..p {
            poll(cx)?;
            for (i, &v) in velocity.iter().enumerate() {
                flow[a] = flow[a] + v.scale(load.force_map[i * p + a]);
            }
        }
        if flow.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "fixed-load parameter flow tangent is not representable",
            ));
        }
        let mut power = PowerSum::default();
        for a in 0..p {
            poll(cx)?;
            for b in 0..p {
                effort[a] = effort[a] + load.impedance[a * p + b] * flow[b];
            }
            if !finite(effort[a]) {
                return Err(invalid(
                    "fixed-load parameter effort tangent is not representable",
                ));
            }
            power.half_inner(flow[a], primal.effort[a], 1.0)?;
            power.half_inner(primal.flow[a], effort[a], 1.0)?;
        }
        let power = power.value()?;
        if !power.is_finite() || flow.iter().chain(&effort).any(|&v| !finite(v)) {
            return Err(invalid("fixed-load parameter tangent is not representable"));
        }
        Ok(Some(HarmonicPortParameterResponse {
            flow_derivative: flow,
            effort_derivative: effort,
            average_absorbed_power_derivative: power,
        }))
    }
}
