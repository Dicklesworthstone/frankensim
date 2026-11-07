//! Prescribed support phasors on the complete original mechanical pencil.
//!
//! A = K - omega^2 M - i omega C. The free solve is
//! A_ff q_f = F_f - A_fc q_c, and the support force ON the structure is
//! R_c = (A q - F)_c. Keeping the cross inertia and velocity terms is essential
//! for base excitation; a static stiffness load does not describe a shaking
//! support. Full-source admission and binary scaling belong to the same owner
//! as the unconstrained harmonic solve.

use super::{
    HarmonicResponseControls, HarmonicResponseError, MatrixLayout, PreparedHarmonicResponse,
    admit_size, finite, invalid, poll, scale, scale_complex, symmetric_entry,
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};

/// Original-coordinate particular motion and its independent cycle-work account.
#[derive(Clone, Debug)]
pub struct SupportedHarmonicResponse {
    /// Full peak generalized displacement, including prescribed supports.
    pub displacement: Vec<C64>,
    /// Full peak velocity, `-i omega displacement`.
    pub velocity: Vec<C64>,
    /// Force exerted BY each support ON the structure, in `support_dofs()` order.
    /// Applied forces at supported coordinates are subtracted exactly once.
    pub support_reaction: Vec<C64>,
    /// Maximum componentwise backward error in the original free rows,
    /// `|Aq-F| / (|A||q|+|F|)`, before support elimination.
    pub relative_residual: f64,
    /// Cycle average `Re(F^H v)/2`, including applied support-coordinate forces.
    pub average_input_power_w: f64,
    /// Cycle average `Re(R_c^H v_c)/2`; signed work supplied by the supports.
    pub average_support_power_w: f64,
    /// Cycle average `Re(v^H D v)/2`, D=(C+C^T)/2, on ALL original coordinates.
    /// Cross-support loss is retained; skew gyroscopic coupling creates no heat.
    pub average_viscous_power_w: f64,
    /// Applied plus support power minus viscous power [W], accumulated without
    /// assigning a missing work term from the balance. Signed roundoff remains.
    pub power_balance_residual_w: f64,
}

/// Reusable supported response at one physical angular frequency.
///
/// The full M/C/K passes the usual harmonic admission: M is SPD, K and the
/// symmetric part of C are numerically PSD, and skew C is retained. The full
/// operator need NOT be invertible: only its free principal block is factored.
/// Thus rigid-body stiffness modes do not prevent supported static response.
/// No penalty stiffness, artificial inertia, damping, or coordinate truncation
/// is introduced. A singular free block still refuses.
///
/// Supports prescribe generalized displacements, including rotations when
/// present; work uses their conjugate forces or moments. Empty and all-support
/// maps are admitted. This is a particular periodic response, not a startup
/// transient, support-actuator dynamics, or a spatial/material certificate.
#[derive(Clone, Debug)]
pub struct PreparedSupportedHarmonicResponse {
    original: PreparedHarmonicResponse,
    free: Option<PreparedHarmonicResponse>,
    free_dofs: Vec<usize>,
    support_dofs: Vec<usize>,
}

impl PreparedSupportedHarmonicResponse {
    /// Prepare complete real row-major M/C/K and a strictly increasing support
    /// index map. All indices refer to the original `n` coordinates. Frequencies
    /// are finite nonnegative [rad/s]; DC is the supported static problem.
    ///
    /// `maximum_coordinates` bounds the FULL source. Numerical payload admission
    /// is `320*n*n + 512*n` bytes, covering full-source admission, both retained
    /// operators, the free LU, and one complete solve. The original size also
    /// passes the shared dense work cap. Caller inputs and allocator metadata are
    /// excluded. Cancellation is polled during assembly, elimination and solves.
    #[allow(clippy::too_many_arguments)] // Full pencil, prescribed map, physical frequency and controls.
    pub fn try_new(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        support_dofs: &[usize],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_size(n, controls)?;
        if support_dofs.iter().any(|&i| i >= n)
            || support_dofs.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "support indices must be strictly increasing original rows",
            ));
        }
        let bytes = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(320))
            .and_then(|x| n.checked_mul(512)?.checked_add(x))
            .ok_or(invalid("supported harmonic payload extent overflow"))?;
        if bytes > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "supported numerical payload bytes",
                requested: bytes,
                cap: controls.maximum_payload_bytes,
            });
        }
        let original = PreparedHarmonicResponse::prepare(
            mass,
            damping,
            stiffness,
            n,
            angular_frequency_rad_s,
            controls,
            cx,
            None,
            None,
            MatrixLayout::Dense,
            false,
        )?;
        let free_dofs: Vec<_> = (0..n)
            .filter(|i| support_dofs.binary_search(i).is_err())
            .collect();
        let free = if free_dofs.is_empty() {
            None
        } else {
            // Principal restriction of the ALREADY admitted/scaled full pencil.
            // Keep exactly the same coordinate and time scales so the support
            // source, solution residual, and reactions share one original A.
            let nf = free_dofs.len();
            let mut operator = Vec::with_capacity(nf * nf);
            let mut scaled_mass = Vec::with_capacity(nf * nf);
            let mut scaled_damping = Vec::with_capacity(nf * nf);
            let mut damping = Vec::with_capacity(nf * nf);
            for &i in &free_dofs {
                poll(cx)?;
                for &j in &free_dofs {
                    let k = i * n + j;
                    operator.push(original.operator[k]);
                    scaled_mass.push(original.scaled_mass[k]);
                    scaled_damping.push(original.scaled_damping[k]);
                    damping.push(original.damping[k]);
                }
            }
            let mut free = PreparedHarmonicResponse {
                n: nf,
                layout: MatrixLayout::Dense,
                angular_frequency_rad_s,
                scaled_frequency: original.scaled_frequency,
                time_exponent: original.time_exponent,
                coordinate_exponents: free_dofs
                    .iter()
                    .map(|&i| original.coordinate_exponents[i])
                    .collect(),
                operator,
                scaled_mass,
                scaled_damping,
                damping,
                factor: None,
                compact_port_factor: None,
                port_load: None,
                maximum_relative_residual: controls.maximum_relative_residual,
            };
            free.factorize(cx)?;
            Some(free)
        };
        poll(cx)?;
        Ok(Self {
            original,
            free,
            free_dofs,
            support_dofs: support_dofs.to_vec(),
        })
    }

    /// Unprescribed original coordinate rows, in increasing order.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Prescribed original coordinate rows; inputs and reactions use this order.
    pub fn support_dofs(&self) -> &[usize] {
        &self.support_dofs
    }

    /// Physical angular frequency [rad/s], including DC.
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.original.angular_frequency_rad_s
    }

    /// Solve one full applied-force vector and one displacement phasor per
    /// support, under `exp(-i omega t)`. Relative phases remain intact. A repeated
    /// load case reuses the free factor; source matrices/state are never changed.
    ///
    /// Reaction and work contractions use compensated double-double arithmetic
    /// on the shared binary-scaled original pencil. Power is restored in source
    /// units only after contraction. No finite-response or near-null relative
    /// guarantee beyond the caller's backward-residual ceiling is implied.
    pub fn solve(
        &self,
        applied_force: &[C64],
        support_displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<SupportedHarmonicResponse, HarmonicResponseError> {
        poll(cx)?;
        let p = &self.original;
        let n = p.n;
        if applied_force.len() != n
            || support_displacement.len() != self.support_dofs.len()
            || applied_force
                .iter()
                .chain(support_displacement)
                .any(|&z| !finite(z))
        {
            return Err(invalid(
                "full finite forces and one displacement per support required",
            ));
        }
        // q = 2^e D y; A_hat = D A D/a^2, b = D F/(2^e a^2),
        // D_ii=2^-coordinate_i and a=2^time. Normalize BOTH independent inputs;
        // a small prescribed motion must survive cancellation of applied forces.
        let mut exponent = i32::MIN;
        for (i, &f) in applied_force.iter().enumerate() {
            for value in [f.re, f.im] {
                if value != 0.0 {
                    exponent = exponent.max(
                        quadratic_fraction(value).1
                            - p.coordinate_exponents[i]
                            - 2 * p.time_exponent,
                    );
                }
            }
        }
        for (&i, &q) in self.support_dofs.iter().zip(support_displacement) {
            for value in [q.re, q.im] {
                if value != 0.0 {
                    exponent =
                        exponent.max(quadratic_fraction(value).1 + p.coordinate_exponents[i]);
                }
            }
        }
        if exponent == i32::MIN {
            exponent = 0;
        }
        let mut rhs = Vec::with_capacity(n);
        for (i, &f) in applied_force.iter().enumerate() {
            rhs.push(scale_complex(
                f,
                -p.coordinate_exponents[i] - 2 * p.time_exponent - exponent,
            )?);
        }
        let mut motion = vec![C64::ZERO; n];
        for (&i, &q) in self.support_dofs.iter().zip(support_displacement) {
            motion[i] = scale_complex(q, p.coordinate_exponents[i] - exponent)?;
        }
        if let Some(free) = &self.free {
            let mut free_rhs = Vec::with_capacity(self.free_dofs.len());
            for &i in &self.free_dofs {
                poll(cx)?;
                let mut sum = parts(rhs[i]);
                for &j in &self.support_dofs {
                    let term = multiply(p.operator[i * n + j], motion[j]);
                    sum[0] = sum[0] - term[0];
                    sum[1] = sum[1] - term[1];
                }
                free_rhs.push(complex(sum)?);
            }
            let (solution, _) = free.solve_normalized(&free_rhs, cx)?;
            for (&i, value) in self.free_dofs.iter().zip(solution) {
                motion[i] = value;
            }
        }
        let mut displacement = Vec::with_capacity(n);
        let mut velocity = Vec::with_capacity(n);
        for (i, &q) in motion.iter().enumerate() {
            poll(cx)?;
            displacement.push(scale_complex(q, exponent - p.coordinate_exponents[i])?);
            velocity.push(scale_complex(
                rotate(q).scale(p.scaled_frequency),
                exponent - p.coordinate_exponents[i] + p.time_exponent,
            )?);
        }
        let mut reactions = Vec::with_capacity(self.support_dofs.len());
        let mut relative_residual = 0.0_f64;
        let mut input = Dd::ZERO;
        let mut support = Dd::ZERO;
        let mut heat = Dd::ZERO;
        for i in 0..n {
            poll(cx)?;
            let mut residual = [-Dd::from_f64(rhs[i].re), -Dd::from_f64(rhs[i].im)];
            let mut denominator = rhs[i].abs();
            let mut viscous = [Dd::ZERO; 2];
            let v = rotate(motion[i]).scale(p.scaled_frequency);
            for j in 0..n {
                let coefficient = p.operator[i * n + j];
                let term = multiply(coefficient, motion[j]);
                residual[0] = residual[0] + term[0];
                residual[1] = residual[1] + term[1];
                denominator += coefficient.abs() * motion[j].abs();
                let d = Dd::from_f64(symmetric_entry(&p.scaled_damping, n, i, j));
                let vj = rotate(motion[j]).scale(p.scaled_frequency);
                viscous[0] = viscous[0] + d * Dd::from_f64(vj.re);
                viscous[1] = viscous[1] + d * Dd::from_f64(vj.im);
            }
            let row_residual = complex(residual)?;
            if self.support_dofs.binary_search(&i).is_ok() {
                reactions.push(scale_complex(
                    row_residual,
                    exponent + 2 * p.time_exponent + p.coordinate_exponents[i],
                )?);
                support = support + dot(residual, v);
            } else {
                let error = if denominator == 0.0 {
                    0.0
                } else {
                    row_residual.abs() / denominator
                };
                if !denominator.is_finite() || !error.is_finite() {
                    return Err(invalid("supported harmonic residual is not representable"));
                }
                relative_residual = relative_residual.max(error);
            }
            input = input + dot(parts(rhs[i]), v);
            heat = heat + dot(viscous, v);
        }
        if relative_residual > p.maximum_relative_residual {
            return Err(HarmonicResponseError::Residual {
                measured: relative_residual,
                limit: p.maximum_relative_residual,
            });
        }
        // All contractions use the same a^3 2^(2e) power scale. Apply the
        // peak-phasor factor 1/2 by an exact binary exponent decrement.
        let power_exponent = 2 * exponent + 3 * p.time_exponent - 1;
        let result = SupportedHarmonicResponse {
            displacement,
            velocity,
            support_reaction: reactions,
            relative_residual,
            average_input_power_w: restore(input, power_exponent)?,
            average_support_power_w: restore(support, power_exponent)?,
            average_viscous_power_w: restore(heat, power_exponent)?,
            power_balance_residual_w: restore(input + support - heat, power_exponent)?,
        };
        poll(cx)?;
        Ok(result)
    }
}

fn rotate(z: C64) -> C64 {
    C64::new(z.im, -z.re)
}
fn parts(z: C64) -> [Dd; 2] {
    [Dd::from_f64(z.re), Dd::from_f64(z.im)]
}
fn multiply(a: C64, b: C64) -> [Dd; 2] {
    let (a, b) = (parts(a), parts(b));
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn dot(a: [Dd; 2], b: C64) -> Dd {
    a[0] * Dd::from_f64(b.re) + a[1] * Dd::from_f64(b.im)
}
fn number(x: Dd) -> Result<f64, HarmonicResponseError> {
    let value = x.to_f64();
    if !x.hi.is_finite() || !x.lo.is_finite() || !value.is_finite() {
        Err(invalid(
            "supported harmonic contraction is not representable",
        ))
    } else {
        Ok(value)
    }
}
fn complex(x: [Dd; 2]) -> Result<C64, HarmonicResponseError> {
    Ok(C64::new(number(x[0])?, number(x[1])?))
}
fn restore(x: Dd, exponent: i32) -> Result<f64, HarmonicResponseError> {
    scale(number(x)?, exponent)
}
