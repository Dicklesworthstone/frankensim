//! Mechanical work/heat tangents through joint velocity/force observations.
//!
//! QR acts on BOTH nominal velocity and force columns. Its same orthogonal
//! rotations carry the two velocity tangents. They are derivatives of the
//! observed functions in a fixed time basis, not derivatives of a QR factor.
//! Keeping force-only rows is essential: zero nominal velocity can still have
//! nonzero first-order source work. No pivot or square root is differentiated.

use super::gradient::{
    ModalTonePowerDerivative, ModalTonePowerGradientBudget, admit_layout, finite, poll,
};
use super::{
    DampedTransition, MAX_COORDINATES, ModalAcousticState, ModalAcousticTimeError, NODES, Sum,
    TonePower, double_delta, double_force, fill_affine_flow, fill_flow, heat_rotation, invalid,
    scaling_with_envelope, small_delta, small_force,
};
use crate::time::forcing::sensitivity::Parameter;
use crate::time::harmonic::integrals::sensitivity::small_action_derivatives;
use crate::time::harmonic::sensitivity::{affine::small_delta_derivative, small_derivative};
use crate::time::harmonic::tones::{ModalToneBudget, PreparedModalTones};
use fs_math::det;

/// Derivatives of coherent work and viscous heat for one physical mode.
/// Each field holds the initial state, prescribed source quadratures and
/// endpoint knots, source clocks, duration and basis/mass normalization fixed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalTonePowerMechanicalGradient {
    /// Independent k = omega^2 [1/s^2], at fixed physical viscous rate c.
    /// Both power derivatives have units J*s^2.
    pub stiffness: ModalTonePowerDerivative,
    /// Physical c [1/s] in q'' + c*q' + k*q = F, at fixed k. Units J*s.
    /// At c=0 this includes the right derivative integral(v^2 dt).
    pub damping_rate: ModalTonePowerDerivative,
}

/// Separately admitted mechanical power tangents of an immutable tone bank.
///
/// Constant and affine envelopes retain every coherent interference term.
/// Cold preparation reuses the original resolved flows, scalar Frechet
/// directions and envelope directions. Joint observation QR preserves signed
/// work and velocity/tangent cross products without energy subtraction or
/// expanded cancelling quadratic forms. No perturbed solve is used.
///
/// Warm queries allocate nothing, evaluate no transcendentals and use O(n^2)
/// contractions, where n=2+2K or 2+4K for K distinct tones. A trajectory reverse
/// sweep adds the existing motion pullback and propagates state cotangents.
/// For a ratio law c=2*zeta*omega, dJ/domega=2*omega*dJ/dk+2*zeta*dJ/dc and
/// dJ/dzeta=2*omega*dJ/dc; those authored chains are not silently applied.
/// These finite-arithmetic derivatives are estimates, not enclosures.
#[derive(Debug)]
pub struct PreparedModalTonePowerMechanicalGradient<'a> {
    bank: &'a PreparedModalTones,
    factors: Box<[MechanicalPower]>,
}

#[derive(Debug)]
struct MechanicalPower {
    // 2n rows of [V, F, V_K, V_C], each block n columns. K=k*h^2, C=c*h.
    // V and F drive the rotations; the tangent blocks follow unchanged.
    rows: Box<[f64]>,
    root_duration: f64,
}

impl PreparedModalTones {
    /// Cold-prepare coherent stiffness and viscous-rate work/heat derivatives.
    /// The original bank, accepted states and source phases are unchanged.
    ///
    /// Additional payload plus peak scratch is admitted before allocation.
    /// The preparation budget's step cap bounds a complete all-mode query
    /// pass; query_budget separately bounds one mode's fixed scratch and work.
    /// Cancellation is polled between modes, short quadrature rows, doubling
    /// rows and QR pivots. Refusal publishes no partial preparation.
    pub fn prepare_power_mechanical_gradients(
        &self,
        preparation_budget: ModalToneBudget,
        query_budget: ModalTonePowerGradientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedModalTonePowerMechanicalGradient<'_>, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        admit_layout(
            self.unique.len(),
            self.tone_count(),
            self.affine,
            query_budget,
        )?;
        admit(self, preparation_budget, &mut cancelled)?;
        let mut factors = Vec::with_capacity(self.modes.len());
        for mode in &self.modes {
            poll(&mut cancelled)?;
            factors.push(MechanicalPower::prepare(
                &mode.power,
                mode.motion[0].transition.drag,
                &self.unique,
                self.affine,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(PreparedModalTonePowerMechanicalGradient {
            bank: self,
            factors: factors.into_boxed_slice(),
        })
    }
}

impl PreparedModalTonePowerMechanicalGradient<'_> {
    /// Derivatives of complete constant-envelope work and heat for one mode.
    pub fn gradients_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalTonePowerMechanicalGradient, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if self.bank.affine {
            return Err(invalid(
                "mechanical power gradients require both affine endpoints",
            ));
        }
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let (c, s) = self.bank.collect_forces(cosine, sine)?;
        let k = self.bank.unique.len();
        let x = mode.power.coordinates(state, &c[..k], &s[..k], None)?;
        self.factors[mode_index].evaluate(&mode.power, state, &x[..2 + 2 * k], &mut cancelled)
    }

    /// Both endpoint pairs use the same interval-start carrier basis. The
    /// derivative holds BOTH authored endpoint values fixed, including their
    /// effect on all cross-frequency transient and heating contributions.
    #[allow(clippy::too_many_arguments)]
    pub fn gradients_mode_affine(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalTonePowerMechanicalGradient, ModalAcousticTimeError> {
        poll(&mut cancelled)?;
        if !self.bank.affine {
            return Err(invalid(
                "affine mechanical power gradients require affine preparation",
            ));
        }
        let mode = self
            .bank
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let [c, s, dc, ds] =
            self.bank
                .collect_affine(cosine_start, sine_start, cosine_end, sine_end)?;
        let k = self.bank.unique.len();
        let x = mode
            .power
            .coordinates(state, &c[..k], &s[..k], Some((&dc[..k], &ds[..k])))?;
        self.factors[mode_index].evaluate(&mode.power, state, &x[..2 + 4 * k], &mut cancelled)
    }
}

impl MechanicalPower {
    fn prepare(
        base: &TonePower,
        original_drag: f64,
        frequencies: &[f64],
        affine: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ModalAcousticTimeError> {
        let (angle, drag, mut time, scale, doublings) =
            scaling_with_envelope(base.omega, original_drag, frequencies, base.h, affine)?;
        let ratio = (base.omega * base.h) / scale;
        let width = if affine { 4 } else { 2 };
        let n = 2 + width * frequencies.len();
        let mut phases = frequencies.iter().map(|f| f * base.h).collect::<Vec<_>>();
        let mut fraction = 1.0;
        for _ in 0..doublings {
            fraction *= 0.5;
            for phase in &mut phases {
                *phase *= 0.5;
            }
        }
        let mut result = Self {
            rows: vec![0.0; 8 * n * n].into_boxed_slice(),
            root_duration: det::sqrt(base.h),
        };
        let mut product = vec![0.0; 8 * n * n];
        let mut e = vec![0.0; n * n];
        let mut directions = [vec![0.0; 2 * n], vec![0.0; 2 * n]];
        let mut next_directions = [vec![0.0; 2 * n], vec![0.0; 2 * n]];
        let mut forces = vec![[0.0; 4]; frequencies.len()];
        let mut deltas = vec![[0.0; 4]; if affine { frequencies.len() } else { 0 }];
        let mut row = vec![0.0; 4 * n];
        for (node, weight) in NODES {
            poll(cancelled)?;
            let flow = DampedTransition::small(angle * node, drag * node);
            for (force, &phase) in forces.iter_mut().zip(&phases) {
                *force = small_force(angle * node, drag * node, phase * node, false);
            }
            for (delta, &phase) in deltas.iter_mut().zip(&phases) {
                *delta = small_delta(angle * node, drag * node, phase * node);
            }
            if affine {
                fill_affine_flow(
                    &mut e,
                    flow,
                    &forces,
                    &deltas,
                    &phases,
                    node,
                    time * node,
                    ratio,
                    fraction * node,
                );
            } else {
                fill_flow(&mut e, flow, &forces, &phases, node, time * node, ratio);
            }
            small_directions(
                &mut directions,
                angle * node,
                drag * node,
                &phases,
                node,
                time * node,
                fraction * node,
                scale,
                &forces,
                &deltas,
            );
            row.fill(0.0);
            let factor = det::sqrt(fraction) * det::sqrt(weight);
            for j in 0..n {
                row[j] = factor * e[j];
                row[2 * n + j] = factor * directions[0][j];
                row[3 * n + j] = factor * directions[1][j];
            }
            for (tone, &phase) in phases.iter().enumerate() {
                let (sine, cosine) = det::sin_cos(phase * node);
                let j = n + 2 + width * tone;
                row[j] = factor * cosine;
                row[j + 1] = factor * sine;
                if affine {
                    row[j + 2] = factor * fraction * node * cosine;
                    row[j + 3] = factor * fraction * node * sine;
                }
            }
            result.append(&mut row, n, cancelled)?;
        }
        let mut flow = DampedTransition::small(angle, drag);
        for (force, &phase) in forces.iter_mut().zip(&phases) {
            *force = small_force(angle, drag, phase, false);
        }
        for (delta, &phase) in deltas.iter_mut().zip(&phases) {
            *delta = small_delta(angle, drag, phase);
        }
        small_directions(
            &mut directions,
            angle,
            drag,
            &phases,
            1.0,
            time,
            fraction,
            scale,
            &forces,
            &deltas,
        );
        for _ in 0..doublings {
            poll(cancelled)?;
            if affine {
                fill_affine_flow(
                    &mut e, flow, &forces, &deltas, &phases, 1.0, time, ratio, fraction,
                );
            } else {
                fill_flow(&mut e, flow, &forces, &phases, 1.0, time, ratio);
            }
            for (before, after) in result
                .rows
                .chunks_exact(4 * n)
                .zip(product.chunks_exact_mut(4 * n))
            {
                poll(cancelled)?;
                for block in 0..4 {
                    for j in 0..n {
                        let mut value =
                            right_column(&before[block * n..(block + 1) * n], &e, j, n, width);
                        if block >= 2 {
                            let direction = &directions[block - 2];
                            value += before[0] * direction[j] + before[1] * direction[n + j];
                        }
                        after[block * n + j] = value;
                    }
                }
            }
            for row in product.chunks_exact_mut(4 * n) {
                result.append(row, n, cancelled)?;
            }
            // The original source blocks do not depend on k/c: each Frechet
            // action has only two dense rows. Full-interval K,C seeds already
            // contain the short fraction^2/fraction; doubling needs no rescale.
            for p in 0..2 {
                for i in 0..2 {
                    for j in 0..n {
                        next_directions[p][i * n + j] =
                            right_column(&directions[p][i * n..(i + 1) * n], &e, j, n, width)
                                + e[i * n] * directions[p][j]
                                + e[i * n + 1] * directions[p][n + j];
                    }
                }
            }
            core::mem::swap(&mut directions, &mut next_directions);
            for (index, (force, phase)) in forces.iter_mut().zip(&mut phases).enumerate() {
                if affine {
                    deltas[index] = double_delta(deltas[index], *force, flow, *phase);
                }
                *force = double_force(*force, flow, *phase, false);
                *phase *= 2.0;
            }
            flow.double();
            time *= 2.0;
            fraction *= 2.0;
        }
        poll(cancelled)?;
        if result.rows.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "mechanical tone power factors must be representable",
            ));
        }
        Ok(result)
    }

    fn append(
        &mut self,
        row: &mut [f64],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ModalAcousticTimeError> {
        let stride = 4 * n;
        for i in 0..2 * n {
            poll(cancelled)?;
            let Some((norm, c, s)) = heat_rotation(self.rows[i * stride + i], row[i]) else {
                continue;
            };
            self.rows[i * stride + i] = norm;
            for j in i + 1..stride {
                let before = self.rows[i * stride + j];
                self.rows[i * stride + j] = c * before + s * row[j];
                row[j] = -s * before + c * row[j];
            }
        }
        // Remaining tangent-only rows are orthogonal to BOTH nominal
        // observations. Their first-order work/heat contributions vanish and
        // remain zero after right composition with any original source flow.
        Ok(())
    }

    fn evaluate(
        &self,
        base: &TonePower,
        state: ModalAcousticState,
        x: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ModalTonePowerMechanicalGradient, ModalAcousticTimeError> {
        let n = x.len();
        let mut work_k = Sum::default();
        let mut work_c = Sum::default();
        let mut heat_k = Sum::default();
        let mut heat_c = Sum::default();
        for row in self.rows.chunks_exact(4 * n) {
            poll(cancelled)?;
            let mut values = [Sum::default(); 4];
            for block in 0..4 {
                for (coefficient, &coordinate) in row[block * n..(block + 1) * n].iter().zip(x) {
                    values[block].add(coefficient * coordinate);
                }
            }
            let [v, f, vk, vc] = values.map(Sum::value);
            // Kappa is a frozen conditioning chart. Its max() is never
            // differentiated. At fixed physical q, x_1=(sum C-k*q)*h/scale
            // contributes -q*h/scale before the bilinear contractions.
            let vk = (vk * base.h).mul_add(
                base.h,
                -((row[1] * (base.h / base.scale)) * state.displacement_m_sqrt_kg),
            );
            let vc = vc * base.h;
            if f != 0.0 {
                let force = base.work_scale * f;
                work_k.add(force * (base.work_scale * vk));
                work_c.add(force * (base.work_scale * vc));
            }
            let direct_heat = self.root_duration * v;
            heat_c.add(direct_heat * direct_heat);
            if base.heat_scale != 0.0 && v != 0.0 {
                let velocity = base.heat_scale * v;
                let k_term = velocity * (base.heat_scale * vk);
                let c_term = velocity * (base.heat_scale * vc);
                heat_k.add(k_term);
                heat_k.add(k_term);
                heat_c.add(c_term);
                heat_c.add(c_term);
            }
        }
        poll(cancelled)?;
        Ok(ModalTonePowerMechanicalGradient {
            stiffness: finite(ModalTonePowerDerivative {
                input_work: work_k.value(),
                viscous_dissipation: heat_k.value(),
            })?,
            damping_rate: finite(ModalTonePowerDerivative {
                input_work: work_c.value(),
                viscous_dissipation: heat_c.value(),
            })?,
        })
    }
}

// Reuse the scalar residual-chart Frechet owner for the common two mechanical
// columns and every constant source. Affine columns differentiate the SAME
// original fill_affine_flow formulas using existing entire force/envelope
// derivatives; there is no additional time integrator or pole formula.
#[allow(clippy::too_many_arguments)]
fn small_directions(
    directions: &mut [Vec<f64>; 2],
    angle: f64,
    drag: f64,
    phases: &[f64],
    node: f64,
    time: f64,
    fraction: f64,
    scale: f64,
    forces: &[[f64; 4]],
    deltas: &[[f64; 4]],
) {
    let affine = !deltas.is_empty();
    let width = if affine { 4 } else { 2 };
    let n = 2 + width * phases.len();
    for row in directions.iter_mut() {
        row.fill(0.0);
    }
    for (tone, &phase) in phases.iter().enumerate() {
        let phase = phase * node;
        let derivative = small_action_derivatives(angle, drag, phase, time, fraction, scale);
        let column = 2 + width * tone;
        for p in 0..2 {
            for i in 0..2 {
                directions[p][i * n] = derivative[p][i][0];
                directions[p][i * n + 1] = derivative[p][i][1];
                directions[p][i * n + column] = derivative[p][i][2];
                directions[p][i * n + column + 1] = derivative[p][i][3];
            }
            if affine {
                let parameter = if p == 0 {
                    Parameter::Stiffness
                } else {
                    Parameter::Damping
                };
                let f = small_derivative(angle, drag, phase, false, parameter);
                let d = small_delta_derivative(angle, drag, phase, parameter);
                let multiplier = if p == 0 {
                    fraction * fraction
                } else {
                    fraction
                };
                directions[p][column + 2] = time * fraction * multiplier * (f[0] - phase * d[1]);
                directions[p][column + 3] = time * fraction * multiplier * (f[1] + phase * d[0]);
                let extra_c = if p == 1 {
                    forces[tone][0] - phase * deltas[tone][1]
                } else {
                    0.0
                };
                let extra_s = if p == 1 {
                    forces[tone][1] + phase * deltas[tone][0]
                } else {
                    0.0
                };
                directions[p][n + column + 2] = fraction
                    * multiplier
                    * (f[2] + drag * f[0] - phase * (d[3] + drag * d[1]) + extra_c);
                directions[p][n + column + 3] = fraction
                    * multiplier
                    * (f[3] + drag * f[1] + phase * (d[2] + drag * d[0]) + extra_s);
            }
        }
    }
}

// Original two-dense-row plus independent source-block right action, applied
// to a joint observation row. Only the storage stride differs from TonePower.
fn right_column(row: &[f64], e: &[f64], j: usize, n: usize, width: usize) -> f64 {
    let mut value = row[0] * e[j] + row[1] * e[n + j];
    if j >= 2 {
        let block = 2 + width * ((j - 2) / width);
        for k in block..block + width {
            value += row[k] * e[k * n + j];
        }
    }
    value
}

fn admit(
    bank: &PreparedModalTones,
    budget: ModalToneBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ModalAcousticTimeError> {
    let k = bank.unique.len();
    let n = (2 + (if bank.affine { 4 } else { 2 }) * k) as u128;
    let square = n * n;
    let modes = bank.modes.len() as u128;
    // Owned 2n-by-4n factors plus one product of that shape, original n-by-n
    // flow, two pairs of two-row Frechet buffers, source/row vectors and fixed
    // scalar-helper scratch. Existing bank and allocator overhead are excluded.
    let storage = modes * (8 * square * 8 + core::mem::size_of::<MechanicalPower>() as u128)
        + (9 * square + 32 * n + 16 * k as u128) * 8
        + 8192;
    let query_work = 64 * square + 256 * bank.tone_count() as u128;
    for (what, value, limit) in [
        (
            "mechanical power input tone count",
            bank.tone_count() as u128,
            budget.maximum_tones,
        ),
        (
            "mechanical power coefficient and preparation storage bytes",
            storage,
            budget.maximum_storage_bytes,
        ),
        (
            "mechanical power all-mode query work score",
            modes * query_work,
            budget.maximum_step_work,
        ),
    ] {
        cap(what, value, limit)?;
    }
    let mut work = 0u128;
    for mode in &bank.modes {
        poll(cancelled)?;
        let (_, _, _, _, doublings) = scaling_with_envelope(
            mode.power.omega,
            mode.motion[0].transition.drag,
            &bank.unique,
            bank.duration_s,
            bank.affine,
        )?;
        // Algebraic loop score, not elapsed time. Joint QR has 2n pivots and
        // four n-column observation blocks; source derivatives reuse fixed
        // 24-term scalar recurrences on eight short Gaussian nodes.
        work += 512 * square
            + 65_536 * k as u128 * if bank.affine { 2 } else { 1 }
            + doublings as u128 * (64 * n * square + 256 * square);
        cap(
            "mechanical power preparation work score",
            work,
            budget.maximum_preparation_work,
        )?;
    }
    // Keep the physical-query fixed-array admission tied to its actual owner.
    debug_assert!((n as usize) <= MAX_COORDINATES);
    Ok(())
}

fn cap(what: &'static str, value: u128, limit: usize) -> Result<(), ModalAcousticTimeError> {
    if value > limit as u128 {
        Err(ModalAcousticTimeError::BudgetExceeded {
            what,
            value: value as f64,
            limit: limit as f64,
        })
    } else {
        Ok(())
    }
}
