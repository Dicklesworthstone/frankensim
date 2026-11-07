//! Independent multi-shaft relations on the original finite rotor inventory.
//!
//! Constant homogeneous rows A enforce A s = 0 on relative shaft rates s.
//! With D = diag(J_i), factor D^(-1/2) A^T = U Sigma V^T after row scaling.
//! The columns of U are orthonormal reaction modes. Write W = D^(1/2) U,
//! E = D^(-1/2) U and k_c = sum_i W_ic a_i. Eliminating only constrained
//! modes gives I_retained = I + sum_c k_c k_c^T, not a prescribed shaft speed.
//! All N-M free shaft momenta remain in the original canonical state.
//!
//! Three-shaft differentials and constant planetary relations fit this model;
//! tooth geometry, mesh compliance, backlash and capacity do not follow from it.

use super::{
    HeldTorque, InertiaTensor3, LieStepError, MAX_AXIAL_ROTORS, Mat3, RotorBank, RotorImpulse,
    RotorState, RotorStep, So3, Vec3, components, dd, dot_dd, increment, invalid, product, scalar,
    subtract3, vector, work_residual,
};
use fs_math::dd::Dd;

#[path = "transmission/passive.rs"]
mod passive;
pub use passive::{TransmissionNetworkStep, TransmissionNetworkStepper};

const ROUNDING: f64 = 8192.0 * f64::EPSILON;

/// One constant dimensionless relation sum_i coefficients[i] * q_i' = 0.
///
/// q_i is the original angle relative to the carrier. A positive row impulse
/// lambda adds coefficients[i]*lambda to axial momentum i. Thus scaling a row
/// rescales its multiplier inversely without changing physical motion.
/// For [input, left, right], [2, -1, -1] is an ideal equal-ratio differential:
/// input rate is the mean of two independently evolving output rates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShaftConstraint<const N: usize> {
    /// Signed coefficients in the original rotor inventory order.
    pub coefficients: [f64; N],
}

/// Prepared independent constant shaft relations, including differential modes.
///
/// Native fs-la Jacobi SVD admits the inertia-weighted constraint rank once.
/// Preparation may allocate; each subsequent operation uses fixed stack work,
/// O(N*M + M*M) arithmetic and one cached native 3x3 gyrostat. There is no step
/// factorization, iteration or allocation. Dependent or unresolved rows refuse;
/// the owner does not choose arbitrary reactions for redundant gear loops.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorTransmission<const N: usize, const M: usize> {
    bank: RotorBank<N>,
    constraints: [ShaftConstraint<N>; M],
    normalized: [[f64; N]; M],
    row_scale: [f64; M],
    weighted_scale: [f64; M],
    fixed_phase: [bool; N],
    reaction_free_components: [bool; 3],
    force_basis: [[f64; M]; N],
    rate_basis: [[f64; M]; N],
    multiplier_basis: [[f64; M]; M],
    recoil_axis: [Vec3; N],
    inertia: InertiaTensor3,
}

/// Explicit perfectly inelastic capture of all declared relations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionEngagement<const N: usize, const M: usize> {
    /// Original accepted transfer; attitude, total momentum and phases remain.
    pub impulse: RotorImpulse<N>,
    /// Actual reciprocal reaction multipliers in authored row units, N m s.
    pub constraint_impulse_nms: [f64; M],
    /// Nonnegative original carrier-plus-rotor kinetic loss, J.
    pub dissipated_energy_j: f64,
    /// Kinetic change plus physical engagement loss, J.
    pub work_energy_residual_j: f64,
}

/// Simultaneous original drives, ideal reactions and external body impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionImpulse<const N: usize, const M: usize> {
    /// Native actual shaft transfers, including drive and constraint reactions.
    pub impulse: RotorImpulse<N>,
    /// Actual reciprocal reaction multipliers in authored row units, N m s.
    pub constraint_impulse_nms: [f64; M],
    /// Work at each authored original shaft drive port, J.
    pub drive_work_j: [f64; N],
    /// Summed ideal reaction work; admitted rounding defect, never heat, J.
    pub reaction_work_j: f64,
    /// Kinetic change minus drives and external work, J.
    pub work_energy_residual_j: f64,
}

/// Free constrained precession with all remaining shaft momenta retained.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionDrift<const N: usize, const M: usize> {
    /// Complete native state after constrained attitude and phase drift.
    pub state_after: RotorState<N>,
    /// Actual absolute axial-momentum changes, N m s.
    pub axial_momentum_change_nms: [f64; N],
    /// Actual reciprocal reaction multipliers in authored row units, N m s.
    pub constraint_impulse_nms: [f64; M],
    /// Change in the original carrier-plus-rotor kinetic form, J.
    pub kinetic_energy_change_j: f64,
    /// Native reduced-gyrostat splitting defect, never transmission loss, J.
    pub drift_energy_defect_j: f64,
    /// Original kinetic change minus native drift defect, J.
    pub work_energy_residual_j: f64,
    /// Spatial angular-momentum change, N m s; no external impulse is applied.
    pub angular_momentum_residual_world: Vec3,
}

/// Driven constrained attitude, shaft spin and relative phase trajectory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionStep<const N: usize, const M: usize> {
    /// Native trajectory; shaft work includes ideal reaction transfers.
    pub motion: RotorStep<N>,
    /// Reactions from both driven endpoints and constrained precession, N m s.
    pub constraint_impulse_nms: [f64; M],
    /// Work at original drive ports, excluding ideal reactions, J.
    pub drive_work_j: [f64; N],
    /// Summed endpoint ideal reaction work, J.
    pub reaction_work_j: f64,
    /// Kinetic change minus drives, external work and native drift defect, J.
    pub work_energy_residual_j: f64,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare independent multi-shaft relations without altering any state.
    ///
    /// Each row removes exactly one relative-rate mode. M=0 retains all native
    /// free rotors; M=N is a complete carrier lock. Rows must be finite, nonzero
    /// and numerically independent after inertia weighting and row equilibration.
    /// A state must explicitly engage or already satisfy the declared relations.
    /// Existing phase offsets survive capture; no tooth-phase reset is inferred.
    /// An explicitly authored singleton row retains that shaft's phase exactly.
    pub fn transmission<const M: usize>(
        &self,
        constraints: [ShaftConstraint<N>; M],
    ) -> Result<RotorTransmission<N, M>, LieStepError> {
        if N > MAX_AXIAL_ROTORS || M > N {
            return Err(invalid("independent bounded shaft transmission rows"));
        }
        let mut normalized = [[0.0; N]; M];
        let mut row_scale = [0.0_f64; M];
        let mut weighted_scale = [0.0_f64; M];
        let mut fixed_phase = [false; N];
        let root_inertia: [f64; N] =
            core::array::from_fn(|i| self.rotors[i].axial_inertia_kg_m2.sqrt());
        let mut weighted = [[0.0; M]; N];
        for e in 0..M {
            let mut terms = 0;
            let mut last_nonzero = 0;
            for (i, &value) in constraints[e].coefficients.iter().enumerate() {
                scalar(value, "finite shaft relation coefficient")?;
                row_scale[e] = row_scale[e].max(value.abs());
                if value != 0.0 {
                    terms += 1;
                    last_nonzero = i;
                }
            }
            if terms == 1 {
                // Exact authored support, not a near-one projector diagonal:
                // this original coordinate has no physical phase freedom.
                fixed_phase[last_nonzero] = true;
            }
            if row_scale[e] == 0.0 {
                return Err(invalid("nonzero shaft relation"));
            }
            for i in 0..N {
                normalized[e][i] = constraints[e].coefficients[i] / row_scale[e];
                if constraints[e].coefficients[i] != 0.0 && normalized[e][i] == 0.0 {
                    return Err(invalid("resolved normalized shaft coefficients"));
                }
                weighted[i][e] = scalar(
                    normalized[e][i] / root_inertia[i],
                    "inertia-weighted shaft coefficient",
                )?;
                weighted_scale[e] = weighted_scale[e].max(weighted[i][e].abs());
            }
            if weighted_scale[e] == 0.0 {
                return Err(invalid("resolved inertia-weighted shaft row"));
            }
            for i in 0..N {
                let value = weighted[i][e] / weighted_scale[e];
                if weighted[i][e] != 0.0 && value == 0.0 {
                    return Err(invalid("resolved equilibrated shaft coefficients"));
                }
                weighted[i][e] = value;
            }
        }
        // A row may exert exactly no carrier torque along a component. Use
        // original represented coefficients and native error-free transforms,
        // with no rank/near-zero tolerance and no discarded residuals.
        let reaction_free_components = core::array::from_fn(|d| {
            constraints.iter().all(|row| {
                exact_zero_sum(
                    0.0,
                    (0..N).map(|i| (row.coefficients[i], components(self.rotors[i].axis_body)[d])),
                )
            })
        });
        let mut force_basis = [[0.0; M]; N];
        let mut rate_basis = [[0.0; M]; N];
        let mut multiplier_basis = [[0.0; M]; M];
        let mut reaction_axis = [Vec3::default(); M];
        if M != 0 {
            let matrix: Vec<f64> = weighted.iter().flatten().copied().collect();
            let factor = fs_la::factor::svd_jacobi(&matrix, N, M);
            let largest = scalar(factor.sigma[0], "transmission largest singular value")?;
            let smallest = scalar(factor.sigma[M - 1], "transmission smallest singular value")?;
            if largest <= 0.0 || smallest <= 512.0 * f64::EPSILON * N as f64 * largest {
                return Err(invalid("independent resolved shaft transmission rank"));
            }
            // Admit the returned numerical basis before using it as mechanics.
            for c in 0..M {
                for d in 0..M {
                    let mut inner = Dd::ZERO;
                    for i in 0..N {
                        inner = inner + dd(factor.u[i * M + c]) * dd(factor.u[i * M + d]);
                    }
                    let wanted = if c == d { 1.0 } else { 0.0 };
                    if scalar((inner - dd(wanted)).to_f64(), "shaft basis orthogonality")?.abs()
                        > ROUNDING
                    {
                        return Err(invalid("resolved orthonormal shaft reaction basis"));
                    }
                }
                for e in 0..M {
                    multiplier_basis[e][c] = scalar(
                        factor.v[e * M + c] / factor.sigma[c],
                        "shaft row multiplier basis",
                    )?;
                }
                let mut axis = [Dd::ZERO; 3];
                for i in 0..N {
                    force_basis[i][c] =
                        product(root_inertia[i], factor.u[i * M + c], "shaft reaction basis")?;
                    rate_basis[i][c] =
                        scalar(factor.u[i * M + c] / root_inertia[i], "shaft rate basis")?;
                    for (value, a) in axis.iter_mut().zip(components(self.rotors[i].axis_body)) {
                        *value = *value + dd(force_basis[i][c]) * dd(a);
                    }
                }
                for d in 0..3 {
                    if reaction_free_components[d] {
                        axis[d] = Dd::ZERO;
                    }
                }
                reaction_axis[c] = vector(
                    Vec3::new(axis[0].to_f64(), axis[1].to_f64(), axis[2].to_f64()),
                    "transmission carrier reaction axis",
                )?;
            }
            for i in 0..N {
                for e in 0..M {
                    let mut rebuilt = Dd::ZERO;
                    let mut scale = dd(weighted[i][e]).abs();
                    for c in 0..M {
                        let term =
                            dd(factor.u[i * M + c]) * dd(factor.sigma[c]) * dd(factor.v[e * M + c]);
                        rebuilt = rebuilt + term;
                        scale = scale + term.abs();
                    }
                    if scalar(
                        (rebuilt - dd(weighted[i][e])).to_f64(),
                        "shaft SVD residual",
                    )?
                    .abs()
                        > ROUNDING * scale.to_f64()
                    {
                        return Err(invalid("resolved original shaft relation basis"));
                    }
                }
            }
        }
        let mut recoil_axis = [Vec3::default(); N];
        for i in 0..N {
            let mut value = [Dd::ZERO; 3];
            for c in 0..M {
                for (entry, k) in value.iter_mut().zip(components(reaction_axis[c])) {
                    *entry = *entry + dd(force_basis[i][c]) * dd(k);
                }
            }
            recoil_axis[i] = if M == N {
                // Full rank has no free axial mode: preserve the native exact
                // lock limit instead of retaining SVD orthogonality roundoff.
                self.rotors[i]
                    .axis_body
                    .scale(self.rotors[i].axial_inertia_kg_m2)
            } else {
                vector(
                    Vec3::new(value[0].to_f64(), value[1].to_f64(), value[2].to_f64()),
                    "transmission original-shaft recoil",
                )?
            };
        }
        let inertia = if M == 0 {
            self.inertia
        } else if M == N {
            self.lock_mode([true; N])?.retained_inertia()
        } else {
            let mut matrix = self.inertia.matrix().m;
            for row in 0..3 {
                for col in row..3 {
                    let mut value = dd(matrix[3 * row + col]);
                    for axis in reaction_axis {
                        let k = components(axis);
                        value = value + dd(k[row]) * dd(k[col]);
                    }
                    matrix[3 * row + col] =
                        scalar(value.to_f64(), "retained transmission inertia")?;
                    matrix[3 * col + row] = matrix[3 * row + col];
                }
            }
            InertiaTensor3::new(Mat3 { m: matrix })?
        };
        Ok(RotorTransmission {
            bank: *self,
            constraints,
            normalized,
            row_scale,
            weighted_scale,
            fixed_phase,
            reaction_free_components,
            force_basis,
            rate_basis,
            multiplier_basis,
            recoil_axis,
            inertia,
        })
    }
}

impl<const N: usize, const M: usize> RotorTransmission<N, M> {
    /// Original row ordering, orientations and coefficient units.
    pub const fn constraints(&self) -> &[ShaftConstraint<N>; M] {
        &self.constraints
    }

    /// Original reduced tensor plus exactly the constrained axial modes.
    pub const fn retained_inertia(&self) -> InertiaTensor3 {
        self.inertia
    }

    /// Capture all declared relations in one perfectly inelastic internal event.
    ///
    /// Total momentum, attitude and original phase offsets are unchanged.
    /// The loss uses the original positive kinetic form of the actual impulse;
    /// no energy rescaling or prescribed shaft velocity is used.
    pub fn engage(
        &self,
        state: &RotorState<N>,
    ) -> Result<TransmissionEngagement<N, M>, LieStepError> {
        self.bank.admit_state(state)?;
        let impulse = if self.admit_transmission(state).is_ok() {
            self.bank.apply_impulse(state, [0.0; N], Vec3::default())?
        } else {
            let rates = self.bank.relative_rates(state)?;
            let mut modes = [Dd::ZERO; M];
            for c in 0..M {
                for i in 0..N {
                    modes[c] = modes[c] + dd(self.force_basis[i][c]) * dd(rates[i]);
                }
            }
            let mut transfer = [0.0; N];
            for i in 0..N {
                let mut value = Dd::ZERO;
                for c in 0..M {
                    value = value - dd(self.force_basis[i][c]) * modes[c];
                }
                transfer[i] = if M == N {
                    product(
                        -self.bank.rotors[i].axial_inertia_kg_m2,
                        rates[i],
                        "full-lock capture transfer",
                    )?
                } else {
                    scalar(value.to_f64(), "transmission capture transfer")?
                };
            }
            let delta_omega = self
                .inertia
                .velocity(self.bank.rotor_momentum(&transfer)?.scale(-1.0))?;
            for (i, value) in transfer.iter_mut().enumerate() {
                *value = scalar(
                    (dd(*value) + dot_dd(self.recoil_axis[i], delta_omega)).to_f64(),
                    "transmission capture carrier recoil",
                )?;
            }
            self.bank.apply_impulse(state, transfer, Vec3::default())?
        };
        self.admit_transmission(&impulse.state_after)?;
        let constraints =
            self.constraint_impulses(&impulse.axial_momentum_change_nms, &[0.0; N])?;
        let delta_carrier = self.bank.carrier_momentum(&impulse.state_after)?
            - self.bank.carrier_momentum(state)?;
        let mut loss = dot_dd(self.bank.inertia.velocity(delta_carrier)?, delta_carrier) * dd(0.5);
        for i in 0..N {
            let delta = impulse.axial_momentum_change_nms[i];
            loss = loss + dd(0.5 * delta / self.bank.rotors[i].axial_inertia_kg_m2) * dd(delta);
        }
        let loss = scalar(loss.to_f64(), "transmission engagement loss")?;
        let residual = scalar(
            (dd(impulse.kinetic_energy_change_j) + dd(loss)).to_f64(),
            "transmission capture work",
        )?;
        if loss < 0.0 {
            return Err(invalid("nonnegative transmission engagement loss"));
        }
        admit_work(residual.abs(), impulse.kinetic_energy_change_j.abs() + loss)?;
        Ok(TransmissionEngagement {
            impulse,
            constraint_impulse_nms: constraints,
            dissipated_energy_j: loss,
            work_energy_residual_j: residual,
        })
    }

    /// Apply all original shaft drives, reactions and external body impulse together.
    ///
    /// The input must already satisfy the relations. Reactions retain all free
    /// generalized drive momenta and have zero work at the common mean rates.
    /// Actual per-shaft transfers are reconstructed from authored row multipliers;
    /// rounding cannot silently become an unreported free-coordinate impulse.
    pub fn apply_impulse(
        &self,
        state: &RotorState<N>,
        drive_impulse_nms: [f64; N],
        external_impulse_body: Vec3,
    ) -> Result<TransmissionImpulse<N, M>, LieStepError> {
        self.admit_transmission(state)?;
        if M == 0 {
            let impulse =
                self.bank
                    .apply_impulse(state, drive_impulse_nms, external_impulse_body)?;
            return Ok(TransmissionImpulse {
                constraint_impulse_nms: [0.0; M],
                drive_work_j: impulse.shaft_work_j,
                reaction_work_j: 0.0,
                work_energy_residual_j: impulse.work_energy_residual_j,
                impulse,
            });
        }
        vector(external_impulse_body, "transmission external impulse")?;
        let delivered = super::super::super::torque_kick(
            state.total_momentum_body,
            external_impulse_body,
            1.0,
        )? - state.total_momentum_body;
        let (transfer, _) = self.impulse_response(&drive_impulse_nms, delivered)?;
        let known_zero_recoil = core::array::from_fn(|d| {
            self.reaction_free_components[d]
                && exact_zero_sum(
                    components(external_impulse_body)[d],
                    (0..N).map(|i| {
                        (
                            -components(self.bank.rotors[i].axis_body)[d],
                            drive_impulse_nms[i],
                        )
                    }),
                )
        });
        let impulse = self.bank.apply_impulse_with_zero_recoil(
            state,
            transfer,
            external_impulse_body,
            known_zero_recoil,
        )?;
        self.admit_transmission(&impulse.state_after)?;
        let constraints =
            self.constraint_impulses(&impulse.axial_momentum_change_nms, &drive_impulse_nms)?;
        let mut drive_work = [0.0; N];
        let mut reaction_work = Dd::ZERO;
        let mut scale = impulse.kinetic_energy_change_j.abs() + impulse.external_work_j.abs();
        for i in 0..N {
            let rate = if M == N || self.fixed_phase[i] {
                Dd::ZERO
            } else {
                dd(impulse.mean_relative_rate_rad_s[i])
            };
            drive_work[i] = scalar(
                (dd(drive_impulse_nms[i]) * rate).to_f64(),
                "transmission drive work",
            )?;
            reaction_work = reaction_work
                + (dd(impulse.axial_momentum_change_nms[i]) - dd(drive_impulse_nms[i])) * rate;
            scale += drive_work[i].abs() + impulse.shaft_work_j[i].abs();
        }
        let reaction_work = scalar(reaction_work.to_f64(), "ideal transmission reaction work")?;
        let residual = work_residual(
            impulse.kinetic_energy_change_j,
            &drive_work,
            impulse.external_work_j,
            0.0,
        )?;
        admit_work(reaction_work.abs() + residual.abs(), scale)?;
        Ok(TransmissionImpulse {
            impulse,
            constraint_impulse_nms: constraints,
            drive_work_j: drive_work,
            reaction_work_j: reaction_work,
            work_energy_residual_j: residual,
        })
    }
    /// Advance free precession with the retained independent shaft modes.
    ///
    /// This is the original second-order gyrostat flow on the exactly reduced
    /// kinetic form. Reactions during precession change the original absolute
    /// rotor momenta while preserving every unconstrained generalized momentum.
    /// Signed nonzero duration permits reversal. No projection loss is added.
    pub fn drift(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<TransmissionDrift<N, M>, LieStepError> {
        let mut drift = self.drift_momenta(state, duration_s)?;
        let first = self.bank.relative_rates(state)?;
        let last = self.bank.relative_rates(&drift.state_after)?;
        let mean = core::array::from_fn(|i| 0.5 * first[i] + 0.5 * last[i]);
        drift.state_after.relative_angle_rad = self.advance_phase(state, &mean, duration_s)?;
        Ok(drift)
    }

    /// Advance held drives and a body- or space-held external torque.
    ///
    /// Two simultaneous constrained half-impulses surround one retained-mode
    /// drift. Original phase increments satisfy all constant rows, with the same
    /// endpoint-rate quadrature as drive work. Spatial torques are transformed
    /// at both actual attitudes. The native drift defect remains numerical
    /// error, separate from drive work and the ideal reactions' zero work.
    /// This fixed-axis/fixed-ratio map is symmetric and second order; timestep
    /// refinement remains necessary. There is no backlash, finite capacity,
    /// tooth contact, parameter-change work or discrete adjoint in this owner.
    pub fn step(
        &self,
        state: &RotorState<N>,
        shaft_torque_nm: [f64; N],
        external_torque: HeldTorque,
        duration_s: f64,
    ) -> Result<TransmissionStep<N, M>, LieStepError> {
        if M == 0 {
            let motion = self
                .bank
                .step(state, shaft_torque_nm, external_torque, duration_s)?;
            return Ok(TransmissionStep {
                constraint_impulse_nms: [0.0; M],
                drive_work_j: motion.shaft_work_j,
                reaction_work_j: 0.0,
                work_energy_residual_j: motion.work_energy_residual_j,
                motion,
            });
        }
        scalar(duration_s, "transmission step duration")?;
        let half = 0.5 * duration_s;
        if duration_s == 0.0 || half == 0.0 {
            return Err(invalid("nonzero resolved transmission duration"));
        }
        let mut drive = [0.0; N];
        for i in 0..N {
            drive[i] = product(
                shaft_torque_nm[i],
                half,
                "transmission shaft torque impulse",
            )?;
        }
        let impulse_at = |rotation: So3| -> Result<Vec3, LieStepError> {
            let value = match external_torque {
                HeldTorque::Body(value) => vector(value, "transmission external body torque")?,
                HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
            };
            Ok(Vec3::new(
                product(value.x, half, "transmission external impulse")?,
                product(value.y, half, "transmission external impulse")?,
                product(value.z, half, "transmission external impulse")?,
            ))
        };
        let first = self.apply_impulse(state, drive, impulse_at(state.rotation)?)?;
        let drift = self.drift_momenta(&first.impulse.state_after, duration_s)?;
        let second = self.apply_impulse(
            &drift.state_after,
            drive,
            impulse_at(drift.state_after.rotation)?,
        )?;
        let mean = core::array::from_fn(|i| {
            0.5 * first.impulse.mean_relative_rate_rad_s[i]
                + 0.5 * second.impulse.mean_relative_rate_rad_s[i]
        });
        let phase = self.advance_phase(state, &mean, duration_s)?;
        let mut motion = self.bank.finish_step(
            state,
            &first.impulse,
            &drift.state_after,
            &second.impulse,
            duration_s,
            Some(phase),
        )?;
        motion.drift_energy_defect_j = drift.drift_energy_defect_j;
        motion.work_energy_residual_j = work_residual(
            motion.kinetic_energy_change_j,
            &motion.shaft_work_j,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        let mut drive_work = [0.0; N];
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work_j[i]) + dd(second.drive_work_j[i])).to_f64(),
                "transmission step drive work",
            )?;
        }
        let mut constraints = [0.0; M];
        for e in 0..M {
            constraints[e] = scalar(
                (dd(first.constraint_impulse_nms[e])
                    + dd(drift.constraint_impulse_nms[e])
                    + dd(second.constraint_impulse_nms[e]))
                .to_f64(),
                "integrated transmission reaction",
            )?;
        }
        let reaction_work = scalar(
            (dd(first.reaction_work_j) + dd(second.reaction_work_j)).to_f64(),
            "transmission endpoint reaction work",
        )?;
        let residual = work_residual(
            motion.kinetic_energy_change_j,
            &drive_work,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        Ok(TransmissionStep {
            motion,
            constraint_impulse_nms: constraints,
            drive_work_j: drive_work,
            reaction_work_j: reaction_work,
            work_energy_residual_j: residual,
        })
    }

    fn drift_momenta(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<TransmissionDrift<N, M>, LieStepError> {
        self.admit_transmission(state)?;
        let omega = self.bank.angular_velocity(state)?;
        // Construct the retained gyrostatic vector from the actual carrier rate:
        // a carrier exactly at rest must not acquire a rounding-sized spin from
        // subtracting two independently accumulated generalized rotor sums.
        let stored = if M == N {
            Vec3::default()
        } else if M == 0 {
            self.bank.rotor_momentum(&state.axial_momentum_nms)?
        } else {
            state.total_momentum_body - self.inertia.momentum(omega)?
        };
        let recovered = self.inertia.velocity(state.total_momentum_body - stored)?;
        admit_vector_error(recovered - omega, omega.norm().max(recovered.norm()))?;
        let mut after = *state;
        (after.rotation, after.total_momentum_body) = self.inertia.momentum_step(
            state.rotation,
            state.total_momentum_body,
            stored,
            duration_s,
        )?;
        let delta_total = after.total_momentum_body - state.total_momentum_body;
        let delta_omega = self.inertia.velocity(delta_total)?;
        let mut requested = [0.0; N];
        let mut change = [0.0; N];
        for i in 0..N {
            requested[i] = scalar(
                dot_dd(self.recoil_axis[i], delta_omega).to_f64(),
                "transmission precession reaction",
            )?;
            after.axial_momentum_nms[i] = increment(
                state.axial_momentum_nms[i],
                requested[i],
                "transmission transported axial momentum",
            )?;
            change[i] = after.axial_momentum_nms[i] - state.axial_momentum_nms[i];
        }
        self.admit_transmission(&after)?;
        self.bank.admit_recoil(
            self.bank.carrier_momentum(state)?,
            self.bank.carrier_momentum(&after)?,
            delta_total,
            &change,
            delta_total,
            &requested,
        )?;
        admit_vector_error(
            self.bank.angular_velocity(&after)? - (omega + delta_omega),
            omega.norm().max((omega + delta_omega).norm()),
        )?;
        let constraints = self.constraint_impulses(&change, &[0.0; N])?;
        let mean = omega + delta_omega.scale(0.5);
        let native_change = scalar(
            dot_dd(mean, delta_total).to_f64(),
            "transmission native drift energy",
        )?;
        let kinetic_change = self.bank.energy_change(state, &after)?;
        let residual = scalar(
            (dd(kinetic_change) - dd(native_change)).to_f64(),
            "transmission kinetic reconstruction",
        )?;
        admit_work(residual.abs(), self.transport_work_scale(state, &after)?)?;
        let momentum_residual = subtract3(
            after.rotation.rotate(after.total_momentum_body)?,
            state.rotation.rotate(state.total_momentum_body)?,
            Vec3::default(),
        );
        vector(momentum_residual, "transmission spatial momentum")?;
        Ok(TransmissionDrift {
            state_after: after,
            axial_momentum_change_nms: change,
            constraint_impulse_nms: constraints,
            kinetic_energy_change_j: kinetic_change,
            drift_energy_defect_j: native_change,
            work_energy_residual_j: residual,
            angular_momentum_residual_world: momentum_residual,
        })
    }

    // Free generalized source impulse before carrier recoil. An explicitly
    // fixed original port is exactly in the reaction space; remove it before
    // the numerical projection so static motor/load efforts do not invent an
    // unresolvable free transfer.
    fn project_drive(&self, drives: &[f64; N]) -> Result<[f64; N], LieStepError> {
        let mut modes = [Dd::ZERO; M];
        for i in 0..N {
            scalar(drives[i], "transmission shaft drive impulse")?;
            if !self.fixed_phase[i] {
                for c in 0..M {
                    modes[c] = modes[c] + dd(self.rate_basis[i][c]) * dd(drives[i]);
                }
            }
        }
        let mut projected = [0.0; N];
        for i in 0..N {
            if M == N || self.fixed_phase[i] {
                continue;
            }
            let mut value = dd(drives[i]);
            for c in 0..M {
                value = value - dd(self.force_basis[i][c]) * modes[c];
            }
            projected[i] = scalar(value.to_f64(), "transmission free drive transfer")?;
        }
        Ok(projected)
    }

    // Preserve structurally torque-free row components before rounding the
    // SVD projection: their carrier action is exactly the authored source.
    fn projected_rotor_momentum(
        &self,
        source: &[f64; N],
        projected: &[f64; N],
    ) -> Result<Vec3, LieStepError> {
        let mut value = components(self.bank.rotor_momentum(projected)?);
        for d in 0..3 {
            if self.reaction_free_components[d] {
                let mut component = Dd::ZERO;
                for i in 0..N {
                    component = component
                        + dd(components(self.bank.rotors[i].axis_body)[d]) * dd(source[i]);
                }
                value[d] = scalar(component.to_f64(), "authored torque-free carrier component")?;
            }
        }
        Ok(Vec3::new(value[0], value[1], value[2]))
    }

    // Affine response for preparation and implicit constitutive stages. State
    // publication, actual rounded transfers and all work/recoil gates still
    // belong to apply_impulse and the original rotor bank.
    fn impulse_response(
        &self,
        drives: &[f64; N],
        external: Vec3,
    ) -> Result<([f64; N], Vec3), LieStepError> {
        let mut transfer = self.project_drive(drives)?;
        let delta_omega = self
            .inertia
            .velocity(external - self.projected_rotor_momentum(drives, &transfer)?)?;
        for (i, value) in transfer.iter_mut().enumerate() {
            *value = scalar(
                (dd(*value) + dot_dd(self.recoil_axis[i], delta_omega)).to_f64(),
                "transmission original shaft impulse",
            )?;
        }
        Ok((transfer, delta_omega))
    }

    fn admit_transmission(&self, state: &RotorState<N>) -> Result<(), LieStepError> {
        let omega = self.bank.angular_velocity(state)?;
        for row in &self.normalized {
            let mut residual = Dd::ZERO;
            let mut scale = Dd::ZERO;
            for i in 0..N {
                let spin =
                    dd(state.axial_momentum_nms[i] / self.bank.rotors[i].axial_inertia_kg_m2);
                let carried = dot_dd(self.bank.rotors[i].axis_body, omega);
                residual = residual + dd(row[i]) * (spin - carried);
                scale = scale + dd(row[i]).abs() * (spin.abs() + carried.abs());
            }
            if scalar(residual.to_f64(), "transmission relative-rate residual")?.abs()
                > ROUNDING * scalar(scale.to_f64(), "transmission rate resolution scale")?
            {
                return Err(invalid(
                    "state does not satisfy shaft transmission relations",
                ));
            }
        }
        // Authored rows alone can conceal a finite forbidden rate when two
        // admitted rows are almost parallel. Admit the orthonormal physical
        // reaction modes too, without multiplying tolerances by row condition.
        for c in 0..M {
            let mut residual = Dd::ZERO;
            let mut scale = Dd::ZERO;
            for i in 0..N {
                let spin =
                    dd(state.axial_momentum_nms[i] / self.bank.rotors[i].axial_inertia_kg_m2);
                let carried = dot_dd(self.bank.rotors[i].axis_body, omega);
                residual = residual + dd(self.force_basis[i][c]) * (spin - carried);
                scale = scale + dd(self.force_basis[i][c]).abs() * (spin.abs() + carried.abs());
            }
            if scalar(residual.to_f64(), "transmission physical rate mode")?.abs()
                > ROUNDING * scalar(scale.to_f64(), "transmission physical rate scale")?
            {
                return Err(invalid(
                    "state does not satisfy physical shaft reaction modes",
                ));
            }
        }
        Ok(())
    }

    fn advance_phase(
        &self,
        state: &RotorState<N>,
        rates: &[f64; N],
        duration: f64,
    ) -> Result<[f64; N], LieStepError> {
        self.advance_phase_with_increment(state, rates, duration)
            .map(|(phase, _)| phase)
    }

    // Retain the actual requested increments separately from ordinary phase
    // addition roundoff for constitutive quadrature and spring-work admission.
    fn advance_phase_with_increment(
        &self,
        state: &RotorState<N>,
        rates: &[f64; N],
        duration: f64,
    ) -> Result<([f64; N], [f64; N]), LieStepError> {
        // An all-locked inventory has no phase degree of freedom. Otherwise
        // remove only numerical row residuals from the admitted mean-rate rule;
        // the orthogonal complement retains every independent differential mode.
        if M == N {
            return Ok((state.relative_angle_rad, [0.0; N]));
        }
        let mut modes = [Dd::ZERO; M];
        for c in 0..M {
            for i in 0..N {
                modes[c] = modes[c] + dd(self.force_basis[i][c]) * dd(rates[i]);
            }
        }
        let mut phase = state.relative_angle_rad;
        let mut actual = [0.0; N];
        let mut requested = [0.0; N];
        for i in 0..N {
            if self.fixed_phase[i] {
                // A singleton row fixes this original phase exactly. Dropping
                // projection roundoff here cannot remove any genuine free mode.
                continue;
            }
            let mut rate = dd(rates[i]);
            for c in 0..M {
                rate = rate - dd(self.rate_basis[i][c]) * modes[c];
            }
            let delta = product(
                duration,
                scalar(rate.to_f64(), "transmission phase rate")?,
                "transmission phase increment",
            )?;
            phase[i] = increment(phase[i], delta, "transmission relative phase")?;
            requested[i] = delta;
            actual[i] = (dd(phase[i]) - dd(state.relative_angle_rad[i])).to_f64();
        }
        for row in &self.normalized {
            let mut residual = Dd::ZERO;
            let mut scale = Dd::ZERO;
            for i in 0..N {
                let term = dd(row[i]) * dd(actual[i]);
                residual = residual + term;
                scale = scale + term.abs();
            }
            if scalar(residual.to_f64(), "transmission phase relation")?.abs()
                > ROUNDING * scalar(scale.to_f64(), "transmission phase scale")?
            {
                return Err(invalid("unresolved transmission phase increments"));
            }
        }
        Ok((phase, requested))
    }

    fn constraint_impulses(
        &self,
        actual: &[f64; N],
        drives: &[f64; N],
    ) -> Result<[f64; M], LieStepError> {
        let mut modes = [Dd::ZERO; M];
        for c in 0..M {
            for i in 0..N {
                modes[c] = modes[c] + dd(self.rate_basis[i][c]) * (dd(actual[i]) - dd(drives[i]));
            }
        }
        let mut result = [0.0; M];
        for e in 0..M {
            let mut scaled = Dd::ZERO;
            for c in 0..M {
                scaled = scaled + dd(self.multiplier_basis[e][c]) * modes[c];
            }
            result[e] = unscale_multiplier(scaled, self.row_scale[e], self.weighted_scale[e])?;
        }
        // Check actual local transfers in original authored units. Neither
        // stored spin nor total energy can hide a lost free-momentum component.
        for i in 0..N {
            let mut rebuilt = Dd::ZERO;
            let mut scale = dd(actual[i]).abs() + dd(drives[i]).abs();
            for e in 0..M {
                let term = dd(self.constraints[e].coefficients[i]) * dd(result[e]);
                rebuilt = rebuilt + term;
                scale = scale + term.abs();
            }
            let residual = scalar(
                (dd(actual[i]) - dd(drives[i]) - rebuilt).to_f64(),
                "transmission original reaction reconstruction",
            )?;
            if residual.abs()
                > ROUNDING * scalar(scale.to_f64(), "transmission impulse resolution scale")?
            {
                return Err(invalid(
                    "unresolved free shaft drive or transmission reaction",
                ));
            }
        }
        Ok(result)
    }
    fn transport_work_scale(
        &self,
        before: &RotorState<N>,
        after: &RotorState<N>,
    ) -> Result<f64, LieStepError> {
        let c0 = self.bank.carrier_momentum(before)?;
        let c1 = self.bank.carrier_momentum(after)?;
        let mean =
            self.bank.inertia.velocity(c0)?.scale(0.5) + self.bank.inertia.velocity(c1)?.scale(0.5);
        let mut scale = components(mean)
            .into_iter()
            .zip(components(c1 - c0))
            .map(|(rate, impulse)| (rate * impulse).abs())
            .sum::<f64>();
        for i in 0..N {
            let mean = (0.5 * before.axial_momentum_nms[i] + 0.5 * after.axial_momentum_nms[i])
                / self.bank.rotors[i].axial_inertia_kg_m2;
            scale += (mean * (after.axial_momentum_nms[i] - before.axial_momentum_nms[i])).abs();
        }
        scalar(scale, "transmission transport work scale")
    }
}

fn admit_work(residual: f64, scale: f64) -> Result<(), LieStepError> {
    scalar(residual, "transmission work residual")?;
    scalar(scale, "transmission work scale")?;
    if residual > ROUNDING * scale {
        return Err(invalid("unresolved ideal transmission work"));
    }
    Ok(())
}

fn admit_vector_error(error: Vec3, scale: f64) -> Result<(), LieStepError> {
    vector(error, "transmission carrier velocity reconstruction")?;
    scalar(scale, "transmission carrier velocity scale")?;
    if error.norm() > ROUNDING * scale {
        return Err(invalid("unresolved retained transmission carrier rate"));
    }
    Ok(())
}

// Oppositely sized factors have a representable product; combining them first
// avoids an avoidable overflowing/underflowing intermediate on row rescaling.
// Same-direction divisions are monotone in magnitude and cannot hide a finite
// answer behind an overflowing intermediate.
fn unscale_multiplier(value: Dd, first: f64, second: f64) -> Result<f64, LieStepError> {
    let result = if (first >= 1.0) != (second >= 1.0) {
        value / (dd(first) * dd(second))
    } else {
        value / dd(first) / dd(second)
    };
    let result = scalar(result.to_f64(), "transmission authored row impulse")?;
    if result == 0.0 && value != Dd::ZERO {
        return Err(invalid("resolved original shaft reaction units"));
    }
    Ok(result)
}

// Conservative exact-expression recognition, not a general sum certificate.
// A product or addition with any residual, overflow or uncertain underflow
// retains the ordinary recoil guard. No near-zero threshold selects zero.
fn exact_zero_sum(mut total: f64, terms: impl Iterator<Item = (f64, f64)>) -> bool {
    if !total.is_finite() {
        return false;
    }
    // Keep all possible 106-bit product tails above the subnormal floor.
    const MIN_PRODUCT: f64 = f64::from_bits(55_u64 << 52);
    for (a, b) in terms {
        if a == 0.0 || b == 0.0 {
            continue;
        }
        let (product, tail) = fs_math::eft::two_prod(a, b);
        if !product.is_finite() || product.abs() < MIN_PRODUCT || tail != 0.0 {
            return false;
        }
        let (sum, tail) = fs_math::eft::two_sum(total, product);
        if !sum.is_finite() || tail != 0.0 {
            return false;
        }
        total = sum;
    }
    total == 0.0
}
