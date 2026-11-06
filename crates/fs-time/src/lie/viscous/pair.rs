//! Reciprocal viscous torque between two finite-inertia rotating rigid bodies.
//!
//! D is fixed in body 0's material axes. With Q = R0^-1 R1, the slip is
//! r = omega0 - Q omega1 and the torques are -D r and Q^T D r. Both bodies
//! recoil; the couple conserves total spatial angular momentum and dissipates
//! r^T D r. External drives, internal transfer work, physical heat and native
//! free-flow energy defects have separate accounts.
//!
//! Two midpoint torque half-stages surround the original compensated free
//! motion of BOTH bodies. At each actual endpoint configuration the relative
//! mobility K = I0^-1 + Q I1^-1 Q^T is factored by native Cholesky. A second
//! 3-by-3 SPD solve determines the common impulse, including semidefinite D.
//! There is no iteration, numerical Jacobian, per-step eigensolve or D inverse.

use super::{
    Dd, HeldTorque, InertiaTensor3, LieStepError, Mat3, So3, Vec3, Vector, add, admit_damping,
    apply, cholesky, compensated, components, dd, dot, drive_at, energy, held, invalid, rounded,
    scalar, subtract, to_space, velocity,
};

/// Two fixed positive inertia tensors, a material damping tensor and interval.
///
/// Preparation caches inertia-normalized principal mobility columns. Each
/// step uses two native 3-by-3 Cholesky factorizations and one solve per torque
/// half-stage. Those four small factors use the native owner's bounded heap
/// storage; no dense six-body-coordinate system or generic nonlinear solver
/// is introduced.
///
/// The coupler's orientation is attached to body 0. This is a smooth viscous
/// rotational element, without a joint constraint, translational forces,
/// shaft elasticity, dry friction, thermal allocation or internal rotors.
#[derive(Debug, Clone)]
pub struct ViscousPairStepper {
    inertia: [InertiaTensor3; 2],
    damping_first_body: Mat3,
    duration_s: f64,
    half_s: f64,
    reference_moment: f64,
    mobility_columns: [[Vector; 3]; 2],
}

/// One body's accepted motion and its separate source/transfer accounts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousPairBodyStep {
    pub rotation_after: So3,
    /// Authoritative material-frame COM angular momentum [kg m^2/s].
    pub momentum_after_body: Vec3,
    pub angular_velocity_after_body: Vec3,
    /// Sum of the actual endpoint external-drive impulses in space [N m s].
    pub drive_angular_impulse_space: Vec3,
    /// Equal and opposite to the other body's coupling impulse [N m s].
    pub coupling_angular_impulse_space: Vec3,
    /// Spatial momentum change minus drive and coupling impulses [N m s].
    pub angular_momentum_residual_space: Vec3,
    /// Signed work of this body's external drive [J].
    pub drive_work_j: f64,
    /// Signed mechanical work transferred to this body by the coupler [J].
    /// The two transfer works sum to minus the pair's relative-slip heat.
    pub coupling_work_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Native torque-free drift energy change alone, never counted as heat.
    pub drift_energy_defect_j: f64,
    /// Kinetic change - drive work - coupling work - drift defect [J].
    /// Includes final public-momentum rounding; no energy projection.
    pub work_energy_residual_j: f64,
}

/// Accepted finite pair, with shared physical heat and total balances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousPairStep {
    pub bodies: [ViscousPairBodyStep; 2],
    /// Sum(s * r_mid^T D r_mid) over both torque half-stages [J].
    /// Nonnegative mechanical dissipation, with no assigned heat destination.
    pub dissipated_energy_j: f64,
    /// Total spatial momentum change minus the two external drives [N m s].
    /// Coupling impulses cancel exactly in the public receipt.
    pub angular_momentum_residual_space: Vec3,
    /// Total kinetic change - external work + heat - drift defects [J].
    pub work_energy_residual_j: f64,
}

impl ViscousPairStepper {
    /// Prepare a viscous element between two independently rotating bodies.
    ///
    /// D [N m s/rad] must be finite, symmetric and positive semidefinite.
    /// The existing compensated principal-minor admission is reused without
    /// clipping negative eigenvalues. The duration must be finite and positive.
    ///
    /// The two actual endpoint configurations each admit s * tr(K D) <= 2,
    /// where s = duration / 2. This numerical upper relaxation-rate estimate
    /// prevents unresolved midpoint sign reversal; it is not an error bound
    /// or a replacement for rotation/decay timestep refinement.
    pub fn new(
        inertia: [InertiaTensor3; 2],
        damping_first_body: Mat3,
        duration_s: f64,
    ) -> Result<Self, LieStepError> {
        scalar(duration_s, "viscous pair duration")?;
        let half_s = 0.5 * duration_s;
        if duration_s <= 0.0 || half_s == 0.0 {
            return Err(invalid("positive resolved viscous pair duration"));
        }
        let damping_first_body = admit_damping(damping_first_body)?;
        let reference_moment = inertia
            .iter()
            .flat_map(|i| components(i.principal_moments()))
            .fold(f64::INFINITY, f64::min);
        let mut mobility_columns = [[[Dd::ZERO; 3]; 3]; 2];
        for body in 0..2 {
            for (axis, moment) in components(inertia[body].principal_moments())
                .into_iter()
                .enumerate()
            {
                // Taking square roots before division avoids underflow in
                // reference_moment / moment across widely separated scales.
                let weight = reference_moment.sqrt() / moment.sqrt();
                let mut column = [Dd::ZERO; 3];
                column[axis] = dd(weight);
                mobility_columns[body][axis] = to_space(inertia[body].principal_to_body(), column)?;
            }
        }
        Ok(Self {
            inertia,
            damping_first_body,
            duration_s,
            half_s,
            reference_moment,
            mobility_columns,
        })
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// The damping tensor carried by body 0's material frame [N m s/rad].
    pub const fn damping_first_body(&self) -> Mat3 {
        self.damping_first_body
    }

    /// Advance both finite bodies under reciprocal drag and held drives.
    ///
    /// Each drive may independently be held in its body's material frame or
    /// in space. Spatial drives are transformed at BOTH actual endpoint
    /// attitudes. Both momenta retain DD impulse tails through native free
    /// drift. Source work and relative-slip heat use the accepted midpoint
    /// rates, with local constitutive/work checks before publication.
    ///
    /// The symmetric composition is second order for this smooth law. Every
    /// input is supplied by value; invalid data, unresolved mobility/decay or
    /// numerical refusal leaves both states unchanged. Numerical receipts
    /// include group arithmetic and final public-state rounding. No finite
    /// thermal reservoir, feedback law, constraints or adjoint is inferred.
    pub fn step(
        &self,
        rotations: [So3; 2],
        momenta_body: [Vec3; 2],
        drives: [HeldTorque; 2],
    ) -> Result<ViscousPairStep, LieStepError> {
        let initial = [
            held::finite_vector(momenta_body[0], "viscous pair momentum")?,
            held::finite_vector(momenta_body[1], "viscous pair momentum")?,
        ];
        let first = self.stage(rotations, initial, drives)?;
        let mut rotations_after = rotations;
        let mut drifted = first.momentum;
        for body in 0..2 {
            (rotations_after[body], drifted[body]) = self.inertia[body].compensated_momentum_step(
                rotations[body],
                first.momentum[body],
                self.duration_s,
            )?;
        }
        let second = self.stage(rotations_after, drifted, drives)?;
        let pair_impulse = add(first.coupling_impulse_space, second.coupling_impulse_space);
        let public_pair_impulse = rounded(pair_impulse, "viscous pair coupling impulse")?;
        let finish = |body: usize| -> Result<_, LieStepError> {
            let momentum_after_body =
                rounded(second.momentum[body], "viscous pair output momentum")?;
            let accepted =
                held::finite_vector(momentum_after_body, "viscous pair accepted momentum")?;
            let angular_velocity_after_body = rounded(
                velocity(self.inertia[body], accepted)?,
                "viscous pair output angular velocity",
            )?;
            let kinetic =
                energy(self.inertia[body], accepted)? - energy(self.inertia[body], initial[body])?;
            let drift = energy(self.inertia[body], drifted[body])?
                - energy(self.inertia[body], first.momentum[body])?;
            let work = first.drive_work[body] + second.drive_work[body];
            let transfer = first.coupling_work[body] + second.coupling_work[body];
            let drive_impulse = add(
                to_space(rotations[body], first.drive_impulse[body])?,
                to_space(rotations_after[body], second.drive_impulse[body])?,
            );
            let coupling = if body == 0 {
                pair_impulse
            } else {
                pair_impulse.map(|v| -v)
            };
            let balance = subtract(
                subtract(
                    to_space(rotations_after[body], accepted)?,
                    to_space(rotations[body], initial[body])?,
                ),
                add(drive_impulse, coupling),
            );
            let receipt = ViscousPairBodyStep {
                rotation_after: rotations_after[body],
                momentum_after_body,
                angular_velocity_after_body,
                drive_angular_impulse_space: rounded(drive_impulse, "viscous pair drive impulse")?,
                coupling_angular_impulse_space: if body == 0 {
                    public_pair_impulse
                } else {
                    public_pair_impulse.scale(-1.0)
                },
                angular_momentum_residual_space: rounded(balance, "viscous pair angular balance")?,
                drive_work_j: scalar(work.to_f64(), "viscous pair drive work")?,
                coupling_work_j: scalar(transfer.to_f64(), "viscous pair transfer work")?,
                kinetic_energy_change_j: scalar(kinetic.to_f64(), "viscous pair kinetic change")?,
                drift_energy_defect_j: scalar(drift.to_f64(), "viscous pair drift defect")?,
                work_energy_residual_j: scalar(
                    (kinetic - work - transfer - drift).to_f64(),
                    "viscous pair body work balance",
                )?,
            };
            Ok((receipt, kinetic - work - drift, balance))
        };
        let (a, account_a, balance_a) = finish(0)?;
        let (b, account_b, balance_b) = finish(1)?;
        let heat = first.heat + second.heat;
        Ok(ViscousPairStep {
            bodies: [a, b],
            dissipated_energy_j: scalar(heat.to_f64(), "viscous pair heat")?,
            angular_momentum_residual_space: rounded(
                add(balance_a, balance_b),
                "viscous pair total angular balance",
            )?,
            work_energy_residual_j: scalar(
                (account_a + account_b + heat).to_f64(),
                "viscous pair total work balance",
            )?,
        })
    }

    fn stage(
        &self,
        rotations: [So3; 2],
        momentum: [Vector; 2],
        drives: [HeldTorque; 2],
    ) -> Result<PairStage, LieStepError> {
        let before = [
            velocity(self.inertia[0], momentum[0])?,
            velocity(self.inertia[1], momentum[1])?,
        ];
        let mut drive_impulse = [[Dd::ZERO; 3]; 2];
        for body in 0..2 {
            let drive = drive_at(rotations[body], drives[body])?;
            for axis in 0..3 {
                drive_impulse[body][axis] = compensated::checked(
                    dd(self.half_s) * drive[axis],
                    "viscous pair drive impulse",
                )?;
            }
        }
        let trial = if self.damping_first_body.max_abs() == 0.0 {
            [Dd::ZERO; 3]
        } else {
            let drive_change = [
                velocity(self.inertia[0], drive_impulse[0])?,
                velocity(self.inertia[1], drive_impulse[1])?,
            ];
            let g = add(
                slip(rotations, before)?,
                slip(rotations, drive_change)?.map(|v| dd(0.5) * v),
            );
            self.coupling_impulse(rotations, g)?
        };
        let trial_space = to_space(rotations[0], trial)?;
        let trial_body = [
            trial,
            compensated::spatial_vector_to_body(rotations[1], trial_space)?.map(|v| -v),
        ];
        let mut after = momentum;
        let mut actual = [[Dd::ZERO; 3]; 2];
        let mut mean = before;
        for body in 0..2 {
            for axis in 0..3 {
                after[body][axis] = compensated::checked(
                    momentum[body][axis] + drive_impulse[body][axis] + trial_body[body][axis],
                    "viscous pair kicked momentum",
                )?;
            }
            actual[body] = subtract(after[body], momentum[body]);
            let final_velocity = velocity(self.inertia[body], after[body])?;
            mean[body] =
                core::array::from_fn(|i| dd(0.5) * before[body][i] + dd(0.5) * final_velocity[i]);
        }
        let relative = slip(rotations, mean)?;
        let drag = apply(self.damping_first_body, relative)?;
        let impulse = drag.map(|v| dd(-self.half_s) * v);
        let coupling_impulse_space = to_space(rotations[0], impulse)?;
        let coupling = [
            impulse,
            compensated::spatial_vector_to_body(rotations[1], coupling_impulse_space)?.map(|v| -v),
        ];
        // Rotation and midpoint products can cancel a transverse component
        // while their roundoff still follows the common physical impulse.
        // Retain componentwise drive/kick scaling, and allow that local
        // vector-sized transport error. No background momentum, kinetic
        // energy or absolute tolerance floor enters this constitutive gate.
        let coupling_scale = trial
            .iter()
            .chain(impulse.iter())
            .map(|v| v.to_f64().abs())
            .fold(0.0_f64, f64::max);
        let mut drive_work = [Dd::ZERO; 2];
        let mut coupling_work = [Dd::ZERO; 2];
        let mut work_scale = 0.0_f64;
        for body in 0..2 {
            let mut body_work_scale = 0.0_f64;
            for axis in 0..3 {
                let residual = scalar(
                    (actual[body][axis] - drive_impulse[body][axis] - coupling[body][axis])
                        .to_f64(),
                    "viscous pair constitutive residual",
                )?;
                let impulses = [
                    actual[body][axis],
                    drive_impulse[body][axis],
                    coupling[body][axis],
                ];
                let scale = impulses
                    .iter()
                    .map(|v| v.to_f64().abs())
                    .fold(coupling_scale, f64::max);
                if residual.abs() > 8192.0 * f64::EPSILON * scale {
                    return Err(invalid("unresolved reciprocal viscous torque"));
                }
                for component in impulses {
                    body_work_scale = body_work_scale.max(
                        scalar(
                            (component * mean[body][axis]).to_f64(),
                            "viscous pair component work",
                        )?
                        .abs(),
                    );
                }
            }
            work_scale = work_scale.max(body_work_scale);
            drive_work[body] = dot(drive_impulse[body], mean[body]);
            coupling_work[body] = dot(coupling[body], mean[body]);
            let residual = scalar(
                (dot(actual[body], mean[body]) - drive_work[body] - coupling_work[body]).to_f64(),
                "viscous pair stage body work",
            )?;
            if residual.abs() > 24576.0 * f64::EPSILON * body_work_scale {
                return Err(invalid("unresolved reciprocal viscous body work"));
            }
        }
        let heat = dd(self.half_s) * dot(relative, drag);
        if scalar(heat.to_f64(), "viscous pair dissipative work")? < 0.0 {
            return Err(invalid("nonnegative viscous pair dissipative work"));
        }
        for axis in 0..3 {
            work_scale = work_scale.max(
                scalar(
                    (impulse[axis] * relative[axis]).to_f64(),
                    "viscous pair slip work",
                )?
                .abs(),
            );
        }
        let residual = scalar(
            (coupling_work[0] + coupling_work[1] + heat).to_f64(),
            "viscous pair reciprocal work",
        )?;
        if residual.abs() > 24576.0 * f64::EPSILON * work_scale {
            return Err(invalid("unresolved reciprocal viscous work"));
        }
        Ok(PairStage {
            momentum: after,
            drive_impulse,
            coupling_impulse_space,
            drive_work,
            coupling_work,
            heat,
        })
    }

    fn coupling_impulse(&self, rotations: [So3; 2], g: Vector) -> Result<Vector, LieStepError> {
        // Columns satisfy sum(c c^T) = I_ref K. Normalizing inertia BEFORE
        // building mobility avoids an overflowing inverse for tiny inertias.
        let mut columns = [[Dd::ZERO; 3]; 6];
        columns[..3].copy_from_slice(&self.mobility_columns[0]);
        for axis in 0..3 {
            columns[3 + axis] = compensated::spatial_vector_to_body(
                rotations[0],
                to_space(rotations[1], self.mobility_columns[1][axis])?,
            )?;
        }
        let mut mobility = [0.0; 9];
        for i in 0..3 {
            for j in 0..=i {
                let mut entry = Dd::ZERO;
                for column in columns {
                    entry = entry + column[i] * column[j];
                }
                let entry = scalar(entry.to_f64(), "viscous pair relative mobility")?;
                mobility[3 * i + j] = entry;
                mobility[3 * j + i] = entry;
            }
        }
        let factor = cholesky(&mobility, 3)
            .map_err(|_| invalid("resolved positive viscous pair mobility"))?;
        let lower: [Vector; 3] =
            core::array::from_fn(|column| core::array::from_fn(|row| dd(factor.l(row, column))));
        let damped = [
            apply(self.damping_first_body, lower[0])?,
            apply(self.damping_first_body, lower[1])?,
            apply(self.damping_first_body, lower[2])?,
        ];
        let mut matrix = [0.0; 9];
        let mut relaxation = Dd::ZERO;
        for i in 0..3 {
            for j in 0..=i {
                let transfer = compensated::divide(
                    dd(self.half_s) * dot(lower[i], damped[j]),
                    dd(self.reference_moment),
                )?;
                if i == j {
                    relaxation = relaxation + transfer;
                }
                let entry = scalar(
                    (dd(if i == j { 1.0 } else { 0.0 }) + dd(0.5) * transfer).to_f64(),
                    "viscous pair midpoint factor",
                )?;
                matrix[3 * i + j] = entry;
                matrix[3 * j + i] = entry;
            }
        }
        let resolved = scalar(
            relaxation.to_f64() * (1.0 + 512.0 * f64::EPSILON),
            "viscous pair relaxation number",
        )?;
        if resolved < 0.0 || resolved > 2.0 {
            return Err(invalid("viscous pair interval does not resolve relaxation"));
        }
        let midpoint_factor =
            cholesky(&matrix, 3).map_err(|_| invalid("viscous pair midpoint factor not SPD"))?;
        let drag = apply(self.damping_first_body, g)?;
        // z = L^T j: (Id + s/(2 I_ref) L^T D L) z = -s L^T D g.
        // This needs neither D^-1 nor a hand-written triangular solver.
        let mut z = [0.0; 3];
        for i in 0..3 {
            z[i] = scalar(
                (dd(-self.half_s) * dot(lower[i], drag)).to_f64(),
                "viscous pair relative impulse right hand side",
            )?;
        }
        midpoint_factor.solve(&mut z);
        let mut relative_midpoint = g;
        for i in 0..3 {
            let change = lower[0][i] * dd(z[0]) + lower[1][i] * dd(z[1]) + lower[2][i] * dd(z[2]);
            relative_midpoint[i] = compensated::checked(
                g[i] + dd(0.5) * compensated::divide(change, dd(self.reference_moment))?,
                "viscous pair relative midpoint",
            )?;
        }
        Ok(apply(self.damping_first_body, relative_midpoint)?.map(|v| dd(-self.half_s) * v))
    }
}

struct PairStage {
    momentum: [Vector; 2],
    drive_impulse: [Vector; 2],
    coupling_impulse_space: Vector,
    drive_work: [Dd; 2],
    coupling_work: [Dd; 2],
    heat: Dd,
}

fn slip(rotations: [So3; 2], rates: [Vector; 2]) -> Result<Vector, LieStepError> {
    Ok(subtract(
        rates[0],
        compensated::spatial_vector_to_body(rotations[0], to_space(rotations[1], rates[1])?)?,
    ))
}
