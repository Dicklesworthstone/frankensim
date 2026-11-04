//! Shared admission, scalar branch/root and actual-state friction checks.
//! Both the one-port and simultaneous-port owners use these exact operations.
use super::*;
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

pub(super) fn admit_port(network: &CoupledModalSystem, port: &ModalFrictionPort,
    config: ModalFrictionConfig) -> Result<Vec<f64>, ModalFrictionError>
{
    // The numerical contact owner is shared with stateful TangentialLaw ports.
    // A law is admitted separately; attachment projection has no dry-law default.
    port.law.regularized_traction_1d(0.0, 0.0, 1.0)?;
    admit_attachments(network, &port.left, port.right.as_ref(), config)
}

pub(super) fn admit_attachments(network: &CoupledModalSystem, left: &ModalAttachment,
    right: Option<&ModalAttachment>, config: ModalFrictionConfig)
    -> Result<Vec<f64>, ModalFrictionError>
{
    if !(1..=128).contains(&config.max_iterations)
        || [config.maximum_force_n, config.maximum_slip_speed_m_s,
            config.force_absolute_tolerance_n, config.force_relative_tolerance,
            config.sticking_tolerance_m].iter().any(|x| !x.is_finite() || *x <= 0.0)
        || config.force_relative_tolerance >= 1.0 {
        return Err(invalid("friction requires explicit positive finite physical and numerical budgets").into());
    }
    let mut column = vec![0.0; network.mode_count()];
    for (attachment, sign) in std::iter::once((left, 1.0))
        .chain(right.map(|right| (right, -1.0))) {
        let model = network.models.get(attachment.component)
            .ok_or_else(|| invalid("friction attachment names an unknown component"))?;
        if attachment.shapes.len() != model.modes().len()
            || attachment.shapes.iter().any(|x| !x.is_finite()) {
            return Err(invalid("friction attachment must match its finite mass-normalized basis").into());
        }
        for (k, b) in attachment.shapes.iter().enumerate() {
            let i = network.offsets[attachment.component] + k;
            column[i] = finite(column[i] + sign*b)?;
        }
    }
    Ok(column)
}

pub(super) fn force_tolerance(config: ModalFrictionConfig, applied: f64, expected: f64)
    -> Result<f64, ModalCouplingError>
{
    finite(config.force_absolute_tolerance_n
        +config.force_relative_tolerance*applied.abs().max(expected.abs()))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve(law: &FrictionLaw, config: ModalFrictionConfig, compliance: f64,
    dt: f64, free_slip: f64, normal: f64, gate: Option<&CancelGate>)
    -> Result<(f64, ModalFrictionRegime, usize), ModalFrictionError>
{
    poll(gate)?;
    if normal == 0.0 { return Ok((0.0, ModalFrictionRegime::Released, 0)); }
    // Delegate only the held load-line solve. Actual modal displacement work
    // and endpoint admission remain with check_within below. The shared owner
    // isolates folded Gaussian branches before applying the force ceiling.
    // Its work is bounded by max_iterations; cancellation brackets that leaf.
    let candidate = law.advance(&law.initial_state(), &TangentialStepInput {
        normal_n: normal,
        slip_increment_m: [free_slip, 0.0],
        dt_s: dt,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: config.max_iterations,
            maximum_force_n: config.maximum_force_n,
            force_absolute_tolerance_n: config.force_absolute_tolerance_n,
            force_relative_tolerance: config.force_relative_tolerance,
        },
    });
    poll(gate)?;
    let step = match candidate {
        Err(TriboError::TangentialSolve { residual_n, tolerance_n, iterations }) => {
            return Err(ModalCouplingError::ContactSolve { residual_n, tolerance_n, iterations }.into());
        }
        Err(error @ TriboError::OutsideApplicability { field: "maximum_force_n", value, .. }) => {
            limit("friction traction", value, config.maximum_force_n)?;
            return Err(error.into());
        }
        result => result?,
    };
    let regime = match step.regime {
        fs_tribo::FrictionRegime::Sticking => ModalFrictionRegime::Sticking,
        fs_tribo::FrictionRegime::Sliding => ModalFrictionRegime::Sliding,
    };
    Ok((step.traction_n[0], regime, step.iterations))
}

/// The actual endpoint slip of an accepted candidate, as a distance:
/// (u - b^T v_end) dt, with v_end the candidate's end velocities.
pub(super) fn endpoint_slip(column: &[f64], models: &[ModalAcousticTimeModel], drive: FrictionDrive, dt: f64)
    -> Result<f64, ModalCouplingError>
{
    let mut velocity = 0.0;
    let mut index = 0;
    for model in models {
        for state in model.states() {
            velocity = finite(velocity + column[index] * state.velocity_m_sqrt_kg_per_s)?;
            index += 1;
        }
    }
    finite(finite(drive.speed_m_s - velocity)? * dt)
}

/// The accepted step's work and branch check (MagentaBadger 8060). The branch
/// is checked on `line_slip`, the slip of the line it was solved on: the end
/// slip for the endpoint line. The work is F d on the actual step displacement
/// `distance`, returned with the kinetic residual as (loss, declared, residual):
/// a sliding step with d along F returns it as loss; a sticking step, and a
/// sliding step whose d runs against F, return it signed as declared work.
#[allow(clippy::too_many_arguments)]
pub(super) fn check(law: &FrictionLaw, config: ModalFrictionConfig, dt: f64,
    drive: FrictionDrive, traction: f64, line_slip: f64, distance: f64, regime: ModalFrictionRegime,
    iterations: usize) -> Result<(f64, f64, f64), ModalFrictionError>
{
    check_within(law, config, dt, drive, traction, line_slip, distance, regime, iterations, 1.0)
}

/// [`check`] with the sliding force residual bounded by `share` of the declared
/// tolerance: all of it (1.0) at the actual endpoint, half of it on a prediction
/// solved with headroom. The sticking-slip bound is not scaled.
#[allow(clippy::too_many_arguments)]
pub(super) fn check_within(law: &FrictionLaw, config: ModalFrictionConfig, dt: f64,
    drive: FrictionDrive, traction: f64, line_slip: f64, distance: f64, regime: ModalFrictionRegime,
    iterations: usize, share: f64) -> Result<(f64, f64, f64), ModalFrictionError>
{
    let speed = finite(line_slip/dt)?;
    let mut constitutive_residual_n = 0.0;
    let work = finite(traction*distance)?;
    let (loss, stuck) = match regime {
        ModalFrictionRegime::Shearing { .. } | ModalFrictionRegime::Microslip => return Err(invalid("material shear must be checked by its stateful constitutive owner").into()),
        ModalFrictionRegime::Released => (0.0, 0.0),
        ModalFrictionRegime::Sticking => {
            if line_slip.abs() > config.sticking_tolerance_m {
                return Err(ModalFrictionError::Sticking { residual_m: line_slip,
                    tolerance_m: config.sticking_tolerance_m });
            }
            (0.0, work)
        }
        ModalFrictionRegime::Sliding => {
            limit("friction slip speed", speed.abs(), config.maximum_slip_speed_m_s)?;
            if speed == 0.0 || (traction != 0.0 && traction.signum() != speed.signum()) {
                return Err(invalid("actual sliding force must oppose nonzero relative slip").into());
            }
            let expected = finite(law.kinetic_coefficient(speed.abs())?
                *drive.normal_force_n*speed.signum())?;
            constitutive_residual_n = finite(traction-expected)?;
            let tolerance = share*force_tolerance(config, traction, expected)?;
            if constitutive_residual_n.abs() > tolerance {
                return Err(ModalCouplingError::ContactSolve { residual_n: constitutive_residual_n,
                    tolerance_n: tolerance, iterations }.into());
            }
            // Applied traction, not a nearby law value: this closes the
            // actual discrete work exchanged with the original stepper.
            if work < 0.0 { (0.0, work) } else { (work, 0.0) }
        }
    };
    Ok((loss, stuck, constitutive_residual_n))
}

#[cfg(test)]
mod trait_tests {
    use super::*;

    /// G1/G5: compare against the actual pre-existing modal branch/root owner,
    /// including its endpoint convention and bisection stopping iteration.
    #[test]
    fn g1_trait_stick_slip_matches_modal_friction_owner_bitwise() {
        let config = ModalFrictionConfig {
            max_iterations: 80,
            maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-11,
            force_relative_tolerance: 1e-10,
            sticking_tolerance_m: 1e-12,
        };
        let dt = 1.0 / 48000.0;
        for law in [
            FrictionLaw::Coulomb {
                static_mu: 0.8,
                kinetic_mu: 0.4,
            },
            FrictionLaw::Stribeck {
                static_mu: 0.8,
                kinetic_mu: 0.4,
                characteristic_speed: 0.04,
                viscous_per_speed: 0.0,
            },
            FrictionLaw::Hyperbolic {
                static_mu: 0.8,
                kinetic_mu: 0.4,
                characteristic_speed: 0.04,
            },
            FrictionLaw::VelocityDependent {
                static_mu: 0.8,
                mu_zero: 0.4,
                slope_per_speed: 0.1,
            },
            FrictionLaw::VelocityDependent {
                static_mu: 0.8,
                mu_zero: 0.4,
                slope_per_speed: -0.1,
            },
            FrictionLaw::VelocityDependent {
                static_mu: 0.8,
                mu_zero: 0.4,
                slope_per_speed: -2.0,
            },
        ] {
            for normal in [0.0, 0.1, 1.0] {
                for free_slip in [-0.3 * dt, -0.001 * dt, 0.0, 0.001 * dt, 0.3 * dt] {
                    for compliance in [1e-8, 1e-6, 1e-4] {
                        let (traction, regime, iterations) =
                            solve(&law, config, compliance, dt, free_slip, normal, None).unwrap();
                        let step = law
                            .advance(
                                &law.initial_state(),
                                &TangentialStepInput {
                                    normal_n: normal,
                                    slip_increment_m: [free_slip, 0.0],
                                    dt_s: dt,
                                    contact_temperature_k: None,
                                    nominal_pressure_pa: None,
                                    compliance_m_per_n: compliance,
                                    solve: TangentialSolveConfig {
                                        max_iterations: config.max_iterations,
                                        maximum_force_n: config.maximum_force_n,
                                        force_absolute_tolerance_n: config
                                            .force_absolute_tolerance_n,
                                        force_relative_tolerance: config.force_relative_tolerance,
                                    },
                                },
                            )
                            .unwrap();
                        assert_eq!(step.traction_n[0].to_bits(), traction.to_bits());
                        assert_eq!(step.iterations, iterations);
                        if normal > 0.0 {
                            assert_eq!(
                                step.regime == fs_tribo::FrictionRegime::Sticking,
                                regime == ModalFrictionRegime::Sticking
                            );
                        }
                    }
                }
            }
        }
    }

    /// A folded Gaussian line has three roots here. The old scalar copy
    /// rejected every feasible branch when the force ceiling lay between its
    /// second and third force roots. The constitutive owner isolates first.
    #[test]
    fn g1_modal_gaussian_branch_survives_a_negative_cap_residual() {
        let law = FrictionLaw::Stribeck {
            static_mu: 1.0,
            kinetic_mu: 0.1,
            characteristic_speed: 1.0,
            viscous_per_speed: 0.0,
        };
        let config = ModalFrictionConfig {
            max_iterations: 80,
            maximum_force_n: 0.95,
            maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-14,
            force_relative_tolerance: 1e-14,
            sticking_tolerance_m: 1e-12,
        };
        for sign in [-1.0, 1.0] {
            let (force, regime, used) = solve(&law, config, 3.0, 1.0, sign * 3.05, 1.0, None).unwrap();
            let slip = sign * 3.05 - 3.0 * force;
            // det-ok: independent direct constitutive expression.
            let expected = sign * (0.1 + 0.9 * (-slip * slip).exp());
            assert_eq!(regime, ModalFrictionRegime::Sliding);
            assert!(force.abs() < 0.11);
            assert!(slip.abs() > 2.7);
            assert!((force - expected).abs() < 3e-14);
            assert!(used <= config.max_iterations);
        }
    }
}
