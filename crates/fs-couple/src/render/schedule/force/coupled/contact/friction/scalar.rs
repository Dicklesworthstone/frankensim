//! Shared admission, scalar branch/root and actual-state friction checks.
//! Both the one-port and simultaneous-port owners use these exact operations.
use super::*;

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
    if normal == 0.0 { return Ok((0.0, ModalFrictionRegime::Released, 0)); }
    let static_mu = law.static_coefficient()?;
    let cap = finite(static_mu*normal)?;
    let required = finite(free_slip/compliance)?;
    if required.abs() <= cap {
        limit("friction traction", required.abs(), config.maximum_force_n)?;
        return Ok((required, ModalFrictionRegime::Sticking, 0));
    }
    let sign = free_slip.signum();
    let evaluate = |magnitude: f64| -> Result<(f64, f64), ModalFrictionError> {
        poll(gate)?;
        let slip = finite(free_slip.abs()-compliance*magnitude)?.max(0.0);
        // Use the one-sided kinetic limit at the sticking endpoint; never
        // feed zero into the owner's strictly-positive kinetic query.
        let mu = if slip == 0.0 {
            law.zero_speed_kinetic_coefficient()?
        } else { law.kinetic_coefficient(finite(slip/dt)?)? };
        let expected = finite(mu*normal)?;
        Ok((finite(magnitude-expected)?, force_tolerance(config, magnitude, expected)?))
    };
    if let Some(force) = law.affine_sliding_traction_n(normal, free_slip.abs(), compliance, dt)? {
        limit("friction traction", force, config.maximum_force_n)?;
        if force < required.abs() {
            let (residual, tolerance) = evaluate(force)?;
            if residual.abs() <= tolerance {
                let traction = if force == 0.0 { 0.0 } else { sign * force };
                return Ok((traction, ModalFrictionRegime::Sliding, 0));
            }
        }
    }
    let (zero, _) = evaluate(0.0)?;
    if zero == 0.0 { return Ok((0.0, ModalFrictionRegime::Sliding, 0)); }
    let mut lo = 0.0;
    let mut hi = required.abs().min(config.maximum_force_n);
    let (mut residual, mut tolerance) = evaluate(hi)?;
    if residual.abs() <= tolerance { return Ok((sign*hi, ModalFrictionRegime::Sliding, 0)); }
    if residual < 0.0 {
        return Err(invalid("friction root is not bracketed within the declared force ceiling").into());
    }
    let mut used = 0;
    for iteration in 1..=config.max_iterations {
        used = iteration;
        let mid = f64::midpoint(lo, hi);
        (residual, tolerance) = evaluate(mid)?;
        if residual.abs() <= tolerance { return Ok((sign*mid, ModalFrictionRegime::Sliding, iteration)); }
        if mid == lo || mid == hi { break; }
        if residual < 0.0 { lo = mid; } else { hi = mid; }
    }
    Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: used }.into())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn check(law: &FrictionLaw, config: ModalFrictionConfig, dt: f64,
    drive: FrictionDrive, traction: f64, slip: f64, regime: ModalFrictionRegime,
    iterations: usize) -> Result<(f64, f64), ModalFrictionError>
{
    let speed = finite(slip/dt)?;
    let mut constitutive_residual_n = 0.0;
    let loss = match regime {
        ModalFrictionRegime::Shearing { .. } | ModalFrictionRegime::Microslip => return Err(invalid("material shear must be checked by its stateful constitutive owner").into()),
        ModalFrictionRegime::Released => 0.0,
        ModalFrictionRegime::Sticking => {
            if slip.abs() > config.sticking_tolerance_m {
                return Err(ModalFrictionError::Sticking { residual_m: slip,
                    tolerance_m: config.sticking_tolerance_m });
            }
            0.0
        }
        ModalFrictionRegime::Sliding => {
            limit("friction slip speed", speed.abs(), config.maximum_slip_speed_m_s)?;
            if speed == 0.0 || (traction != 0.0 && traction.signum() != speed.signum()) {
                return Err(invalid("actual sliding force must oppose nonzero relative slip").into());
            }
            let expected = finite(law.kinetic_coefficient(speed.abs())?
                *drive.normal_force_n*speed.signum())?;
            constitutive_residual_n = finite(traction-expected)?;
            let tolerance = force_tolerance(config, traction, expected)?;
            if constitutive_residual_n.abs() > tolerance {
                return Err(ModalCouplingError::ContactSolve { residual_n: constitutive_residual_n,
                    tolerance_n: tolerance, iterations }.into());
            }
            // Applied traction, not a nearby law value: this closes the
            // actual discrete work exchanged with the original stepper.
            finite(traction*slip)?
        }
    };
    Ok((loss, constitutive_residual_n))
}

#[cfg(test)]
mod trait_tests {
    use super::*;
    use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

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
}
