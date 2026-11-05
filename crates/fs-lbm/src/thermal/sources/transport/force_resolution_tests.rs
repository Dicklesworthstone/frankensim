//! Actual native heat/flow under a nonlinear force query after transport.
use super::*;
use crate::d3q19::E3;
use crate::thermal::ThermalForcingError;
use crate::thermal::boundary::ThermalWall;
use crate::thermal::sources::Slab;

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: 128,
        temperature_range: [-2.0, 3.0],
        balance_tolerance: 1e-10,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn density_departure(temperature: f64) -> f64 {
    let rise = temperature - 0.9;
    -0.02 * rise - 0.04 * rise * rise
}

// Equal heat flux at both walls of a two-layer periodic channel produces the
// exact thermal ramp T(t)=0.9+0.1t. A quadratic density law gives a quadratic
// force in time, whose trapezoidal quadrature error is known independently.
// dx=rho0=duration=1; force and velocity scale as dt² and dt, respectively.
fn nonlinear_impulse_error(steps: usize) -> [f64; 2] {
    let dt = 1.0 / steps as f64;
    let mut two = ThermalLbm::slab(1, 2, 0.5 + 0.6 * dt, 0.5 + 0.3 * dt, 0.0);
    let mut three = ThermalLbm3::slab(1, 1, 2, 0.5 + 0.6 * dt, 0.5 + 0.4 * dt, 0.0);
    two.gpop.fill([0.3, 0.15, 0.15, 0.15, 0.15]);
    three
        .gpop
        .fill([0.225, 0.1125, 0.1125, 0.1125, 0.1125, 0.1125, 0.1125]);
    let walls = [ThermalWall::HeatFlux(0.06 * dt); 2];
    two.set_thermal_walls(walls).unwrap();
    three.set_thermal_walls(walls).unwrap();
    let gravity2 = [-0.1, 0.0];
    let gravity3 = [0.06, 0.08, 0.0];
    let external2 = [0.0002, 0.0];
    let external3 = [-0.00015, 0.00005, 0.0];
    let heat = [0.0, 0.04 * dt, 0.04 * dt, 0.0];
    let momentum_tau = two.grid.tau.clone();
    let thermal_tau2 = [two.tau_g; 4];
    let thermal_tau3 = [three.tau_g; 4];
    let mut wall2 = [0.0; 2];
    let mut wall3 = [0.0; 3];
    for step in 0..steps {
        let initial2: Vec<_> = two
            .gpop
            .iter()
            .map(|p| density_departure(p.iter().sum()))
            .collect();
        let initial3: Vec<_> = three
            .gpop
            .iter()
            .map(|p| density_departure(p.iter().sum()))
            .collect();
        let forces2: Vec<_> = initial2
            .iter()
            .map(|&rho| {
                core::array::from_fn(|axis| (external2[axis] + rho * gravity2[axis]) * dt * dt)
            })
            .collect();
        let forces3: Vec<_> = initial3
            .iter()
            .map(|&rho| {
                core::array::from_fn(|axis| (external3[axis] + rho * gravity3[axis]) * dt * dt)
            })
            .collect();
        let before2 = (two.grid.f.clone(), two.grid.fext.clone());
        let before3 = (three.grid.f.clone(), three.grid.fext.clone());
        let mut calls2 = 0;
        let mut calls3 = 0;
        let next_temperature = 0.9 + 0.1 * (step + 1) as f64 * dt;
        let balance2 = source_step_resolved_mapped::<_, 2, 5, _, ThermalSourceError>(
            &mut two,
            &heat,
            &forces2,
            limits(),
            RelaxationField {
                momentum: &momentum_tau,
                thermal: &thermal_tau2,
            },
            |transported, predicted| {
                calls2 += 1;
                assert_eq!(
                    (&transported.grid.f, &transported.grid.fext),
                    (&before2.0, &before2.1)
                );
                let mut resolved = forces2.clone();
                for cell in 1..3 {
                    close(predicted[cell], next_temperature, 2e-13);
                    let mean = 0.5 * initial2[cell] + 0.5 * density_departure(predicted[cell]);
                    resolved[cell] = core::array::from_fn(|axis| {
                        (external2[axis] + mean * gravity2[axis]) * dt * dt
                    });
                }
                Ok(resolved)
            },
            |_, balance| Ok(balance),
        )
        .unwrap();
        let balance3 = source_step_resolved_mapped::<_, 3, 7, _, ThermalSourceError>(
            &mut three,
            &heat,
            &forces3,
            limits(),
            RelaxationField {
                momentum: &momentum_tau,
                thermal: &thermal_tau3,
            },
            |transported, predicted| {
                calls3 += 1;
                assert_eq!(
                    (&transported.grid.f, &transported.grid.fext),
                    (&before3.0, &before3.1)
                );
                let mut resolved = forces3.clone();
                for cell in 1..3 {
                    close(predicted[cell], next_temperature, 2e-13);
                    let mean = 0.5 * initial3[cell] + 0.5 * density_departure(predicted[cell]);
                    resolved[cell] = core::array::from_fn(|axis| {
                        (external3[axis] + mean * gravity3[axis]) * dt * dt
                    });
                }
                Ok(resolved)
            },
            |_, balance| Ok(balance),
        )
        .unwrap();
        assert_eq!((calls2, calls3), (1, 1));
        for balance in [balance2, balance3] {
            close(balance.source, 0.08 * dt, 2e-13);
            close(balance.bottom_wall, 0.06 * dt, 2e-13);
            close(balance.top_wall, 0.06 * dt, 2e-13);
            close(balance.stored_change, 0.2 * dt, 2e-13);
        }
        for cell in 1..3 {
            close(two.gpop[cell].iter().sum(), next_temperature, 2e-13);
            close(three.gpop[cell].iter().sum(), next_temperature, 2e-13);
        }
        // Independently measure actual reflected-link momentum, including
        // wall drag. No resolved force or heat residual is used as its oracle.
        for (cell, inward) in [(1, 1), (2, -1)] {
            for (q, &(ex, ey)) in crate::E.iter().enumerate() {
                if i32::from(ey) == inward {
                    for (axis, component) in [ex, ey].iter().enumerate() {
                        wall2[axis] -= 2.0 * two.grid.f[cell][q] * f64::from(*component);
                    }
                }
            }
            for (q, &(ex, ey, ez)) in E3.iter().enumerate() {
                if i32::from(ez) == inward {
                    for (axis, component) in [ex, ey, ez].iter().enumerate() {
                        wall3[axis] -= 2.0 * three.grid.f[cell][q] * f64::from(*component);
                    }
                }
            }
        }
    }
    let density_integral = -0.02 * 0.1 / 2.0 - 0.04 * 0.1 * 0.1 / 3.0;
    let mut momentum2 = [0.0; 2];
    let mut momentum3 = [0.0; 3];
    for cell in 1..3 {
        let moments = two.grid.moments(cell);
        let (rho, velocity) = three.grid.moments3(cell);
        for axis in 0..2 {
            momentum2[axis] += moments.rho * moments.u[axis];
        }
        for axis in 0..3 {
            momentum3[axis] += rho * velocity[axis];
        }
    }
    let mut errors = [0.0_f64; 2];
    for axis in 0..2 {
        let measured = (momentum2[axis] + wall2[axis]) / dt;
        let expected = 2.0 * (density_integral * gravity2[axis] + external2[axis]);
        errors[0] = errors[0].max((measured - expected).abs());
    }
    for axis in 0..3 {
        let measured = (momentum3[axis] + wall3[axis]) / dt;
        let expected = 2.0 * (density_integral * gravity3[axis] + external3[axis]);
        errors[1] = errors[1].max((measured - expected).abs());
    }
    errors
}

#[test]
fn nonlinear_density_force_uses_one_thermal_endpoint_and_one_flow_impulse() {
    let mut previous = [0.0; 2];
    for steps in [4, 8, 16] {
        let errors = nonlinear_impulse_error(steps);
        println!("nonlinear density quadrature, steps={steps}: {errors:?}");
        close(
            errors[0],
            2.0 * 0.1 * 0.04 * 0.1 * 0.1 / (6.0 * (steps * steps) as f64),
            2e-12,
        );
        close(errors[1], 0.8 * errors[0], 2e-12);
        if steps > 4 {
            for axis in 0..2 {
                close(previous[axis] / errors[axis], 4.0, 3e-4);
            }
        }
        previous = errors;
    }
}

#[derive(Debug, PartialEq)]
enum QueryError {
    Thermal(ThermalSourceError),
    Owner(&'static str),
}
impl From<ThermalSourceError> for QueryError {
    fn from(error: ThermalSourceError) -> Self {
        Self::Thermal(error)
    }
}
impl From<ThermalForcingError> for QueryError {
    fn from(error: ThermalForcingError) -> Self {
        Self::Thermal(error.into())
    }
}

fn rollback<S: Slab<D, Q>, const D: usize, const Q: usize>(
    mut state: S,
    mut reference: S,
    force: [f64; D],
    thermal_tau: f64,
    flow: impl Fn(&S) -> Vec<f64>,
) {
    let (cells, plane) = state.layout().unwrap();
    let forces = vec![force; cells];
    for _ in 0..4 {
        state.advance(&forces).unwrap();
        reference.advance(&forces).unwrap();
    }
    let before = (
        flow(&state),
        state.populations().to_vec(),
        state.wall_exchange(),
    );
    let momentum = state.momentum_relaxation().to_vec();
    let thermal = vec![thermal_tau; cells];
    let rates = RelaxationField {
        momentum: &momentum,
        thermal: &thermal,
    };
    let mut heat = vec![0.0; cells];
    heat[plane..cells - plane].fill(0.002);
    for stage in 0..3 {
        let mut queries = 0;
        let mut finishes = 0;
        let refusal = source_step_resolved_mapped::<S, D, Q, _, QueryError>(
            &mut state,
            &heat,
            &forces,
            limits(),
            rates,
            |_, _| {
                queries += 1;
                if stage == 0 {
                    return Err(QueryError::Owner("force source rejected"));
                }
                let mut resolved = forces.clone();
                if stage == 2 {
                    resolved[plane][0] = f64::NAN;
                }
                Ok(resolved)
            },
            |candidate, _| {
                finishes += 1;
                assert_ne!(candidate.populations(), before.1.as_slice());
                Err::<(), _>(QueryError::Owner("actual endpoint source rejected"))
            },
        );
        assert_eq!(queries, 1);
        match stage {
            0 => assert_eq!(refusal, Err(QueryError::Owner("force source rejected"))),
            1 => assert_eq!(
                refusal,
                Err(QueryError::Owner("actual endpoint source rejected"))
            ),
            _ => assert!(matches!(
                refusal,
                Err(QueryError::Thermal(ThermalSourceError::Force(
                    ThermalForcingError::NonFinite { .. }
                )))
            )),
        }
        assert_eq!(finishes, usize::from(stage == 1));
        assert_eq!(
            (
                flow(&state),
                state.populations().to_vec(),
                state.wall_exchange()
            ),
            before
        );
    }
    // Continue against an independently warmed ordinary solver. Both the
    // early and late refused candidates must leave accepted history intact.
    let zero = vec![0.0; cells];
    for _ in 0..3 {
        reference.advance_local(&forces, &thermal).unwrap();
        source_step_resolved_mapped::<S, D, Q, _, QueryError>(
            &mut state,
            &zero,
            &forces,
            limits(),
            rates,
            |_, _| Ok(forces.clone()),
            |_, balance| Ok(balance),
        )
        .unwrap();
        assert_eq!(flow(&state), flow(&reference));
        assert_eq!(state.populations(), reference.populations());
        assert_eq!(state.wall_exchange(), reference.wall_exchange());
    }
}

#[test]
fn typed_force_and_actual_endpoint_refusals_preserve_native_state_and_history() {
    let mut two = ThermalLbm::slab(5, 6, 0.9, 0.8, 0.0);
    let mut reference2 = ThermalLbm::slab(5, 6, 0.9, 0.8, 0.0);
    two.perturb(0.003);
    reference2.perturb(0.003);
    rollback(two, reference2, [0.0002, -0.0001], 0.8, |s| {
        s.grid
            .f
            .iter()
            .flatten()
            .chain(s.grid.fext.iter().flatten())
            .chain(&s.grid.tau)
            .copied()
            .collect()
    });
    let mut three = ThermalLbm3::slab(3, 2, 6, 0.9, 0.9, 0.0);
    let mut reference3 = ThermalLbm3::slab(3, 2, 6, 0.9, 0.9, 0.0);
    let bulk: Vec<_> = (0..three.grid.f.len())
        .map(|cell| 0.7 + 0.2 * (cell % 5) as f64)
        .collect();
    three
        .grid
        .set_central_moment_bulk_field(&bulk, 1.2)
        .unwrap();
    reference3
        .grid
        .set_central_moment_bulk_field(&bulk, 1.2)
        .unwrap();
    three.perturb(0.003);
    reference3.perturb(0.003);
    rollback(three, reference3, [0.0001, -0.0002, 0.00015], 0.9, |s| {
        // Spatial bulk relaxation is part of the retained source snapshot,
        // including after force-query and actual-endpoint refusals.
        assert!(s.grid.has_independent_bulk_collision());
        for (cell, &rate) in bulk.iter().enumerate() {
            assert!(matches!(s.grid.collision_model_at(cell),
                crate::d3q19::CollisionModel3::CentralMomentBulk { bulk_rate, higher_order_rate, .. }
                    if bulk_rate == rate && higher_order_rate == 1.2));
        }
        s.grid
            .f
            .iter()
            .flatten()
            .chain(s.grid.fext.iter().flatten())
            .chain(&s.grid.tau)
            .copied()
            .collect()
    });
}
