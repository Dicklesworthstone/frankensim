//! Native buoyancy quadrature under actual wall transport and checked sources.
use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::{E3, ThermalLbm3};
use fs_lbm::thermal::boundary::ThermalWall;
use fs_lbm::thermal::sources::transport::RelaxationField;
use fs_lbm::thermal::sources::{ThermalSourceError, ThermalSourceLimits};
use fs_lbm::thermal::{BuoyancyQuadrature, ThermalForcingError};

const W5: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];
const W7: [f64; 7] = [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125];
const E2: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: 1024,
        temperature_range: [-2.0, 3.0],
        balance_tolerance: 1e-10,
    }
}

// A two-layer channel with equal imposed flux at both walls stays spatially
// isothermal by reflection symmetry, while actual scalar streaming raises T.
// With dx=rho=duration=1, forces scale with dt^2 and velocity with dt. This
// permits time refinement at fixed viscosity/diffusivity and physical loads.
fn heating_impulse_error(steps: usize, quadrature: BuoyancyQuadrature) -> [f64; 2] {
    let dt = 1.0 / steps as f64;
    let (initial, reference, wall_rate, volume_rate, gravity_beta) = (0.9, 0.5, 0.06, 0.04, 0.001);
    let tau_f = 0.5 + 3.0 * 0.2 * dt;
    let tau_g2 = 0.5 + 3.0 * 0.1 * dt;
    let tau_g3 = 0.5 + 4.0 * 0.1 * dt;
    let mut two = ThermalLbm::slab(1, 2, tau_f, tau_g2, gravity_beta * dt * dt);
    let mut three = ThermalLbm3::slab(1, 1, 2, tau_f, tau_g3, gravity_beta * dt * dt);
    assert_eq!(
        two.buoyancy_quadrature(),
        BuoyancyQuadrature::OpeningTemperature
    );
    assert_eq!(
        three.buoyancy_quadrature(),
        BuoyancyQuadrature::OpeningTemperature
    );
    let direction2 = [-1.0, 0.0];
    let direction3 = [0.6, 0.8, 0.0];
    two.set_buoyancy_direction(direction2).unwrap();
    three.set_buoyancy_direction(direction3).unwrap();
    two.set_buoyancy_quadrature(quadrature);
    three.set_buoyancy_quadrature(quadrature);
    two.gpop.fill(W5.map(|weight| initial * weight));
    three.gpop.fill(W7.map(|weight| initial * weight));
    two.t_ref = reference;
    three.t_ref = reference;
    let walls = [ThermalWall::HeatFlux(wall_rate * dt); 2];
    two.set_thermal_walls(walls).unwrap();
    three.set_thermal_walls(walls).unwrap();
    let external2 = [0.0002, 0.0];
    let external3 = [-0.00015, 0.00005, 0.0];
    let forces2 = [external2.map(|value| value * dt * dt); 4];
    let forces3 = [external3.map(|value| value * dt * dt); 4];
    let heat = [0.0, volume_rate * dt, volume_rate * dt, 0.0];
    let momentum_tau = [tau_f; 4];
    let thermal_tau2 = [tau_g2; 4];
    let thermal_tau3 = [tau_g3; 4];
    let mut wall_impulse2 = [0.0; 2];
    let mut wall_impulse3 = [0.0; 3];
    for step in 0..steps {
        let (balance2, balance3) = if step % 2 == 0 {
            (
                two.step_with_sources(&heat, &forces2, limits()).unwrap(),
                three.step_with_sources(&heat, &forces3, limits()).unwrap(),
            )
        } else {
            (
                two.step_with_sources_and_transport(
                    &heat,
                    &forces2,
                    RelaxationField {
                        momentum: &momentum_tau,
                        thermal: &thermal_tau2,
                    },
                    limits(),
                )
                .unwrap(),
                three
                    .step_with_sources_and_transport(
                        &heat,
                        &forces3,
                        RelaxationField {
                            momentum: &momentum_tau,
                            thermal: &thermal_tau3,
                        },
                        limits(),
                    )
                    .unwrap(),
            )
        };
        assert_eq!(two.buoyancy_quadrature(), quadrature);
        assert_eq!(three.buoyancy_quadrature(), quadrature);
        for balance in [balance2, balance3] {
            close(balance.source, 2.0 * volume_rate * dt, 2e-13);
            close(balance.bottom_wall, wall_rate * dt, 2e-13);
            close(balance.top_wall, wall_rate * dt, 2e-13);
        }
        for cell in 1..3 {
            let temperature = initial + (wall_rate + volume_rate) * (step + 1) as f64 * dt;
            close(two.gpop[cell].iter().sum(), temperature, 2e-13);
            close(three.gpop[cell].iter().sum(), temperature, 2e-13);
        }
        // Independently measure the stationary-wall impulse from the actual
        // reflected populations: outgoing momentum minus incoming momentum.
        // No force field or volume-balance residual enters this measurement.
        for (cell, inward) in [(1, 1), (2, -1)] {
            for (q, &(ex, ey)) in E2.iter().enumerate() {
                if ey == inward {
                    wall_impulse2[0] -= 2.0 * two.grid.f[cell][q] * f64::from(ex);
                    wall_impulse2[1] -= 2.0 * two.grid.f[cell][q] * f64::from(ey);
                }
            }
            for (q, &(ex, ey, ez)) in E3.iter().enumerate() {
                if i32::from(ez) == inward {
                    for (axis, component) in [ex, ey, ez].iter().enumerate() {
                        wall_impulse3[axis] -= 2.0 * three.grid.f[cell][q] * f64::from(*component);
                    }
                }
            }
        }
    }
    let integral = gravity_beta * (initial - reference + 0.5 * (wall_rate + volume_rate));
    let mut momentum2 = [0.0; 2];
    let mut momentum3 = [0.0; 3];
    for cell in 1..3 {
        let moments = two.grid.moments(cell);
        let (rho3, velocity3) = three.grid.moments3(cell);
        for axis in 0..2 {
            momentum2[axis] += moments.rho * moments.u[axis];
        }
        for axis in 0..3 {
            momentum3[axis] += rho3 * velocity3[axis];
        }
    }
    let mut error = [0.0_f64; 2];
    for axis in 0..2 {
        let measured = (momentum2[axis] + wall_impulse2[axis]) / dt;
        let expected = 2.0 * (integral * direction2[axis] + external2[axis]);
        error[0] = error[0].max((measured - expected).abs());
    }
    for axis in 0..3 {
        let measured = (momentum3[axis] + wall_impulse3[axis]) / dt;
        let expected = 2.0 * (integral * direction3[axis] + external3[axis]);
        error[1] = error[1].max((measured - expected).abs());
    }
    error
}

#[test]
fn transported_wall_and_volume_heat_deliver_the_endpoint_averaged_impulse() {
    let mut previous = [0.0; 2];
    for steps in [4, 8, 16] {
        let opening = heating_impulse_error(steps, BuoyancyQuadrature::OpeningTemperature);
        let averaged = heating_impulse_error(steps, BuoyancyQuadrature::EndpointAverage);
        println!("steps={steps}: opening impulse error {opening:?}; endpoint average {averaged:?}");
        // Source half-kicks already center volume heat. Opening sampling
        // misses half of this interval's WALL heat, producing an O(dt) lag.
        close(opening[0], 0.001 * 0.06 / steps as f64, 2e-12);
        close(opening[1], 0.8 * 0.001 * 0.06 / steps as f64, 2e-12);
        for dimension in 0..2 {
            assert!(averaged[dimension] < 2e-12, "{averaged:?}");
            if steps > 4 {
                close(previous[dimension] / opening[dimension], 2.0, 2e-6);
            }
        }
        previous = opening;
    }
}

#[test]
fn endpoint_force_changes_retain_one_heat_step_at_the_incoming_physical_velocity() {
    let mut two = ThermalLbm::slab(5, 6, 0.9, 0.8, 2e-4);
    let mut reference2 = ThermalLbm::slab(5, 6, 0.9, 0.8, 2e-4);
    let mut three = ThermalLbm3::slab(5, 2, 6, 0.9, 0.9, 2e-4);
    let mut reference3 = ThermalLbm3::slab(5, 2, 6, 0.9, 0.9, 2e-4);
    for slab in [&mut two, &mut reference2] {
        slab.perturb(0.003);
    }
    for slab in [&mut three, &mut reference3] {
        slab.perturb(0.003);
    }
    for _ in 0..2 {
        for slab in [&mut two, &mut reference2] {
            slab.step_with_external_forces(&vec![[0.0002, -0.0001]; slab.gpop.len()])
                .unwrap();
        }
        for slab in [&mut three, &mut reference3] {
            slab.step_with_external_forces(&vec![[-0.0001, 0.00015, 0.00005]; slab.gpop.len()])
                .unwrap();
        }
    }
    two.set_buoyancy_quadrature(BuoyancyQuadrature::EndpointAverage);
    three.set_buoyancy_quadrature(BuoyancyQuadrature::EndpointAverage);
    reference2.gbeta = 0.0;
    reference3.gbeta = 0.0;
    let zero2 = vec![[0.0; 2]; two.gpop.len()];
    let zero3 = vec![[0.0; 3]; three.gpop.len()];
    let heat2 = vec![0.0; two.gpop.len()];
    let heat3 = vec![0.0; three.gpop.len()];
    for (direction2, direction3, sign) in [
        ([1.0, 0.0], [0.6, 0.8, 0.0], 1.0),
        ([-0.6, 0.8], [0.0, -0.6, 0.8], -1.0),
        ([0.0, -1.0], [-1.0, 0.0, 0.0], 0.0),
    ] {
        // Give the ordinary scalar reference the same evolved input flow,
        // but preserve both independently evolved thermal populations/history.
        reference2.grid = two.grid.clone();
        reference3.grid.f.clone_from(&three.grid.f);
        reference3.grid.fext.clone_from(&three.grid.fext);
        reference3.grid.tau.clone_from(&three.grid.tau);
        two.set_buoyancy_direction(direction2).unwrap();
        three.set_buoyancy_direction(direction3).unwrap();
        reference2.step_with_external_forces(&zero2).unwrap();
        reference3.step_with_external_forces(&zero3).unwrap();
        two.step_with_sources(
            &heat2,
            &vec![[sign * 0.0003, -sign * 0.0001]; heat2.len()],
            limits(),
        )
        .unwrap();
        three
            .step_with_sources(
                &heat3,
                &vec![[sign * 0.0001, -sign * 0.0002, sign * 0.0003]; heat3.len()],
                limits(),
            )
            .unwrap();
        for (actual, expected) in two
            .gpop
            .iter()
            .flatten()
            .zip(reference2.gpop.iter().flatten())
        {
            close(*actual, *expected, 2e-13);
        }
        for (actual, expected) in three
            .gpop
            .iter()
            .flatten()
            .zip(reference3.gpop.iter().flatten())
        {
            close(*actual, *expected, 2e-13);
        }
        for (actual, expected) in two
            .last_wall_exchange()
            .unwrap()
            .iter()
            .zip(reference2.last_wall_exchange().unwrap())
        {
            close(actual.heat_in, expected.heat_in, 2e-13);
        }
        for (actual, expected) in three
            .last_wall_exchange()
            .unwrap()
            .iter()
            .zip(reference3.last_wall_exchange().unwrap())
        {
            close(actual.heat_in, expected.heat_in, 2e-13);
        }
    }
}

#[test]
fn endpoint_force_overflow_refuses_the_complete_source_candidate() {
    let mut two = ThermalLbm::slab(1, 2, 0.9, 0.8, f64::MAX);
    let mut three = ThermalLbm3::slab(1, 1, 2, 0.9, 0.9, f64::MAX);
    two.gpop.fill([0.0; 5]);
    three.gpop.fill([0.0; 7]);
    two.t_ref = 0.0;
    three.t_ref = 0.0;
    two.t_bottom = 12.0;
    three.t_bottom = 12.0;
    two.set_buoyancy_direction([1.0, 0.0]).unwrap();
    three.set_buoyancy_direction([0.0, -1.0, 0.0]).unwrap();
    two.set_buoyancy_quadrature(BuoyancyQuadrature::EndpointAverage);
    three.set_buoyancy_quadrature(BuoyancyQuadrature::EndpointAverage);
    let before2 = (
        two.grid.f.clone(),
        two.gpop.clone(),
        two.grid.fext.clone(),
        two.last_wall_exchange(),
    );
    let before3 = (
        three.grid.f.clone(),
        three.gpop.clone(),
        three.grid.fext.clone(),
        three.last_wall_exchange(),
    );
    let wide = ThermalSourceLimits {
        temperature_range: [-20.0, 20.0],
        ..limits()
    };
    // Opening temperature and force are exactly zero. Only the transported
    // endpoint makes gbeta * mean(T) overflow; refusal must discard that heat,
    // its wall receipt and history as well as all flow/force changes.
    assert!(matches!(
        two.step_with_sources(&[0.0; 4], &[[0.0; 2]; 4], wide),
        Err(ThermalSourceError::Force(
            ThermalForcingError::NonFinite { .. }
        ))
    ));
    assert!(matches!(
        three.step_with_sources(&[0.0; 4], &[[0.0; 3]; 4], wide),
        Err(ThermalSourceError::Force(
            ThermalForcingError::NonFinite { .. }
        ))
    ));
    assert_eq!(
        (
            two.grid.f.clone(),
            two.gpop.clone(),
            two.grid.fext.clone(),
            two.last_wall_exchange()
        ),
        before2
    );
    assert_eq!(
        (
            three.grid.f.clone(),
            three.gpop.clone(),
            three.grid.fext.clone(),
            three.last_wall_exchange()
        ),
        before3
    );
    assert_eq!(
        two.buoyancy_quadrature(),
        BuoyancyQuadrature::EndpointAverage
    );
    assert_eq!(
        three.buoyancy_quadrature(),
        BuoyancyQuadrature::EndpointAverage
    );
    // The admitted follow-up uses the original thermal state, including a
    // fresh derivative history; no refused wall heating can leak into it.
    two.gbeta = 0.0;
    three.gbeta = 0.0;
    let balance2 = two
        .step_with_sources(&[0.0; 4], &[[0.0; 2]; 4], wide)
        .unwrap();
    let balance3 = three
        .step_with_sources(&[0.0; 4], &[[0.0; 3]; 4], wide)
        .unwrap();
    close(balance2.stored_change, 4.0, 2e-13);
    close(balance3.stored_change, 3.0, 2e-13);
}
