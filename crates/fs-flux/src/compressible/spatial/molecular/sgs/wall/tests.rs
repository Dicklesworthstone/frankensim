use super::*;

fn config() -> EquilibriumWallConfig {
    EquilibriumWallConfig {
        maximum_wall_distance_plus: 1000.0,
        maximum_mach_number: 0.3,
    }
}

fn wall(velocity: [f64; 3]) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity,
        thermal: WallThermalCondition::Adiabatic,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn equilibrium_wall_recovers_manufactured_spalding_friction_across_inner_layers() {
    // Manufacture physical samples from the defining profile, independently
    // of the production inversion. A fixed friction velocity must be recovered
    // at every sample height, from the viscous to the logarithmic layer.
    let (rho, mu, friction_velocity) = (1.2, 1.8e-5, 0.01);
    for u_plus in [0.1, 1.0, 5.0, 10.0, 20.0] {
        let ku = 0.41_f64 * u_plus;
        let y_plus = u_plus + (ku.exp() - 1.0 - ku - ku.powi(2) / 2.0 - ku.powi(3) / 6.0) / 8.6;
        let distance = y_plus * mu / (rho * friction_velocity);
        let speed = u_plus * friction_velocity;
        let viscosity = config()
            .viscosity(
                wall([0.0; 3]),
                1,
                [speed, 0.0, 0.0],
                rho,
                300.0,
                distance,
                mu,
            )
            .unwrap();
        close(
            viscosity.secant * speed / distance,
            rho * friction_velocity.powi(2),
            5e-12,
        );
        assert!(viscosity.tangent >= viscosity.secant);
        assert!(viscosity.tangent <= 2.0 * viscosity.secant);

        // A numerical traction derivative checks the nonlinear stability
        // coefficient, not just the secant used to compute the actual flux.
        let epsilon = 1e-5 * speed;
        let stress = |velocity| {
            let value = config()
                .viscosity(
                    wall([0.0; 3]),
                    1,
                    [velocity, 0.0, 0.0],
                    rho,
                    300.0,
                    distance,
                    mu,
                )
                .unwrap();
            value.secant * velocity / distance
        };
        let derivative =
            distance * (stress(speed + epsilon) - stress(speed - epsilon)) / (2.0 * epsilon);
        close(viscosity.tangent, derivative, 2e-8);
    }
}

#[test]
fn equilibrium_wall_preserves_rest_viscous_limit_and_wall_relative_frame() {
    let mu = 1e-4;
    let rest = config()
        .viscosity(
            wall([0.25, 0.0, 0.0]),
            1,
            [0.25, 0.0, 0.0],
            1.0,
            10.0,
            0.1,
            mu,
        )
        .unwrap();
    assert_eq!(rest.secant.to_bits(), mu.to_bits());
    assert_eq!(rest.tangent.to_bits(), mu.to_bits());
    let tiny = config()
        .viscosity(wall([0.0; 3]), 1, [1e-12, 0.0, 0.0], 1.0, 10.0, 0.1, mu)
        .unwrap();
    close(tiny.secant, mu, 1e-12);
    let base = config()
        .viscosity(wall([0.0; 3]), 1, [0.25, 0.0, -0.125], 1.0, 10.0, 0.1, mu)
        .unwrap();
    let shifted = config()
        .viscosity(
            wall([8.0, 0.0, -4.0]),
            1,
            [8.25, 0.0, -4.125],
            1.0,
            10.0,
            0.1,
            mu,
        )
        .unwrap();
    assert_eq!(base.secant.to_bits(), shifted.secant.to_bits());
    assert_eq!(base.tangent.to_bits(), shifted.tangent.to_bits());
}

#[test]
fn equilibrium_wall_refuses_incompatible_wall_and_authored_domain_exit() {
    assert!(config().admit_wall(MolecularWall3d::SlipAdiabatic).is_err());
    assert!(
        config()
            .admit_wall(MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
            })
            .is_err()
    );
    for value in [0.0, -1.0, 1001.0, f64::NAN, f64::INFINITY] {
        assert!(
            EquilibriumWallConfig {
                maximum_wall_distance_plus: value,
                ..config()
            }
            .validate()
            .is_err()
        );
    }
    assert!(
        config()
            .viscosity(wall([0.0; 3]), 1, [0.0, 0.31, 0.0], 1.0, 1.0, 0.1, 1e-4)
            .unwrap_err()
            .contains("Mach")
    );
    let limited = EquilibriumWallConfig {
        maximum_wall_distance_plus: 1.0,
        ..config()
    };
    assert!(
        limited
            .viscosity(wall([0.0; 3]), 1, [0.2, 0.0, 0.0], 1.0, 10.0, 0.1, 1e-4)
            .unwrap_err()
            .contains("distance-plus")
    );
    assert!(
        config()
            .viscosity(wall([0.0; 3]), 1, [0.0; 3], 1.0, 10.0, 0.1, 0.0)
            .is_err()
    );
}
