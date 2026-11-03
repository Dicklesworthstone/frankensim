use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn law(a: f64, b: f64) -> RecoilLaw {
    RecoilLaw {
        linear_rate_per_s: a,
        quadratic_per_m: b,
        gas_heat_fraction: 1.0,
    }
}

#[test]
fn reduced_mass_recoil_matches_independent_slip_oracle_and_subdivision() {
    let (density, volume, dt, velocity) = (1.5, 0.4, 0.2, 1.1);
    let support = SlidingSupport {
        mass_kg: 0.8,
        velocity_m_s: -0.3,
    };
    for (a, b) in [(0.0, 0.0), (1.7, 0.0), (0.0, 2.3), (1.7, 2.3)] {
        let law = law(a, b);
        let mass = density * volume;
        let reduced = mass * support.mass_kg / (mass + support.mass_kg);
        let elapsed = dt * (1.0 + mass / support.mass_kg);
        let slip: f64 = velocity - support.velocity_m_s;
        let denominator = if a == 0.0 {
            1.0 + b * slip.abs() * elapsed
        } else {
            (a * elapsed).exp() + b * slip.abs() / a * (a * elapsed).exp_m1()
        };
        let final_slip = slip / denominator;
        let center =
            (mass * velocity + support.mass_kg * support.velocity_m_s) / (mass + support.mass_kg);
        let result = advance(
            density * velocity,
            density,
            velocity,
            volume,
            dt,
            support,
            law,
        )
        .unwrap();
        close(
            result.momentum / density,
            center + reduced / mass * final_slip,
            5e-15,
        );
        close(
            result.support_velocity_m_s,
            center - reduced / support.mass_kg * final_slip,
            5e-15,
        );
        close(
            result.dissipation_j,
            0.5 * reduced * (slip * slip - final_slip * final_slip),
            3e-15,
        );
        close(result.momentum_residual_ns, 0.0, 3e-15);
        close(
            result.work_to_gas_j + result.support_kinetic_energy_change_j,
            0.0,
            3e-15,
        );
        let half = advance(
            density * velocity,
            density,
            velocity,
            volume,
            dt / 2.0,
            support,
            law,
        )
        .unwrap();
        let split = advance(
            half.momentum,
            density,
            half.momentum / density,
            volume,
            dt / 2.0,
            SlidingSupport {
                velocity_m_s: half.support_velocity_m_s,
                ..support
            },
            law,
        )
        .unwrap();
        close(split.momentum, result.momentum, 5e-15);
        close(
            split.support_velocity_m_s,
            result.support_velocity_m_s,
            5e-15,
        );
        close(
            half.dissipation_j + split.dissipation_j,
            result.dissipation_j,
            4e-15,
        );
        close(
            half.work_to_gas_j + split.work_to_gas_j,
            result.work_to_gas_j,
            4e-15,
        );
    }
}

#[test]
fn strong_recoil_and_large_support_limit_use_actual_mass_and_initial_velocity() {
    for mass in [1e-6, 1.0, 1e6, 1e12] {
        let support = SlidingSupport {
            mass_kg: mass,
            velocity_m_s: 0.0,
        };
        let result = advance(2.0, 1.0, 2.0, 1.0, 1.0, support, law(1e4, 2.0)).unwrap();
        let center = 2.0 / (1.0 + mass);
        close(result.momentum, center, 3e-15);
        close(result.support_velocity_m_s, center, 3e-10);
        close(result.momentum_residual_ns, 0.0, 2e-15);
        close(result.dissipation_j, 2.0 * mass / (1.0 + mass), 5e-15);
    }
    let support = SlidingSupport {
        mass_kg: 1e12,
        velocity_m_s: 0.0,
    };
    let result = advance(2.0, 1.0, 2.0, 1.0, 0.2, support, law(1.7, 2.3)).unwrap();
    let (fixed, _) = super::super::super::drag::combined_decay(1.7, 2.3, 2.0, 0.2);
    close(result.momentum, 2.0 * fixed, 3e-12);
    assert!(result.support_velocity_m_s > 0.0);
    assert!(
        advance(
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            SlidingSupport {
                mass_kg: f64::from_bits(1),
                velocity_m_s: 0.0
            },
            law(1.0, 1.0)
        )
        .is_err()
    );
}
