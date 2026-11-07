//! Source-independent moving-face and extensive RK2 arithmetic.
//!
//! Original constitutive owners admit every input/output state. Geometry uses
//! the same face velocities as these fluxes, so volume motion is not a second
//! independently rounded transport model.

pub(super) fn laboratory_flux(relative: [f64; 3], velocity: f64) -> [f64; 3] {
    [
        relative[0],
        relative[1] + velocity * relative[0],
        relative[2] + velocity * (relative[1] + (0.5 * velocity) * relative[0]),
    ]
}

pub(super) fn predictor(
    state: [f64; 3],
    before_volume: f64,
    after_volume: f64,
    dt_area: f64,
    fluxes: [[f64; 3]; 2],
) -> [f64; 3] {
    // A purely translating uniform state needs no multiply/divide round trip.
    if before_volume == after_volume && fluxes[0] == fluxes[1] {
        return state;
    }
    core::array::from_fn(|k| {
        (before_volume * state[k] + dt_area * (fluxes[0][k] - fluxes[1][k])) / after_volume
    })
}

pub(super) fn corrector(
    before: [f64; 3],
    stage: [f64; 3],
    before_volume: f64,
    after_volume: f64,
    dt_area: f64,
    fluxes: [[f64; 3]; 2],
) -> [f64; 3] {
    if before == stage && before_volume == after_volume && fluxes[0] == fluxes[1] {
        return before;
    }
    core::array::from_fn(|k| {
        (0.5 * before_volume * before[k]
            + 0.5 * (after_volume * stage[k] + dt_area * (fluxes[0][k] - fluxes[1][k])))
            / after_volume
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64, scale: f64) {
        assert!(a.is_finite() && b.is_finite());
        assert!(
            (a - b).abs() <= 128.0 * f64::EPSILON * scale.max(1.0),
            "{a:.16e} != {b:.16e}"
        );
    }

    #[test]
    fn face_frame_flux_equals_independent_laboratory_f_minus_wu() {
        let rho = 2.0;
        let thermal = -19.0; // The common caloric datum may be negative.
        let pressure = 3.0;
        let velocity = 3.0;
        let face = 7.0;
        let relative_velocity = velocity - face;
        let relative_energy = rho * (thermal + 0.5 * relative_velocity * relative_velocity);
        let f = [
            rho * relative_velocity,
            rho * relative_velocity * relative_velocity + pressure,
            relative_velocity * (relative_energy + pressure),
        ];
        let physical_energy = rho * (thermal + 0.5 * velocity * velocity);
        let independent = [
            rho * velocity - face * rho,
            rho * velocity * velocity + pressure - face * rho * velocity,
            velocity * (physical_energy + pressure) - face * physical_energy,
        ];
        let actual = laboratory_flux(f, face);
        for k in 0..3 {
            near(actual[k], independent[k], independent[k].abs());
        }
    }

    #[test]
    fn wall_has_no_mass_transport_and_carries_exact_pressure_work() {
        for (p, w) in [(3.0, -7.0), (1.405e6, 25.0), (1e-200, 1e160)] {
            let flux = laboratory_flux([0.0, p, 0.0], w);
            assert_eq!(flux, [0.0, p, p * w]);
        }
    }

    #[test]
    fn galilean_boosts_compose_for_mass_momentum_and_total_energy_flux() {
        let relative = [-2.3, 8.7, -11.2];
        for (a, b) in [(-13.0, 9.0), (2000.0, -1997.0), (0.0, 3.0)] {
            let once = laboratory_flux(relative, a + b);
            let twice = laboratory_flux(laboratory_flux(relative, a), b);
            for k in 0..3 {
                near(once[k], twice[k], 1e7);
            }
        }
    }

    #[test]
    fn translating_uniform_state_keeps_the_original_payload_bits() {
        let q = [10.3, -31.713, -72903.111];
        let flux = [0.0, 1.405e6, -3.7e6];
        let stage = predictor(q, 0.037, 0.037, 1e-6, [flux; 2]);
        let after = corrector(q, stage, 0.037, 0.037, 1e-6, [flux; 2]);
        assert_eq!(stage.map(f64::to_bits), q.map(f64::to_bits));
        assert_eq!(after.map(f64::to_bits), q.map(f64::to_bits));
    }

    #[test]
    fn moving_volume_gcl_preserves_an_arbitrary_uniform_conservative_state() {
        let q = [2.0, 6.0, 19.0];
        let physical_flux = [6.0, 21.0, 66.0];
        let w = [-0.5, 0.25];
        let fluxes = w.map(|speed| core::array::from_fn(|k| physical_flux[k] - speed * q[k]));
        let dt_area = 0.125;
        let v0 = 0.5;
        let v1 = v0 + dt_area * (w[1] - w[0]);
        let stage = predictor(q, v0, v1, dt_area, fluxes);
        let after = corrector(q, stage, v0, v1, dt_area, fluxes);
        for k in 0..3 {
            near(stage[k], q[k], q[k].abs());
            near(after[k], q[k], q[k].abs());
        }
    }

    #[test]
    fn extensive_rk2_closes_nonuniform_moving_cells_with_shared_boundary_fluxes() {
        let q = [[2.0, 6.0, 19.0], [3.0, -4.0, 29.0], [1.5, 2.0, 11.0]];
        let volumes = [0.3, 0.5, 0.2];
        let velocities = [-0.3, -0.1, 0.4, 0.7];
        let dt_area = 0.001;
        let faces_a = [
            [0.0, 10.0, -3.0],
            [2.0, 8.0, 9.0],
            [-1.0, 7.0, -2.0],
            [0.0, 5.0, 3.5],
        ];
        let faces_b = [
            [0.0, 11.0, -3.3],
            [1.7, 8.5, 8.0],
            [-0.8, 6.0, -1.0],
            [0.0, 5.4, 3.78],
        ];
        let mut change = [0.0; 3];
        for i in 0..3 {
            let next_volume = volumes[i] + dt_area * (velocities[i + 1] - velocities[i]);
            let stage = predictor(
                q[i],
                volumes[i],
                next_volume,
                dt_area,
                [faces_a[i], faces_a[i + 1]],
            );
            let next = corrector(
                q[i],
                stage,
                volumes[i],
                next_volume,
                dt_area,
                [faces_b[i], faces_b[i + 1]],
            );
            for k in 0..3 {
                change[k] += next_volume * next[k] - volumes[i] * q[i][k];
            }
        }
        for k in 0..3 {
            let expected =
                0.5 * dt_area * (faces_a[0][k] + faces_b[0][k] - faces_a[3][k] - faces_b[3][k]);
            near(change[k], expected, 30.0);
        }
    }

    #[test]
    fn pressure_work_rk2_refines_against_the_independent_adiabatic_volume_ode() {
        // Manufactured uniform chamber energy ODE dQ/dt=-p*dV/dt,
        // p=(gamma-1)*Q/V. This tests extensive time integration, not a
        // homogeneous approximation to the resolved moving-wall PDE.
        let gamma = 1.4;
        let finish = 0.4;
        let rate = -0.3;
        let initial_energy = 2.5;
        let final_volume: f64 = 1.0 + rate * finish;
        let expected = initial_energy * final_volume.powf(1.0 - gamma);
        let mut errors = Vec::new();
        for count in [10, 20, 40, 80] {
            let dt = finish / count as f64;
            let mut q = [1.0, 0.0, initial_energy];
            for i in 0..count {
                let v0 = 1.0 + rate * (i as f64 * dt);
                let v1 = 1.0 + rate * ((i + 1) as f64 * dt);
                let p0 = (gamma - 1.0) * q[2];
                let f0 = [[0.0, p0, 0.0], [0.0, p0, p0 * rate]];
                let stage = predictor(q, v0, v1, dt, f0);
                let p1 = (gamma - 1.0) * stage[2];
                q = corrector(q, stage, v0, v1, dt, [[0.0, p1, 0.0], [0.0, p1, p1 * rate]]);
            }
            errors.push((q[2] * final_volume - expected).abs());
        }
        for pair in errors.windows(2) {
            assert!(pair[0] / pair[1] > 3.9, "{errors:?}");
        }
    }
}
