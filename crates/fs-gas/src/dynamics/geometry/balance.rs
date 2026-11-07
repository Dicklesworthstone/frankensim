//! Shared stationary quasi-1D momentum and metric arithmetic.
//!
//! Pressure traction and the side-wall source are cancelled before area
//! multiplication. The stationary side wall changes momentum but performs no
//! energy work. Constitutive admission stays with the original flow owner.

/// One forward finite-volume stage with a cell pressure source p dA.
pub(in crate::dynamics) fn pressure_balanced_update(
    state: [f64; 3],
    dt_over_volume: f64,
    areas: [f64; 2],
    fluxes: [[f64; 3]; 2],
    pressure: f64,
) -> [f64; 3] {
    core::array::from_fn(|k| {
        let traction = if k == 1 { pressure } else { 0.0 };
        state[k]
            + dt_over_volume
                * (areas[0] * (fluxes[0][k] - traction) - areas[1] * (fluxes[1][k] - traction))
    })
}

/// Face-weighted wave crossing rate for admitted positive geometry/speeds.
/// The owning mesh checks the result before it controls an accepted step.
#[cfg(any(feature = "material", test))]
pub(super) fn face_signal_rate(areas: [f64; 2], volume: f64, speeds: [f64; 2]) -> f64 {
    let left = areas[0] / volume;
    let right = areas[1] / volume;
    if left == right && speeds[0] == speeds[1] {
        speeds[0] * left
    } else {
        (0.5 * left) * speeds[0] + (0.5 * right) * speeds[1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64, scale: f64) {
        assert!(a.is_finite() && b.is_finite());
        assert!(
            (a - b).abs() <= 64.0 * f64::EPSILON * scale.max(1.0),
            "{a:.16e} != {b:.16e}"
        );
    }

    #[test]
    fn resting_pressure_cancels_before_large_unequal_area_products() {
        let q = [10.0, 0.0, -713.0];
        let p = 1_405_488.503_926_915;
        for areas in [[0.003, 0.09], [2e299, 8e299], [8e299, 2e299]] {
            let next = pressure_balanced_update(q, 0.125, areas, [[0.0, p, 0.0]; 2], p);
            assert_eq!(next, q);
        }
    }

    #[test]
    fn uniform_moving_state_uses_geometric_advection_without_pressure_acceleration() {
        let rho = 10.0;
        let u = -23.0;
        let p = 1_405_488.503_926_915;
        let q = [rho, rho * u, 0.5 * rho * u * u - 900.0];
        let flux = [rho * u, rho * u * u + p, u * (q[2] + p)];
        let areas = [0.01, 0.037];
        let dt_over_volume = 1e-7;
        let got = pressure_balanced_update(q, dt_over_volume, areas, [flux; 2], p);
        let da = areas[1] - areas[0];
        let expected = [
            rho - dt_over_volume * da * rho * u,
            rho * u - dt_over_volume * da * rho * u * u,
            q[2] - dt_over_volume * da * u * (q[2] + p),
        ];
        for k in 0..3 {
            near(got[k], expected[k], q[k].abs());
        }
    }

    #[test]
    fn integrated_momentum_includes_side_pressure_and_mass_energy_only_boundary_fluxes() {
        let areas = [0.03, 0.05, 0.08, 0.04, 0.09];
        let volumes = [0.003, 0.012, 0.018, 0.007];
        let pressures = [1.4e6, 1.2e6, 1.7e6, 1.1e6];
        let fluxes = [
            [3.0, 1.5e6, 8e6],
            [7.0, 1.3e6, 6e6],
            [-2.0, 1.8e6, -2e6],
            [5.0, 1.6e6, 4e6],
            [4.0, 1.2e6, 1e6],
        ];
        let q = [10.0, 40.0, 5e4];
        let dt = 3e-7;
        let mut change = [0.0; 3];
        let mut side = 0.0;
        for i in 0..4 {
            let next = pressure_balanced_update(
                q,
                dt / volumes[i],
                [areas[i], areas[i + 1]],
                [fluxes[i], fluxes[i + 1]],
                pressures[i],
            );
            for k in 0..3 {
                change[k] += (next[k] - q[k]) * volumes[i];
            }
            side += pressures[i] * (areas[i + 1] - areas[i]);
        }
        for k in 0..3 {
            let expected = dt
                * (areas[0] * fluxes[0][k] - areas[4] * fluxes[4][k]
                    + if k == 1 { side } else { 0.0 });
            near(
                change[k],
                expected,
                q[k].abs() * volumes.iter().sum::<f64>(),
            );
        }
    }

    #[test]
    fn pressure_source_quadrature_refines_at_second_order_for_smooth_nozzle() {
        // A=1+x^2, p=2+x^2. The independent integrated pressure force is
        // integral_0^1 p A' dx = 5/2. The runtime source samples each cell.
        let mut errors = Vec::new();
        for n in [8, 16, 32, 64] {
            let mut momentum = 0.0;
            for i in 0..n {
                let l = i as f64 / n as f64;
                let r = (i + 1) as f64 / n as f64;
                let c = 0.5 * (l + r);
                let volume = (r - l) + (r * r * r - l * l * l) / 3.0;
                let result = pressure_balanced_update(
                    [0.0; 3],
                    1.0 / volume,
                    [1.0 + l * l, 1.0 + r * r],
                    [[0.0; 3]; 2],
                    2.0 + c * c,
                );
                momentum += result[1] * volume;
            }
            errors.push((momentum - 2.5).abs());
        }
        for w in errors.windows(2) {
            assert!(w[0] / w[1] > 3.99 && w[0] / w[1] < 4.01, "{errors:?}");
        }
    }

    #[test]
    fn face_metric_uses_both_areas_and_survives_extensive_unit_scaling() {
        near(
            face_signal_rate([2.0, 8.0], 3.0, [5.0, 11.0]),
            49.0 / 3.0,
            17.0,
        );
        for scale in [1e-280, 1.0, 1e280] {
            near(
                face_signal_rate([2.0 * scale, 8.0 * scale], 3.0 * scale, [5e10, 11e10]),
                49e10 / 3.0,
                17e10,
            );
        }
        assert_eq!(
            face_signal_rate([0.07; 2], 0.014, [400.0; 2]),
            400.0 * (0.07 / 0.014)
        );
    }
}
