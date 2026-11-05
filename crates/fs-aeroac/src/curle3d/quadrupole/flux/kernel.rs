//! Instantaneous inviscid source arithmetic, independent of the flow integrator.

use super::super::kernel::reynolds_moment;
use fs_math::dd::Dd;

pub(super) fn rotate(rotation: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    rotation.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

pub(super) fn proper_rotation(r: [[f64; 3]; 3]) -> bool {
    if r.iter().flatten().any(|x| !x.is_finite()) {
        return false;
    }
    let tolerance = 512.0 * f64::EPSILON;
    for i in 0..3 {
        for j in 0..3 {
            let dot = r[i][0] * r[j][0] + r[i][1] * r[j][1] + r[i][2] * r[j][2];
            if (dot - if i == j { 1.0 } else { 0.0 }).abs() > tolerance {
                return false;
            }
        }
    }
    let determinant = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        + r[0][1] * (r[1][2] * r[2][0] - r[1][0] * r[2][2])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    (determinant - 1.0).abs() <= 3.0 * tolerance
}

/// Actual rho ui uj and pressure-density difference; molecular/SGS stresses absent.
pub(super) fn inviscid_moment(
    reference_density: f64,
    sound_speed: f64,
    reference_pressure: f64,
    density: f64,
    pressure: f64,
    volume: f64,
    velocity: [f64; 3],
) -> Option<[f64; 6]> {
    if [
        reference_density,
        sound_speed,
        reference_pressure,
        density,
        pressure,
        volume,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return None;
    }
    let mut moment = reynolds_moment(density, volume, velocity)?;
    // Differences precede products; DD retains cancellation in a near-isentropic
    // state without subtracting two large absolute pressure-like fields first.
    let density_difference = Dd::from_f64(density) - Dd::from_f64(reference_density);
    let pressure_difference = Dd::from_f64(pressure) - Dd::from_f64(reference_pressure);
    let c = Dd::from_f64(sound_speed);
    let isotropic = (pressure_difference - (density_difference * c) * c) * Dd::from_f64(volume);
    for value in &mut moment[..3] {
        *value = (Dd::from_f64(*value) + isotropic).to_f64();
    }
    moment.iter().all(|x| x.is_finite()).then_some(moment)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 4e-14 * (1.0 + expected.abs()),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g0_isentropic_difference_reduces_to_physical_reynolds_moment() {
        let q = inviscid_moment(1.0, 4.0, 100.0, 2.0, 116.0, 0.5, [3.0, -2.0, 1.0]).unwrap();
        assert_eq!(q, [9.0, 4.0, 1.0, -6.0, 3.0, -2.0]);
        // Ambient pressure can be large without contaminating an exactly
        // representable isentropic difference by a spurious constant trace.
        let q = inviscid_moment(1.0, 4.0, 1e12, 1.0078125, 1e12 + 0.125, 2.0, [0.0; 3]).unwrap();
        assert_eq!(q, [0.0; 6]);
    }

    #[test]
    fn g0_entropy_pressure_contribution_keeps_its_trace_and_sign() {
        let q = inviscid_moment(1.0, 4.0, 100.0, 2.0, 120.0, 0.5, [3.0, -2.0, 1.0]).unwrap();
        assert_eq!(q, [11.0, 6.0, 3.0, -6.0, 3.0, -2.0]);
        let negative = inviscid_moment(1.0, 4.0, 100.0, 2.0, 108.0, 0.5, [0.0; 3]).unwrap();
        assert_eq!(negative, [-4.0, -4.0, -4.0, 0.0, 0.0, 0.0]);
        assert_eq!(
            inviscid_moment(1.0, 4.0, 100.0, 1.0, 100.0, 1.0, [0.0; 3]).unwrap(),
            [0.0; 6]
        );
    }

    #[test]
    fn g3_full_stress_rotates_as_a_second_order_tensor() {
        let r = [[0.48, -0.6, 0.64], [0.36, 0.8, 0.48], [-0.8, 0.0, 0.6]];
        assert!(proper_rotation(r));
        let u = [2.0, -3.0, 0.5];
        let q = inviscid_moment(1.2, 343.0, 101325.0, 1.3, 112000.0, 0.002, u).unwrap();
        let full = [[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]];
        let rotated =
            inviscid_moment(1.2, 343.0, 101325.0, 1.3, 112000.0, 0.002, rotate(r, u)).unwrap();
        for (entry, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .enumerate()
        {
            let mut expected = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    expected += r[i][a] * full[a][b] * r[j][b];
                }
            }
            close(rotated[entry], expected);
        }
    }

    #[test]
    fn g3_cell_volume_and_consistent_density_pressure_rescaling() {
        let q = inviscid_moment(1.0, 3.0, 12.0, 1.2, 13.0, 0.1, [0.5, -0.25, 1.0]).unwrap();
        let volume = inviscid_moment(1.0, 3.0, 12.0, 1.2, 13.0, 0.4, [0.5, -0.25, 1.0]).unwrap();
        let density = inviscid_moment(2.0, 3.0, 24.0, 2.4, 26.0, 0.1, [0.5, -0.25, 1.0]).unwrap();
        for i in 0..6 {
            close(volume[i], 4.0 * q[i]);
            close(density[i], 2.0 * q[i]);
        }
    }

    #[test]
    fn g0_invalid_rotations_domains_and_unrepresentable_stress_refuse() {
        assert!(!proper_rotation([
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]));
        assert!(!proper_rotation([
            [1.0, 0.1, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]));
        assert!(!proper_rotation([[f64::NAN; 3]; 3]));
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(inviscid_moment(1.0, 3.0, 12.0, bad, 13.0, 0.1, [0.0; 3]).is_none());
            assert!(inviscid_moment(1.0, 3.0, bad, 1.2, 13.0, 0.1, [0.0; 3]).is_none());
        }
        assert!(
            inviscid_moment(1.0, 3.0, 12.0, 1.2, 13.0, 0.1, [f64::INFINITY, 0.0, 0.0]).is_none()
        );
        assert!(inviscid_moment(1.0, 3.0, 12.0, 1e308, 13.0, 1e308, [0.0; 3]).is_none());
    }
}
