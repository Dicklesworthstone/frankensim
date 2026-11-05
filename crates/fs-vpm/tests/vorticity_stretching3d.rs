//! Local 3-D vortex stretching and enstrophy production from analytic gradients.

use fs_vpm::filament3d::WakeKinematics;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn rigid_rotation_has_vorticity_but_no_stretching_or_enstrophy_production() {
    let rate = 3.25;
    let sample = WakeKinematics {
        velocity: [0.0; 3],
        gradient: [[0.0, -rate, 0.0], [rate, 0.0, 0.0], [0.0, 0.0, 0.0]],
    };
    assert_eq!(sample.vorticity(), [0.0, 0.0, 2.0 * rate]);
    assert_eq!(sample.strain_rate(), [[0.0; 3]; 3]);
    assert_eq!(sample.vorticity_stretching(), [0.0; 3]);
    assert_eq!(sample.enstrophy_production_rate(), 0.0);
}

#[test]
fn incompressible_strain_stretches_and_tilts_the_actual_curl() {
    let sample = WakeKinematics {
        velocity: [0.0; 3],
        gradient: [[1.0, -2.0, 0.0], [2.0, -1.0, 3.0], [0.0, -3.0, 0.0]],
    };
    assert_eq!(sample.divergence(), 0.0);
    let omega = sample.vorticity();
    assert_eq!(omega, [-6.0, 0.0, 4.0]);
    assert_eq!(sample.vorticity_stretching(), [-6.0, 0.0, 0.0]);
    assert_eq!(sample.enstrophy_production_rate(), 36.0);
    let strain = sample.strain_rate();
    for i in 0..3 { for j in 0..3 { assert_eq!(strain[i][j], strain[j][i]); } }
    let s_omega: [f64; 3] = std::array::from_fn(|i| {
        strain[i][0] * omega[0] + strain[i][1] * omega[1] + strain[i][2] * omega[2]
    });
    let production = omega[0] * s_omega[0] + omega[1] * s_omega[1] + omega[2] * s_omega[2];
    assert_eq!(production, sample.enstrophy_production_rate());
}

#[test]
fn gradient_scaling_gives_physical_stretching_and_enstrophy_units() {
    let base = WakeKinematics {
        velocity: [0.0; 3],
        gradient: [[0.7, -1.1, 0.4], [1.6, -0.2, -0.8], [0.3, 0.9, -0.5]],
    };
    let base_stretch = base.vorticity_stretching();
    let base_prod = base.enstrophy_production_rate();
    for rate_scale in [1e-8, 1e-3, 1.0, 1e3, 1e8] {
        let scaled = WakeKinematics {
            velocity: [0.0; 3],
            gradient: base.gradient.map(|row| row.map(|value| value * rate_scale)),
        };
        let stretch = scaled.vorticity_stretching();
        for i in 0..3 {
            near(stretch[i], base_stretch[i] * rate_scale * rate_scale,
                2e-12 * (1.0 + base_stretch[i].abs() * rate_scale * rate_scale));
        }
        near(scaled.enstrophy_production_rate(), base_prod * rate_scale * rate_scale * rate_scale,
            5e-12 * (1.0 + base_prod.abs() * rate_scale * rate_scale * rate_scale));
    }
}
