//! G1/G3 local coupled acoustic modes on explicitly authored SI tensor data.
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::piezoelectric::tensor::{TensorPiezoelectric, waves::PiezoelectricWaveCondition};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::tensor::{mandel_basis, rotation};

fn material() -> TensorPiezoelectric {
    let mut claims = ClaimSet::new();
    for role in TensorPiezoelectric::requirements().unwrap() {
        let name = role.name();
        let suffix = name.rsplit('_').next().unwrap().as_bytes();
        let value = if name == "density" {
            7500.0
        } else {
            let (i, j) = ((suffix[0] - b'1') as usize, (suffix[1] - b'1') as usize);
            if name.starts_with("piezoelectric_stiffness") {
                if i == j {
                    [100e9, 80e9, 120e9, 60e9, 70e9, 50e9][i]
                } else {
                    0.0
                }
            } else if name.starts_with("piezoelectric_permittivity") {
                if i == j { [2e-8, 1.5e-8, 1e-8][i] } else { 0.0 }
            } else {
                match (i, j) {
                    (2, 2) => 15.0,
                    (0, 5) => 8.0,
                    _ => 0.0,
                }
            }
        };
        let dims = role.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored piezoacoustic tensor fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic piezo".into(),
            phase: "solid".into(),
            process: "authored stress-charge tensor".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    TensorPiezoelectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        ContentHash([31; 32]),
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-10 * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn apply(q: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    q.map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum())
}

#[test]
fn charge_free_wave_uses_gauss_constraint_and_does_not_clamp_transverse_displacement() {
    let material = material();
    let fixed = material
        .directional_waves(
            [0.0, 0.0, 1.0],
            PiezoelectricWaveCondition::FixedField,
            0.01,
        )
        .unwrap();
    let free = material
        .directional_waves(
            [0.0, 0.0, 1.0],
            PiezoelectricWaveCondition::ChargeFree,
            0.01,
        )
        .unwrap();
    for (i, expected) in [25e9, 35e9, 120e9].into_iter().enumerate() {
        near(fixed.mechanical.modes[i].directional_modulus_pa, expected);
        near(
            free.mechanical.modes[i].directional_modulus_pa,
            expected + if i == 2 { 225.0 / 1e-8 } else { 0.0 },
        );
        near(
            free.mechanical.modes[i].phase_speed_m_s,
            (free.mechanical.modes[i].directional_modulus_pa / 7500.0).sqrt(),
        );
    }
    assert!(
        material.open_circuit_stiffness_mandel_pa()[5][5] > material.stiffness_mandel_pa()[5][5]
    );
    // That local fixed-D shear stiffening MUST NOT enter the charge-free bulk wave.
    near(
        free.mechanical.modes[0].directional_modulus_pa,
        fixed.mechanical.modes[0].directional_modulus_pa,
    );
    near(free.field_per_displacement_gradient_v_m[2][2], -15.0 / 1e-8);
    assert_eq!(fixed.field_per_displacement_gradient_v_m, [[0.0; 3]; 3]);
}

#[test]
fn rotated_piezoacoustic_modes_preserve_speeds_and_satisfy_gauss_law() {
    let material = material();
    let q = rotation([0.0, 1.0, 0.0], 0.61);
    let rotated = material.rotated(ContentHash([32; 32]), q).unwrap();
    let n = [0.3, -0.4, 0.5];
    let before = material
        .directional_waves(n, PiezoelectricWaveCondition::ChargeFree, 0.01)
        .unwrap();
    let after = rotated
        .directional_waves(apply(q, n), PiezoelectricWaveCondition::ChargeFree, 0.01)
        .unwrap();
    for i in 0..3 {
        near(
            before.mechanical.modes[i].phase_speed_m_s,
            after.mechanical.modes[i].phase_speed_m_s,
        );
        let expected_velocity = apply(q, before.mechanical.modes[i].energy_velocity_m_s);
        for (actual, expected) in after.mechanical.modes[i]
            .energy_velocity_m_s
            .into_iter()
            .zip(expected_velocity)
        {
            near(actual, expected);
        }
        let a = before.mechanical.modes[i].polarization;
        let n = before.mechanical.direction;
        let basis = mandel_basis();
        let strain = core::array::from_fn(|k| {
            (0..3)
                .map(|i| (0..3).map(|j| basis[k][i][j] * n[j] * a[i]).sum::<f64>())
                .sum()
        });
        let response = material
            .field_response(strain, before.field_per_displacement_gradient_v_m[i])
            .unwrap();
        near(
            (0..3).map(|j| n[j] * response.displacement_c_m2[j]).sum(),
            0.0,
        );
        let field = before.field_per_displacement_gradient_v_m[i];
        // These fields are per UNIT strain, typically 10^9 V/m. Test the
        // angular constraint after scaling so cancellation in large products
        // cannot be mistaken for a nonzero transverse electric field.
        let field_scale = field
            .iter()
            .fold(0.0_f64, |scale, value| scale.max(value.abs()));
        let field = field.map(|value| value / field_scale.max(f64::MIN_POSITIVE));
        for j in 0..3 {
            near(field[j] * n[(j + 1) % 3] - field[(j + 1) % 3] * n[j], 0.0);
        }
    }
    assert!(
        material
            .directional_waves([0.0; 3], PiezoelectricWaveCondition::ChargeFree, 0.01)
            .is_err()
    );
    assert!(
        material
            .directional_waves(n, PiezoelectricWaveCondition::FixedField, 1.0)
            .is_err()
    );
}

#[test]
fn g0_charge_free_energy_velocity_includes_directional_electric_transport() {
    let material = material();
    let n: [f64; 3] = [0.6, 0.8, 0.0];
    let waves = material
        .directional_waves(n, PiezoelectricWaveCondition::ChargeFree, 0.01)
        .unwrap();
    let mode = waves.mechanical.modes[0];
    near(mode.polarization[2].abs(), 1.0);
    // For this decoupled z-polarized mode:
    // rho*omega(k)^2 = G_zx*kx^2 + G_zy*ky^2
    //                 + e_x,zx^2*kx^4/(2*(eps_x*kx^2+eps_y*ky^2)).
    // Differentiate this scalar dispersion exactly, not a frozen C_eff(n).
    let (gx, gy, e2, ex, ey, rho) = (25e9, 35e9, 64.0, 2e-8, 1.5e-8, 7500.0);
    let eps = ex * n[0].powi(2) + ey * n[1].powi(2);
    let modulus = gx * n[0].powi(2) + gy * n[1].powi(2) + 0.5 * e2 * n[0].powi(4) / eps;
    let speed = (modulus / rho).sqrt();
    near(mode.phase_speed_m_s, speed);
    let expected = [
        (gx * n[0] + e2 * n[0].powi(3) / eps - 0.5 * e2 * n[0].powi(5) * ex / eps.powi(2))
            / (rho * speed),
        (gy * n[1] - 0.5 * e2 * n[0].powi(4) * ey * n[1] / eps.powi(2)) / (rho * speed),
        0.0,
    ];
    for (actual, expected) in mode.energy_velocity_m_s.into_iter().zip(expected) {
        near(actual, expected);
    }
    let frozen_tensor_y_velocity = gy * n[1] / (rho * speed);
    assert!((mode.energy_velocity_m_s[1] - frozen_tensor_y_velocity).abs() > 10.0);
    near(
        (0..3).map(|i| n[i] * mode.energy_velocity_m_s[i]).sum(),
        speed,
    );
}

#[test]
fn g3_coupled_energy_velocities_match_the_full_wave_vector_dispersion_gradient() {
    let material = material();
    let k = [0.3, -0.4, 0.5];
    for condition in [
        PiezoelectricWaveCondition::FixedField,
        PiezoelectricWaveCondition::ChargeFree,
    ] {
        let waves = material.directional_waves(k, condition, 0.01).unwrap();
        for mode in 0..3 {
            let speed = waves.mechanical.modes[mode].phase_speed_m_s;
            for axis in 0..3 {
                let h = 1e-5;
                let mut plus = k;
                let mut minus = k;
                plus[axis] += h;
                minus[axis] -= h;
                let frequency = |k: [f64; 3]| {
                    k.iter().map(|v| v * v).sum::<f64>().sqrt()
                        * material
                            .directional_waves(k, condition, 0.01)
                            .unwrap()
                            .mechanical
                            .modes[mode]
                            .phase_speed_m_s
                };
                let derivative = (frequency(plus) - frequency(minus)) / (2.0 * h);
                assert!(
                    (waves.mechanical.modes[mode].energy_velocity_m_s[axis] - derivative).abs()
                        < 2e-8 * speed
                );
            }
            near(
                (0..3)
                    .map(|i| {
                        waves.mechanical.direction[i]
                            * waves.mechanical.modes[mode].energy_velocity_m_s[i]
                    })
                    .sum(),
                speed,
            );
        }
    }
}
