//! Actual augmented production pencil versus independent physical equations.

use super::*;

fn matrix(values: &[f64], n: usize) -> Csr {
    let mut a = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[i * n + j] != 0.0 {
                a.push(i, j, values[i * n + j]);
            }
        }
    }
    a.assemble()
}

fn scalar_wall(mass: f64, stiffness: f64, resistance: f64) -> Wall {
    Wall {
        dofs: vec![0],
        mass: matrix(&[1.0], 1),
        area_m2: 1.0,
        surface_density: mass,
        stiffness,
        resistance,
        arms: Vec::new(),
    }
}

fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.0),
        "{a:?} != {b:?}, error {:e}",
        (a - b).abs()
    );
}

#[test]
fn lossless_local_wall_pole_has_finite_coupled_motion_and_analytic_tangent() {
    // omega=1, sigma=K=1 is EXACTLY Z=0. The complete cavity/wall
    // equations instead give p=0, v=u, w=-i u without an impedance inverse.
    let fluid = [
        matrix(&[0.0], 1),
        matrix(&[1.0], 1),
        matrix(&[0.0], 1),
        matrix(&[0.0], 1),
    ];
    let pencil = CoupledPencil::prepare(
        1.0,
        fluid.each_ref(),
        &[0],
        vec![scalar_wall(1.0, 1.0, 0.0)],
        MechanicalHarmonicLimits::default(),
        &mut || false,
    )
    .unwrap();
    let u = C64::new(0.3, -0.2);
    let du = C64::new(-0.1, 0.05);
    let x = pencil.solve(&[u], &mut || false).unwrap();
    let dx = pencil
        .frequency_derivative(&x.pressure, &[u], &[du], &mut || false)
        .unwrap();
    let (wall, account) = pencil.wall_fields(&x.pressure, &mut || false).unwrap();
    let tangent = pencil
        .wall_derivatives(&x.pressure, &dx.pressure, &mut || false)
        .unwrap();
    near(x.pressure[0], C64::ZERO, 3e-15);
    near(wall[0].normal_velocity_m_s[0], u, 3e-15);
    near(wall[0].displacement_m[0], C64::new(0.0, -1.0) * u, 3e-15);
    near(dx.pressure[0], C64::new(0.0, 2.0) * u, 5e-15);
    near(
        tangent[0].normal_velocity_m_s_per_rad_s[0],
        u.scale(2.0) + du,
        5e-15,
    );
    near(
        tangent[0].displacement_m_per_rad_s[0],
        C64::new(0.0, -1.0) * (u + du),
        5e-15,
    );
    assert!((account.kinetic_energy_j - 0.25 * u.norm_sq()).abs() < 1e-15);
    assert!((account.elastic_energy_j - 0.25 * u.norm_sq()).abs() < 1e-15);
    assert_eq!(account.dissipation_w(), 0.0);
    assert!(account.fluid_to_wall_power_w.abs() < 1e-15);
    assert!(x.backward_error < 5e-13 && dx.backward_error < 5e-13);
    assert!(super::super::mechanical_admittance(1.0, 1.0, 1.0, 0.0, [], &mut || false).is_err());
}

#[test]
fn genuinely_singular_coupled_pole_is_not_regularized() {
    // At omega=1, sigma=4,K=3,C=A=1 the COUPLED resonance is exact,
    // whereas the local wall impedance is i and is perfectly invertible.
    let fluid = [
        matrix(&[0.0], 1),
        matrix(&[1.0], 1),
        matrix(&[0.0], 1),
        matrix(&[0.0], 1),
    ];
    assert!(matches!(
        CoupledPencil::prepare(
            1.0,
            fluid.each_ref(),
            &[0],
            vec![scalar_wall(4.0, 3.0, 0.0)],
            MechanicalHarmonicLimits::default(),
            &mut || false
        ),
        Err(PressureError::Factor(LuComplexError::Singular { .. }))
    ));
}

const W: [f64; 4] = [0.8, 0.2, 0.2, 0.5];
const K: [f64; 4] = [5.0, -1.0, -1.0, 3.0];
const M: [f64; 4] = [0.4, 0.1, 0.1, 0.3];
const G: [f64; 4] = [0.15, 0.02, 0.02, 0.12];
const ARMS: [(f64, f64); 2] = [(1.2, 0.15), (0.4, 1.0)]; // k,tau.

// Independent ORIGINAL-z formulation: p,w,z0,z1; v=i omega w.
// No production extension coordinates or assembled block matrices are used.
fn original_z_system(omega: f64, input: &[C64; 2], slope: &[C64; 2]) -> (Vec<C64>, Vec<C64>) {
    let n = 8;
    let mut a = vec![C64::ZERO; n * n];
    let mut da = vec![C64::ZERO; n * n];
    for i in 0..2 {
        for j in 0..2 {
            let ij = i * 2 + j;
            a[i * n + j] = C64::new(K[ij] - omega * omega * M[ij], omega * G[ij]);
            da[i * n + j] = C64::new(-2.0 * omega * M[ij], G[ij]);
            a[i * n + 2 + j] = C64::from_re(-omega * omega * W[ij]);
            da[i * n + 2 + j] = C64::from_re(-2.0 * omega * W[ij]);
            a[(2 + i) * n + j] = C64::from_re(-W[ij]);
            a[(2 + i) * n + 2 + j] = C64::new(
                (2.5 + 1.2 + 0.4 - 0.7 * omega * omega) * W[ij],
                omega * 0.03 * W[ij],
            );
            da[(2 + i) * n + 2 + j] = C64::new(-2.0 * 0.7 * omega * W[ij], 0.03 * W[ij]);
            for (arm, &(k, _)) in ARMS.iter().enumerate() {
                a[(2 + i) * n + 4 + 2 * arm + j] = C64::from_re(-det::sqrt(k) * W[ij]);
            }
        }
    }
    for (arm, &(k, tau)) in ARMS.iter().enumerate() {
        for i in 0..2 {
            let row = 4 + 2 * arm + i;
            a[row * n + 2 + i] = C64::from_re(-det::sqrt(k));
            a[row * n + row] = C64::new(1.0, omega * tau);
            da[row * n + row] = C64::new(0.0, tau);
        }
    }
    let factor = lu_complex_with_cancel(&a, n, &mut || false).unwrap();
    let mut x = vec![C64::ZERO; n];
    for i in 0..2 {
        x[i] = C64::new(0.0, omega) * input[i];
    }
    factor.solve_with_cancel(&mut x, &mut || false).unwrap();
    let mut dx = vec![C64::ZERO; n];
    for i in 0..2 {
        dx[i] = C64::new(0.0, 1.0) * input[i] + C64::new(0.0, omega) * slope[i];
    }
    for i in 0..n {
        for j in 0..n {
            dx[i] = dx[i] - da[i * n + j] * x[j];
        }
    }
    factor.solve_with_cancel(&mut dx, &mut || false).unwrap();
    (x, dx)
}

fn direct_form(a: &[C64], b: &[C64], matrix: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            s += matrix[i * 2 + j] * (a[i].re * b[j].re + a[i].im * b[j].im);
        }
    }
    s
}

#[test]
fn complete_original_maxwell_states_slopes_storage_and_reciprocal_work_agree() {
    let fluid = [
        matrix(&K, 2),
        matrix(&M, 2),
        matrix(&G, 2),
        matrix(&[0.0; 4], 2),
    ];
    let input = [C64::new(0.1, 0.2), C64::new(-0.06, 0.01)];
    let slope = [C64::new(-0.02, 0.01), C64::new(-0.001, -0.006)];
    for omega in [0.3, 1.1, 4.0, 10.0] {
        let wall = Wall {
            dofs: vec![0, 1],
            mass: matrix(&W, 2),
            area_m2: W.iter().sum(),
            surface_density: 0.7,
            stiffness: 2.5,
            resistance: 0.03,
            arms: ARMS
                .iter()
                .map(|&(k, tau)| (det::sqrt(k), 1.0 / tau))
                .collect(),
        };
        let pencil = CoupledPencil::prepare(
            omega,
            fluid.each_ref(),
            &[0, 1],
            vec![wall],
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
        let x = pencil.solve(&input, &mut || false).unwrap();
        let dx = pencil
            .frequency_derivative(&x.pressure, &input, &slope, &mut || false)
            .unwrap();
        let (walls, account) = pencil.wall_fields(&x.pressure, &mut || false).unwrap();
        let powers = pencil.wall_dissipation(&walls, &mut || false).unwrap();
        assert_eq!(powers.len(), 1);
        assert_eq!(powers[0].area_m2, W.iter().sum::<f64>());
        assert_eq!(powers[0].mass_kg, 0.7 * W.iter().sum::<f64>());
        assert_eq!(
            powers[0].resistance_dissipation_w,
            account.resistance_dissipation_w
        );
        assert_eq!(
            powers[0].relaxation_dissipation_w,
            account.relaxation_dissipation_w
        );
        let derivatives = pencil
            .wall_derivatives(&x.pressure, &dx.pressure, &mut || false)
            .unwrap();
        let (exact, tangent) = original_z_system(omega, &input, &slope);
        let mut velocity = [C64::ZERO; 2];
        for i in 0..2 {
            near(x.pressure[i], exact[i], 2e-12);
            near(dx.pressure[i], tangent[i], 3e-11);
            near(walls[0].displacement_m[i], exact[2 + i], 2e-12);
            velocity[i] = C64::new(0.0, omega) * exact[2 + i];
            near(walls[0].normal_velocity_m_s[i], velocity[i], 2e-12);
            near(
                derivatives[0].displacement_m_per_rad_s[i],
                tangent[2 + i],
                3e-11,
            );
            near(
                derivatives[0].normal_velocity_m_s_per_rad_s[i],
                C64::new(0.0, 1.0) * exact[2 + i] + C64::new(0.0, omega) * tangent[2 + i],
                3e-11,
            );
            for arm in 0..2 {
                assert_eq!(
                    walls[0].relaxation_extension[arm][i],
                    x.pressure[6 + 2 * arm + i]
                );
                assert_eq!(
                    derivatives[0].relaxation_extension_per_rad_s[arm][i],
                    dx.pressure[6 + 2 * arm + i]
                );
                near(
                    walls[0].relaxation_memory[arm][i],
                    exact[4 + 2 * arm + i],
                    2e-12,
                );
                near(
                    derivatives[0].relaxation_memory_per_rad_s[arm][i],
                    tangent[4 + 2 * arm + i],
                    3e-11,
                );
            }
        }
        let vv = direct_form(&velocity, &velocity, &W);
        let ww = direct_form(&exact[2..4], &exact[2..4], &W);
        assert!((account.kinetic_energy_j - 0.25 * 0.7 * vv).abs() < 2e-12);
        assert!((account.elastic_energy_j - 0.25 * 2.5 * ww).abs() < 2e-12);
        assert!((account.resistance_dissipation_w - 0.5 * 0.03 * vv).abs() < 2e-12);
        let (mut memory_storage, mut memory_heat) = (0.0, 0.0);
        for (arm, &(k, tau)) in ARMS.iter().enumerate() {
            let e: Vec<_> = (0..2)
                .map(|i| exact[2 + i].scale(det::sqrt(k)) - exact[4 + 2 * arm + i])
                .collect();
            let ee = direct_form(&e, &e, &W);
            memory_storage += 0.25 * ee;
            memory_heat += 0.5 * ee / tau;
        }
        assert!(account.relaxation_energy_j > 0.0);
        assert!((account.relaxation_energy_j - memory_storage).abs() < 2e-12);
        assert!((account.relaxation_dissipation_w - memory_heat).abs() < 2e-12);
        let power = 0.5 * direct_form(&exact[..2], &velocity, &W);
        let supplied = 0.5 * (exact[0].conj() * input[0] + exact[1].conj() * input[1]).re;
        let robin = 0.5 * direct_form(&exact[..2], &exact[..2], &G);
        assert!((account.fluid_to_wall_power_w - power).abs() < 2e-12);
        assert!(account.power_balance_residual_w().abs() < 2e-12);
        assert!((supplied - robin - account.dissipation_w()).abs() < 3e-12);
        assert!(x.backward_error < 5e-13 && dx.backward_error < 5e-13);
    }

    // A nearly frozen arm has a representable z much smaller than c*w and e.
    // Check against the original memory equation, without subtracting those
    // nearly equal phasors; its derivative must retain the explicit omega term.
    let fluid = [
        matrix(&[0.0], 1),
        matrix(&[1.0], 1),
        matrix(&[0.0], 1),
        matrix(&[0.0], 1),
    ];
    let mut wall = scalar_wall(1.0, 1.0, 1.0);
    wall.arms.push((1.0, 1e-20));
    let pencil = CoupledPencil::prepare(
        1.0,
        fluid.each_ref(),
        &[0],
        vec![wall],
        MechanicalHarmonicLimits::default(),
        &mut || false,
    )
    .unwrap();
    let input = [C64::new(0.2, -0.1)];
    let x = pencil.solve(&input, &mut || false).unwrap();
    let dx = pencil
        .frequency_derivative(&x.pressure, &input, &[C64::ZERO], &mut || false)
        .unwrap();
    let (fields, _) = pencil.wall_fields(&x.pressure, &mut || false).unwrap();
    let derivatives = pencil
        .wall_derivatives(&x.pressure, &dx.pressure, &mut || false)
        .unwrap();
    let denominator = C64::new(1.0, 1e20);
    let expected = fields[0].displacement_m[0] / denominator;
    let expected_slope = derivatives[0].displacement_m_per_rad_s[0] / denominator
        - fields[0].displacement_m[0] * C64::new(0.0, 1e20) / (denominator * denominator);
    assert!((fields[0].relaxation_memory[0][0] - expected).abs() < 5e-13 * expected.abs());
    assert!(
        (derivatives[0].relaxation_memory_per_rad_s[0][0] - expected_slope).abs()
            < 5e-13 * expected_slope.abs()
    );
}

#[test]
fn exact_pressure_release_retains_all_wall_coordinates_and_complete_resource_limits() {
    let fluid = [
        matrix(&[4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 3.0], 3),
        matrix(&[0.3, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.4], 3),
        matrix(&[0.0; 9], 3),
        matrix(&[0.0; 9], 3),
    ];
    let wall = || Wall {
        dofs: vec![0, 1],
        mass: matrix(&W, 2),
        area_m2: W.iter().sum(),
        surface_density: 0.7,
        stiffness: 2.5,
        resistance: 0.1,
        arms: vec![(1.0, 2.0)],
    };
    let limits = MechanicalHarmonicLimits {
        max_states: 8,
        max_dense_bytes: 32 * 8 * 8,
        residual_tolerance: 5e-13,
    };
    let pencil = CoupledPencil::prepare(
        1.1,
        fluid.each_ref(),
        &[0, 2],
        vec![wall()],
        limits,
        &mut || false,
    )
    .unwrap();
    assert_eq!(pencil.wall_dofs(0), Some([0, 1].as_slice()));
    for bad in [
        MechanicalHarmonicLimits {
            max_states: 7,
            ..limits
        },
        MechanicalHarmonicLimits {
            max_dense_bytes: 32 * 8 * 8 - 1,
            ..limits
        },
    ] {
        assert!(matches!(
            CoupledPencil::prepare(
                1.1,
                fluid.each_ref(),
                &[0, 2],
                vec![wall()],
                bad,
                &mut || false
            ),
            Err(PressureError::Limit { .. })
        ));
    }
    let input = [C64::new(0.1, 0.2), C64::ZERO, C64::new(-0.02, 0.01)];
    let primary = pencil.solve(&input, &mut || false).unwrap();
    let tangent = pencil
        .frequency_derivative(&primary.pressure, &input, &[C64::ZERO; 3], &mut || false)
        .unwrap();
    assert_eq!(primary.pressure[1], C64::ZERO);
    assert_eq!(tangent.pressure[1], C64::ZERO);
    let (walls, _) = pencil
        .wall_fields(&primary.pressure, &mut || false)
        .unwrap();
    assert_eq!(walls[0].displacement_m.len(), 2);
    assert_eq!(walls[0].normal_velocity_m_s.len(), 2);
    assert_eq!(walls[0].relaxation_memory[0].len(), 2);
    let mut calls = 0;
    assert!(matches!(
        pencil.frequency_derivative(&primary.pressure, &input, &[C64::ZERO; 3], &mut || {
            calls += 1;
            calls == 7
        }),
        Err(PressureError::Cancelled)
    ));
    let replay = pencil
        .frequency_derivative(&primary.pressure, &input, &[C64::ZERO; 3], &mut || false)
        .unwrap();
    assert_eq!(tangent.pressure, replay.pressure);
    assert_eq!(
        tangent.backward_error.to_bits(),
        replay.backward_error.to_bits()
    );
}

#[test]
fn resolved_fast_relaxation_heat_survives_rounded_displacement_memory_cancellation() {
    let fluid = [
        matrix(&[0.0], 1),
        matrix(&[1.0], 1),
        matrix(&[0.0], 1),
        matrix(&[0.0], 1),
    ];
    let mut wall = scalar_wall(1.0, 1.0, 0.0);
    let coupling = 1e20;
    let rate = 1e20;
    wall.arms.push((coupling, rate));
    let pencil = CoupledPencil::prepare(
        1.0,
        fluid.each_ref(),
        &[0],
        vec![wall],
        MechanicalHarmonicLimits::default(),
        &mut || false,
    )
    .unwrap();
    let solution = pencil.solve(&[C64::new(0.3, -0.2)], &mut || false).unwrap();
    let (fields, account) = pencil
        .wall_fields(&solution.pressure, &mut || false)
        .unwrap();
    let extension = fields[0].relaxation_extension[0][0];
    assert_eq!(extension, solution.pressure[3]);
    assert!(extension.abs() > 0.0);
    // The physical extension is twenty orders smaller than c*w and z, so
    // subtracting those rounded OUTPUT phasors discards its nonzero heat.
    let recovered = fields[0].displacement_m[0].scale(coupling) - fields[0].relaxation_memory[0][0];
    assert_eq!(recovered, C64::ZERO);
    let power = pencil.wall_dissipation(&fields, &mut || false).unwrap()[0];
    let expected = 0.5 * rate * extension.norm_sq();
    assert!(expected > 0.0);
    assert!((power.relaxation_dissipation_w - expected).abs() < 2e-12 * expected);
    assert_eq!(
        power.relaxation_dissipation_w,
        account.relaxation_dissipation_w
    );
    assert_eq!(power.resistance_dissipation_w, 0.0);
    let mut broken = fields.clone();
    broken[0].relaxation_extension.clear();
    assert!(matches!(
        pencil.wall_dissipation(&broken, &mut || false),
        Err(PressureError::Invalid(_))
    ));
    broken = fields.clone();
    broken[0].relaxation_extension[0][0].re = f64::NAN;
    assert!(matches!(
        pencil.wall_dissipation(&broken, &mut || false),
        Err(PressureError::Invalid(_))
    ));
    assert!(matches!(
        pencil.wall_dissipation(&fields, &mut || true),
        Err(PressureError::Cancelled)
    ));
    assert_eq!(
        pencil.wall_dissipation(&fields, &mut || false).unwrap()[0],
        power
    );
}
