use super::{Viscosity, no_slip_traction, traction};

#[derive(Clone, Copy, PartialEq)]
enum End {
    Periodic,
    Slip,
    Zero,
    NoSlip,
}

fn native<const D: usize>(
    axis: usize,
    distance: f64,
    velocity: [[f64; D]; 2],
    face: Viscosity,
    gradient: impl Fn(usize, usize, usize) -> f64,
) -> [f64; D] {
    std::array::from_fn(|a| {
        let normal = (velocity[1][a] - velocity[0][a]) / distance;
        if a == axis {
            let divergence: f64 = (0..D)
                .filter(|&j| j != axis)
                .map(|j| 0.5 * gradient(0, j, j) + 0.5 * gradient(1, j, j))
                .sum();
            ((4.0 / 3.0) * face.shear + face.bulk) * normal
                + (face.bulk - (2.0 / 3.0) * face.shear) * divergence
        } else {
            face.shear * (normal + 0.5 * gradient(0, axis, a) + 0.5 * gradient(1, axis, a))
        }
    })
}

/// Independent finite-volume assembly and cell quadratic form. The production
/// stress is compared to the original direct-normal/centered-cross stress.
fn balance<const D: usize>(
    shape: [usize; D],
    spacing: [f64; D],
    ends: [[End; 2]; D],
    velocity: &[[f64; D]],
    material: &[Viscosity],
    corrected: bool,
) -> (f64, f64, Vec<[f64; D]>) {
    let (production, dissipation, force, _) = balance_with_walls(
        shape,
        spacing,
        ends,
        [[[0.0; D]; 2]; D],
        velocity,
        material,
        corrected,
    );
    (production, dissipation, force)
}

fn balance_with_walls<const D: usize>(
    shape: [usize; D],
    spacing: [f64; D],
    ends: [[End; 2]; D],
    wall_velocity: [[[f64; D]; 2]; D],
    velocity: &[[f64; D]],
    material: &[Viscosity],
    corrected: bool,
) -> (f64, f64, Vec<[f64; D]>, f64) {
    let n = velocity.len();
    let mut c = vec![[[0.0; D]; D]; n];
    let mut h = c.clone();
    let mut force = vec![[0.0; D]; n];
    for i in 0..n {
        for axis in 0..D {
            let stride: usize = shape[..axis].iter().product();
            let position = (i / stride) % shape[axis];
            let mut jump = [[0.0; D]; 2];
            for (side, sign) in [(0, -1.0), (1, 1.0)] {
                let outside = position == if side == 0 { 0 } else { shape[axis] - 1 };
                let other = if !outside || ends[axis][side] == End::Periodic {
                    let index = if outside {
                        if side == 0 {
                            i + (shape[axis] - 1) * stride
                        } else {
                            i - (shape[axis] - 1) * stride
                        }
                    } else if side == 0 {
                        i - stride
                    } else {
                        i + stride
                    };
                    velocity[index]
                } else {
                    let mut ghost = velocity[i];
                    if ends[axis][side] == End::Slip {
                        ghost[axis] = -ghost[axis];
                    } else if ends[axis][side] == End::NoSlip {
                        ghost = std::array::from_fn(|a| {
                            2.0 * wall_velocity[axis][side][a] - velocity[i][a]
                        });
                    }
                    ghost
                };
                for a in 0..D {
                    jump[side][a] = sign * (other[a] - velocity[i][a]) / spacing[axis];
                }
            }
            for a in 0..D {
                c[i][a][axis] = 0.5 * (jump[0][a] + jump[1][a]);
                h[i][a][axis] = 0.5 * (jump[1][a] - jump[0][a]);
            }
        }
    }
    let mut wall_power = 0.0;
    for i in 0..n {
        for axis in 0..D {
            let stride: usize = shape[..axis].iter().product();
            let position = (i / stride) % shape[axis];
            if position + 1 < shape[axis] || ends[axis][1] == End::Periodic {
                let k = if position + 1 < shape[axis] {
                    i + stride
                } else {
                    i - (shape[axis] - 1) * stride
                };
                let cells = [material[i], material[k]];
                let face = Viscosity {
                    shear: 0.5 * cells[0].shear + 0.5 * cells[1].shear,
                    bulk: 0.5 * cells[0].bulk + 0.5 * cells[1].bulk,
                };
                let read = |side: usize, a: usize, j: usize| c[[i, k][side]][a][j];
                let pair = [velocity[i], velocity[k]];
                let stress = if corrected {
                    traction(axis, spacing[axis], pair, cells, read)
                } else {
                    native(axis, spacing[axis], pair, face, read)
                };
                for a in 0..D {
                    force[i][a] += stress[a] / spacing[axis];
                    force[k][a] -= stress[a] / spacing[axis];
                }
            }
            for (side, sign) in [(0, 1.0), (1, -1.0)] {
                let outside = position == if side == 0 { 0 } else { shape[axis] - 1 };
                if outside && ends[axis][side] == End::Slip {
                    let mu = material[i].shear;
                    let bulk = material[i].bulk;
                    let normal = sign * 2.0 * velocity[i][axis] / spacing[axis];
                    let transverse: f64 = (0..D).filter(|&j| j != axis).map(|j| c[i][j][j]).sum();
                    let traction =
                        ((4.0 / 3.0) * mu + bulk) * normal + (bulk - (2.0 / 3.0) * mu) * transverse;
                    force[i][axis] -= sign * traction / spacing[axis];
                } else if outside && ends[axis][side] == End::NoSlip {
                    let wall = wall_velocity[axis][side];
                    let stress = no_slip_traction(
                        axis,
                        -sign,
                        0.5 * spacing[axis],
                        velocity[i],
                        wall,
                        material[i],
                        |a, j| c[i][a][j],
                    );
                    for a in 0..D {
                        force[i][a] -= sign * stress[a] / spacing[axis];
                        wall_power -= sign * wall[a] * stress[a] / spacing[axis];
                    }
                }
            }
        }
    }
    let mut production = 0.0;
    let mut dissipation = 0.0;
    for i in 0..n {
        let mu = material[i].shear;
        let bulk = material[i].bulk;
        let mut divergence = 0.0;
        let mut sym2 = 0.0;
        let mut h2 = 0.0;
        let mut diagonal_h2 = 0.0;
        for a in 0..D {
            production += velocity[i][a] * force[i][a];
            divergence += c[i][a][a];
            diagonal_h2 += h[i][a][a].powi(2);
            for j in 0..D {
                sym2 += (0.5 * (c[i][a][j] + c[i][j][a])).powi(2);
                h2 += h[i][a][j].powi(2);
            }
        }
        dissipation += 2.0 * mu * sym2
            + (bulk - (2.0 / 3.0) * mu) * divergence.powi(2)
            + mu * h2
            + (bulk + mu / 3.0) * diagonal_h2;
    }
    (production, dissipation, force, wall_power)
}

fn supported_ends<const D: usize>(shape: [usize; D], spacing: [f64; D]) {
    let n = shape.iter().product();
    let velocity: Vec<[f64; D]> = (0..n)
        .map(|i| std::array::from_fn(|a| 2.0 * ((13 * i + 7 * a + 1) as f64).sin()))
        .collect();
    let material: Vec<_> = (0..n)
        .map(|i| Viscosity {
            shear: 10.0_f64.powf(4.0 * ((3 * i + 1) as f64).sin()),
            bulk: if i % 3 == 0 {
                0.0
            } else {
                10.0_f64.powf(2.0 * ((5 * i + 4) as f64).sin())
            },
        })
        .collect();
    let choices = [
        [End::Periodic; 2],
        [End::Slip; 2],
        [End::Zero; 2],
        [End::Slip, End::Zero],
        [End::Zero, End::Slip],
    ];
    for mut encoded in 0..5_usize.pow(D as u32) {
        let ends = std::array::from_fn(|_| {
            let pair = choices[encoded % 5];
            encoded /= 5;
            pair
        });
        let (production, dissipation, _) =
            balance(shape, spacing, ends, &velocity, &material, true);
        assert!(dissipation >= 0.0 && production <= 0.0);
        assert!(
            (production + dissipation).abs() <= 3e-14 * dissipation.max(1.0),
            "energy rate {production}, independent quadratic form {dissipation}"
        );
    }
}

#[test]
fn arbitrary_positive_materials_dissipate_with_periodic_slip_and_zero_flux_ends() {
    supported_ends([4, 5], [0.7, 1.3]);
    supported_ends([4, 3, 2], [0.7, 1.3, 0.9]);
    supported_ends([4, 4, 1], [0.7, 1.3, 0.9]);
}

#[test]
fn concentrated_viscosity_no_longer_amplifies_the_periodic_four_by_four_mode() {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    let velocity: Vec<_> = (0..16)
        .map(|i| [-a[i % 4] * b[i / 4], a[i / 4] * b[i % 4]])
        .collect();
    let material: Vec<_> = (0..16)
        .map(|i| Viscosity {
            shear: if i == 0 { 300.0 } else { 1.0 },
            bulk: 0.0,
        })
        .collect();
    let ends = [[End::Periodic; 2]; 2];
    let (original, _, _) = balance([4, 4], [1.0; 2], ends, &velocity, &material, false);
    let (corrected, dissipation, _) = balance([4, 4], [1.0; 2], ends, &velocity, &material, true);
    assert!((original - 3.613333333333333).abs() < 2e-13);
    assert!((corrected + 2.366666666666667).abs() < 2e-13);
    assert!((corrected + dissipation).abs() < 2e-13);
}

fn uniform_coefficients<const D: usize>() {
    let velocity: [[f64; D]; 2] =
        std::array::from_fn(|side| std::array::from_fn(|a| ((side * 7 + a + 1) as f64).sin()));
    let gradient = |side: usize, a: usize, j: usize| ((side * 7 + a * 3 + j + 1) as f64).cos();
    for viscosity in [
        Viscosity::default(),
        Viscosity {
            shear: 1.7,
            bulk: 0.3,
        },
    ] {
        for axis in 0..D {
            let original = native(axis, 0.7, velocity, viscosity, gradient);
            let compatible = traction(axis, 0.7, velocity, [viscosity; 2], gradient);
            for a in 0..D {
                assert!((original[a] - compatible[a]).abs() < 3e-15);
            }
        }
    }
}

#[test]
fn uniform_coefficients_match_original_newtonian_stress() {
    uniform_coefficients::<2>();
    uniform_coefficients::<3>();
}

#[test]
fn extreme_midpoint_coefficient_does_not_erase_finite_cell_stress() {
    let cells = [
        Viscosity {
            shear: 1.0,
            bulk: 0.0,
        },
        Viscosity {
            shear: 2.0,
            bulk: 0.0,
        },
    ];
    let face = Viscosity {
        shear: 1e30,
        bulk: 0.0,
    };
    let velocity = [[0.0; 3], [0.0, 1.0, 0.0]];
    let gradient = |_: usize, _: usize, _: usize| 0.0;
    let native = native(0, 1.0, velocity, face, gradient);
    let compatible = traction(0, 1.0, velocity, cells, gradient);
    assert_eq!(compatible, [0.0, 1.5, 0.0]);
    // A subtraction-based correction loses all of this representable stress.
    assert_eq!(native[1] + (compatible[1] - native[1]), 0.0);
    let source_nonmechanical_energy = 7.0;
    assert_eq!(source_nonmechanical_energy - 0.5 * compatible[1], 6.25);
    let native_energy = source_nonmechanical_energy - 0.5 * native[1];
    assert_eq!(native_energy + 0.5 * (native[1] - compatible[1]), 0.0);
}

#[test]
fn cell_quadrature_recovers_affine_stress() {
    let gradient = [[0.3, -0.2, 0.7], [0.8, -0.1, 0.4], [-0.6, 0.9, 0.2]];
    let cells = [
        Viscosity {
            shear: 2.0,
            bulk: 0.3,
        },
        Viscosity {
            shear: 5.0,
            bulk: 0.8,
        },
    ];
    let mean = Viscosity {
        shear: 3.5,
        bulk: 0.55,
    };
    for axis in 0..3 {
        let left = [0.1, -0.2, 0.5];
        let right = std::array::from_fn(|a| left[a] + 0.7 * gradient[a][axis]);
        let read = |_: usize, a: usize, j: usize| gradient[a][j];
        let compatible = traction(axis, 0.7, [left, right], cells, read);
        let expected = native(axis, 0.7, [left, right], mean, read);
        for a in 0..3 {
            assert!((compatible[a] - expected[a]).abs() < 4e-14);
        }
    }
}

#[test]
fn coupled_variable_shear_and_bulk_stress_converges_at_second_order() {
    let mut previous: Option<f64> = None;
    for n in [16, 32, 64] {
        let h = 2.0 * std::f64::consts::PI / n as f64;
        let mut velocity = Vec::new();
        let mut material = Vec::new();
        let mut exact = Vec::new();
        for i in 0..n * n {
            let x = (i % n) as f64 * h + 0.5 * h;
            let y = (i / n) as f64 * h + 0.5 * h;
            let (sx, cx, sy, cy) = (x.sin(), x.cos(), y.sin(), y.cos());
            let u = [sx * sy, cx * sy];
            let mu = 2.0 + 0.25 * cx + 0.15 * sy;
            let bulk = 0.2 + 0.1 * cy;
            let lambda = bulk - (2.0 / 3.0) * mu;
            let divergence = cx * (sy + cy);
            let grad_divergence = [-sx * (sy + cy), cx * (cy - sy)];
            let grad_mu = [-0.25 * sx, 0.15 * cy];
            let grad_lambda = [sx / 6.0, -0.1 * sy - 0.1 * cy];
            let gradient = [[cx * sy, sx * cy], [-sx * sy, cx * cy]];
            velocity.push(u);
            material.push(Viscosity { shear: mu, bulk });
            exact.push(std::array::from_fn::<_, 2, _>(|a| {
                -2.0 * mu * u[a]
                    + (mu + lambda) * grad_divergence[a]
                    + (0..2)
                        .map(|j| grad_mu[j] * (gradient[a][j] + gradient[j][a]))
                        .sum::<f64>()
                    + grad_lambda[a] * divergence
            }));
        }
        let (_, _, force) = balance(
            [n; 2],
            [h; 2],
            [[End::Periodic; 2]; 2],
            &velocity,
            &material,
            true,
        );
        let error = (force
            .iter()
            .zip(&exact)
            .flat_map(|(got, want)| (0..2).map(move |a| (got[a] - want[a]).powi(2)))
            .sum::<f64>()
            / (2 * n * n) as f64)
            .sqrt();
        if let Some(coarse) = previous {
            let ratio = coarse / error;
            assert!((3.8..4.1).contains(&ratio), "refinement ratio {ratio}");
        }
        previous = Some(error);
    }
}

fn prescribed_wall_balance<const D: usize>(shape: [usize; D], spacing: [f64; D]) {
    let n = shape.iter().product();
    let velocity: Vec<[f64; D]> = (0..n)
        .map(|i| std::array::from_fn(|a| ((11 * i + 5 * a + 1) as f64).sin()))
        .collect();
    let material: Vec<_> = (0..n)
        .map(|i| Viscosity {
            shear: 10.0_f64.powf(4.0 * ((3 * i + 1) as f64).sin()),
            bulk: if i % 3 == 0 {
                0.0
            } else {
                10.0_f64.powf(2.0 * ((7 * i + 2) as f64).sin())
            },
        })
        .collect();
    let choices = [
        [End::NoSlip; 2],
        [End::NoSlip, End::Slip],
        [End::Slip, End::NoSlip],
        [End::Periodic; 2],
        [End::Zero; 2],
    ];
    for motion in [0.0, 1.0] {
        let walls = std::array::from_fn(|axis| {
            std::array::from_fn(|side| {
                std::array::from_fn(|a| {
                    if a == axis {
                        0.0
                    } else {
                        motion * ((axis * 11 + side * 7 + a + 1) as f64).cos()
                    }
                })
            })
        });
        for mut encoded in 0..5_usize.pow(D as u32) {
            let ends = std::array::from_fn(|_| {
                let pair = choices[encoded % 5];
                encoded /= 5;
                pair
            });
            let (rate, dissipation, _, power) =
                balance_with_walls(shape, spacing, ends, walls, &velocity, &material, true);
            assert!(dissipation >= 0.0);
            let scale = (power.abs() + rate.abs() + dissipation).max(1.0);
            assert!(
                (power - rate - dissipation).abs() < 3e-14 * scale,
                "wall power {power}, kinetic rate {rate}, quadratic dissipation {dissipation}"
            );
            if motion == 0.0 {
                assert_eq!(power, 0.0);
                assert!(rate <= 0.0);
            }
        }
    }
}

#[test]
fn prescribed_walls_close_power_at_corners_and_singleton_axes() {
    // Every plane orientation, mixed ends, distinct corner velocities and
    // eight-decade positive viscosity; no corner velocity averaging is needed.
    prescribed_wall_balance([4, 5], [0.7, 1.3]);
    prescribed_wall_balance([1, 1], [0.7, 1.3]);
    prescribed_wall_balance([4, 3, 2], [0.7, 1.3, 0.9]);
    prescribed_wall_balance([1, 3, 1], [0.7, 1.3, 0.9]);
    prescribed_wall_balance([1, 1, 1], [0.7, 1.3, 0.9]);
}

#[test]
fn moving_wall_can_supply_kinetic_energy_with_positive_viscous_dissipation() {
    let mut walls = [[[0.0; 2]; 2]; 2];
    walls[0][1][1] = 1.0;
    let (rate, dissipation, force, power) = balance_with_walls(
        [1, 1],
        [1.0; 2],
        [[End::NoSlip; 2], [End::Periodic; 2]],
        walls,
        &[[0.0, 0.25]],
        &[Viscosity {
            shear: 2.0,
            bulk: 0.0,
        }],
        true,
    );
    assert_eq!(force, vec![[0.0, 2.0]]);
    assert_eq!(rate, 0.5);
    assert_eq!(power, 3.0);
    assert_eq!(dissipation, 2.5);
    assert_eq!(power - rate, dissipation);
}

#[test]
fn variable_viscosity_wall_driven_couette_converges_to_the_resistance_integral() {
    // mu(x)=1+x, U(0)=0, U(1)=1: tau=1/ln(2), U(x)=ln(1+x)/ln(2).
    // Steady cell values and wall load converge quadratically in this example;
    // inserting the continuum profile in the half-cell wall traction is only
    // first-order locally. Do not infer a general second-order wall theorem.
    let exact_traction = 1.0 / 2.0_f64.ln();
    let mut previous: Option<[f64; 3]> = None;
    for n in [8, 16, 32, 64] {
        let h = 1.0 / n as f64;
        let material: Vec<_> = (0..n)
            .map(|i| Viscosity {
                shear: 1.0 + (i as f64 + 0.5) * h,
                bulk: 0.0,
            })
            .collect();
        let mut resistance = 0.5 * h / material[0].shear + 0.5 * h / material[n - 1].shear;
        for pair in material.windows(2) {
            resistance += h / (0.5 * pair[0].shear + 0.5 * pair[1].shear);
        }
        let stress = 1.0 / resistance;
        let mut partial = 0.5 * h / material[0].shear;
        let mut velocity = Vec::new();
        let mut profile_error = 0.0_f64;
        for i in 0..n {
            let u = stress * partial;
            velocity.push([0.0, u]);
            let exact = (1.0 + (i as f64 + 0.5) * h).ln() * exact_traction;
            profile_error = profile_error.max((u - exact).abs());
            if i + 1 < n {
                partial += h / (0.5 * material[i].shear + 0.5 * material[i + 1].shear);
            }
        }
        let mut walls = [[[0.0; 2]; 2]; 2];
        walls[0][1][1] = 1.0;
        let (rate, dissipation, force, power) = balance_with_walls(
            [n, 1],
            [h, 1.0],
            [[End::NoSlip; 2], [End::Periodic; 2]],
            walls,
            &velocity,
            &material,
            true,
        );
        assert!(force.iter().flatten().all(|v| v.abs() < 1e-9));
        assert!((power * h - stress).abs() < 1e-12);
        assert!((power - rate - dissipation).abs() < 1e-10);
        let local = no_slip_traction(
            0,
            -1.0,
            0.5 * h,
            [0.0, (1.0 + 0.5 * h).ln() * exact_traction],
            [0.0; 2],
            material[0],
            |_, _| 0.0,
        )[1];
        let errors = [
            profile_error,
            (stress - exact_traction).abs(),
            (local - exact_traction).abs(),
        ];
        if let Some(coarse) = previous {
            let ratios = std::array::from_fn::<_, 3, _>(|i| coarse[i] / errors[i]);
            assert!((3.5..4.2).contains(&ratios[0]));
            assert!((3.5..4.2).contains(&ratios[1]));
            assert!((1.9..2.1).contains(&ratios[2]));
            eprintln!("Couette n={n}: profile/load/local-wall refinement {ratios:?}");
        }
        previous = Some(errors);
    }
}
