use super::*;

fn limits(cells: [usize; 2]) -> RectangularThermalLimits {
    RectangularThermalLimits {
        cells,
        maximum_cell_diffusion_phase: 1.0,
        maximum_relative_residual: 2e-8,
        maximum_payload_bytes: 512 * 1024 * 1024,
        sparse: fs_sparse::direct::DirectLimits {
            max_symbolic_entries: 8_000_000,
            max_analysis_work: 500_000_000,
            max_numeric_work: 1_000_000_000_000,
            max_front_dim: 4096,
        },
    }
}
fn layer(
    thickness: [f64; 4],
    cells: [usize; 4],
    conductivity: [[f64; 2]; 2],
    capacity: f64,
) -> Layer {
    Layer {
        thickness,
        cells,
        conductivity,
        capacity,
    }
}
fn solve(
    dim: [f64; 2],
    kg: f64,
    cap: C64,
    omega: f64,
    layers: &[Layer],
    baths: [Bath; 4],
    cells: [usize; 2],
) -> JacketField {
    prepare(
        dim,
        kg,
        cap,
        omega,
        layers,
        baths,
        limits(cells),
        &mut || false,
    )
    .unwrap()
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-25), "{a:e} != {b:e}");
}

// Independent continuous homogeneous rectangle, driven ONLY inside the gas
// subrectangle. Its double-sine forcing coefficients are integrated exactly.
fn sine_reference(
    dim: [f64; 2],
    thickness: [f64; 4],
    k: f64,
    c: f64,
    omega: f64,
    point: [f64; 2],
) -> (C64, [C64; 3]) {
    let size = [
        dim[0] + thickness[0] + thickness[1],
        dim[1] + thickness[2] + thickness[3],
    ];
    let lower = [thickness[0], thickness[2]];
    let point = [point[0] + lower[0], point[1] + lower[1]];
    let mut mean = C64::ZERO;
    let mut field = [C64::ZERO; 3];
    for m in 1..=180 {
        let a = m as f64 * std::f64::consts::PI / size[0];
        let ix = ((a * lower[0]).cos() - (a * (lower[0] + dim[0])).cos()) / a;
        for n in 1..=180 {
            let b = n as f64 * std::f64::consts::PI / size[1];
            let iy = ((b * lower[1]).cos() - (b * (lower[1] + dim[1])).cos()) / b;
            let response = C64::new(0.0, omega * 4.0 * ix * iy / (size[0] * size[1]))
                / C64::new(k * (a * a + b * b), omega * c);
            mean = mean + response.scale(ix * iy / (dim[0] * dim[1]));
            field[0] = field[0] + response.scale((a * point[0]).sin() * (b * point[1]).sin());
            field[1] = field[1] + response.scale(a * (a * point[0]).cos() * (b * point[1]).sin());
            field[2] = field[2] + response.scale(b * (a * point[0]).sin() * (b * point[1]).cos());
        }
    }
    (mean, field)
}

#[test]
fn original_gas_and_solid_domains_converge_to_independent_forced_rectangle() {
    let dim = [1.0, 1.5];
    let thickness = [0.25, 0.25, 0.375, 0.375];
    let point = [8.5 / 32.0, 19.5 / 32.0 * 1.5];
    let reference = sine_reference(dim, thickness, 0.7, 2.0, 0.9, point);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let value = solve(
            dim,
            0.7,
            C64::from_re(2.0),
            0.9,
            &[layer(thickness, [n / 4; 4], [[0.7, 0.0], [0.0, 0.7]], 2.0)],
            [Bath::Isothermal; 4],
            [n, n],
        );
        errors.push(relative(value.mean_gas, reference.0));
        assert!(value.forms[0].conduction > 0.0 && value.forms[1].conduction > 0.0);
        if n == 32 {
            let actual = value.at(0, point).unwrap();
            assert!(relative(actual[0], reference.1[0]) < 0.002);
            assert!(relative(actual[1], reference.1[1]) < 0.006);
            assert!(relative(actual[2], reference.1[2]) < 0.006);
            let gas = value.at(0, [0.0, 0.63]).unwrap();
            let solid = value.at(1, [0.0, 0.63]).unwrap();
            assert!(relative(gas[0], solid[0]) < 1e-13);
        }
    }
    eprintln!("conjugate sine-reference mean errors: {errors:?}");
    assert!(errors[0] / errors[1] > 3.7 && errors[1] / errors[2] > 3.7);
    assert!(errors[2] < 0.002);
}

#[test]
fn original_full_tensor_corner_fields_and_each_heat_channel_are_retained() {
    let baths = [
        Bath::Conductance {
            conductance_w_m2_k: 0.3,
        },
        Bath::Conductance {
            conductance_w_m2_k: 0.7,
        },
        Bath::Adiabatic,
        Bath::Isothermal,
    ];
    let layers = [
        layer([0.25; 4], [4; 4], [[2.0, 0.6], [0.6, 0.9]], 3.0),
        layer([0.125; 4], [2; 4], [[0.6, -0.12], [-0.12, 1.1]], 1.3),
    ];
    let cap = C64::new(2.0, -0.35);
    let value = solve([1.0, 1.5], 0.7, cap, 0.9, &layers, baths, [16, 24]);
    let corner = value.at(1, [-0.13, -0.17]).unwrap();
    assert!(corner[1].abs() > 1e-4 && corner[2].abs() > 1e-4);
    assert!(
        value
            .forms
            .iter()
            .all(|f| f.conduction > 0.0 && f.field > 0.0)
    );
    assert!(value.relative_heat_defect < 1e-11 && value.relative_work_defect < 1e-11);
    let molecular = -0.9 * cap.im * value.forms[0].field;
    let bath = value.bath_forms[0] * 0.3 + value.bath_forms[1] * 0.7;
    assert!(molecular > 0.0 && bath > 0.0);
    let output = value.forms.iter().map(|f| f.conduction).sum::<f64>() + molecular + bath;
    near(output, 0.9 * value.forms[0].temperature_integral.im, 1e-11);
    let mass = value.forms[0].temperature_integral * cap
        + value.forms[1].temperature_integral.scale(3.0)
        + value.forms[2].temperature_integral.scale(1.3);
    assert!((C64::new(0.0, 0.9) * mass + value.bath_flux - C64::new(0.0, 1.35)).abs() < 1e-11);
    // The isotropic normal-only replacement measurably changes the corner response.
    let mut diagonal = layers;
    for v in &mut diagonal {
        v.conductivity[0][1] = 0.0;
        v.conductivity[1][0] = 0.0;
    }
    let dropped = solve([1.0, 1.5], 0.7, cap, 0.9, &diagonal, baths, [16, 24]);
    assert!(relative(value.mean_gas, dropped.mean_gas) > 1e-3);
}

#[test]
fn reflection_preserves_original_rotated_tensors_and_asymmetric_baths() {
    let layer = layer(
        [0.25, 0.5, 0.125, 0.25],
        [4, 8, 2, 4],
        [[4.0, 1.2], [1.2, 2.0]],
        3.0,
    );
    let baths = [0.3, 1.1, 0.8, 0.5].map(|g| Bath::Conductance {
        conductance_w_m2_k: g,
    });
    let a = solve(
        [1.0, 1.5],
        0.7,
        C64::from_re(2.0),
        0.9,
        &[layer],
        baths,
        [16, 24],
    );
    let mut mirror = layer;
    mirror.thickness.swap(0, 1);
    mirror.cells.swap(0, 1);
    mirror.conductivity[0][1] *= -1.0;
    mirror.conductivity[1][0] *= -1.0;
    let mut mb = baths;
    mb.swap(0, 1);
    let b = solve(
        [1.0, 1.5],
        0.7,
        C64::from_re(2.0),
        0.9,
        &[mirror],
        mb,
        [16, 24],
    );
    assert!(relative(a.mean_gas, b.mean_gas) < 2e-10);
    for (region, p) in [(0, [0.27, 0.71]), (1, [-0.13, 0.27]), (1, [-0.17, -0.07])] {
        let x = a.at(region, p).unwrap();
        let y = b.at(region, [1.0 - p[0], p[1]]).unwrap();
        assert!(relative(x[0], y[0]) < 2e-10);
        assert!(relative(x[1], -y[1]) < 2e-10);
        assert!(relative(x[2], y[2]) < 2e-10);
    }
}

#[test]
fn splitting_an_original_homogeneous_solid_preserves_every_discrete_field() {
    let one = layer([0.5; 4], [8; 4], [[2.0, 0.4], [0.4, 1.0]], 3.0);
    let half = layer([0.25; 4], [4; 4], one.conductivity, one.capacity);
    let baths = [Bath::Conductance {
        conductance_w_m2_k: 0.7,
    }; 4];
    let a = solve(
        [1.0, 1.5],
        0.7,
        C64::new(2.0, -0.2),
        0.9,
        &[one],
        baths,
        [16, 24],
    );
    let b = solve(
        [1.0, 1.5],
        0.7,
        C64::new(2.0, -0.2),
        0.9,
        &[half, half],
        baths,
        [16, 24],
    );
    assert_eq!(a.axes, b.axes);
    assert_eq!(a.stored, b.stored);
    assert_eq!(a.mean_gas, b.mean_gas);
    assert_eq!(a.bath_flux, b.bath_flux);
    near(
        a.forms[1].conduction,
        b.forms[1].conduction + b.forms[2].conduction,
        1e-14,
    );
    near(a.forms[1].field, b.forms[1].field + b.forms[2].field, 1e-14);
}

#[test]
fn resolved_thin_jacket_converges_to_existing_finite_capacity_walls() {
    let g = 0.8;
    let capacity = 0.3;
    let omega = 0.9;
    let kg = 0.7;
    let old = super::super::super::physics::prepare(
        [1.0, 1.5],
        C64::new(0.0, omega * 2.0 / kg),
        [C64::new(g, omega * capacity).scale(1.0 / kg); 4],
        [super::super::super::RectangularThermalWall::Finite {
            conductance_w_m2_k: g,
            areal_heat_capacity_j_m2_k: capacity,
        }; 4],
        limits([24, 32]),
        &mut || false,
    )
    .unwrap();
    let mut errors = Vec::new();
    for d in [0.02, 0.01, 0.005] {
        let value = solve(
            [1.0, 1.5],
            kg,
            C64::from_re(2.0),
            omega,
            &[layer(
                [d; 4],
                [2; 4],
                [[0.7, 0.0], [0.0, 0.7]],
                capacity / d,
            )],
            [Bath::Conductance {
                conductance_w_m2_k: g,
            }; 4],
            [24, 32],
        );
        errors.push(relative(value.mean_gas, old.mean.scale(0.5)));
    }
    eprintln!("resolved thin-jacket limit errors: {errors:?}");
    assert!(errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8);
    assert!(errors[2] < 0.02);
}

#[test]
fn physical_length_and_conductivity_scaling_preserves_temperature_and_work() {
    let original = layer([0.25; 4], [4; 4], [[2.0, 0.6], [0.6, 0.9]], 3.0);
    let baths = [Bath::Conductance {
        conductance_w_m2_k: 0.7,
    }; 4];
    let a = solve(
        [1.0, 1.5],
        0.7,
        C64::new(2.0, -0.2),
        0.9,
        &[original],
        baths,
        [16, 24],
    );
    let scale = 0.125;
    let mut small = original;
    small.thickness = small.thickness.map(|t| t * scale);
    small.conductivity = small.conductivity.map(|r| r.map(|k| k * scale * scale));
    let b = solve(
        [scale, 1.5 * scale],
        0.7 * scale * scale,
        C64::new(2.0, -0.2),
        0.9,
        &[small],
        [Bath::Conductance {
            conductance_w_m2_k: 0.7 * scale,
        }; 4],
        [16, 24],
    );
    assert!(relative(a.mean_gas, b.mean_gas) < 1e-11);
    let x = a.at(1, [-0.13, 0.71]).unwrap();
    let y = b.at(1, [-0.13 * scale, 0.71 * scale]).unwrap();
    assert!(relative(x[0], y[0]) < 1e-11);
    for d in 1..3 {
        assert!(relative(x[d], y[d].scale(scale)) < 1e-11);
    }
    for region in 0..2 {
        near(
            a.forms[region].conduction * scale * scale,
            b.forms[region].conduction,
            1e-11,
        );
    }
}

#[test]
fn opaque_gas_and_solid_keep_original_gradients_and_nonzero_heat() {
    let value = solve(
        [0.2, 0.3],
        0.7,
        C64::new(2.0, -0.05),
        2000.0,
        &[layer([0.02; 4], [4; 4], [[2.0, 0.4], [0.4, 1.0]], 3.0)],
        [Bath::Isothermal; 4],
        [24, 32],
    );
    assert!(value.lift.is_some());
    let p = [0.079, 0.118];
    let f = value.at(0, p).unwrap();
    assert!(f[1].abs() > 1e-5 && f[2].abs() > 1e-5);
    for axis in 0..2 {
        let mut low = p;
        let mut high = p;
        low[axis] -= 1e-7;
        high[axis] += 1e-7;
        let difference =
            (value.at(0, high).unwrap()[0] - value.at(0, low).unwrap()[0]).scale(0.5e7);
        assert!(relative(difference, f[1 + axis]) < 1e-5);
    }
    assert!(value.forms.iter().all(|f| f.conduction > 0.0));
}

#[test]
fn original_material_geometry_limits_and_cancelled_preparation_refuse() {
    let good = layer([0.25; 4], [4; 4], [[2.0, 0.4], [0.4, 1.0]], 3.0);
    let build = |layers: &[Layer], limits, cancelled: &mut dyn FnMut() -> bool| {
        prepare(
            [1.0, 1.5],
            0.7,
            C64::from_re(2.0),
            0.9,
            layers,
            [Bath::Isothermal; 4],
            limits,
            &mut || cancelled(),
        )
    };
    assert!(build(&[], limits([16, 24]), &mut || false).is_err());
    let mut invalid = good;
    invalid.conductivity = [[1.0, 2.0], [2.0, 1.0]];
    assert!(build(&[invalid], limits([16, 24]), &mut || false).is_err());
    invalid = good;
    invalid.thickness[0] = 0.0;
    assert!(build(&[invalid], limits([16, 24]), &mut || false).is_err());
    let mut limited = limits([16, 24]);
    limited.maximum_payload_bytes = 32;
    assert!(build(&[good], limited, &mut || false).is_err());
    limited = limits([16, 24]);
    limited.maximum_cell_diffusion_phase = 1e-5;
    assert!(build(&[good], limited, &mut || false).is_err());
    assert!(build(&[good], limits([128, 128]), &mut || false).is_err());
    for stop in [0, 9, 100, 700] {
        let mut count = 0;
        assert!(
            build(&[good], limits([16, 24]), &mut || {
                count += 1;
                count > stop
            })
            .is_err()
        );
    }
    let a = build(&[good], limits([16, 24]), &mut || false).unwrap();
    let b = build(&[good], limits([16, 24]), &mut || false).unwrap();
    assert_eq!(a.stored, b.stored);
    assert!(a.at(0, [-0.1, 0.5]).is_err());
    assert!(a.at(1, [0.4, 0.5]).is_err());
    assert!(a.at(2, [-0.1, 0.5]).is_err());
}
