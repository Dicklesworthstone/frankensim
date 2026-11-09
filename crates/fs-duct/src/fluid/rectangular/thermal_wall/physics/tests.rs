//! G1/G3: independent diffusion solutions, physical work and geometry refinement.
use super::*;
use fs_sparse::direct::DirectLimits;

fn limits(cells: [usize; 2]) -> RectangularThermalLimits {
    RectangularThermalLimits {
        cells,
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-8,
        maximum_payload_bytes: 512 * 1024 * 1024,
        sparse: DirectLimits {
            max_symbolic_entries: 8_000_000,
            max_analysis_work: 500_000_000,
            max_numeric_work: 100_000_000_000,
            max_front_dim: 4096,
        },
    }
}
fn solve(
    dim: [f64; 2],
    z: C64,
    beta: [C64; 4],
    walls: [Wall; 4],
    cells: [usize; 2],
) -> ThermalField {
    prepare(dim, z, beta, walls, limits(cells), &mut || false).unwrap()
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn wall(beta: C64) -> Wall {
    Wall::Finite {
        conductance_w_m2_k: beta.re,
        areal_heat_capacity_j_m2_k: beta.im,
    }
}
fn exp(z: C64) -> C64 {
    C64::new(z.re.exp() * z.im.cos(), z.re.exp() * z.im.sin())
}
fn sinh(z: C64) -> C64 {
    (exp(z) - exp(z.scale(-1.0))).scale(0.5)
}
fn cosh(z: C64) -> C64 {
    (exp(z) + exp(z.scale(-1.0))).scale(0.5)
}

// Independent dimensional 1-D Robin ODE: h=1+A cosh(qx)+B sinh(qx).
fn robin(l: f64, z: C64, beta: [C64; 2], x: f64) -> [C64; 3] {
    let q = z.sqrt();
    let s = sinh(q.scale(l));
    let c = cosh(q.scale(l));
    let a = beta[0];
    let b = q.scale(-1.0);
    let d = q * s + beta[1] * c;
    let e = q * c + beta[1] * s;
    let determinant = a * e - b * d;
    let aa = (beta[0].scale(-1.0) * e + b * beta[1]) / determinant;
    let bb = (a * beta[1].scale(-1.0) + beta[0] * d) / determinant;
    [
        C64::ONE + aa * cosh(q.scale(x)) + bb * sinh(q.scale(x)),
        q * (aa * sinh(q.scale(x)) + bb * cosh(q.scale(x))),
        C64::ONE + aa * s / (q.scale(l)) + bb * (c - C64::ONE) / (q.scale(l)),
    ]
}
#[test]
fn insulated_constant_mode_and_zero_capacity_faces_are_exactly_retained() {
    let dims = [0.8, 1.3];
    for z in [
        C64::new(0.0, 1e-20),
        C64::new(0.7, 3.0),
        C64::new(9.0, 30.0),
    ] {
        let a = solve(dims, z, [C64::ZERO; 4], [Wall::Adiabatic; 4], [12, 20]);
        let b = solve(dims, z, [C64::ZERO; 4], [wall(C64::ZERO); 4], [12, 20]);
        assert_eq!(a.stored, b.stored);
        assert_eq!(a.factor_nonzeros, 0);
        assert_eq!(a.gradient_form, 0.0);
        for p in [[0.0, 0.0], [0.13, 0.79], dims] {
            assert_eq!(a.at(p).unwrap(), [C64::ONE, C64::ZERO, C64::ZERO]);
        }
        assert!(relative(a.mean, C64::ONE) < 1e-14);
        near(a.field_form, dims[0] * dims[1], 1e-14);
        for (i, l) in [dims[1], dims[1], dims[0], dims[0]].into_iter().enumerate() {
            near(a.wall_forms[i], l, 1e-14);
        }
    }
}
// A double infinite sine expansion independent of the Q1 assembly/factor.
// 255 odd modes in each direction make the remaining mean tail negligible
// relative to the independently refined bilinear mesh errors below.
fn dirichlet_mean(dim: [f64; 2], z: C64) -> C64 {
    let mut sum = C64::ZERO;
    let pi = core::f64::consts::PI;
    for i in 0..255 {
        for j in 0..255 {
            let m = (2 * i + 1) as f64;
            let n = (2 * j + 1) as f64;
            let lambda = (pi * m / dim[0]).powi(2) + (pi * n / dim[1]).powi(2);
            sum = sum + (z / (z + C64::from_re(lambda))).scale(64.0 / (pi.powi(4) * m * m * n * n));
        }
    }
    sum
}
#[test]
fn four_isothermal_faces_converge_to_independent_rectangular_series() {
    let dim = [1.0, 1.3];
    let z = C64::new(0.5, 6.0);
    let expected = dirichlet_mean(dim, z);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let result = solve(dim, z, [C64::ZERO; 4], [Wall::Isothermal; 4], [n, n]);
        let error = relative(result.mean, expected);
        errors.push(error);
        for p in [[0.0, 0.37], [dim[0], 0.79], [0.34, 0.0], [0.81, dim[1]]] {
            assert!(result.at(p).unwrap()[0].abs() < 2e-15);
        }
        assert!(result.gradient_form > 0.0 && result.factor_nonzeros > 0);
    }
    eprintln!("rectangular all-isothermal refinement: {errors:?}");
    assert!(errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5);
    assert!(errors[2] < 0.002);
}
#[test]
fn independently_stored_unequal_walls_converge_to_robin_ode_fields_and_work() {
    let dim = [1.0, 1.3];
    let z = C64::new(0.7, 3.0);
    let beta = [C64::new(0.5, 0.8), C64::new(2.0, 0.1), C64::ZERO, C64::ZERO];
    let walls = [
        wall(beta[0]),
        wall(beta[1]),
        Wall::Adiabatic,
        Wall::Adiabatic,
    ];
    let expected = robin(dim[0], z, [beta[0], beta[1]], 0.37);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let result = solve(dim, z, beta, walls, [n, n]);
        let value = result.at([0.37, 0.47]).unwrap();
        errors.push(relative(result.mean, expected[2]));
        assert!(value[2].abs() < 1e-11);
        if n == 32 {
            assert!(relative(value[0], expected[0]) < 0.001);
            let midpoint = 11.5 / 32.0;
            let gradient = result.at([midpoint, 0.47]).unwrap()[1];
            let gradient_reference = robin(dim[0], z, [beta[0], beta[1]], midpoint)[1];
            assert!(relative(gradient, gradient_reference) < 0.002);
            let mut gradient = 0.0;
            let mut field = 0.0;
            let segments = 10000;
            for i in 0..=segments {
                let x = dim[0] * i as f64 / segments as f64;
                let f = robin(dim[0], z, [beta[0], beta[1]], x);
                let weight = if i == 0 || i == segments {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                gradient += weight * f[1].norm_sq();
                field += weight * f[0].norm_sq();
            }
            gradient *= dim[0] * dim[1] / (3.0 * segments as f64);
            field *= dim[0] * dim[1] / (3.0 * segments as f64);
            near(result.gradient_form, gradient, 0.001);
            near(result.field_form, field, 0.001);
        }
    }
    eprintln!("rectangular unequal Robin refinement: {errors:?}");
    assert!(errors[0] / errors[1] > 3.6 && errors[1] / errors[2] > 3.6);
    assert!(errors[2] < 0.001);
}
#[test]
fn two_dimensional_molecular_and_four_bath_work_close_original_acoustic_port() {
    // rho=omega=kappa=R=1, T0=2, Cp=3+0.7i under -i omega.
    let dim = [1.0, 1.3];
    let cp = C64::new(3.0, 0.7);
    let z = C64::new(cp.im, cp.re);
    let beta = [
        C64::new(0.5, 0.8),
        C64::new(2.0, 0.1),
        C64::new(0.2, 1.1),
        C64::new(0.9, 0.3),
    ];
    let fields = solve(dim, z, beta, beta.map(wall), [16, 20]);
    let positive = fields.gradient_form
        + z.re * fields.field_form
        + beta
            .iter()
            .zip(fields.wall_forms)
            .map(|(b, f)| b.re * f)
            .sum::<f64>();
    let rhs = (z * fields.mean.conj()).re * dim[0] * dim[1];
    near(positive, rhs, 1e-11);
    let physical_heat = positive / (4.0 * cp.norm_sq());
    let compliance = (C64::ONE - fields.mean / cp.conj()).scale(0.5);
    let input = 0.5 * (C64::new(0.0, dim[0] * dim[1]) * compliance).re;
    near(physical_heat, input, 1e-12);
    assert!(fields.gradient_form > 0.0 && fields.field_form > 0.0);
    assert!(fields.wall_forms.into_iter().all(|v| v > 0.0));
    // Interior derivatives agree with independently differenced local values.
    let p = [0.331, 0.477];
    let eps = 1e-6;
    let v = fields.at(p).unwrap();
    for axis in 0..2 {
        let mut a = p;
        let mut b = p;
        a[axis] += eps;
        b[axis] -= eps;
        let slope = (fields.at(a).unwrap()[0] - fields.at(b).unwrap()[0]).scale(0.5 / eps);
        assert!(relative(slope, v[1 + axis]) < 1e-8);
    }
}
#[test]
fn axis_exchange_and_physical_length_units_preserve_original_solution() {
    let dim = [0.8, 1.3];
    let z = C64::new(0.7, 3.0);
    let beta = [
        C64::new(0.5, 0.8),
        C64::new(2.0, 0.1),
        C64::new(0.2, 1.1),
        C64::new(0.9, 0.3),
    ];
    let a = solve(dim, z, beta, beta.map(wall), [12, 16]);
    let swapped = [beta[2], beta[3], beta[0], beta[1]];
    let b = solve([dim[1], dim[0]], z, swapped, swapped.map(wall), [16, 12]);
    for p in [[0.13, 0.79], [0.8, 0.3], [0.0, 0.0]] {
        let x = a.at(p).unwrap();
        let y = b.at([p[1], p[0]]).unwrap();
        for (i, j) in [(0, 0), (1, 2), (2, 1)] {
            assert!(relative(x[i], y[j]) < 2e-10);
        }
    }
    let scale = 0.001;
    let scaled_beta = beta.map(|v| v.scale(1.0 / scale));
    let c = solve(
        dim.map(|v| v * scale),
        z.scale(1.0 / (scale * scale)),
        scaled_beta,
        scaled_beta.map(wall),
        [12, 16],
    );
    assert!(relative(a.mean, c.mean) < 2e-10);
    near(a.gradient_form, c.gradient_form, 2e-10);
    near(a.field_form, c.field_form / (scale * scale), 2e-10);
    for i in 0..4 {
        near(a.wall_forms[i], c.wall_forms[i] / scale, 2e-10);
    }
}
#[test]
fn finite_wall_limit_recovers_the_original_isothermal_boundary() {
    let dim = [1.0, 1.3];
    let z = C64::new(0.0, 3.0);
    let reference = solve(dim, z, [C64::ZERO; 4], [Wall::Isothermal; 4], [12, 16]);
    let mut previous = f64::INFINITY;
    for g in [1e2, 1e4, 1e6] {
        let beta = [C64::from_re(g); 4];
        let actual = solve(dim, z, beta, beta.map(wall), [12, 16]);
        let error = relative(actual.mean, reference.mean);
        assert!(error < previous * 0.02);
        previous = error;
    }
    assert!(previous < 1e-5);
}
#[test]
fn opaque_thermal_complement_preserves_nonzero_physical_gradients() {
    let dim = [0.1, 0.13];
    let z = C64::new(0.0, 30_000.0);
    let fields = solve(dim, z, [C64::ZERO; 4], [Wall::Isothermal; 4], [32, 40]);
    assert!(fields.complement);
    assert!(fields.mean.re > 0.7 && fields.mean.re < 1.0);
    assert!(fields.gradient_form > 0.0);
    near(
        fields.gradient_form,
        z.im * fields.mean.im * dim[0] * dim[1],
        1e-10,
    );
}
#[test]
fn original_scope_budget_and_cancellation_refuse_without_partial_fields() {
    let dim = [1.0, 1.3];
    let z = C64::new(0.0, 3.0);
    let walls = [Wall::Isothermal; 4];
    let beta = [C64::ZERO; 4];
    assert!(prepare(dim, z, beta, walls, limits([1, 8]), &mut || false).is_err());
    assert!(prepare(dim, z, beta, walls, limits([128, 128]), &mut || false).is_err());
    assert!(
        prepare(
            dim,
            C64::new(-1.0, 3.0),
            beta,
            walls,
            limits([8, 8]),
            &mut || false
        )
        .is_err()
    );
    let mut invalid = limits([8, 8]);
    invalid.maximum_payload_bytes = 1;
    assert!(prepare(dim, z, beta, walls, invalid, &mut || false).is_err());
    let mut bounded = limits([8, 8]);
    bounded.sparse.max_analysis_work = 1;
    assert!(matches!(
        prepare(dim, z, beta, walls, bounded, &mut || false),
        Err(Error::Sparse(_))
    ));
    let mut count = 0;
    assert!(
        prepare(dim, z, beta, walls, limits([8, 8]), &mut || {
            count += 1;
            count > 4
        })
        .is_err()
    );
    let a = solve(dim, z, beta, walls, [8, 8]);
    let b = solve(dim, z, beta, walls, [8, 8]);
    assert_eq!(a.stored, b.stored);
    assert!(a.at([-1e-9, 0.2]).is_err());
    assert!(a.at([0.2, f64::NAN]).is_err());
    let mut phase = limits([8, 8]);
    phase.maximum_cell_diffusion_phase = 1e-9;
    assert!(prepare(dim, z, beta, walls, phase, &mut || false).is_err());
}
