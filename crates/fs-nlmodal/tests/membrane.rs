//! Tensor membrane elasticity must survive modal reduction with its physical
//! shear, Poisson and orientation dependence and without changing the isotropic limit.

use fs_material::sheet::SheetMatrix;
use fs_nlmodal::{
    SampledPlateMode, SineMode, SosModalStorage, VkModel, VkPlateParams, VkSampledPlateParams,
    von_karman_sampled_plate, von_karman_sampled_plate_with_membrane,
};

fn geometry() -> VkPlateParams {
    VkPlateParams {
        lx: 0.37,
        ly: 0.21,
        h: 0.004,
        young: 12.0e9,
        nu: 0.3,
        rho: 600.0,
        pretension_n_m: 0.0,
    }
}

fn modes() -> Vec<SampledPlateMode> {
    let side = 33;
    [(1.0, 1.0), (2.0, 2.0)]
        .into_iter()
        .map(|(m, n)| SampledPlateMode {
            omega: 300.0 * m,
            nx: side,
            ny: side,
            w: (0..side * side)
                .map(|k| {
                    let x = (k % side) as f64 / (side - 1) as f64;
                    let y = (k / side) as f64 / (side - 1) as f64;
                    (m * core::f64::consts::PI * x).sin() * (n * core::f64::consts::PI * y).sin()
                })
                .collect(),
        })
        .collect()
}

fn tensor_params(p: &VkPlateParams, membrane: SheetMatrix) -> VkSampledPlateParams {
    VkSampledPlateParams {
        lx: p.lx,
        ly: p.ly,
        areal_density: p.rho * p.h,
        membrane,
    }
}

fn principal(e1: f64, e2: f64, nu12: f64, g12: f64, h: f64) -> SheetMatrix {
    let denom = 1.0 - nu12 * nu12 * e2 / e1;
    SheetMatrix::orthotropic(
        e1 * h / denom,
        e2 * h / denom,
        nu12 * e2 * h / denom,
        g12 * h,
    )
    .unwrap()
}

fn brackets(storage: &SosModalStorage, q: &[f64]) -> Vec<f64> {
    storage
        .channels
        .iter()
        .map(|channel| {
            let mut value = 0.0;
            for (i, qi) in q.iter().enumerate() {
                for (j, qj) in q.iter().enumerate() {
                    value += qi * channel.coupling[i * q.len() + j] * qj;
                }
            }
            value
        })
        .collect()
}

fn quartic(model: &VkModel, q: &[f64]) -> f64 {
    brackets(&model.storage, q)
        .iter()
        .zip(&model.storage.channels)
        .map(|(e, channel)| 0.25 * channel.coefficient * e * e)
        .sum()
}

#[test]
fn tensor_sampled_plate_recovers_the_isotropic_membrane_energy() {
    let stress = [
        SineMode { m: 1, n: 1 },
        SineMode { m: 2, n: 2 },
        SineMode { m: 3, n: 1 },
    ];
    for nu in [0.0, 0.3, 0.49] {
        let mut p = geometry();
        p.nu = nu;
        let legacy = von_karman_sampled_plate(&p, &modes(), &stress).unwrap();
        let tensor = principal(p.young, p.young, nu, p.young / (2.0 * (1.0 + nu)), p.h);
        let full =
            von_karman_sampled_plate_with_membrane(&tensor_params(&p, tensor), &modes(), &stress)
                .unwrap();
        assert_eq!(full.storage.omegas, legacy.storage.omegas);
        for q in [[0.003, -0.001], [0.0, 0.002], [0.001, 0.004]] {
            let expected = quartic(&legacy, &q);
            let actual = quartic(&full, &q);
            assert!((actual / expected - 1.0).abs() < 2.0e-13, "nu={nu}");
        }
    }
}

#[test]
fn orthotropic_hardening_retains_independent_shear_and_poisson_elasticity() {
    let p = geometry();
    let (e1, e2, nu12) = (14.0e9, 1.1e9, 0.65);
    // A major orthotropic Poisson ratio above 1/2 is admissible here:
    // nu12*nu21 < 1. Applying an isotropic scalar admission would reject it.
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 2, n: 2 }];
    let raw = von_karman_sampled_plate(&p, &modes(), &stress).unwrap();
    let q = [0.003, -0.001];
    let e = brackets(&raw.storage, &q);
    let mut energies = Vec::new();
    for g12 in [0.3e9, 1.5e9] {
        let tensor = principal(e1, e2, nu12, g12, p.h);
        let full =
            von_karman_sampled_plate_with_membrane(&tensor_params(&p, tensor), &modes(), &stress)
                .unwrap();
        let expected: f64 = stress
            .iter()
            .zip(&e)
            .map(|(mode, e)| {
                let kx = mode.m as f64 * core::f64::consts::PI / p.lx;
                let ky = mode.n as f64 * core::f64::consts::PI / p.ly;
                let compliance = (ky.powi(4) / e1
                    + kx.powi(4) / e2
                    + kx.powi(2) * ky.powi(2) * (1.0 / g12 - 2.0 * nu12 / e1))
                    / p.h;
                e * e / (8.0 * compliance)
            })
            .sum();
        let actual = quartic(&full, &q);
        assert!((actual / expected - 1.0).abs() < 3.0e-13);
        energies.push(actual);
    }
    assert!(
        energies[1] > 1.15 * energies[0],
        "G12 must affect hardening"
    );
}

// Independent area integration in material axes: no stiffness inversion or
// analytic sine/cosine overlap from the production implementation is reused.
fn integrated_compliance(p: &VkPlateParams, angle: f64, stress: &[SineMode; 2]) -> [f64; 4] {
    let nodes: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683_1,
        0.0,
        0.538_469_310_105_683_1,
        0.906_179_845_938_664,
    ];
    let weights = [
        0.236_926_885_056_189_1,
        0.478_628_670_499_366_5,
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
    ];
    let points: Vec<_> = (0..8)
        .flat_map(|segment| {
            nodes.into_iter().zip(weights).map(move |(node, weight)| {
                ((segment as f64 + (node + 1.0) * 0.5) / 8.0, weight / 16.0)
            })
        })
        .collect();
    let (s, c) = angle.sin_cos();
    let mut matrix = [0.0; 4];
    let norm = (4.0 / (p.lx * p.ly)).sqrt();
    for &(x, wx) in &points {
        for &(y, wy) in &points {
            let local: Vec<_> = stress
                .iter()
                .map(|mode| {
                    let kx = mode.m as f64 * core::f64::consts::PI / p.lx;
                    let ky = mode.n as f64 * core::f64::consts::PI / p.ly;
                    let ss = norm * (kx * p.lx * x).sin() * (ky * p.ly * y).sin();
                    let cc = norm * (kx * p.lx * x).cos() * (ky * p.ly * y).cos();
                    let (nx, ny, nxy) = (-ky * ky * ss, -kx * kx * ss, -kx * ky * cc);
                    [
                        c * c * nx + s * s * ny + 2.0 * c * s * nxy,
                        s * s * nx + c * c * ny - 2.0 * c * s * nxy,
                        -c * s * nx + c * s * ny + (c * c - s * s) * nxy,
                    ]
                })
                .collect();
            for i in 0..2 {
                for j in 0..2 {
                    let (a, b) = (local[i], local[j]);
                    let energy = a[0] * b[0] / 14.0e9 + a[1] * b[1] / 1.1e9
                        - 0.65 / 14.0e9 * (a[0] * b[1] + a[1] * b[0])
                        + a[2] * b[2] / 0.3e9;
                    matrix[i * 2 + j] += wx * wy * p.lx * p.ly * energy / p.h;
                }
            }
        }
    }
    matrix
}

#[test]
fn rotated_membrane_matches_complementary_energy_and_retains_airy_mode_coupling() {
    let p = geometry();
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 2, n: 2 }];
    let q = [0.003, -0.001];
    let raw = von_karman_sampled_plate(&p, &modes(), &stress).unwrap();
    let e = brackets(&raw.storage, &q);
    let angle = 0.37;
    let tensor = principal(14.0e9, 1.1e9, 0.65, 0.3e9, p.h)
        .rotated(angle)
        .unwrap();
    let params = tensor_params(&p, tensor);
    let full = von_karman_sampled_plate_with_membrane(&params, &modes(), &stress).unwrap();
    let c = integrated_compliance(&p, angle, &stress);
    let determinant = c[0] * c[3] - c[1] * c[2];
    let expected =
        (c[3] * e[0] * e[0] - 2.0 * c[1] * e[0] * e[1] + c[0] * e[1] * e[1]) / (8.0 * determinant);
    let actual = quartic(&full, &q);
    assert!(
        (actual / expected - 1.0).abs() < 2.0e-11,
        "{actual} vs {expected}"
    );
    let diagonal_only = (e[0] * e[0] / c[0] + e[1] * e[1] / c[3]) / 8.0;
    assert!((diagonal_only / expected - 1.0).abs() > 0.01);
    let reversed =
        von_karman_sampled_plate_with_membrane(&params, &modes(), &[stress[1], stress[0]]).unwrap();
    assert!((quartic(&reversed, &q) / expected - 1.0).abs() < 2.0e-11);
    let mut scaled = params;
    scaled.membrane =
        SheetMatrix::from_coefficients(tensor.coefficients().map(|v| 2.5 * v)).unwrap();
    let stiffer = von_karman_sampled_plate_with_membrane(&scaled, &modes(), &stress).unwrap();
    assert!((quartic(&stiffer, &q) / actual - 2.5).abs() < 2.0e-12);
}

#[test]
fn sampled_tensor_mass_normalization_preserves_physical_displacement_energy() {
    let p = geometry();
    let tensor = principal(14.0e9, 1.1e9, 0.65, 0.3e9, p.h)
        .rotated(0.37)
        .unwrap();
    let mut params = tensor_params(&p, tensor);
    let stress = [SineMode { m: 1, n: 1 }, SineMode { m: 2, n: 2 }];
    let baseline = von_karman_sampled_plate_with_membrane(&params, &modes(), &stress).unwrap();
    let q = [0.003, -0.001];
    params.areal_density *= 3.0;
    let heavier = von_karman_sampled_plate_with_membrane(&params, &modes(), &stress).unwrap();
    let physical_q = q.map(|v| v * 3.0f64.sqrt());
    assert!((quartic(&heavier, &physical_q) / quartic(&baseline, &q) - 1.0).abs() < 3.0e-13);
    params.areal_density = f64::INFINITY;
    assert!(von_karman_sampled_plate_with_membrane(&params, &modes(), &stress).is_err());
}
