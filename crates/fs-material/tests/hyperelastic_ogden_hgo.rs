//! M0.3-W1a: independent material-point oracles, spectral coalescence,
//! objectivity, dimensioned cards, and physical-domain refusals.
//! Published parameters are never fitted here. Prescribed J=1 tests remove
//! the lateral pressure; finite-bulk traction-free tests are separate.

use fs_material::hyper::{
    FibreInvariants, Hyperelastic, HyperelasticModel, MAX_FIBRE_FAMILIES, MAX_OGDEN_TERMS,
    OgdenTerm, Volumetric,
};
use fs_material::tensor::{mandel_basis, matmul3, rotate_mandel, rotation};
use fs_material::{IsotropicElastic, MaterialError};
use fs_math::det;
use std::time::Instant;

const I: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
// Original mu/alpha convention, not the alternative 2*mu/alpha² convention.
// Ogden (1972), Proc R Soc A 326:565–584, doi:10.1098/rspa.1972.0026;
// reproduced in Holzapfel (2000), Nonlinear Solid Mechanics §6.5.
const TRELOAR: [OgdenTerm; 3] = [
    OgdenTerm {
        mu_pa: 0.63e6,
        alpha: 1.3,
    },
    OgdenTerm {
        mu_pa: 0.0012e6,
        alpha: 5.0,
    },
    OgdenTerm {
        mu_pa: -0.01e6,
        alpha: -2.0,
    },
];
const MU: f64 = 422_500.0;
// HGO (2000), doi:10.1023/A:1010835316564, Fig.14 media parameters:
// c=3 kPa, k1=2.3632 kPa, k2=.8393; our c10=c/2, not c.
const C10: f64 = 1500.0;
const K1: f64 = 2363.2;
const K2: f64 = 0.8393;

fn ogden(ratio: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: TRELOAR.to_vec(),
            volumetric: Volumetric::Quadratic { kappa: MU * ratio },
        },
        32.0,
    )
    .unwrap()
}

fn hgo(
    fibres: Vec<[f64; 3]>,
    dispersion: f64,
    invariants: FibreInvariants,
    ratio: f64,
) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::HolzapfelGasserOgden {
            c10: C10,
            k1: K1,
            k2: K2,
            dispersion,
            fibres,
            invariants,
            volumetric: Volumetric::Quadratic { kappa: C10 * ratio },
        },
        8.0,
    )
    .unwrap()
}

fn diag(x: [f64; 3]) -> [f64; 9] {
    [x[0], 0., 0., 0., x[1], 0., 0., 0., x[2]]
}
fn matrix(x: &[f64; 9]) -> [[f64; 3]; 3] {
    core::array::from_fn(|i| core::array::from_fn(|j| x[3 * i + j]))
}
fn flat(x: [[f64; 3]; 3]) -> [f64; 9] {
    core::array::from_fn(|i| x[i / 3][i % 3])
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn power(x: f64, a: f64) -> f64 {
    det::exp(a * det::ln(x))
}
fn close(a: f64, b: f64, scale: f64, tol: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tol * scale,
        "{a:e} vs {b:e}, scale {scale:e}, tolerance {tol:e}"
    );
}
fn stream(name: &str) -> fs_rand::Stream {
    // Test identity, not worker order, is the counter-stream key.
    let seed = name.bytes().fold(0xcbf29ce484222325u64, |s, b| {
        (s ^ u64::from(b)).wrapping_mul(0x100000001b3)
    });
    fs_rand::StreamKey {
        seed,
        kernel: 26,
        tile: 0,
    }
    .stream()
}

fn equal_laws(a: &Hyperelastic, b: &Hyperelastic, f: &[f64; 9], scale: f64, tol: f64) {
    close(a.energy(f), b.energy(f), scale, tol);
    for (x, y) in a.piola(f).unwrap().into_iter().zip(b.piola(f).unwrap()) {
        close(x, y, scale, tol);
    }
    for (x, y) in a
        .tangent(f)
        .unwrap()
        .into_iter()
        .flatten()
        .zip(b.tangent(f).unwrap().into_iter().flatten())
    {
        close(x, y, scale, tol);
    }
}

fn fd_tangent(law: &Hyperelastic, f: [f64; 9]) -> f64 {
    let a = law.tangent(&f).unwrap();
    let scale = a.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
    let mut worst = 0.0f64;
    for j in 0..9 {
        let (mut up, mut dn) = (f, f);
        up[j] += 1e-6;
        dn[j] -= 1e-6;
        let (p, m) = (law.piola(&up).unwrap(), law.piola(&dn).unwrap());
        for i in 0..9 {
            close(a[i][j], (p[i] - m[i]) / 2e-6, scale, 1e-6);
            worst = worst.max((a[i][j] - (p[i] - m[i]) / 2e-6).abs() / scale);
            close(a[i][j], a[j][i], scale, 1e-12);
        }
    }
    worst
}

#[test]
fn g0_ogden_one_term_alpha_two_is_decoupled_neo_hookean() {
    let (mu, kappa) = (2e5, 5e6);
    let law = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: mu,
                alpha: 2.0,
            }],
            volumetric: Volumetric::Quadratic { kappa },
        },
        4.0,
    )
    .unwrap();
    let nh = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: mu / 2.0,
            c01: 0.0,
            kappa,
        },
        4.0,
    )
    .unwrap();
    let mut rng = stream("g0_ogden_one_term_alpha_two_is_decoupled_neo_hookean");
    for _ in 0..1000 {
        let f = core::array::from_fn(|i| I[i] + 0.3 * (2.0 * rng.next_f64() - 1.0));
        assert!(determinant(&f) > 0.0);
        equal_laws(&law, &nh, &f, kappa, 1e-13);
    }
}

#[test]
fn g0_hgo_without_fibres_is_neo_hookean() {
    let f = [1.2, 0.1, 0., 0., 0.85, 0.05, 0.03, 0., 1.05];
    let nh = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: C10,
            c01: 0.,
            kappa: 5e6,
        },
        8.,
    )
    .unwrap();
    let empty = hgo(vec![], 0., FibreInvariants::Full, 5e6 / C10);
    equal_laws(&empty, &nh, &f, 5e6, 1e-14);
    let mut zero = hgo(
        vec![[1., 0., 0.]],
        0.2,
        FibreInvariants::Isochoric,
        5e6 / C10,
    );
    if let HyperelasticModel::HolzapfelGasserOgden { k1, .. } = &mut zero.model {
        *k1 = 0.0;
    }
    equal_laws(&zero, &nh, &f, 5e6, 1e-14);
    for invariants in [FibreInvariants::Full, FibreInvariants::Isochoric] {
        let a = hgo(vec![[1., 0., 0.]], 1. / 3., invariants, 5e6 / C10);
        let b = hgo(vec![[0., 1., 0.]], 1. / 3., invariants, 5e6 / C10);
        equal_laws(&a, &b, &f, 5e6, 1e-12);
    }
}

#[test]
fn g0_reference_state_and_coalescent_stretches() {
    let started = Instant::now();
    let law = ogden(20.0);
    let kappa = 20.0 * MU;
    let lambda = kappa - 2. * MU / 3.;
    // Reuse the elastic-law and Mandel tensor owners for the independent
    // small-strain limit. Contract dP/dF with symmetric Mandel directions.
    let young = 9. * kappa * MU / (3. * kappa + MU);
    let poisson = (3. * kappa - 2. * MU) / (2. * (3. * kappa + MU));
    let elastic = IsotropicElastic::new(young, poisson, 0.1)
        .unwrap()
        .stiffness();
    let a = law.tangent(&I).unwrap();
    let mut max_reference_error = 0.0f64;
    for p in law.piola(&I).unwrap() {
        close(p, 0., MU, 1e-14);
    }
    for i in 0..9 {
        for j in 0..9 {
            let (r, s, u, v) = (i / 3, i % 3, j / 3, j % 3);
            let delta = |a, b| f64::from(u8::from(a == b));
            let expected = lambda * delta(r, s) * delta(u, v)
                + MU * (delta(r, u) * delta(s, v) + delta(r, v) * delta(s, u));
            close(a[i][j], expected, kappa, 1e-10);
            max_reference_error = max_reference_error.max((a[i][j] - expected).abs() / kappa);
        }
    }
    let basis = mandel_basis();
    let mandel: [[f64; 6]; 6] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..9)
                .flat_map(|r| {
                    (0..9).map(move |s| basis[i][r / 3][r % 3] * a[r][s] * basis[j][s / 3][s % 3])
                })
                .sum()
        })
    });
    for i in 0..6 {
        for j in 0..6 {
            close(mandel[i][j], elastic[i][j], kappa, 1e-10);
        }
    }
    let rotated = rotate_mandel(mandel, rotation([0., 0., 1.], 0.73));
    for i in 0..6 {
        for j in 0..6 {
            close(rotated[i][j], mandel[i][j], kappa, 1e-10);
        }
    }
    let mut max_fd_error = 0.0f64;
    for f in [
        I,
        diag([1.4, power(1.4, -0.5), power(1.4, -0.5)]),
        diag([1.4, 1.4, power(1.4, -2.)]),
        diag([1.1; 3]),
    ] {
        max_fd_error = max_fd_error.max(fd_tangent(&law, f));
    }
    // The HGO energy is C1 at E=0; the declared inactive-side Hessian
    // there is the matrix tangent. FD checks use E strictly away from zero.
    let fibre = hgo(vec![[1., 0., 0.]], 0., FibreInvariants::Full, 20.0);
    for p in fibre.piola(&I).unwrap() {
        close(p, 0., C10, 1e-14);
    }
    assert!(
        fibre
            .tangent(&I)
            .unwrap()
            .iter()
            .flatten()
            .all(|v| v.is_finite())
    );
    fd_tangent(&fibre, diag([1.2, 1.1, 1.1]));
    fd_tangent(&fibre, diag([0.8, 1.1, 1.1]));
    println!(
        "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g0_reference_state_and_coalescent_stretches\",\"tier\":\"G0\",\"inputs\":{{\"parameter_set\":\"ogden1972-treloar\",\"citation\":\"10.1098/rspa.1972.0026\",\"kappa_over_mu\":20}},\"provenance\":{{\"source_ids\":[\"ogden1972\"],\"M8_W1_ledger_row\":\"ogden-1972-treloar-three-term-parameters\"}},\"rungs\":{{\"states\":4}},\"observed\":{{\"max_reference_tangent_relative_error\":{max_reference_error},\"max_fd_relative_error\":{max_fd_error}}},\"tolerance\":{{\"reference\":1e-10,\"fd\":1e-6}},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
        started.elapsed().as_secs_f64()
    );
}

#[test]
fn g0_frame_indifference_and_material_symmetry() {
    let mut rng = stream("g0_frame_indifference_and_material_symmetry");
    let laws = [
        ogden(20.),
        hgo(vec![[1., 0., 0.]], 0.1, FibreInvariants::Full, 20.),
    ];
    for law in &laws {
        for _ in 0..100 {
            let mut axis = core::array::from_fn(|_| 2. * rng.next_f64() - 1.);
            let norm = det::sqrt(axis.iter().map(|x| x * x).sum());
            for x in &mut axis {
                *x /= norm;
            }
            let q = rotation(axis, 6. * rng.next_f64() - 3.);
            let f = core::array::from_fn(|i| I[i] + 0.2 * (2. * rng.next_f64() - 1.));
            let qf = flat(matmul3(&q, &matrix(&f)));
            let p = law.piola(&f).unwrap();
            let expected = flat(matmul3(&q, &matrix(&p)));
            let scale = p.iter().map(|v| v.abs()).fold(MU, f64::max);
            for (a, b) in law.piola(&qf).unwrap().into_iter().zip(expected) {
                close(a, b, scale, 1e-12);
            }
            close(law.energy(&qf), law.energy(&f), scale, 1e-12);
        }
    }
    let f = diag([1.3, 0.9, 1.05]);
    let transverse = rotation([1., 0., 0.], 0.7);
    let arbitrary = rotation([0., 0., 1.], 0.7);
    for law in &laws {
        let fr = flat(matmul3(&matrix(&f), &transverse));
        close(law.energy(&fr), law.energy(&f), MU, 1e-12);
    }
    let fr = flat(matmul3(&matrix(&f), &arbitrary));
    assert!((laws[1].energy(&fr) - laws[1].energy(&f)).abs() > 1.0);
    close(laws[0].energy(&fr), laws[0].energy(&f), MU, 1e-12);
}

fn hgo_incompressible(lambda: f64, families: usize) -> f64 {
    let e = lambda * lambda - 1.;
    2. * C10 * (lambda * lambda - 1. / lambda)
        + if e > 0. {
            families as f64 * 2. * K1 * e * det::exp(K2 * e * e) * lambda * lambda
        } else {
            0.
        }
}

#[test]
fn g1_hgo_fibre_uniaxial_analytic() {
    let started = Instant::now();
    for invariants in [FibreInvariants::Isochoric, FibreInvariants::Full] {
        for families in [1, 2] {
            let law = hgo(vec![[1., 0., 0.]; families], 0., invariants, 1e5);
            let matrix_only = hgo(vec![], 0., invariants, 1e5);
            let mut worst = 0.0f64;
            for n in 0..=180 {
                let l = 0.7 + 0.9 * n as f64 / 180.;
                let t = power(l, -0.5);
                let f = diag([l, t, t]);
                let p = law.piola(&f).unwrap();
                // A Lagrange pressure enforces traction-free transverse
                // faces at prescribed J=1: sigma11-sigma22, not raw P11.
                let sigma = p[0] * l - p[4] * t;
                let oracle = hgo_incompressible(l, families);
                let error = if oracle == 0.0 {
                    close(sigma, 0.0, C10, 1e-14);
                    0.0
                } else {
                    (sigma - oracle).abs() / oracle.abs()
                };
                worst = worst.max(error);
                assert!(error < 1e-6, "l={l}, {invariants:?}, {sigma} vs {oracle}");
                if l <= 1. {
                    assert_eq!(law.energy(&f).to_bits(), matrix_only.energy(&f).to_bits());
                    assert_eq!(p, matrix_only.piola(&f).unwrap());
                }
            }
            println!(
                "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g1_hgo_fibre_uniaxial_analytic\",\"tier\":\"G1\",\"inputs\":{{\"parameter_set\":\"hgo2000-media\",\"citation\":\"10.1023/A:1010835316564\",\"kappa_over_c10\":1e5,\"families\":{families},\"invariants\":\"{invariants:?}\",\"kinematics\":\"J=1 with lateral-pressure projection\"}},\"provenance\":{{\"source_ids\":[\"hgo2000-fig14\"],\"M8_W1_ledger_row\":\"hgo-2000-fig14-media-parameters\"}},\"rungs\":{{\"stretch_count\":181,\"min\":0.7,\"max\":1.6}},\"observed\":{{\"max_relative_error\":{worst}}},\"tolerance\":1e-6,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
                started.elapsed().as_secs_f64()
            );
        }
    }
}

fn bisect(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let low = f(lo);
    let high = f(hi);
    assert!(
        low <= 0. && high >= 0.,
        "root not bracketed: {low:e}, {high:e}"
    );
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if mid == lo || mid == hi {
            break;
        }
        if f(mid) > 0. {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

// Independent diagonal compressible Ogden stress (no production AD/eigen).
fn ogden_diagonal(stretches: [f64; 3], kappa: f64) -> [f64; 3] {
    let j = stretches.iter().product::<f64>();
    let tau = stretches.map(|l| {
        TRELOAR
            .iter()
            .map(|p| p.mu_pa * power(l * power(j, -1. / 3.), p.alpha))
            .sum::<f64>()
    });
    let mean = tau.iter().sum::<f64>() / 3.;
    core::array::from_fn(|i| (tau[i] - mean + kappa * j * (j - 1.)) / stretches[i])
}

#[derive(Clone, Copy, Debug)]
enum Mode {
    Uniaxial,
    Equibiaxial,
    PureShear,
}
impl Mode {
    fn max(self) -> f64 {
        match self {
            Self::Uniaxial => 7.5,
            Self::Equibiaxial => 4.5,
            Self::PureShear => 5.,
        }
    }
    fn lateral(self) -> usize {
        if matches!(self, Self::Uniaxial) { 1 } else { 2 }
    }
    fn stretches(self, l: f64, j: f64) -> [f64; 3] {
        match self {
            Self::Uniaxial => [l, det::sqrt(j / l), det::sqrt(j / l)],
            Self::Equibiaxial => [l, l, j / (l * l)],
            Self::PureShear => [l, 1., j / l],
        }
    }
    fn oracle(self, l: f64) -> f64 {
        let exponent = match self {
            Self::Uniaxial => 0.5,
            Self::Equibiaxial => 2.,
            Self::PureShear => 1.,
        };
        TRELOAR
            .iter()
            .map(|p| p.mu_pa * (power(l, p.alpha - 1.) - power(l, -exponent * p.alpha - 1.)))
            .sum()
    }
}

#[test]
fn g2_treloar_ogden1972_reference_curves() {
    let started = Instant::now();
    let law = ogden(1e4);
    for mode in [Mode::Uniaxial, Mode::Equibiaxial, Mode::PureShear] {
        let mut worst = 0.0f64;
        for count in [33, 65, 129] {
            for n in 0..count {
                let l = 1. + (mode.max() - 1.) * n as f64 / (count - 1) as f64;
                let s = mode.stretches(l, 1.);
                let p = law.piola(&diag(s)).unwrap();
                let k = mode.lateral();
                let nominal = p[0] - p[3 * k + k] * s[k] / s[0];
                let expected = mode.oracle(l);
                let error = if expected == 0.0 {
                    close(nominal, 0.0, MU, 1e-14);
                    0.0
                } else {
                    (nominal - expected).abs() / expected.abs()
                };
                worst = worst.max(error);
                assert!(error < 1e-3, "{mode:?} l={l}: {nominal} vs {expected}");
            }
        }
        println!(
            "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g2_treloar_ogden1972_reference_curves\",\"tier\":\"G2\",\"inputs\":{{\"parameter_set\":\"ogden1972-treloar\",\"citation\":\"10.1098/rspa.1972.0026\",\"kappa_over_mu\":1e4,\"kinematics\":\"J=1 with lateral-pressure projection\"}},\"provenance\":{{\"source_ids\":[\"ogden1972\",\"treloar1944\"],\"M8_W1_ledger_row\":\"ogden-1972-treloar-three-term-parameters\"}},\"rungs\":{{\"stretch_counts\":[33,65,129],\"min\":1,\"max\":{}}},\"mode\":\"{mode:?}\",\"observed\":{{\"max_relative_error\":{worst}}},\"tolerance\":1e-3,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            mode.max(),
            started.elapsed().as_secs_f64()
        );
    }
    log_treloar_points(&law);
}

#[test]
fn g1_compressible_ogden_matches_diagonal_root_and_bulk_ladder() {
    for mode in [Mode::Uniaxial, Mode::Equibiaxial, Mode::PureShear] {
        let mut previous = None;
        for ratio in [1e4, 2e4, 4e4] {
            let law = ogden(ratio);
            let kappa = MU * ratio;
            let l = mode.max();
            let k = mode.lateral();
            let j = bisect(1., 1.1, |j| ogden_diagonal(mode.stretches(l, j), kappa)[k]);
            let stretches = mode.stretches(l, j);
            let p = law.piola(&diag(stretches)).unwrap();
            let exact = ogden_diagonal(stretches, kappa);
            for i in 0..3 {
                close(p[4 * i], exact[i], MU, 1e-10);
            }
            close(p[4 * k], 0., kappa, 1e-14);
            let error = (p[0] - mode.oracle(l)).abs() / mode.oracle(l).abs();
            if let Some(old) = previous {
                assert!(error / old > 0.49 && error / old < 0.51);
            }
            if ratio == 1e4 && !matches!(mode, Mode::PureShear) {
                // Counterexample to the literal finite-bulk interpretation
                // of the bead's 1e-3 band: no law or gate is changed to hide it.
                assert!(error > 1e-3);
            }
            previous = Some(error);
            println!(
                "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g1_compressible_ogden_matches_diagonal_root_and_bulk_ladder\",\"mode\":\"{mode:?}\",\"kappa_over_mu\":{ratio},\"J\":{j},\"incompressible_relative_residual\":{error},\"oracle\":\"traction-free compressible root\",\"verdict\":\"pass\"}}"
            );
        }
    }
}

#[test]
fn g1_hgo_finite_bulk_residual_has_derived_first_order_coefficient() {
    let l = 1.6;
    let d = 2. * K1 * (l * l - 1.) * det::exp(K2 * (l * l - 1.) * (l * l - 1.)) * l * l;
    let limit = hgo_incompressible(l, 1);
    // Full fibres depend only on prescribed lambda. Zero lateral stress:
    // kappa(J-1) = (2c10/3) J^(-5/3) (lambda²-J/lambda).
    // Hence J=1+A/kappa+O(kappa^-2), A=(2c10/3)(lambda²-1/lambda).
    // sigma = 2c10 J^(-5/3)(lambda²-J/lambda)+D/J, so its first
    // derivative at J=1 is B=2c10(-5lambda²/3+2/(3lambda))-D.
    let aa = 2. * C10 / 3. * (l * l - 1. / l);
    let bb = 2. * C10 * (-5. * l * l / 3. + 2. / (3. * l)) - d;
    let mut previous = None;
    for ratio in [1e5, 2e5, 4e5] {
        let kappa = C10 * ratio;
        let j = bisect(0.99, 1.01, |j| {
            kappa * (j - 1.) - 2. * C10 / 3. * power(j, -5. / 3.) * (l * l - j / l)
        });
        let t = det::sqrt(j / l);
        let law = hgo(vec![[1., 0., 0.]], 0., FibreInvariants::Full, ratio);
        let p = law.piola(&diag([l, t, t])).unwrap();
        let sigma = p[0] * l / j;
        let exact = 2. * C10 * power(j, -5. / 3.) * (l * l - j / l) + d / j;
        close(sigma, exact, limit, 1e-10);
        close(p[4], 0., kappa, 1e-14);
        let residual = sigma - limit;
        close(residual * kappa, aa * bb, (aa * bb).abs(), 1e-4);
        assert!(
            residual.abs() / limit > 1e-6,
            "the literal finite-bulk band is not satisfied"
        );
        if let Some(old) = previous {
            assert!(residual / old > 0.499 && residual / old < 0.501);
        }
        previous = Some(residual);
        println!(
            "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g1_hgo_finite_bulk_residual_has_derived_first_order_coefficient\",\"kappa_over_c10\":{ratio},\"J\":{j},\"observed_relative_residual\":{},\"literal_band\":1e-6,\"literal_band_satisfied\":false,\"derived_leading_relative_residual\":{},\"verdict\":\"pass\"}}",
            residual.abs() / limit,
            (aa * bb / kappa / limit).abs()
        );
    }
}

#[test]
fn g1_hgo_two_symmetric_families_match_independent_lateral_root() {
    let beta = 29. * std::f64::consts::PI / 180.;
    let (c, s) = (det::cos(beta), det::sin(beta));
    let law = hgo(
        vec![[c, s, 0.], [c, -s, 0.]],
        0.,
        FibreInvariants::Full,
        1e5,
    );
    for n in 1..=30 {
        let l = 1. + 0.6 * n as f64 / 30.;
        // Incompressible two-family uniaxial state diag(l,t,1/(l*t)).
        // sigma22-sigma33=0 fixes t in one dimension; no production stress
        // is used by the oracle. Symmetry cancels shear exactly.
        let diff = |t: f64| {
            let z = 1. / (l * t);
            let i4 = l * l * c * c + t * t * s * s;
            let e = i4 - 1.;
            let response = if e > 0. {
                4. * K1 * e * det::exp(K2 * e * e)
            } else {
                0.
            };
            2. * C10 * (t * t - z * z) + response * t * t * s * s
        };
        let t = bisect(0.2, 2., diff);
        let z = 1. / (l * t);
        close(diff(t), 0., C10, 1e-14);
        let e = l * l * c * c + t * t * s * s - 1.;
        let expected =
            2. * C10 * (l * l - z * z) + 4. * K1 * e * det::exp(K2 * e * e) * l * l * c * c;
        let p = law.piola(&diag([l, t, z])).unwrap();
        close(p[0] * l - p[8] * z, expected, expected.abs().max(C10), 1e-6);
        close(p[4] * t - p[8] * z, 0., C10, 1e-6);
        for i in [1, 2, 3, 5, 6, 7] {
            close(p[i], 0., C10, 1e-12);
        }
    }
}

fn log_treloar_points(law: &Hyperelastic) {
    let started = Instant::now();
    // Bucchi & Hearn (2013), Predictions of aneurysm formation in
    // distensible tubes, Part A, Appendix A, Table A1 (author manuscript):
    // https://researchportal.port.ac.uk/files/202484/Predictions_of_Aneurysm_Formation_in_Distensible_Tubes_PART_A%5B2%5D.pdf
    // These are published retrieved/digitized observations from Treloar
    // (1944), not a newly measured or refitted dataset. Uncertainty/fit band
    // and a redistribution licence are unstated: diagnostic logging only.
    // Original table gives engineering strain; lambda=1+strain.
    let uniaxial = [
        (0.144, 140110.),
        (0.247, 230196.),
        (0.410, 320368.),
        (0.615, 410600.),
        (0.901, 500946.),
        (1.186, 581299.),
        (1.431, 671589.),
        (2.023, 852311.),
        (2.595, 1052996.),
        (3.024, 1223493.),
        (3.761, 1594292.),
        (4.375, 1954917.),
        (4.765, 2305242.),
        (5.176, 2695562.),
        (5.404, 3045652.),
        (5.631, 3415724.),
        (5.879, 3775830.),
        (6.066, 4145853.),
        (6.172, 4495766.),
        (6.277, 4875658.),
        (6.444, 5235646.),
        (6.509, 5595496.),
    ];
    let equibiaxial = [
        (0.02, 92057.04),
        (0.06, 155979.),
        (0.11, 236293.47),
        (0.14, 257218.20),
        (0.20, 326084.40),
        (0.31, 434367.18),
        (0.42, 508452.30),
        (0.68, 647695.44),
        (0.94, 763159.14),
        (1.49, 959977.17),
        (2.03, 1239503.31),
        (2.43, 1440147.24),
        (2.75, 1706940.),
        (3.07, 1972378.98),
        (3.26, 2202364.62),
        (3.45, 2418459.30),
    ];
    let shear = [
        (0.03802, 72176.5845),
        (0.12388, 167213.412),
        (0.19062, 250829.928),
        (0.3053, 338332.185),
        (0.44878, 422092.908),
        (0.85081, 593629.587),
        (1.38715, 765421.326),
        (1.97142, 941098.806),
        (2.46932, 1124206.38),
        (2.95773, 1288308.06),
        (3.35973, 1467438.66),
        (3.68509, 1627449.57),
        (3.95282, 1802518.83),
    ];
    for (mode, points) in [
        (Mode::Uniaxial, uniaxial.as_slice()),
        (Mode::Equibiaxial, equibiaxial.as_slice()),
        (Mode::PureShear, shear.as_slice()),
    ] {
        for &(strain, observed) in points {
            let l = 1. + strain;
            let s = mode.stretches(l, 1.);
            let k = mode.lateral();
            let p = law.piola(&diag(s)).unwrap();
            let nominal = p[0] - p[4 * k] * s[k] / s[0];
            println!(
                "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W1a\",\"test\":\"g2_treloar_ogden1972_reference_curves\",\"tier\":\"G2\",\"mode\":\"{mode:?}\",\"inputs\":{{\"parameter_set\":\"ogden1972-treloar\",\"citation\":\"10.1098/rspa.1972.0026\",\"kappa_over_mu\":1e4}},\"provenance\":{{\"source_ids\":[\"Bucchi-Hearn-2013-Appendix-A-Treloar-1944\"],\"M8_W1_ledger_row\":\"treloar-1944-points-bucchi-hearn-2013\",\"uncertainty\":\"unstated\",\"redistribution_license\":\"unstated\"}},\"rungs\":{{\"stretch\":{l}}},\"observed\":{{\"nominal_pa\":{observed},\"computed_nominal_pa\":{nominal},\"relative_deviation\":{}}},\"tolerance\":null,\"verdict\":\"diagnostic-only-unpinned-source-band\",\"elapsed_s\":{}}}",
                (nominal - observed) / observed,
                started.elapsed().as_secs_f64()
            );
        }
    }
}

#[test]
fn ogden_term_cap_and_cap_plus_one() {
    for n in [MAX_OGDEN_TERMS, MAX_OGDEN_TERMS + 1] {
        let law = Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![
                    OgdenTerm {
                        mu_pa: 1.,
                        alpha: 2.
                    };
                    n
                ],
                volumetric: Volumetric::Quadratic { kappa: 10. },
            },
            3.,
        );
        if n == MAX_OGDEN_TERMS {
            assert!(law.unwrap().tangent(&I).is_ok());
        } else {
            assert!(matches!(law, Err(MaterialError::Parameters { .. })));
        }
    }
}
#[test]
fn fibre_family_cap_and_cap_plus_one() {
    for n in [MAX_FIBRE_FAMILIES, MAX_FIBRE_FAMILIES + 1] {
        let mut model = hgo(vec![], 0., FibreInvariants::Full, 10.).model;
        if let HyperelasticModel::HolzapfelGasserOgden { fibres, .. } = &mut model {
            *fibres = vec![[1., 0., 0.]; n];
        }
        let law = Hyperelastic::new(model, 3.);
        if n == MAX_FIBRE_FAMILIES {
            assert!(law.unwrap().tangent(&I).is_ok());
        } else {
            assert!(matches!(law, Err(MaterialError::Parameters { .. })));
        }
    }
}
#[test]
fn ogden_refuses_nonpositive_mu_alpha_product() {
    for (mu, alpha) in [
        (0., 2.),
        (1., 0.),
        (-1., 2.),
        (1., -2.),
        (f64::INFINITY, 2.),
        (1., f64::NAN),
        (1e308, 2.),
        (1e308, 1e-308),
    ] {
        assert!(matches!(
            Hyperelastic::new(
                HyperelasticModel::Ogden {
                    terms: vec![OgdenTerm { mu_pa: mu, alpha }],
                    volumetric: Volumetric::Quadratic { kappa: 10. },
                },
                3.
            ),
            Err(MaterialError::Parameters { .. })
        ));
    }
    assert!(ogden(10.).piola(&I).is_ok()); // published negative/negative term admitted
}
#[test]
fn hgo_refuses_non_unit_or_zero_fibre() {
    for a in [[0.; 3], [1. + 2e-12, 0., 0.], [f64::NAN, 0., 0.]] {
        let mut model = hgo(vec![], 0., FibreInvariants::Full, 10.).model;
        if let HyperelasticModel::HolzapfelGasserOgden { fibres, .. } = &mut model {
            *fibres = vec![a];
        }
        assert!(matches!(
            Hyperelastic::new(model, 3.),
            Err(MaterialError::Parameters { .. })
        ));
    }
}
#[test]
fn hgo_refuses_dispersion_outside_zero_to_one_third() {
    for dispersion in [-1e-12, 1. / 3. + 1e-12, f64::NAN] {
        let mut model = hgo(vec![], 0., FibreInvariants::Full, 10.).model;
        if let HyperelasticModel::HolzapfelGasserOgden { dispersion: d, .. } = &mut model {
            *d = dispersion;
        }
        assert!(matches!(
            Hyperelastic::new(model, 3.),
            Err(MaterialError::Parameters { .. })
        ));
    }
}
#[test]
fn principal_stretch_at_limit_accepted_beyond_refused() {
    for mut law in [
        ogden(10.),
        hgo(vec![[1., 0., 0.]], 0., FibreInvariants::Full, 10.),
    ] {
        law.stretch_limit = 3.;
        for l in [3., 1. / 3.] {
            let boundary = diag([l, 1., 1.]);
            assert!(law.piola(&boundary).is_ok());
            assert!(law.tangent(&boundary).is_ok());
            let beyond = diag([
                if l > 1. {
                    l * (1. + 1e-12)
                } else {
                    l / (1. + 1e-12)
                },
                1.,
                1.,
            ]);
            assert!(matches!(
                law.piola(&beyond),
                Err(MaterialError::State { .. })
            ));
            assert!(matches!(
                law.tangent(&beyond),
                Err(MaterialError::State { .. })
            ));
        }
        assert_eq!(
            law.principal_stretches(&diag([2., 0.5, 1.])).unwrap(),
            [0.5, 1., 2.]
        );
    }
}
#[test]
fn det_f_nonpositive_refused() {
    for law in [ogden(10.), hgo(vec![], 0., FibreInvariants::Full, 10.)] {
        for f in [
            diag([0., 1., 1.]),
            diag([-1., 1., 1.]),
            diag([f64::NAN, 1., 1.]),
        ] {
            assert!(matches!(law.piola(&f), Err(MaterialError::State { .. })));
            assert!(matches!(law.tangent(&f), Err(MaterialError::State { .. })));
        }
    }
}

#[test]
fn g0_bulk_and_stretch_parameter_domains_refuse() {
    for law in [ogden(10.), hgo(vec![], 0., FibreInvariants::Full, 10.)] {
        for limit in [0.0, 0.99, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                Hyperelastic::new(law.model.clone(), limit),
                Err(MaterialError::Parameters { .. })
            ));
        }
        for kappa in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            for bulk in [
                Volumetric::Quadratic { kappa },
                Volumetric::LogSquared { kappa },
            ] {
                let mut model = law.model.clone();
                match &mut model {
                    HyperelasticModel::Ogden { volumetric, .. }
                    | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => {
                        *volumetric = bulk
                    }
                    _ => unreachable!("new-law fixtures"),
                }
                assert!(matches!(
                    Hyperelastic::new(model, 3.),
                    Err(MaterialError::Parameters { .. })
                ));
            }
        }
    }
    assert!(matches!(
        Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![],
                volumetric: Volumetric::Quadratic { kappa: 10. },
            },
            3.
        ),
        Err(MaterialError::Parameters { .. })
    ));
}

/// G0/G3: a declared material domain is enforced by the original NH/MR
/// consumers as well as the spectral laws, independently of rigid rotations.
#[test]
fn g0_neo_hookean_and_mooney_rivlin_stretch_domains() {
    let q = rotation([0.0, 0.0, 1.0], 0.73);
    for model in [
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 5.0,
        },
        HyperelasticModel::MooneyRivlin {
            c10: 0.8,
            c01: 0.2,
            kappa: 7.0,
        },
    ] {
        let law = Hyperelastic::new(model, 3.0).unwrap();
        for stretch in [3.0, 1.0 / 3.0] {
            let at = diag([stretch, 1.0, 1.0]);
            let outside = diag([
                if stretch > 1.0 {
                    stretch * (1.0 + 1e-12)
                } else {
                    stretch / (1.0 + 1e-12)
                },
                1.0,
                1.0,
            ]);
            for f in [at, flat(matmul3(&q, &matrix(&at)))] {
                assert!(law.piola(&f).is_ok(), "inclusive limit: {f:?}");
                assert!(law.tangent(&f).is_ok(), "inclusive limit: {f:?}");
            }
            for f in [outside, flat(matmul3(&q, &matrix(&outside)))] {
                assert!(matches!(law.piola(&f), Err(MaterialError::State { .. })));
                assert!(matches!(law.tangent(&f), Err(MaterialError::State { .. })));
            }
        }
        // Positive determinant alone is insufficient: an isochoric stretch
        // can leave BOTH the tensile and compressive material domains.
        let isochoric = diag([50.0, 1.0 / 50.0, 1.0]);
        assert!(matches!(
            law.piola(&isochoric),
            Err(MaterialError::State { .. })
        ));
        assert!(matches!(
            law.tangent(&isochoric),
            Err(MaterialError::State { .. })
        ));
    }
}

#[test]
fn g0_neo_hookean_and_mooney_rivlin_invalid_stretch_limits_refuse() {
    for model in [
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 5.0,
        },
        HyperelasticModel::MooneyRivlin {
            c10: 0.8,
            c01: 0.2,
            kappa: 7.0,
        },
    ] {
        let identity_only = Hyperelastic::new(model.clone(), 1.0).unwrap();
        assert!(identity_only.piola(&I).is_ok());
        assert!(identity_only.tangent(&I).is_ok());
        for limit in [-1.0, 0.0, 0.999999, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                Hyperelastic::new(model.clone(), limit),
                Err(MaterialError::Parameters { .. })
            ));
            let mut changed = identity_only.clone();
            changed.stretch_limit = limit;
            assert!(matches!(
                changed.piola(&I),
                Err(MaterialError::Parameters { .. })
            ));
            assert!(matches!(
                changed.tangent(&I),
                Err(MaterialError::Parameters { .. })
            ));
        }
    }
}

#[test]
fn g0_log_bulk_split_and_parameter_derived_admissibility() {
    for volumetric in [
        Volumetric::Quadratic { kappa: 5e6 },
        Volumetric::LogSquared { kappa: 5e6 },
    ] {
        let law = Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 1e5,
                    alpha: 2.,
                }],
                volumetric,
            },
            3.,
        )
        .unwrap();
        let f = diag([1.2, 0.9, 1.1]);
        close(
            law.energy(&f),
            law.isochoric_energy(&f) + law.volumetric_energy(determinant(&f)),
            5e6,
            1e-14,
        );
        fd_tangent(&law, f);
        let expected = matches!(volumetric, Volumetric::Quadratic { .. }).then_some(true);
        assert_eq!(law.admissibility().polyconvex, expected);
    }
    assert_eq!(ogden(10.).admissibility().polyconvex, None);
    assert_eq!(
        hgo(vec![[1., 0., 0.]], 0.1, FibreInvariants::Full, 10.)
            .admissibility()
            .polyconvex,
        Some(true)
    );
    assert_eq!(
        hgo(vec![[1., 0., 0.]], 0.1, FibreInvariants::Isochoric, 10.)
            .admissibility()
            .polyconvex,
        None
    );
    // Hydrostatic expansion exposes the documented full/isochoric difference.
    let full = hgo(vec![[1., 0., 0.]], 0., FibreInvariants::Full, 10.);
    let iso = hgo(vec![[1., 0., 0.]], 0., FibreInvariants::Isochoric, 10.);
    let f = diag([1.1; 3]);
    assert!(full.isochoric_energy(&f) > 1.);
    close(iso.isochoric_energy(&f), 0., C10, 1e-14);
}

fn model_card(hgo_model: bool) -> fs_matdb::ConstitutiveModelCard {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use fs_qty::Dims;
    let one = Dims([0, 0, 0, 0, 0, 0]);
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let mut parameters = std::collections::BTreeMap::new();
    for (name, value, dims) in [
        ("stretch_limit", 8., one),
        ("volumetric_kind", 0., one),
        ("kappa", 5e6, pa),
    ] {
        parameters.insert(name.to_string(), LawParameter { value, dims });
    }
    if hgo_model {
        for (name, value, dims) in [
            ("c10", C10, pa),
            ("k1", K1, pa),
            ("k2", K2, one),
            ("dispersion", 0., one),
            ("fibre_invariants", 1., one),
            ("fibre_count", 1., one),
            ("fibre_0_x", 1., one),
            ("fibre_0_y", 0., one),
            ("fibre_0_z", 0., one),
        ] {
            parameters.insert(name.to_string(), LawParameter { value, dims });
        }
    } else {
        parameters.insert(
            "term_count".to_string(),
            LawParameter {
                value: 3.,
                dims: one,
            },
        );
        for (i, p) in TRELOAR.iter().enumerate() {
            parameters.insert(
                format!("mu_{i}"),
                LawParameter {
                    value: p.mu_pa,
                    dims: pa,
                },
            );
            parameters.insert(
                format!("alpha_{i}"),
                LawParameter {
                    value: p.alpha,
                    dims: one,
                },
            );
        }
    }
    let source = fs_blake3::hash_bytes(
        b"authored software transport fixture; numerical literature coefficients only",
    );
    ConstitutiveModelCard {
        law: LawId(
            if hgo_model {
                "hgo-dispersed"
            } else {
                "ogden-decoupled"
            }
            .to_string(),
        ),
        law_version: 1,
        parameters,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: fs_evidence::ValidityDomain::unconstrained(),
        sources: vec![source],
        provenance: Provenance {
            source: "authored law-card transport test".to_string(),
            license: "test fixture".to_string(),
            artifact: Some(source),
        },
    }
}

#[test]
fn g0_dimensioned_model_pack_reaches_real_hyperelastic_registry() {
    use fs_matdb::{
        MODEL_PACK_TARGET_BASIS, ModelNormalizationReceipt, ModelNormalizationTarget,
        NormalizedModelPack,
    };
    use fs_material::graph::{
        LawRegistry, check_consistent_tangent, check_free_energy_consistency,
    };
    use fs_material::hyper::cards::register_hyperelastic_laws;
    for is_hgo in [false, true] {
        let card = model_card(is_hgo);
        let source = card.provenance.artifact.unwrap();
        let receipts = card
            .parameters
            .iter()
            .map(|(name, p)| {
                ModelNormalizationReceipt::new(
                    ModelNormalizationTarget::Parameter {
                        model: card.content_hash(),
                        parameter: name.clone(),
                    },
                    fs_blake3::hash_bytes(format!("{name}={}", p.value).as_bytes()),
                    p.dims,
                    1.,
                    0.,
                    "SI",
                    MODEL_PACK_TARGET_BASIS,
                    None,
                    None,
                )
            })
            .collect();
        let pack = NormalizedModelPack::new(
            "hyperelastic-fixture",
            "test-v1",
            source,
            "authored fixture",
            vec![card.clone()],
            receipts,
        )
        .unwrap();
        let loaded =
            NormalizedModelPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes())
                .unwrap();
        let transported = &loaded.models()[0];
        assert_eq!(transported, &card);
        let material = Hyperelastic::from_model_card(transported).unwrap();
        let mut registry = LawRegistry::new();
        register_hyperelastic_laws(&mut registry);
        let node = registry.instantiate("solid/material", transported).unwrap();
        let f = [1.2, 0.1, 0., 0., 0.9, 0.04, 0., 0., 1.05];
        let h: Vec<_> = f.iter().zip(I).map(|(a, b)| a - b).collect();
        let output = node.evaluate(&[], &h).unwrap();
        for (a, b) in output.outputs.iter().zip(material.piola(&f).unwrap()) {
            close(*a, b, 5e6, 1e-13);
        }
        check_consistent_tangent("solid/material", node.as_ref(), &[], &h, 1e-3).unwrap();
        check_free_energy_consistency("solid/material", node.as_ref(), &[], &h, 1e-3).unwrap();
        assert!(node.evaluate(&[], &h[..8]).is_err());
        assert!(node.evaluate(&[0.], &h).is_err());
        assert!(node.evaluate(&[], &[-2.; 9]).is_err());
        let mut wrong = card.clone();
        wrong.parameters.get_mut("kappa").unwrap().dims = fs_qty::Dims([0; 6]);
        assert!(Hyperelastic::from_model_card(&wrong).is_err());
        wrong = card.clone();
        wrong.parameters.remove("volumetric_kind");
        assert!(Hyperelastic::from_model_card(&wrong).is_err());
        wrong = card.clone();
        wrong.parameters.insert(
            "unexpected".to_string(),
            fs_matdb::LawParameter {
                value: 1.,
                dims: fs_qty::Dims([0; 6]),
            },
        );
        assert!(Hyperelastic::from_model_card(&wrong).is_err());
        wrong = card.clone();
        wrong.law_version = 2;
        assert!(Hyperelastic::from_model_card(&wrong).is_err());
        wrong = card.clone();
        wrong
            .parameters
            .get_mut(if is_hgo { "fibre_count" } else { "term_count" })
            .unwrap()
            .value = 1.5;
        assert!(Hyperelastic::from_model_card(&wrong).is_err());
    }
}
