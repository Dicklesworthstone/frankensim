//! M0.3-W3: the prescribed active fibre stress on top of the passive Ogden/HGO laws. Analytic oracles are written
//! out here and share no code with the law: the isometric Cauchy stress, the traction-free free-contraction root
//! at prescribed J = 1 (a Lagrange pressure removes the lateral stress), objectivity, the finite-difference tangent,
//! and the active-work ledger. No parameter is fitted; the force-length families are test fixtures.

use fs_material::active_tissue::{
    ActiveFibreStress, ActiveTissue, ForceLength, MAX_FORCE_LENGTH_KNOTS,
};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, Volumetric};
use fs_material::tensor::{matmul3, rotation};
use fs_material::MaterialError;
use fs_math::det;
use std::time::Instant;

const I: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
const T_MAX: f64 = 1.0e5;
const MU: f64 = 2.0e5;
const DOMAIN: (f64, f64) = (0.5, 1.5);

fn families() -> Vec<ForceLength> {
    vec![
        ForceLength::Gaussian {
            stretch_opt: 1.0,
            width: 0.4,
        },
        ForceLength::Parabola {
            stretch_opt: 1.0,
            half_width: 0.6,
        },
        ForceLength::Table {
            knots: vec![(0.5, 0.0), (0.8, 0.7), (1.0, 1.0), (1.2, 0.6), (1.5, 0.0)],
        },
    ]
}

// Closed forms of the three fixtures, written out here (no call into the law): family 0 Gaussian, 1 parabola, 2 table.
fn closed_form(family: usize, l: f64) -> f64 {
    match family {
        0 => det::exp(-((l - 1.0) / 0.4) * ((l - 1.0) / 0.4)),
        1 => (1.0 - ((l - 1.0) / 0.6) * ((l - 1.0) / 0.6)).max(0.0),
        _ => {
            if l <= 0.5 || l >= 1.5 {
                0.0
            } else if l < 0.8 {
                0.7 * (l - 0.5) / 0.3
            } else if l < 1.0 {
                0.7 + 0.3 * (l - 0.8) / 0.2
            } else if l < 1.2 {
                1.0 - 0.4 * (l - 1.0) / 0.2
            } else {
                0.6 - 0.6 * (l - 1.2) / 0.3
            }
        }
    }
}

fn active(force_length: ForceLength) -> ActiveFibreStress {
    ActiveFibreStress::new(T_MAX, force_length, [1., 0., 0.], DOMAIN).unwrap()
}

fn neo_hookean() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: MU,
            lambda: 5.0e6,
        },
        4.0,
    )
    .unwrap()
}

fn hgo(fibres: Vec<[f64; 3]>) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::HolzapfelGasserOgden {
            c10: 1.5e3,
            k1: 2.3632e3,
            k2: 0.8393,
            dispersion: 0.0,
            fibres,
            invariants: FibreInvariants::Full,
            volumetric: Volumetric::Quadratic { kappa: 1.5e8 },
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

fn close(a: f64, b: f64, scale: f64, tol: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tol * scale,
        "{a:e} vs {b:e}, scale {scale:e}, tolerance {tol:e}"
    );
}

fn stream(name: &str) -> fs_rand::Stream {
    let seed = name.bytes().fold(0xcbf29ce484222325u64, |s, b| {
        (s ^ u64::from(b)).wrapping_mul(0x100000001b3)
    });
    fs_rand::StreamKey {
        seed,
        kernel: 28,
        tile: 0,
    }
    .stream()
}

fn bisect(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lo) <= 0. && f(hi) >= 0., "root not bracketed");
    for _ in 0..200 {
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

#[test]
fn g0_zero_activation_is_passive_law_bitwise() {
    let law = ActiveTissue::new(hgo(vec![[1., 0., 0.]]), active(families().remove(0))).unwrap();
    let mut rng = stream("g0_zero_activation_is_passive_law_bitwise");
    for _ in 0..1000 {
        let f: [f64; 9] = core::array::from_fn(|i| I[i] + 0.2 * (2. * rng.next_f64() - 1.));
        assert!(determinant(&f) > 0.0);
        let p = law.piola(&f, 0.0).unwrap();
        let q = law.passive.piola(&f).unwrap();
        assert!(p.iter().zip(&q).all(|(a, b)| a.to_bits() == b.to_bits()));
        let a = law.tangent(&f, 0.0).unwrap();
        let b = law.passive.tangent(&f).unwrap();
        assert!(
            a.iter()
                .flatten()
                .zip(b.iter().flatten())
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
    }
}

#[test]
fn g1_isometric_active_stress_analytic() {
    let started = Instant::now();
    let mut worst = 0.0f64;
    for (n, family) in families().into_iter().enumerate() {
        let law = ActiveTissue::new(hgo(vec![[1., 0., 0.]]), active(family.clone())).unwrap();
        for step in 0..50 {
            let l = 0.6 + 0.8 * step as f64 / 49.;
            let t = 0.9;
            let f = diag([l, t, t]);
            let j = l * t * t;
            for activation in [0.25, 0.5, 1.0] {
                let p = law.active_piola(&f, activation).unwrap();
                // Cauchy stress sigma = P F^T / J for a diagonal F.
                let sigma = [p[0] * l / j, p[4] * t / j, p[8] * t / j];
                let expected = T_MAX * activation * closed_form(n, l);
                close(family.factor(l), closed_form(n, l), 1.0, 1e-13);
                close(sigma[0], expected, T_MAX, 1e-12);
                worst = worst.max((sigma[0] - expected).abs() / T_MAX);
                // The orthogonal components are exactly zero.
                assert!(p.iter().enumerate().all(|(i, v)| i == 0 || *v == 0.0));
                assert!(sigma[1] == 0.0 && sigma[2] == 0.0);
            }
        }
        println!(
            "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W3\",\"test\":\"g1_isometric_active_stress_analytic\",\"tier\":\"G1\",\"inputs\":{{\"passive\":\"hgo-fixture\",\"active\":\"force-length-family-{n}\",\"activation\":[0.25,0.5,1.0]}},\"provenance\":{{\"dataset\":\"none (analytic fixture)\",\"M8_W1_ledger_row\":null}},\"rungs\":{{\"stretch_count\":50,\"min\":0.6,\"max\":1.4}},\"observed\":{{\"max_relative_error\":{worst}}},\"tolerance\":1e-12,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            started.elapsed().as_secs_f64()
        );
    }
}

#[test]
fn g1_free_contraction_equilibrium() {
    let started = Instant::now();
    let passive = neo_hookean();
    for activation in [0.25, 0.5, 1.0] {
        let family = ForceLength::Gaussian {
            stretch_opt: 1.0,
            width: 0.4,
        };
        let law = ActiveTissue::new(passive.clone(), active(family.clone())).unwrap();
        // Prescribed J = 1 with the lateral pressure eliminated: the traction-free block balances
        // sigma11 - sigma22 = 0. For the compressible Neo-Hookean at J = 1, sigma11 - sigma22 = mu (l^2 - 1/l);
        // the active Cauchy stress acts along the fibre only. This oracle uses no production code.
        let analytic = |l: f64| MU * (l * l - 1. / l) + T_MAX * activation * family.factor(l);
        let produced = |l: f64| {
            let t = det::exp(-0.5 * det::ln(l));
            let f = diag([l, t, t]);
            let p = law.piola(&f, activation).unwrap();
            let j = determinant(&f);
            p[0] * l / j - p[4] * t / j
        };
        let root_analytic = bisect(0.5, 1.0, analytic);
        let root_produced = bisect(0.5, 1.0, produced);
        close(root_produced, root_analytic, 1.0, 1e-10);
        assert!(root_analytic < 1.0, "an activated block shortens");
        println!(
            "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W3\",\"test\":\"g1_free_contraction_equilibrium\",\"tier\":\"G1\",\"inputs\":{{\"passive\":\"neo-hookean-fixture\",\"active\":\"gaussian-fixture\",\"activation\":{activation}}},\"provenance\":{{\"dataset\":\"none (analytic fixture)\",\"M8_W1_ledger_row\":null}},\"rungs\":{{\"bracket\":[0.5,1.0]}},\"observed\":{{\"root\":{root_produced},\"root_difference\":{}}},\"tolerance\":1e-10,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            (root_produced - root_analytic).abs(),
            started.elapsed().as_secs_f64()
        );
    }
}

#[test]
fn g0_active_power_ledger() {
    // Over a prescribed stretch history with time-varying activation, W_total = Delta Psi_passive + W_active.
    let law = ActiveTissue::new(
        hgo(vec![[1., 0., 0.]]),
        active(ForceLength::Gaussian {
            stretch_opt: 1.0,
            width: 0.4,
        }),
    )
    .unwrap();
    let g = [0.9, 0.2, -0.1, 0.1, -0.7, 0.15, -0.2, 0.05, 0.4];
    let f_at = |t: f64| -> [f64; 9] { core::array::from_fn(|i| I[i] + 0.3 * t * g[i]) };
    let f_dot: [f64; 9] = core::array::from_fn(|i| 0.3 * g[i]);
    let activation_at = |t: f64| 3. * t * t - 2. * t * t * t;
    // Composite 5-point Gauss-Legendre over 200 subintervals: exact to rounding for this smooth integrand.
    let nodes = [
        0.0,
        0.538_469_310_105_683_1,
        -0.538_469_310_105_683_1,
        0.906_179_845_938_664,
        -0.906_179_845_938_664,
    ];
    let weights = [
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_5,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
        0.236_926_885_056_189_1,
    ];
    let (mut w_total, mut w_active) = (0.0f64, 0.0f64);
    let pieces = 200;
    for k in 0..pieces {
        let (a, b) = (k as f64 / pieces as f64, (k + 1) as f64 / pieces as f64);
        for (x, w) in nodes.iter().zip(weights) {
            let t = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let (f, act) = (f_at(t), activation_at(t));
            let p = law.piola(&f, act).unwrap();
            let total: f64 = p.iter().zip(&f_dot).map(|(u, v)| u * v).sum();
            w_total += 0.5 * (b - a) * w * total;
            w_active += 0.5 * (b - a) * w * law.active_power(&f, &f_dot, act).unwrap();
        }
    }
    let delta_psi = law.passive.energy(&f_at(1.0)) - law.passive.energy(&f_at(0.0));
    let scale = w_total.abs().max(delta_psi.abs()).max(w_active.abs());
    assert!(w_active.abs() > 1e-3 * scale, "the active work must be a visible part of the ledger");
    close(w_total, delta_psi + w_active, scale, 1e-10);
    println!(
        "{{\"suite\":\"irp-s0-materials\",\"bead\":\"M0.3-W3\",\"test\":\"g0_active_power_ledger\",\"tier\":\"G0\",\"observed\":{{\"w_total\":{w_total},\"delta_psi\":{delta_psi},\"w_active\":{w_active},\"closure\":{}}},\"tolerance\":1e-10,\"verdict\":\"pass\"}}",
        (w_total - delta_psi - w_active).abs() / scale
    );
}

#[test]
fn g0_objectivity_and_the_active_stress_follows_the_rotated_fibre() {
    let law = ActiveTissue::new(
        hgo(vec![[1., 0., 0.]]),
        active(ForceLength::Parabola {
            stretch_opt: 1.0,
            half_width: 0.6,
        }),
    )
    .unwrap();
    let mut rng = stream("g0_objectivity");
    for _ in 0..100 {
        let mut axis: [f64; 3] = core::array::from_fn(|_| 2. * rng.next_f64() - 1.);
        let norm = det::sqrt(axis.iter().map(|x| x * x).sum());
        for x in &mut axis {
            *x /= norm;
        }
        let q = rotation(axis, 6. * rng.next_f64() - 3.);
        let f: [f64; 9] = core::array::from_fn(|i| I[i] + 0.15 * (2. * rng.next_f64() - 1.));
        let qf = flat(matmul3(&q, &matrix(&f)));
        let activation = rng.next_f64();
        let p = law.piola(&f, activation).unwrap();
        let expected = flat(matmul3(&q, &matrix(&p)));
        // The floor is the passive bulk modulus: J rounds at 1e-16 and the volumetric stress is kappa (J - 1).
        let scale = p.iter().map(|v| v.abs()).fold(1.0e8, f64::max);
        for (a, b) in law.piola(&qf, activation).unwrap().into_iter().zip(expected) {
            close(a, b, scale, 1e-12);
        }
    }
}

#[test]
fn g0_tangent_matches_fd_and_is_declared_non_symmetric() {
    for family in families() {
        let law = ActiveTissue::new(hgo(vec![[1., 0., 0.]]), active(family)).unwrap();
        // The fibre stretch stays away from 1, where the passive HGO fibre switches on.
        let f = [0.92, 0.05, -0.02, 0.03, 1.04, 0.06, -0.04, 0.02, 0.97];
        let activation = 0.6;
        let total = law.tangent(&f, activation).unwrap();
        let passive = law.passive.tangent(&f).unwrap();
        let part: Vec<Vec<f64>> = (0..9)
            .map(|i| (0..9).map(|j| total[i][j] - passive[i][j]).collect())
            .collect();
        let active_scale = part.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        let total_scale = total.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(active_scale > 1e3, "the active tangent is a visible part");
        for j in 0..9 {
            let (mut up, mut dn) = (f, f);
            up[j] += 1e-6;
            dn[j] -= 1e-6;
            let (pa, ma) = (
                law.active_piola(&up, activation).unwrap(),
                law.active_piola(&dn, activation).unwrap(),
            );
            let (pt, mt) = (
                law.piola(&up, activation).unwrap(),
                law.piola(&dn, activation).unwrap(),
            );
            for i in 0..9 {
                // The active stress alone, against its own scale, then the total against the total scale.
                close(part[i][j], (pa[i] - ma[i]) / 2e-6, active_scale, 1e-6);
                close(total[i][j], (pt[i] - mt[i]) / 2e-6, total_scale, 1e-6);
            }
        }
        let asymmetry = (0..9)
            .flat_map(|i| (0..9).map(move |j| (i, j)))
            .map(|(i, j)| (part[i][j] - part[j][i]).abs())
            .fold(0.0, f64::max);
        assert!(asymmetry > 1e-3 * active_scale, "the active tangent is not symmetric");
        assert!(!law.admissibility().tangent_symmetric);
    }
}

#[test]
fn refuses_activation_outside_unit_interval() {
    let law = ActiveTissue::new(neo_hookean(), active(families().remove(0))).unwrap();
    let f = diag([1.0, 1.0, 1.0]);
    assert!(law.piola(&f, 1.0).is_ok());
    assert!(law.piola(&f, 0.0).is_ok());
    for bad in [1.0 + 1e-12, -1e-12, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            law.piola(&f, bad),
            Err(MaterialError::State { .. })
        ));
        assert!(matches!(
            law.tangent(&f, bad),
            Err(MaterialError::State { .. })
        ));
    }
}

#[test]
fn refuses_fibre_stretch_outside_force_length_domain() {
    let law = ActiveTissue::new(neo_hookean(), active(families().remove(0))).unwrap();
    let at = |l: f64| diag([l, 1.0, 1.0]);
    for bound in [DOMAIN.0, DOMAIN.1] {
        assert!(law.piola(&at(bound), 0.5).is_ok(), "the bound itself is accepted");
    }
    assert!(matches!(
        law.piola(&at(DOMAIN.1 * (1.0 + 1e-12)), 0.5),
        Err(MaterialError::State { .. })
    ));
    assert!(matches!(
        law.piola(&at(DOMAIN.0 * (1.0 - 1e-12)), 0.5),
        Err(MaterialError::State { .. })
    ));
    // At zero activation no active stress is evaluated, so the domain is not consulted.
    assert!(law.piola(&at(DOMAIN.1 * 1.01), 0.0).is_ok());
}

#[test]
fn refuses_non_unit_or_mismatched_fibre() {
    let family = || families().remove(0);
    for fibre in [[1., 1., 0.], [0., 0., 0.], [1.0 + 1e-9, 0., 0.], [f64::NAN, 0., 0.]] {
        assert!(matches!(
            ActiveFibreStress::new(T_MAX, family(), fibre, DOMAIN),
            Err(MaterialError::Parameters { .. })
        ));
    }
    // The active fibre must be one of the passive HGO card's fibres.
    let transverse = ActiveFibreStress::new(T_MAX, family(), [0., 1., 0.], DOMAIN).unwrap();
    assert!(matches!(
        ActiveTissue::new(hgo(vec![[1., 0., 0.]]), transverse.clone()),
        Err(MaterialError::Parameters { .. })
    ));
    assert!(ActiveTissue::new(hgo(vec![[1., 0., 0.], [0., 1., 0.]]), transverse.clone()).is_ok());
    assert!(ActiveTissue::new(neo_hookean(), transverse).is_ok());
}

#[test]
fn refuses_nonpositive_t_max_and_invalid_force_length() {
    let family = || families().remove(0);
    for t_max in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ActiveFibreStress::new(t_max, family(), [1., 0., 0.], DOMAIN),
            Err(MaterialError::Parameters { .. })
        ));
    }
    for domain in [(0.0, 1.5), (1.5, 0.5), (1.0, 1.0), (0.5, f64::INFINITY)] {
        assert!(matches!(
            ActiveFibreStress::new(T_MAX, family(), [1., 0., 0.], domain),
            Err(MaterialError::Parameters { .. })
        ));
    }
    let knots = |n: usize| ForceLength::Table {
        knots: (0..n).map(|i| (0.5 + i as f64 * 0.01, 0.5)).collect(),
    };
    assert!(ActiveFibreStress::new(T_MAX, knots(MAX_FORCE_LENGTH_KNOTS), [1., 0., 0.], DOMAIN).is_ok());
    assert!(matches!(
        ActiveFibreStress::new(T_MAX, knots(MAX_FORCE_LENGTH_KNOTS + 1), [1., 0., 0.], DOMAIN),
        Err(MaterialError::Parameters { .. })
    ));
    for bad in [
        ForceLength::Gaussian { stretch_opt: 1.0, width: 0.0 },
        ForceLength::Gaussian { stretch_opt: -1.0, width: 0.4 },
        ForceLength::Parabola { stretch_opt: 1.0, half_width: f64::NAN },
        ForceLength::Table { knots: vec![(1.0, 1.0)] },
        ForceLength::Table { knots: vec![(1.0, 1.0), (1.0, 0.5)] },
        ForceLength::Table { knots: vec![(0.8, 1.0), (1.0, 1.2)] },
    ] {
        assert!(matches!(
            ActiveFibreStress::new(T_MAX, bad, [1., 0., 0.], DOMAIN),
            Err(MaterialError::Parameters { .. })
        ));
    }
}
