//! fs-nlmodal's joint loaded modal basis, `loaded::LoadedModalBasis`, against cases solved by
//! hand: one mode under a known added mass, two equal modes coupled by a mutual mass, two unequal
//! ones by the quadratic of their pencil; the basis's mass normalization, the attachments carried
//! with it, the structure it builds, and its refusals.
//! Run with `cargo test -p fs-nlmodal --test loaded_basis`.

// Exact float comparisons here assert values carried unchanged.
#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_nlmodal::component::{ModalAcousticTimeBudget, ModalAttachment, ModalStructureRefusal};
use fs_nlmodal::loaded::LoadedModalBasis;
use fs_ports::component::{Component, PortRole};
use fs_ports::{PortKind, StableId};

const TAU: f64 = core::f64::consts::TAU;
const EPS: f64 = f64::EPSILON;

#[test]
fn one_mode_under_added_mass_mu_sounds_at_f0_over_root_one_plus_mu() {
    for (f0, mu) in [(329.63, 0.0), (329.63, 0.012), (1_318.5, 0.3), (88.0, 2.5)] {
        let w0 = TAU * f0;
        let basis = LoadedModalBasis::new(&[w0], &[mu]).expect("definite");
        let w = basis.angular_frequencies_rad_s()[0];
        let exact = w0 / (1.0 + mu).sqrt();
        // The reduction divides w0^2 twice by sqrt(1 + mu) and each side takes a square root:
        // five roundings, each of at most eps / 2 of the value.
        assert!(
            (w - exact).abs() <= 4.0 * EPS * exact,
            "f0 {f0} mu {mu}: {w} against {exact}"
        );
        // Phi = 1 / sqrt(1 + mu), normalized in M + M_a.
        let phi = basis.basis()[0];
        assert!(
            (phi * phi * (1.0 + mu) - 1.0).abs() <= 8.0 * EPS,
            "phi {phi}"
        );
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 11,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// `w'^2` against `lambda` within fs-modal's certified residual, widened by what rounds outside
/// it. The residual is itself computed in floating point: each of its `n` entries is a sum of
/// `2 n` products of entries up to the stiffness scale `largest` (with `M + M_a` near the identity
/// here, its inverse norm is about the Euclidean one), so it can understate the true residual by
/// up to `(2 n + 2) eps` of twice that scale. `w'` is a rounded square root, squared again here,
/// and the formula by hand rounds a few times: eight of `lambda` for both.
fn within_the_certificate(basis: &LoadedModalBasis, j: usize, lambda: f64, largest: f64) {
    let w = basis.angular_frequencies_rad_s()[j];
    let n = basis.modes() as f64;
    let bound = basis.residuals()[j] + 8.0 * EPS * lambda + 2.0 * (2.0 * n + 2.0) * EPS * largest;
    assert!(
        (w * w - lambda).abs() <= bound,
        "mode {j}: {} against {lambda}, bound {bound}",
        w * w
    );
}

#[test]
fn two_equal_modes_split_by_their_mutual_mass() {
    let (w0, mu, m) = (TAU * 329.63, 0.05, 0.02);
    let basis = LoadedModalBasis::new(&[w0, w0], &[mu, m, m, mu]).expect("definite");
    // (1, 1) carries mass 1 + mu + m, (1, -1) carries 1 + mu - m: the in-phase pair is lower.
    within_the_certificate(&basis, 0, w0 * w0 / (1.0 + mu + m), w0 * w0);
    within_the_certificate(&basis, 1, w0 * w0 / (1.0 + mu - m), w0 * w0);
    let phi = basis.basis();
    // Column 0 is (1, 1) / sqrt(2 (1 + mu + m)) up to sign, column 1 is (1, -1) / sqrt(2 (1 + mu - m)).
    let (a, b) = (
        (2.0 * (1.0 + mu + m)).sqrt().recip(),
        (2.0 * (1.0 + mu - m)).sqrt().recip(),
    );
    // Each computed vector is within an angle of about eta lambda_1 / (lambda_1 - lambda_0) of
    // the exact one (Davis and Kahan), with eta the backward error of the Cholesky reduction and
    // the Jacobi solve, taken as 8 eps for two modes; lambda_1 / (lambda_1 - lambda_0) is
    // (1 + mu + m) / (2 m). Four roundings of the normalization add 4 eps.
    let tolerance = (8.0 * (1.0 + mu + m) / (2.0 * m) + 4.0) * EPS;
    assert!(
        (phi[0].abs() - a).abs() <= tolerance * a && (phi[2] - phi[0]).abs() <= 2.0 * tolerance * a
    );
    assert!(
        (phi[1].abs() - b).abs() <= tolerance * b && (phi[3] + phi[1]).abs() <= 2.0 * tolerance * b
    );
}

#[test]
fn two_unequal_modes_follow_the_quadratic_of_their_pencil() {
    let (w1, w2, m) = (TAU * 329.63, TAU * 352.74, 0.08);
    let basis = LoadedModalBasis::new(&[w1, w2], &[0.0, m, m, 0.0]).expect("definite");
    // det(K - lambda M) = (1 - m^2) lambda^2 - (k1 + k2) lambda + k1 k2 = 0 with k = w^2, the
    // same products the basis forms. The discriminant (k1 + k2)^2 - 4 (1 - m^2) k1 k2 is about a
    // hundredth of its first term here, so it is formed as (k1 - k2)^2 + 4 m^2 k1 k2, free of that
    // cancellation, and the smaller root in the form free of cancellation too.
    let (k1, k2) = (w1 * w1, w2 * w2);
    let (a, b, c) = (1.0 - m * m, k1 + k2, k1 * k2);
    let root = ((k1 - k2) * (k1 - k2) + 4.0 * m * m * c).sqrt();
    within_the_certificate(&basis, 0, 2.0 * c / (b + root), k2);
    within_the_certificate(&basis, 1, (b + root) / (2.0 * a), k2);
    // The mutual mass pushes the pair apart.
    let w = basis.angular_frequencies_rad_s();
    assert!(w[0] < w1 && w[1] > w2, "{w:?}");
}

/// `gamma(k) = k u / (1 - k u)` with the unit roundoff `u = eps / 2`: the relative error bound of
/// `k` successive roundings (N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd
/// ed., SIAM 2002, Lemma 3.1). An inner product of length `k` is within `gamma(k)` of the sum of
/// its terms' magnitudes.
fn gamma(k: usize) -> f64 {
    let ku = k as f64 * 0.5 * EPS;
    ku / (1.0 - ku)
}

#[test]
fn the_basis_is_normalized_in_the_loaded_mass_and_diagonalizes_the_stiffness() {
    let w = [TAU * 329.63, TAU * 331.2, TAU * 1_318.0];
    let added = [
        0.012, 0.004, -0.001, //
        0.004, 0.030, 0.002, //
        -0.001, 0.002, 0.006,
    ];
    let n = 3;
    let basis = LoadedModalBasis::new(&w, &added).expect("definite");
    let mut mass = added.to_vec();
    let mut stiffness = vec![0.0; 9];
    for k in 0..3 {
        mass[4 * k] += 1.0;
        stiffness[4 * k] = w[k] * w[k];
    }
    // Every Gershgorin disc of M lies right of 1, so its smallest eigenvalue is at least 1 and the
    // inverse-mass norm of a vector is at most its Euclidean norm.
    let gershgorin = (0..n)
        .map(|k| {
            mass[k * n + k]
                - (0..n)
                    .filter(|&l| l != k)
                    .map(|l| mass[k * n + l].abs())
                    .sum::<f64>()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(gershgorin >= 1.0, "{gershgorin}");
    let unit = basis.congruent(&mass).expect("symmetric");
    let diagonal = basis.congruent(&stiffness).expect("symmetric");
    let phi = basis.basis();
    let loaded = basis.angular_frequencies_rad_s();
    let lambda = |j: usize| loaded[j] * loaded[j];
    // `|phi_i|^T |A| |phi_j|`: the scale of the rounding of `phi_i^T A phi_j`. `congruent` forms
    // `A Phi` and then `Phi^T (A Phi)`, two inner products of length n: within `gamma(2 n)` of it.
    let absolute = |a: &[f64], i: usize, j: usize| -> f64 {
        (0..n)
            .map(|k| {
                (0..n)
                    .map(|l| phi[k * n + i].abs() * a[k * n + l].abs() * phi[l * n + j].abs())
                    .sum::<f64>()
            })
            .sum()
    };
    // fs-modal normalizes each vector itself: `phi^T (M phi)` in 2 n roundings, a square root and
    // a division of each entry, so the exact `phi_i^T M phi_i` is 1 within this.
    let normalized = |i: usize| gamma(2 * n + 2) * absolute(&mass, i, i) + gamma(2);
    // The exact residual `r_j = K phi_j - lambda_j M phi_j` in the inverse-mass norm. fs-modal
    // states it rounded: each entry of `K phi` and `M phi` is n fused steps and the difference
    // one more, `gamma(n + 1)` of `|K| |phi_j| + lambda_j |M| |phi_j|`; the norm itself (a Cholesky
    // solve, an inner product and a square root) is within `gamma(4 n + 2)` of its value. The
    // stated `lambda_j` is `w'_j^2` within `gamma(3)`, which the last factor carries.
    let certified = |j: usize| -> f64 {
        let spread = (0..n)
            .map(|k| {
                let row: f64 = (0..n)
                    .map(|l| {
                        (stiffness[k * n + l].abs() + lambda(j) * mass[k * n + l].abs())
                            * phi[l * n + j].abs()
                    })
                    .sum();
                row * row
            })
            .sum::<f64>()
            .sqrt();
        basis.residuals()[j] * (1.0 + gamma(4 * n + 2)) + gamma(n + 1) * spread * (1.0 + gamma(3))
    };
    let norm = |i: usize| (1.0 + normalized(i)).sqrt();
    // For i != j, `(lambda_i - lambda_j) phi_i^T M phi_j = phi_i^T r_j - phi_j^T r_i` exactly (K and
    // M symmetric), and `|phi_i^T r_j| <= ||phi_i||_M ||r_j||_{M^-1}`: M-orthogonality to what the
    // residuals certify, over the gap between the two modes.
    let orthogonal = |i: usize, j: usize| -> f64 {
        let gap = (lambda(i) - lambda(j)).abs() - gamma(4) * (lambda(i) + lambda(j));
        assert!(gap > 0.0, "modes {i} and {j} are not separated");
        (norm(i) * certified(j) + norm(j) * certified(i)) / gap
    };
    for i in 0..3 {
        for j in 0..3 {
            let rounding = gamma(2 * n);
            let (m_bound, k_bound) = if i == j {
                // `phi_i^T K phi_i = lambda_i phi_i^T M phi_i + phi_i^T r_i`, against `w'_i^2`.
                (
                    normalized(i),
                    lambda(i) * (normalized(i) + gamma(3)) + norm(i) * certified(i),
                )
            } else {
                // `phi_i^T K phi_j = lambda_j phi_i^T M phi_j + phi_i^T r_j`, or the same with i and
                // j exchanged; the computed entry is mirrored, so either bound holds.
                let o = orthogonal(i, j);
                (
                    o,
                    (lambda(j) * o + norm(i) * certified(j))
                        .min(lambda(i) * o + norm(j) * certified(i)),
                )
            };
            let m_bound = m_bound + rounding * absolute(&mass, i, j);
            let k_bound = k_bound + rounding * absolute(&stiffness, i, j);
            let want = if i == j { 1.0 } else { 0.0 };
            assert!(
                (unit[i * 3 + j] - want).abs() <= m_bound,
                "M ({i}, {j}): {} against {want}, bound {m_bound}",
                unit[i * 3 + j]
            );
            let want = if i == j { lambda(i) } else { 0.0 };
            assert!(
                (diagonal[i * 3 + j] - want).abs() <= k_bound,
                "K ({i}, {j}): {} against {want}, bound {k_bound}",
                diagonal[i * 3 + j]
            );
            assert_eq!(unit[i * 3 + j], unit[j * 3 + i], "mirrored bit for bit");
        }
    }
    assert!(
        loaded.windows(2).all(|pair| pair[0] <= pair[1]),
        "ascending"
    );
}

#[test]
fn the_held_loss_factor_is_sinc_squared_and_nan_for_a_step_that_is_not_positive_and_finite() {
    let w = [TAU * 329.63, TAU * 1_318.5];
    let basis = LoadedModalBasis::new(&w, &[0.0; 4]).expect("definite");
    let dt = 1.0 / 44_100.0;
    let factors = basis.held_loss_factors(dt);
    assert_eq!(factors.len(), 2);
    for (j, &factor) in factors.iter().enumerate() {
        let x = 0.5 * basis.angular_frequencies_rad_s()[j] * dt;
        let sinc = x.sin() / x;
        assert_eq!(factor, sinc * sinc, "mode {j}");
    }
    for bad in [0.0, -dt, f64::NAN, f64::INFINITY] {
        let factors = basis.held_loss_factors(bad);
        assert_eq!(factors.len(), 2);
        assert!(factors.iter().all(|f| f.is_nan()), "dt {bad}: {factors:?}");
    }
}

#[test]
fn an_attachment_carried_to_the_loaded_modes_keeps_its_power() {
    let w = [TAU * 329.63, TAU * 352.74];
    let basis = LoadedModalBasis::new(&w, &[0.01, 0.003, 0.003, 0.02]).expect("definite");
    let strike = ModalAttachment {
        id: StableId::new("strike").expect("id"),
        kind: PortKind::MechanicalForceVelocity,
        width: NonZeroUsize::new(2).expect("two"),
        weights: vec![2.1, -1.4, 0.3, 0.9],
    };
    let carried = basis.loaded_attachment(strike.clone()).expect("carried");
    assert_eq!(
        (carried.id.clone(), carried.width),
        (strike.id.clone(), strike.width)
    );
    // A row's flow is the same whether formed in the loaded coordinates or in the unloaded ones
    // they stand for: s' . y = s . (Phi y).
    let phi = basis.basis();
    for y in [[1.0, 0.0], [0.0, 1.0], [0.37, -1.2]] {
        let x = basis.unloaded(&y);
        for row in 0..2 {
            let loaded: f64 = (0..2).map(|j| carried.weights[row * 2 + j] * y[j]).sum();
            let unloaded: f64 = (0..2).map(|k| strike.weights[row * 2 + k] * x[k]).sum();
            // Each side is two nested inner products of length n = 2 (the carried weights `s Phi`
            // and then `. y`; `Phi y` and then `s .`), each within `gamma(2 n)` of the sum of the
            // magnitudes of the products `s_k Phi_kj y_j`.
            let magnitudes: f64 = (0..2)
                .flat_map(|k| (0..2).map(move |j| (k, j)))
                .map(|(k, j)| (strike.weights[row * 2 + k] * phi[k * 2 + j] * y[j]).abs())
                .sum();
            let bound = 2.0 * gamma(4) * magnitudes;
            assert!(
                (loaded - unloaded).abs() <= bound,
                "row {row} at {y:?}: {loaded} against {unloaded}, bound {bound}"
            );
        }
    }
}

#[test]
fn the_structure_runs_the_loaded_modes_undamped_behind_the_modal_port() {
    let w = [TAU * 329.63, TAU * 352.74];
    let basis = LoadedModalBasis::new(&w, &[0.01, 0.003, 0.003, 0.02]).expect("definite");
    let strike = ModalAttachment {
        id: StableId::new("strike").expect("id"),
        kind: PortKind::MechanicalForceVelocity,
        width: NonZeroUsize::MIN,
        weights: vec![2.1, -1.4],
    };
    let structure = basis
        .structure(
            44_100,
            ModalAcousticTimeBudget::audible_reference(),
            StableId::new("modal").expect("id"),
            PortKind::MechanicalForceVelocity,
            vec![strike.clone()],
        )
        .expect("structure");
    let ports = structure.ports();
    assert_eq!(ports.len(), 3);
    assert_eq!(
        (ports[0].id().as_str(), ports[0].width(), ports[0].role()),
        ("modal", 2, PortRole::EffortIn)
    );
    assert_eq!((ports[1].id().as_str(), ports[1].width()), ("strike", 1));
    assert_eq!(ports[2].id().as_str(), "heat");
    let modes = structure.model().modes();
    for (mode, &loaded) in modes.iter().zip(basis.angular_frequencies_rad_s()) {
        assert_eq!(mode.angular_frequency_rad_s, loaded);
        assert_eq!(mode.damping_ratio, 0.0);
    }
    // One step from rest at 300 K with a unit effort on one input at a time: the flows out of the
    // modal port (two), the strike row, and the heat port.
    let carried = basis.loaded_attachment(strike).expect("carried").weights;
    let step = |efforts: [f64; 3]| {
        let inputs = [efforts[0], efforts[1], efforts[2], 300.0];
        let (start, mut next, mut outputs) = (structure.at_rest(), structure.at_rest(), [0.0; 4]);
        with_cx(|cx| {
            structure.propose(&start, 1.0 / 44_100.0, &inputs, &mut next, &mut outputs, cx)
        })
        .expect("step");
        outputs
    };
    let alone = [step([1.0, 0.0, 0.0]), step([0.0, 1.0, 0.0])];
    for (j, out) in alone.iter().enumerate() {
        // The modal port is the loaded coordinates themselves: a force on one moves that one
        // alone, and the strike row reads it through its carried weight (two roundings of at
        // most eps / 2 each side).
        assert!(out[j] != 0.0 && out[1 - j] == 0.0, "mode {j}: {out:?}");
        assert!(
            (out[2] - carried[j] * out[j]).abs() <= 4.0 * EPS * out[2].abs(),
            "mode {j}: {} against {}",
            out[2],
            carried[j] * out[j]
        );
        assert_eq!(out[3], 0.0, "undamped: no heat");
    }
    // A unit force on the strike row pushes each loaded mode with its carried weight. From rest
    // fs-modal moves an undamped mode by (F / w^2) (1 - cos w dt), so a flow is linear in the
    // force to seven roundings of at most eps / 2; the strike row's flow is their carried sum.
    let struck = step([0.0, 0.0, 1.0]);
    let (mut flow, mut size) = (0.0_f64, 0.0_f64);
    for j in 0..2 {
        assert!(
            (struck[j] - carried[j] * alone[j][j]).abs() <= 4.0 * EPS * struck[j].abs(),
            "mode {j}: {} against {}",
            struck[j],
            carried[j] * alone[j][j]
        );
        flow += carried[j] * struck[j];
        size += (carried[j] * struck[j]).abs();
    }
    assert!(
        (struck[2] - flow).abs() <= 4.0 * EPS * size,
        "{} against {flow}",
        struck[2]
    );
    assert_eq!(struck[3], 0.0);
}

#[test]
fn the_basis_refuses_what_it_cannot_represent() {
    let w = TAU * 329.63;
    let refused = |w: &[f64], added: &[f64]| LoadedModalBasis::new(w, added).unwrap_err();
    assert_eq!(refused(&[], &[]), ModalStructureRefusal::NoModes);
    for value in [0.0, -w, f64::NAN, f64::INFINITY] {
        match refused(&[w, value], &[0.0; 4]) {
            ModalStructureRefusal::Frequency {
                mode: 1,
                value: got,
            } => {
                assert_eq!(got.to_bits(), value.to_bits());
            }
            other => panic!("{value}: {other:?}"),
        }
    }
    let next = |x: f64| f64::from_bits(x.to_bits() + 1);
    assert_eq!(
        refused(&[w, w], &[0.0; 3]),
        ModalStructureRefusal::MatrixShape {
            expected: 4,
            found: 3
        }
    );
    assert_eq!(
        refused(&[w, w], &[0.0, f64::NAN, f64::NAN, 0.0]),
        ModalStructureRefusal::NonFiniteEntry { row: 0, column: 1 }
    );
    assert_eq!(
        refused(&[w, w], &[0.0, 0.01, next(0.01), 0.0]),
        ModalStructureRefusal::Asymmetric { row: 1, column: 0 }
    );
    // 1 + M_a with a pivot at zero, and one indefinite by its mutual mass.
    for added in [vec![-1.0], vec![0.0, 1.5, 1.5, 0.0]] {
        let n = if added.len() == 1 { 1 } else { 2 };
        assert!(matches!(
            refused(&vec![w; n], &added),
            ModalStructureRefusal::LoadedMass(_)
        ));
    }
    let basis = LoadedModalBasis::new(&[w, w], &[0.0; 4]).expect("unloaded");
    assert_eq!(
        basis.congruent(&[1.0, 2.0, next(2.0), 1.0]),
        Err(ModalStructureRefusal::Asymmetric { row: 1, column: 0 })
    );
    // A signed zero equals its mirror.
    LoadedModalBasis::new(&[w, w], &[0.0, 0.0, -0.0, 0.0]).expect("symmetric");
    // More modes than fs-la's dense eigensolver admits (4 n^2 + 3 n above 64 Mi at n = 4096),
    // refused before the matrix is read.
    assert!(matches!(
        refused(&vec![w; 4_096], &[]),
        ModalStructureRefusal::LoadedSize(_)
    ));
    // Finite weights that overflow when carried to the in-phase mode (about 0.7 of each, summed)
    // are refused with the attachment's index in the list given.
    let coupled = LoadedModalBasis::new(&[w, w], &[0.0, 0.02, 0.02, 0.0]).expect("definite");
    let row = |id: &str, weights: Vec<f64>| ModalAttachment {
        id: StableId::new(id).expect("id"),
        kind: PortKind::MechanicalForceVelocity,
        width: NonZeroUsize::MIN,
        weights,
    };
    let huge = vec![1.5e308, 1.5e308];
    assert_eq!(
        coupled.loaded_attachment(row("huge", huge.clone())),
        Err(ModalStructureRefusal::NonFiniteWeight { attachment: 0 })
    );
    let refusal = coupled
        .structure(
            44_100,
            ModalAcousticTimeBudget::audible_reference(),
            StableId::new("modal").expect("id"),
            PortKind::MechanicalForceVelocity,
            vec![row("fine", vec![1.0, 0.5]), row("huge", huge)],
        )
        .map(|_| ())
        .unwrap_err();
    assert_eq!(
        refusal,
        ModalStructureRefusal::NonFiniteWeight { attachment: 1 }
    );
}
