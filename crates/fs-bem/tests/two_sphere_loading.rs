//! The two-sphere oracle of the multi-port radiation law: two pulsating spheres of radius 1 m with
//! centres 3 m apart, each sphere's uniform normal velocity one mode. The closed-body loading
//! bake's mutual impedance `Z_12` (the force on one sphere per unit velocity of the other), its
//! self impedance in the other's presence, and the pair-mean law's matrices sampled from it are
//! quantities no input contains; they are compared with an independent series solution.
//!
//! The reference: the field is a sum of outgoing partial waves `A_n h_n(k r_1) P_n` about one
//! centre and `B_n h_n(k r_2) P_n` about the other (exp(-i w t)). On each sphere the normal
//! derivative of the whole field is projected onto `P_n` by Gauss-Legendre quadrature, which
//! needs no addition theorem, and set to `i w rho` times the sphere's velocity; a small pivoted
//! elimination solves for `A` and `B`. `Z_rs` is `4 pi a^2` times the mean pressure on sphere `r`
//! for unit velocity of sphere `s`. No BEM kernel, panel or production solver is used in it. It is
//! checked for its own convergence in the number of partial waves and quadrature nodes, and
//! against a quantity it does not form from the impedance: the power the same partial waves carry
//! to the far field, `pi / (rho c)` times the integral of `|D(mu)|^2` over the axis angle, must
//! equal the power the spheres deliver, `Re(v^T Z v) / 2`. With equal spheres the mirror between
//! them makes `Z_12 = Z_21` and `Z_11 = Z_22` hold for any projection formula, right or wrong, so
//! those equalities are not used as a check. A dropped term in the projections breaks the power
//! balance by at least 5e-5 of it at these ka, a flipped sign by 0.2 to 38 percent; the balance
//! is required to 2e-11.
//!
//! The tolerance comes from panel refinement: three meshes at a constant refinement ratio
//! `r = 4/3` (18, 24 and 32 latitude rings, 2448 to 7936 panels, all on the dense route) give the
//! observed order `p` and the grid convergence index of the finest,
//! `GCI = 1.25 |Z_3 - Z_2| / (r^p - 1)` (the index from P. J. Roache, J. Fluids Eng. 116 (1994)
//! 405-413; the factor of safety 1.25 for three meshes with an observed order from I. B. Celik et
//! al., J. Fluids Eng. 130 (2008) 078001). Each entry must converge monotonically with an observed
//! order within `[1, 3]` around the formal order 2 of flat panels with one collocation point
//! each (outside it the ladder is not in the asymptotic range the factor of safety rests on), and
//! the finest must lie within its GCI of the reference.
//!
//! The law: a `PairMeanLaw` of two modes, sphere 1 at ka 0.5 and sphere 2 at ka 1 (half-widths
//! zero), asks the bake for ka 0.5, 0.75 and 1. Its added mass and resistance from the finest
//! bake must equal `-Im Z / w` and `Re Z` of the reference at the pair-mean frequencies, within
//! the GCI of the entries they are taken from: this pins the law's sign convention and its choice
//! of frequency against a solution it never saw.
//!
//! Sphere 2's mesh is sphere 1's turned half a revolution about the x axis (`y` and `z` negated,
//! exact in floating point), so exchanging the spheres maps the discrete operator onto itself. A
//! collocation bake of two meshes that are not images of each other is not reciprocal beyond its
//! discretization error, and the bake refuses a row above its dense reciprocity bound.
//!
//! The panel ladder is ignored by default (minutes in release, about 3 GB at 7936 panels): run it
//! with `cargo test --release -p fs-bem --test two_sphere_loading -- --ignored`. The series test
//! runs with the crate's tests.

// The series is written index by index in the notation of its equations.
#![allow(clippy::needless_range_loop, clippy::many_single_char_names)]

use std::f64::consts::{PI, TAU};

use fs_bem::helmholtz::{self, Medium};
use fs_bem::modal_loading::{
    LoadingRoute, ModalLoadingBake, PairMeanLaw, SurfaceModeTraces, modal_radiation_loading,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::hash_domain;
use fs_math::c64::C64;

const RADIUS_M: f64 = 1.0;
const SEPARATION_M: f64 = 3.0;
/// The refinement ratio of the panel ladder.
const RATIO: f64 = 4.0 / 3.0;
/// The ka of the two modes of the law and of their pair mean.
const KA: [f64; 3] = [0.5, 0.75, 1.0];
/// The series' own convergence, as a fraction of `|Z_11|`, checked by the series test.
const SERIES: f64 = 1.0e-11;

/// Outgoing spherical Hankel functions `h_n(x)`, `n < count`, by upward recurrence:
/// `h_0 = -i e^{ix} / x`, `h_1 = -e^{ix} (x + i) / x^2`. The complex value is accurate (the
/// recurrence is stable for the growing `y_n`); its real part `j_n` alone loses digits at high
/// order, and the series never uses it alone.
fn hankel(x: f64, count: usize) -> Vec<C64> {
    let e = C64::new(x.cos(), x.sin());
    let mut h = vec![
        e * C64::new(0.0, -1.0 / x),
        -(e * C64::new(x, 1.0)).scale(1.0 / (x * x)),
    ];
    for n in 1..count {
        let next = h[n].scale((2 * n + 1) as f64 / x) - h[n - 1];
        h.push(next);
    }
    h.truncate(count.max(1));
    h
}

/// `h_n'(x)` from `h_0' = -h_1` and `h_n' = h_{n-1} - (n + 1) h_n / x`; needs `count + 1` values.
fn hankel_derivative(h: &[C64], x: f64, n: usize) -> C64 {
    if n == 0 {
        -h[1]
    } else {
        h[n - 1] - h[n].scale((n + 1) as f64 / x)
    }
}

/// `P_n(c)` and `P_n'(c)`, `n < count`, with `P'_{n+1} = P'_{n-1} + (2n + 1) P_n` (no division by
/// `1 - c^2`, which vanishes on the axis).
fn legendre(c: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
    let mut p = vec![1.0, c];
    let mut d = vec![0.0, 1.0];
    for n in 1..count {
        p.push(((2 * n + 1) as f64 * c * p[n] - n as f64 * p[n - 1]) / (n + 1) as f64);
        d.push(d[n - 1] + (2 * n + 1) as f64 * p[n]);
    }
    p.truncate(count);
    d.truncate(count);
    (p, d)
}

/// Gauss-Legendre nodes and weights on [-1, 1] by Newton's method on `P_q`.
fn gauss(q: usize) -> Vec<(f64, f64)> {
    (0..q)
        .map(|i| {
            let mut x = (PI * (i as f64 + 0.75) / (q as f64 + 0.5)).cos();
            for _ in 0..100 {
                let (p, d) = legendre(x, q + 1);
                let step = p[q] / d[q];
                x -= step;
                if step.abs() <= 4.0 * f64::EPSILON {
                    break;
                }
            }
            let (_, d) = legendre(x, q + 1);
            (x, 2.0 / ((1.0 - x * x) * d[q] * d[q]))
        })
        .collect()
}

/// For a sphere of radius `a` at the origin and outgoing waves about a centre at `s` on the axis:
/// the projections onto `P_n` (`n < count`) of wave `m`'s normal derivative (`g`) and value (`f`)
/// on the sphere, row-major `[m * count + n]`.
fn projections(k: f64, a: f64, s: f64, count: usize, nodes: &[(f64, f64)]) -> (Vec<C64>, Vec<C64>) {
    let mut g = vec![C64::ZERO; count * count];
    let mut f = vec![C64::ZERO; count * count];
    for &(mu, weight) in nodes {
        let rho = (a * a + s * s - 2.0 * a * s * mu).sqrt();
        let c = (a * mu - s) / rho;
        let d_rho = (a - s * mu) / rho;
        let d_c = mu / rho - (a * mu - s) * (a - s * mu) / (rho * rho * rho);
        let h = hankel(k * rho, count + 1);
        let (p_other, d_other) = legendre(c, count);
        let (p_own, _) = legendre(mu, count);
        for m in 0..count {
            let value = h[m].scale(p_other[m]);
            let derivative = hankel_derivative(&h, k * rho, m).scale(k * d_rho * p_other[m])
                + h[m].scale(d_other[m] * d_c);
            for n in 0..count {
                let w = weight * 0.5 * (2 * n + 1) as f64 * p_own[n];
                g[m * count + n] = g[m * count + n] + derivative.scale(w);
                f[m * count + n] = f[m * count + n] + value.scale(w);
            }
        }
    }
    (g, f)
}

/// Gaussian elimination with partial pivoting on a dense complex system.
fn solve(mut m: Vec<C64>, mut b: Vec<C64>) -> Vec<C64> {
    let n = b.len();
    for p in 0..n {
        let pivot = (p..n)
            .max_by(|&i, &j| m[i * n + p].abs().total_cmp(&m[j * n + p].abs()))
            .expect("rows");
        for c in 0..n {
            m.swap(p * n + c, pivot * n + c);
        }
        b.swap(p, pivot);
        for r in p + 1..n {
            let factor = m[r * n + p] / m[p * n + p];
            for c in p..n {
                m[r * n + c] = m[r * n + c] - factor * m[p * n + c];
            }
            b[r] = b[r] - factor * b[p];
        }
    }
    for p in (0..n).rev() {
        let mut sum = b[p];
        for c in p + 1..n {
            sum = sum - m[p * n + c] * b[c];
        }
        b[p] = sum / m[p * n + p];
    }
    b
}

/// The series with `count` partial waves about each centre and `q` quadrature nodes. Sphere 1 is
/// at the origin, sphere 2 at `+SEPARATION_M` on the axis; both expansions use the same axis.
struct Series {
    medium: Medium,
    k: f64,
    count: usize,
    q: usize,
    /// `h_n(ka)`, `n <= count`.
    h: Vec<C64>,
    /// The `2 count` boundary rows: sphere 1's projections, then sphere 2's.
    matrix: Vec<C64>,
    /// The other sphere's waves' values projected on sphere 1 (`f1`) and on sphere 2 (`f2`).
    f1: Vec<C64>,
    f2: Vec<C64>,
}

impl Series {
    fn new(medium: Medium, ka: f64, count: usize, q: usize) -> Self {
        let (a, d) = (RADIUS_M, SEPARATION_M);
        let k = ka / a;
        let nodes = gauss(q);
        let (g1, f1) = projections(k, a, d, count, &nodes);
        let (g2, f2) = projections(k, a, -d, count, &nodes);
        let h = hankel(ka, count + 1);
        let size = 2 * count;
        let mut matrix = vec![C64::ZERO; size * size];
        for n in 0..count {
            let own = hankel_derivative(&h, ka, n).scale(k);
            matrix[n * size + n] = own;
            matrix[(count + n) * size + count + n] = own;
            for m in 0..count {
                // Sphere 1's row n takes sphere 2's waves B_m; sphere 2's row takes sphere 1's A_m.
                matrix[n * size + count + m] = g1[m * count + n];
                matrix[(count + n) * size + m] = g2[m * count + n];
            }
        }
        Self {
            medium,
            k,
            count,
            q,
            h,
            matrix,
            f1,
            f2,
        }
    }

    /// `[A_0 .. A_{count-1}, B_0 .. B_{count-1}]` for the spheres' normal velocities.
    fn waves(&self, velocity: [f64; 2]) -> Vec<C64> {
        let omega = self.k * self.medium.sound_speed;
        let mut rhs = vec![C64::ZERO; 2 * self.count];
        for (sphere, &v) in velocity.iter().enumerate() {
            rhs[sphere * self.count] = C64::new(0.0, omega * self.medium.density * v);
        }
        solve(self.matrix.clone(), rhs)
    }

    /// Row-major `[Z_11, Z_12, Z_21, Z_22]`: column 0 is sphere 1 moving.
    fn impedance(&self) -> [C64; 4] {
        let surface = 4.0 * PI * RADIUS_M * RADIUS_M;
        let count = self.count;
        let mut z = [C64::ZERO; 4];
        for column in 0..2 {
            let mut velocity = [0.0; 2];
            velocity[column] = 1.0;
            let x = self.waves(velocity);
            let (a_n, b_n) = x.split_at(count);
            let mut mean = [a_n[0] * self.h[0], b_n[0] * self.h[0]];
            for m in 0..count {
                mean[0] = mean[0] + b_n[m] * self.f1[m * count];
                mean[1] = mean[1] + a_n[m] * self.f2[m * count];
            }
            z[column] = mean[0].scale(surface);
            z[2 + column] = mean[1].scale(surface);
        }
        z
    }

    /// The power the partial waves carry to the far field [W]: with `h_n(kr) -> (-i)^{n+1}
    /// e^{ikr} / (kr)`, the pressure is `e^{ikr} D(mu) / r` with
    /// `D = sum (-i)^{n+1} (A_n + B_n e^{-ikd mu}) P_n(mu) / k`, and the power is
    /// `pi / (rho c)` times the integral of `|D|^2` over `mu`.
    fn far_field_power(&self, velocity: [f64; 2]) -> f64 {
        let x = self.waves(velocity);
        let (a_n, b_n) = x.split_at(self.count);
        let mut total = 0.0;
        for (mu, weight) in gauss(self.q) {
            let (p, _) = legendre(mu, self.count);
            let phase = -self.k * SEPARATION_M * mu;
            let shift = C64::new(phase.cos(), phase.sin());
            let mut turn = C64::new(0.0, -1.0);
            let mut d = C64::ZERO;
            for n in 0..self.count {
                d = d + (turn * (a_n[n] + b_n[n] * shift)).scale(p[n] / self.k);
                turn = turn * C64::new(0.0, -1.0);
            }
            total += weight * (d.re * d.re + d.im * d.im);
        }
        PI / (self.medium.density * self.medium.sound_speed) * total
    }
}

/// The reference `[Z_11, Z_12, Z_21, Z_22]` of the converged series.
fn reference(medium: Medium, ka: f64) -> [C64; 4] {
    Series::new(medium, ka, 16, 64).impedance()
}

/// A closed sphere of radius `a` centred at `z0` on the axis, from latitude rings that share their
/// vertices exactly (`4 L (L - 1)` triangles), outward.
fn ring_sphere(latitude: usize, a: f64, z0: f64) -> Vec<[[f64; 3]; 3]> {
    let longitude = 2 * latitude;
    let mut vertices = vec![[0.0, 0.0, z0 + a]];
    for j in 1..latitude {
        let theta = PI * j as f64 / latitude as f64;
        for i in 0..longitude {
            let phi = TAU * i as f64 / longitude as f64;
            vertices.push([
                a * theta.sin() * phi.cos(),
                a * theta.sin() * phi.sin(),
                z0 + a * theta.cos(),
            ]);
        }
    }
    let south = vertices.len();
    vertices.push([0.0, 0.0, z0 - a]);
    let ring = |j: usize, i: usize| 1 + j * longitude + i % longitude;
    let mut triangles = Vec::new();
    for i in 0..longitude {
        triangles.push([vertices[0], vertices[ring(0, i)], vertices[ring(0, i + 1)]]);
        for j in 0..latitude - 2 {
            let (p, q, r, s) = (
                ring(j, i),
                ring(j + 1, i),
                ring(j + 1, i + 1),
                ring(j, i + 1),
            );
            triangles.push([vertices[p], vertices[q], vertices[r]]);
            triangles.push([vertices[p], vertices[r], vertices[s]]);
        }
        triangles.push([
            vertices[ring(latitude - 2, i)],
            vertices[south],
            vertices[ring(latitude - 2, i + 1)],
        ]);
    }
    triangles
}

/// The bake at `omegas` on two ring spheres of `latitude` rings each, centred at `z = -d/2` and,
/// turned half a revolution about the x axis, at `z = +d/2`; with its panel count. A refusal (a
/// row above the dense reciprocity bound among them) fails with the rung and the error's own
/// numbers.
#[expect(
    clippy::result_large_err,
    reason = "ModalLoadingError is the bake's refusal type; boxing it is outside this test"
)]
fn bem(medium: Medium, omegas: &[f64], latitude: usize) -> (ModalLoadingBake, usize) {
    let mut triangles = ring_sphere(latitude, RADIUS_M, -0.5 * SEPARATION_M);
    let half = triangles.len();
    let turned: Vec<[[f64; 3]; 3]> = triangles
        .iter()
        .map(|triangle| triangle.map(|[x, y, z]| [x, -y, -z]))
        .collect();
    triangles.extend(turned);
    let surface = SpherePanels::from_triangles(triangles).expect("two closed spheres");
    let panels = surface.areas().len();
    let shapes: Vec<Vec<f64>> = (0..2)
        .map(|sphere| {
            (0..panels)
                .map(|i| {
                    if (i < half) == (sphere == 0) {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
    let bake = helmholtz::with_dense_workers(workers, || {
        modal_radiation_loading(
            SurfaceModeTraces::ClosedBody {
                surface: &surface,
                shapes: &shapes,
            },
            LoadingRoute::ClosedBodyGeometryPolicy,
            omegas,
            medium,
            hash_domain("analytic-fixture", "two pulsating spheres".as_bytes()),
            4 * omegas.len(),
            || false,
        )
    })
    .unwrap_or_else(|error| {
        panic!("bake refused at {latitude} rings ({panels} panels): {error:?}")
    });
    (bake, panels)
}

#[test]
fn the_series_reference_converges_and_balances_the_power_it_radiates() {
    let medium = Medium::air();
    for ka in KA {
        let coarse = Series::new(medium, ka, 12, 48).impedance();
        let series = Series::new(medium, ka, 16, 64);
        let fine = series.impedance();
        let scale = fine[0].abs();
        for (entry, (c, f)) in coarse.iter().zip(&fine).enumerate() {
            assert!(
                (*c - *f).abs() <= SERIES * scale,
                "ka {ka} entry {entry}: {c:?} against {f:?}"
            );
        }
        // Each side of the balance carries the series' truncation, within SERIES of |Z_11| per
        // entry (above); a velocity pair weights the four entries by at most v1^2 + v2^2 in all.
        for velocity in [[1.0, 0.0], [1.0, 1.0], [1.0, -1.0]] {
            let delivered = 0.5
                * (0..2)
                    .flat_map(|r| (0..2).map(move |s| (r, s)))
                    .map(|(r, s)| velocity[r] * velocity[s] * fine[2 * r + s].re)
                    .sum::<f64>();
            let radiated = series.far_field_power(velocity);
            let weight = velocity[0] * velocity[0] + velocity[1] * velocity[1];
            assert!(
                (delivered - radiated).abs() <= 2.0 * SERIES * scale * weight,
                "ka {ka} velocity {velocity:?}: delivered {delivered}, radiated {radiated}"
            );
        }
        // The mutual impedance is not small beside the self term, so a missing or misplaced mutual
        // term would show: the pulsating sphere's field at the other centre gives
        // |Z_12| / |Z_alone| = a / d, a third, at low ka.
        let alone = C64::new(ka * ka, -ka).scale(
            4.0 * PI * RADIUS_M * RADIUS_M * medium.density * medium.sound_speed / (1.0 + ka * ka),
        );
        assert!(fine[1].abs() >= 0.1 * alone.abs(), "ka {ka}");
    }
}

#[test]
#[ignore = "two-sphere panel ladder to 7936 panels: minutes in release, about 3 GB"]
fn the_bake_and_the_pair_mean_law_converge_to_the_two_sphere_reference_within_the_gci() {
    let medium = Medium::air();
    let omega = |ka: f64| ka / RADIUS_M * medium.sound_speed;
    // Sphere 1 a mode at ka 0.5, sphere 2 a mode at ka 1; the law asks for ka 0.5, 0.75 and 1.
    let law =
        PairMeanLaw::new(&[omega(KA[0]), omega(KA[2])], &[0.0, 0.0], &[None, None]).expect("law");
    let omegas = law.frequencies();
    assert_eq!(omegas.len(), 3);
    let latitudes = [18, 24, 32];
    let rungs: Vec<(ModalLoadingBake, usize)> =
        latitudes.iter().map(|&l| bem(medium, &omegas, l)).collect();
    let panels: Vec<usize> = rungs.iter().map(|(_, p)| *p).collect();
    // gci[ka][entry], with the reference's own convergence added to each.
    let mut gci = [[0.0; 4]; 3];
    let mut truth = [[C64::ZERO; 4]; 3];
    for (i, &w) in omegas.iter().enumerate() {
        let ka = KA[i];
        truth[i] = reference(medium, ka);
        let slack = SERIES * truth[i][0].abs();
        for entry in 0..4 {
            let z: Vec<C64> = rungs
                .iter()
                .map(|(bake, _)| bake.at_omega(w).expect("declared").impedance()[entry])
                .collect();
            let (e21, e32) = ((z[1] - z[0]).abs(), (z[2] - z[1]).abs());
            let order = (e21 / e32).ln() / RATIO.ln();
            let index = 1.25 * e32 / (RATIO.powf(order) - 1.0);
            let errors: Vec<f64> = z.iter().map(|z| (*z - truth[i][entry]).abs()).collect();
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"two-sphere-loading\",\"ka\":{ka},\"entry\":{entry},\"panels\":[{},{},{}],\"observed\":[{},{}],\"oracle\":[{},{}],\"errors\":[{},{},{}],\"order\":{order},\"gci\":{index}}}",
                panels[0],
                panels[1],
                panels[2],
                z[2].re,
                z[2].im,
                truth[i][entry].re,
                truth[i][entry].im,
                errors[0],
                errors[1],
                errors[2],
            );
            assert!(
                e32 < e21,
                "ka {ka} entry {entry}: the ladder does not converge"
            );
            assert!(
                (1.0..=3.0).contains(&order),
                "ka {ka} entry {entry}: observed order {order} outside [1, 3]"
            );
            assert!(
                errors[0] > errors[1] && errors[1] > errors[2],
                "ka {ka} entry {entry}: errors {errors:?}"
            );
            assert!(
                errors[2] <= index + slack,
                "ka {ka} entry {entry}: finest error {} above its GCI {index} (order {order})",
                errors[2]
            );
            gci[i][entry] = index + slack;
        }
    }
    // The law on the finest bake. Its own entries come from Z_11 at ka 0.5 and Z_22 at ka 1, its
    // mutual entry from (Z_12 + Z_21) / 2 at ka 0.75; |Im| and |Re| of an entry's error are each at
    // most its modulus, so within the GCI (over w for the added mass).
    let loading = law.sample(&rungs[2].0).expect("sample");
    let expected = [
        (0, 0, truth[0][0], gci[0][0]),
        (0, 1, truth[1][1], 0.5 * (gci[1][1] + gci[1][2])),
        (1, 1, truth[2][3], gci[2][3]),
    ];
    for (r, s, z, tolerance) in expected {
        let w = law.pair_mean(r, s);
        let (mass, real) = (loading.added_mass[r * 2 + s], loading.resistance[r * 2 + s]);
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"two-sphere-pair-mean-law\",\"entry\":[{r},{s}],\"added_mass\":[{mass},{}],\"resistance\":[{real},{}],\"gci\":{tolerance}}}",
            -z.im / w,
            z.re
        );
        assert!(
            (mass - -z.im / w).abs() <= tolerance / w,
            "M_a ({r}, {s}): {mass} against {}",
            -z.im / w
        );
        assert!(
            (real - z.re).abs() <= tolerance,
            "R ({r}, {s}): {real} against {}",
            z.re
        );
    }
    // A pulsating sphere's own added mass and resistance are positive: -Im Z / w and Re Z of
    // rho c S (x^2 - i x) / (1 + x^2).
    assert!(loading.added_mass[0] > 0.0 && loading.resistance[0] > 0.0);
    assert!(loading.added_mass[3] > 0.0 && loading.resistance[3] > 0.0);
}
