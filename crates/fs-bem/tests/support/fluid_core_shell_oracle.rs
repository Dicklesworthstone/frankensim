//! Independent concentric three-fluid oracle: four continuity equations per
//! spherical harmonic, with no BEM kernels, matrices, or fs-la factorization.
//! Numerical coefficients are explicit test fixtures, not measured fluids.
use fs_math::c64::C64;

pub struct CoreShellOracle {
    omega: f64,
    radii: [f64; 2],
    rho: [f64; 3],
    k: [f64; 3],
    modes: Vec<[C64; 4]>,
}

impl CoreShellOracle {
    pub fn new(omega: f64, radii: [f64; 2], rho: [f64; 3], c: [f64; 3], order: usize) -> Self {
        let k = c.map(|speed| omega / speed);
        let mut modes = Vec::new();
        for l in 0..=order {
            let (jo, djo, ho, dho) = radial(l, k[0] * radii[0]);
            let (js, djs, hs, dhs) = radial(l, k[1] * radii[0]);
            let (ja, dja, ha, dha) = radial(l, k[1] * radii[1]);
            let (jc, djc, _, _) = radial(l, k[2] * radii[1]);
            let momentum = std::array::from_fn::<_, 3, _>(|r| k[r] / rho[r]);
            // p0=j0+a*h0, p1=b*j1+c*y1, p2=d*j2. Match pressure
            // and (1/rho)*dp/dr on the two original spherical boundaries.
            let matrix = [
                [ho, -C64::from_re(js), -C64::from_re(hs.im), C64::ZERO],
                [
                    dho.scale(momentum[0]),
                    C64::from_re(-momentum[1] * djs),
                    C64::from_re(-momentum[1] * dhs.im),
                    C64::ZERO,
                ],
                [
                    C64::ZERO,
                    C64::from_re(ja),
                    C64::from_re(ha.im),
                    C64::from_re(-jc),
                ],
                [
                    C64::ZERO,
                    C64::from_re(momentum[1] * dja),
                    C64::from_re(momentum[1] * dha.im),
                    C64::from_re(-momentum[2] * djc),
                ],
            ];
            let rhs = [
                C64::from_re(-jo),
                C64::from_re(-momentum[0] * djo),
                C64::ZERO,
                C64::ZERO,
            ];
            modes.push(solve_four(matrix, rhs));
        }
        Self {
            omega,
            radii,
            rho,
            k,
            modes,
        }
    }

    /// Outgoing amplitude for a unit plane wave traveling along positive z:
    /// p_scattered(r*d) ~ exp(i*k0*r)/r * far(d_z).
    pub fn far(&self, cosine: f64) -> C64 {
        let mut result = C64::ZERO;
        let mut previous = 0.0;
        let mut legendre = 1.0;
        for (l, coefficients) in self.modes.iter().enumerate() {
            result = result
                + coefficients[0] * C64::new(0.0, -((2 * l + 1) as f64) * legendre / self.k[0]);
            let next =
                ((2 * l + 1) as f64 * cosine * legendre - l as f64 * previous) / (l + 1) as f64;
            previous = legendre;
            legendre = next;
        }
        result
    }

    /// Total pressure and positive-z particle velocity on the positive axis.
    /// The test points avoid the origin and either material boundary.
    pub fn axial(&self, region: usize, radius: f64) -> (C64, C64) {
        let mut pressure = C64::ZERO;
        let mut derivative = C64::ZERO;
        let mut phase = C64::ONE;
        for l in 0..self.modes.len() {
            let (p, dp) = self.radial_field(l, region, radius);
            let weight = phase.scale((2 * l + 1) as f64);
            pressure = pressure + weight * p;
            derivative = derivative + weight * dp;
            phase = phase * C64::new(0.0, 1.0);
        }
        (
            pressure,
            derivative / C64::new(0.0, self.omega * self.rho[region]),
        )
    }

    fn radial_field(&self, l: usize, region: usize, radius: f64) -> (C64, C64) {
        let [a, b, c, d] = self.modes[l];
        let (j, dj, h, dh) = radial(l, self.k[region] * radius);
        let (p, dp) = match region {
            0 => (C64::from_re(j) + a * h, C64::from_re(dj) + a * dh),
            1 => (b.scale(j) + c.scale(h.im), b.scale(dj) + c.scale(dh.im)),
            2 => (d.scale(j), d.scale(dj)),
            _ => panic!("unknown fluid region"),
        };
        (p, dp.scale(self.k[region]))
    }

    pub fn maximum_continuity_defect(&self) -> f64 {
        let mut error = 0.0_f64;
        for l in 0..self.modes.len() {
            for boundary in 0..2 {
                let (a, da) = self.radial_field(l, boundary, self.radii[boundary]);
                let (b, db) = self.radial_field(l, boundary + 1, self.radii[boundary]);
                let va = da.scale(1.0 / self.rho[boundary]);
                let vb = db.scale(1.0 / self.rho[boundary + 1]);
                error = error.max((a - b).abs() / (1.0 + a.abs() + b.abs()));
                error = error.max((va - vb).abs() / (1.0 + va.abs() + vb.abs()));
            }
        }
        error
    }

    pub fn maximum_lossless_defect(&self) -> f64 {
        self.modes
            .iter()
            .map(|mode| {
                let a = mode[0];
                (a.re + a.abs() * a.abs()).abs()
            })
            .fold(0.0, f64::max)
    }
}

fn spherical_j(l: usize, x: f64) -> f64 {
    // Power series is stable for the low-frequency test domain, including
    // harmonics far above x where upward recurrence of j_l is unstable.
    let mut term = 1.0;
    for n in 1..=l {
        term *= x / (2 * n + 1) as f64;
    }
    let mut result = term;
    for n in 1..=80 {
        term *= -x * x / (2 * n * (2 * l + 2 * n + 1)) as f64;
        result += term;
        if term.abs() <= 2e-17 * result.abs() {
            break;
        }
    }
    result
}

fn spherical_y(l: usize, x: f64) -> f64 {
    let mut previous = -x.cos() / x;
    if l == 0 {
        return previous;
    }
    let mut current = -x.cos() / (x * x) - x.sin() / x;
    for n in 1..l {
        let next = (2 * n + 1) as f64 * current / x - previous;
        previous = current;
        current = next;
    }
    current
}

fn radial(l: usize, x: f64) -> (f64, f64, C64, C64) {
    let j = spherical_j(l, x);
    let y = spherical_y(l, x);
    let dj = l as f64 * j / x - spherical_j(l + 1, x);
    let dy = l as f64 * y / x - spherical_y(l + 1, x);
    (j, dj, C64::new(j, y), C64::new(dj, dy))
}

fn solve_four(mut a: [[C64; 4]; 4], mut b: [C64; 4]) -> [C64; 4] {
    for row in 0..4 {
        let scale = a[row].iter().map(|z| z.abs()).fold(0.0, f64::max);
        for value in &mut a[row] {
            *value = value.scale(1.0 / scale);
        }
        b[row] = b[row].scale(1.0 / scale);
    }
    for column in 0..4 {
        let pivot = (column..4)
            .max_by(|&a_row, &b_row| a[a_row][column].abs().total_cmp(&a[b_row][column].abs()))
            .unwrap();
        assert!(
            a[pivot][column].abs() > 1e-250,
            "singular independent shell oracle"
        );
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..4 {
            let factor = a[row][column] / a[column][column];
            for col in column + 1..4 {
                a[row][col] = a[row][col] - factor * a[column][col];
            }
            b[row] = b[row] - factor * b[column];
        }
    }
    let mut x = [C64::ZERO; 4];
    for row in (0..4).rev() {
        let mut value = b[row];
        for col in row + 1..4 {
            value = value - a[row][col] * x[col];
        }
        x[row] = value / a[row][row];
    }
    x
}

#[test]
fn independent_shell_oracle_conserves_work_and_converges_in_order() {
    let make =
        |order| CoreShellOracle::new(0.6, [1.0, 0.45], [1.0, 1.6, 0.65], [1.0, 0.75, 1.3], order);
    let (low, high) = (make(14), make(18));
    assert!(high.maximum_continuity_defect() < 2e-14);
    assert!(high.maximum_lossless_defect() < 2e-14);
    for direction in [-1.0, 0.0, 1.0] {
        assert!((low.far(direction) - high.far(direction)).abs() < 2e-14);
    }
    for (region, r) in [(0, 1.8), (1, 0.65), (2, 0.15)] {
        let (pa, va) = low.axial(region, r);
        let (pb, vb) = high.axial(region, r);
        assert!(
            (pa - pb).abs() < 2e-13,
            "region {region}: pressure truncation {}",
            (pa - pb).abs()
        );
        assert!(
            (va - vb).abs() < 2e-13,
            "region {region}: velocity truncation {}",
            (va - vb).abs()
        );
    }
}

#[test]
fn independent_shell_oracle_recovers_unchanged_plane_wave() {
    let oracle = CoreShellOracle::new(0.6, [1.0, 0.45], [1.0; 3], [1.0; 3], 14);
    for direction in [-1.0, 0.0, 1.0] {
        assert!(oracle.far(direction).abs() < 2e-14);
    }
    for (region, r) in [(0, 1.8), (1, 0.65), (2, 0.15)] {
        let expected = C64::new((0.6_f64 * r).cos(), (0.6_f64 * r).sin());
        let (p, v) = oracle.axial(region, r);
        assert!((p - expected).abs() < 2e-13);
        assert!((v - expected).abs() < 2e-13);
    }
}
