//! Independent four-equation spherical-mode continuity with complex k and
//! real density. Manufactured test media, no BEM kernels or fs-la factorization.
use fs_math::c64::C64;

pub struct MolecularShellOracle {
    omega: f64,
    rho: [f64; 3],
    k: [C64; 3],
    modes: Vec<[C64; 4]>,
}
impl MolecularShellOracle {
    pub fn new(omega: f64, radii: [f64; 2], rho: [f64; 3], k: [C64; 3], order: usize) -> Self {
        let mut modes = Vec::new();
        for l in 0..=order {
            let [jo, yo, djo, dyo] = radial(l, k[0].scale(radii[0]));
            let [js, ys, djs, dys] = radial(l, k[1].scale(radii[0]));
            let [ja, ya, dja, dya] = radial(l, k[1].scale(radii[1]));
            let [jc, _, djc, _] = radial(l, k[2].scale(radii[1]));
            let imaginary = C64::new(0.0, 1.0);
            let (ho, dho) = (jo + imaginary * yo, djo + imaginary * dyo);
            let m = std::array::from_fn::<_, 3, _>(|r| k[r].scale(1.0 / rho[r]));
            // p0=j0+a*h0, p1=b*j1+c*y1, p2=d*j2. Match p and
            // dp/dr divided by REAL rho on both physical interfaces.
            modes.push(solve_four([
                [ho, -js, -ys, C64::ZERO],
                [m[0] * dho, -m[1] * djs, -m[1] * dys, C64::ZERO],
                [C64::ZERO, ja, ya, -jc],
                [C64::ZERO, m[1] * dja, m[1] * dya, -m[2] * djc],
            ], [-jo, -m[0] * djo, C64::ZERO, C64::ZERO]));
        }
        Self { omega, rho, k, modes }
    }

    pub fn axial(&self, region: usize, radius: f64) -> (C64, C64) {
        let (mut p, mut dp, mut phase) = (C64::ZERO, C64::ZERO, C64::ONE);
        for l in 0..self.modes.len() {
            let (value, derivative) = self.field(l, region, radius);
            let weight = phase.scale((2 * l + 1) as f64);
            p = p + weight * value;
            dp = dp + weight * derivative;
            phase = phase * C64::new(0.0, 1.0);
        }
        (p, dp / C64::new(0.0, self.omega * self.rho[region]))
    }

    pub fn outward_work(&self, region: usize, radius: f64) -> f64 {
        // Orthogonality integrates the complete plane-wave modal field on
        // the sphere, with the same peak-phasor one-half as p*conj(v).
        let mut power = 0.0;
        for l in 0..self.modes.len() {
            let (p, dp) = self.field(l, region, radius);
            let v = dp / C64::new(0.0, self.omega * self.rho[region]);
            power += (2 * l + 1) as f64 * (p * v.conj()).re;
        }
        2.0 * core::f64::consts::PI * radius * radius * power
    }

    pub fn continuity_defect(&self, radii: [f64; 2]) -> f64 {
        let mut error = 0.0_f64;
        for l in 0..self.modes.len() {
            for boundary in 0..2 {
                let (a, da) = self.field(l, boundary, radii[boundary]);
                let (b, db) = self.field(l, boundary + 1, radii[boundary]);
                error = error.max((a - b).abs() / (1.0 + a.abs() + b.abs()));
                let (va, vb) = (da.scale(1.0 / self.rho[boundary]), db.scale(1.0 / self.rho[boundary + 1]));
                error = error.max((va - vb).abs() / (1.0 + va.abs() + vb.abs()));
            }
        }
        error
    }

    fn field(&self, l: usize, region: usize, radius: f64) -> (C64, C64) {
        let [a, b, c, d] = self.modes[l];
        let [j, y, dj, dy] = radial(l, self.k[region].scale(radius));
        let imaginary = C64::new(0.0, 1.0);
        let (p, dp) = match region {
            0 => (j + a * (j + imaginary * y), dj + a * (dj + imaginary * dy)),
            1 => (b * j + c * y, b * dj + c * dy),
            2 => (d * j, d * dj),
            _ => panic!("unknown region"),
        };
        (p, dp * self.k[region])
    }
}

fn spherical_j(l: usize, z: C64) -> C64 {
    let mut term = C64::ONE;
    for n in 1..=l { term = (term * z).scale(1.0 / (2 * n + 1) as f64); }
    let mut sum = term;
    for n in 1..=80 {
        term = -(term * z * z).scale(1.0 / (2 * n * (2 * l + 2 * n + 1)) as f64);
        sum = sum + term;
        if term.abs() <= 2e-17 * sum.abs() { break; }
    }
    sum
}
fn spherical_y(l: usize, z: C64) -> C64 {
    let sine = C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh());
    let cosine = C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh());
    let mut previous = -cosine / z;
    if l == 0 { return previous; }
    let mut current = -cosine / (z * z) - sine / z;
    for n in 1..l {
        let next = current.scale((2 * n + 1) as f64) / z - previous;
        previous = current;
        current = next;
    }
    current
}
fn radial(l: usize, z: C64) -> [C64; 4] {
    let (j, y) = (spherical_j(l, z), spherical_y(l, z));
    [j, y, j.scale(l as f64) / z - spherical_j(l + 1, z), y.scale(l as f64) / z - spherical_y(l + 1, z)]
}
fn solve_four(mut a: [[C64; 4]; 4], mut b: [C64; 4]) -> [C64; 4] {
    for row in 0..4 {
        let scale = a[row].iter().map(|z| z.abs()).fold(0.0, f64::max);
        for value in &mut a[row] { *value = value.scale(1.0 / scale); }
        b[row] = b[row].scale(1.0 / scale);
    }
    for column in 0..4 {
        let pivot = (column..4).max_by(|&a_row, &b_row| a[a_row][column].abs().total_cmp(&a[b_row][column].abs())).unwrap();
        assert!(a[pivot][column].abs() > 1e-250);
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..4 {
            let factor = a[row][column] / a[column][column];
            for col in column + 1..4 { a[row][col] = a[row][col] - factor * a[column][col]; }
            b[row] = b[row] - factor * b[column];
        }
    }
    let mut x = [C64::ZERO; 4];
    for row in (0..4).rev() {
        let mut value = b[row];
        for col in row + 1..4 { value = value - a[row][col] * x[col]; }
        x[row] = value / a[row][row];
    }
    x
}
