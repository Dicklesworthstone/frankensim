//! Independent axisymmetric two-sphere continuum solution.
//! Exact pressure/superficial-velocity matching gives local partial waves.
//! The axial addition theorem translates outgoing waves (checked separately
//! against Gauss-Legendre projection), and a separate small
//! pivoted elimination sums all interactions. No BEM kernels, panel matrices,
//! material formulas, production LU or receiver rows are used here.
use super::*;
use crate::helmholtz::bc::IncidentSample;

fn regular(z: C64, count: usize) -> Vec<C64> {
    let mut result = Vec::with_capacity(count);
    let mut leading = C64::ONE;
    for n in 0..count {
        if n > 0 {
            leading = leading * z.scale(1.0 / (2 * n + 1) as f64);
        }
        let (mut term, mut sum) = (leading, leading);
        for m in 1..48 {
            term = term * (z * z).scale(-1.0 / (2 * m * (2 * n + 2 * m + 1)) as f64);
            sum = sum + term;
        }
        result.push(sum);
    }
    result
}

fn exp_i(z: C64) -> C64 {
    C64::new(z.re.cos(), z.re.sin()).scale((-z.im).exp())
}

fn outgoing(z: C64, count: usize) -> Vec<C64> {
    let e = exp_i(z) / z;
    let mut result = vec![
        C64::new(0.0, -1.0) * e,
        -e * (C64::ONE + C64::new(0.0, 1.0) / z),
    ];
    for n in 1..count - 1 {
        result.push(result[n].scale((2 * n + 1) as f64) / z - result[n - 1]);
    }
    result
}

fn derivative(f: &[C64], z: C64, n: usize) -> C64 {
    f[n].scale(n as f64) / z - f[n + 1]
}

fn legendre(x: f64, n: usize) -> Vec<f64> {
    let mut p = vec![1.0; n];
    if n > 1 {
        p[1] = x;
    }
    for l in 1..n - 1 {
        p[l + 1] = ((2 * l + 1) as f64 * x * p[l] - l as f64 * p[l - 1]) / (l + 1) as f64;
    }
    p
}

fn gauss(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            for _ in 0..20 {
                let p = legendre(x, n + 1);
                let dp = n as f64 * (x * p[n] - p[n - 1]) / (x * x - 1.0);
                let step = p[n] / dp;
                x -= step;
                if step.abs() < 2e-16 {
                    break;
                }
            }
            let p = legendre(x, n + 1);
            let dp = n as f64 * (x * p[n] - p[n - 1]) / (x * x - 1.0);
            (x, 2.0 / ((1.0 - x * x) * dp * dp))
        })
        .collect()
}

fn solve(mut a: Vec<C64>, mut b: Vec<C64>) -> Vec<C64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        for j in 0..n {
            a.swap(k * n + j, p * n + j);
        }
        b.swap(k, p);
        for i in k + 1..n {
            let f = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                a[i * n + j] = a[i * n + j] - f * a[k * n + j];
            }
            b[i] = b[i] - f * b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] = b[i] - a[i * n + j] * b[j];
        }
        b[i] = b[i] / a[i * n + i];
    }
    b
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |a, i| a * i as f64)
}

// Axial spherical-wave addition theorem. The squared zero-order Wigner 3-j
// symbol has a finite factorial expression; no oscillatory projection loses
// digits when the target regular wave is very small at high order.
fn translate(target: usize, source: usize, separation_k: f64) -> C64 {
    let h = outgoing(C64::from_re(separation_k.abs()), target + source + 2);
    let mut result = C64::ZERO;
    for order in target.abs_diff(source)..=target + source {
        if (target + source + order) % 2 != 0 {
            continue;
        }
        let g = (target + source + order) / 2;
        let ratio =
            factorial(g) / (factorial(g - target) * factorial(g - source) * factorial(g - order));
        let wigner2 =
            factorial(2 * (g - target)) * factorial(2 * (g - source)) * factorial(2 * (g - order))
                / factorial(2 * g + 1)
                * ratio
                * ratio;
        let phase_power = target as isize - source as isize + order as isize;
        let mut sign = if phase_power.rem_euclid(4) == 0 {
            1.0
        } else {
            -1.0
        };
        if separation_k < 0.0 && order % 2 == 1 {
            sign = -sign;
        }
        result =
            result + h[order].scale(sign * ((2 * target + 1) * (2 * order + 1)) as f64 * wigner2);
    }
    result
}

#[test]
fn axial_addition_matches_independent_gauss_projection() {
    for separation in [-1.5, 1.5] {
        let radius_k = 0.5;
        let j = regular(C64::from_re(radius_k), 6);
        for target in 0..5 {
            for source in 0..5 {
                let mut projected = C64::ZERO;
                for (mu, w) in gauss(80) {
                    let dz = separation + radius_k * mu;
                    let r = (radius_k * radius_k * (1.0 - mu * mu) + dz * dz).sqrt();
                    let h = outgoing(C64::from_re(r), 7);
                    projected = projected
                        + h[source].scale(
                            legendre(dz / r, 6)[source]
                                * legendre(mu, 6)[target]
                                * w
                                * (2 * target + 1) as f64
                                / 2.0,
                        ) / j[target];
                }
                let exact = translate(target, source, separation);
                assert!(
                    (exact - projected).abs() < 3e-8 * exact.abs(),
                    "target={target} source={source}: {exact:?} != {projected:?}"
                );
            }
        }
    }
}

pub(super) struct Pair {
    radius: f64,
    centers: [f64; 2],
    ko: C64,
    ki: [C64; 2],
    rho_o: f64,
    rho_i: [C64; 2],
    bulk: [C64; 2],
    omega: f64,
    order: usize,
    outgoing: Vec<C64>,
    transmitted: Vec<C64>,
    pub(super) inward_power: [f64; 2],
}
impl Pair {
    pub(super) fn new(
        gas: &GasState,
        response: [PorousFluidResponse; 2],
        radius: f64,
        centers: [f64; 2],
        order: usize,
        coupled: bool,
    ) -> Self {
        let omega = response[0].angular_frequency_rad_s();
        assert_eq!(omega, response[1].angular_frequency_rad_s());
        let ko = C64::from_re(omega / gas.sound_speed);
        let ki = response.map(|r| r.characteristic().1.conj());
        let rho_i = response.map(|r| r.dynamic_density_kg_m3().conj());
        let bulk = response.map(|r| r.dynamic_bulk_modulus_pa().conj());
        let n = 2 * order;
        let mut translation = vec![C64::ZERO; n * n];
        let mut scattering = vec![C64::ZERO; n];
        let mut transmission = vec![C64::ZERO; n];
        let mut incoming = vec![C64::ZERO; n];
        let zo = ko.scale(radius);
        let jo = regular(zo, order + 1);
        let ho = outgoing(zo, order + 1);
        for body in 0..2 {
            let zi = ki[body].scale(radius);
            let ji = regular(zi, order + 1);
            let eta = rho_i[body].scale(1.0 / gas.density);
            let mut phase = exp_i(ko.scale(centers[body]));
            for l in 0..order {
                let index = body * order + l;
                let djo = derivative(&jo, zo, l);
                let dho = derivative(&ho, zo, l);
                let dji = derivative(&ji, zi, l);
                let denominator = ko * ji[l] * dho - ki[body] * dji * ho[l] / eta;
                scattering[index] = (ki[body] * dji * jo[l] / eta - ko * ji[l] * djo) / denominator;
                transmission[index] =
                    (C64::new(0.0, 1.0) / ko.scale(radius * radius)) / denominator;
                incoming[index] = phase.scale((2 * l + 1) as f64);
                phase = phase * C64::new(0.0, 1.0);
            }
            if coupled {
                let source = 1 - body;
                let distance = centers[body] - centers[source];
                for l in 0..order {
                    for m in 0..order {
                        translation[(body * order + l) * n + source * order + m] =
                            translate(l, m, ko.re * distance);
                    }
                }
            }
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut drive = vec![C64::ZERO; n];
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] =
                    -scattering[i] * ho[i % order] * translation[i * n + j] / ho[j % order];
            }
            matrix[i * n + i] = C64::ONE;
            drive[i] = scattering[i] * ho[i % order] * incoming[i];
        }
        // Scale unknowns by their outgoing boundary amplitude, avoiding the
        // large high-order Hankel factors in an unscaled multipole solve.
        let outgoing = solve(matrix, drive)
            .into_iter()
            .enumerate()
            .map(|(i, value)| value / ho[i % order])
            .collect::<Vec<_>>();
        let transmitted = (0..n)
            .map(|i| {
                let mut exciting = incoming[i];
                for j in 0..n {
                    exciting = exciting + translation[i * n + j] * outgoing[j];
                }
                transmission[i] * exciting
            })
            .collect::<Vec<_>>();
        let mut inward_power = [0.0; 2];
        for body in 0..2 {
            let z = ki[body].scale(radius);
            let j = regular(z, order + 1);
            for l in 0..order {
                let b = transmitted[body * order + l];
                let p = b * j[l];
                let v = b * ki[body] * derivative(&j, z, l) / (C64::new(0.0, omega) * rho_i[body]);
                inward_power[body] -=
                    2.0 * core::f64::consts::PI * radius * radius * (p * v.conj()).re
                        / (2 * l + 1) as f64;
            }
        }
        Self {
            radius,
            centers,
            ko,
            ki,
            rho_o: gas.density,
            rho_i,
            bulk,
            omega,
            order,
            outgoing,
            transmitted,
            inward_power,
        }
    }

    pub(super) fn exterior_axis(&self, z: f64) -> IncidentSample {
        let incident = exp_i(self.ko.scale(z));
        let (mut p, mut gradient) = (incident, C64::new(0.0, 1.0) * self.ko * incident);
        for body in 0..2 {
            let dz = z - self.centers[body];
            let (r, direction) = (dz.abs(), dz.signum());
            let kr = self.ko.scale(r);
            let h = outgoing(kr, self.order + 1);
            let angular = legendre(direction, self.order);
            for l in 0..self.order {
                let b = self.outgoing[body * self.order + l].scale(angular[l]);
                p = p + b * h[l];
                gradient = gradient + b * self.ko * derivative(&h, kr, l).scale(direction);
            }
        }
        IncidentSample {
            pressure_pa: p,
            velocity_m_s: [
                C64::ZERO,
                C64::ZERO,
                gradient / C64::new(0.0, self.omega * self.rho_o),
            ],
        }
    }
    pub(super) fn interior_axis(&self, body: usize, r: f64) -> IncidentSample {
        assert!(r > 0.0 && r < self.radius);
        let z = self.ki[body].scale(r);
        let j = regular(z, self.order + 1);
        let (mut p, mut gradient) = (C64::ZERO, C64::ZERO);
        for l in 0..self.order {
            let b = self.transmitted[body * self.order + l];
            p = p + b * j[l];
            gradient = gradient + b * self.ki[body] * derivative(&j, z, l);
        }
        IncidentSample {
            pressure_pa: p,
            velocity_m_s: [
                C64::ZERO,
                C64::ZERO,
                gradient / (C64::new(0.0, self.omega) * self.rho_i[body]),
            ],
        }
    }
    pub(super) fn volume_losses(&self, body: usize, count: usize) -> [f64; 2] {
        let mut result = [0.0; 2];
        for (x, w) in gauss(count) {
            let r = 0.5 * self.radius * (x + 1.0);
            let z = self.ki[body].scale(r);
            let j = regular(z, self.order + 1);
            let (mut p2, mut gradient2) = (0.0, 0.0);
            for l in 0..self.order {
                let weight = self.transmitted[body * self.order + l].norm_sq() / (2 * l + 1) as f64;
                p2 += weight * j[l].norm_sq();
                gradient2 += weight
                    * ((self.ki[body] * derivative(&j, z, l)).norm_sq()
                        + (l * (l + 1)) as f64 * j[l].norm_sq() / (r * r));
            }
            let volume = 2.0 * core::f64::consts::PI * self.radius * w * r * r;
            result[0] += volume * 0.5 * self.omega * self.rho_i[body].im * gradient2
                / self.rho_i[body].scale(self.omega).norm_sq();
            result[1] += volume * 0.5 * self.omega * self.bulk[body].recip().im * p2;
        }
        result
    }
}
