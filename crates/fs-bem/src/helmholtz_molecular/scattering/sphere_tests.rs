//! G1/G2: relaxing impedance-sphere scattering against independent partial waves.
use super::MolecularScatteringOptions;
use crate::helmholtz::{
    Formulation,
    bc::{Admittance, IncidentField, PlaneWave},
    laplace::complex::ComplexOptions,
    molecular::MolecularExterior,
};
use crate::panel3d::SpherePanels;
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::{c64::C64, det};
use std::f64::consts::PI;

fn material() -> MolecularExterior {
    // Deliberately synthetic slow reservoir, not an atmospheric-loss fit.
    let data = "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

fn exp_i(z: C64) -> C64 {
    C64::new(det::cos(z.re), det::sin(z.re)).scale(det::exp(-z.im))
}

// NIST DLMF 10.53.1: a power series, avoiding unstable upward j_l recurrence
// for small complex ka. https://dlmf.nist.gov/10.53.E1
fn regular(z: C64, count: usize) -> Vec<C64> {
    let mut leading = C64::ONE;
    let mut values = Vec::with_capacity(count);
    for n in 0..count {
        if n > 0 {
            leading = leading * z.scale(1.0 / (2 * n + 1) as f64);
        }
        let (mut term, mut sum) = (leading, leading);
        for m in 1..48 {
            term = term * (z * z).scale(-1.0 / (2 * m * (2 * n + 2 * m + 1)) as f64);
            sum = sum + term;
        }
        values.push(sum);
    }
    values
}

// Outgoing exp(ikr) seeds plus NIST DLMF 10.51.1, valid for complex z.
// This does not call the BEM Green function or its derivatives.
fn outgoing(z: C64, count: usize) -> Vec<C64> {
    let iz = C64::new(0.0, 1.0) / z;
    let ez = exp_i(z) / z;
    let mut values = vec![-C64::new(0.0, 1.0) * ez, -ez * (C64::ONE + iz)];
    for n in 1..count - 1 {
        values.push(values[n].scale((2 * n + 1) as f64) / z - values[n - 1]);
    }
    values
}

fn derivative(values: &[C64], z: C64, n: usize) -> C64 {
    // NIST DLMF 10.51.2, with f'_0=-f_1.
    if n == 0 {
        -values[1]
    } else {
        values[n - 1] - values[n].scale((n + 1) as f64) / z
    }
}

struct SphereOracle {
    k: C64,
    rho_omega: f64,
    scattered: Vec<C64>,
    absorbed_power_w: f64,
}

impl SphereOracle {
    #[allow(clippy::too_many_arguments)]
    fn new(
        k: C64,
        radius: f64,
        rho_omega: f64,
        rho_c: f64,
        omega_over_c: f64,
        beta: C64,
        amplitude: C64,
        count: usize,
    ) -> Self {
        let ka = k.scale(radius);
        let j = regular(ka, count);
        let h = outgoing(ka, count);
        let mut phase = C64::ONE;
        let mut scattered = Vec::with_capacity(count);
        let mut surface_norm = 0.0;
        // DLMF 10.60.7: exp(ikr cos(theta)) = sum (2l+1)i^l j_l P_l.
        // Enforce d_r p_total + i*rho*omega*Y*p_total = 0 independently.
        // rho*omega*Y = (omega/c_eq)*beta, NOT k*beta in relaxing gas.
        let wall = C64::new(0.0, omega_over_c) * beta;
        for n in 0..count {
            let incident = amplitude * phase.scale((2 * n + 1) as f64);
            let a = -incident * (k * derivative(&j, ka, n) + wall * j[n])
                / (k * derivative(&h, ka, n) + wall * h[n]);
            scattered.push(a);
            let total = incident * j[n] + a * h[n];
            surface_norm += total.abs().powi(2) / (2 * n + 1) as f64;
            phase = phase * C64::new(0.0, 1.0);
        }
        // Legendre orthogonality gives the exact spherical surface integral;
        // no panel quadrature or BEM trace is used to predict wall dissipation.
        let absorbed_power_w = 2.0 * PI * radius * radius * beta.re / rho_c * surface_norm;
        Self {
            k,
            rho_omega,
            scattered,
            absorbed_power_w,
        }
    }

    fn sample(&self, radius: f64, mu: f64) -> (C64, C64) {
        let kr = self.k.scale(radius);
        let h = outgoing(kr, self.scattered.len());
        let (mut previous, mut legendre) = (1.0, 1.0);
        let (mut pressure, mut gradient) = (C64::ZERO, C64::ZERO);
        for (n, &a) in self.scattered.iter().enumerate() {
            if n == 1 {
                legendre = mu;
            } else if n > 1 {
                let next =
                    ((2 * n - 1) as f64 * mu * legendre - (n - 1) as f64 * previous) / n as f64;
                previous = legendre;
                legendre = next;
            }
            pressure = pressure + (a * h[n]).scale(legendre);
            gradient = gradient + (a * self.k * derivative(&h, kr, n)).scale(legendre);
        }
        (pressure, gradient * C64::new(0.0, -1.0 / self.rho_omega))
    }
}

#[test]
fn molecular_robin_scattering_converges_in_pressure_velocity_and_absorption() {
    let model = material();
    let omega = 1000.0;
    let radius = 0.2;
    let distance = 0.6;
    let gas = model.source().gas();
    let k = model.source().relaxing_wavenumber(omega).unwrap();
    let rho_omega = gas.density * omega;
    let rho_c = gas.density * gas.sound_speed;
    let beta = C64::new(0.45, 0.2);
    let amplitude = C64::new(1.0, -0.3);
    assert!(
        k.im > 0.01 * k.re,
        "fixture must resolve molecular absorption"
    );
    assert!((k.scale(gas.sound_speed / omega) - C64::ONE).abs() > 0.02);
    let oracle = |count| {
        SphereOracle::new(
            k,
            radius,
            rho_omega,
            rho_c,
            omega / gas.sound_speed,
            beta,
            amplitude,
            count,
        )
    };
    let exact = oracle(24);
    let fewer = oracle(16);
    let points = [
        [0.0, 0.0, distance],
        [distance, 0.0, 0.0],
        [0.0, 0.0, -distance],
    ];
    let expected = points.map(|p| exact.sample(distance, p[2] / distance));
    for (point, &(p, u)) in points.iter().zip(&expected) {
        let (p16, u16) = fewer.sample(distance, point[2] / distance);
        assert!((p16 - p).abs() < 1e-12 * p.abs());
        assert!((u16 - u).abs() < 1e-12 * u.abs());
    }
    assert!((fewer.absorbed_power_w / exact.absorbed_power_w - 1.0).abs() < 1e-12);
    let incident = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: amplitude,
    });
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let n = surface.areas().len();
        let betas = vec![beta; n];
        let options = MolecularScatteringOptions {
            boundary: ComplexOptions {
                max_panels: n,
                max_bytes: 32 * n * n,
                max_work: 512 * (n * n) as u64,
            },
            minimum_source_clearance_m: 1e-4,
            maximum_source_point_products: 100_000,
        };
        let field = model
            .solve_scattering_batch(
                &surface,
                omega,
                &[incident],
                Admittance::Local {
                    beta: &betas,
                    declared_active: false,
                },
                Formulation::PlainCbie,
                &options,
                || false,
            )
            .unwrap()
            .remove(0);
        let observed = field
            .total_pressure_velocity_at(&points, 0.05, || false)
            .unwrap();
        let (mut dp2, mut p2, mut du2, mut u2) = (0.0, 0.0, 0.0, 0.0);
        for ((sample, point), &(p, u)) in observed.iter().zip(&points).zip(&expected) {
            let incident_p = amplitude * exp_i(k.scale(point[2]));
            let scattered_p = sample.pressure_pa - incident_p;
            let total_radial_u = sample
                .velocity_m_s
                .iter()
                .zip(point)
                .fold(C64::ZERO, |sum, (&v, &x)| sum + v.scale(x / distance));
            let incident_radial_u = (incident_p * k).scale(point[2] / (distance * rho_omega));
            dp2 += (scattered_p - p).abs().powi(2);
            p2 += p.abs().powi(2);
            du2 += (total_radial_u - incident_radial_u - u).abs().powi(2);
            u2 += u.abs().powi(2);
        }
        let solution = field.solution();
        let absorbed = solution.absorbed_power_w;
        assert!(absorbed > 0.0);
        assert!((absorbed + solution.total_boundary_work_w).abs() < 1e-11 * absorbed);
        errors.push([
            (dp2 / p2).sqrt(),
            (du2 / u2).sqrt(),
            (absorbed / exact.absorbed_power_w - 1.0).abs(),
        ]);
    }
    eprintln!("molecular sphere [pressure, radial velocity, absorbed power] errors: {errors:?}");
    for quantity in 0..3 {
        assert!(
            errors[1][quantity] < 0.8 * errors[0][quantity],
            "{errors:?}"
        );
        assert!(
            errors[2][quantity] < 0.8 * errors[1][quantity],
            "{errors:?}"
        );
        assert!(errors[2][quantity] < 0.08, "{errors:?}");
    }
}
