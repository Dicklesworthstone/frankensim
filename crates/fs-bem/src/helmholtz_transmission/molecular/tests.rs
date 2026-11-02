use super::*;
use crate::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

pub(super) fn state(
    temperature: f64,
    pressure: f64,
    relaxation_ms: Option<f64>,
) -> RelaxationState {
    // Manufactured slow caloric reservoir; not measured atmospheric data.
    let text = relaxation_ms.map_or_else(
        || "declared 1 synthetic 0".into(),
        |ms| {
            format!(
                "declared 1 synthetic 1 vibration 1 100 1 power {} 293.15 0 200 500",
                101.325 * ms
            )
        },
    );
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}
pub(super) fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    })
}
fn wave(k: C64) -> DampedWaveNumber {
    DampedWaveNumber {
        k,
        phase_origin_m: 0.0,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-30),
        "{a:?} != {b:?}"
    );
}

#[test]
fn complex_contrast_keeps_weak_static_and_extreme_loss_limits() {
    let (x, y, nx, ny) = ([0.0, 0.0, 1.0], [0.0; 3], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);
    for (ko, ki) in [(0.5, 0.7), (1e-8, 1.000001e-8), (0.5, 0.5 + 1e-12)] {
        let actual =
            kernels::difference(wave(C64::from_re(ko)), wave(C64::from_re(ki)), x, nx, y, ny);
        let expected = super::super::kernel_difference(ko, ki, x, nx, y, ny);
        for (a, b) in actual.into_iter().zip(expected) {
            close(a, b, 2e-12);
        }
        let (a, b, c) = kernels::self_contrast(wave(C64::from_re(ko)), wave(C64::from_re(ki)), 0.1);
        let (d, e, f) = super::super::self_contrast(ko, ki, 0.1);
        for (actual, expected) in [(a, d), (b, e), (c, f)] {
            close(actual, expected, 2e-12);
        }
    }
    let ko = wave(C64::new(0.4, 0.03));
    assert_eq!(kernels::difference(ko, ko, x, nx, y, ny), [C64::ZERO; 4]);
    // The interior exponential is below underflow. Each difference must be
    // minus the original exterior kernel, even with enormous midpoint terms.
    let (s, d, dp, n) = crate::helmholtz::kernels(ko, x, nx, y, ny);
    for loss in [800.0, 20_000.0, 1e6] {
        let actual = kernels::difference(ko, wave(C64::new(0.7, loss)), x, nx, y, ny);
        for (a, b) in actual.into_iter().zip([-s, -d, -dp, -n]) {
            close(a, b, 2e-14);
        }
    }
    let ki = wave(ko.k + C64::new(1e-12, -2e-13));
    let step = 1e-4;
    let plus = crate::helmholtz::kernels(wave(ko.k + C64::from_re(step)), x, nx, y, ny);
    let minus = crate::helmholtz::kernels(wave(ko.k - C64::from_re(step)), x, nx, y, ny);
    let expected = [
        plus.0 - minus.0,
        plus.1 - minus.1,
        plus.2 - minus.2,
        plus.3 - minus.3,
    ]
    .map(|z| z.scale(0.5 / step) * (ki.k - ko.k));
    for (a, b) in kernels::difference(ko, ki, x, nx, y, ny)
        .into_iter()
        .zip(expected)
    {
        close(a, b, 1e-7);
    }
}

#[test]
fn empty_reservoirs_recover_native_transmission_and_molecular_identity_is_transparent() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let outside = state(293.15, 101325.0, None);
    let inside = state(353.15, 101325.0, None);
    let omega = 200.0;
    let native = FluidInterface::new(
        &surface,
        outside.gas(),
        inside.gas(),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let molecular = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let a = native.scatter_field(plane(), 0.01, &gate).unwrap();
    let b = molecular.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    for (&actual, &expected) in b
        .total_pressure_pa()
        .iter()
        .zip(a.solution().total_pressure_pa())
    {
        close(actual, expected, 2e-12);
    }
    let velocity_scale = a
        .solution()
        .total_normal_velocity_m_s()
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    for (&actual, &expected) in b
        .total_normal_velocity_m_s()
        .iter()
        .zip(a.solution().total_normal_velocity_m_s())
    {
        assert!((actual - expected).abs() <= 2e-11 * velocity_scale);
    }
    let relaxing = state(293.15, 101325.0, Some(5.0));
    let owner = MolecularFluidInterface::new(
        &surface,
        &relaxing,
        &relaxing,
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert!(owner.k_o.k.im > 0.0);
    assert!(std::ptr::eq(owner.exterior_source(), &relaxing));
    assert!(std::ptr::eq(owner.interior_source(), &relaxing));
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.5, -0.2),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.2, 0.0, 3.0],
            volume_velocity_m3_s: C64::new(0.002, 0.001),
        }),
    ];
    let incident = IncidentField::Sources(&sources);
    let field = owner.scatter_field(incident, 0.01, 1000, &gate).unwrap();
    assert!(
        field
            .scattered_pressure_pa()
            .iter()
            .chain(field.scattered_normal_velocity_m_s())
            .all(|&z| z == C64::ZERO)
    );
    for (point, interior) in [([0.0, 0.0, 2.0], false), ([0.0, 0.0, 0.3], true)] {
        let expected = incident
            .sample_wave(
                owner.k_o,
                relaxing.gas().density * omega,
                None,
                point,
                &mut || false,
            )
            .unwrap();
        let actual = if interior {
            field.interior_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        } else {
            field.total_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        }
        .unwrap()[0];
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.velocity_m_s, expected.velocity_m_s);
    }
}

// Entire regular-Bessel series and outgoing seeds/recurrence, independent of
// the panel kernels. The continuum coefficients enforce p and (1/rho)dp/dr
// continuity mode by mode, with physical SI rho and omega.
pub(super) fn regular(z: C64, count: usize) -> Vec<C64> {
    let mut out = Vec::with_capacity(count);
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
        out.push(sum);
    }
    out
}
pub(super) fn exp_i(z: C64) -> C64 {
    C64::new(z.re.cos(), z.re.sin()).scale((-z.im).exp())
}
pub(super) fn outgoing(z: C64, count: usize) -> Vec<C64> {
    let ez = exp_i(z) / z;
    let mut out = vec![
        C64::new(0.0, -1.0) * ez,
        -ez * (C64::ONE + C64::new(0.0, 1.0) / z),
    ];
    for n in 1..count - 1 {
        out.push(out[n].scale((2 * n + 1) as f64) / z - out[n - 1]);
    }
    out
}
pub(super) fn derivative(f: &[C64], z: C64, n: usize) -> C64 {
    f[n].scale(n as f64) / z - f[n + 1]
}
struct Sphere {
    ko: C64,
    ki: C64,
    rho_o: f64,
    rho_i: f64,
    omega: f64,
    scattered: Vec<C64>,
    transmitted: Vec<C64>,
    inward_power: f64,
}
impl Sphere {
    fn new(outside: &RelaxationState, inside: &RelaxationState, omega: f64) -> Self {
        let (ko, ki) = (
            outside.relaxing_wavenumber(omega).unwrap(),
            inside.relaxing_wavenumber(omega).unwrap(),
        );
        let eta = inside.gas().density / outside.gas().density;
        let (jo, ji, ho) = (regular(ko, 17), regular(ki, 17), outgoing(ko, 17));
        let (mut phase, mut scattered, mut transmitted, mut outward) =
            (C64::ONE, Vec::new(), Vec::new(), 0.0);
        for n in 0..16 {
            let (djo, dji, dho) = (
                derivative(&jo, ko, n),
                derivative(&ji, ki, n),
                derivative(&ho, ko, n),
            );
            let denominator = ko * ji[n] * dho - (ki * dji * ho[n]).scale(1.0 / eta);
            let a = ((ki * dji * jo[n]).scale(1.0 / eta) - ko * ji[n] * djo) / denominator;
            let b = (C64::new(0.0, 1.0) / ko) / denominator;
            let amplitude = phase.scale((2 * n + 1) as f64);
            scattered.push(amplitude * a);
            transmitted.push(amplitude * b);
            let p = amplitude * b * ji[n];
            let v = amplitude * b * ki * dji * C64::new(0.0, -1.0 / (inside.gas().density * omega));
            outward += (p * v.conj()).re / (2 * n + 1) as f64;
            phase = phase * C64::new(0.0, 1.0);
        }
        Self {
            ko,
            ki,
            rho_o: outside.gas().density,
            rho_i: inside.gas().density,
            omega,
            scattered,
            transmitted,
            inward_power: -2.0 * std::f64::consts::PI * outward,
        }
    }
    fn axis(&self, r: f64, interior: bool) -> IncidentSample {
        let (k, rho, coefficients) = if interior {
            (self.ki, self.rho_i, &self.transmitted)
        } else {
            (self.ko, self.rho_o, &self.scattered)
        };
        let z = k.scale(r);
        let f = if interior {
            regular(z, 17)
        } else {
            outgoing(z, 17)
        };
        let (mut p, mut gradient) = (C64::ZERO, C64::ZERO);
        for (n, &c) in coefficients.iter().enumerate() {
            p = p + c * f[n];
            gradient = gradient + c * k * derivative(&f, z, n);
        }
        if !interior {
            p = p + exp_i(z);
            gradient = gradient + C64::new(0.0, 1.0) * k * exp_i(z);
        }
        IncidentSample {
            pressure_pa: p,
            velocity_m_s: [
                C64::ZERO,
                C64::ZERO,
                gradient * C64::new(0.0, -1.0 / (rho * self.omega)),
            ],
        }
    }
    fn volume_loss(&self, inside: &RelaxationState) -> f64 {
        let mut integral = 0.0;
        for i in 0..=256 {
            let r = i as f64 / 256.0;
            let j = regular(self.ki.scale(r), 16);
            let norm: f64 = self
                .transmitted
                .iter()
                .zip(j)
                .enumerate()
                .map(|(n, (&b, j))| (b * j).abs().powi(2) / (2 * n + 1) as f64)
                .sum();
            let weight = if i == 0 || i == 256 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * r * r * norm;
        }
        let bulk = inside.caloric_response(self.omega).unwrap().bulk_modulus_pa;
        2.0 * std::f64::consts::PI * self.omega * bulk.recip().im * integral / (3.0 * 256.0)
    }
}

#[test]
fn molecular_transmission_converges_to_complex_sphere_and_bulk_dissipation() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(293.15, 101325.0, Some(5.0));
    let omega = 200.0;
    let exact = Sphere::new(&outside, &inside, omega);
    assert!(exact.inward_power > 0.0);
    assert!((exact.volume_loss(&inside) - exact.inward_power).abs() < 1e-8 * exact.inward_power);
    let (mut previous_field, mut previous_power) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let owner = MolecularFluidInterface::new(
            &surface,
            &outside,
            &inside,
            omega,
            Options::default(),
            &gate,
        )
        .unwrap();
        let field = owner.scatter_field(plane(), 0.01, 10_000, &gate).unwrap();
        let mut error = 0.0_f64;
        for (radius, interior) in [(0.35, true), (2.0, false)] {
            let point = [0.0, 0.0, radius];
            let actual = if interior {
                field.interior_pressure_velocity_at(&[point], 0.01, 10_000, &gate)
            } else {
                field.total_pressure_velocity_at(&[point], 0.01, 10_000, &gate)
            }
            .unwrap()[0];
            let expected = exact.axis(radius, interior);
            error = error.max(
                (actual.pressure_pa - expected.pressure_pa).abs() / expected.pressure_pa.abs(),
            );
            error = error.max(
                (actual.velocity_m_s[2] - expected.velocity_m_s[2]).abs()
                    / expected.velocity_m_s[2].abs(),
            );
        }
        let power = -field.interface_work.surface_power.re;
        let power_error = (power - exact.inward_power).abs() / exact.inward_power;
        assert!(
            error < previous_field,
            "mesh={level}: field={error}, previous={previous_field}"
        );
        assert!(
            power_error < previous_power,
            "mesh={level}: power={power_error}, previous={previous_power}"
        );
        (previous_field, previous_power) = (error, power_error);
        if level == 2 {
            assert!(error < 0.08, "field={error}");
            assert!(power_error < 0.15, "power={power_error}");
            assert!(power > 0.0);
        }
    }
}

#[test]
fn molecular_interface_refuses_unresolved_sources_resources_and_cancellation() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(293.15, 101325.0, Some(5.0));
    assert!(
        MolecularFluidInterface::new(
            &surface,
            &outside,
            &inside,
            200.0,
            Options {
                maximum_bytes: 0,
                ..Options::default()
            },
            &gate
        )
        .is_err()
    );
    let mismatch = state(293.15, 102325.0, Some(5.0));
    assert!(
        MolecularFluidInterface::new(
            &surface,
            &outside,
            &mismatch,
            200.0,
            Options::default(),
            &gate
        )
        .is_err()
    );
    let owner = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert!(owner.scatter_field(plane(), 0.01, 19, &gate).is_err());
    let inside_source = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&inside_source), 0.01, 1000, &gate)
            .is_err()
    );
    let solved = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(
        solved
            .total_pressure_velocity_at(&[[0.0; 3]], 0.01, 1000, &gate)
            .is_err()
    );
    assert!(
        solved
            .interior_pressure_velocity_at(&[[0.0, 0.0, 2.0]], 0.01, 1000, &gate)
            .is_err()
    );
    assert!(
        solved
            .total_pressure_velocity_at(&[[0.0, 0.0, 2.0]], 0.01, 19, &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        owner.scatter_field(plane(), 0.01, 1000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        solved.interior_pressure_velocity_at(&[[0.0; 3]], 0.01, 1000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
}
