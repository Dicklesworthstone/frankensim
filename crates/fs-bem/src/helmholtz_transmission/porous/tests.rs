use super::*;
use crate::helmholtz::bc::incident::IncidentSample;
use crate::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{Geometry, complex::Prepared};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn law() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}
fn plane(amplitude: C64) -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: amplitude,
    })
}
fn controls() -> MolecularNearFieldOptions {
    MolecularNearFieldOptions {
        quadrature: near_field::Options {
            relative_tolerance: 2e-7,
            maximum_kernel_evaluations: 2_000_000,
            ..near_field::Options::default()
        },
        ..MolecularNearFieldOptions::default()
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-25),
        "{a:?} != {b:?}"
    );
}
fn wave(k: C64) -> DampedWaveNumber {
    DampedWaveNumber {
        k,
        phase_origin_m: 0.0,
    }
}

#[test]
fn complex_density_reduces_to_real_transmission_and_receiver_momentum() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let expected = assemble(&surface, 0.5, 0.7, 1.4, 0.4, &gate).unwrap();
    let actual = assemble_wave(
        &surface,
        wave(C64::from_re(0.5)),
        wave(C64::from_re(0.7)),
        0.5,
        C64::from_re(1.4),
        C64::from_re(0.4),
        &gate,
        molecular::kernels::difference,
        molecular::kernels::self_contrast,
    )
    .unwrap();
    for (&a, &b) in actual.iter().zip(&expected) {
        assert!(
            (a - b).abs() <= 3e-14 * b.abs().max(1e-25),
            "{a:?} != {b:?}"
        );
    }
    let points = [[0.2, 0.1, 2.0]];
    let geometry = Geometry::new(&surface, &points, 0.01).unwrap();
    let k = wave(C64::new(0.5, 0.1));
    let legacy = Prepared::new(&geometry, k, 413.0, controls().quadrature, &mut || false).unwrap();
    let complex = Prepared::new_equivalent_fluid(
        &geometry,
        k,
        C64::from_re(413.0),
        controls().quadrature,
        &mut || false,
    )
    .unwrap();
    let p = vec![C64::new(0.7, -0.2); surface.areas().len()];
    let v = vec![C64::new(0.002, 0.001); p.len()];
    let a = legacy.evaluate(&p, &v, &mut || false).unwrap();
    let b = complex.evaluate(&p, &v, &mut || false).unwrap();
    assert_eq!(a.scalar.pressure, b.scalar.pressure);
    assert_eq!(a.particle_velocity_m_s, b.particle_velocity_m_s);
    assert_eq!(
        a.scalar.quadrature_error_estimate_pa,
        b.scalar.quadrature_error_estimate_pa
    );
    assert_eq!(
        a.quadrature_error_estimate_m_s,
        b.quadrature_error_estimate_m_s
    );
}

// Independent entire spherical-Bessel series and outgoing recurrence. No BEM
// kernels or production transmission algebra supply the continuum oracle.
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
fn gauss(count: usize) -> Vec<(f64, f64)> {
    (0..count)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (count as f64 + 0.5)).cos();
            let mut d = 0.0;
            for _ in 0..32 {
                let (mut a, mut b) = (1.0, x);
                for n in 2..=count {
                    let next = ((2 * n - 1) as f64 * x * b - (n - 1) as f64 * a) / n as f64;
                    (a, b) = (b, next);
                }
                d = count as f64 * (x * b - a) / (x * x - 1.0);
                let step = b / d;
                x -= step;
                if step.abs() < 2e-16 {
                    break;
                }
            }
            (0.5 * (x + 1.0), 1.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}
struct Sphere {
    radius: f64,
    ko: C64,
    ki: C64,
    rho_o: f64,
    rho_i: C64,
    bulk: C64,
    omega: f64,
    scattered: Vec<C64>,
    transmitted: Vec<C64>,
    inward_power: f64,
}
impl Sphere {
    fn new(gas: &GasState, response: PorousFluidResponse, radius: f64) -> Self {
        let omega = response.angular_frequency_rad_s();
        let (ko, ki, rho_i, bulk) = (
            C64::from_re(omega / gas.sound_speed),
            response.characteristic().1.conj(),
            response.dynamic_density_kg_m3().conj(),
            response.dynamic_bulk_modulus_pa().conj(),
        );
        let eta = rho_i.scale(1.0 / gas.density);
        let (zo, zi) = (ko.scale(radius), ki.scale(radius));
        let (jo, ji, ho) = (regular(zo, 17), regular(zi, 17), outgoing(zo, 17));
        let (mut phase, mut scattered, mut transmitted, mut outward) =
            (C64::ONE, Vec::new(), Vec::new(), 0.0);
        for n in 0..16 {
            let (djo, dji, dho) = (
                derivative(&jo, zo, n),
                derivative(&ji, zi, n),
                derivative(&ho, zo, n),
            );
            let denominator = ko * ji[n] * dho - ki * dji * ho[n] / eta;
            let a = (ki * dji * jo[n] / eta - ko * ji[n] * djo) / denominator;
            let b = (C64::new(0.0, 1.0) / (ko.scale(radius * radius))) / denominator;
            let amplitude = phase.scale((2 * n + 1) as f64);
            scattered.push(amplitude * a);
            transmitted.push(amplitude * b);
            let p = amplitude * b * ji[n];
            let v = amplitude * b * ki * dji / (C64::new(0.0, omega) * rho_i);
            outward += (p * v.conj()).re / (2 * n + 1) as f64;
            phase = phase * C64::new(0.0, 1.0);
        }
        Self {
            radius,
            ko,
            ki,
            rho_o: gas.density,
            rho_i,
            bulk,
            omega,
            scattered,
            transmitted,
            inward_power: -2.0 * std::f64::consts::PI * radius * radius * outward,
        }
    }
    fn axis(&self, r: f64, interior: bool) -> IncidentSample {
        let (k, rho, c) = if interior {
            (self.ki, self.rho_i, &self.transmitted)
        } else {
            (self.ko, C64::from_re(self.rho_o), &self.scattered)
        };
        let z = k.scale(r);
        let f = if interior {
            regular(z, 17)
        } else {
            outgoing(z, 17)
        };
        let (mut p, mut gradient) = (C64::ZERO, C64::ZERO);
        for (n, &coefficient) in c.iter().enumerate() {
            p = p + coefficient * f[n];
            gradient = gradient + coefficient * k * derivative(&f, z, n);
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
                gradient / (C64::new(0.0, self.omega) * rho),
            ],
        }
    }
    fn volume_losses(&self, count: usize) -> [f64; 2] {
        let mut result = [0.0; 2];
        for (x, w) in gauss(count) {
            let r = self.radius * x;
            let z = self.ki.scale(r);
            let j = regular(z, 17);
            let (mut p2, mut gradient2) = (0.0, 0.0);
            for (n, &b) in self.transmitted.iter().enumerate() {
                let weight = b.norm_sq() / (2 * n + 1) as f64;
                p2 += weight * j[n].norm_sq();
                gradient2 += weight
                    * ((self.ki * derivative(&j, z, n)).norm_sq()
                        + (n * (n + 1)) as f64 * j[n].norm_sq() / (r * r));
            }
            let volume = 4.0 * std::f64::consts::PI * self.radius * w * r * r;
            result[0] += volume * 0.5 * self.omega * self.rho_i.im * gradient2
                / self.rho_i.scale(self.omega).norm_sq();
            result[1] += volume * 0.5 * self.omega * self.bulk.recip().im * p2;
        }
        result
    }
}

#[test]
fn finite_porous_body_converges_to_partial_waves_and_positive_volume_heat() {
    let gate = CancelGate::new();
    let gas = gas();
    let law = law();
    let omega = core::f64::consts::TAU * 500.0;
    let radius = 0.02;
    let exact = Sphere::new(&gas, law.effective_fluid(&gas, omega).unwrap(), radius);
    let loss = exact.volume_losses(48);
    let refined = exact.volume_losses(64);
    assert!(loss.iter().all(|&value| value > 0.0));
    for i in 0..2 {
        assert!((loss[i] - refined[i]).abs() < 2e-12 * refined[i]);
    }
    assert!((loss[0] + loss[1] - exact.inward_power).abs() < 2e-11 * exact.inward_power);
    eprintln!(
        "porous sphere continuum: viscous={} W, compression={} W, inward={} W",
        loss[0], loss[1], exact.inward_power
    );
    let (mut previous_field, mut previous_power) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let owner = PorousFluidInterface::new(
            &surface,
            &gas,
            &law,
            PoreFluid::Fluid(&gas),
            omega,
            Options::default(),
            &gate,
        )
        .unwrap();
        assert!(std::ptr::eq(owner.porous_source(), &law));
        assert_eq!(
            *owner.material_response(),
            law.effective_fluid(&gas, omega).unwrap()
        );
        let field = owner
            .scatter_field(plane(C64::ONE), 1e-4, 10000, &gate)
            .unwrap();
        let inside = [0.0, 0.0, 0.37 * radius];
        let outside = [0.0, 0.0, 2.0 * radius];
        let interior = field
            .interior_absorption_adaptive_at(&[inside], 1e-4, 10000, controls(), &gate)
            .unwrap();
        let exterior = field
            .total_pressure_velocity_adaptive_at(&[outside], 1e-4, 10000, controls(), &gate)
            .unwrap();
        let mut error = 0.0_f64;
        for (actual, r, interior_flag) in [
            (&interior.fields, inside[2], true),
            (&exterior, outside[2], false),
        ] {
            let expected = exact.axis(r, interior_flag);
            error = error.max(
                (actual.scalar.pressure[0] - expected.pressure_pa).abs()
                    / expected.pressure_pa.abs(),
            );
            error = error.max(
                (actual.particle_velocity_m_s[0][2] - expected.velocity_m_s[2]).abs()
                    / expected.velocity_m_s[2].abs(),
            );
        }
        let expected = exact.axis(inside[2], true);
        let viscous = 0.5 * omega * exact.rho_i.im * expected.velocity_m_s[2].norm_sq();
        let thermal = 0.5 * omega * exact.bulk.recip().im * expected.pressure_pa.norm_sq();
        assert!(interior.viscous_loss_w_m3[0] > 0.0 && interior.compression_loss_w_m3[0] > 0.0);
        let power = -field.interface_work.surface_power.re;
        let power_error = (power - exact.inward_power).abs() / exact.inward_power;
        eprintln!(
            "porous sphere level={level}: field_error={error:.6e}, boundary_power_error={power_error:.6e}"
        );
        assert!(
            error < previous_field,
            "level={level}: field={error}, previous={previous_field}"
        );
        assert!(
            power_error < previous_power,
            "level={level}: power={power_error}, previous={previous_power}"
        );
        (previous_field, previous_power) = (error, power_error);
        if level == 2 {
            assert!(error < 0.08, "field={error}");
            assert!(power_error < 0.15, "power={power_error}");
            assert!((interior.viscous_loss_w_m3[0] - viscous).abs() < 0.16 * viscous);
            assert!((interior.compression_loss_w_m3[0] - thermal).abs() < 0.16 * thermal);
        }
    }
}

#[test]
fn porous_incident_superposition_keeps_coherent_heat_and_phase() {
    let gate = CancelGate::new();
    let (gas, law) = (gas(), law());
    let surface = SpherePanels::icosphere(0.02, 0).unwrap();
    let owner = PorousFluidInterface::new(
        &surface,
        &gas,
        &law,
        PoreFluid::Fluid(&gas),
        core::f64::consts::TAU * 500.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::ONE,
        }),
        IncidentSource::PlaneWave(PlaneWave {
            direction: [1.0, 0.0, 0.0],
            amplitude_pa: C64::new(0.4, 0.6),
        }),
    ];
    let a = owner
        .scatter_field(IncidentField::Sources(&sources[..1]), 1e-4, 1000, &gate)
        .unwrap();
    let b = owner
        .scatter_field(IncidentField::Sources(&sources[1..]), 1e-4, 1000, &gate)
        .unwrap();
    let combined = owner
        .scatter_field(IncidentField::Sources(&sources), 1e-4, 1000, &gate)
        .unwrap();
    let point = [[0.002, 0.001, 0.004]];
    let fa = a
        .interior_absorption_adaptive_at(&point, 1e-4, 1000, controls(), &gate)
        .unwrap();
    let fb = b
        .interior_absorption_adaptive_at(&point, 1e-4, 1000, controls(), &gate)
        .unwrap();
    let fc = combined
        .interior_absorption_adaptive_at(&point, 1e-4, 1000, controls(), &gate)
        .unwrap();
    close(
        fc.fields.scalar.pressure[0],
        fa.fields.scalar.pressure[0] + fb.fields.scalar.pressure[0],
        3e-13,
    );
    for axis in 0..3 {
        close(
            fc.fields.particle_velocity_m_s[0][axis],
            fa.fields.particle_velocity_m_s[0][axis] + fb.fields.particle_velocity_m_s[0][axis],
            3e-12,
        );
    }
    let heat = fc.viscous_loss_w_m3[0] + fc.compression_loss_w_m3[0];
    let incoherent = fa.viscous_loss_w_m3[0]
        + fa.compression_loss_w_m3[0]
        + fb.viscous_loss_w_m3[0]
        + fb.compression_loss_w_m3[0];
    assert!((heat - incoherent).abs() > 0.01 * heat);
    let rotated_sources = sources.map(|source| match source {
        IncidentSource::PlaneWave(p) => IncidentSource::PlaneWave(PlaneWave {
            direction: p.direction,
            amplitude_pa: p.amplitude_pa * C64::new(0.0, 1.0),
        }),
        _ => unreachable!(),
    });
    let rotated = owner
        .scatter_field(IncidentField::Sources(&rotated_sources), 1e-4, 1000, &gate)
        .unwrap();
    let fr = rotated
        .interior_absorption_adaptive_at(&point, 1e-4, 1000, controls(), &gate)
        .unwrap();
    assert!(
        (fr.viscous_loss_w_m3[0] - fc.viscous_loss_w_m3[0]).abs() < 3e-13 * fc.viscous_loss_w_m3[0]
    );
    assert!(
        (fr.compression_loss_w_m3[0] - fc.compression_loss_w_m3[0]).abs()
            < 3e-13 * fc.compression_loss_w_m3[0]
    );
}

fn relaxing_gas(with_reservoir: bool) -> RelaxationState {
    let text = if with_reservoir {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas(), vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}
#[test]
fn original_relaxing_pore_law_and_refusals_survive_finite_body_binding() {
    let gate = CancelGate::new();
    let (gas, law) = (gas(), law());
    let surface = SpherePanels::icosphere(0.02, 0).unwrap();
    let omega = core::f64::consts::TAU * 500.0;
    let empty = relaxing_gas(false);
    let relaxing = relaxing_gas(true);
    let plain = PorousFluidInterface::new(
        &surface,
        &gas,
        &law,
        PoreFluid::Fluid(&gas),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let identity = PorousFluidInterface::new(
        &surface,
        &gas,
        &law,
        PoreFluid::RelaxingGas(&empty),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(plain.contrast, identity.contrast);
    let owner = PorousFluidInterface::new(
        &surface,
        &gas,
        &law,
        PoreFluid::RelaxingGas(&relaxing),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        *owner.material_response(),
        law.effective_relaxing(&relaxing, omega).unwrap()
    );
    assert_ne!(owner.bulk, plain.bulk);
    match owner.pore_fluid_source() {
        PoreFluid::RelaxingGas(original) => assert!(std::ptr::eq(original, &relaxing)),
        _ => panic!("original pore source lost"),
    }
    assert!(
        PorousFluidInterface::new(
            &surface,
            &gas,
            &law,
            PoreFluid::Fluid(&gas),
            1.0,
            Options::default(),
            &gate,
        )
        .is_err()
    );
    assert!(
        PorousFluidInterface::new(
            &surface,
            &gas,
            &law,
            PoreFluid::Fluid(&gas),
            omega,
            Options {
                maximum_bytes: 0,
                ..Options::default()
            },
            &gate,
        )
        .is_err()
    );
    let internal = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&internal), 1e-4, 1000, &gate)
            .is_err()
    );
    let field = owner
        .scatter_field(plane(C64::ONE), 1e-4, 1000, &gate)
        .unwrap();
    assert!(
        field
            .interior_pressure_velocity_adaptive_at(
                &[[0.0, 0.0, 0.05]],
                1e-4,
                1000,
                controls(),
                &gate,
            )
            .is_err()
    );
    assert!(
        field
            .total_pressure_velocity_adaptive_at(&[[0.0; 3]], 1e-4, 1000, controls(), &gate,)
            .is_err()
    );
    assert!(
        field
            .interior_pressure_velocity_adaptive_at(
                &[[0.0; 3]],
                1e-4,
                1000,
                MolecularNearFieldOptions {
                    maximum_row_bytes: 0,
                    ..controls()
                },
                &gate,
            )
            .is_err()
    );
    gate.request();
    assert!(matches!(
        field.interior_absorption_adaptive_at(&[[0.0; 3]], 1e-4, 1000, controls(), &gate,),
        Err(HelmholtzError::Cancelled)
    ));
}
