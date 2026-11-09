use super::*;
use crate::helmholtz::bc::{
    ForceDipole, IncidentSource, Monopole, PlaneWave, StressQuadrupole,
};
use crate::helmholtz::transmission::FluidInterface;
use super::super::tests::state;
use std::f64::consts::PI;

const OMEGA: f64 = 200.0;

fn controls() -> InternalSourceControls {
    InternalSourceControls {
        minimum_source_clearance_m: 0.01,
        maximum_relative_residual: 1e-10,
    }
}
fn receivers() -> MolecularNearFieldOptions {
    MolecularNearFieldOptions {
        quadrature: near_field::Options {
            relative_tolerance: 1e-6,
            maximum_depth: 12,
            maximum_kernel_evaluations: 4_000_000,
        },
        ..MolecularNearFieldOptions::default()
    }
}
fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
    for (index, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        assert!((a - b).abs() <= tolerance * scale, "index {index}: {a:?} != {b:?}");
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-30)
}
fn sources() -> [IncidentSource; 3] {
    [
        IncidentSource::Monopole(Monopole {
            position_m: [0.1, -0.12, 0.05],
            volume_velocity_m3_s: C64::new(0.003, -0.001),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [-0.15, 0.1, 0.07],
            force_on_fluid_n: [
                C64::new(0.001, 0.002),
                C64::new(-0.0007, 0.001),
                C64::new(0.0003, -0.0002),
            ],
        }),
        IncidentSource::StressQuadrupole(StressQuadrupole {
            position_m: [0.03, 0.02, -0.1],
            stress_moment_n_m: [
                C64::new(0.0003, -0.0002),
                C64::new(-0.0002, 0.0001),
                C64::new(0.0001, 0.0004),
                C64::new(0.00005, 0.0001),
                C64::new(-0.00003, 0.00007),
                C64::new(0.00004, -0.00009),
            ],
        }),
    ]
}

#[test]
fn empty_reservoirs_recover_original_internal_source_fields_and_work() {
    let outside = state(293.15, 101325.0, None);
    let inside = state(353.15, 101325.0, None);
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let gate = CancelGate::new();
    let native = FluidInterface::new(
        &surface, outside.gas(), inside.gas(), OMEGA, Options::default(), &gate,
    ).unwrap();
    let molecular = MolecularFluidInterface::new(
        &surface, &outside, &inside, OMEGA, Options::default(), &gate,
    ).unwrap();
    let terms = sources();
    let incident = IncidentField::Sources(&terms);
    let a = native.radiate_internal_field(incident, controls(), &gate).unwrap();
    let b = molecular.radiate_internal_field(incident, controls(), 1000, &gate).unwrap();
    assert!(std::ptr::eq(molecular.exterior_source(), &outside));
    assert!(std::ptr::eq(molecular.interior_source(), &inside));
    close(b.total_pressure_pa(), &a.transmitted().pressure, 3e-12);
    close(b.total_normal_velocity_m_s(), &a.transmitted().velocity, 3e-11);
    close(b.reflected_pressure_pa(), a.reflected_pressure_pa(), 3e-11);
    close(b.reflected_normal_velocity_m_s(), a.reflected_normal_velocity_m_s(), 3e-11);
    close(
        &[b.interface_work.surface_power],
        &[a.transmitted().power_diagnostics.surface_power],
        3e-11,
    );
    for (region, point) in [
        (InternalSourceRegion::Interior, [0.37, -0.25, 0.2]),
        (InternalSourceRegion::Exterior, [2.1, 0.7, -0.4]),
    ] {
        let va = a.pressure_velocity_at(region, &[point], 0.01, receivers().quadrature, &gate).unwrap();
        let vb = b.pressure_velocity_at(region, &[point], 0.01, receivers(), &gate).unwrap();
        close(&vb.scalar.pressure, &va.scalar.pressure, 2e-9);
        close(&vb.particle_velocity_m_s[0], &va.particle_velocity_m_s[0], 2e-9);
    }
}

fn sine(z: C64) -> C64 {
    C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh())
}
fn cosine(z: C64) -> C64 {
    C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh())
}
fn phase(k: C64, radius: f64) -> C64 {
    C64::new((k.re * radius).cos(), (k.re * radius).sin()).scale((-k.im * radius).exp())
}

// Independent spherical l=0 continuity for radius one. p_i=C exp(i ki r)/r
// + A j0(ki r); p_o=B exp(i ko r)/r. Real rho supplies Euler momentum.
fn central_coefficients(ko: C64, ki: C64, rho_o: f64, rho_i: f64, q: C64) -> [C64; 3] {
    let imaginary = C64::new(0.0, 1.0);
    let c = C64::new(0.0, -OMEGA * rho_i / (4.0 * PI)) * q;
    let p = c * phase(ki, 1.0);
    let dp = p * (imaginary * ki - C64::ONE);
    let j = sine(ki) / ki;
    let dj = cosine(ki) - j;
    let h = phase(ko, 1.0);
    let dh = h * (imaginary * ko - C64::ONE);
    let b = (dp * j - p * dj) / ((dh * j).scale(rho_i / rho_o) - h * dj);
    let a = (b * h - p) / j;
    [c, a, b]
}
fn central_field(k: C64, rho: f64, r: f64, outgoing: C64, regular: C64) -> (C64, C64) {
    let z = k.scale(r);
    let j = sine(z) / z;
    let wave = (outgoing * phase(k, r)).scale(1.0 / r);
    let p = wave + regular * j;
    let dp = wave * (C64::new(0.0, 1.0) * k - C64::from_re(1.0 / r))
        + (regular * (cosine(z) - j)).scale(1.0 / r);
    (p, dp / C64::new(0.0, OMEGA * rho))
}
fn outward_power(k: C64, rho: f64, r: f64, outgoing: C64, regular: C64) -> f64 {
    let (p, v) = central_field(k, rho, r, outgoing, regular);
    2.0 * PI * r * r * (p * v.conj()).re
}

#[test]
fn internal_molecular_emitter_refines_to_complex_spherical_transmission_and_outward_work() {
    // These are manufactured caloric reservoirs, not measured air absorption.
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(353.15, 101325.0, Some(5.0));
    let (ko, ki) = (
        outside.relaxing_wavenumber(OMEGA).unwrap(),
        inside.relaxing_wavenumber(OMEGA).unwrap(),
    );
    let (rho_o, rho_i) = (outside.gas().density, inside.gas().density);
    let q = C64::new(0.004, -0.001);
    let terms = [Monopole { position_m: [0.0; 3], volume_velocity_m3_s: q }];
    let [c, a, b] = central_coefficients(ko, ki, rho_o, rho_i, q);
    let (po, vo) = central_field(ko, rho_o, 1.0, b, C64::ZERO);
    let (pi, vi) = central_field(ki, rho_i, 1.0, c, a);
    assert!(relative(pi, po) < 2e-13 && relative(vi, vo) < 2e-13);
    let exact_work = outward_power(ko, rho_o, 1.0, b, C64::ZERO);
    assert!(exact_work > 0.0);
    let farther_work = outward_power(ko, rho_o, 3.0, b, C64::ZERO);
    assert!(farther_work > 0.0 && farther_work < exact_work);
    let compliance = outside.caloric_response(OMEGA).unwrap().bulk_modulus_pa.recip();
    let volume_loss = 2.0 * PI * OMEGA * compliance.im * b.norm_sq()
        * ((-2.0 * ko.im).exp() - (-6.0 * ko.im).exp()) / (2.0 * ko.im);
    assert!(((exact_work - farther_work) / volume_loss - 1.0).abs() < 2e-10);
    let gate = CancelGate::new();
    let mut previous_error = f64::INFINITY;
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let owner = MolecularFluidInterface::new(
            &surface, &outside, &inside, OMEGA, Options::default(), &gate,
        ).unwrap();
        let field = owner.radiate_internal_field(
            IncidentField::Monopoles(&terms), controls(), 1000, &gate,
        ).unwrap();
        let mut error = 0.0;
        for (region, radius, k, rho, outgoing, regular) in [
            (InternalSourceRegion::Interior, 0.35, ki, rho_i, c, a),
            (InternalSourceRegion::Exterior, 3.0, ko, rho_o, b, C64::ZERO),
        ] {
            let actual = field.pressure_velocity_at(
                region, &[[0.0, 0.0, radius]], 0.01, receivers(), &gate,
            ).unwrap();
            let (p, v) = central_field(k, rho, radius, outgoing, regular);
            error += relative(actual.scalar.pressure[0], p);
            error += ((actual.particle_velocity_m_s[0][2] - v).norm_sq()
                + actual.particle_velocity_m_s[0][0].norm_sq()
                + actual.particle_velocity_m_s[0][1].norm_sq()).sqrt() / v.abs();
        }
        let reflected = field.reflected_pressure_velocity_at(
            &[[0.0, 0.0, 0.35]], 0.01, receivers(), &gate,
        ).unwrap();
        let (pr, _) = central_field(ki, rho_i, 0.35, C64::ZERO, a);
        error += relative(reflected.scalar.pressure[0], pr);
        assert!(error < previous_error, "level {level}: {error} >= {previous_error}");
        previous_error = error;
        let work = field.interface_work.surface_power.re;
        assert!(work > 0.0);
        assert!(field.relative_backward_residual < controls().maximum_relative_residual);
        if level == 2 {
            assert!(error < 0.15, "full pressure/velocity/reflection error {error}");
            assert!((work / exact_work - 1.0).abs() < 0.06);
        }
    }
}

#[test]
fn matched_molecular_media_recover_exact_mixed_emitters_and_weak_reflection_survives() {
    let gas = state(293.15, 101325.0, Some(5.0));
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new();
    let owner = MolecularFluidInterface::new(
        &surface, &gas, &gas, OMEGA, Options::default(), &gate,
    ).unwrap();
    let terms = sources();
    let sources = IncidentField::Sources(&terms);
    let field = owner.radiate_internal_field(sources, controls(), 1000, &gate).unwrap();
    assert_eq!(field.relative_backward_residual, 0.0);
    assert!(field.reflected_pressure_pa().iter()
        .chain(field.reflected_normal_velocity_m_s()).all(|&z| z == C64::ZERO));
    for (region, point) in [
        (InternalSourceRegion::Interior, [0.37, -0.25, 0.2]),
        (InternalSourceRegion::Exterior, [2.1, 0.7, -0.4]),
    ] {
        let actual = field.pressure_velocity_at(region, &[point], 0.01, receivers(), &gate).unwrap();
        let expected = sources.sample_wave(
            owner.k_i, gas.gas().density * OMEGA, None, point, &mut || false,
        ).unwrap();
        assert_eq!(actual.scalar.pressure[0], expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s[0], expected.velocity_m_s);
        assert_eq!(actual.scalar.quadrature_error_estimate_pa[0], 0.0);
        assert_eq!(actual.quadrature_error_estimate_m_s[0], [0.0; 3]);
    }
    let emitter = [Monopole {
        position_m: [0.1, -0.05, 0.07],
        volume_velocity_m3_s: C64::new(0.003, -0.001),
    }];
    let mut first = None;
    for jump in [2.0_f64.powi(-18), 2.0_f64.powi(-30)] {
        let inside = state(293.15, 101325.0, Some(5.0 * (1.0 + jump)));
        let contrast = MolecularFluidInterface::new(
            &surface, &gas, &inside, OMEGA, Options::default(), &gate,
        ).unwrap();
        let response = contrast.radiate_internal_field(
            IncidentField::Monopoles(&emitter), controls(), 1000, &gate,
        ).unwrap();
        let reflected = response.reflected_pressure_velocity_at(
            &[[0.3, -0.15, 0.2]], 0.01, receivers(), &gate,
        ).unwrap();
        assert!(reflected.scalar.pressure[0].abs() > 1e-22);
        let normalized = reflected.scalar.pressure[0].scale(1.0 / jump);
        if let Some(reference) = first {
            close(&[normalized], &[reference], 2e-4);
        }
        first = Some(normalized);
    }
}

#[test]
fn molecular_internal_volume_force_and_stress_sources_interfere_coherently() {
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(293.15, 101325.0, Some(5.0));
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new();
    let owner = MolecularFluidInterface::new(
        &surface, &outside, &inside, OMEGA, Options::default(), &gate,
    ).unwrap();
    let terms = sources();
    let point = [[2.1, 0.7, -0.4]];
    let mixed = owner.radiate_internal_field(
        IncidentField::Sources(&terms), controls(), 1000, &gate,
    ).unwrap();
    let actual = mixed.pressure_velocity_at(
        InternalSourceRegion::Exterior, &point, 0.01, receivers(), &gate,
    ).unwrap();
    let mut pressure = C64::ZERO;
    let mut velocity = [C64::ZERO; 3];
    let mut separate_intensity = [0.0; 3];
    for term in &terms {
        let one = owner.radiate_internal_field(
            IncidentField::Sources(std::slice::from_ref(term)), controls(), 1000, &gate,
        ).unwrap();
        let observed = one.pressure_velocity_at(
            InternalSourceRegion::Exterior, &point, 0.01, receivers(), &gate,
        ).unwrap();
        pressure = pressure + observed.scalar.pressure[0];
        for axis in 0..3 {
            velocity[axis] = velocity[axis] + observed.particle_velocity_m_s[0][axis];
            separate_intensity[axis] += 0.5
                * (observed.scalar.pressure[0] * observed.particle_velocity_m_s[0][axis].conj()).re;
        }
    }
    close(&actual.scalar.pressure, &[pressure], 1e-10);
    close(&actual.particle_velocity_m_s[0], &velocity, 1e-10);
    let coherent = 0.5 * (pressure * velocity[0].conj()).re;
    assert!((coherent - separate_intensity[0]).abs() > 1e-4 * coherent.abs());
}

#[test]
fn internal_molecular_sources_keep_native_geometry_and_explicit_work_refusals() {
    let gas = state(293.15, 101325.0, Some(5.0));
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let gate = CancelGate::new();
    let owner = MolecularFluidInterface::new(
        &surface, &gas, &gas, OMEGA, Options::default(), &gate,
    ).unwrap();
    let emitter = [Monopole { position_m: [0.0; 3], volume_velocity_m3_s: C64::ONE }];
    let incident = IncidentField::Monopoles(&emitter);
    let field = owner.radiate_internal_field(incident, controls(), 20, &gate).unwrap();
    assert!(owner.radiate_internal_field(incident, controls(), 19, &gate).is_err());
    for position_m in [[2.0, 0.0, 0.0], surface.centroids()[0]] {
        let wrong = [Monopole { position_m, ..emitter[0] }];
        assert!(owner.radiate_internal_field(
            IncidentField::Monopoles(&wrong), controls(), 1000, &gate,
        ).is_err());
    }
    let plane = PlaneWave { direction: [0.0, 0.0, 1.0], amplitude_pa: C64::ONE };
    assert!(owner.radiate_internal_field(IncidentField::PlaneWave(plane), controls(), 1000, &gate).is_err());
    let mixed = [IncidentSource::Monopole(emitter[0]), IncidentSource::PlaneWave(plane)];
    assert!(owner.radiate_internal_field(IncidentField::Sources(&mixed), controls(), 1000, &gate).is_err());
    assert!(owner.radiate_internal_field(IncidentField::Monopoles(&[]), controls(), 1000, &gate).is_err());
    for limits in [
        InternalSourceControls { minimum_source_clearance_m: 0.0, ..controls() },
        InternalSourceControls { maximum_relative_residual: f64::NAN, ..controls() },
    ] {
        assert!(owner.radiate_internal_field(incident, limits, 1000, &gate).is_err());
    }
    for (region, point) in [
        (InternalSourceRegion::Interior, [2.0, 0.0, 0.0]),
        (InternalSourceRegion::Exterior, [0.2, 0.0, 0.0]),
        (InternalSourceRegion::Interior, [0.0; 3]),
        (InternalSourceRegion::Exterior, surface.centroids()[0]),
    ] {
        assert!(field.pressure_velocity_at(region, &[point], 0.01, receivers(), &gate).is_err());
    }
    assert!(field.pressure_velocity_at(
        InternalSourceRegion::Exterior, &[[2.0, 0.0, 0.0]; 21], 0.01, receivers(), &gate,
    ).is_err());
    let no_rows = MolecularNearFieldOptions { maximum_row_bytes: 1, ..receivers() };
    assert!(field.pressure_velocity_at(
        InternalSourceRegion::Exterior, &[[2.0, 0.0, 0.0]], 0.01, no_rows, &gate,
    ).is_err());
    let no_work = MolecularNearFieldOptions {
        quadrature: near_field::Options { maximum_kernel_evaluations: 80 * 20 - 1, ..receivers().quadrature },
        ..receivers()
    };
    assert!(field.pressure_velocity_at(
        InternalSourceRegion::Exterior, &[[2.0, 0.0, 0.0]], 0.01, no_work, &gate,
    ).is_err());
    gate.request();
    assert!(matches!(
        owner.radiate_internal_field(incident, controls(), 1000, &gate),
        Err(HelmholtzError::Cancelled),
    ));
    assert!(matches!(
        field.pressure_velocity_at(
            InternalSourceRegion::Exterior, &[[2.0, 0.0, 0.0]], 0.01, receivers(), &gate,
        ),
        Err(HelmholtzError::Cancelled),
    ));
}
