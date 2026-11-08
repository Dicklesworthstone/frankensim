//! Internal physical sources crossing a penetrable fluid boundary.
//! The homogeneous media are explicit numerical conformance fixtures.

use fs_bem::helmholtz::bc::{
    ForceDipole, IncidentField, IncidentSource, Monopole, PlaneWave, StressQuadrupole,
};
use fs_bem::helmholtz::transmission::internal::{InternalSourceControls, InternalSourceRegion};
use fs_bem::helmholtz::transmission::{FluidInterface, Options};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::{near_field, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;
use std::cell::Cell;
use std::f64::consts::PI;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn controls() -> InternalSourceControls {
    InternalSourceControls {
        minimum_source_clearance_m: 0.01,
        maximum_relative_residual: 1e-10,
    }
}
fn receivers() -> near_field::Options {
    near_field::Options {
        relative_tolerance: 1e-6,
        maximum_depth: 12,
        maximum_kernel_evaluations: 4_000_000,
    }
}
fn phase(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-25)
}
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(relative(a, b) < tolerance, "{a:?} != {b:?}");
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
fn matched_media_keep_exact_internal_and_external_mixed_source_fields() {
    struct Counted {
        queries: Cell<usize>,
    }
    impl FluidAcoustics for Counted {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            self.queries.set(self.queries.get() + 1);
            Ok(fluid(2.0, 4.0))
        }
    }
    let material = Counted {
        queries: Cell::new(0),
    };
    let source_view: &dyn FluidAcoustics = &material;
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let gate = CancelGate::new();
    let model = FluidInterface::new(
        &surface,
        source_view,
        source_view,
        1.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(material.queries.get(), 1);
    let terms = sources();
    let field = IncidentField::Sources(&terms);
    for _ in 0..2 {
        let response = model
            .radiate_internal_field(field, controls(), &gate)
            .unwrap();
        assert_eq!(response.relative_backward_residual, 0.0);
        assert!(
            response
                .reflected_pressure_pa()
                .iter()
                .chain(response.reflected_normal_velocity_m_s())
                .all(|&z| z == C64::ZERO)
        );
        for (region, point) in [
            (InternalSourceRegion::Interior, [0.37, -0.25, 0.2]),
            (InternalSourceRegion::Exterior, [2.1, 0.7, -0.4]),
        ] {
            let result = response
                .pressure_velocity_at(region, &[point], 0.01, receivers(), &gate)
                .unwrap();
            let expected = field
                .pressure_velocity_at(
                    0.25,
                    Medium {
                        density: 2.0,
                        sound_speed: 4.0,
                    },
                    point,
                )
                .unwrap();
            assert_eq!(result.scalar.pressure[0], expected.pressure_pa);
            assert_eq!(result.particle_velocity_m_s[0], expected.velocity_m_s);
            assert_eq!(result.scalar.quadrature_error_estimate_pa[0], 0.0);
            assert_eq!(result.quadrature_error_estimate_m_s[0], [0.0; 3]);
        }
    }
    assert_eq!(material.queries.get(), 1);
}

// Independent spherical l=0 continuity solution, with R=1. The interior is
// C exp(i ki r)/r + A j0(ki r), and the exterior is B exp(i ko r)/r.
// Match pressure and p'/(i omega rho) directly; no boundary integral enters.
fn central_coefficients(
    omega: f64,
    rho_o: f64,
    c_o: f64,
    rho_i: f64,
    c_i: f64,
    q: C64,
) -> (C64, C64, C64) {
    let (ko, ki) = (omega / c_o, omega / c_i);
    let c = C64::new(0.0, -omega * rho_i / (4.0 * PI)) * q;
    let p = c * phase(ki);
    let dp = p * C64::new(-1.0, ki);
    let j = ki.sin() / ki;
    let dj = ki.cos() - ki.sin() / ki;
    let h = phase(ko);
    let dh = h * C64::new(-1.0, ko);
    let b = (dp.scale(j) - p.scale(dj)) / (dh.scale(rho_i * j / rho_o) - h.scale(dj));
    let a = (b * h - p).scale(1.0 / j);
    (c, a, b)
}
fn central_field(
    omega: f64,
    rho: f64,
    c: f64,
    radius: f64,
    outgoing: C64,
    regular: C64,
) -> (C64, C64) {
    let k = omega / c;
    let kr = k * radius;
    let wave = outgoing * phase(kr).scale(1.0 / radius);
    let pressure = wave + regular.scale(kr.sin() / kr);
    let derivative = wave * C64::new(-1.0 / radius, k)
        + regular.scale(k * (kr.cos() / kr - kr.sin() / (kr * kr)));
    (pressure, derivative * C64::new(0.0, -1.0 / (omega * rho)))
}

#[test]
fn internal_monopole_refines_to_independent_spherical_transmission_and_power() {
    let (omega, rho_o, c_o, rho_i, c_i) = (0.6, 1.2, 1.0, 2.7, 0.8);
    let exterior = fluid(rho_o, c_o);
    let interior = fluid(rho_i, c_i);
    let q = C64::new(0.004, -0.001);
    let terms = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: q,
    }];
    let (c, a, b) = central_coefficients(omega, rho_o, c_o, rho_i, c_i, q);
    let expected_power = 2.0 * PI * b.norm_sq() / (rho_o * c_o);
    let gate = CancelGate::new();
    let mut previous = f64::INFINITY;
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let model = FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            omega,
            Options::default(),
            &gate,
        )
        .unwrap();
        let response = model
            .radiate_internal_field(IncidentField::Monopoles(&terms), controls(), &gate)
            .unwrap();
        let mut error = 0.0;
        for (region, radius, rho, sound_speed, outgoing, regular) in [
            (InternalSourceRegion::Interior, 0.35, rho_i, c_i, c, a),
            (
                InternalSourceRegion::Exterior,
                3.0,
                rho_o,
                c_o,
                b,
                C64::ZERO,
            ),
        ] {
            let actual = response
                .pressure_velocity_at(region, &[[0.0, 0.0, radius]], 0.01, receivers(), &gate)
                .unwrap();
            let (pressure, velocity) =
                central_field(omega, rho, sound_speed, radius, outgoing, regular);
            error += relative(actual.scalar.pressure[0], pressure);
            error += ((actual.particle_velocity_m_s[0][2] - velocity).norm_sq()
                + actual.particle_velocity_m_s[0][0].norm_sq()
                + actual.particle_velocity_m_s[0][1].norm_sq())
            .sqrt()
                / velocity.abs();
        }
        assert!(error < previous, "level {level}: {error} >= {previous}");
        previous = error;
        let power = response.transmitted().power_diagnostics.surface_power.re;
        assert!(power > 0.0);
        assert!(response.relative_backward_residual < controls().maximum_relative_residual);
        if level == 2 {
            assert!(error < 0.08, "full pressure/velocity error {error}");
            assert!((power / expected_power - 1.0).abs() < 0.06);
        }
    }
}

#[test]
fn weak_interior_reflection_survives_nearly_matched_media() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let exterior = fluid(1.0, 4.0);
    let terms = [Monopole {
        position_m: [0.1, -0.05, 0.07],
        volume_velocity_m3_s: C64::new(0.003, -0.001),
    }];
    let gate = CancelGate::new();
    let mut first = None;
    for jump in [2.0_f64.powi(-20), 2.0_f64.powi(-40)] {
        let interior = fluid(1.0 + jump, 4.0);
        let model = FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            1.0,
            Options::default(),
            &gate,
        )
        .unwrap();
        let response = model
            .radiate_internal_field(IncidentField::Monopoles(&terms), controls(), &gate)
            .unwrap();
        let reflected = response
            .reflected_pressure_velocity_at(&[[0.3, -0.15, 0.2]], 0.01, receivers(), &gate)
            .unwrap();
        assert!(reflected.scalar.pressure[0].abs() > 1e-20);
        let normalized = reflected.scalar.pressure[0].scale(1.0 / jump);
        if let Some(reference) = first {
            near(normalized, reference, 5e-6);
        }
        first = Some(normalized);
    }
}

#[test]
fn all_three_internal_multipoles_superpose_before_exterior_intensity() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let exterior = fluid(1.0, 3.0);
    let interior = fluid(1.7, 2.5);
    let gate = CancelGate::new();
    let model = FluidInterface::new(
        &surface,
        &exterior,
        &interior,
        1.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let terms = sources();
    let point = [[2.1, 0.7, -0.4]];
    let mixed = model
        .radiate_internal_field(IncidentField::Sources(&terms), controls(), &gate)
        .unwrap()
        .pressure_velocity_at(
            InternalSourceRegion::Exterior,
            &point,
            0.01,
            receivers(),
            &gate,
        )
        .unwrap();
    let mut pressure = C64::ZERO;
    let mut velocity = [C64::ZERO; 3];
    let mut separate_intensity = [0.0; 3];
    for term in &terms {
        let result = model
            .radiate_internal_field(
                IncidentField::Sources(std::slice::from_ref(term)),
                controls(),
                &gate,
            )
            .unwrap()
            .pressure_velocity_at(
                InternalSourceRegion::Exterior,
                &point,
                0.01,
                receivers(),
                &gate,
            )
            .unwrap();
        pressure = pressure + result.scalar.pressure[0];
        for axis in 0..3 {
            velocity[axis] = velocity[axis] + result.particle_velocity_m_s[0][axis];
            separate_intensity[axis] +=
                0.5 * (result.scalar.pressure[0] * result.particle_velocity_m_s[0][axis].conj()).re;
        }
    }
    near(mixed.scalar.pressure[0], pressure, 1e-10);
    for axis in 0..3 {
        near(mixed.particle_velocity_m_s[0][axis], velocity[axis], 1e-10);
    }
    let intensity = 0.5 * (pressure * velocity[0].conj()).re;
    assert!((intensity - separate_intensity[0]).abs() > 1e-4 * intensity.abs());
}

#[test]
fn internal_source_and_receiver_regions_clearance_and_cancellation_refuse() {
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let material = fluid(1.0, 4.0);
    let gate = CancelGate::new();
    let model = FluidInterface::new(
        &surface,
        &material,
        &material,
        1.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let terms = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    let field = IncidentField::Monopoles(&terms);
    let response = model
        .radiate_internal_field(field, controls(), &gate)
        .unwrap();
    for (region, point) in [
        (InternalSourceRegion::Interior, [2.0, 0.0, 0.0]),
        (InternalSourceRegion::Exterior, [0.2, 0.0, 0.0]),
        (InternalSourceRegion::Interior, [0.0; 3]),
    ] {
        assert!(
            response
                .pressure_velocity_at(region, &[point], 0.01, receivers(), &gate)
                .is_err()
        );
    }
    for position_m in [[2.0, 0.0, 0.0], surface.centroids()[0]] {
        let wrong = [Monopole {
            position_m,
            ..terms[0]
        }];
        assert!(
            model
                .radiate_internal_field(IncidentField::Monopoles(&wrong), controls(), &gate)
                .is_err()
        );
    }
    let plane = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    });
    assert!(
        model
            .radiate_internal_field(plane, controls(), &gate)
            .is_err()
    );
    for limits in [
        InternalSourceControls {
            minimum_source_clearance_m: 0.0,
            ..controls()
        },
        InternalSourceControls {
            maximum_relative_residual: f64::NAN,
            ..controls()
        },
    ] {
        assert!(model.radiate_internal_field(field, limits, &gate).is_err());
    }
    assert!(
        model
            .radiate_internal_field(IncidentField::Monopoles(&[]), controls(), &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        model.radiate_internal_field(field, controls(), &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        response.pressure_velocity_at(
            InternalSourceRegion::Exterior,
            &[[2.0, 0.0, 0.0]],
            0.01,
            receivers(),
            &gate,
        ),
        Err(HelmholtzError::Cancelled)
    ));
}
