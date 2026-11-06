//! Three adjacent fluids, two shared curved interfaces, one coupled response.
//! The concentric reference solves independent spherical-wave continuity;
//! offset-core and transparent-interface cases exercise the general owner.
use fs_bem::helmholtz::bc::{IncidentField, IncidentSource, Monopole, PlaneWave};
use fs_bem::helmholtz::transmission::nested::{
    CoreShellFluidInterface, CoreShellInterface, CoreShellRegion,
};
use fs_bem::helmholtz::transmission::{FluidInterface, Options, TransmissionError};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::{near_field, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;

#[path = "support/fluid_core_shell_oracle.rs"]
mod oracle;
use oracle::CoreShellOracle;

const OMEGA: f64 = 0.6;
const DIRECTIONS: [[f64; 3]; 3] = [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn fluids() -> [FluidAcousticState; 3] {
    [fluid(1.0, 1.0), fluid(1.6, 0.75), fluid(0.65, 1.3)]
}
fn views(states: &[FluidAcousticState; 3]) -> [&dyn FluidAcoustics; 3] {
    [&states[0], &states[1], &states[2]]
}
fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    })
}
fn sphere(radius: f64, level: u32, center: [f64; 3]) -> SpherePanels {
    let original = SpherePanels::icosphere(radius, level).unwrap();
    SpherePanels::from_triangles(
        original
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| {
                triangle.map(|point| std::array::from_fn(|axis| point[axis] + center[axis]))
            })
            .collect(),
    )
    .unwrap()
}
fn exp(phase: f64) -> C64 {
    C64::new(phase.cos(), phase.sin())
}
fn norm(values: &[C64]) -> f64 {
    values.iter().map(|z| z.abs().powi(2)).sum::<f64>().sqrt()
}
fn relative_error(actual: &[C64], expected: &[C64]) -> f64 {
    let difference: Vec<_> = actual.iter().zip(expected).map(|(&a, &b)| a - b).collect();
    norm(&difference) / norm(expected)
}
fn assert_close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(1.0, f64::max);
    for (index, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - b).abs() < tolerance * scale,
            "index {index}: {a:?} != {b:?}"
        );
    }
}

#[test]
fn concentric_refraction_converges_to_independent_three_fluid_waves() {
    let reference =
        CoreShellOracle::new(OMEGA, [1.0, 0.45], [1.0, 1.6, 0.65], [1.0, 0.75, 1.3], 18);
    assert!(reference.maximum_continuity_defect() < 2e-14);
    assert!(reference.maximum_lossless_defect() < 2e-14);
    let exact_far: Vec<_> = DIRECTIONS.iter().map(|d| reference.far(d[2])).collect();
    let receivers = [
        (CoreShellRegion::Exterior, 0, 1.8),
        (CoreShellRegion::Shell, 1, 0.65),
        (CoreShellRegion::Core, 2, 0.15),
    ];
    let mut errors = Vec::new();
    let gate = CancelGate::new();
    let states = fluids();
    for level in [1, 2] {
        let outer = sphere(1.0, level, [0.0; 3]);
        let core = sphere(0.45, level, [0.0; 3]);
        let owner = CoreShellFluidInterface::new(
            &outer,
            &core,
            views(&states),
            OMEGA,
            0.01,
            Options::default(),
            &gate,
        )
        .unwrap();
        assert_eq!(owner.active_interface_count(), 2);
        let field = owner.scatter_field(plane(), 0.01, &gate).unwrap();
        let actual_far = field.far_field(&DIRECTIONS, 1_000_000, &gate).unwrap();
        let mut mesh_errors = vec![relative_error(&actual_far, &exact_far)];
        for &(region, index, r) in &receivers {
            let sample = field
                .pressure_velocity_at(
                    region,
                    &[[0.0, 0.0, r]],
                    0.01,
                    near_field::Options::default(),
                    &gate,
                )
                .unwrap();
            let (pressure, velocity) = reference.axial(index, r);
            mesh_errors.push(relative_error(&sample.scalar.pressure, &[pressure]));
            mesh_errors.push(relative_error(
                &sample.particle_velocity_m_s[0],
                &[C64::ZERO, C64::ZERO, velocity],
            ));
            assert!(sample.scalar.quadrature_error_estimate_pa[0].is_finite());
            assert!(
                sample.quadrature_error_estimate_m_s[0]
                    .iter()
                    .all(|v| v.is_finite() && *v >= 0.0)
            );
        }
        for boundary in [CoreShellInterface::Outer, CoreShellInterface::Core] {
            let trace = field.interface_trace(boundary).unwrap();
            assert_eq!(
                trace.pressure_pa().len(),
                owner.surface(boundary).areas().len()
            );
            assert!(trace.interface_work().surface_power.abs().is_finite());
        }
        errors.push(mesh_errors);
    }
    for (quantity, (&coarse, &fine)) in errors[0].iter().zip(&errors[1]).enumerate() {
        assert!(
            fine < coarse,
            "quantity {quantity}: refinement {coarse} -> {fine}"
        );
        // Constant panels and the original equivalent-disc self terms have
        // mesh error. This gate concerns that error, not quadrature estimates.
        assert!(
            fine < 0.08,
            "quantity {quantity}: fine relative error {fine}"
        );
    }
}

#[test]
fn matched_adjacent_fluids_reduce_to_the_original_single_interface() {
    let outside = fluid(1.0, 1.0);
    let inside = fluid(1.6, 0.75);
    let outer = sphere(1.0, 1, [0.0; 3]);
    let core = sphere(0.35, 1, [0.12, 0.0, 0.08]);
    let gate = CancelGate::new();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.8, -0.3),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.2, -0.1, 2.0],
            volume_velocity_m3_s: C64::new(0.4, 0.7),
        }),
    ];
    let incident = IncidentField::Sources(&sources);
    for transparent_outer in [true, false] {
        let materials: [&dyn FluidAcoustics; 3] = if transparent_outer {
            [&outside, &outside, &inside]
        } else {
            [&outside, &inside, &inside]
        };
        let (surface, active, eliminated) = if transparent_outer {
            (&core, CoreShellInterface::Core, CoreShellInterface::Outer)
        } else {
            (&outer, CoreShellInterface::Outer, CoreShellInterface::Core)
        };
        let nested = CoreShellFluidInterface::new(
            &outer,
            &core,
            materials,
            OMEGA,
            0.01,
            Options::default(),
            &gate,
        )
        .unwrap();
        let original =
            FluidInterface::new(surface, &outside, &inside, OMEGA, Options::default(), &gate)
                .unwrap();
        assert_eq!(nested.active_interface_count(), 1);
        let actual = nested.scatter_field(incident, 0.01, &gate).unwrap();
        let expected = original.scatter_field(incident, 0.01, &gate).unwrap();
        assert!(std::ptr::eq(actual.radiating_surface(), surface));
        assert!(actual.interface_trace(eliminated).is_none());
        assert_close(
            &actual.scattered().pressure,
            &expected.solution().scattered().pressure,
            3e-12,
        );
        assert_close(
            &actual.scattered().velocity,
            &expected.solution().scattered().velocity,
            3e-12,
        );
        let trace = actual.interface_trace(active).unwrap();
        assert_close(
            trace.pressure_pa(),
            expected.solution().total_pressure_pa(),
            3e-12,
        );
        assert_close(
            trace.normal_velocity_m_s(),
            expected.solution().total_normal_velocity_m_s(),
            3e-12,
        );
        // The core sample is inside either surviving inclusion. Eliminating
        // a matched boundary must retain the original interior propagation.
        let point = [[0.12, 0.0, 0.08]];
        let actual_inside = actual
            .pressure_velocity_at(
                CoreShellRegion::Core,
                &point,
                0.01,
                near_field::Options::default(),
                &gate,
            )
            .unwrap();
        let expected_inside = expected
            .solution()
            .interior_pressure_at(&point, 0.01, near_field::Options::default())
            .unwrap();
        assert_close(
            &actual_inside.scalar.pressure,
            &expected_inside.pressure,
            2e-7,
        );
    }
}

#[test]
fn completely_matched_nested_geometry_retains_exact_incident_fields() {
    let outer = sphere(1.0, 0, [0.0; 3]);
    let core = sphere(0.35, 0, [0.12, 0.0, 0.08]);
    let material = fluid(1.0, 1.0);
    let view: &dyn FluidAcoustics = &material;
    let gate = CancelGate::new();
    let owner = CoreShellFluidInterface::new(
        &outer,
        &core,
        [view; 3],
        OMEGA,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(owner.active_interface_count(), 0);
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.8, -0.3),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.2, -0.1, 2.0],
            volume_velocity_m3_s: C64::new(0.4, 0.7),
        }),
    ];
    let incident = IncidentField::Sources(&sources);
    let field = owner.scatter_field(incident, 0.01, &gate).unwrap();
    for value in field
        .scattered()
        .pressure
        .iter()
        .chain(&field.scattered().velocity)
    {
        assert_eq!(*value, C64::ZERO);
    }
    for boundary in [CoreShellInterface::Outer, CoreShellInterface::Core] {
        assert!(field.interface_trace(boundary).is_none());
    }
    for (region, point) in [
        (CoreShellRegion::Exterior, [0.0, 0.0, 1.6]),
        (CoreShellRegion::Shell, [0.0, 0.0, 0.65]),
        (CoreShellRegion::Core, [0.12, 0.0, 0.08]),
    ] {
        let expected = incident
            .pressure_velocity_at(
                OMEGA,
                Medium {
                    density: 1.0,
                    sound_speed: 1.0,
                },
                point,
            )
            .unwrap();
        let actual = field
            .pressure_velocity_at(
                region,
                &[point],
                0.01,
                near_field::Options::default(),
                &gate,
            )
            .unwrap();
        assert_eq!(actual.scalar.pressure, vec![expected.pressure_pa]);
        assert_eq!(actual.particle_velocity_m_s, vec![expected.velocity_m_s]);
        assert_eq!(actual.scalar.quadrature_error_estimate_pa, vec![0.0]);
        assert_eq!(actual.quadrature_error_estimate_m_s, vec![[0.0; 3]]);
    }
}

#[test]
fn offset_core_retains_translation_phase_in_all_three_regions() {
    let states = fluids();
    let shift = [2.0, -3.0, 4.0];
    let core_center = [0.12, 0.0, 0.08];
    let translated_core = std::array::from_fn(|axis| core_center[axis] + shift[axis]);
    let outer = sphere(1.0, 1, [0.0; 3]);
    let core = sphere(0.35, 1, core_center);
    let moved_outer = sphere(1.0, 1, shift);
    let moved_core = sphere(0.35, 1, translated_core);
    let gate = CancelGate::new();
    let original = CoreShellFluidInterface::new(
        &outer,
        &core,
        views(&states),
        OMEGA,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let translated = CoreShellFluidInterface::new(
        &moved_outer,
        &moved_core,
        views(&states),
        OMEGA,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let a = original.scatter_field(plane(), 0.01, &gate).unwrap();
    let b = translated.scatter_field(plane(), 0.01, &gate).unwrap();
    let far_a = a.far_field(&DIRECTIONS, 1_000_000, &gate).unwrap();
    let far_b = b.far_field(&DIRECTIONS, 1_000_000, &gate).unwrap();
    for (i, d) in DIRECTIONS.iter().enumerate() {
        let outgoing_phase =
            OMEGA * (shift[2] - (0..3).map(|axis| d[axis] * shift[axis]).sum::<f64>());
        assert_close(&[far_b[i]], &[far_a[i] * exp(outgoing_phase)], 3e-11);
    }
    for (region, point) in [
        (CoreShellRegion::Exterior, [0.0, 0.0, 1.8]),
        (CoreShellRegion::Shell, [0.0, 0.0, 0.65]),
        (CoreShellRegion::Core, core_center),
    ] {
        let shifted = std::array::from_fn(|axis| point[axis] + shift[axis]);
        let observed_a = a
            .pressure_velocity_at(
                region,
                &[point],
                0.01,
                near_field::Options::default(),
                &gate,
            )
            .unwrap();
        let observed_b = b
            .pressure_velocity_at(
                region,
                &[shifted],
                0.01,
                near_field::Options::default(),
                &gate,
            )
            .unwrap();
        let phase = exp(OMEGA * shift[2]);
        assert_close(
            &observed_b.scalar.pressure,
            &[observed_a.scalar.pressure[0] * phase],
            2e-8,
        );
        assert_close(
            &observed_b.particle_velocity_m_s[0],
            &observed_a.particle_velocity_m_s[0].map(|v| v * phase),
            2e-8,
        );
    }
}

#[test]
fn topology_region_source_and_resource_refusals_preserve_reusable_owner() {
    let states = fluids();
    let outer = sphere(1.0, 0, [0.0; 3]);
    let core = sphere(0.35, 0, [0.0; 3]);
    let crossing = sphere(0.35, 0, [0.8, 0.0, 0.0]);
    let gate = CancelGate::new();
    assert!(
        CoreShellFluidInterface::new(
            &outer,
            &crossing,
            views(&states),
            OMEGA,
            0.01,
            Options::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        CoreShellFluidInterface::new(
            &outer,
            &core,
            views(&states),
            OMEGA,
            0.6,
            Options::default(),
            &gate
        )
        .is_err()
    );
    assert!(matches!(
        CoreShellFluidInterface::new(
            &outer,
            &core,
            views(&states),
            OMEGA,
            0.01,
            Options {
                maximum_panels: 39,
                ..Options::default()
            },
            &gate
        ),
        Err(TransmissionError::Helmholtz(HelmholtzError::WorkCap { .. }))
    ));
    assert!(matches!(
        CoreShellFluidInterface::new(
            &outer,
            &core,
            views(&states),
            OMEGA,
            0.01,
            Options {
                maximum_bytes: 1,
                ..Options::default()
            },
            &gate
        ),
        Err(TransmissionError::Helmholtz(
            HelmholtzError::CacheBudget { .. }
        ))
    ));
    // A closed triangulation with an inward dent must not pass the deliberately
    // conservative convex outer-shell admission, even if the core fits.
    let original = sphere(1.0, 1, [0.0; 3]);
    let mut triangles = original.triangles().unwrap().to_vec();
    let vertex = triangles[0][0];
    for point in triangles.iter_mut().flatten() {
        if *point == vertex {
            *point = vertex.map(|coordinate| 0.35 * coordinate);
        }
    }
    let concave = SpherePanels::from_triangles(triangles).unwrap();
    assert!(
        CoreShellFluidInterface::new(
            &concave,
            &core,
            views(&states),
            OMEGA,
            0.01,
            Options::default(),
            &gate
        )
        .is_err()
    );
    let owner = CoreShellFluidInterface::new(
        &outer,
        &core,
        views(&states),
        OMEGA,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    // Being outside the core is insufficient: incident source singularities
    // must also remain outside the outer boundary of the complete medium.
    let shell_source = [Monopole {
        position_m: [0.0, 0.0, 0.6],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&shell_source), 0.01, &gate)
            .is_err()
    );
    let response = owner.scatter_field(plane(), 0.01, &gate).unwrap();
    assert!(
        response
            .pressure_velocity_at(
                CoreShellRegion::Shell,
                &[[0.0; 3]],
                0.01,
                near_field::Options::default(),
                &gate
            )
            .is_err()
    );
    assert!(
        response
            .pressure_velocity_at(
                CoreShellRegion::Core,
                &[[0.0, 0.0, 0.6]],
                0.01,
                near_field::Options::default(),
                &gate
            )
            .is_err()
    );
    assert!(
        response
            .pressure_velocity_at(
                CoreShellRegion::Shell,
                &[[0.0, 0.0, 0.6]],
                0.01,
                near_field::Options {
                    maximum_kernel_evaluations: 80,
                    ..near_field::Options::default()
                },
                &gate
            )
            .is_err()
    );
    assert!(response.far_field(&DIRECTIONS, 1, &gate).is_err());
    assert!(
        response
            .pressure_velocity_at(
                CoreShellRegion::Shell,
                &[[0.0, 0.0, 0.6]],
                0.01,
                near_field::Options::default(),
                &gate
            )
            .is_ok()
    );
    assert!(owner.scatter_field(plane(), 0.01, &gate).is_ok());
    gate.request();
    assert!(matches!(
        owner.scatter_field(plane(), 0.01, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        response.pressure_velocity_at(
            CoreShellRegion::Shell,
            &[[0.0, 0.0, 0.6]],
            0.01,
            near_field::Options::default(),
            &gate
        ),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        response.far_field(&DIRECTIONS, 1_000_000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        CoreShellFluidInterface::new(
            &outer,
            &core,
            views(&states),
            OMEGA,
            0.01,
            Options::default(),
            &gate
        ),
        Err(TransmissionError::Helmholtz(HelmholtzError::Cancelled))
    ));
}
