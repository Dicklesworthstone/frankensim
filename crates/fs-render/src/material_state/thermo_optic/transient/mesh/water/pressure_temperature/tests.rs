use super::*;
use crate::charts::Ray;
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, propagate};
use crate::gradient_index::bundle::{BundleOptions, propagate_bundle};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::water_refraction::IapwsWaterRefraction;
use fs_qty::Temperature;
use fs_rep_mesh::TetComplex;
use fs_thermochem::water::if97::If97State;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 111, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx);
    });
}
fn mesh(reverse: bool) -> ConductionMesh {
    let mut tet = [0, 1, 2, 3];
    if reverse { tet.swap(0, 1); }
    ConductionMesh::new(TetComplex::from_tets(4, vec![tet]), vec![
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    ]).unwrap()
}
fn point(x: [f64; 3]) -> Point3 { Point3::new(x[0], x[1], x[2]) }
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}; tolerance {tolerance}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.002, max_step_m: 0.01,
        position_tolerance_m: 1e-12, direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-12, optical_depth_tolerance: 1e-12,
        ..GradientIndexOptions::default()
    }
}

#[test]
fn g1_typed_and_pascal_sources_preserve_original_vertices_and_affine_jets() {
    with_cx(|_, cx| {
        for reverse in [false, true] {
            let mesh = mesh(reverse);
            let temperature = [300.0, 305.0, 297.0, 302.0];
            let pressure_pa = [30e6, 38e6, 24e6, 33e6];
            let pressure = pressure_pa.map(Pressure::new);
            let thermal = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &temperature).unwrap();
            let typed = P1PressureSnapshot::try_new(
                cx, &mesh, NodalPressureValues::Quantities(&pressure), ContentHash([111; 32]),
            ).unwrap();
            let native = P1PressureSnapshot::try_new(
                cx, &mesh, NodalPressureValues::Pascals(&pressure_pa), ContentHash([111; 32]),
            ).unwrap();
            let NodalPressureValues::Quantities(original) = typed.nodal_pressure_values() else { panic!("lost typed source"); };
            assert!(core::ptr::eq(original, &pressure[..]));
            let NodalPressureValues::Pascals(original) = native.nodal_pressure_values() else { panic!("lost Pa source"); };
            assert!(core::ptr::eq(original, &pressure_pa[..]));
            assert_eq!(native.nodal_pressure_values().get(2), Some(Pressure::new(24e6)));
            let a = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([112; 32]), typed).unwrap();
            let b = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([112; 32]), native).unwrap();
            assert_eq!(a.content_identity(), b.content_identity());
            let (a, b) = (a.element(0).unwrap(), b.element(0).unwrap());
            for (vertex, position) in mesh.positions().iter().enumerate() {
                let sample = a.sample_jet(point(*position)).unwrap();
                assert_eq!(sample.pressure_pa, pressure_pa[vertex]);
                assert_eq!(sample.temperature.sample.temperature_k, temperature[vertex]);
                assert_eq!(sample, b.sample_jet(point(*position)).unwrap());
            }
            let p = point([0.2, 0.1, 0.15]);
            let jet = a.sample_jet(p).unwrap();
            assert_eq!(jet, b.sample_jet(p).unwrap());
            near(jet.temperature.sample.temperature_k, 301.0, 1e-12);
            near(jet.pressure_pa, 31.45e6, 1e-7);
            vector_near(jet.temperature.sample.gradient_k_per_m, Vec3::new(5.0, -3.0, 2.0), 1e-12);
            vector_near(jet.pressure_gradient_pa_per_m, Vec3::new(8e6, -6e6, 3e6), 1e-8);
            assert_eq!(jet.temperature.diagonal_k_per_m2, [0.0; 3]);
            assert_eq!(jet.pressure_diagonal_pa_per_m2, [0.0; 3]);
            assert_eq!(jet.pressure_off_diagonal_pa_per_m2, [0.0; 3]);
            assert!(!a.is_spatially_uniform());
        }
    });
}

#[test]
fn g1_joint_p1_optical_hessian_matches_original_scalar_sources_with_mixed_gradients() {
    with_cx(|_, cx| {
        let mesh = mesh(false);
        for (p0, t0, gp, region) in [
            (30e6, 300.0, [8e6, -6e6, 3e6], If97Region::Region1),
            (8e4, 400.0, [2e4, -1.3e4, 0.7e4], If97Region::Region2),
        ] {
            let gt = [5.0, -3.0, 2.0];
            let temperature = [t0, t0 + gt[0], t0 + gt[1], t0 + gt[2]];
            let pressure = [p0, p0 + gp[0], p0 + gp[1], p0 + gp[2]];
            let thermal = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &temperature).unwrap();
            let pressure = P1PressureSnapshot::try_new(
                cx, &mesh, NodalPressureValues::Pascals(&pressure), ContentHash([113; 32]),
            ).unwrap();
            let source = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([114; 32]), pressure).unwrap();
            let element = source.element(0).unwrap();
            let state = element.bind_water(region).unwrap();
            assert_eq!(state.snapshot_identity(), source.content_identity());
            let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
            let optics = state.bind_refraction(&absorption).unwrap();
            let x = [0.18, 0.23, 0.16];
            assert_eq!(state.sample_gradient(point(x)).unwrap(), state.sample_jet(point(x)).unwrap().sample());
            for wavelength in [226.5, 589.0, 1013.98] {
                let original = IapwsWaterRefraction::prepare(wavelength).unwrap();
                // This oracle evaluates the original scalar equations only:
                // no pressure-element or spatial derivative implementation.
                let phase = |y: [f64; 3]| {
                    let t = t0 + gt[0] * y[0] + gt[1] * y[1] + gt[2] * y[2];
                    let p = p0 + gp[0] * y[0] + gp[1] * y[1] + gp[2] * y[2];
                    let native = If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region).unwrap();
                    original.sample_if97_state(&native).unwrap().phase_index
                };
                let jet = optics.prepare(wavelength).unwrap().optical_jet(point(x)).unwrap();
                let g = jet.optical().phase_gradient_per_m();
                let gradient = [g.x, g.y, g.z];
                let (n, h) = (phase(x), 0.002);
                for i in 0..3 {
                    let (mut a, mut b) = (x, x);
                    a[i] += h; b[i] -= h;
                    near(gradient[i], (phase(a) - phase(b)) / (2.0 * h), 2e-9);
                    near(jet.hessian().diagonal()[i], (phase(a) - 2.0 * n + phase(b)) / (h * h), 3e-9);
                }
                for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                    let (mut pp, mut pm, mut mp, mut mm) = (x, x, x, x);
                    pp[i] += h; pp[j] += h; pm[i] += h; pm[j] -= h;
                    mp[i] -= h; mp[j] += h; mm[i] -= h; mm[j] -= h;
                    near(jet.hessian().off_diagonal()[k],
                        (phase(pp) - phase(pm) - phase(mp) + phase(mm)) / (4.0 * h * h), 3e-9);
                }
                assert!(jet.hessian().diagonal().into_iter().any(|h| h.abs() > 1e-10));
            }
        }
    });
}

#[test]
fn g1_original_joint_element_drives_native_rays_and_angular_bundles() {
    with_cx(|_, cx| {
        let mesh = mesh(false);
        let temperature = [300.0, 305.0, 297.0, 302.0];
        let pressure = [30e6, 38e6, 24e6, 33e6];
        let thermal = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &temperature).unwrap();
        let pressure = P1PressureSnapshot::try_new(
            cx, &mesh, NodalPressureValues::Pascals(&pressure), ContentHash([115; 32]),
        ).unwrap();
        let source = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([116; 32]), pressure).unwrap();
        let element = source.element(0).unwrap();
        let state = element.bind_water(If97Region::Region1).unwrap();
        let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
        let optics = state.bind_refraction(&absorption).unwrap();
        let prepared = optics.prepare(589.0).unwrap();
        let ray = Ray { origin: point([0.15, 0.15, 0.1]), dir: Vec3::new(0.0, 0.0, 1.0) };
        let length = 0.2;
        let axis = Vec3::new(1.0, 0.0, 0.0);
        let central = propagate(&prepared, ray, axis, 589.0, length, options()).unwrap();
        let bundle = propagate_bundle(&prepared, ray, axis, 589.0, length, BundleOptions {
            integration: options(),
            position_derivative_tolerance_m: 1e-12,
            direction_derivative_tolerance: 1e-12,
        }).unwrap();
        let p = bundle.propagation();
        vector_near(p.ray.origin.delta_from(central.ray.origin), Vec3::new(0.0, 0.0, 0.0), 3e-10);
        vector_near(p.ray.dir, central.ray.dir, 3e-10);
        near(p.optical_path_length_m, central.optical_path_length_m, 3e-10);
        near(p.group_optical_path_length_m.unwrap(), central.group_optical_path_length_m.unwrap(), 3e-10);
        near(p.power_transmittance(), (-0.02 * length).exp(), 1e-12);
        assert!(p.ray.dir.x > 1e-6 && p.ray.dir.y < -1e-6);
        assert!((p.group_optical_path_length_m.unwrap() - p.optical_path_length_m).abs() > 1e-5);
        let epsilon = 2e-5;
        for (basis, actual) in bundle.launch_basis().into_iter().zip(bundle.angular_derivatives()) {
            let neighbour = |sign: f64| {
                let d = Vec3::new(ray.dir.x + sign * epsilon * basis.x,
                    ray.dir.y + sign * epsilon * basis.y, ray.dir.z + sign * epsilon * basis.z);
                let r = Ray { dir: d.scale(1.0 / d.dot(d).sqrt()), ..ray };
                propagate(&prepared, r, axis, 589.0, length, options()).unwrap().ray
            };
            let (a, b) = (neighbour(1.0), neighbour(-1.0));
            vector_near(actual.position_m_per_rad, a.origin.delta_from(b.origin).scale(0.5 / epsilon), 3e-7);
            vector_near(actual.direction_per_rad,
                Vec3::new(a.dir.x - b.dir.x, a.dir.y - b.dir.y, a.dir.z - b.dir.z).scale(0.5 / epsilon), 3e-7);
        }
    });
}

#[test]
fn g3_original_pressure_geometry_identities_and_admission_remain_explicit() {
    with_cx(|gate, cx| {
        let mesh = mesh(false);
        let temperature = [300.0; 4];
        let values = [3e6; 4];
        let thermal = P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &temperature).unwrap();
        let pressure = |id| P1PressureSnapshot::try_new(
            cx, &mesh, NodalPressureValues::Pascals(&values), ContentHash([id; 32]),
        );
        let a = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([117; 32]), pressure(118).unwrap()).unwrap();
        let b = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([119; 32]), pressure(118).unwrap()).unwrap();
        let c = P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([117; 32]), pressure(120).unwrap()).unwrap();
        assert_ne!(a.content_identity(), b.content_identity());
        assert_ne!(a.content_identity(), c.content_identity());
        assert_eq!(a.temperature_snapshot_identity(), ContentHash([117; 32]));
        assert_eq!(a.pressure_snapshot().snapshot_identity(), ContentHash([118; 32]));
        assert!(matches!(pressure(0), Err(PressureTemperatureSnapshotError::MissingPressureIdentity)));
        assert!(matches!(P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([0; 32]), pressure(118).unwrap()),
            Err(PressureTemperatureSnapshotError::MissingTemperatureIdentity)));
        let foreign_mesh = super::tests::mesh(false);
        let foreign = P1PressureSnapshot::try_new(
            cx, &foreign_mesh, NodalPressureValues::Pascals(&values), ContentHash([118; 32]),
        ).unwrap();
        assert!(matches!(P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([117; 32]), foreign),
            Err(PressureTemperatureSnapshotError::ForeignGeometry)));
        assert!(matches!(P1PressureSnapshot::try_new(
            cx, &mesh, NodalPressureValues::Pascals(&values[..3]), ContentHash([118; 32]),
        ), Err(PressureTemperatureSnapshotError::PressureLength { .. })));
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let bad = [invalid, 3e6, 3e6, 3e6];
            assert!(matches!(P1PressureSnapshot::try_new(
                cx, &mesh, NodalPressureValues::Pascals(&bad), ContentHash([118; 32]),
            ), Err(PressureTemperatureSnapshotError::InvalidPressure { vertex: 0 })));
        }
        let element = a.element(0).unwrap();
        assert!(matches!(element.pressure_temperature_gradient(point([1.1, 0.0, 0.0])),
            Err(GradientIndexError::OutsideFieldDomain)));
        assert!(matches!(a.element(1), Err(PressureTemperatureSnapshotError::Temperature(
            TransientTemperatureFieldError::InvalidElement { element: 1 }))));
        gate.request();
        assert!(matches!(pressure(118), Err(PressureTemperatureSnapshotError::Cancelled)));
    });
}
