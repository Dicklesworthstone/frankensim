use super::super::{P1ContactError, P1ContactInterface};
use super::{P1CurvedContactInterface, P1CurvedContactOptions};
use crate::charts::Ray;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{
    BeerLambertAbsorption, DielectricCoatingStack, DielectricGlass, GlassProvenance, SellmeierIor,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{GradientIndexError, GradientIndexOptions};
use crate::material_state::thermo_optic::transient::mesh::P1MeshError;
use crate::material_state::thermo_optic::transient::{
    P1ConductionElement, P1ConductionTemperatureSnapshot,
};
use crate::material_state::thermo_optic::{MaterialThermoOpticField, ThermoOpticFieldError};
use crate::polarization::Stokes;
use fs_blake3::ContentHash;
use fs_conduction::ConductionMesh;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};
use fs_rep_mesh::TetComplex;

type Binding<'a> = MaterialThermoOpticField<'a, P1ConductionElement<'a>>;
const AXIS: Vec3 = Vec3::new(0.0, 1.0, 0.0);

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn same_ray(a: Ray, b: Ray, tolerance: f64) {
    vector(
        a.origin.delta_from(b.origin),
        Vec3::new(0.0, 0.0, 0.0),
        tolerance,
    );
    vector(a.dir, b.dir, tolerance);
}
fn options() -> P1CurvedContactOptions {
    P1CurvedContactOptions {
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                initial_step_m: 0.01,
                max_step_m: 0.04,
                min_step_m: 1e-13,
                relative_tolerance: 2e-12,
                position_tolerance_m: 2e-12,
                direction_tolerance: 2e-12,
                optical_path_tolerance_m: 2e-12,
                optical_depth_tolerance: 2e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
        max_incident_distance_m: 1.0,
        max_outgoing_distance_m: 0.09,
    }
}
fn launch(sine: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, -0.1, -0.2),
        dir: Vec3::new(sine, 0.0, (1.0 - sine * sine).sqrt()),
    }
}
fn coefficient(n: f64) -> f64 {
    0.5 * (n - 1.0 / n) * 0.003
}
fn index(n: f64, temperature: f64) -> f64 {
    n + coefficient(n) * (temperature - 300.0)
}
fn fixture(
    reverse: bool,
    values: [f64; 8],
    indices: [f64; 2],
    absorption: [f64; 2],
    f: impl FnOnce(&CancelGate, &Cx<'_>, [&Binding<'_>; 2]),
) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 116,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let positions = vec![
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut tets = vec![[0, 1, 2, 3], [4, 5, 6, 7]];
        if reverse {
            for tet in &mut tets {
                tet.swap(0, 1);
            }
        }
        let mesh = ConductionMesh::new(TetComplex::from_tets(8, tets), positions).unwrap();
        let snapshot = P1ConductionTemperatureSnapshot::try_new(&cx, &mesh, &values).unwrap();
        let elements = [snapshot.element(0).unwrap(), snapshot.element(1).unwrap()];
        // Nondispersive reference laws make both phase and group integrals
        // independently available from the affine-index analytical solution.
        let glasses = core::array::from_fn::<_, 2, _>(|i| {
            DielectricGlass::new(
                SellmeierIor::try_new_in_band(
                    [indices[i] * indices[i] - 1.0, 0.0, 0.0],
                    [0.0; 3],
                    (1000.0, 2000.0),
                )
                .unwrap(),
                BeerLambertAbsorption::try_constant(absorption[i]).unwrap(),
                GlassProvenance::Custom,
            )
        });
        let law = SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                d: [0.003, 0.0, 0.0],
                e: [0.0; 2],
                resonance_um: 0.0,
            },
            300.0,
            [200.0, 450.0],
            [1000.0, 2000.0],
            ContentHash([116; 32]),
        )
        .unwrap();
        let a = MaterialThermoOpticField::try_new(
            &glasses[0],
            &law,
            &elements[0],
            ContentHash([117; 32]),
        )
        .unwrap();
        let b = MaterialThermoOpticField::try_new(
            &glasses[1],
            &law,
            &elements[1],
            ContentHash([117; 32]),
        )
        .unwrap();
        f(&gate, &cx, [&a, &b]);
    });
}

// Independent closed form for n(z)=n0+g*(z-z0): d(n*u)/ds=(0,0,g).
// Horizontal momentum k is constant and vertical momentum q grows linearly.
fn affine_z(ray: Ray, n0: f64, g: f64, distance: f64) -> (Ray, f64) {
    let k = n0 * ray.dir.x;
    let q0 = n0 * ray.dir.z;
    let q = q0 + g * distance;
    let n = k.hypot(q);
    let primitive = |v: f64| 0.5 * (v * k.hypot(v) + k * k * (v / k).asinh());
    (
        Ray {
            origin: ray.origin.offset(Vec3::new(
                k / g * ((q / k).asinh() - (q0 / k).asinh()),
                0.0,
                distance * (2.0 * q0 + g * distance) / (n + n0),
            )),
            dir: Vec3::new(k / n, 0.0, q / n),
        },
        (primitive(q) - primitive(q0)) / g,
    )
}

#[test]
fn g1_curved_contact_preserves_affine_rays_snell_and_material_integrals() {
    for reverse in [false, true] {
        fixture(
            reverse,
            [310.0, 310.0, 310.0, 280.0, 290.0, 290.0, 290.0, 270.0],
            [1.3, 1.7],
            [0.2, 0.6],
            |_, cx, bindings| {
                let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
                let ray = launch(0.55);
                let split = contact.split(cx, 0, ray, AXIS, 1500.0, options()).unwrap();
                let incident = split.incident;
                let n0 = index(1.3, 304.0);
                let gi = coefficient(1.3) * 30.0;
                let gt = coefficient(1.7) * -20.0;
                let (expected_i, phase_i) = affine_z(ray, n0, gi, incident.distance_m);
                same_ray(incident.ray, expected_i, 5e-11);
                near(incident.optical_path_length_m, phase_i, 5e-11);
                assert!(incident.ray.origin.z < 0.0);
                assert!(incident.ray.origin.z.abs() <= options().event.surface_tolerance_m);
                let ni = index(1.3, 310.0 + 30.0 * incident.ray.origin.z);
                let nt = index(1.7, 290.0 - 20.0 * incident.ray.origin.z);
                let k = n0 * ray.dir.x;
                near(ni * incident.ray.dir.x, k, 5e-11);
                assert!((incident.ray.dir.x - ray.dir.x).abs() > 1e-4);
                let si = incident.ray.dir.x;
                let ci = incident.ray.dir.z;
                let st = ni / nt * si;
                let ct = (1.0 - st * st).sqrt();
                let start_t = Ray {
                    origin: incident.ray.origin,
                    dir: Vec3::new(st, 0.0, ct),
                };
                let (expected_t, phase_t) = affine_z(start_t, nt, gt, 0.09);
                let transmitted = split.transmitted.unwrap();
                let t = transmitted.propagation();
                same_ray(t.ray, expected_t, 8e-11);
                near(t.final_sample.phase_index() * t.ray.dir.x, k, 5e-11);
                near(t.distance_m, incident.distance_m + 0.09, 3e-15);
                near(t.optical_path_length_m, phase_i + phase_t, 8e-11);
                near(
                    t.group_optical_path_length_m.unwrap(),
                    phase_i + phase_t,
                    8e-11,
                );
                near(
                    t.optical_depth,
                    0.2 * incident.distance_m + 0.6 * 0.09,
                    2e-13,
                );
                let rs = (ni * ci - nt * ct) / (ni * ci + nt * ct);
                let rp = (nt * ci - ni * ct) / (nt * ci + ni * ct);
                let reflectance = 0.5 * (rs * rs + rp * rp);
                near(split.fresnel.reflectance, reflectance, 2e-13);
                near(
                    transmitted.power_fraction(),
                    (1.0 - reflectance) * (-t.optical_depth).exp(),
                    2e-13,
                );
                let reflected = split.reflected.unwrap();
                let start_r = Ray {
                    origin: incident.ray.origin,
                    dir: Vec3::new(si, 0.0, -ci),
                };
                let (expected_r, phase_r) = affine_z(start_r, ni, gi, 0.09);
                same_ray(reflected.propagation().ray, expected_r, 8e-11);
                near(
                    reflected.propagation().optical_path_length_m,
                    phase_i + phase_r,
                    8e-11,
                );
                near(
                    reflected.power_fraction(),
                    reflectance * (-0.2 * (incident.distance_m + 0.09)).exp(),
                    2e-13,
                );
                assert_eq!(
                    (reflected.element_index(), transmitted.element_index()),
                    (0, 1)
                );
                assert_eq!(transmitted.arrival_face(), None);
                assert_eq!(split.coating_absorbed_power_fraction, 0.0);
                vector(transmitted.frame().horizontal(), AXIS, 2e-12);
            },
        );
    }
}

#[test]
fn g1_uniform_curved_contact_converges_to_the_exact_straight_owner() {
    fixture(
        false,
        [300.0; 8],
        [1.2, 1.8],
        [0.2, 0.6],
        |_, cx, bindings| {
            let exact = P1ContactInterface::try_new(cx, bindings, [3, 3])
                .unwrap()
                .split(cx, 0, launch(0.6), AXIS, 1500.0, 0.09)
                .unwrap();
            let curved = P1CurvedContactInterface::try_new(cx, bindings, [3, 3])
                .unwrap()
                .split(cx, 0, launch(0.6), AXIS, 1500.0, options())
                .unwrap();
            near(curved.incident.distance_m, exact.incident.distance_m, 2e-8);
            assert!(
                curved.incident.ray.origin.z < 0.0,
                "arrival must retain its integrated position"
            );
            for (a, b) in [
                (curved.reflected.unwrap(), exact.reflected.unwrap()),
                (curved.transmitted.unwrap(), exact.transmitted.unwrap()),
            ] {
                same_ray(a.propagation().ray, b.propagation().ray, 3e-8);
                near(
                    a.propagation().optical_path_length_m,
                    b.propagation().optical_path_length_m,
                    4e-8,
                );
                near(
                    a.propagation().group_optical_path_length_m.unwrap(),
                    b.propagation().group_optical_path_length_m.unwrap(),
                    4e-8,
                );
                near(a.power_fraction(), b.power_fraction(), 2e-8);
            }
            assert!(curved.boundary_evaluations > 0);
            assert!(curved.attempted_steps > 0);
            assert!(curved.contact_resolution_m > 0.0 && curved.contact_resolution_m < 3e-8);
        },
    );
}

#[test]
fn g1_tightening_contact_tolerance_converges_across_a_thermal_jump() {
    fixture(
        false,
        [310.0, 310.0, 310.0, 280.0, 290.0, 290.0, 290.0, 270.0],
        [1.3, 1.7],
        [0.2, 0.6],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            let ray = launch(0.55);
            let n0 = index(1.3, 304.0);
            let ni = index(1.3, 310.0);
            let nt = index(1.7, 290.0);
            let gi = coefficient(1.3) * 30.0;
            let k = n0 * ray.dir.x;
            let exact_distance = ((ni * ni - k * k).sqrt() - n0 * ray.dir.z) / gi;
            let (arrival, phase_i) = affine_z(ray, n0, gi, exact_distance);
            let (_, phase_t) = affine_z(
                Ray {
                    origin: arrival.origin,
                    dir: Vec3::new(k / nt, 0.0, (1.0 - (k / nt).powi(2)).sqrt()),
                },
                nt,
                coefficient(1.7) * -20.0,
                0.09,
            );
            let mut loose = options();
            loose.event.surface_tolerance_m = 1e-4;
            let a = contact.split(cx, 0, ray, AXIS, 1500.0, loose).unwrap();
            let b = contact.split(cx, 0, ray, AXIS, 1500.0, options()).unwrap();
            let error = |p: f64| (p - phase_i - phase_t).abs();
            let loose_error = error(a.transmitted.unwrap().propagation().optical_path_length_m);
            let tight_error = error(b.transmitted.unwrap().propagation().optical_path_length_m);
            assert!(loose_error > 1e-7);
            assert!(
                tight_error < loose_error * 0.01,
                "{tight_error} versus {loose_error}"
            );
            assert!(tight_error < 3e-8);
            assert!(b.contact_resolution_m < a.contact_resolution_m * 0.01);
            near(b.incident.distance_m, exact_distance, 2e-8);
        },
    );
}

fn film_source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([118; 32]), ConductorDataStatus::Representative).unwrap()
}
fn absorbing_film(n: f64, k: f64, thickness_m: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(1000.0 + 125.0 * i as f64, n, k).unwrap()
        }),
        thickness_m,
        125.0,
        film_source(),
    )
    .unwrap()
}

#[test]
fn g1_curved_coating_retains_independent_polarized_power_and_frame_covariance() {
    let transparent = CoherentDielectricCoating::try_new(
        SellmeierIor::try_new_in_band([1.4 * 1.4 - 1.0, 0.0, 0.0], [0.0; 3], (1000.0, 2000.0))
            .unwrap(),
        225e-9,
        film_source(),
    )
    .unwrap();
    let absorbing = absorbing_film(2.1, 0.12, 180e-9);
    let films = [
        CoatingLayer::Dielectric(&transparent),
        CoatingLayer::Absorbing(&absorbing),
    ];
    let stack = DielectricCoatingStack::try_mixed_in_band(&films, (1000.0, 2000.0)).unwrap();
    fixture(
        false,
        [300.0, 300.0, 300.0, 280.0, 300.0, 300.0, 300.0, 340.0],
        [1.2, 1.8],
        [0.0; 2],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3])
                .unwrap()
                .with_coating_stack(stack);
            // Set conserved tangential momentum so the mathematical contact has
            // sin(theta_i)=0.6. The numerical collar contributes O(surface tolerance).
            let ray = launch(0.72 / index(1.2, 296.0));
            let mut opts = options();
            opts.max_outgoing_distance_m = 0.0;
            let input = Stokes::try_new(1.0, 0.3, 0.4, 0.2).unwrap();
            let a = contact
                .split_polarized(cx, 0, ray, AXIS, input, 1500.0, opts)
                .unwrap();
            // Independent characteristic matrices at n_i=1.2,n_t=1.8,k_parallel=.72.
            // Columns are s/p; entries are R,T,A. The same physical films are used
            // by the exact-owner controls, without using that owner's output here.
            let s = [0.06570139486987908, 0.7729547986894787, 0.16134380644064228];
            let p = [0.0343781678055786, 0.7980701299636407, 0.16755170223078075];
            near(
                a.reflected.unwrap().power_fraction(),
                0.65 * s[0] + 0.35 * p[0],
                3e-9,
            );
            near(
                a.transmitted.unwrap().power_fraction(),
                0.65 * s[1] + 0.35 * p[1],
                3e-9,
            );
            near(
                a.coating_absorbed_power_fraction,
                0.65 * s[2] + 0.35 * p[2],
                3e-9,
            );
            near(
                a.reflected.unwrap().power_fraction()
                    + a.transmitted.unwrap().power_fraction()
                    + a.coating_absorbed_power_fraction,
                1.0,
                2e-13,
            );
            let b = contact
                .split_polarized(
                    cx,
                    0,
                    ray,
                    Vec3::new(-ray.dir.z, 0.0, ray.dir.x),
                    Stokes::try_new(1.0, -0.3, -0.4, 0.2).unwrap(),
                    1500.0,
                    opts,
                )
                .unwrap();
            for (a, b) in [
                (a.reflected.unwrap(), b.reflected.unwrap()),
                (a.transmitted.unwrap(), b.transmitted.unwrap()),
            ] {
                let a_stokes = a.stokes();
                let b_stokes = b.stokes();
                for (a, b) in [a_stokes.i, a_stokes.q, a_stokes.u, a_stokes.v]
                    .into_iter()
                    .zip([b_stokes.i, b_stokes.q, b_stokes.u, b_stokes.v])
                {
                    near(a, b, 2e-13);
                }
                vector(a.frame().horizontal(), b.frame().horizontal(), 2e-13);
            }
        },
    );
}

#[test]
fn g1_curved_bulk_keeps_stokes_in_its_parallel_transported_frame() {
    // The upper gradient has a y component, so its curved outgoing trajectory
    // leaves the original incidence plane and rotates its transverse frame.
    fixture(
        false,
        [300.0, 300.0, 300.0, 280.0, 240.0, 280.0, 340.0, 290.0],
        [1.2, 1.8],
        [0.0, 0.4],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            let input = Stokes::try_new(1.0, 0.2, 0.4, 0.5).unwrap();
            let mut opts = options();
            opts.max_outgoing_distance_m = 0.0;
            let start = contact
                .split_polarized(cx, 0, launch(0.5), AXIS, input, 1500.0, opts)
                .unwrap()
                .transmitted
                .unwrap();
            opts.max_outgoing_distance_m = 0.2;
            let end = contact
                .split_polarized(cx, 0, launch(0.5), AXIS, input, 1500.0, opts)
                .unwrap()
                .transmitted
                .unwrap();
            assert!(end.propagation().ray.dir.y.abs() > 1e-3);
            assert!(start.frame().horizontal().dot(end.frame().horizontal()) < 1.0 - 1e-6);
            let beer = (-0.4_f64 * 0.2).exp();
            for (a, b) in [
                start.stokes().i,
                start.stokes().q,
                start.stokes().u,
                start.stokes().v,
            ]
            .into_iter()
            .zip([
                end.stokes().i,
                end.stokes().q,
                end.stokes().u,
                end.stokes().v,
            ]) {
                near(a * beer, b, 2e-12);
            }
            near(start.stokes().degree(), end.stokes().degree(), 2e-13);
        },
    );
}

#[test]
fn g1_a_curved_ray_can_return_to_its_departure_face_and_cross_again() {
    fixture(
        false,
        [300.0, 300.0, 300.0, 300.0, 440.0, 440.0, 440.0, 210.0],
        [1.675, 1.5],
        [0.0; 2],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            let mut ray = launch(0.9995);
            ray.origin = Point3::new(-0.2, -0.1, -0.001);
            let mut opts = options();
            opts.max_outgoing_distance_m = 1.0;
            let split = contact.split(cx, 0, ray, AXIS, 1500.0, opts).unwrap();
            let returned = split.transmitted.unwrap();
            assert_eq!(returned.arrival_face(), Some(3));
            assert_eq!(returned.element_index(), 1);
            assert!(returned.propagation().ray.dir.z < -0.02);
            assert!(returned.propagation().ray.origin.z > 0.0);
            assert!(returned.propagation().ray.origin.z <= opts.event.surface_tolerance_m);
            // n(z)=1.675-.2875z turns at z≈.00291 m; its full return arc
            // is 2*n(0)*cos(theta)/.2875, safely inside this finite face.
            let return_arc = 2.0 * 1.675 * ray.dir.z / 0.2875;
            near(
                returned.propagation().distance_m - split.incident.distance_m,
                return_arc,
                2e-6,
            );
            opts.max_outgoing_distance_m = 0.01;
            let next = contact.split_arrival(cx, &returned, opts).unwrap();
            assert_eq!(next.incident, *returned.propagation());
            let exited = next.transmitted.unwrap();
            assert_eq!(exited.element_index(), 0);
            vector(
                exited.propagation().ray.dir,
                Vec3::new(ray.dir.x, 0.0, -ray.dir.z),
                2e-9,
            );
            near(
                exited.propagation().distance_m,
                returned.propagation().distance_m + 0.01,
                3e-15,
            );
            assert!(exited.contact_resolution_m() >= returned.contact_resolution_m());
        },
    );
}

#[test]
fn g3_curved_contact_shares_work_limits_and_preserves_cancellation() {
    fixture(
        false,
        [310.0, 310.0, 310.0, 280.0, 290.0, 290.0, 290.0, 270.0],
        [1.3, 1.7],
        [0.0; 2],
        |gate, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            let ray = launch(0.55);
            let baseline = contact.split(cx, 0, ray, AXIS, 1500.0, options()).unwrap();
            assert!(baseline.attempted_steps > 1 && baseline.boundary_evaluations > 1);
            assert!(matches!(
                contact.split_arrival(cx, &baseline.transmitted.unwrap(), options()),
                Err(P1ContactError::ForeignArrival)
            ));
            assert!(matches!(
                contact.split(
                    cx,
                    0,
                    Ray {
                        origin: ray.origin,
                        dir: ray.dir.scale(-1.0)
                    },
                    AXIS,
                    1500.0,
                    options(),
                ),
                Err(P1ContactError::ContactNotFirstFace { .. })
            ));
            let mut opts = options();
            opts.event.integration.max_steps = baseline.attempted_steps;
            opts.event.max_boundary_queries = baseline.boundary_evaluations;
            assert!(contact.split(cx, 0, ray, AXIS, 1500.0, opts).is_ok());
            opts.event.integration.max_steps -= 1;
            assert!(matches!(
                contact.split(cx, 0, ray, AXIS, 1500.0, opts),
                Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::StepBudgetExceeded
                )))
            ));
            opts = options();
            opts.event.max_boundary_queries = baseline.boundary_evaluations - 1;
            assert!(matches!(
                contact.split(cx, 0, ray, AXIS, 1500.0, opts),
                Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::BoundaryQueryBudgetExceeded
                )))
            ));
            opts = options();
            opts.max_incident_distance_m = 0.01;
            assert!(matches!(
                contact.split(cx, 0, ray, AXIS, 1500.0, opts),
                Err(P1ContactError::IncidentDistanceLimit)
            ));
            assert!(contact.split(cx, 0, ray, AXIS, 900.0, options()).is_err());
            gate.request();
            assert!(matches!(
                contact.split(cx, 0, ray, AXIS, 1500.0, options()),
                Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::Cancelled
                )))
            ));
        },
    );
}

#[test]
fn g3_curved_contact_refuses_unresolved_edges_and_out_of_law_collar_temperatures() {
    fixture(
        false,
        [300.0; 8],
        [1.2, 1.8],
        [0.0; 2],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            let mut ray = launch(0.0);
            let length = (0.9_f64 * 0.9 + 0.2 * 0.2).sqrt();
            ray.dir = Vec3::new(0.0, -0.9 / length, 0.2 / length);
            // Both incident coordinates are strictly interior; only the intended
            // contact lies on the shared edge y=-1 of two original faces.
            assert!(matches!(
                contact.split(cx, 0, ray, Vec3::new(1.0, 0.0, 0.0), 1500.0, options()),
                Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::BoundaryResolutionLimit
                ))) | Err(P1ContactError::UnresolvedFace)
            ));
        },
    );
    fixture(
        false,
        [300.0, 300.0, 300.0, 300.0, 200.0, 200.0, 200.0, 450.0],
        [1.2, 1.8],
        [0.0; 2],
        |_, cx, bindings| {
            let contact = P1CurvedContactInterface::try_new(cx, bindings, [3, 3]).unwrap();
            // Original target vertices are in the accepted law domain. Its affine
            // continuation at the still-incident numerical endpoint falls below
            // 200 K and must refuse, rather than clamp or average temperatures.
            assert!(matches!(
                contact.split(cx, 0, launch(0.4), AXIS, 1500.0, options()),
                Err(P1ContactError::Mesh(P1MeshError::Material(
                    ThermoOpticFieldError::Field(GradientIndexError::OutsideFieldDomain)
                ))) | Err(P1ContactError::Mesh(P1MeshError::Propagation(
                    GradientIndexError::OutsideFieldDomain
                )))
            ));
        },
    );
}
