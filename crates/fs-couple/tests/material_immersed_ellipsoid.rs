//! Source-derived nonspherical inertia, geometry and torque on real D3Q19 steps.
use fs_couple::material_immersed::{
    AdditionalMaterialEllipsoid, BodyStateSi, ImmersedDomain, RetainedFluidEllipsoids,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{Se3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{E3, Q3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

const AXES: [f64; 3] = [0.014, 0.009, 0.006];
const DX: f64 = 0.01;
const DT: f64 = 0.1;
const DENSITY: f64 = 3000.0;
type Model<'a> = RetainedFluidEllipsoids<'a, ResolvedLiquidState, 1>;

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn rotation(v: [f64; 3]) -> So3 {
    So3::exp(So3Tangent::new(vector(v))).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, delta={:.4e}",
        a - b
    );
}
fn card(properties: &[(&str, Dims, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: if name == "thermal_conductivity" {
                    PropertyKey::with_quantity(
                        name,
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity,
                            ValueForm::Static,
                        )),
                    )
                } else {
                    PropertyKey::new(name, dims)
                },
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic ellipsoid coupling fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic source".into(),
            phase: "fixture".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn fluid(density: f64, viscosity: f64) -> ResolvedLiquidState {
    resolve_liquid_state(
        &card(&[
            ("density", Density::DIMS, density),
            ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, viscosity),
            (
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                1000.0,
            ),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        ]),
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn source(body_from_shape: So3) -> ResolvedBody {
    let geometry = VolumeMoments::ellipsoid(vector(AXES))
        .unwrap()
        .transformed(Se3::from_parts(body_from_shape, Vec3::new(0.0, 0.0, 0.0)).unwrap())
        .unwrap();
    resolve_body(
        &[RegionBinding {
            key: "ellipsoid",
            geometry,
            material: &card(&[("density", Density::DIMS, DENSITY)]),
            point: &point(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn model<'a>(
    fluid: &'a ResolvedLiquidState,
    material: &'a ResolvedBody,
    body_from_shape: So3,
) -> Model<'a> {
    Model::new(
        fluid,
        ImmersedDomain {
            cells: [8; 3],
            periodic: [true; 3],
            cell_size_m: DX,
            timestep_s: DT,
            maximum_speed_m_s: 0.001,
            maximum_cells: 512,
            samples_per_axis: 4,
        },
        [AdditionalMaterialEllipsoid::new(material, AXES, body_from_shape).unwrap()],
        [BodyStateSi {
            center_m: [0.0415, 0.0387, 0.0431],
            orientation: rotation([0.4, -0.2, 0.3])
                .compose(body_from_shape.inverse())
                .unwrap(),
            velocity_m_s: [0.0001, -0.00007, 0.00003],
            angular_momentum_body_kg_m2_s: array(
                body_from_shape.rotate(vector([2e-9, -1e-9, 3e-9])).unwrap(),
            ),
        }],
        &CancelGate::new(),
    )
    .unwrap()
}
fn populations(model: &Model<'_>) -> Vec<[f64; Q3]> {
    (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| (0..8).map(move |x| model.grid().populations(x, y, z)))
        })
        .collect()
}
fn momentum_si(model: &Model<'_>, origin: [f64; 3]) -> (Vec3, Vec3) {
    let mut linear = Vec3::new(0.0, 0.0, 0.0);
    let mut angular = linear;
    let origin = vector(origin);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let position = Vec3::new(x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5).scale(DX);
                let f = model.grid().populations(x, y, z);
                for q in 0..Q3 {
                    let e = [E3[q].0, E3[q].1, E3[q].2];
                    let impulse = vector(e.map(|c| model.impulse_ns(f[q] * f64::from(c)).unwrap()));
                    linear = linear + impulse;
                    angular = angular + (position - origin).cross(impulse);
                }
            }
        }
    }
    let body = model.state_si(0).unwrap();
    let impulse = vector(body.velocity_m_s).scale(model.sources()[0].material().inertia().mass());
    (
        linear + impulse,
        angular
            + (vector(body.center_m) - origin).cross(impulse)
            + body
                .orientation
                .rotate(vector(body.angular_momentum_body_kg_m2_s))
                .unwrap(),
    )
}

#[test]
fn ellipsoid_moments_and_sampled_geometry_share_the_complete_material_frame() {
    let volume = 4.0 * core::f64::consts::PI / 3.0 * AXES[0] * AXES[1] * AXES[2];
    let principal = [
        (AXES[1].powi(2) + AXES[2].powi(2)) * volume / 5.0,
        (AXES[0].powi(2) + AXES[2].powi(2)) * volume / 5.0,
        (AXES[0].powi(2) + AXES[1].powi(2)) * volume / 5.0,
    ];
    let basis = rotation([0.2, 0.4, -0.3]);
    let body = source(basis);
    close(body.inertia().mass(), DENSITY * volume, 1e-17);
    let r = basis.matrix();
    let actual = body.inertia().inertia_com().m;
    for i in 0..3 {
        for j in 0..3 {
            let expected: f64 = (0..3)
                .map(|k| DENSITY * principal[k] * r.m[3 * i + k] * r.m[3 * j + k])
                .sum();
            close(actual[3 * i + j], expected, 3e-22);
        }
    }
    assert!(actual[1].abs() > 1e-9);
    assert!(AdditionalMaterialEllipsoid::new(&body, AXES, basis).is_ok());
    assert!(AdditionalMaterialEllipsoid::new(&body, AXES, So3::identity()).is_err());
    assert!(
        AdditionalMaterialEllipsoid::new(&body, [AXES[0] * 2.0, AXES[1], AXES[2]], basis).is_err()
    );
    for invalid in [[0.0, 1.0, 1.0], [-1.0, 1.0, 1.0], [f64::INFINITY, 1.0, 1.0]] {
        assert!(VolumeMoments::ellipsoid(vector(invalid)).is_err());
    }
    let sphere = VolumeMoments::sphere(0.01).unwrap();
    let equal_axes = VolumeMoments::ellipsoid(Vec3::new(0.01, 0.01, 0.01)).unwrap();
    close(sphere.volume_m3(), equal_axes.volume_m3(), 1e-20);
    for (a, b) in sphere
        .inertia_m5()
        .m
        .into_iter()
        .zip(equal_axes.inertia_m5().m)
    {
        close(a, b, 1e-24);
    }
}

#[test]
fn rotating_full_tensor_ellipsoid_obeys_si_fluid_and_body_angular_balance() {
    let fluid = fluid(1000.0, 0.11);
    let basis = rotation([0.2, 0.4, -0.3]);
    let material = source(basis);
    let mut model = model(&fluid, &material, basis);
    let origin = [-0.011, 0.023, -0.037];
    let force = [2e-6, -1e-6, 3e-7];
    let torque = [3e-9, -2e-9, 1e-9];
    let held = [[
        force[0], force[1], force[2], torque[0], torque[1], torque[2],
    ]];
    assert!(core::ptr::eq(model.sources()[0].material(), &material));
    for (lattice, si) in model.bodies_lattice()[0]
        .inertia
        .matrix()
        .m
        .into_iter()
        .zip(material.inertia().inertia_com().m)
    {
        close(model.inertia_kg_m2(lattice).unwrap(), si, 3e-22);
    }
    let mut total_hydrodynamic_couple = 0.0;
    for _ in 0..4 {
        let before = momentum_si(&model, origin);
        let report = model
            .step_forces_si(held, origin, &CancelGate::new())
            .unwrap();
        let after = momentum_si(&model, origin);
        let arm = vector(report.bodies[0].midpoint_center).scale(DX) - vector(origin);
        let imposed = (arm.cross(vector(force)) + vector(torque)).scale(DT);
        let seam = vector(
            report
                .coupling
                .periodic
                .angular_momentum_change
                .map(|h| model.angular_impulse_nms(h).unwrap()),
        );
        for axis in 0..3 {
            close(array(after.0 - before.0)[axis], force[axis] * DT, 3e-16);
            close(
                array(after.1 - before.1 - seam)[axis],
                array(imposed)[axis],
                3e-17,
            );
            close(
                model
                    .angular_impulse_nms(report.bodies[0].momentum_residual[axis + 3])
                    .unwrap(),
                0.0,
                1e-19,
            );
        }
        assert_eq!(report.coupling.lubrication.active_pairs, 0);
        total_hydrodynamic_couple += report.coupling.bodies[0]
            .hydrodynamic
            .body_angular_impulse_about_center
            .into_iter()
            .map(f64::abs)
            .sum::<f64>();
    }
    assert!(total_hydrodynamic_couple > 1e-6);
    assert_eq!(model.accepted_steps(), 4);
}

#[test]
fn body_basis_relabeling_preserves_motion_and_rebind_or_refusal_preserves_state() {
    let original = fluid(1000.0, 0.11);
    let updated = fluid(1000.0, 0.22);
    let incompatible = fluid(1100.0, 0.22);
    let basis = rotation([0.2, 0.4, -0.3]);
    let principal = source(So3::identity());
    let reexpressed = source(basis);
    let mut a = model(&original, &principal, So3::identity());
    let mut b = model(&original, &reexpressed, basis);
    let input = [[2e-6, -1e-6, 3e-7, 3e-9, -2e-9, 1e-9]];
    for _ in 0..3 {
        let _ = a
            .step_forces_si(input, [0.0; 3], &CancelGate::new())
            .unwrap();
        let _ = b
            .step_forces_si(input, [0.0; 3], &CancelGate::new())
            .unwrap();
        for (fa, fb) in populations(&a).into_iter().zip(populations(&b)) {
            for (x, y) in fa.into_iter().zip(fb) {
                close(x, y, 3e-14);
            }
        }
        let sa = a.state_si(0).unwrap();
        let sb = b.state_si(0).unwrap();
        let actual_shape_pose = sb.orientation.compose(basis).unwrap().matrix();
        for (x, y) in sa
            .orientation
            .matrix()
            .m
            .into_iter()
            .zip(actual_shape_pose.m)
        {
            close(x, y, 2e-12);
        }
        for axis in 0..3 {
            close(sa.center_m[axis], sb.center_m[axis], 1e-14);
            close(sa.velocity_m_s[axis], sb.velocity_m_s[axis], 1e-14);
        }
        let ha = sa
            .orientation
            .rotate(vector(sa.angular_momentum_body_kg_m2_s))
            .unwrap();
        let hb = sb
            .orientation
            .rotate(vector(sb.angular_momentum_body_kg_m2_s))
            .unwrap();
        for (x, y) in array(ha).into_iter().zip(array(hb)) {
            close(x, y, 1e-20);
        }
    }
    let fields = populations(&b);
    let bodies = *b.bodies_lattice();
    b.rebind_fluid(&updated, &CancelGate::new()).unwrap();
    assert_eq!(fields, populations(&b));
    assert_eq!(&bodies, b.bodies_lattice());
    assert_eq!(b.accepted_steps(), 3);
    assert!(core::ptr::eq(b.scaling().material(), &updated));
    close(b.grid().viscosity(), 0.22 / 1000.0 * DT / (DX * DX), 1e-16);
    assert!(b.rebind_fluid(&incompatible, &CancelGate::new()).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(b.rebind_fluid(&original, &cancelled).is_err());
    assert!(b.step_forces_si(input, [0.0; 3], &cancelled).is_err());
    assert!(
        b.step_forces_si([[1.0; 6]], [0.0; 3], &CancelGate::new())
            .is_err()
    );
    assert_eq!(fields, populations(&b));
    assert_eq!(&bodies, b.bodies_lattice());
    assert_eq!(b.accepted_steps(), 3);
    assert!(core::ptr::eq(b.scaling().material(), &updated));
    a.rebind_fluid(&updated, &CancelGate::new()).unwrap();
    let _ = a
        .step_forces_si(input, [0.0; 3], &CancelGate::new())
        .unwrap();
    let _ = b
        .step_forces_si(input, [0.0; 3], &CancelGate::new())
        .unwrap();
    for (fa, fb) in populations(&a).into_iter().zip(populations(&b)) {
        for (x, y) in fa.into_iter().zip(fb) {
            close(x, y, 3e-14);
        }
    }
}
