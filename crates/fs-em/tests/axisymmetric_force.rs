//! Independent Lorentz/dipole force references and native mechanical loads.
//! Synthetic materials exercise the numerical model, not measured device data.
mod common;

use common::{config, limits, material, near};
use fs_em::{
    force::MaxwellForceLimits,
    planar::{
        MagneticError,
        axisymmetric::{
            AxisymmetricMagnetostatics,
            nonlinear::{
                AxisymmetricNonlinearAssembly, AxisymmetricNonlinearControls,
                AxisymmetricNonlinearMagnetostatics,
            },
        },
        nonlinear::IsotropicMagneticMaterial,
        permanent::PlanarRecoilMaterial,
    },
    transfer::PlanarFrame,
};
use fs_evidence::ValidityDomain;
use fs_ga::{Mat3, Se3, Vec3, Wrench};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::MaterialPropertySelection,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    free_floating_forward_dynamics,
};
use fs_qty::QuantitySpec;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::NewtonKrylovConfig;

fn force_limits() -> MaxwellForceLimits {
    MaxwellForceLimits {
        max_cells: 10000,
        max_edges: 20000,
    }
}
fn frame() -> PlanarFrame {
    // Meridional coordinates [r,z] map to global [x,z].
    PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap()
}
fn mesh(n: usize, reverse: bool, angle: f64) -> TriComplex2 {
    let points: Vec<_> = (0..=2 * n)
        .flat_map(|j| (0..=n).map(move |i| [i as f64 / n as f64, -1.0 + j as f64 / n as f64]))
        .collect();
    let mut faces = Vec::new();
    for j in 0..2 * n {
        for i in 0..n {
            let a = (j * (n + 1) + i) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            // Reflection about z=0 preserves the discretized current/magnet,
            // so its own axial force cancels without subtracting two solves.
            faces.extend(if j < n {
                [[a, b, c], [a, c, d]]
            } else {
                [[a, b, d], [b, c, d]]
            });
        }
    }
    if reverse {
        for face in &mut faces {
            face.swap(1, 2);
        }
    }
    TriComplex2::from_indexed_triangles(
        "axisymmetric force fixture",
        points,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}
fn cells(mesh: &TriComplex2, include: impl Fn([f64; 2]) -> bool) -> Vec<usize> {
    mesh.faces()
        .iter()
        .enumerate()
        .filter_map(|(cell, face)| {
            let center = [0, 1].map(|axis| {
                face.iter()
                    .map(|&v| mesh.vertices()[v as usize][axis] / 3.0)
                    .sum()
            });
            include(center).then_some(cell)
        })
        .collect()
}
fn collar(mesh: &TriComplex2) -> Vec<usize> {
    cells(mesh, |[r, z]| r < 0.75 && z.abs() < 0.5)
}
fn coil(mesh: &TriComplex2) -> Vec<usize> {
    cells(mesh, |[r, z]| r > 0.25 && r < 0.5 && z.abs() < 0.25)
}
fn core(mesh: &TriComplex2) -> Vec<usize> {
    cells(mesh, |[r, z]| r < 0.25 && z.abs() < 0.25)
}
fn boundary(mesh: &TriComplex2, offset: f64, gradient: f64) -> Vec<(usize, f64)> {
    mesh.vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 1.0 || p[1].abs() == 1.0)
        .map(|(v, p)| (v, offset + gradient * p[1]))
        .collect()
}
fn linear(n: usize, reverse: bool, angle: f64) -> AxisymmetricMagnetostatics {
    let mesh = mesh(n, reverse, angle);
    let materials = vec![material(1.0); mesh.faces().len()];
    AxisymmetricMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap()
}
fn axial(wrench: Wrench) -> f64 {
    near(wrench.force.x, 0.0);
    near(wrench.force.y, 0.0);
    for value in [wrench.torque.x, wrench.torque.y, wrench.torque.z] {
        near(value, 0.0);
    }
    wrench.force.z
}
fn nonlinear_controls() -> AxisymmetricNonlinearControls {
    AxisymmetricNonlinearControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-11,
            linear_restart: 64,
            max_linear_cycles: 80,
            forcing_minimum: 1e-12,
            forcing_maximum: 0.1,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 40,
        max_krylov_basis_entries: 1000000,
    }
}
fn nonlinear(
    mesh: TriComplex2,
    models: Vec<IsotropicMagneticMaterial>,
) -> AxisymmetricNonlinearMagnetostatics {
    AxisymmetricNonlinearMagnetostatics::new(
        mesh,
        models,
        AxisymmetricNonlinearAssembly {
            limits: limits(),
            quadrature_order: 3,
            max_quadrature_points: 100000,
        },
        &mut || false,
    )
    .unwrap()
}
fn source(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic axisymmetric force law".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "axisymmetric force fixture".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn recoil() -> IsotropicMagneticMaterial {
    let claims = [
        (
            RECOIL_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        ),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, 0.5),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            -100.0,
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            100.0,
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            100.0,
        ),
    ]
    .into_iter()
    .map(|(key, dims, value)| PropertyClaim {
        key: PropertyKey::new(key, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    })
    .collect();
    let resolved = ResolvedRecoilMaterial::resolve(
        &source(claims),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    IsotropicMagneticMaterial::Recoil(PlanarRecoilMaterial::new(resolved, [0.0, 1.0]).unwrap())
}
fn saturation() -> IsotropicMagneticMaterial {
    let card = source(vec![PropertyClaim {
        key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
        value: PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: FLUX_DENSITY_DIMS,
            knots: vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)],
            dims: FIELD_STRENGTH_DIMS,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "B",
            QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
            0.0,
            4.0,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }]);
    let point = QueryPoint::new()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap();
    IsotropicMagneticMaterial::Reversible(
        ReversibleMagneticMaterial::resolve(
            &card,
            &point,
            "B",
            4.0,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
}

#[test]
fn source_free_stress_is_zero_on_axis_sloping_and_toroidal_surfaces() {
    for reverse in [false, true] {
        let p = linear(8, reverse, core::f64::consts::TAU);
        let field = p
            .solve(
                &vec![0.0; p.mesh().faces().len()],
                &boundary(p.mesh(), 0.3, 0.2),
                config(),
                &mut || false,
            )
            .unwrap();
        // u=.3+.2z is an exact Maxwell quadrupole. Closed contours enclose
        // no source, even when their sides slope or their rotation is a torus.
        for mut region in [
            collar(p.mesh()),
            cells(p.mesh(), |[r, z]| r + 1.3 * (z - 0.1).abs() < 0.75),
            cells(p.mesh(), |[r, z]| r > 0.25 && r < 0.75 && z.abs() < 0.5),
        ] {
            let force = field
                .maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || false)
                .unwrap();
            near(axial(force), 0.0);
            region.reverse();
            assert_eq!(
                force,
                field
                    .maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || false)
                    .unwrap()
            );
        }
    }
}

#[test]
fn coil_stress_converges_to_independent_lorentz_and_fixed_current_work() {
    let beta = 0.2;
    let height = 0.5;
    let radial_moment = (0.5_f64.powi(3) - 0.25_f64.powi(3)) / 3.0;
    // Uniform Jphi=1, imposed Br=-beta*r. The centered coil's self-force
    // vanishes by symmetry. This integral uses no production force function.
    let exact = core::f64::consts::TAU * beta * height * radial_moment;
    // Independent coenergy interaction for a rigid axial shift at FIXED J:
    // integral Jphi*Aphi_external dV. The external quadrupole stays fixed.
    let interaction =
        |shift: f64| core::f64::consts::TAU * height * radial_moment * (0.3 + beta * shift);
    let delta = 1e-5;
    near(
        (interaction(delta) - interaction(-delta)) / (2.0 * delta),
        exact,
    );
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let p = linear(n, false, core::f64::consts::TAU);
        let conductor = coil(p.mesh());
        let mut current = vec![0.0; p.mesh().faces().len()];
        for &cell in &conductor {
            current[cell] = 1.0;
        }
        let field = p
            .solve(
                &current,
                &boundary(p.mesh(), 0.3, beta),
                config(),
                &mut || false,
            )
            .unwrap();
        let force = field
            .maxwell_wrench(
                &collar(p.mesh()),
                &frame(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let error = (axial(force) - exact).abs();
        assert!(
            error < 0.6 * previous,
            "n={n}: {error}, previous={previous}"
        );
        previous = error;
        let lorentz = field
            .impressed_current_wrench(&conductor, &frame(), [0.0; 3], force_limits(), &mut || {
                false
            })
            .unwrap();
        near(axial(lorentz), exact);
        let mut reordered = conductor;
        reordered.reverse();
        assert_eq!(
            lorentz,
            field
                .impressed_current_wrench(
                    &reordered,
                    &frame(),
                    [0.0; 3],
                    force_limits(),
                    &mut || false
                )
                .unwrap()
        );
    }
    assert!(previous < 0.04 * exact);
}

#[test]
fn coil_polarity_frame_origin_and_native_mechanics_keep_the_same_load() {
    let p = linear(8, false, core::f64::consts::TAU);
    let mut current = vec![0.0; p.mesh().faces().len()];
    for cell in coil(p.mesh()) {
        current[cell] = 1.0;
    }
    let region = collar(p.mesh());
    let mut reference = 0.0;
    for sign in [1.0, -1.0] {
        let source: Vec<_> = current.iter().map(|j| sign * j).collect();
        let field = p
            .solve(
                &source,
                &boundary(p.mesh(), 0.3, 0.2),
                config(),
                &mut || false,
            )
            .unwrap();
        let local = field
            .maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || false)
            .unwrap();
        if sign > 0.0 {
            reference = axial(local);
        } else {
            near(axial(local), -reference);
        }
        // r maps to global y, and the magnetic axis maps to global x.
        let placement =
            PlanarFrame::new([2.0, 3.0, 4.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        let origin = [2.2, 3.3, 4.4];
        let force = Vec3::new(local.force.z, 0.0, 0.0);
        let expected = Wrench::new(Vec3::new(-0.2, -0.3, -0.4).cross(force), force);
        let actual = field
            .maxwell_wrench(&region, &placement, origin, force_limits(), &mut || false)
            .unwrap();
        for (got, want) in actual.to_array().into_iter().zip(expected.to_array()) {
            near(got, want);
        }
        let inertia = SpatialInertia::new(
            2.0,
            Vec3::new(0.0, 0.0, 0.0),
            Mat3 {
                m: [0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2],
            },
        )
        .unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "coil",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia,
        )])
        .unwrap();
        let acceleration = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[],
            &[],
            &[],
            Vec3::new(0.0, 0.0, 0.0),
            &[local],
        )
        .unwrap();
        near(
            acceleration.base_spatial_acceleration_body.linear.z,
            local.force.z / 2.0,
        );
    }
}

#[test]
fn permanent_magnet_force_converges_without_an_impressed_current_surrogate() {
    // Recoil mu equals ambient mu=1, so fixed axial M=Br/mu=.5.
    // dBz_external/dz=.4. The centered magnet's self-force is zero.
    let volume = core::f64::consts::PI * 0.25_f64.powi(2) * 0.5;
    let expected = volume * 0.5 * 0.4;
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let mesh = mesh(n, false, core::f64::consts::TAU);
        let body = core(&mesh);
        let mut models = vec![IsotropicMagneticMaterial::Linear(material(1.0)); mesh.faces().len()];
        for &cell in &body {
            models[cell] = recoil();
        }
        let p = nonlinear(mesh, models);
        let initial: Vec<_> = p
            .mesh()
            .vertices()
            .iter()
            .map(|p| 0.3 + 0.2 * p[1])
            .collect();
        let field = p
            .solve(
                &vec![0.0; p.mesh().faces().len()],
                &boundary(p.mesh(), 0.3, 0.2),
                &initial,
                nonlinear_controls(),
                &mut || false,
            )
            .unwrap();
        let force = field
            .maxwell_wrench(
                &collar(p.mesh()),
                &frame(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let error = (axial(force) - expected).abs();
        assert!(
            error < 0.6 * previous,
            "n={n}: {error}, previous={previous}"
        );
        previous = error;
        let free_current = field
            .impressed_current_wrench(&body, &frame(), [0.0; 3], force_limits(), &mut || false)
            .unwrap();
        assert_eq!(free_current.to_array(), [0.0; 6]);
        assert!(field.recoil_operating_points().iter().any(Option::is_some));
    }
    assert!(previous < 0.04 * expected);
}

#[test]
fn saturated_core_has_an_exterior_load_but_its_nonlinear_surface_refuses() {
    let mesh = mesh(8, false, core::f64::consts::TAU);
    let body = core(&mesh);
    let mut models = vec![IsotropicMagneticMaterial::Linear(material(0.25)); mesh.faces().len()];
    let law = saturation();
    for &cell in &body {
        models[cell] = law.clone();
    }
    let p = nonlinear(mesh, models);
    let initial: Vec<_> = p
        .mesh()
        .vertices()
        .iter()
        .map(|p| 0.7 + 0.2 * p[1])
        .collect();
    let field = p
        .solve(
            &vec![0.0; p.mesh().faces().len()],
            &boundary(p.mesh(), 0.7, 0.2),
            &initial,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        body.iter()
            .any(|&cell| field.sample(cell, [1.0 / 3.0; 3]).unwrap().flux_density_t[1] > 1.0)
    );
    let force = field
        .maxwell_wrench(
            &collar(p.mesh()),
            &frame(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert!(axial(force) > 0.001);
    assert!(
        field
            .maxwell_wrench(&body, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
    assert_eq!(
        field
            .impressed_current_wrench(&body, &frame(), [0.0; 3], force_limits(), &mut || false)
            .unwrap()
            .to_array(),
        [0.0; 6]
    );
}

#[test]
fn contour_sources_full_revolution_limits_and_cancellation_are_required() {
    let p = linear(8, false, core::f64::consts::TAU);
    let mut current = vec![0.0; p.mesh().faces().len()];
    let conductor = coil(p.mesh());
    for &cell in &conductor {
        current[cell] = 1.0;
    }
    let region = collar(p.mesh());
    let field = p
        .solve(
            &current,
            &boundary(p.mesh(), 0.3, 0.2),
            config(),
            &mut || false,
        )
        .unwrap();
    for invalid in [
        vec![],
        vec![region[0], region[0]],
        vec![p.mesh().faces().len()],
    ] {
        assert!(
            field
                .maxwell_wrench(&invalid, &frame(), [0.0; 3], force_limits(), &mut || false)
                .is_err()
        );
        assert!(
            field
                .impressed_current_wrench(&invalid, &frame(), [0.0; 3], force_limits(), &mut || {
                    false
                })
                .is_err()
        );
    }
    assert!(
        field
            .maxwell_wrench(&conductor, &frame(), [0.0; 3], force_limits(), &mut || {
                false
            })
            .is_err()
    );
    let all: Vec<_> = (0..p.mesh().faces().len()).collect();
    assert!(
        field
            .maxwell_wrench(&all, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
    for insufficient in [
        MaxwellForceLimits {
            max_cells: p.mesh().faces().len() - 1,
            ..force_limits()
        },
        MaxwellForceLimits {
            max_edges: p.mesh().edges().len() - 1,
            ..force_limits()
        },
    ] {
        assert!(
            field
                .maxwell_wrench(&region, &frame(), [0.0; 3], insufficient, &mut || false)
                .is_err()
        );
    }
    assert!(
        field
            .maxwell_wrench(
                &region,
                &frame(),
                [f64::NAN, 0.0, 0.0],
                force_limits(),
                &mut || false
            )
            .is_err()
    );
    let before = field
        .maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        field.maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || {
            polls += 1;
            polls > region.len() + p.mesh().faces().len() + 7
        }),
        Err(MagneticError::Cancelled)
    ));
    assert_eq!(
        before,
        field
            .maxwell_wrench(&region, &frame(), [0.0; 3], force_limits(), &mut || false)
            .unwrap()
    );
    assert!(matches!(
        field
            .impressed_current_wrench(&conductor, &frame(), [0.0; 3], force_limits(), &mut || true),
        Err(MagneticError::Cancelled)
    ));
    let sector = linear(4, false, core::f64::consts::PI);
    let sf = sector
        .solve(
            &vec![0.0; sector.mesh().faces().len()],
            &boundary(sector.mesh(), 0.3, 0.0),
            config(),
            &mut || false,
        )
        .unwrap();
    assert!(
        sf.maxwell_wrench(
            &collar(sector.mesh()),
            &frame(),
            [0.0; 3],
            force_limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        sf.impressed_current_wrench(
            &coil(sector.mesh()),
            &frame(),
            [0.0; 3],
            force_limits(),
            &mut || false
        )
        .is_err()
    );

    // Matching permeability alone does not make two original material states
    // one homogeneous exterior source snapshot.
    let mesh = mesh(8, false, core::f64::consts::TAU);
    let alternate = source(vec![PropertyClaim {
        key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
        value: PropertyValue::Scalar {
            value: 1.0,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }]);
    let other = common::resolve(&alternate, &QueryPoint::new());
    let mut materials = vec![material(1.0); mesh.faces().len()];
    for cell in cells(&mesh, |[_, z]| z > 0.0) {
        materials[cell] = other.clone();
    }
    let p = AxisymmetricMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let field = p
        .solve(
            &vec![0.0; p.mesh().faces().len()],
            &boundary(p.mesh(), 0.3, 0.0),
            config(),
            &mut || false,
        )
        .unwrap();
    assert!(
        field
            .maxwell_wrench(
                &collar(p.mesh()),
                &frame(),
                [0.0; 3],
                force_limits(),
                &mut || false
            )
            .is_err()
    );
}
