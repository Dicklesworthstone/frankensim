//! Independent cylindrical saturation/recoil fields from synthetic source cards.
//! These fixtures establish numerical behavior, not measured material validity.
mod common;

use common::{config, limits, material};
use fs_em::planar::{
    MagneticError,
    axisymmetric::{
        AxisymmetricMagnetostatics,
        nonlinear::{
            AxisymmetricNonlinearAssembly, AxisymmetricNonlinearControls,
            AxisymmetricNonlinearError, AxisymmetricNonlinearMagnetostatics,
        },
    },
    nonlinear::IsotropicMagneticMaterial,
    permanent::{MagneticEnergyReference, PlanarRecoilMaterial},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::QuantitySpec;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};

fn source(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic cylindrical magnetics".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured axisymmetric magnetic law".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn saturation(maximum: f64) -> IsotropicMagneticMaterial {
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
            maximum,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
}
fn recoil(
    mu: f64,
    br: f64,
    parallel: [f64; 2],
    transverse: f64,
    direction: [f64; 2],
) -> IsotropicMagneticMaterial {
    let claims = [
        (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            parallel[0],
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            parallel[1],
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            transverse,
        ),
    ]
    .into_iter()
    .map(|(name, dims, value)| PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    })
    .collect();
    let law = ResolvedRecoilMaterial::resolve(
        &source(claims),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    IsotropicMagneticMaterial::Recoil(PlanarRecoilMaterial::new(law, direction).unwrap())
}
fn assembly(order: usize) -> AxisymmetricNonlinearAssembly {
    AxisymmetricNonlinearAssembly {
        limits: limits(),
        quadrature_order: order,
        max_quadrature_points: 100000,
    }
}
fn controls() -> AxisymmetricNonlinearControls {
    AxisymmetricNonlinearControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-11,
            linear_restart: 24,
            max_linear_cycles: 40,
            forcing_minimum: 1e-12,
            forcing_maximum: 0.1,
            forcing_gamma: 0.9,
            forcing_exponent: 1.5,
            globalization: Globalization::LineSearch(LineSearchConfig {
                armijo: 1e-4,
                contraction: 0.5,
                minimum_step: 1.0 / 1048576.0,
            }),
        },
        max_newton_iterations: 40,
        max_krylov_basis_entries: 100000,
    }
}
fn rectangle(r0: f64, radius: f64, height: f64, nr: usize, nz: usize, angle: f64) -> TriComplex2 {
    let mut points = Vec::new();
    let mut faces = Vec::new();
    for j in 0..=nz {
        for i in 0..=nr {
            let r = if i == nr {
                radius
            } else {
                r0 + (radius - r0) * i as f64 / nr as f64
            };
            let z = if j == nz {
                height
            } else {
                height * j as f64 / nz as f64
            };
            points.push([r, z]);
        }
    }
    for j in 0..nz {
        for i in 0..nr {
            let a = (j * (nr + 1) + i) as u32;
            let b = a + (nr + 1) as u32;
            faces.push([a, a + 1, b + 1]);
            faces.push([a, b + 1, b]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "nonlinear meridional fixture",
        points,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}
fn outer(mesh: &TriComplex2, radius: f64, value: f64) -> Vec<(usize, f64)> {
    mesh.vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == radius)
        .map(|(v, _)| (v, value))
        .collect()
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-8 * scale,
        "{actual} != {expected}; scale={scale}"
    );
}

#[test]
fn saturated_current_solution_has_exact_source_energy_and_legendre_work() {
    // u=1-r gives Bz=2-3r and Br=0. The source curve has H'=4 above
    // B=1, H'=1 below it. Ampere therefore gives Jphi=12 for r<1/3,
    // Jphi=3 outside. The mesh resolves that actual saturation interface.
    let height = 0.5;
    for angle in [core::f64::consts::TAU, core::f64::consts::TAU / 7.0] {
        let mesh = rectangle(0.0, 1.0, height, 6, 3, angle);
        let currents: Vec<_> = mesh
            .faces()
            .iter()
            .map(|face| {
                let r = face
                    .iter()
                    .map(|&v| mesh.vertices()[v as usize][0])
                    .sum::<f64>()
                    / 3.0;
                if r < 1.0 / 3.0 { 12.0 } else { 3.0 }
            })
            .collect();
        let p = AxisymmetricNonlinearMagnetostatics::new(
            mesh,
            vec![saturation(4.0); currents.len()],
            assembly(3),
            &mut || false,
        )
        .unwrap();
        let boundary = outer(p.mesh(), 1.0, 0.0);
        for sign in [1.0, -1.0] {
            let j: Vec<_> = currents.iter().map(|v| sign * v).collect();
            let f = p
                .solve(
                    &j,
                    &boundary,
                    &vec![0.0; p.mesh().vertices().len()],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for (&[r, _], &u) in p.mesh().vertices().iter().zip(f.potential_over_radius_t()) {
                near(u, sign * (1.0 - r), 1.0);
            }
            let scale = angle * height;
            near(f.energy_j(), scale * 5.0 / 36.0, scale);
            near(f.coenergy_j(), scale * 7.0 / 36.0, scale);
            near(f.impressed_current_pairing_j(), scale / 3.0, scale);
            near(f.boundary_pairing_j(), 0.0, scale);
            near(f.work_balance_error_j(), 0.0, scale);
            assert!((2.0 * f.energy_j() - f.impressed_current_pairing_j()).abs() > 0.05 * scale);
            assert!(f.report().converged && f.report().iterations > 0);
            assert_eq!(f.energy_reference(), MagneticEnergyReference::ZeroFlux);
        }
    }
}

#[test]
fn linear_limit_reuses_original_cylindrical_fields_in_both_components() {
    let mu = 0.5;
    let exact = |r: f64, z: f64| 1.0 - r + 0.3 + 0.2 * z;
    let mesh = rectangle(0.0, 1.0, 0.7, 6, 4, 1.2);
    let boundary: Vec<_> = mesh
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 1.0 || p[1] == 0.0 || p[1] == 0.7)
        .map(|(v, p)| (v, exact(p[0], p[1])))
        .collect();
    let mut initial = vec![0.0; mesh.vertices().len()];
    for &(v, u) in &boundary {
        initial[v] = u;
    }
    let currents = vec![3.0 / mu; mesh.faces().len()];
    let linear = AxisymmetricMagnetostatics::new(
        mesh.clone(),
        vec![material(mu); currents.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let nonlinear = AxisymmetricNonlinearMagnetostatics::new(
        mesh,
        vec![IsotropicMagneticMaterial::Linear(material(mu)); currents.len()],
        assembly(3),
        &mut || false,
    )
    .unwrap();
    let reference = linear
        .solve(&currents, &boundary, config(), &mut || false)
        .unwrap();
    let f = nonlinear
        .solve(&currents, &boundary, &initial, controls(), &mut || false)
        .unwrap();
    for (&u, &v) in f
        .potential_over_radius_t()
        .iter()
        .zip(reference.potential_over_radius_t())
    {
        near(u, v, 2.0);
    }
    for (cell, face) in nonlinear.mesh().faces().iter().enumerate() {
        let r = face
            .iter()
            .map(|&v| nonlinear.mesh().vertices()[v as usize][0])
            .sum::<f64>()
            / 3.0;
        let z = face
            .iter()
            .map(|&v| nonlinear.mesh().vertices()[v as usize][1])
            .sum::<f64>()
            / 3.0;
        let sample = f.sample(cell, [1.0 / 3.0; 3]).unwrap();
        near(sample.flux_density_t[0], -0.2 * r, 2.0);
        near(sample.flux_density_t[1], 2.6 - 3.0 * r + 0.4 * z, 3.0);
        near(
            sample.magnetic_field_a_m[1],
            sample.flux_density_t[1] / mu,
            6.0,
        );
        near(sample.azimuthal_potential_wb_m, r * exact(r, z), 2.0);
    }
    near(f.energy_j(), reference.energy_j(), reference.energy_j());
    near(f.work_balance_error_j(), 0.0, reference.energy_j());
}

#[test]
fn recoil_keeps_shifted_storage_negative_coenergy_and_original_margins() {
    let mu = 0.25;
    let br = 0.8;
    let mesh = rectangle(0.0, 1.0, 2.0, 3, 2, 1.2);
    let p = AxisymmetricNonlinearMagnetostatics::new(
        mesh,
        vec![recoil(mu, br, [-50.0, 50.0], 50.0, [0.0, 1.0]); 12],
        assembly(3),
        &mut || false,
    )
    .unwrap();
    let zero = vec![0.0; p.mesh().vertices().len()];
    let currents = vec![0.0; p.mesh().faces().len()];
    let f = p
        .solve(
            &currents,
            &outer(p.mesh(), 1.0, 0.0),
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    let volume = 1.2;
    let energy = volume * br * br / (2.0 * mu);
    near(f.energy_j(), energy, energy);
    near(f.coenergy_j(), -energy, energy);
    near(f.work_balance_error_j(), 0.0, energy);
    assert_eq!(
        f.energy_reference(),
        MagneticEnergyReference::FixedRemanence
    );
    for states in f.recoil_operating_points() {
        for state in states.unwrap() {
            near(state.parallel_field_a_m, -br / mu, br / mu);
            near(state.transverse_field_a_m, 0.0, br / mu);
            assert!(state.lower_parallel_margin_a_m > 0.0);
        }
    }
    let boundary = outer(p.mesh(), 1.0, br / 2.0);
    let mut initial = zero;
    for &(v, u) in &boundary {
        initial[v] = u;
    }
    let released = p
        .solve(&currents, &boundary, &initial, controls(), &mut || false)
        .unwrap();
    for &u in released.potential_over_radius_t() {
        near(u, br / 2.0, br);
    }
    near(released.energy_j(), 0.0, energy);
    assert!(released.report().iterations > 0);
}

#[test]
fn full_cell_source_limits_include_unsampled_vertices_and_axis_orientation() {
    let mesh = TriComplex2::from_indexed_triangles(
        "source bound",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::axisymmetric(1.0).unwrap(),
    )
    .unwrap();
    // u=r: Bz=3r reaches 3 at the vertex; every order-three interior
    // quadrature sample is below 2.6. Sampling alone would admit a bad field.
    let p = AxisymmetricNonlinearMagnetostatics::new(
        mesh.clone(),
        vec![saturation(2.6)],
        assembly(3),
        &mut || false,
    )
    .unwrap();
    let boundary = [(0, 0.0), (1, 1.0), (2, 0.0)];
    assert!(matches!(
        p.solve(&[0.0], &boundary, &[0.0, 1.0, 0.0], controls(), &mut || {
            false
        }),
        Err(AxisymmetricNonlinearError::FieldRange { .. })
    ));
    let p = AxisymmetricNonlinearMagnetostatics::new(
        mesh.clone(),
        vec![recoil(1.0, 0.5, [-1.0, 2.1], 1.0, [0.0, 1.0])],
        assembly(3),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        p.solve(&[0.0], &boundary, &[0.0, 1.0, 0.0], controls(), &mut || {
            false
        }),
        Err(AxisymmetricNonlinearError::RecoilRange { .. })
    ));
    assert!(
        AxisymmetricNonlinearMagnetostatics::new(
            mesh,
            vec![recoil(1.0, 0.5, [-2.0, 2.0], 2.0, [1.0, 0.0])],
            assembly(3),
            &mut || false
        )
        .is_err()
    );
    let annulus = rectangle(0.5, 1.0, 1.0, 1, 1, 1.0);
    assert!(
        AxisymmetricNonlinearMagnetostatics::new(
            annulus,
            vec![recoil(1.0, 0.5, [-2.0, 2.0], 2.0, [1.0, 0.0]); 2],
            assembly(3),
            &mut || false
        )
        .is_ok()
    );
}

#[test]
fn saturated_solenoid_core_converges_to_independent_ampere_energy() {
    let a = 0.02_f64;
    let b = 0.04_f64;
    let radius = 0.06_f64;
    let height = 0.08;
    let j = 100.0;
    let mu = 0.5;
    // Hcore=J(b-a)=2 and the explicit source has H=4B-3 above B=1.
    let core_b = 1.25;
    let core_energy = 0.875;
    let angle = core::f64::consts::TAU;
    let flux_integral =
        core_b * a * a / 2.0 + mu * j * (b * (b * b - a * a) / 2.0 - (b.powi(3) - a.powi(3)) / 3.0);
    let primitive = |r: f64| b * b * r * r / 2.0 - 2.0 * b * r.powi(3) / 3.0 + r.powi(4) / 4.0;
    let expected = angle
        * height
        * (core_energy * a * a / 2.0 + mu * j * j / 2.0 * (primitive(b) - primitive(a)));
    let mut previous = f64::INFINITY;
    for nr in [12, 24, 48] {
        let mesh = rectangle(0.0, radius, height, nr, 4, angle);
        let mut models = Vec::new();
        let mut currents = Vec::new();
        for face in mesh.faces() {
            let r = face
                .iter()
                .map(|&v| mesh.vertices()[v as usize][0])
                .sum::<f64>()
                / 3.0;
            models.push(if r < a {
                saturation(4.0)
            } else {
                IsotropicMagneticMaterial::Linear(material(mu))
            });
            currents.push(if r > a && r < b { j } else { 0.0 });
        }
        let p = AxisymmetricNonlinearMagnetostatics::new(mesh, models, assembly(3), &mut || false)
            .unwrap();
        let outer_u = flux_integral / radius.powi(2);
        // This millimetre-scale fixture has a much smaller current residual
        // [A m2] than the unit-domain manufactured fields above.
        let mut solve_controls = controls();
        solve_controls.solver.absolute_tolerance = 1e-14;
        let f = p
            .solve(
                &currents,
                &outer(p.mesh(), radius, outer_u),
                &vec![outer_u; p.mesh().vertices().len()],
                solve_controls,
                &mut || false,
            )
            .unwrap();
        let error = (f.energy_j() - expected).abs() / expected;
        assert!(
            error < 0.45 * previous,
            "nr={nr}: {error}, previous {previous}"
        );
        previous = error;
        near(f.work_balance_error_j(), 0.0, expected);
        assert!(f.report().converged);
    }
    assert!(previous < 0.002, "solenoid energy error={previous}");
}

#[test]
fn bounded_quadrature_solver_and_cancellation_refuse_without_a_field() {
    let mesh = rectangle(0.0, 1.0, 1.0, 3, 2, 1.0);
    let count = mesh.faces().len();
    let models = vec![saturation(4.0); count];
    assert!(
        AxisymmetricNonlinearMagnetostatics::new(
            mesh.clone(),
            models.clone(),
            assembly(2),
            &mut || false
        )
        .is_err()
    );
    let mut insufficient = assembly(3);
    insufficient.max_quadrature_points = count * 9 - 1;
    assert!(
        AxisymmetricNonlinearMagnetostatics::new(
            mesh.clone(),
            models.clone(),
            insufficient,
            &mut || false
        )
        .is_err()
    );
    let p =
        AxisymmetricNonlinearMagnetostatics::new(mesh, models, assembly(3), &mut || false).unwrap();
    let boundary = outer(p.mesh(), 1.0, 0.0);
    let initial = vec![0.0; p.mesh().vertices().len()];
    let currents = vec![1.0; count];
    let mut exhausted = controls();
    exhausted.max_newton_iterations = 0;
    exhausted.max_krylov_basis_entries = 0;
    assert!(matches!(
        p.solve(&currents, &boundary, &initial, exhausted, &mut || false),
        Err(AxisymmetricNonlinearError::Nonconvergence(_))
    ));
    let mut insufficient = controls();
    insufficient.max_krylov_basis_entries = 1;
    assert!(
        p.solve(&currents, &boundary, &initial, insufficient, &mut || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        p.solve(&currents, &boundary, &initial, controls(), &mut || {
            polls += 1;
            polls > 30
        }),
        Err(AxisymmetricNonlinearError::Magnetic(
            MagneticError::Cancelled
        ))
    ));
    assert!(
        p.solve(
            &[f64::NAN; 12],
            &boundary,
            &initial,
            controls(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        p.solve(
            &currents,
            &[(0, 0.0), (4, 0.0), (8, 0.0)],
            &initial,
            controls(),
            &mut || false
        )
        .is_err()
    );
    let field = p
        .solve(
            &vec![0.0; count],
            &boundary,
            &initial,
            exhausted,
            &mut || false,
        )
        .unwrap();
    assert!(field.sample(0, [0.5, 0.5, 0.5]).is_err());
    near(field.energy_j(), 0.0, 1.0);
}
