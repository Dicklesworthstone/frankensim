//! Analytic G0/G1 fields from manufactured recoil claims. These source-declared
//! test models are NOT measured material or seed-catalog validation.
mod common;
use common::*;
use fs_em::{
    force::MaxwellForceLimits,
    planar::{
        MagneticError, PlanarMagnetostatics,
        eddy::{EddyControls, EddyError, PlanarEddy},
        permanent::{MagneticEnergyReference, PlanarMagneticMaterial, PlanarRecoilMaterial},
    },
    transfer::{CellSample, PlanarFrame},
    winding::{PlanarWindingBasis, WindingConfig},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};

fn claims(mu: f64, br: f64, parallel: [f64; 2], transverse: f64) -> Vec<PropertyClaim> {
    [
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
        provenance: Provenance {
            source: "manufactured reversible recoil law".into(),
            license: "test-only".into(),
            artifact: None,
        },
    })
    .collect()
}
fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured-magnet".into(),
            phase: "solid".into(),
            process: "recoil-regression".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn recoil(mu: f64, br: f64) -> ResolvedRecoilMaterial {
    ResolvedRecoilMaterial::resolve(
        &card(claims(mu, br, [-1e6, 1e6], 1e6)),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn magnet(law: ResolvedRecoilMaterial, direction: [f64; 2]) -> PlanarMagneticMaterial {
    PlanarMagneticMaterial::Recoil(PlanarRecoilMaterial::new(law, direction).unwrap())
}
fn triangle() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "natural-recoil-source",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap()
}
fn layered_problem(law: ResolvedRecoilMaterial, direction: [f64; 2]) -> PlanarMagnetostatics {
    PlanarMagnetostatics::new_with_material_models(
        layered(1.0, 2.0, 0.3),
        vec![
            magnet(law.clone(), direction),
            magnet(law, direction),
            PlanarMagneticMaterial::Linear(material(1e-6)),
            PlanarMagneticMaterial::Linear(material(1e-6)),
        ],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn enclosed_recoil_magnets_drive_static_wrenches_but_remanence_on_the_contour_refuses() {
    let n = 16;
    let depth = 0.2;
    let mu = 1e-6;
    let br = 0.01;
    for (direction, magnetized_exterior) in [
        ([1.0, 0.0], false),
        ([-1.0, 0.0], false),
        ([1.0, 0.0], true),
    ] {
        let mut vertices = Vec::new();
        for y in 0..=n {
            for x in 0..=n {
                vertices.push([
                    -1.0 + 2.0 * x as f64 / n as f64,
                    -1.0 + 2.0 * y as f64 / n as f64,
                ]);
            }
        }
        let mut faces = Vec::new();
        for y in 0..n {
            for x in 0..n {
                let a = (y * (n + 1) + x) as u32;
                let b = a + 1;
                let d = a + (n + 1) as u32;
                faces.extend([[a, b, d + 1], [a, d + 1, d]]);
            }
        }
        let ambient = PlanarMagneticMaterial::Linear(material(mu));
        let permanent = magnet(recoil(mu, br), direction);
        let mut models = Vec::new();
        let mut enclosed = Vec::new();
        for (cell, face) in faces.iter().enumerate() {
            let center: [f64; 2] = [0, 1].map(|axis| {
                face.iter()
                    .map(|&v| vertices[v as usize][axis])
                    .sum::<f64>()
                    / 3.0
            });
            let radius = center[0].abs().max(center[1].abs());
            if radius < 0.5 {
                enclosed.push(cell);
            }
            models.push(if magnetized_exterior || radius < 0.25 {
                permanent.clone()
            } else {
                ambient.clone()
            });
        }
        let mesh = TriComplex2::from_indexed_triangles(
            "recoil-exterior-force",
            vertices,
            faces,
            Metric2::planar(depth).unwrap(),
        )
        .unwrap();
        let boundary: Vec<_> = mesh
            .boundary_trace()
            .unwrap()
            .vertices
            .into_iter()
            .map(|v| (v as usize, 1e-3 * mesh.vertices()[v as usize][0]))
            .collect();
        let problem =
            PlanarMagnetostatics::new_with_material_models(mesh, models, limits(), &mut || false)
                .unwrap();
        let currents = vec![0.0; problem.mesh().faces().len()];
        let field = problem
            .solve(&currents, &boundary, config(), &mut || false)
            .unwrap();
        assert!(field.problem().has_permanent_magnets());
        assert_eq!(
            field.impressed_current_density_a_m2(),
            Some(currents.as_slice())
        );
        let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let load = field.maxwell_wrench(
            &enclosed,
            &frame,
            [0.0, 0.0, depth / 2.0],
            MaxwellForceLimits {
                max_cells: 10000,
                max_edges: 20000,
            },
            &mut || false,
        );
        if magnetized_exterior {
            assert!(matches!(
                load,
                Err(MagneticError::Invalid(
                    "Maxwell force contour requires an induced-linear exterior without remanence"
                ))
            ));
        } else {
            let load = load.unwrap();
            // With equal recoil/ambient slopes, M=Br/mu. The uniform applied
            // B=(0,-0.001,0) gives torque_z=(M cross B)_z * magnet volume.
            // The original P1 exterior stress remains a finite-mesh estimate.
            let exact_torque = -direction[0] * (br / mu) * 0.25 * depth * 1e-3;
            assert!(
                (load.torque.z - exact_torque).abs() < 0.2 * exact_torque.abs(),
                "load={load:?}, exact torque={exact_torque}"
            );
            assert!(load.force.x.abs() < 1e-8 && load.force.y.abs() < 1e-8);
            let profile: Vec<_> = problem
                .mesh()
                .faces()
                .iter()
                .enumerate()
                .map(|(cell, face)| {
                    if problem.recoil_materials()[cell].is_some() {
                        let x_sum: f64 = face
                            .iter()
                            .map(|&v| problem.mesh().vertices()[v as usize][0])
                            .sum();
                        1000.0 * x_sum.signum()
                    } else {
                        0.0
                    }
                })
                .collect();
            let prescribed: Vec<_> = boundary.iter().map(|&(v, _)| v).collect();
            let winding = PlanarWindingBasis::new_with_fixed_remanence(
                &problem,
                vec![profile.clone()],
                &prescribed,
                vec![0.0],
                WindingConfig {
                    max_windings: 1,
                    max_basis_values: 20000,
                    pairing_tolerance_h: 1e-8,
                    current_balance_tolerance_turns: 1e-8,
                    inverse_residual_tolerance: 1e-8,
                },
                1e-8,
                config(),
                &mut || false,
            )
            .unwrap();
            let actual = winding.field_for_currents(&[0.7], &mut || false).unwrap();
            let currents: Vec<_> = profile.iter().map(|j| 0.7 * j).collect();
            assert!(!actual.is_incremental_response());
            assert_eq!(
                actual.impressed_current_density_a_m2(),
                Some(currents.as_slice())
            );
            let homogeneous: Vec<_> = prescribed.iter().map(|&v| (v, 0.0)).collect();
            let direct = problem
                .solve(&currents, &homogeneous, config(), &mut || false)
                .unwrap();
            let load = |field: &fs_em::planar::PlanarField<'_>| {
                field
                    .maxwell_wrench(
                        &enclosed,
                        &frame,
                        [0.0, 0.0, depth / 2.0],
                        MaxwellForceLimits {
                            max_cells: 10000,
                            max_edges: 20000,
                        },
                        &mut || false,
                    )
                    .unwrap()
            };
            for (got, expected) in load(&actual)
                .to_array()
                .into_iter()
                .zip(load(&direct).to_array())
            {
                near(got, expected);
            }
        }
    }
}

#[test]
fn constant_remanence_needs_its_natural_boundary_load_and_reuses_field_transfer() {
    let law = recoil(1.5e-6, 0.9);
    let problem = PlanarMagnetostatics::new_with_material_models(
        triangle(),
        vec![magnet(law.clone(), [1.0, 0.0])],
        limits(),
        &mut || false,
    )
    .unwrap();
    // Only the bottom edge has prescribed A=0. On the other edges H_t=0.
    // Exact solution A=Br*y: B=Br*x_hat, H=0. Curl(Br)=0 inside the cell,
    // so using only an interior curl source would wrongly return A=B=0.
    let field = problem
        .solve(&[0.0], &[(0, 0.0), (1, 0.0)], config(), &mut || false)
        .unwrap();
    near(field.potential_wb_m()[2], 0.9);
    near(field.magnetic_flux_density_t()[0][0], 0.9);
    near(field.magnetic_flux_density_t()[0][1], 0.0);
    near(field.magnetic_field_a_m()[0][0], 0.0);
    near(field.magnetic_field_a_m()[0][1], 0.0);
    near(field.energy_j(), 0.0);
    assert!(field.report().unwrap().converged_euclidean());
    assert_eq!(
        field.energy_reference(),
        MagneticEnergyReference::FixedRemanence
    );
    assert_eq!(&problem.materials()[0], law.material());
    assert_eq!(
        problem.recoil_materials()[0].as_ref().unwrap().material(),
        &law
    );
    let operating = field.recoil_operating_points()[0].unwrap();
    near(operating.parallel_field_a_m, 0.0);
    near(operating.lower_parallel_margin_a_m, 1e6);
    let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    let samples = field
        .sample_cells(
            &frame,
            &[CellSample {
                triangle: 0,
                position_m: [0.25, 0.25, 0.1],
            }],
            1e-12,
            1,
            &mut || false,
        )
        .unwrap();
    for (got, expected) in samples[0].into_iter().zip([0.9, 0.0, 0.0]) {
        near(got, expected);
    }
}

#[test]
fn impressed_current_and_remanence_add_without_confusing_physical_h_or_energy() {
    let mu = 2e-6;
    let br = 0.8;
    let j = 3e5;
    for sign in [-1.0, 1.0] {
        let problem = PlanarMagnetostatics::new_with_material_models(
            triangle(),
            vec![magnet(recoil(mu, br), [sign, 0.0])],
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = problem
            .solve(&[j], &[(0, 0.0), (1, 0.0)], config(), &mut || false)
            .unwrap();
        // K_22=volume/mu; b_2=volume*(J/3+signed_Br/mu).
        near(
            field.magnetic_flux_density_t()[0][0],
            sign * br + mu * j / 3.0,
        );
        near(field.magnetic_field_a_m()[0][0], j / 3.0);
        near(
            field.recoil_operating_points()[0]
                .unwrap()
                .parallel_field_a_m,
            sign * j / 3.0,
        );
        near(field.energy_j(), 0.5 * 0.1 * mu * (j / 3.0) * (j / 3.0));
    }
}

#[test]
fn material_swaps_and_orientation_match_the_independent_layered_magnet_air_solution() {
    for (mu, br) in [(2e-6, 0.8), (4e-6, 1.2)] {
        for sign in [-1.0, 1.0] {
            let problem = layered_problem(recoil(mu, br), [0.0, sign]);
            let bc = ends(0.0);
            let field = problem
                .solve(&[0.0; 4], &bc, config(), &mut || false)
                .unwrap();
            // With equal layer widths and equal A on the two end boundaries,
            // integral B_y dx=0. Interface tangential H continuity then gives
            // H_y=-signed_Br/(mu_magnet+mu_air), independent of FE assembly.
            let h = -sign * br / (mu + 1e-6);
            for cell in 0..4 {
                near(field.magnetic_field_a_m()[cell][0], 0.0);
                near(field.magnetic_field_a_m()[cell][1], h);
                near(field.magnetic_flux_density_t()[cell][0], 0.0);
                near(
                    field.magnetic_flux_density_t()[cell][1],
                    if cell < 2 {
                        mu * h + sign * br
                    } else {
                        1e-6 * h
                    },
                );
                assert_eq!(field.recoil_operating_points()[cell].is_some(), cell < 2);
            }
            near(field.energy_j(), 0.5 * 0.6 * (mu + 1e-6) * h * h);
            let shifted: Vec<_> = bc.iter().map(|&(v, a)| (v, a + 7.0)).collect();
            let gauge = problem
                .solve(&[0.0; 4], &shifted, config(), &mut || false)
                .unwrap();
            near(gauge.energy_j(), field.energy_j());
            for (a, b) in gauge
                .magnetic_flux_density_t()
                .iter()
                .zip(field.magnetic_flux_density_t())
            {
                near(a[0], b[0]);
                near(a[1], b[1]);
            }
        }
    }
}

#[test]
fn source_pinned_temperature_refresh_updates_recoil_slope_remanence_and_limits_together() {
    let mut temperature_claims = claims(2e-6, 0.9, [-4e5, 5e5], 1e4);
    for (claim, hot) in temperature_claims
        .iter_mut()
        .zip([4e-6, 0.6, -2e5, 3e5, 1e4])
    {
        let PropertyValue::Scalar { value, dims } = claim.value.clone() else {
            unreachable!()
        };
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, value), (400.0, hot)],
            dims,
        };
        claim.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        claim.interpolation = InterpolationPolicy::LinearInside;
    }
    let card = card(temperature_claims);
    let cold = ResolvedRecoilMaterial::resolve(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let hot = cold
        .at_query_point(&card, &QueryPoint::new().with("T", 400.0).unwrap())
        .unwrap();
    assert_eq!(
        cold.material().card_identity(),
        hot.material().card_identity()
    );
    assert_ne!(cold.material().identity(), hot.material().identity());
    for (a, b) in cold
        .material()
        .properties()
        .iter()
        .zip(hot.material().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    for (law, expected_b, expected_margin) in [(cold.clone(), 0.3, 1e5), (hot, 0.12, 8e4)] {
        let problem = layered_problem(law, [0.0, 1.0]);
        let field = problem
            .solve(&[0.0; 4], &ends(0.0), config(), &mut || false)
            .unwrap();
        near(field.magnetic_flux_density_t()[0][1], expected_b);
        near(
            field.recoil_operating_points()[0]
                .unwrap()
                .lower_parallel_margin_a_m,
            expected_margin,
        );
    }
    assert!(
        cold.at_query_point(&card, &QueryPoint::new().with("T", 450.0).unwrap())
            .is_err()
    );
    let replacement = self::card(claims(2e-6, 0.7, [-4e5, 5e5], 1e4));
    assert!(
        cold.at_query_point(&replacement, &QueryPoint::new().with("T", 300.0).unwrap())
            .is_err()
    );
}

#[test]
fn computed_parallel_and_transverse_recoil_exits_refuse_before_publishing_a_field() {
    let law = ResolvedRecoilMaterial::resolve(
        &card(claims(2e-6, 0.8, [-1e5, 2e5], 1e4)),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let problem = PlanarMagnetostatics::new_with_material_models(
        triangle(),
        vec![magnet(law, [1.0, 0.0])],
        limits(),
        &mut || false,
    )
    .unwrap();
    for (at_x, at_y) in [(0.0, 0.4), (0.0, 1.4), (-0.04, 0.8)] {
        // All DOFs are prescribed: bounds must also be checked when there is
        // no algebraic solve or residual report to drive downstream admission.
        assert!(matches!(
            problem.solve(
                &[0.0],
                &[(0, 0.0), (1, at_x), (2, at_y)],
                config(),
                &mut || false
            ),
            Err(MagneticError::Material(_))
        ));
    }
    let valid = problem
        .solve(
            &[0.0],
            &[(0, 0.0), (1, 0.0), (2, 0.9)],
            config(),
            &mut || false,
        )
        .unwrap();
    near(
        valid.recoil_operating_points()[0]
            .unwrap()
            .parallel_field_a_m,
        5e4,
    );
    assert!(valid.report().is_none());
}

#[test]
fn recoil_requires_original_units_explicit_source_roles_and_fixed_condition_domains() {
    let base = claims(2e-6, 0.8, [-1e6, 1e6], 1e6);
    let admit = |claims, point| {
        ResolvedRecoilMaterial::resolve(
            &card(claims),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    let mut missing = base.clone();
    missing.remove(0);
    assert!(admit(missing, QueryPoint::new()).is_err());
    let mut generic = base.clone();
    generic[0].key = PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS);
    assert!(admit(generic, QueryPoint::new()).is_err());
    let mut wrong_units = base.clone();
    wrong_units[1].key = PropertyKey::new(REMANENT_FLUX_DENSITY_PROPERTY, Dims([0; 6]));
    wrong_units[1].value = PropertyValue::Scalar {
        value: 0.8,
        dims: Dims([0; 6]),
    };
    assert!(admit(wrong_units, QueryPoint::new()).is_err());
    let mut field_frozen = base.clone();
    field_frozen[0].value = PropertyValue::Curve {
        abscissa: "B".into(),
        abscissa_dims: FLUX_DENSITY_DIMS,
        knots: vec![(0.0, 2e-6), (1.0, 2e-6)],
        dims: MAGNETIC_PERMEABILITY_DIMS,
    };
    field_frozen[0].interpolation = InterpolationPolicy::LinearInside;
    field_frozen[0].validity = ValidityDomain::unconstrained().with("B", 0.0, 1.0);
    assert!(admit(field_frozen, QueryPoint::new().with("B", 0.5).unwrap()).is_err());
    for (br, parallel, transverse) in [
        (0.0, [-1e6, 1e6], 1e6),
        (0.8, [1.0, 10.0], 1e6),
        (0.8, [-1e6, 1e6], -1.0),
        (0.8, [1e6, -1e6], 1e6),
    ] {
        assert!(admit(claims(2e-6, br, parallel, transverse), QueryPoint::new()).is_err());
    }
    for direction in [[0.0, 0.0], [2.0, 0.0], [f64::NAN, 0.0]] {
        assert!(PlanarRecoilMaterial::new(recoil(2e-6, 0.8), direction).is_err());
    }
}

#[test]
fn affine_windings_admit_fixed_magnets_while_dynamic_induction_still_refuses() {
    let problem = layered_problem(recoil(2e-6, 0.8), [0.0, 1.0]);
    let field = problem
        .solve(&[0.0; 4], &ends(0.0), config(), &mut || false)
        .unwrap();
    assert!(matches!(
        PlanarEddy::new(
            &field,
            vec![],
            &[0, 1, 2, 3],
            &ends(0.0),
            0.0,
            EddyControls {
                linear: config(),
                energy_tolerance_j: 1e-8,
                current_balance_tolerance_a: 1e-8
            },
            &mut || false
        ),
        Err(EddyError::Invalid(
            "permanent-magnet induction requires explicit remanence source and parameter-work ownership"
        ))
    ));
    let winding = PlanarWindingBasis::new_with_fixed_remanence(
        &problem,
        vec![vec![1.0, 1.0, -1.0, -1.0]],
        &[0, 3, 2, 5],
        vec![0.0],
        WindingConfig {
            max_windings: 1,
            max_basis_values: 1000,
            pairing_tolerance_h: 1e-8,
            current_balance_tolerance_turns: 1e-8,
            inverse_residual_tolerance: 1e-8,
        },
        1e-8,
        config(),
        &mut || false,
    )
    .unwrap();
    assert!(winding.bias_field().is_some());
    assert!(winding.unit_fields()[0].is_incremental_response());
    let reconstructed = winding.field_for_currents(&[0.0], &mut || false).unwrap();
    near(reconstructed.energy_j(), field.energy_j());
    assert!(matches!(
        problem.solve(&[0.0; 4], &ends(0.0), config(), &mut || true),
        Err(MagneticError::Cancelled)
    ));
}
