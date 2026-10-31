//! Tensor source, orientation and energy oracles on actual native 3D fields.
use super::*;
use fs_em::spatial::{
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        induction::nonlinear::{
            SpatialNonlinearInduction, SpatialNonlinearInductionControls,
        },
        nonlinear::{NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics},
    },
};
use fs_material::conductor::{
    ElectricalPropertyBasis,
    spatial::{ResolvedOrthotropicConductor3, ResolvedSpatialConductor, SpatialConductorMaterial},
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
fn rotated(sign: f64) -> [[f64; 3]; 3] {
    let c = 0.5_f64.sqrt();
    [[c, sign * c, 0.0], [-sign * c, c, 0.0], [0.0, 0.0, 1.0]]
}
fn tensor_source(
    basis: ElectricalPropertyBasis,
    axes: [[f64; 3]; 3],
) -> (SpatialMagneticMaterial, ResolvedOrthotropicConductor3) {
    let names = ["principal-electrical-1", "principal-electrical-2", "principal-electrical-3"];
    let mut values = vec![
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 2.0),
        (ELECTRICAL_CONDUCTIVITY_PROPERTY, ELECTRICAL_CONDUCTIVITY_DIMS, 5.0),
    ];
    for (name, sigma) in names.into_iter().zip([2.0_f64, 5.0, 11.0]) {
        let value = match basis {
            ElectricalPropertyBasis::Conductivity => sigma,
            ElectricalPropertyBasis::Resistivity => sigma.recip(),
        };
        values.push((name, basis.dimensions(), value));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = vec![];
    for (name, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured spatial orthotropic conductor".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).unwrap());
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic tensor".into(), phase: "solid".into(),
            process: "fixed axes".into(), revision: 0,
        },
        claims, vec![],
    ).unwrap();
    let material = resolve_material_state_point(
        &card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let conductor = ResolvedOrthotropicConductor3::try_from_material(
        material.clone(), basis, names.map(|name| PropertyKey::new(name, basis.dimensions())), axes,
    ).unwrap();
    (SpatialMagneticMaterial::Linear(material), conductor)
}
fn tensor_cube(axes: [[f64; 3]; 3]) -> (SpatialVectorMagnetostatics, ResolvedOrthotropicConductor3) {
    let (material, conductor) = tensor_source(ElectricalPropertyBasis::Conductivity, axes);
    let (positions, tets) = grid(1);
    let problem = SpatialVectorMagnetostatics::new(
        positions, tets, vec![material; 6], assembly(), &mut || false,
    ).unwrap();
    (problem, conductor)
}

#[test]
fn exact_principal_sources_resistivity_and_fixed_frame_are_retained() {
    let (_, sigma) = tensor_source(ElectricalPropertyBasis::Conductivity, rotated(1.0));
    let (_, rho) = tensor_source(ElectricalPropertyBasis::Resistivity, rotated(1.0));
    // Independent rotation: sigma_xx=sigma_yy=3.5, sigma_xy=-1.5, sigma_zz=11.
    for conductor in [&sigma, &rho] {
        let expected = [[3.5, -1.5, 0.0], [-1.5, 3.5, 0.0], [0.0, 0.0, 11.0]];
        let actual = conductor.conductivity().matrix_s_m();
        for i in 0..3 {
            for j in 0..3 { near(actual[i][j], expected[i][j], 2e-14); }
        }
        assert_eq!(conductor.principal_axes(), Some(rotated(1.0)));
        assert_eq!(conductor.with_material(conductor.material().clone()).unwrap(), *conductor);
        assert!(conductor.has_admitted_static_sources());
    }
    assert!(!sigma.has_same_electrical_source(&rho));
    let keys: [PropertyKey; 3] = sigma.source_keys().cloned().collect::<Vec<_>>().try_into().unwrap();
    for axes in [
        [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ] {
        assert!(ResolvedOrthotropicConductor3::try_from_material(
            sigma.material().clone(), sigma.source_basis(), keys.clone(), axes,
        ).is_err());
    }
    let mut absent = keys.clone();
    absent[2] = PropertyKey::new("missing-third-source", ELECTRICAL_CONDUCTIVITY_DIMS);
    assert!(ResolvedOrthotropicConductor3::try_from_material(
        sigma.material().clone(), sigma.source_basis(), absent, IDENTITY,
    ).is_err());
    let mut wrong_units = keys;
    wrong_units[1] = PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS);
    assert!(ResolvedOrthotropicConductor3::try_from_material(
        sigma.material().clone(), sigma.source_basis(), wrong_units, IDENTITY,
    ).is_err());
}

#[test]
fn rotated_native_modes_match_independent_tensor_decay_and_lorentz_work() {
    // Exact P1 integration and charge elimination on the six cube tetrahedra:
    // K=4/mu=2, M=tr(sigma)/40-(sigma_xy+sigma_xz+sigma_yz)/120.
    // The minimizing nodal gradient is independent of a constant SPD tensor.
    for (axes, m) in [(IDENTITY, 0.45), (rotated(1.0), 0.4625), (rotated(-1.0), 0.4375)] {
        let (problem, conductor) = tensor_cube(axes);
        let old = [[4.0; 3]; 6];
        let initial = problem.solve(
            &old, &trace(&problem, [0.0; 3], false), static_controls(), &mut || false,
        ).unwrap();
        let state = SpatialInduction::new_with_conductors(
            &initial, vec![Some(conductor.clone()); 6], preparation(), 2e-9, &mut || false,
        ).unwrap();
        for alpha in [0.0, 0.6] {
            let (dt, k): (f64, f64) = (0.03, 2.0);
            let velocity: Vec<_> = problem.positions_m().iter().map(|p| p.map(|x| alpha * x)).collect();
            let endpoint = if alpha == 0.0 { [[0.0; 3]; 6] } else { old };
            let supply = if alpha == 0.0 { 0.5 * k } else { k };
            let effective = k + 2.0 * alpha * m;
            let a1 = (m - 0.5 * dt * effective + dt * supply) / (m + 0.5 * dt * effective);
            let amid = 0.5 * (1.0 + a1);
            let rate = (a1 - 1.0) / dt;
            let conductor_rate = rate + 2.0 * alpha * amid;
            let candidate = state.step_candidate(
                &endpoint, &velocity, dt, controls(), &mut || false,
            ).unwrap();
            let step = candidate.step();
            near(candidate.state().potential_wb()[diagonal(&problem)], a1, 3e-9);
            near(step.heat_j, dt * m * conductor_rate.powi(2), 3e-9);
            near(step.mechanical_work_j, -2.0 * alpha * amid * m * conductor_rate * dt, 3e-9);
            near(step.impressed_work_j, supply * (a1 - 1.0), 3e-9);
            near(step.stored_after_j - step.stored_before_j, 0.5 * k * (a1 * a1 - 1.0), 3e-9);
            near(step.energy_residual_j, 0.0, 3e-9);
            assert!(step.heat_j > 0.0);
            assert!(step.maximum_ampere_residual_a < 2e-8 && step.maximum_charge_residual_a < 2e-8);
            let sigma = conductor.conductivity().matrix_s_m();
            for cell in &step.cells {
                for q in 0..4 {
                    let e = cell.conductor_electric_field_v_m[q];
                    let j = std::array::from_fn(|i| (0..3).map(|k| sigma[i][k] * e[k]).sum::<f64>());
                    let force = cross(j, step.stage_flux_density_t[cell.cell]);
                    for i in 0..3 {
                        near(cell.current_density_a_m2[q][i], j[i], 3e-9);
                        near(cell.lorentz_force_density_n_m3[q][i], force[i], 3e-9);
                    }
                }
            }
            assert_eq!(state.time_s(), 0.0);
            assert_eq!(candidate.state().conductors()[0].as_ref().unwrap().principal_axes(), Some(axes));
        }
    }
}

#[test]
fn original_scalar_constructor_and_mixed_source_owner_agree_exactly() {
    let (problem, tensor) = tensor_cube(rotated(1.0));
    let scalar = ResolvedConductorState::try_from_material(tensor.material().clone()).unwrap();
    let initial = problem.solve(
        &[[4.0; 3]; 6], &trace(&problem, [0.0; 3], false), static_controls(), &mut || false,
    ).unwrap();
    let scalar_state = SpatialInduction::new(
        &initial, vec![Some(scalar.clone()); 6], preparation(), 2e-9, &mut || false,
    ).unwrap();
    let mixed_state = SpatialInduction::new_with_conductors(
        &initial, vec![Some(ResolvedSpatialConductor::from(scalar)); 6],
        preparation(), 2e-9, &mut || false,
    ).unwrap();
    let a = scalar_state.step_candidate(&[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.04, controls(), &mut || false).unwrap();
    let b = mixed_state.step_candidate(&[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.04, controls(), &mut || false).unwrap();
    assert_eq!(a.state().potential_wb(), b.state().potential_wb());
    assert_eq!(a.step().heat_j, b.step().heat_j);
    let mut heterogeneous = mixed_state.conductors().to_vec();
    heterogeneous[0] = Some(tensor.into());
    let mixed = SpatialInduction::new_with_conductors(
        &initial, heterogeneous, preparation(), 2e-9, &mut || false,
    ).unwrap().step_candidate(
        &[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.04, controls(), &mut || false,
    ).unwrap();
    assert!(mixed.step().heat_j > 0.0);
    near(mixed.step().energy_residual_j, 0.0, 3e-9);
}

#[test]
fn nonlinear_owner_retains_rotated_tensor_and_separates_numerical_loss() {
    let (material, conductor) = tensor_source(ElectricalPropertyBasis::Conductivity, rotated(1.0));
    let SpatialMagneticMaterial::Linear(material) = material else { unreachable!() };
    let (positions, tets) = grid(1);
    let problem = SpatialNonlinearVectorMagnetostatics::new(
        positions, tets, vec![SpatialNonlinearMaterial::Linear(material); 6],
        assembly(), &mut || false,
    ).unwrap();
    let boundary = VectorMagneticBoundary {
        potential: problem.exterior_faces().into_iter().map(|face| VectorPotentialFace {
            face, edge_integrals_wb: [0.0; 3],
        }).collect(),
    };
    let initial = problem.solve(
        &[[4.0; 3]; 6], &boundary, &vec![0.0; problem.mesh().edges.len()],
        NonlinearVectorControls {
            field: static_controls(), relative_residual_tolerance: 2e-10,
            maximum_newton_iterations: 40, maximum_backtracks: 35,
        }, &mut || false,
    ).unwrap();
    let state = SpatialNonlinearInduction::new_with_conductors(
        &initial, vec![Some(conductor); 6], preparation(), 2e-9, &mut || false,
    ).unwrap();
    let (m, k, dt): (f64, f64, f64) = (0.4625, 2.0, 0.04);
    let a1 = m / (m + k * dt);
    let candidate = state.step_candidate(
        &[[0.0; 3]; 6], &[[0.0; 3]; 8], dt,
        SpatialNonlinearInductionControls {
            induction: controls(), maximum_newton_iterations: 40, maximum_backtracks: 35,
        }, &mut || false,
    ).unwrap();
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], a1, 3e-9);
    near(candidate.step().heat_j, m * (a1 - 1.0).powi(2) / dt, 3e-9);
    near(candidate.step().numerical_dissipation_j, 0.5 * k * (a1 - 1.0).powi(2), 3e-9);
    near(candidate.step().energy_residual_j, 0.0, 3e-9);
    assert_eq!(candidate.state().conductors()[0].as_ref().unwrap().principal_axes(), Some(rotated(1.0)));
}
