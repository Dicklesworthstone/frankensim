//! Manufactured physical interfaces, not measured material data.
//! Independent two-body balances check the native full-tensor source owner.

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConductivityTensor, ContactResistance, EnergyConfig, FaceContactResistance,
    SolidMaterial, SourceFaceContact, TENSOR_SDIRK2_GAMMA, TensorConductionSolution,
    TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_nonlinear_tensor_conduction, advance_nonlinear_tensor_conduction_sdirk2,
    solve_nonlinear_tensor_conduction, solve_tensor_conduction,
};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

#[derive(Clone, Copy)]
enum Response {
    Constant,
    Curve,
    Surface,
}

fn resistance([a, b]: [f64; 2]) -> f64 {
    0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0
}

fn kelvin() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn pressure() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn separation() -> QuantitySpec {
    QuantitySpec::dimensional(fs_qty::Length::DIMS)
}

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("TA", kelvin(), 350.0)
        .unwrap()
        .with_quantity("TB", kelvin(), 320.0)
        .unwrap()
        .with_quantity("normal_pressure", pressure(), 1e6)
        .unwrap()
        .with_quantity("normal_separation", separation(), 1e-6)
        .unwrap()
}

fn contact(a: [f64; 2], b: [f64; 2], reverse: bool, response: Response) -> SourceFaceContact {
    let value = match response {
        Response::Constant => PropertyValue::Scalar {
            value: 0.5,
            dims: RD,
        },
        Response::Curve => PropertyValue::Curve {
            abscissa: "TA".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [a[0], a[1]]
                .map(|t| (t, 0.5 + (t - 320.0) / 1024.0))
                .to_vec(),
            dims: RD,
        },
        Response::Surface => PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()],
            axis_dims: [Temperature::DIMS; 2],
            samples: [[a[0], b[0]], [a[1], b[0]], [a[0], b[1]], [a[1], b[1]]]
                .map(|t| (t, resistance(t)))
                .to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: RD,
        },
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value,
            // Each validity axis carries the quantity that point() gives it: a query
            // axis typed differently from its claim's axis is outside the claim's domain.
            validity: ValidityDomain::unconstrained()
                .with_quantity("TA", kelvin(), a[0], a[1])
                .with_quantity("TB", kelvin(), b[0], b[1])
                .with_quantity("normal_pressure", pressure(), 1e6, 1e6)
                .with_quantity("normal_separation", separation(), 1e-6, 1e-6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if matches!(response, Response::Constant) {
                InterpolationPolicy::ConstantWithinValidity
            } else {
                InterpolationPolicy::LinearInside
            },
            observations: vec![],
            provenance: Provenance {
                source: "manufactured ordered contact for native tensor regression".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "fixture".into(),
            process: "ground finish".into(),
            revision: 0,
        },
        texture_frame: "fixed reference normal and surface finish".into(),
    };
    let [sa, sb] = if reverse { ["B", "A"] } else { ["A", "B"] };
    let card = Arc::new(
        InterfaceSystemCard::assemble(
            surface(sa),
            surface(sb),
            SystemContext {
                medium: "dry gap".into(),
                third_body: None,
                environment: "test".into(),
                history: "unaged".into(),
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let mut query = point();
    query = query
        .with("TA", 350.0_f64.clamp(a[0], a[1]))
        .unwrap()
        .with("TB", 320.0_f64.clamp(b[0], b[1]))
        .unwrap();
    SourceFaceContact::new(
        if reverse { [1, 0] } else { [0, 1] },
        card,
        &query,
        if reverse { ["TB", "TA"] } else { ["TA", "TB"] },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn setup(response: Response) -> ThermalSetup {
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup
        .source_face_contacts
        .push(contact([256.0, 384.0], [256.0, 384.0], false, response));
    setup
}

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        ..TensorNonlinearConfig::default()
    }
}

fn domain() -> VoxelDomain {
    VoxelDomain::from_fn(2, 1, 1, 1.0, |p| Voxel::Solid(u16::from(p[0] > 1.0))).unwrap()
}

fn tensors() -> [ConductivityTensor; 2] {
    [2.0, 4.0]
        .map(|k| ConductivityTensor::new([[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]).unwrap())
}

fn solids() -> [SolidMaterial; 2] {
    [SolidMaterial::new("A", 2.0), SolidMaterial::new("B", 4.0)]
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-7,
        "{actual:.16e} != {expected:.16e}"
    );
}

fn check_source(solution: &TensorConductionSolution, surfaces: [f64; 2], heat: f64) {
    let state = &solution.source_contact_states[0];
    assert_eq!(solution.source_contact_states.len(), 1);
    let reverse = state.cells == [1, 0];
    let expected = if reverse {
        [surfaces[1], surfaces[0]]
    } else {
        surfaces
    };
    for i in 0..2 {
        near(state.surface_temperature_k[i], expected[i]);
    }
    near(
        state.heat_rate_first_to_second_w,
        if reverse { -heat } else { heat },
    );
    near(state.source.state_coordinate("TA").unwrap(), surfaces[0]);
    near(state.source.state_coordinate("TB").unwrap(), surfaces[1]);
    assert_eq!(state.source.state_coordinate("normal_pressure"), Some(1e6));
    assert_eq!(
        state.source.state_coordinate("normal_separation"),
        Some(1e-6)
    );
    assert_eq!(state.resistance_relative_residual, 0.0);
    assert_eq!(state.resistance_m2_k_w, state.applied_resistance_m2_k_w);
    near(state.contact_temperature_residual_k, 0.0);
    near(solution.outward_heat_w[0][1], heat);
    near(solution.outward_heat_w[1][0], -heat);
    assert!(solution.report.max_heat_residual <= controls().linear.tolerance);
    assert!(heat * (surfaces[1].recip() - surfaces[0].recip()) >= 0.0);
}

#[test]
fn full_tensor_contact_retains_supported_surfaces_and_original_order() {
    // Adiabatic tangential faces imply K_eff=1/(K^-1)_xx, not K_xx.
    let tensors = [
        ConductivityTensor::new([[3.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0]]).unwrap(),
        ConductivityTensor::new([[4.0, -0.5, 0.0], [-0.5, 1.5, 0.0], [0.0, 0.0, 2.0]]).unwrap(),
    ];
    let surfaces = [350.0, 320.0];
    let q = 30.0 / resistance(surfaces);
    let centers = [surfaces[0] + 0.2 * q, surfaces[1] - 3.0 * q / 23.0];
    assert!(centers[0] > 352.0 && centers[1] < 318.0);
    let mut previous = None::<TensorConductionSolution>;
    for reverse in [false, true] {
        let mut setup = setup(Response::Surface);
        setup.fixed_temperature = vec![(0, centers[0]), (1, centers[1])];
        setup.source_face_contacts[0] =
            contact([348.0, 352.0], [318.0, 322.0], reverse, Response::Surface);
        let before = setup.clone();
        let solved = solve_nonlinear_tensor_conduction(
            &domain(),
            &solids(),
            &setup,
            &[&tensors[0], &tensors[1]],
            &centers,
            &controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let solution = solved.solution;
        check_source(&solution, surfaces, q);
        for cell in 0..2 {
            near(solution.temperature[cell], centers[cell]);
            for f in 2..6 {
                near(solution.outward_heat_w[cell][f], 0.0);
            }
        }
        near(solution.report.balance.residual_w, 0.0);
        if let Some(old) = previous {
            for cell in 0..2 {
                for face in 0..6 {
                    near(
                        old.face_temperature_k[cell][face],
                        solution.face_temperature_k[cell][face],
                    );
                    near(
                        old.outward_heat_w[cell][face],
                        solution.outward_heat_w[cell][face],
                    );
                }
            }
        }
        previous = Some(solution);
        assert_eq!(setup, before);
    }
}

// The only root unknown is actual contact heat. Both storage equations and
// each bulk half-cell drop are eliminated independently of the spatial solver.
fn sensible_oracle(interval: f64, carry: f64) -> ([f64; 2], [f64; 2], f64) {
    let values = |q: f64| {
        let energy = carry + interval * q;
        let t = [360.0 - energy / 2.0, 300.0 + energy / 5.0];
        (t, [t[0] - 0.25 * q, t[1] + 0.125 * q])
    };
    let (mut lo, mut hi) = (0.0, 100.0);
    for _ in 0..80 {
        let q = 0.5 * (lo + hi);
        let (_, s) = values(q);
        if s[0] - s[1] - resistance(s) * q > 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let (t, s) = values(q);
    (t, s, q)
}

#[test]
fn both_sensible_time_drivers_solve_original_source_without_double_counting_heat() {
    let mut setup = setup(Response::Surface);
    // Native Newton owns the budget; scalar outer Picard controls are unrelated.
    setup.source_contact_controls.max_iterations = 0;
    let before = setup.clone();
    let tensors = tensors();
    let previous = [360.0, 300.0];
    let dt = 0.35;
    let check = |solution: &TensorConductionSolution, expected: ([f64; 2], [f64; 2], f64)| {
        for cell in 0..2 {
            near(solution.temperature[cell], expected.0[cell]);
        }
        check_source(solution, expected.1, expected.2);
        near(
            2.0 * (solution.temperature[0] - previous[0])
                + 5.0 * (solution.temperature[1] - previous[1]),
            0.0,
        );
    };
    let step = advance_nonlinear_tensor_conduction(
        &domain(),
        &solids(),
        &setup,
        &[&tensors[0], &tensors[1]],
        &[2.0, 5.0],
        &previous,
        2.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    check(&step.solution.solution, sensible_oracle(dt, 0.0));
    near(step.record.stored_energy_change_j, 0.0);
    near(step.record.closure_j, 0.0);
    assert_eq!(step.record.source_j, 0.0);
    assert_eq!(step.record.boundary_outflow_j, 0.0);

    let step = advance_nonlinear_tensor_conduction_sdirk2(
        &domain(),
        &solids(),
        [&setup, &setup],
        &[&tensors[0], &tensors[1]],
        &[2.0, 5.0],
        &previous,
        2.0,
        dt,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let first = sensible_oracle(a, 0.0);
    check(&step.stages[0].solution, first);
    check(
        &step.stages[1].solution,
        sensible_oracle(a, (dt - a) * first.2),
    );
    near(step.record.stored_energy_change_j, 0.0);
    near(step.record.closure_j, 0.0);
    assert_eq!(step.record.source_j, 0.0);
    assert_eq!(step.record.boundary_outflow_j, 0.0);
    assert_eq!(setup, before);
    assert_eq!(previous, [360.0, 300.0]);
}

#[test]
fn original_constant_and_curve_contacts_replay_their_physical_laws() {
    let tensors = tensors();
    for response in [Response::Constant, Response::Curve] {
        let mut setup = setup(response);
        let surfaces = [350.0, 320.0];
        let r = if matches!(response, Response::Constant) {
            0.5
        } else {
            0.5 + 30.0 / 1024.0
        };
        let q = 30.0 / r;
        let centers = [surfaces[0] + 0.25 * q, surfaces[1] - 0.125 * q];
        setup.fixed_temperature = vec![(0, centers[0]), (1, centers[1])];
        let source = solve_nonlinear_tensor_conduction(
            &domain(),
            &solids(),
            &setup,
            &[&tensors[0], &tensors[1]],
            &centers,
            &controls(),
            &CancelGate::new(),
        )
        .unwrap();
        check_source(&source.solution, surfaces, q);
        near(
            source.solution.source_contact_states[0].resistance_m2_k_w,
            r,
        );
        if matches!(response, Response::Constant) {
            assert!(
                solve_tensor_conduction(
                    &domain(),
                    &solids(),
                    &setup,
                    &tensors,
                    &controls().linear,
                    &CancelGate::new()
                )
                .is_err()
            );
            setup.source_face_contacts.clear();
            setup.face_contacts.push(FaceContactResistance {
                cells: [0, 1],
                resistance_m2_k_w: r,
            });
            let fixed = solve_nonlinear_tensor_conduction(
                &domain(),
                &solids(),
                &setup,
                &[&tensors[0], &tensors[1]],
                &centers,
                &controls(),
                &CancelGate::new(),
            )
            .unwrap();
            assert!(fixed.solution.source_contact_states.is_empty());
            for cell in 0..2 {
                for face in 0..6 {
                    near(
                        source.solution.outward_heat_w[cell][face],
                        fixed.solution.outward_heat_w[cell][face],
                    );
                    near(
                        source.solution.face_temperature_k[cell][face],
                        fixed.solution.face_temperature_k[cell][face],
                    );
                }
            }
        }
    }
}

#[test]
fn ambiguous_unsupported_cancelled_and_exhausted_interfaces_refuse_atomically() {
    let setup = setup(Response::Surface);
    let before = setup.clone();
    let tensors = tensors();
    let previous = [360.0, 300.0];
    let solve = |setup: &ThermalSetup, config: &TensorNonlinearConfig, gate: &CancelGate| {
        advance_nonlinear_tensor_conduction(
            &domain(),
            &solids(),
            setup,
            &[&tensors[0], &tensors[1]],
            &[2.0, 5.0],
            &previous,
            0.0,
            0.35,
            config,
            gate,
        )
    };
    let mut face = setup.clone();
    face.face_contacts.push(FaceContactResistance {
        cells: [1, 0],
        resistance_m2_k_w: 0.0,
    });
    assert!(solve(&face, &controls(), &CancelGate::new()).is_err());
    let mut pair = setup.clone();
    pair.contacts.push(ContactResistance {
        materials: (1, 0),
        resistance_m2_k_w: 0.0,
    });
    assert!(solve(&pair, &controls(), &CancelGate::new()).is_err());
    let mut duplicate = setup.clone();
    duplicate.source_face_contacts.push(contact(
        [256.0, 384.0],
        [256.0, 384.0],
        true,
        Response::Surface,
    ));
    assert!(solve(&duplicate, &controls(), &CancelGate::new()).is_err());
    let mut unsupported = setup.clone();
    unsupported.source_face_contacts[0] =
        contact([319.0, 321.0], [319.0, 321.0], false, Response::Surface);
    assert!(solve(&unsupported, &controls(), &CancelGate::new()).is_err());
    let exhausted = TensorNonlinearConfig {
        max_newton_iterations: 0,
        ..controls()
    };
    assert!(matches!(
        solve(&setup, &exhausted, &CancelGate::new()),
        Err(ChtError::SolverNotConverged { .. })
    ));
    assert!(
        advance_nonlinear_tensor_conduction_sdirk2(
            &domain(),
            &solids(),
            [&setup, &unsupported],
            &[&tensors[0], &tensors[1]],
            &[2.0, 5.0],
            &previous,
            0.0,
            0.35,
            &controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(solve(&setup, &controls(), &gate), Err(ChtError::Cancelled));
    assert_eq!(setup, before);
    assert_eq!(previous, [360.0, 300.0]);
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn nonlinear_contact_and_original_latent_storage_share_one_conservative_stage() {
    use fs_blake3::ContentHash;
    use fs_lbm::conjugate::{
        TensorCaloricLaw, TensorCaloricStorage, advance_caloric_tensor_conduction,
        advance_caloric_tensor_conduction_sdirk2,
    };
    use fs_material::phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
    };
    let cold = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xC1; 32]),
        [
            (-100.0, 250.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (300.0, 400.0, 1.0),
        ]
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 2.0,
        })
        .to_vec(),
    )
    .unwrap();
    let hot = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        ContentHash([0xC2; 32]),
        SolidLiquidPhase::Solid,
        -100.0,
        &[250.0, 450.0].map(|t| HeatCapacityKnot {
            temperature_k: t,
            specific_heat_capacity_j_kg_k: 2.0,
            bulk_density_kg_m3: 3.0,
        }),
    )
    .unwrap();
    let caloric: [&dyn TensorCaloricLaw; 2] = [&cold, &hot];
    let previous = [20.0, 120.0];
    let mass = [2.0, 3.0];
    let tensor =
        ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
    let setup = setup(Response::Surface);
    let before = setup.clone();
    let oracle = |interval: f64, carry: f64| {
        let values = |q: f64| {
            let energy = carry + interval * q;
            let h = [20.0 + energy / 2.0, 120.0 - energy / 3.0];
            let t = [300.0 + (h[0] - 100.0).max(0.0) / 2.0, 300.0 + h[1] / 2.0];
            (h, [t[0] + 0.5 * q, t[1] - 0.5 * q])
        };
        let (mut lo, mut hi) = (0.0, 40.0);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            let (_, s) = values(q);
            if s[1] - s[0] - resistance(s) * q > 0.0 {
                lo = q;
            } else {
                hi = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let (h, s) = values(q);
        (h, s, q)
    };
    for dt in [1.0, 8.0] {
        for second_order in [false, true] {
            let storage = TensorCaloricStorage {
                laws: &caloric,
                reference_mass_kg: &mass,
                previous_specific_energy_j_kg: &previous,
            };
            let (endpoint, record, expected) = if second_order {
                let step = advance_caloric_tensor_conduction_sdirk2(
                    &domain(),
                    &solids(),
                    [&setup, &setup],
                    &[&tensor, &tensor],
                    storage,
                    [&[], &[]],
                    4.0,
                    dt,
                    &controls(),
                    &CancelGate::new(),
                )
                .unwrap();
                let a = TENSOR_SDIRK2_GAMMA * dt;
                let first = oracle(a, 0.0);
                for cell in 0..2 {
                    near(step.stages[0].specific_energy_j_kg[cell], first.0[cell]);
                }
                check_source(&step.stages[0].solution.solution, first.1, -first.2);
                let [_, last] = step.stages;
                (last, step.record, oracle(a, (dt - a) * first.2))
            } else {
                let step = advance_caloric_tensor_conduction(
                    &domain(),
                    &solids(),
                    &setup,
                    &[&tensor, &tensor],
                    storage,
                    &[],
                    4.0,
                    dt,
                    &controls(),
                    &CancelGate::new(),
                )
                .unwrap();
                (step.endpoint, step.record, oracle(dt, 0.0))
            };
            for cell in 0..2 {
                near(endpoint.specific_energy_j_kg[cell], expected.0[cell]);
            }
            check_source(&endpoint.solution.solution, expected.1, -expected.2);
            let phase = cold
                .values_at_specific_enthalpy(endpoint.specific_energy_j_kg[0])
                .unwrap();
            if dt == 1.0 {
                assert_eq!(phase.temperature_k, 300.0);
                assert_eq!(
                    cold.temperature_slope_kg_k_j(endpoint.specific_energy_j_kg[0])
                        .unwrap(),
                    0.0
                );
                assert!(phase.liquid_mass_fraction > 0.2 && phase.liquid_mass_fraction < 1.0);
            } else {
                assert!(phase.temperature_k > 300.0);
            }
            near(
                2.0 * (endpoint.specific_energy_j_kg[0] - previous[0])
                    + 3.0 * (endpoint.specific_energy_j_kg[1] - previous[1]),
                0.0,
            );
            near(record.stored_energy_change_j, 0.0);
            near(record.closure_j, 0.0);
            assert_eq!(record.source_j, 0.0);
            assert_eq!(record.boundary_outflow_j, 0.0);
        }
    }
    assert_eq!(setup, before);
    assert_eq!(previous, [20.0, 120.0]);
}
