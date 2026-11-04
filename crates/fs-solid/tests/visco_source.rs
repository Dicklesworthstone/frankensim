//! G1/G3: original source changes, physical Maxwell memories and separate work.
//! These manufactured cards are analytic software fixtures, not measured data.

#[path = "visco_source/dynamics.rs"]
mod dynamics;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::visco_tensor::source::{
    ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity,
    TensorKelvinViscosityProperties, TensorPronyProperties,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Temperature};
use fs_solid::linear::Formulation;
use fs_solid::linear3::visco::source::{ViscoMaterialSource, ViscoSourceError};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticError, TetStaticSolveConfig,
    TetThermalStrainState,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const HEX: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const HEXES: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 926,
                    kernel_id: 926,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite()
            && (actual - expected).abs() <= 4e-10 * actual.abs().max(expected.abs()).max(1e-10),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("omega", 0.0)
        .unwrap()
}
fn claim(name: &str, dims: Dims, values: [f64; 3]) -> PropertyClaim {
    let (value, interpolation) = if values[0] == values[1] && values[1] == values[2] {
        (
            PropertyValue::Scalar {
                value: values[0],
                dims,
            },
            InterpolationPolicy::ConstantWithinValidity,
        )
    } else {
        (
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [300.0, 350.0, 400.0].into_iter().zip(values).collect(),
                dims,
            },
            InterpolationPolicy::LinearInside,
        )
    };
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        interpolation,
        validity: ValidityDomain::unconstrained()
            .with("T", 300.0, 400.0)
            .with("omega", 0.0, 1000.0),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured source-continuation fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn properties() -> (TensorPronyProperties, TensorKelvinViscosityProperties) {
    (
        TensorPronyProperties {
            elasticity: TensorElasticity::Isotropic,
            strain_limit: "strain_limit".into(),
            branches: vec![TensorBranchProperties {
                stiffness_upper: core::array::from_fn(|i| format!("branch_{i:02}")),
                relaxation_time: "tau".into(),
            }],
        },
        TensorKelvinViscosityProperties {
            viscosity_upper: core::array::from_fn(|i| format!("kelvin_{i:02}")),
        },
    )
}
fn card(changing_stiffness: bool) -> MaterialCard {
    let (properties, kelvin) = properties();
    let mut claims = ClaimSet::new();
    for entry in [
        claim(
            "young_modulus",
            Pressure::DIMS,
            if changing_stiffness {
                [1200.0, 600.0, 1200.0]
            } else {
                [1200.0; 3]
            },
        ),
        claim("poisson_ratio", Dims::NONE, [0.0; 3]),
        claim("density", Density::DIMS, [1000.0, 900.0, 800.0]),
        claim("strain_limit", Dims::NONE, [0.03, 0.03, 0.0003]),
        claim("tau", Dims([0, 0, 1, 0, 0, 0]), [0.07, 0.14, 0.035]),
    ] {
        claims.insert_claim(entry).unwrap();
    }
    let mut index = 0;
    for i in 0..6 {
        for j in i..6 {
            let values = match (i, j, changing_stiffness) {
                (0, 0, true) => [600.0, 300.0, 600.0],
                (0, 0, false) => [600.0; 3],
                (1, 1, true) => [0.0, 200.0, 0.0],
                _ => [0.0; 3],
            };
            claims
                .insert_claim(claim(
                    &properties.branches[0].stiffness_upper[index],
                    Pressure::DIMS,
                    values,
                ))
                .unwrap();
            claims
                .insert_claim(claim(
                    &kelvin.viscosity_upper[index],
                    DynViscosity::DIMS,
                    if i == 0 && j == 0 {
                        [8.0, 4.0, 8.0]
                    } else {
                        [0.0; 3]
                    },
                ))
                .unwrap();
            index += 1;
        }
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "analytic coefficient continuation".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard) -> ViscoMaterialSource {
    let (properties, kelvin) = properties();
    ViscoMaterialSource {
        material: ResolvedTensorMaterial::resolve(
            card,
            &point(300.0),
            properties,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .with_kelvin_viscosity(card, kelvin, MaterialPropertySelection::SingleClaimOnly)
        .unwrap(),
        material_to_mesh: IDENTITY,
    }
}
fn prepare(
    source: &ViscoMaterialSource,
    mesh: SolidMeshRef<'_>,
    fixed: &[usize],
    cx: &Cx<'_>,
) -> ViscoElements {
    let material =
        TetViscoMaterial::from_source(&source.material, source.material_to_mesh).unwrap();
    ViscoElements::new(
        mesh,
        TetViscoField::Uniform(&material),
        fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 10_000,
        },
        cx,
    )
    .unwrap()
    .with_material_sources(vec![source.clone()], cx)
    .unwrap()
}
fn tet(source: &ViscoMaterialSource, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    prepare(
        source,
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &TETS,
        },
        &fixed,
        cx,
    )
}
fn phi(h: f64, tau: f64) -> f64 {
    tau / h * (1.0 - (-h / tau).exp())
}
fn solve_config() -> TetStaticSolveConfig {
    TetStaticSolveConfig {
        maximum_iterations: 100,
        relative_residual_tolerance: 2e-11,
    }
}

#[test]
fn original_sparse_history_preserves_singular_branch_memory_and_reference_mass() {
    with_cx(|cx, _| {
        let card = card(true);
        let source = source(&card);
        let original = tet(&source, &[3, 7], cx);
        let next = original
            .at_query_points(&[&card], &[point(350.0)], cx)
            .unwrap();
        let mut state = original.state(cx).unwrap();
        let (a, b, h) = (0.002, -0.001, 0.02);
        let mut u = [0.0; 12];
        u[3] = a;
        u[7] = b;
        state.step_internal_force(&u, h, cx).unwrap();
        let before = state.energy();
        let q = [phi(h, 0.07) * a, phi(h, 0.07) * b];
        let changed = state.rebind_materials(&next, cx).unwrap();
        let work =
            0.5 * VOLUME * (-600.0 * (a * a + b * b) - 300.0 * q[0] * q[0] + 200.0 * q[1] * q[1]);
        close(changed.parameter_work_j, work);
        close(changed.history.parameter_work_j(), work);
        assert_eq!(changed.history.displacement_m(), &u);
        assert_eq!(changed.history.energy().work_j, before.work_j);
        assert_eq!(changed.history.energy().dissipated_j, before.dissipated_j);
        // The branch was initially rank one. Its previously unstressed y spring
        // strain must survive when the original source activates that direction.
        let stored =
            0.5 * VOLUME * (600.0 * (a * a + b * b) + 300.0 * q[0] * q[0] + 200.0 * q[1] * q[1]);
        close(changed.history.energy().stored_j, stored);
        for i in 0..original.equilibrium().mass.nrows() {
            assert_eq!(
                original.equilibrium().mass.row(i),
                next.equilibrium().mass.row(i)
            );
        }
        assert_eq!(
            next.equilibrium().total_mass_kg,
            original.equilibrium().total_mass_kg
        );
        close(next.reference_material(0).unwrap().density_kg_m3, 1000.0);
        close(
            next.material_sources().unwrap()[0].material.density_kg_m3(),
            900.0,
        );
        close(
            next.reference_material_sources().unwrap()[0]
                .material
                .density_kg_m3(),
            1000.0,
        );
        close(next.equilibrium().stiffness.get(0, 0), 600.0 * VOLUME);
        assert_eq!(
            next.material_sources().unwrap()[0].material_to_mesh,
            IDENTITY
        );
        assert_eq!(state.energy(), before);
        let mut continued = changed.history;
        let duration: f64 = 0.05;
        let decay = (-duration / 0.14).exp();
        let held = continued.step_internal_force(&u, duration, cx).unwrap();
        close(
            held.full_force_n[3],
            VOLUME * (600.0 * a + 300.0 * q[0] * decay),
        );
        close(
            held.full_force_n[7],
            VOLUME * (600.0 * b + 200.0 * q[1] * decay),
        );
        close(
            held.step_dissipated_j,
            0.5 * VOLUME * (300.0 * q[0] * q[0] + 200.0 * q[1] * q[1]) * (1.0 - decay * decay),
        );
        assert_eq!(held.step_work_j, 0.0);
        close(
            held.energy.work_j + continued.parameter_work_j(),
            held.energy.stored_j + held.energy.dissipated_j,
        );
    });
}

#[test]
fn changed_relaxation_and_kelvin_rates_have_zero_parameter_storage_work() {
    with_cx(|cx, _| {
        let card = card(false);
        let source = source(&card);
        let original = tet(&source, &[3], cx);
        let next = original
            .at_query_points(&[&card], &[point(350.0)], cx)
            .unwrap();
        let mut state = original.state(cx).unwrap();
        let mut u = [0.0; 12];
        u[3] = 0.002;
        state.step_internal_force(&u, 0.02, cx).unwrap();
        let old = state.energy();
        let changed = state.rebind_materials(&next, cx).unwrap();
        assert_eq!(changed.parameter_work_j, 0.0);
        assert_eq!(changed.history.energy(), old);
        let q0 = phi(0.02, 0.07) * u[3];
        let mut history = changed.history;
        let (du, h, tau, stiffness, viscosity): (f64, f64, f64, f64, f64) =
            (0.0003, 0.03, 0.14, 600.0, 4.0);
        let decay = (-h / tau).exp();
        let q1 = decay * q0 + phi(h, tau) * du;
        let stress0 = stiffness * q0;
        let constant = stiffness * tau * du / h;
        let transient = stress0 - constant;
        let heat = VOLUME
            * (constant * constant * h / (stiffness * tau)
                + 2.0 * constant * transient * (1.0 - decay) / stiffness
                + transient * transient * (1.0 - decay * decay) / (2.0 * stiffness)
                + viscosity * du * du / h);
        u[3] += du;
        let response = history.step_internal_force(&u, h, cx).unwrap();
        close(
            response.full_force_n[3],
            VOLUME * (1200.0 * u[3] + stiffness * q1 + viscosity * du / h),
        );
        close(response.step_dissipated_j, heat);
        assert_eq!(history.parameter_work_j(), 0.0);
    });
}

#[test]
fn source_softening_continues_native_quasistatic_force_control() {
    with_cx(|cx, _| {
        let card = card(true);
        let source = source(&card);
        let original = tet(&source, &[3], cx);
        let next = original
            .at_query_points(&[&card], &[point(350.0)], cx)
            .unwrap();
        let h = 0.02;
        let old_operator = original.quasistatic_operator(h, cx).unwrap();
        let new_operator = next.quasistatic_operator(h, cx).unwrap();
        let mut state = original.state(cx).unwrap();
        let mut force = [0.0; 12];
        force[3] = 0.1;
        let first = state
            .step_quasistatic(&old_operator, &force, &[0.0; 12], solve_config(), cx)
            .unwrap();
        let u0 = first.displacement_m[3];
        close(
            u0,
            (force[3] / VOLUME) / (1200.0 + 600.0 * phi(h, 0.07) + 8.0 / h),
        );
        let q0 = phi(h, 0.07) * u0;
        let mut continued = state.rebind_materials(&next, cx).unwrap().history;
        assert!(
            continued
                .step_quasistatic(&old_operator, &force, &[0.0; 12], solve_config(), cx,)
                .is_err()
        );
        let decay = (-h / 0.14).exp();
        let increment = (force[3] / VOLUME - 600.0 * u0 - 300.0 * decay * q0)
            / (600.0 + 300.0 * phi(h, 0.14) + 4.0 / h);
        let response = continued
            .step_quasistatic(&new_operator, &force, &[0.0; 12], solve_config(), cx)
            .unwrap();
        close(response.displacement_m[3], u0 + increment);
        close(response.constitutive.full_force_n[3], force[3]);
        close(
            response.constitutive.energy.work_j + continued.parameter_work_j(),
            response.constitutive.energy.stored_j + response.constitutive.energy.dissipated_j,
        );
        assert!(response.true_relative_residual <= solve_config().relative_residual_tolerance);
    });
}

fn thermal(identity: [u8; 32], free: f64, index: usize) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        350.0,
        348.0 - 2.0 * index as f64,
        [free, 0.0, 0.0, 0.0, 0.0, 0.0],
        identity,
        ContentHash([73; 32]),
        0.03,
    )
    .unwrap()
}

#[test]
fn source_rebind_preserves_every_thermal_sample_and_original_bbar_projection() {
    with_cx(|cx, _| {
        let card = card(true);
        let source = source(&card);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let original = prepare(
                &source,
                SolidMeshRef::Hex8 {
                    nodes_m: &HEX,
                    elements: &HEXES,
                    formulation,
                    stiffness_points_per_dir: 2,
                },
                &[],
                cx,
            );
            let next = original
                .at_query_points(&[&card], &[point(350.0)], cx)
                .unwrap();
            let old_id = original
                .reference_material(0)
                .unwrap()
                .material_state_identity;
            let new_id = next.reference_material(0).unwrap().material_state_identity;
            let free: Vec<_> = (0..8).map(|q| 0.0002 * (q + 1) as f64).collect();
            let old: Vec<_> = free
                .iter()
                .enumerate()
                .map(|(q, &e)| thermal(old_id, e, q))
                .collect();
            let replacement: Vec<_> = free
                .iter()
                .enumerate()
                .map(|(q, &e)| thermal(new_id, e, q))
                .collect();
            fn field(states: &[TetThermalStrainState]) -> QuadraticThermalStrainField<'_> {
                QuadraticThermalStrainField::AtQuadrature {
                    states,
                    points_per_element: 8,
                }
            }
            let mut state = original.state(cx).unwrap();
            state
                .step_internal_force_with_thermal_strain(&[0.0; 24], 0.02, field(&old), cx)
                .unwrap();
            assert!(state.rebind_materials(&next, cx).is_err());
            let wrong_temperature = original
                .at_query_points(&[&card], &[point(325.0)], cx)
                .unwrap();
            let wrong_id = wrong_temperature
                .reference_material(0)
                .unwrap()
                .material_state_identity;
            let wrong_receipts: Vec<_> = free
                .iter()
                .enumerate()
                .map(|(q, &e)| thermal(wrong_id, e, q))
                .collect();
            assert!(
                state
                    .rebind_materials_with_thermal_strain(
                        &wrong_temperature,
                        field(&wrong_receipts),
                        cx,
                    )
                    .is_err()
            );
            assert!(
                state
                    .rebind_materials_with_thermal_strain(&next, field(&old), cx)
                    .is_err()
            );
            let mut altered = replacement.clone();
            altered[3] = thermal(new_id, free[3] + 1e-6, 3);
            assert!(
                state
                    .rebind_materials_with_thermal_strain(&next, field(&altered), cx)
                    .is_err()
            );
            // The held-coefficient history admits a spatial temperature field.
            // One elementwise source query cannot replace it with an average.
            let mut gradient = old.clone();
            gradient[7] = TetThermalStrainState::try_new(
                351.0,
                old[7].reference_temperature_k(),
                [free[7] + 1e-4, 0.0, 0.0, 0.0, 0.0, 0.0],
                old_id,
                old[7].thermal_law_identity(),
                old[7].maximum_absolute_mandel_strain(),
            )
            .unwrap();
            let mut nonuniform = state.clone();
            nonuniform
                .step_internal_force_with_thermal_strain(&[0.0; 24], 0.01, field(&gradient), cx)
                .unwrap();
            let rebound_gradient: Vec<_> = gradient
                .iter()
                .map(|s| {
                    TetThermalStrainState::try_new(
                        s.temperature_k(),
                        s.reference_temperature_k(),
                        s.free_strain_mandel(),
                        new_id,
                        s.thermal_law_identity(),
                        s.maximum_absolute_mandel_strain(),
                    )
                    .unwrap()
                })
                .collect();
            assert!(matches!(
                nonuniform.rebind_materials_with_thermal_strain(
                    &next,
                    field(&rebound_gradient),
                    cx,
                ),
                Err(ViscoSourceError::Invalid(
                    "element coefficient temperature must match every held thermal quadrature state"
                ))
            ));
            let changed = state
                .rebind_materials_with_thermal_strain(&next, field(&replacement), cx)
                .unwrap();
            assert_eq!(
                changed.history.thermal_strains().unwrap(),
                replacement.as_slice()
            );
            assert!(matches!(
                changed.history.thermal_field(),
                Some(QuadraticThermalStrainField::AtQuadrature {
                    points_per_element: 8,
                    ..
                })
            ));
            assert_eq!(
                changed.history.eigenstrain_work_j(),
                state.eigenstrain_work_j()
            );
            let average = free.iter().sum::<f64>() / 8.0;
            let ramp = phi(0.02, 0.07);
            let work = free
                .iter()
                .map(|&e| {
                    let strain = if formulation == Formulation::BBar {
                        [
                            (2.0 * e + average) / 3.0,
                            (average - e) / 3.0,
                            (average - e) / 3.0,
                        ]
                    } else {
                        [e, 0.0, 0.0]
                    };
                    0.5 / 8.0
                        * (-600.0 * strain.iter().map(|e| e * e).sum::<f64>()
                            + ramp
                                * ramp
                                * (-300.0 * strain[0].powi(2) + 200.0 * strain[1].powi(2)))
                })
                .sum::<f64>();
            close(changed.parameter_work_j, work);
            let mut history = changed.history;
            let held = history.step_internal_force(&[0.0; 24], 0.03, cx).unwrap();
            assert_eq!(held.step_eigenstrain_work_j, 0.0);
            assert_eq!(held.step_work_j, 0.0);
            close(
                held.energy.work_j + history.parameter_work_j(),
                held.energy.stored_j + held.energy.dissipated_j,
            );
        }
    });
}

#[test]
fn foreign_sources_geometry_domains_and_cancellation_refuse_without_losing_history() {
    with_cx(|cx, gate| {
        let card = card(true);
        let source = source(&card);
        let original = tet(&source, &[3], cx);
        let independent = tet(&source, &[3], cx);
        let mut state = original.state(cx).unwrap();
        let mut u = [0.0; 12];
        u[3] = 0.002;
        state.step_internal_force(&u, 0.02, cx).unwrap();
        let before = state.energy();
        assert!(state.rebind_materials(&independent, cx).is_err());
        assert!(
            original
                .at_query_points(&[&card], &[point(450.0)], cx)
                .is_err()
        );
        assert!(original.at_query_points(&[], &[], cx).is_err());
        let foreign = self::card(false);
        assert!(
            original
                .at_query_points(&[&foreign], &[point(350.0)], cx)
                .is_err()
        );
        let narrow = original
            .at_query_points(&[&card], &[point(400.0)], cx)
            .unwrap();
        assert!(state.rebind_materials(&narrow, cx).is_err());
        assert_eq!(state.energy(), before);
        assert_eq!(state.parameter_work_j(), 0.0);
        assert_eq!(state.displacement_m(), &u);
        gate.request();
        assert!(matches!(
            original.at_query_points(&[&card], &[point(350.0)], cx),
            Err(ViscoSourceError::Native(ViscoElementError::Elastic(
                TetElasticError::Cancelled
            )))
        ));
        assert_eq!(state.energy(), before);
    });
}
