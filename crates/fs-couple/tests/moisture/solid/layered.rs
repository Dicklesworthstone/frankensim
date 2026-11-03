//! G1 heterogeneous compatible swelling and G3 conservative physical coupling.
//! Authored coefficients verify the implementation, not a measured laminate.
use super::*;
use fs_couple::moisture::layered::{MoistureLayer, MoistureLayeredFace, MoistureLayeredSlab};
use fs_material::state_point::{
    MonotoneMoistureSorption, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

#[path = "layered/deforming.rs"]
mod deforming;

fn material(rho: f64, beta: [f64; 3], upper: f64, saturation: f64) -> MaterialCard {
    let original = card_with_swelling([0.01; 2], upper, beta);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Scalar {
                value: 1.1 * rho,
                dims: Density::DIMS,
            };
        }
        if matches!(claim.value, PropertyValue::Scalar { .. }) {
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(
                "sorption",
                qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ),
            value: PropertyValue::Curve {
                abscissa: "relative_humidity".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 0.0), (1.0, saturation)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    293.15,
                    293.15,
                )
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.0,
                    1.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored layered solid sorption".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn layer(card: MaterialCard, length: f64, moisture: Vec<f64>, cx: &Cx<'_>) -> MoistureLayer {
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    let state = resolve_material_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_with_quantity(
            "sorption",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let sorption =
        MonotoneMoistureSorption::from_binding(&card, &state, "sorption", [0.0, 1.0]).unwrap();
    MoistureLayer::new(slab(card, length, moisture, cx), sorption, cx).unwrap()
}

fn sources(count: usize) -> Vec<MoistureLayerSolidSource> {
    vec![
        MoistureLayerSolidSource {
            stiffness: SolidStiffnessModel::Isotropic,
            elastic_selection: MaterialPropertySelection::SingleClaimOnly,
            expansion_selection: MaterialPropertySelection::SingleClaimOnly,
            axes: axes()
        };
        count
    ]
}
fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}
fn bind_stack(stack: &MoistureLayeredSlab, fixed: Vec<usize>, cx: &Cx<'_>) -> MoistureLayeredSolid {
    MoistureLayeredSolid::new(
        stack,
        0.01,
        sources(stack.layers().len()),
        fixed,
        controls(),
        cx,
    )
    .unwrap()
}

#[test]
fn g1_identical_layer_split_reproduces_the_original_native_solid() {
    with_cx(|cx| {
        let card = material(500.0, [0.01, 0.04, 0.02], 0.3, 0.3);
        let whole = slab(card.clone(), 0.2, vec![0.2, 0.15], cx);
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(card.clone(), 0.1, vec![0.2], cx),
                layer(card, 0.1, vec![0.15], cx),
            ],
            256,
            cx,
        )
        .unwrap();
        let original = bind(
            &whole,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        let layered = bind_stack(&stack, free_supports(), cx);
        let a = original.resolve(&whole, cx).unwrap();
        let b = layered.resolve(&stack, cx).unwrap();
        assert_eq!(original.reference_nodes_m(), layered.reference_nodes_m());
        assert_eq!(original.tetrahedra(), layered.tetrahedra());
        assert_eq!(layered.layer_cells(), &[0..1, 1..2]);
        for (u, v) in a
            .solution()
            .displacement_m()
            .iter()
            .zip(b.solution().displacement_m())
        {
            for k in 0..3 {
                near(u[k], v[k], 1e-12);
            }
        }
        near(
            a.stress().elastic_energy_j(),
            b.stress().elastic_energy_j(),
            1e-9,
        );
        near(
            a.tangent().unwrap().pencil().total_mass_kg,
            b.wet_mass_kg(),
            1e-13,
        );
        for (left, right) in a.materials().iter().zip(b.materials()) {
            near(left.density_kg_m3, right.density_kg_m3, 1e-10);
        }
    });
}

#[test]
fn g1_unlike_materials_and_rotated_directions_share_compatible_axial_displacement() {
    with_cx(|cx| {
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(500.0, [0.01, 0.0, 0.0], 0.3, 0.3),
                    0.1,
                    vec![0.2],
                    cx,
                ),
                layer(
                    material(800.0, [0.0, 0.04, 0.0], 0.3, 0.2),
                    0.2,
                    vec![0.15, 0.15],
                    cx,
                ),
            ],
            256,
            cx,
        )
        .unwrap();
        let mut selected = sources(2);
        selected[1].axes.principal_to_mesh = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        selected[1].axes.principal_frame = ContentHash([8; 32]);
        selected[1].stiffness = SolidStiffnessModel::OrientedOrthotropic {
            principal_to_mesh: selected[1].axes.principal_to_mesh,
            orientation_identity: ContentHash([9; 32]),
        };
        let solid =
            MoistureLayeredSolid::new(&stack, 0.01, selected, free_supports(), controls(), cx)
                .unwrap();
        let state = solid.resolve(&stack, cx).unwrap();
        assert_eq!(solid.reference_nodes_m().len(), 16);
        assert_eq!(solid.tetrahedra().len(), 18);
        assert_eq!(solid.layer_cells(), &[0..1, 1..3]);
        for (x, u) in solid
            .reference_nodes_m()
            .iter()
            .zip(state.solution().displacement_m())
        {
            let expected = 0.001 * x[0].min(0.1) + 0.002 * (x[0] - 0.1).max(0.0);
            near(u[0], expected, 2e-10);
            near(u[1], 0.0, 2e-10);
            near(u[2], 0.0, 2e-10);
        }
        for element in state.stress().elements() {
            for &value in element.stress_mandel_pa() {
                near(value, 0.0, 2.0);
            }
        }
        for (k, range) in solid.layer_cells().iter().enumerate() {
            for i in range.clone() {
                assert_eq!(
                    state.elastic_states()[i].card_identity(),
                    stack.layers()[k].slab().card().content_hash()
                );
            }
        }
        near(
            state.wet_mass_kg(),
            0.5 * 0.1 * 1.2 + 0.8 * 0.2 * 1.15,
            1e-13,
        );
        near(
            state.tangent().unwrap().pencil().total_mass_kg,
            stack.wet_mass_kg(),
            1e-13,
        );
        near(
            state
                .current_density_kg_m3()
                .iter()
                .zip(state.geometry().element_volumes_m3())
                .map(|(rho, volume)| rho * volume)
                .sum(),
            stack.wet_mass_kg(),
            1e-13,
        );
    });
}

#[test]
fn g1_each_constrained_layer_retains_its_own_stress_energy_and_wet_mass() {
    with_cx(|cx| {
        let betas = [[0.01, 0.04, 0.02], [0.03, -0.01, 0.01]];
        let moisture = [0.2, 0.15];
        let lengths = [0.1, 0.2];
        let stack = MoistureLayeredSlab::new(
            (0..2)
                .map(|k| {
                    layer(
                        material(500.0 + 300.0 * k as f64, betas[k], 0.3, 0.3),
                        lengths[k],
                        vec![moisture[k]],
                        cx,
                    )
                })
                .collect(),
            256,
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, (0..36).collect(), cx);
        let state = solid.resolve(&stack, cx).unwrap();
        let mut energy = 0.0;
        for k in 0..2 {
            let e = state.elastic_states()[k]
                .property("young_modulus")
                .unwrap()
                .value_si();
            let free = betas[k].map(|beta| beta * (moisture[k] - 0.1));
            let trace = free.iter().sum::<f64>();
            let stress = free.map(|value| -e * (0.4 * trace + 0.8 * value));
            for element in &state.stress().elements()[6 * k..6 * (k + 1)] {
                for j in 0..3 {
                    near(element.stress_mandel_pa()[j], stress[j], 1e-5);
                }
            }
            energy += -0.5 * lengths[k] * 0.001 * (0..3).map(|j| free[j] * stress[j]).sum::<f64>();
        }
        near(state.stress().elastic_energy_j(), energy, 1e-9);
        assert!(state.tangent().is_none());
        assert!(
            state
                .solution()
                .displacement_m()
                .iter()
                .all(|u| *u == [0.0; 3])
        );
    });
}

#[test]
fn g1_layered_water_uptake_changes_the_native_generalized_eigenvalue() {
    with_cx(|cx| {
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(material(500.0, [0.01; 3], 0.3, 0.3), 0.1, vec![0.1], cx),
                layer(material(800.0, [0.02; 3], 0.3, 0.2), 0.1, vec![0.1], cx),
            ],
            256,
            cx,
        )
        .unwrap();
        // The one free interface displacement belongs to BOTH actual layers.
        let binding = bind_stack(&stack, (0..36).filter(|d| *d != 12).collect(), cx);
        let before = binding.resolve(&stack, cx).unwrap();
        let after = binding
            .advance(
                &stack,
                0.25,
                &[MoistureLayeredFace::PrescribedRelativeHumidity(0.9); 2],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        let eigenvalue = |state: &MoistureSolidState| {
            let pencil = state.tangent().unwrap().pencil();
            assert_eq!(pencil.free_dofs, &[12]);
            pencil.stiffness.get(0, 0) / pencil.mass.get(0, 0)
        };
        assert!(eigenvalue(&after.solid) > 0.0);
        assert!(eigenvalue(&after.solid) < eigenvalue(&before));
        assert!(after.solid.wet_mass_kg() > before.wet_mass_kg());
        assert!(after.solid.stress().elastic_energy_j() > 0.0);
        assert_ne!(after.solid.solution().displacement_m()[4][0], 0.0);
    });
}

#[test]
fn g3_finite_vapor_transport_changes_native_response_and_late_source_failure_rolls_back() {
    with_cx(|cx| {
        let make = |upper| {
            MoistureLayeredSlab::new(
                vec![
                    layer(material(500.0, [0.01; 3], upper, 0.3), 0.1, vec![0.1], cx),
                    layer(material(800.0, [0.02; 3], 0.3, 0.2), 0.1, vec![0.1], cx),
                ],
                256,
                cx,
            )
            .unwrap()
        };
        let stack = make(0.3);
        let solid = bind_stack(&stack, free_supports(), cx);
        let before = solid.resolve(&stack, cx).unwrap();
        let gas = atmosphere::atmosphere(0.9, 1.0);
        let (candidate, after) = solid
            .advance_finite_atmosphere(
                &gas,
                &stack,
                0.5,
                [true; 2],
                config(),
                nonlinear(),
                atmosphere::controls(),
                cx,
            )
            .unwrap();
        assert!(candidate.atmosphere().vapor_mass_kg().unwrap() < gas.vapor_mass_kg().unwrap());
        assert!(after.wet_mass_kg() > before.wet_mass_kg());
        assert_eq!(after.time_s(), candidate.moisture().time_s());
        near(
            after.wet_mass_kg(),
            candidate.moisture().wet_mass_kg(),
            1e-13,
        );
        assert_ne!(
            before.tangent().unwrap().identity(),
            after.tangent().unwrap().identity()
        );
        assert!(
            after.geometry().element_volumes_m3().iter().sum::<f64>()
                > before.geometry().element_volumes_m3().iter().sum::<f64>()
        );

        let limited = make(0.1001);
        let binding = bind_stack(&limited, free_supports(), cx);
        let water = gas
            .advance_layered(
                &limited,
                0.5,
                [true; 2],
                config(),
                nonlinear(),
                atmosphere::controls(),
                cx,
            )
            .unwrap();
        assert!(water.moisture().layers()[0].slab().moisture()[0] > 0.1001);
        assert!(
            binding
                .advance_finite_atmosphere(
                    &gas,
                    &limited,
                    0.5,
                    [true; 2],
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(limited.time_s(), 0.0);
        assert_eq!(limited.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(limited.layers()[1].slab().moisture(), &[0.1]);
        near(gas.relative_humidity().unwrap(), 0.9, 1e-13);
        assert!(binding.resolve(&stack, cx).is_err());
        let mut wrong = sources(2);
        wrong[1].axes.mesh_frame = ContentHash([11; 32]);
        assert!(
            MoistureLayeredSolid::new(&stack, 0.01, wrong, free_supports(), controls(), cx)
                .is_err()
        );
    });
}

#[test]
fn g3_adaptive_windows_publish_exact_time_mechanics_and_discard_late_failure() {
    use fs_couple::moisture::adaptive::{MoistureAdvanceError, MoistureTimeControls};
    with_cx(|cx| {
        let time = MoistureTimeControls {
            initial_step_s: 0.025,
            minimum_step_s: 1e-8,
            maximum_step_s: 0.025,
            maximum_evaluations: 1000,
            moisture_absolute_tolerance: 1e-4,
            activity_absolute_tolerance: 1e-4,
            reservoir_water_absolute_tolerance_kg: 1e-7,
            relative_tolerance: 0.0,
            water_balance_tolerance_kg: 1e-9,
        };
        let make = |upper| {
            MoistureLayeredSlab::new(
                vec![
                    layer(material(500.0, [0.01; 3], upper, 0.3), 0.1, vec![0.1], cx),
                    layer(material(800.0, [0.02; 3], 0.3, 0.2), 0.1, vec![0.1], cx),
                ],
                256,
                cx,
            )
            .unwrap()
        };
        let stack = make(0.3);
        let binding = bind_stack(&stack, free_supports(), cx);
        let fixed = binding
            .advance_adaptive_to(
                &stack,
                0.1,
                &[MoistureLayeredFace::PrescribedRelativeHumidity(0.9); 2],
                config(),
                nonlinear(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(fixed.solid.time_s(), 0.1);
        assert!(fixed.transport.summary.accepted_substeps >= 2);
        near(
            fixed.solid.wet_mass_kg(),
            fixed.moisture.wet_mass_kg(),
            1e-13,
        );
        let gas = atmosphere::atmosphere(0.9, 0.1);
        let (next_gas, finite) = binding
            .advance_finite_atmosphere_adaptive_to(
                &gas,
                &stack,
                0.1,
                [true; 2],
                config(),
                nonlinear(),
                atmosphere::controls(),
                time,
                cx,
            )
            .unwrap();
        assert_eq!(finite.solid.time_s(), 0.1);
        near(
            next_gas.vapor_mass_kg().unwrap() + finite.moisture.water_mass_kg(),
            gas.vapor_mass_kg().unwrap() + stack.water_mass_kg(),
            1e-9,
        );
        let restricted = make(0.1001);
        let constrained = bind_stack(&restricted, free_supports(), cx);
        assert!(matches!(
            constrained.advance_finite_atmosphere_adaptive_to(
                &gas,
                &restricted,
                0.1,
                [true; 2],
                config(),
                nonlinear(),
                atmosphere::controls(),
                time,
                cx
            ),
            Err(MoistureAdvanceError::FinalInventory(_))
        ));
        assert_eq!(restricted.time_s(), 0.0);
        assert_eq!(restricted.layers()[0].slab().moisture(), &[0.1]);
        near(gas.relative_humidity().unwrap(), 0.9, 1e-13);
    });
}
