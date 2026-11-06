//! Source natural strain reaches the actual conservative quadratic solid.
use super::*;
use fs_couple::finite_solid::thermal::MeanThermalReference;
use fs_couple::finite_solid::thermal::conservative::{
    ConservativeThermalError, ConservativeThermalField, ConservativeThermalPreparation,
};
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::conservative::ConservativeFormulation;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn preparation(formulation: ConservativeFormulation) -> ConservativeThermalPreparation {
    ConservativeThermalPreparation {
        formulation,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 64,
        maximum_mass_assembly_bytes: 128 * 1024,
        material_settings: DiscreteGradientSettings {
            minimum_relative_strain_increment: 1e-6,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-12,
        },
    }
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 15.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn mean_material(card: &MaterialCard, endpoint_k: f64, reverse_roles: bool) -> FiniteSolidMaterial {
    let (reference_axis, endpoint_axis) = if reverse_roles {
        ("T_endpoint", "T_reference")
    } else {
        ("T_reference", "T_endpoint")
    };
    let quantity = QuantitySpec::dimensional(Temperature::DIMS);
    let query = QueryPoint::new()
        .with_quantity(reference_axis, quantity, 300.0)
        .unwrap()
        .with_quantity(endpoint_axis, quantity, endpoint_k)
        .unwrap();
    let thermal = MeanThermalReference::resolve(
        card,
        &PropertyKey::new(
            "mean_engineering_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
        ),
        &query,
        reference_axis,
        endpoint_axis,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(thermal.reference_axis(), reference_axis);
    assert_eq!(thermal.endpoint_axis(), endpoint_axis);
    FiniteSolidMaterial::resolve_mean_thermal(
        card,
        0,
        thermal,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g0_mean_reference_changes_preserve_motion_inventory_and_signed_source_work() {
    with_cx(|_, cx| {
        let card = super::super::thermal::thermal_card(0.01, 1.0, false);
        let cold = ConservativeThermalField::uniform(
            mean_material(&card, 300.0, false),
            ReferenceDeformation::default(),
            1,
            cx,
        )
        .unwrap();
        let cold_material = cold.material_at(0).unwrap();
        let refreshed = cold_material.at_temperature(&card, 310.0).unwrap();
        let independently_prepared = mean_material(&card, 310.0, false);
        assert_eq!(refreshed.density(), cold_material.density());
        assert_eq!(refreshed.law(), independently_prepared.law());
        assert_eq!(
            refreshed.thermal_reference().unwrap().principal_stretches(),
            independently_prepared
                .thermal_reference()
                .unwrap()
                .principal_stretches()
        );
        for refused in [311.0, -1.0, f64::NAN] {
            assert!(cold_material.at_temperature(&card, refused).is_err());
        }
        let foreign = super::super::thermal::thermal_card(0.02, 1.0, false);
        assert!(cold_material.at_temperature(&foreign, 310.0).is_err());
        let nonthermal = FiniteSolidMaterial::resolve(
            &card,
            0,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(nonthermal.at_temperature(&card, 305.0).is_err());
        let hot =
            ConservativeThermalField::uniform(refreshed, ReferenceDeformation::default(), 1, cx)
                .unwrap();
        let endpoint_density = FiniteSolidMaterial::resolve(
            &card,
            0,
            &point(310.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .reference_density_kg_m3();
        let natural_volume = 1.1_f64.powi(3);
        near(endpoint_density, 6.0 / natural_volume, 1e-13);
        assert_ne!(hot.reference_density_kg_m3()[0], endpoint_density);

        let (nodes, elements) = mesh(cx);
        let speed = [0.25, -0.1, 0.05];
        let original = cold
            .prepare(
                &nodes,
                &elements,
                &vec![[0.0; 3]; nodes.len()],
                &vec![speed; nodes.len()],
                &[],
                preparation(ConservativeFormulation::FullMaterial),
                cx,
            )
            .unwrap();
        // Begin with an accepted moving trajectory at nonzero physical time.
        let moving = original
            .advance(0.02, &[], &[], [0.0; 3], settings(), 1e-8, cx)
            .unwrap()
            .solid;
        let old_u = moving.native().displacement_m().to_vec();
        let old_v = moving.native().velocity_mps().to_vec();
        let old_time = moving.native().time_s();
        let kinetic = moving.native().kinetic_energy_j();
        let heating = moving.change_reference(&hot, &[], 1e-8, cx).unwrap();
        // Authored Ogden alpha=2 has U(Je)=20*(Je-1)^2. A fixed F=I
        // carries no deviatoric energy, and storage is per original volume.
        let expected = natural_volume * 20.0 * (1.0 / natural_volume - 1.0).powi(2) / 6.0;
        near(heating.energy.reference_work_j, expected, 1e-10);
        near(heating.energy.stored_after_j, expected, 1e-10);
        near(heating.energy.cumulative_balance_error_j, 0.0, 1e-8);
        assert!(std::ptr::eq(heating.solid.source_field(), &hot));
        assert_eq!(heating.solid.native().displacement_m(), old_u);
        assert_eq!(heating.solid.native().velocity_mps(), old_v);
        assert_eq!(heating.solid.native().time_s(), old_time);
        near(heating.solid.native().kinetic_energy_j(), kinetic, 0.0);
        near(heating.solid.native().reference_mass_kg(), 1.0, 2e-12);
        assert_eq!(heating.solid.native().reference_density_kg_m3(), &[6.0]);

        let cooling = heating
            .solid
            .change_reference(&cold, &[], 1e-8, cx)
            .unwrap();
        near(cooling.energy.reference_work_j, -expected, 1e-10);
        near(cooling.energy.cumulative_balance_error_j, 0.0, 1e-8);
        let continued = cooling
            .solid
            .advance(0.01, &[], &[], [0.0; 3], settings(), 1e-8, cx)
            .unwrap();
        near(continued.solid.native().time_s(), 0.03, 1e-15);
        near(continued.energy.cumulative_balance_error_j, 0.0, 1e-8);
        near(continued.solid.native().kinetic_energy_j(), kinetic, 1e-10);
        assert_eq!(moving.native().displacement_m(), old_u);
        assert_eq!(moving.native().velocity_mps(), old_v);
        assert_eq!(moving.native().time_s(), old_time);

        // Both endpoints have the same card, datum, numerical G and law
        // context. Interchanging the source's named temperature roles must
        // still refuse: a numerical endpoint does not define its convention.
        let reversed = ConservativeThermalField::uniform(
            mean_material(&card, 310.0, true),
            ReferenceDeformation::default(),
            1,
            cx,
        )
        .unwrap();
        assert_eq!(reversed.reference_at(0), hot.reference_at(0));
        assert!(matches!(
            moving.change_reference(&reversed, &[], 1e-8, cx),
            Err(ConservativeThermalError::IncompatibleSource { .. })
        ));

        // The next mean coefficient is a fresh source evaluation. A second
        // claim becomes valid only at the new endpoint, so the original claim
        // must remain pinned instead of re-running a broad selection policy.
        let ThermalReference::Mean(mean) = cold_material.thermal_reference().unwrap() else {
            panic!("expected original mean source");
        };
        let selected = mean.source().resolved().properties()[0]
            .answer()
            .receipt
            .selected;
        let density_selected = cold_material.density().properties()[0]
            .answer()
            .receipt
            .selected;
        let mut varying = card.claims().claim(selected).unwrap().clone();
        let mut competing = varying.clone();
        varying.value = PropertyValue::Curve {
            abscissa: "T_endpoint".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 0.005), (310.0, 0.01)],
            dims: INVERSE_TEMPERATURE_DIMS,
        };
        varying.interpolation = InterpolationPolicy::LinearInside;
        competing.value = PropertyValue::Scalar {
            value: 0.02,
            dims: INVERSE_TEMPERATURE_DIMS,
        };
        let quantity = QuantitySpec::dimensional(Temperature::DIMS);
        competing.validity = ValidityDomain::unconstrained()
            .with_quantity("T_reference", quantity, 300.0, 310.0)
            .with_quantity("T_endpoint", quantity, 305.0, 310.0);
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(card.claims().claim(density_selected).unwrap().clone())
            .unwrap();
        claims.insert_claim(varying).unwrap();
        claims.insert_claim(competing).unwrap();
        let varying_card =
            MaterialCard::assemble(card.id().clone(), claims, card.models().to_vec()).unwrap();
        let varying_cold = mean_material(&varying_card, 300.0, false);
        let varying_hot = varying_cold.at_temperature(&varying_card, 310.0).unwrap();
        for stretch in varying_hot
            .thermal_reference()
            .unwrap()
            .principal_stretches()
        {
            near(stretch, 1.1, 1e-14);
        }
        assert_eq!(varying_hot.density(), varying_cold.density());
        let ThermalReference::Mean(before) = varying_cold.thermal_reference().unwrap() else {
            panic!("expected mean source");
        };
        let ThermalReference::Mean(after) = varying_hot.thermal_reference().unwrap() else {
            panic!("expected refreshed mean source");
        };
        near(
            before.source().resolved().properties()[0].value_si(),
            0.005,
            0.0,
        );
        near(
            after.source().resolved().properties()[0].value_si(),
            0.01,
            0.0,
        );
        assert_eq!(
            before.source().resolved().properties()[0]
                .answer()
                .receipt
                .selected,
            after.source().resolved().properties()[0]
                .answer()
                .receipt
                .selected
        );
    });
}

#[test]
fn g0_thermal_reference_events_refuse_changed_source_frame_datum_and_manufacturing_map() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let identity = ReferenceDeformation::default();
        let original =
            ConservativeThermalField::uniform(bind_expansion(&source, Q, 300.0), identity, 1, cx)
                .unwrap();
        let other_source = card([1.08, 0.97, 1.03], "T", false);
        let source_changed = ConservativeThermalField::uniform(
            bind_expansion(&other_source, Q, 300.0),
            identity,
            1,
            cx,
        )
        .unwrap();
        let frame_changed = ConservativeThermalField::uniform(
            bind_expansion(
                &source,
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                300.0,
            ),
            identity,
            1,
            cx,
        )
        .unwrap();
        let different_datum = ExpansionThermalReference::try_new(
            &source,
            path(&source, false, "T", point(305.0)),
            305.0,
            frame(Q),
        )
        .unwrap();
        let datum_changed = ConservativeThermalField::uniform(
            FiniteSolidMaterial::resolve_thermal(
                &source,
                0,
                ThermalReference::Expansion(different_datum),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
            identity,
            1,
            cx,
        )
        .unwrap();
        // At each own datum all three candidates have the same numerical
        // natural strain and density. Their physical source bindings differ.
        for candidate in [&source_changed, &frame_changed, &datum_changed] {
            assert_eq!(candidate.reference_at(0), original.reference_at(0));
            assert_eq!(candidate.reference_density_kg_m3(), &[6.0]);
        }
        let manufacturing_changed = ConservativeThermalField::uniform(
            bind_expansion(&source, Q, 300.0),
            ReferenceDeformation::try_new([1.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap(),
            1,
            cx,
        )
        .unwrap();
        let (nodes, elements) = mesh(cx);
        let body = original
            .prepare(
                &nodes,
                &elements,
                &vec![[0.0; 3]; nodes.len()],
                &vec![[0.0; 3]; nodes.len()],
                &[],
                preparation(ConservativeFormulation::FullMaterial),
                cx,
            )
            .unwrap();
        for candidate in [
            &source_changed,
            &frame_changed,
            &datum_changed,
            &manufacturing_changed,
        ] {
            assert!(matches!(
                body.change_reference(candidate, &[], 1e-8, cx),
                Err(ConservativeThermalError::IncompatibleSource { .. })
            ));
            assert!(matches!(
                body.change_reference_and_advance(
                    candidate,
                    &[],
                    0.01,
                    &[],
                    &[],
                    [0.0; 3],
                    settings(),
                    1e-8,
                    cx
                ),
                Err(ConservativeThermalError::IncompatibleSource { .. })
            ));
        }
        assert_eq!(body.native().displacement_m(), vec![[0.0; 3]; nodes.len()]);
        assert_eq!(body.native().velocity_mps(), vec![[0.0; 3]; nodes.len()]);
        assert_eq!(body.native().time_s(), 0.0);
        assert!(std::ptr::eq(body.source_field(), &original));
    });
}

#[test]
fn g1_source_expansion_and_contraction_supply_exact_constraint_work_to_native_motion() {
    with_cx(|_, cx| {
        let source = card([1.01; 3], "T", false);
        let frame = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let datum = ConservativeThermalField::uniform(
            bind_expansion(&source, frame, 300.0),
            ReferenceDeformation::default(),
            1,
            cx,
        )
        .unwrap();
        let (nodes, elements) = mesh(cx);
        let formulation = ConservativeFormulation::ExactIncompressible {
            maximum_relative_volume_residual: 1e-9,
            maximum_initial_relative_volume_rate_per_s: 1e-9,
            maximum_point_volume_work_defect: 1e-10,
        };
        let body = datum
            .prepare(
                &nodes,
                &elements,
                &vec![[0.0; 3]; nodes.len()],
                &vec![[0.0; 3]; nodes.len()],
                &[],
                preparation(formulation),
                cx,
            )
            .unwrap();
        for (temperature, stretch) in [(310.0, 1.01_f64), (290.0, 0.99_f64)] {
            let refreshed = datum
                .material_at(0)
                .unwrap()
                .at_temperature(&source, temperature)
                .unwrap();
            assert_eq!(refreshed, bind_expansion(&source, frame, temperature));
            let target = ConservativeThermalField::uniform(
                refreshed,
                ReferenceDeformation::default(),
                1,
                cx,
            )
            .unwrap();
            near(
                target.reference_at(0).unwrap().volume_ratio(),
                stretch.powi(3),
                1e-14,
            );
            // A fixed-F event cannot manufacture the new natural volume.
            assert!(body.change_reference(&target, &[], 1e-8, cx).is_err());
            let changed = body
                .change_reference_and_advance(
                    &target,
                    &[],
                    0.05,
                    &[],
                    &[],
                    [0.0; 3],
                    settings(),
                    1e-8,
                    cx,
                )
                .unwrap();
            let next = &changed.step.solid;
            // Each P1 pressure shape integrates to V/4=1/24 on this
            // original tetrahedron. The native multiplier is tensile
            // positive, so the source supplies -p_interval * Delta J_g.
            let pressure_sum: f64 = next.native().algorithmic_pressure_pa().iter().sum();
            let reference_work = -(stretch.powi(3) - 1.0) * pressure_sum / 24.0;
            assert!(reference_work > 0.0);
            near(changed.reference.reference_work_j, 0.0, 1e-11);
            near(
                changed.step.energy.constraint_reference_work_j,
                reference_work,
                1e-10,
            );
            near(next.native().total_energy_j(), reference_work, 1e-8);
            near(changed.step.energy.cumulative_balance_error_j, 0.0, 1e-8);
            near(next.native().reference_mass_kg(), 1.0, 2e-12);
            assert!(std::ptr::eq(next.source_field(), &target));
            let held = next
                .advance(0.01, &[], &[], [0.0; 3], settings(), 1e-8, cx)
                .unwrap();
            near(held.energy.constraint_reference_work_j, 0.0, 0.0);
            near(held.energy.cumulative_balance_error_j, 0.0, 1e-8);
            near(held.solid.native().total_energy_j(), reference_work, 1e-8);
            assert!(std::ptr::eq(held.solid.source_field(), &target));
            assert_eq!(body.native().displacement_m(), vec![[0.0; 3]; nodes.len()]);
            assert_eq!(body.native().velocity_mps(), vec![[0.0; 3]; nodes.len()]);
            assert_eq!(body.native().time_s(), 0.0);
        }
        assert!(
            datum
                .material_at(0)
                .unwrap()
                .at_temperature(&source, 311.0)
                .is_err()
        );
    });
}

#[test]
fn source_reference_field_retains_order_original_density_and_exact_element_coverage() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let hot = bind_expansion(&source, Q, 310.0);
        let warm = bind_expansion(&source, Q, 305.0);
        let manufacture =
            ReferenceDeformation::try_new([1.03, 0.16, 0.02, 0.0, 0.97, -0.01, 0.0, 0.0, 1.02])
                .unwrap();
        let thermal = match hot.thermal_reference().unwrap() {
            ThermalReference::Expansion(reference) => *reference.reference_deformation(),
            _ => panic!("expected authored directional source"),
        };
        let expected = manufacture.followed_by(thermal).unwrap();
        let reversed = thermal.followed_by(manufacture).unwrap();
        assert_ne!(
            expected.deformation_gradient(),
            reversed.deformation_gradient()
        );
        let uniform = ConservativeThermalField::uniform(hot.clone(), manufacture, 2, cx).unwrap();
        assert_eq!(uniform.element_count(), 2);
        assert_eq!(uniform.reference_at(0), Some(&expected));
        assert_eq!(uniform.reference_at(1), Some(&expected));
        assert_eq!(uniform.reference_before_thermal_at(1), Some(&manufacture));
        assert!(uniform.material_at(2).is_none());
        assert!(uniform.reference_at(2).is_none());
        assert!(uniform.reference_before_thermal_at(2).is_none());
        assert_eq!(uniform.material_at(0).unwrap(), &hot);
        assert_eq!(uniform.reference_density_kg_m3(), &[6.0, 6.0]);

        let references = [manufacture, ReferenceDeformation::default()];
        let field =
            ConservativeThermalField::per_element(vec![hot.clone(), warm.clone()], &references, cx)
                .unwrap();
        assert_eq!(field.material_at(0).unwrap(), &hot);
        assert_eq!(field.material_at(1).unwrap(), &warm);
        assert_eq!(field.reference_at(0), Some(&expected));
        assert_ne!(field.reference_at(0), field.reference_at(1));
        assert_eq!(field.reference_density_kg_m3(), &[6.0, 6.0]);
        assert!(ConservativeThermalField::uniform(hot.clone(), manufacture, 0, cx).is_err());
        assert!(ConservativeThermalField::per_element(vec![hot], &[], cx).is_err());
        let (nodes, elements) = mesh(cx);
        assert!(matches!(
            field.prepare(
                &nodes,
                &elements,
                &vec![[0.0; 3]; nodes.len()],
                &vec![[0.0; 3]; nodes.len()],
                &[],
                preparation(ConservativeFormulation::FullMaterial),
                cx
            ),
            Err(ConservativeThermalError::Invalid(_))
        ));
    });
}

#[test]
fn large_sourced_natural_expansion_moves_conservatively_with_original_mass_in_all_formulations() {
    with_cx(|_, cx| {
        // Existing authored source: 300->310 K has s=4 and rho0=6,
        // while endpoint density is 6/64. The original inventory must win.
        let source = super::super::thermal::thermal_card(0.3, 1.0, false);
        let material = super::super::thermal::prepared(&source);
        let reference = ReferenceDeformation::default();
        let field = ConservativeThermalField::uniform(material, reference, 1, cx).unwrap();
        let g = field.reference_at(0).unwrap();
        near(g.volume_ratio(), 64.0, 1e-12);
        assert!(
            field
                .material_at(0)
                .unwrap()
                .law()
                .law()
                .tangent(g.deformation_gradient())
                .is_err()
        );
        let (nodes, elements) = mesh(cx);
        let displacement: Vec<_> = nodes.iter().map(|x| x.map(|x| 3.0 * x)).collect();
        let speed = [0.25, -0.1, 0.05];
        let velocities = vec![speed; nodes.len()];
        let mass = 6.0 / 6.0;
        let kinetic = 0.5 * mass * speed.iter().map(|v| v * v).sum::<f64>();
        let exact = ConservativeFormulation::ExactIncompressible {
            maximum_relative_volume_residual: 1e-9,
            maximum_initial_relative_volume_rate_per_s: 1e-9,
            maximum_point_volume_work_defect: 1e-10,
        };
        let bulk = ConservativeFormulation::OriginalMaterial {
            bulk: DiscreteBulkSettings {
                minimum_volume_ratio: 0.5,
                maximum_volume_ratio: 1.5,
                maximum_iterations: 96,
                absolute_pressure_tolerance_pa: 1e-10,
                relative_pressure_tolerance: 1e-12,
                absolute_work_tolerance_pa: 1e-10,
                relative_work_tolerance: 1e-12,
            },
            maximum_relative_volume_residual: 1e-9,
            maximum_point_volume_work_defect: 1e-10,
        };
        for formulation in [ConservativeFormulation::FullMaterial, exact, bulk] {
            let mut body = field
                .prepare(
                    &nodes,
                    &elements,
                    &displacement,
                    &velocities,
                    &[],
                    preparation(formulation),
                    cx,
                )
                .unwrap();
            assert!(std::ptr::eq(body.source_field(), &field));
            near(body.native().reference_mass_kg(), mass, 2e-12);
            near(body.native().stored_energy_j(), 0.0, 1e-9);
            near(body.native().kinetic_energy_j(), kinetic, 1e-12);
            let acceleration = vec![[0.3, -0.2, 0.1]; nodes.len()];
            let load = body
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap();
            for axis in 0..3 {
                let total: f64 = load
                    .iter()
                    .filter(|(dof, _)| dof % 3 == axis)
                    .map(|(_, f)| f)
                    .sum();
                near(total, mass * acceleration[0][axis], 2e-12);
            }
            let mut time = 0.0;
            for dt in [0.01, 0.015] {
                let before_u = body.native().displacement_m().to_vec();
                let before_v = body.native().velocity_mps().to_vec();
                let before_t = body.native().time_s();
                let step = body
                    .advance(dt, &[], &[], [0.0; 3], settings(), 1e-8, cx)
                    .unwrap();
                assert_eq!(body.native().displacement_m(), before_u);
                assert_eq!(body.native().velocity_mps(), before_v);
                assert_eq!(body.native().time_s(), before_t);
                time += dt;
                for (node, (u, v)) in step
                    .solid
                    .native()
                    .displacement_m()
                    .iter()
                    .zip(step.solid.native().velocity_mps())
                    .enumerate()
                {
                    for axis in 0..3 {
                        near(u[axis], displacement[node][axis] + time * speed[axis], 2e-9);
                        near(v[axis], speed[axis], 2e-8);
                    }
                }
                near(step.inertia.numerical_loss_j, 0.0, 0.0);
                near(step.inertia.kinetic_after_j, kinetic, 1e-9);
                near(step.energy.balance_error_j, 0.0, 1e-8);
                near(step.energy.cumulative_balance_error_j, 0.0, 1e-8);
                assert!(std::ptr::eq(step.solid.source_field(), &field));
                assert_eq!(
                    step.solid
                        .reference_acceleration_load_n(&acceleration, cx)
                        .unwrap(),
                    load
                );
                near(step.solid.native().reference_mass_kg(), mass, 2e-12);
                body = step.solid;
            }
        }
    });
}
