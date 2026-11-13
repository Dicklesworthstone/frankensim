//! Synthetic spatial source/history regressions, not measured hardness accuracy.
use super::*;
use fs_couple::solid_enthalpy::plastic::hardness::ThermoplasticHardnessFieldError;
use fs_matdb::{HardnessLoadStep, HardnessTestContext, ObservationDataset};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_qty::{QtyAny, QuantitySpec};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;

fn response(scale: f64, temperature: f64, alpha: f64) -> f64 {
    scale * (1.0 - 0.001 * (temperature - 300.0) + 8.0 * alpha)
}

fn with_hardness(original: Fixture, scale: f64, maximum_alpha: f64) -> (Fixture, PropertyKey) {
    let mut claims = original.card.claims().clone();
    let provenance = Provenance {
        source: "synthetic independent spatial indentation test".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic original spatial material".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", &scale.to_le_bytes()),
            caveats: "numerical regression only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic indentation revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), ALPHA.into()],
                axis_dims: [Temperature::DIMS, Dims::NONE],
                samples: [
                    [280.0, 0.0],
                    [320.0, 0.0],
                    [280.0, maximum_alpha],
                    [320.0, maximum_alpha],
                ]
                .map(|point| (point, response(scale, point[0], point[1])))
                .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with_quantity(
                    ALPHA,
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                    maximum_alpha,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(original.card.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &QueryPoint::new().with("T", 280.0).unwrap(),
            upper: &QueryPoint::new().with("T", 320.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        original.caloric.curve().knots()[0].specific_enthalpy_j_kg,
        1e-5,
    )
    .unwrap();
    (
        Fixture {
            card,
            caloric,
            expansion: original.expansion,
        },
        key,
    )
}

fn pair_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>, Vec<(usize, f64)>) {
    let mut nodes = NODES.to_vec();
    nodes.extend(NODES.map(|[x, y, z]| [x + 3.0 * L, y, z]));
    let mut prescribed = ROLLERS.to_vec();
    prescribed.extend(ROLLERS.map(|(dof, value)| (dof + 12, value)));
    (nodes, vec![[0, 1, 2, 3], [4, 5, 6, 7]], prescribed)
}

#[test]
fn spatial_hardness_uses_each_original_temperature_and_plastic_history() {
    context(|cx| {
        let cp = [10.0, 100.0];
        let scale = [2e8, 4e8];
        let stress = [3.4e6, 3.8e6];
        let (a, key_a) = with_hardness(Fixture::new(cp[0]), scale[0], 0.04);
        let (b, key_b) = with_hardness(Fixture::new(cp[1]), scale[1], 0.04);
        let selections = [(key_a.clone(), SELECT), (key_b.clone(), SELECT)];
        let thermal = initial(&[&a, &b], &[], &[0.0; 2]);
        let (nodes, tets, prescribed) = pair_mesh();
        let solid = bind(&thermal, &[&a, &b], &nodes, &tets, cx);
        let binding = solid
            .bind_hardness_field(&thermal, &selections, 2, cx)
            .unwrap();
        let original = binding.at_accepted(&solid, &thermal, cx).unwrap();
        let native = solid
            .advance(
                &thermal,
                &(),
                0.001,
                &prescribed,
                &[(3, stress[0] * V / L), (15, stress[1] * V / L)],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let native_solution = native.solid().solution().clone();
        let native_report = format!("{:?}", native.report());
        let accepted = binding.check_candidate(native, cx).unwrap();
        assert_eq!(accepted.mechanics.solid().solution(), &native_solution);
        assert_eq!(format!("{:?}", accepted.mechanics.report()), native_report);
        assert_eq!(accepted.hardness.time_s(), 0.001);
        assert_eq!(accepted.hardness.points().len(), 2);
        for (i, point) in accepted.hardness.points().iter().enumerate() {
            let cell = CellId(100 + i as u64);
            let temperature = accepted.mechanics.thermal().temperature_k(cell).unwrap();
            let alpha = native_solution.element_states[i].alpha;
            assert_eq!((point.element, point.cell), (i, cell));
            assert_eq!(point.hardness.key(), &selections[i].0);
            assert_eq!(
                point.hardness.resolved().state_coordinate("T"),
                Some(temperature)
            );
            assert_eq!(
                point.hardness.resolved().state_coordinate(ALPHA),
                Some(alpha)
            );
            assert_eq!(
                point.hardness.property().answer().receipt.selected,
                original.points()[i]
                    .hardness
                    .property()
                    .answer()
                    .receipt
                    .selected,
            );
            close(
                point.hardness.hardness_pa(),
                response(scale[i], temperature, alpha),
                1e-5,
            );
            // Independent uniaxial plastic/heat balance, not a call to the source adapter.
            let b = RHO * cp[i] * H - SOFTENING * (2.0 * YIELD - stress[i]);
            let c = YIELD * (stress[i] - YIELD);
            let theta = 2.0 * c / (b + (b * b + 4.0 * SOFTENING.powi(2) * c).sqrt());
            let expected_alpha = (stress[i] - YIELD + SOFTENING * theta) / H;
            close(temperature, 300.0 + theta, 2e-7);
            close(alpha, expected_alpha, 2e-9);
            close(
                point.hardness.hardness_pa(),
                response(scale[i], 300.0 + theta, expected_alpha),
                8.0,
            );
        }
        assert_eq!(
            accepted.hardness.hardness_pa().collect::<Vec<_>>(),
            accepted
                .hardness
                .points()
                .iter()
                .map(|p| p.hardness.hardness_pa())
                .collect::<Vec<_>>(),
        );
        assert!(
            accepted.hardness.points()[1].hardness.hardness_pa()
                > accepted.hardness.points()[0].hardness.hardness_pa()
        );
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(solid.solution().element_states[0].alpha, 0.0);
    });
}

#[test]
fn spatial_hardness_retains_complete_cyclic_branches_and_original_work() {
    context(|cx| {
        let scale = 3e8;
        let (fixture, key) = with_hardness(chaboche::fixture(None, true), scale, 0.04);
        let mut thermal = initial(&[&fixture], &[], &[0.0]);
        let mut solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx)
            .with_voce_hardening(&[Some(SELECT)], cx)
            .unwrap()
            .with_chaboche_branches(&[(SELECT, chaboche::schema())], cx)
            .unwrap();
        let binding = solid
            .bind_hardness_field(&thermal, &[(key, SELECT)], 1, cx)
            .unwrap();
        let mut previous_alpha = 0.0;
        for (index, strain) in [0.012, -0.012].into_iter().enumerate() {
            let mut prescribed = ROLLERS.to_vec();
            prescribed.push((3, L * strain));
            let native = solid
                .advance(
                    &thermal,
                    &(),
                    (index + 1) as f64 * 0.001,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let complete_history = native.solid().solution().clone();
            let report = format!("{:?}", native.report());
            let accepted = binding.check_candidate(native, cx).unwrap();
            assert_eq!(accepted.mechanics.solid().solution(), &complete_history);
            assert_eq!(format!("{:?}", accepted.mechanics.report()), report);
            let history = &complete_history.element_states[0];
            assert!(history.plastic.alpha > previous_alpha);
            assert!(
                history
                    .backstress_pa
                    .iter()
                    .all(|branch| branch.iter().any(|x| *x != 0.0))
            );
            let temperature = accepted
                .mechanics
                .thermal()
                .temperature_k(CellId(100))
                .unwrap();
            close(
                accepted.hardness.points()[0].hardness.hardness_pa(),
                response(scale, temperature, history.plastic.alpha),
                1e-5,
            );
            assert_eq!(
                accepted.hardness.points()[0]
                    .hardness
                    .resolved()
                    .state_coordinate(ALPHA),
                Some(history.plastic.alpha),
            );
            previous_alpha = history.plastic.alpha;
            (thermal, solid, _, _) = accepted.mechanics.into_parts();
        }
    });
}

#[test]
fn spatial_hardness_late_coverage_refusal_preserves_accepted_owners() {
    context(|cx| {
        let (a, key_a) = with_hardness(Fixture::new(100.0), 2e8, 0.04);
        let (b, key_b) = with_hardness(Fixture::new(100.0), 4e8, 1e-5);
        let thermal = initial(&[&a, &b], &[], &[0.0; 2]);
        let (nodes, tets, prescribed) = pair_mesh();
        let solid = bind(&thermal, &[&a, &b], &nodes, &tets, cx);
        let selections = [(key_a, SELECT), (key_b, SELECT)];
        assert!(
            solid
                .bind_hardness_field(&thermal, &selections, 1, cx)
                .is_err()
        );
        let binding = solid
            .bind_hardness_field(&thermal, &selections, 2, cx)
            .unwrap();
        let old_field = binding.at_accepted(&solid, &thermal, cx).unwrap();
        let original_solution = solid.solution().clone();
        let foreign = bind(&thermal, &[&a, &b], &nodes, &tets, cx);
        assert!(matches!(
            binding.at_accepted(&foreign, &thermal, cx),
            Err(ThermoplasticHardnessFieldError::Invalid(_))
        ));
        let stale = thermal
            .advance_to(0.001, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert!(matches!(
            binding.at_accepted(&solid, &stale, cx),
            Err(ThermoplasticHardnessFieldError::Invalid(_))
        ));
        let native = solid
            .advance(
                &thermal,
                &(),
                0.001,
                &prescribed,
                &[(3, 3.4e6 * V / L), (15, 3.4e6 * V / L)],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert!(native.solid().solution().element_states[1].alpha > 1e-5);
        assert!(matches!(
            binding.check_candidate(native, cx),
            Err(ThermoplasticHardnessFieldError::Source {
                element: 1,
                cell: CellId(101),
                ..
            }),
        ));
        assert_eq!(solid.solution(), &original_solution);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(
            binding.at_accepted(&solid, &thermal, cx).unwrap(),
            old_field
        );
    });
}
