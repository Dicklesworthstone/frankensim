//! Original coupled heat/history drives H; synthetic numerical evidence only.
use super::*;
use fs_couple::solid_enthalpy::hex::hardness::HexThermoplasticHardnessError;
use fs_matdb::{HardnessLoadStep, HardnessTestContext, ObservationDataset};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::state_point::YIELD_STRESS_PROPERTY;
use fs_qty::QtyAny;

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;

fn response(temperature: f64, alpha: f64) -> f64 {
    2e8 * (1.0 - 0.01 * (temperature - 300.0) + 10.0 * alpha)
}

fn fixture(maximum_alpha: f64) -> (Fixture, PropertyKey) {
    let original = Fixture::new();
    let mut claims = original.card.claims().clone();
    claims
        .insert_claim(property(
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            Y,
            -SOFTENING,
        ))
        .unwrap();
    claims
        .insert_claim(property(
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            H,
            0.0,
        ))
        .unwrap();
    let provenance = Provenance {
        source: "synthetic original coupled Hex8 indentation".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic coupled Hex8 material".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"coupled Hex8 hardness"),
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
                    [308.0, 0.0],
                    [280.0, maximum_alpha],
                    [308.0, maximum_alpha],
                ]
                .map(|p| (p, response(p[0], p[1])))
                .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_spec(), 280.0, 308.0)
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
            heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(308.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-8,
    )
    .unwrap();
    (Fixture { card, caloric }, key)
}

fn bind_j2<'a>(
    fixture: &'a Fixture,
    heat: &DepositionSimulation<'_>,
    formulation: Formulation,
    cx: &Cx<'_>,
) -> ThermoplasticHexSolid<'a, ThermoplasticSource<'a>> {
    let source = ThermoplasticSource::try_new(
        SourceQuery {
            card: &fixture.card,
            caloric: &fixture.caloric,
            context: point(300.0),
            temperature_axis: "T".into(),
            selection: SELECT,
            hardening_basis: HardeningPropertyBasis::PlasticStrainModulus,
            thermal_expansion: None,
            strain_limit: 0.04,
        },
        0.0,
    )
    .unwrap();
    ThermoplasticHexSolid::try_new(
        heat,
        &NODES,
        &HEX,
        vec![source],
        &(0..8).map(cell).collect::<Vec<_>>(),
        discretization(formulation),
        cx,
    )
    .unwrap()
}

fn imposed(scale: f64) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .flat_map(|x| {
            [
                scale * x[0] * (1.0 + 0.5 * x[1] / L + 0.3 * x[2] / L),
                0.0,
                0.0,
            ]
        })
        .enumerate()
        .collect()
}

#[test]
fn accepted_hardness_follows_distinct_coupled_point_heat_and_reversal_without_changing_accounts() {
    context(|cx| {
        let (fixture, key) = fixture(0.04);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let mut heat = thermal(&fixture, &[V / 8.0; 8], 0.5, 0.002);
            let mut solid = bind_j2(&fixture, &heat, formulation, cx);
            let mut observer = Probe {
                steps: 0,
                fail: false,
            };
            let binding = solid
                .bind_hardness_field(&heat, &[(key.clone(), SELECT)], cx)
                .unwrap();
            let initial = binding.at_accepted(&solid, &heat, cx).unwrap();
            for (interval, scale) in [0.006, -0.006].into_iter().enumerate() {
                let target = (interval + 1) as f64 * DT;
                let candidate = solid
                    .advance(
                        &heat,
                        &observer,
                        target,
                        &imposed(scale),
                        &[],
                        settings(),
                        tolerances(),
                        thermal_config(),
                        config(),
                        CouplingMethod::RelaxedPicard,
                        cx,
                    )
                    .unwrap();
                let solution = candidate.solid().solution().clone();
                let cells = candidate.thermal().cells().clone();
                let report = format!("{:?}", candidate.report());
                let energy = candidate.quadrature_energy().to_vec();
                let observer_steps = candidate.observer().steps;
                let accepted = binding.check_candidate(candidate, cx).unwrap();
                assert_eq!(accepted.mechanics.solid().solution(), &solution);
                assert_eq!(accepted.mechanics.thermal().cells(), &cells);
                assert_eq!(format!("{:?}", accepted.mechanics.report()), report);
                assert_eq!(accepted.mechanics.quadrature_energy(), energy.as_slice());
                assert_eq!(accepted.mechanics.observer().steps, observer_steps);
                assert_eq!(accepted.hardness.time_s(), target);
                assert_eq!(accepted.hardness.points().len(), 8);
                for (q, point) in accepted.hardness.points().iter().enumerate() {
                    let temperature = accepted.mechanics.thermal().temperature_k(cell(q)).unwrap();
                    let alpha = solution.quadrature_states[q].alpha;
                    assert!(temperature > 300.0 + 0.5 * (interval + 1) as f64);
                    assert!(alpha > solid.solution().quadrature_states[q].alpha);
                    assert_eq!(
                        (point.element, point.local_point, point.cell),
                        (0, q, cell(q))
                    );
                    assert_eq!(point.hardness.key(), &key);
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
                        initial.points()[q]
                            .hardness
                            .property()
                            .answer()
                            .receipt
                            .selected
                    );
                    near(
                        point.hardness.hardness_pa(),
                        response(temperature, alpha),
                        1e-5,
                    );
                }
                assert!(solution.quadrature_states[7].alpha > solution.quadrature_states[0].alpha);
                assert_ne!(
                    accepted.hardness.points()[0].hardness.hardness_pa(),
                    accepted.hardness.points()[7].hardness.hardness_pa()
                );
                let report = accepted.mechanics.report();
                near(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                near(
                    report.thermal_audit.stored_energy_change_j,
                    report.external_heat_j
                        + report.material_energy.material.physical_heat_j
                        + report.material_energy.material.reversible_heat_j,
                    3e-9,
                );
                assert_eq!(observer.steps, 2 * interval);
                (heat, solid, observer, _, _) = accepted.mechanics.into_parts();
            }
        }
    });
}

#[test]
fn late_point_coverage_and_original_inventory_identity_gate_the_whole_candidate() {
    context(|cx| {
        let (wide, _) = fixture(0.04);
        let wide_heat = thermal(&wide, &[V / 8.0; 8], 0.5, 0.002);
        let wide_solid = bind_j2(&wide, &wide_heat, Formulation::BBar, cx);
        let wide_candidate = wide_solid
            .advance(
                &wide_heat,
                &(),
                DT,
                &imposed(0.006),
                &[],
                settings(),
                tolerances(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                cx,
            )
            .unwrap();
        // Set a source boundary strictly between the last native point and all
        // earlier points. The independent mechanical/card roles are unchanged.
        let history = &wide_candidate.solid().solution().quadrature_states;
        let earlier = history[..7]
            .iter()
            .map(|state| state.alpha)
            .fold(0.0_f64, f64::max);
        assert!(history[7].alpha > earlier);
        let (fixture, key) = fixture(0.5 * (earlier + history[7].alpha));
        let heat = thermal(&fixture, &[V / 8.0; 8], 0.5, 0.002);
        let solid = bind_j2(&fixture, &heat, Formulation::BBar, cx);
        let binding = solid
            .bind_hardness_field(&heat, &[(key, SELECT)], cx)
            .unwrap();
        let original_field = binding.at_accepted(&solid, &heat, cx).unwrap();
        let original_solution = solid.solution().clone();
        let original_cells = heat.cells().clone();

        let foreign = bind_j2(&fixture, &heat, Formulation::BBar, cx);
        assert!(matches!(
            binding.at_accepted(&foreign, &heat, cx),
            Err(HexThermoplasticHardnessError::Invalid(_))
        ));
        let mut volumes = [V / 8.0; 8];
        volumes[7] *= 1.1;
        let changed_inventory = thermal(&fixture, &volumes, 0.5, 0.002);
        assert!(matches!(
            binding.at_accepted(&solid, &changed_inventory, cx),
            Err(HexThermoplasticHardnessError::Native(_))
        ));
        let stale = heat.advance_to(DT, thermal_config(), || Ok(())).unwrap().0;
        assert!(matches!(
            binding.at_accepted(&solid, &stale, cx),
            Err(HexThermoplasticHardnessError::Invalid(_))
        ));

        let candidate = solid
            .advance(
                &heat,
                &(),
                DT,
                &imposed(0.006),
                &[],
                settings(),
                tolerances(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                cx,
            )
            .unwrap();
        assert!(matches!(binding.check_candidate(candidate, cx),
            Err(HexThermoplasticHardnessError::Source { point: 7, cell: id, .. }) if id == cell(7)));
        assert_eq!(solid.solution(), &original_solution);
        assert_eq!(heat.cells(), &original_cells);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(
            binding.at_accepted(&solid, &heat, cx).unwrap(),
            original_field
        );
    });
}
