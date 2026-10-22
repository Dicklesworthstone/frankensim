//! Pressure-conditioned source conductivity reaches the existing thermal
//! transport and lumped enthalpy owner without losing source-cell extrema.

mod support;

use fs_conduction::lumped::{
    BiotGate, LUMPED_BIOT_CEILING, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig,
    LumpedThermalTransport, MAX_LUMPED_THERMAL_TRANSPORT_SAMPLES, solve_lumped_enthalpy,
};
use fs_conduction::material::{CONDUCTIVITY_DIMS, ConductivityTable};
use fs_conduction::{ConductionError, TEMPERATURE_DIMS};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PINNED_CLAIM_POLICY_TAG,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use support::with_cx;

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal_conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    )
}

fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature(), t)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            pressure,
        )
        .unwrap()
}

fn source() -> PropertyClaim {
    PropertyClaim {
        key: key(),
        value: PropertyValue::Surface2D {
            axes: ["temperature".into(), "pressure".into()],
            axis_dims: [TEMPERATURE_DIMS, fs_qty::Pressure::DIMS],
            // A four-cell fan with a low-conductivity interior vertex.
            // The enclosing rectangle's corners alone cannot bound its minimum.
            samples: vec![
                ([300.0, 100_000.0], 35.0),
                ([800.0, 100_000.0], 35.0),
                ([300.0, 300_000.0], 35.0),
                ([800.0, 300_000.0], 35.0),
                ([550.0, 200_000.0], 0.05),
            ],
            triangles: vec![[0, 1, 4], [1, 3, 4], [3, 2, 4], [2, 0, 4]],
            dims: CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", temperature(), 300.0, 800.0)
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                100_000.0,
                300_000.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "manufactured conductivity(T, p) fixture; not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn card(source: PropertyClaim) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(source).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured pressure-conditioned conductor".into(),
            phase: "solid".into(),
            process: "numerical regression fixture".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn table(claims: &ClaimSet, points: &[QueryPoint]) -> Result<ConductivityTable, ConductionError> {
    ConductivityTable::from_claims_at_query_points(
        claims,
        &key(),
        "temperature",
        points,
        SelectionPolicy::SingleClaimOnly,
    )
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g3_surface_slice_retains_every_cell_edge_and_pins_added_samples() {
    let card = card(source());
    let coarse = table(
        card.claims(),
        &[point(300.0, 150_000.0), point(800.0, 150_000.0)],
    )
    .unwrap();
    assert_eq!(
        coarse.knots().iter().map(|(t, _)| *t).collect::<Vec<_>>(),
        vec![300.0, 425.0, 675.0, 800.0],
        "both crossed edges survive even when they have equal conductivity"
    );
    for ((_, actual), expected) in coarse.knots().iter().zip([35.0, 17.525, 17.525, 35.0]) {
        assert_close(*actual, expected);
    }
    let explicit = table(
        card.claims(),
        &[300.0, 425.0, 675.0, 800.0].map(|t| point(t, 150_000.0)),
    )
    .unwrap();
    assert_eq!(coarse, explicit);
    let refined = table(
        card.claims(),
        &[300.0, 400.0, 550.0, 700.0, 800.0].map(|t| point(t, 150_000.0)),
    )
    .unwrap();
    let selected = coarse.receipts()[0].selected;
    for ((t, value), receipt) in refined.knots().iter().zip(refined.receipts()) {
        card.claims().verify_receipt(receipt).unwrap();
        assert_eq!(receipt.selected, selected);
        assert_eq!(
            &receipt.axis_quantities,
            point(*t, 150_000.0).axis_quantities()
        );
        assert!(
            receipt
                .query_point
                .contains(&("pressure".into(), 150_000.0))
        );
        let canonical = card
            .claims()
            .query_pinned_typed(&key(), &point(*t, 150_000.0), selected)
            .unwrap();
        assert_close(*value, canonical.evidence.value.value);
        if [400.0, 550.0, 700.0].contains(t) {
            assert_eq!(receipt.policy, PINNED_CLAIM_POLICY_TAG);
        }
    }
    for t in [350.0, 500.0, 750.0] {
        let canonical = card
            .claims()
            .query_pinned_typed(&key(), &point(t, 150_000.0), selected)
            .unwrap();
        assert_close(coarse.eval(t).unwrap(), canonical.evidence.value.value);
        assert_close(refined.eval(t).unwrap(), coarse.eval(t).unwrap());
    }
}

#[test]
fn g1_pressure_conditioned_transport_drives_enthalpy_and_gates_hidden_minimum() {
    let card = card(source());
    let transport_at = |pressure| {
        LumpedThermalTransport::from_material_card_without_radiation(
            &card,
            &key(),
            "temperature",
            &[point(300.0, pressure), point(800.0, pressure)],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap()
    };
    // cp = 100 J/(kg K), with density, mass and area consistent with Lc = 1 mm.
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [(0.0, 300.0), (50_000.0, 800.0)]
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 5_000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let body_at = |pressure| {
        LumpedEnthalpyBody::try_new_with_transport(
            "pressure-conditioned specimen",
            0.1,
            0.02,
            20.0,
            0.001,
            transport_at(pressure),
            &curve,
        )
        .unwrap()
    };
    let admitted = body_at(150_000.0);
    let nonuniform = body_at(200_000.0);
    assert_ne!(
        admitted.transport().identity(),
        nonuniform.transport().identity()
    );
    assert_eq!(
        admitted.transport().material_card_identity(),
        Some(card.content_hash())
    );
    assert_close(admitted.maximum_biot(300.0), 0.02 / 17.525);
    assert_close(nonuniform.maximum_biot(300.0), 0.02 / 0.05);
    assert!(nonuniform.maximum_biot(300.0) > LUMPED_BIOT_CEILING);
    let config = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 20_000.0,
        ambient_temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        internal_power_w: 0.0,
        duration_s: 10.0,
        maximum_step_s: 1.0,
        maximum_steps: 10,
        enthalpy_tolerance_j_kg: 1e-8,
    };
    let march = with_cx(|cx| {
        assert!(
            matches!(solve_lumped_enthalpy(cx, &nonuniform, BiotGate::corpus_default(), config),
                Err(ConductionError::ScenarioRow { what, .. }) if what.contains("Biot")),
            "a caller endpoint grid must not bypass the interior Biot limit"
        );
        solve_lumped_enthalpy(cx, &admitted, BiotGate::corpus_default(), config).unwrap()
    });
    // Backward Euler for C dT/dt = h A (Tamb - T), C = m cp = 10 J/K.
    let expected = 300.0 + 200.0 / 1.04_f64.powi(10);
    let last = march.samples().last().unwrap();
    assert!((last.phase_state.temperature_k() - expected).abs() < 1e-7);
    assert!(march.cumulative_absolute_energy_residual_j() < 1e-6);
    for pair in march.samples().windows(2) {
        let end = &pair[1];
        assert_eq!(end.radiation_into_body_w, 0.0);
        assert_close(
            end.convection_into_body_w,
            0.4 * (300.0 - end.phase_state.temperature_k()),
        );
        let stored = admitted.mass_kg()
            * (end.phase_state.specific_enthalpy_j_kg()
                - pair[0].phase_state.specific_enthalpy_j_kg());
        assert!((stored - (end.time_s - pair[0].time_s) * end.net_power_into_body_w).abs() < 1e-7);
    }
}

#[test]
fn g1_pressure_data_changes_spatial_temperature_and_heat_flux() {
    use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
    use fs_conduction::field::ScalarField;
    use fs_conduction::fixtures::{on_box_face, unit_cube};
    use fs_conduction::material::ConductivityModel;
    use fs_conduction::mesh::ConductionMesh;
    use fs_conduction::solve::{
        ConductionProblem, InitialGuess, Nonlinearity, SolveConfig, element_heat_flux, solve,
    };

    let card = card(source());
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold plate",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(450.0).unwrap(),
        )
        .unwrap()
        .region(
            "hot fluid",
            |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::robin(20.0, 650.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let source = ScalarField::Uniform(0.0);
    let mut responses = Vec::new();
    for (pressure, conductivity) in [(100_000.0, 35.0), (150_000.0, 17.525)] {
        let material = ConductivityModel::isotropic(
            table(
                card.claims(),
                &[point(300.0, pressure), point(800.0, pressure)],
            )
            .unwrap(),
        );
        let solution = with_cx(|cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                SolveConfig {
                    nonlinearity: Nonlinearity::FixedPoint {
                        relaxation: 1.0,
                        max_backtracks: 8,
                    },
                    initial: InitialGuess::Uniform(550.0),
                    ..SolveConfig::default()
                },
            )
            .unwrap()
        });
        // The whole solution stays in the source's 425..675 K plateau at
        // 150 kPa. Fourier conduction in series with the Robin film gives
        // T(x) = 450 + slope*x, slope = (650-450)/(L + k/h), L = 1 m.
        let slope = 200.0 / (1.0 + conductivity / 20.0);
        let flux = -conductivity * slope;
        assert!(solution.report.free_dofs > 0);
        assert_eq!(solution.report.material_receipts, material.receipts().len());
        for (position, temperature) in mesh.positions().iter().zip(&solution.temperature) {
            assert!((temperature - (450.0 + slope * position[0])).abs() < 1e-8);
        }
        for q in element_heat_flux(&mesh, &material, &solution.temperature).unwrap() {
            assert!((q[0] - flux).abs() < 1e-7);
            assert!(q[1].abs() < 1e-7 && q[2].abs() < 1e-7);
        }
        assert!((solution.report.energy.robin_out_w - flux).abs() < 1e-7);
        assert!((solution.report.energy.dirichlet_in_w - flux).abs() < 1e-7);
        assert!(solution.report.energy.closure_w.abs() < 1e-7);
        let fluid_trace = &solution.report.robin_fluxes[0];
        responses.push((fluid_trace.mean_wall_temperature_k, fluid_trace.heat_rate_w));
    }
    assert!(responses[0].0 < responses[1].0);
    assert!(responses[0].1 < responses[1].1);
}

#[test]
fn g0_surface_gaps_discrete_support_and_source_handoffs_refuse() {
    let endpoints = [point(300.0, 200_000.0), point(800.0, 200_000.0)];
    let mut gap = source();
    if let PropertyValue::Surface2D {
        samples, triangles, ..
    } = &mut gap.value
    {
        *samples = [
            [300.0, 100_000.0],
            [400.0, 100_000.0],
            [300.0, 300_000.0],
            [700.0, 100_000.0],
            [800.0, 100_000.0],
            [800.0, 300_000.0],
        ]
        .into_iter()
        .map(|at| (at, 35.0))
        .collect();
        *triangles = vec![[0, 1, 2], [3, 4, 5]];
    }
    let gap = card(gap);
    for point in &endpoints {
        gap.claims()
            .query_typed(&key(), point, SelectionPolicy::SingleClaimOnly)
            .unwrap();
    }
    assert!(table(gap.claims(), &endpoints).is_err());

    let mut discrete = source();
    discrete.interpolation = InterpolationPolicy::TabulatedOnly;
    if let PropertyValue::Surface2D { triangles, .. } = &mut discrete.value {
        triangles.clear();
    }
    let discrete = card(discrete);
    let endpoints = [point(300.0, 100_000.0), point(800.0, 100_000.0)];
    for point in &endpoints {
        discrete
            .claims()
            .query_typed(&key(), point, SelectionPolicy::SingleClaimOnly)
            .unwrap();
    }
    assert!(table(discrete.claims(), &endpoints).is_err());

    let mut claims = ClaimSet::new();
    let selected = claims.insert_claim(source()).unwrap();
    let mut competitor = source();
    competitor.value = PropertyValue::Scalar {
        value: 50.0,
        dims: CONDUCTIVITY_DIMS,
    };
    competitor.interpolation = InterpolationPolicy::ConstantWithinValidity;
    competitor.validity =
        competitor
            .validity
            .with_quantity("temperature", temperature(), 400.0, 450.0);
    claims.insert_claim(competitor).unwrap();
    let endpoints = [point(300.0, 150_000.0), point(800.0, 150_000.0)];
    assert!(
        table(&claims, &endpoints).is_err(),
        "hidden source competition cannot become one constitutive law"
    );
    ConductivityTable::from_claims_pinned_at_query_points(
        &claims,
        &key(),
        "temperature",
        &endpoints,
        selected,
    )
    .unwrap();
}

#[test]
fn g0_surface_transport_keeps_quantity_context_positivity_and_sample_limits() {
    let card = card(source());
    let endpoints = [point(300.0, 150_000.0), point(800.0, 150_000.0)];
    assert!(
        ConductivityTable::from_claims_at_query_points(
            card.claims(),
            &PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
            "temperature",
            &endpoints,
            SelectionPolicy::SingleClaimOnly,
        )
        .is_err(),
        "dimensional equality cannot erase the source quantity kind"
    );
    let difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    let wrong_temperature = endpoints.clone().map(|p| {
        let t = p.axes()["temperature"];
        p.with_quantity("temperature", difference, t).unwrap()
    });
    assert!(table(card.claims(), &wrong_temperature).is_err());
    assert!(
        table(
            card.claims(),
            &[point(300.0, 150_000.0), point(800.0, 200_000.0)]
        )
        .is_err()
    );
    let oversized = vec![point(300.0, 150_000.0); MAX_LUMPED_THERMAL_TRANSPORT_SAMPLES + 1];
    assert!(
        LumpedThermalTransport::from_material_card_without_radiation(
            &card,
            &key(),
            "temperature",
            &oversized,
            SelectionPolicy::SingleClaimOnly,
        )
        .is_err()
    );
    let mut nonphysical = source();
    if let PropertyValue::Surface2D { samples, .. } = &mut nonphysical.value {
        samples[4].1 = 0.0;
    }
    let nonphysical = self::card(nonphysical);
    assert!(
        matches!(
            table(
                nonphysical.claims(),
                &[point(300.0, 200_000.0), point(800.0, 200_000.0)]
            ),
            Err(ConductionError::Conductivity { .. })
        ),
        "an interior zero cannot hide behind positive endpoint samples"
    );
}

#[test]
fn g3_legacy_temperature_curve_keeps_source_knots_and_axis_dimensions() {
    let mut source = source();
    source.key = PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS);
    source.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: TEMPERATURE_DIMS,
        knots: vec![(300.0, 35.0), (550.0, 0.05), (800.0, 35.0)],
        dims: CONDUCTIVITY_DIMS,
    };
    source.validity = ValidityDomain::unconstrained().with("T", 300.0, 800.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(source.clone()).unwrap();
    let table = ConductivityTable::from_claims(
        &claims,
        "thermal_conductivity",
        &[300.0, 800.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        table.knots(),
        &[(300.0, 35.0), (550.0, 0.05), (800.0, 35.0)]
    );
    for receipt in table.receipts() {
        claims.verify_receipt(receipt).unwrap();
    }
    if let PropertyValue::Curve { abscissa_dims, .. } = &mut source.value {
        *abscissa_dims = fs_qty::Pressure::DIMS;
    }
    let mut wrong_axis = ClaimSet::new();
    wrong_axis.insert_claim(source).unwrap();
    assert!(
        matches!(
            ConductivityTable::from_claims(
                &wrong_axis,
                "thermal_conductivity",
                &[300.0, 800.0],
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(ConductionError::Conductivity { .. })
        ),
        "the name T cannot turn pressure coordinates into temperature"
    );
}
