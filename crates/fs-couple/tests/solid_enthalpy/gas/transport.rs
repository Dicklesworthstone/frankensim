//! Nonlinear internal material transport uses the same accepted gas/solid state.
use super::*;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};

#[path = "transport/heat.rs"]
mod heat;

fn nonlinear_card(slope: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = transport_material(1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, slope * 280.0), (600.0, slope * 600.0)],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = caloric_for(&card);
    (card, caloric)
}

fn thermal_pair<'a>(
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    powers: [f64; 2],
) -> DepositionSimulation<'a> {
    let events = [30, 10]
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 77,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 10.0 * ([340.0, 360.0][i] - 280.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64],
                    end_m: [1.0, 0.0, i as f64],
                    volume_m3: 1.0 / 6.0,
                    surface_area_m2: 3.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: powers[i],
                },
            }),
        })
        .chain(std::iter::once(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(30),
                b: CellId(10),
                area_m2: 0.5,
                // Deliberately distinct from the material value: the accepted source
                // callback must override this coefficient in the ORIGINAL residual.
                conductance_w_k: 0.5,
            }),
        }))
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0] as &dyn ThermalLaw),
            (MaterialId(1), curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn material_cells<'a>(
    cards: [&'a MaterialCard; 2],
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    [30, 10].map(|id| {
        let i = usize::from(id == 10);
        MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(i as u64),
            curve: curves[i],
            point: ThermalMaterialPoint {
                material: cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
        }
    })
}

const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.3,
    distance_b_m: 0.7,
}];

// Manufacture source powers from three independently written endpoint energy
// equations. The expected temperature, source coefficients and displacement do
// not come from a second invocation of the production transport resolver.
fn manufactured_heat() -> [f64; 3] {
    let [a, b, g] = [400.0, 440.0, 410.0];
    let conductance = 0.5 / (0.3 / (0.005 * a) + 0.7 / (0.01 * b));
    let internal = 2.0 * conductance * (a - b);
    let ga = 0.1 / (0.001 / 0.03 + 0.02 / (0.005 * a));
    let gb = 0.2 / (0.001 / 0.03 + 0.02 / (0.01 * b));
    let qa = 2.0 * ga * (a - g);
    let qb = 2.0 * gb * (b - g);
    let gas_capacity = 2.5 * 100000.0 * 0.1 / 390.0;
    [
        10.0 * (a - 340.0) + internal + qa,
        10.0 * (b - 360.0) - internal + qb,
        gas_capacity * (g - 390.0) - qa - qb,
    ]
}

#[test]
fn nonlinear_heterogeneous_internal_transport_closes_independent_energy_and_deformation() {
    with_cx(|cx| {
        let (a, ca) = nonlinear_card(0.005);
        let (b, cb) = nonlinear_card(0.01);
        let curves = [ca.curve(), cb.curve()];
        let heat = manufactured_heat();
        let thermal = thermal_pair(curves, [heat[0] / 2.0, heat[1] / 2.0]);
        let context = point(340.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::super::source(30, &a, ca.curve(), &context),
                super::super::source(10, &b, cb.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = material_cells([&a, &b], curves, &context, &key);
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &LINKS,
        };
        let transport = network.prepare(&thermal, || false).unwrap();
        let gas = gas(390.0, 0.1);
        let patches = contacts([0.1, 0.2]);
        let next = source
            .step_with_isochoric_gas_transport(
                &thermal,
                &gas,
                2.0,
                heat[2],
                &patches,
                &transport,
                solve(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            400.0,
            2e-7,
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            440.0,
            2e-7,
        );
        near(next.gas().state().state().temperature, 410.0, 2e-7);
        let sample = next.report().transport.as_ref().unwrap();
        assert_eq!(sample.links.len(), 1);
        assert_eq!(sample.links[0].cells, [CellId(30), CellId(10)]);
        near(sample.links[0].conductivity_w_m_k[0], 2.0, 1e-8);
        near(sample.links[0].conductivity_w_m_k[1], 4.4, 1e-8);
        near(
            sample.links[0].conductance_w_k,
            0.5 / (0.3 / 2.0 + 0.7 / 4.4),
            1e-8,
        );
        for patch in &next.report().contacts {
            let internal = &sample
                .cells
                .iter()
                .find(|cell| cell.cell == patch.cell)
                .unwrap()
                .conductivity
                .as_ref()
                .unwrap();
            assert_eq!(internal.query(), patch.source.query());
            assert_eq!(
                internal
                    .property_by_key(&key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                patch
                    .source
                    .property_by_key(&key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        let stored = next
            .thermal()
            .cells()
            .iter()
            .map(|(id, cell)| {
                cell.mass_kg * (cell.specific_h_j_kg - thermal.cells()[id].specific_h_j_kg)
            })
            .sum::<f64>()
            + gas.mass_kg()
                * (next.gas().state().specific_internal_energy_j_kg()
                    - gas.state().specific_internal_energy_j_kg());
        near(stored, heat.iter().sum(), 2e-7);
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        let equilibrium = next.state().solve(&loads, &FREE, config(), cx).unwrap();
        let u = equilibrium.solution.displacement.displacement_m();
        near(u[1][0], 0.0007, 1e-10);
        near(u[2][1], 0.0007, 1e-10);
        near(u[3][2], 0.0006, 1e-10);
        near(u[4][2], -0.0008, 1e-10);
        let frozen = source
            .step_with_isochoric_gas(
                &thermal,
                &gas,
                2.0,
                heat[2],
                &patches,
                solve(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        assert!(frozen.report().transport.is_none());
        assert!((frozen.thermal().temperature_k(CellId(30)).unwrap() - 400.0).abs() > 0.1);
        // A second interval reuses the SAME prepared claims at fresh states.
        let again = source
            .step_with_isochoric_gas_transport(
                next.thermal(),
                next.gas(),
                2.1,
                0.0,
                &patches,
                &transport,
                solve(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        assert_ne!(
            again.report().transport.as_ref().unwrap().links[0].conductance_w_k,
            sample.links[0].conductance_w_k
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas.state().state().temperature, 390.0);
    });
}

#[test]
fn changed_material_curve_changes_the_same_finite_gas_and_spatial_motion() {
    with_cx(|cx| {
        let heat = manufactured_heat();
        let mut results = Vec::new();
        for slope in [0.01, 0.02] {
            let (a, ca) = nonlinear_card(0.005);
            let (b, cb) = nonlinear_card(slope);
            let context = point(340.0);
            let curves = [ca.curve(), cb.curve()];
            let thermal = thermal_pair(curves, [heat[0] / 2.0, heat[1] / 2.0]);
            let source = SolidEnthalpy::try_new(
                &thermal,
                &PAIR_NODES,
                &PAIR_TETS,
                vec![
                    super::super::source(30, &a, ca.curve(), &context),
                    super::super::source(10, &b, cb.curve(), &context),
                ],
                TetAssemblyBudget::standard(),
                0.1,
                config().mesh_frame,
                cx,
            )
            .unwrap();
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let cells = material_cells([&a, &b], curves, &context, &key);
            let network = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &LINKS,
            };
            let prepared = network.prepare(&thermal, || false).unwrap();
            let next = source
                .step_with_isochoric_gas_transport(
                    &thermal,
                    &gas(390.0, 0.1),
                    2.0,
                    heat[2],
                    &contacts([0.1, 0.2]),
                    &prepared,
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let equilibrium = next.state().solve(&loads, &FREE, config(), cx).unwrap();
            results.push((
                next.gas().state().state().temperature,
                equilibrium.solution.displacement.displacement_m()[4][2],
            ));
            assert!(next.report().combined_energy_residual_j.abs() < 2e-7);
        }
        assert!((results[0].0 - results[1].0).abs() > 0.01);
        assert!((results[0].1 - results[1].1).abs() > 1e-7);
    });
}

#[test]
fn mismatched_source_context_and_failed_nonlinear_trial_preserve_both_inventories() {
    with_cx(|cx| {
        let (a, ca) = nonlinear_card(0.005);
        let (b, cb) = nonlinear_card(0.01);
        let curves = [ca.curve(), cb.curve()];
        let thermal = thermal_pair(curves, [0.0, 0.0]);
        let context = point(340.0);
        let other = context
            .clone()
            .with("different-fixed-condition", 1.0)
            .unwrap();
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::super::source(30, &a, ca.curve(), &context),
                super::super::source(10, &b, cb.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let gas = gas(390.0, 0.1);
        let saved = gas.clone();
        let old = thermal.cells().clone();
        let bad_cells = material_cells([&a, &b], curves, &other, &key);
        let bad_network = MaterialNetwork {
            cells: &bad_cells,
            resistors: &[],
            links: &LINKS,
        };
        let bad = bad_network.prepare(&thermal, || false).unwrap();
        assert!(matches!(
            source.step_with_isochoric_gas_transport(
                &thermal,
                &gas,
                1.0,
                0.0,
                &contacts([0.1, 0.2]),
                &bad,
                solve(),
                &[0.0; 3],
                cx
            ),
            Err(SolidGasError::Invalid(_))
        ));
        let cells = material_cells([&a, &b], curves, &context, &key);
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &LINKS,
        };
        let prepared = network.prepare(&thermal, || false).unwrap();
        // A rejected source-domain probe cannot leak either tentative history.
        assert!(
            source
                .step_with_isochoric_gas_transport(
                    &thermal,
                    &gas,
                    1.0,
                    0.0,
                    &contacts([0.1, 0.2]),
                    &prepared,
                    solve(),
                    &[1e6, 0.0, 0.0],
                    cx
                )
                .is_err()
        );
        let clone = ca.curve().clone();
        let wrong_law = material_cells([&a, &b], [&clone, cb.curve()], &context, &key);
        let wrong = MaterialNetwork {
            cells: &wrong_law,
            resistors: &[],
            links: &LINKS,
        };
        assert!(wrong.prepare(&thermal, || false).is_err());
        assert_eq!(thermal.cells(), &old);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, saved);
    });
}
