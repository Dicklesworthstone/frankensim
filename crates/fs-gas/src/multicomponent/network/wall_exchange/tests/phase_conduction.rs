//! G1 original latent inventory and G3 phase-source changes in real gas exchange.
use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, PhaseWallConductionSource,
};
use crate::multicomponent::wall_exchange::MaterialWallCaloric;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Area, Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell, MaterialId};
use std::sync::Arc;

// Manufactured chart: cp=1, Tm=5, latent heat=10, rho=20. On its
// isothermal plateau phi=(h-4)/10, and the independently authored source is
// k=scale*(0.1+1.9*phi). These numbers are not measured material data.
fn source(scale: f64, maximum_fraction: f64, anchor: f64) -> Arc<PhaseWallConductionSource> {
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut claims = ClaimSet::new();
    let claim = |end: f64, multiplier: f64| PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Curve {
            abscissa: "phi".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![
                (0.0, multiplier * 0.1),
                (end, multiplier * (0.1 + 1.9 * end)),
            ],
            dims: key.dims(),
        },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured phase wall".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    let selected = claims.insert_claim(claim(maximum_fraction, scale)).unwrap();
    if maximum_fraction < 1.0 {
        // A tempting replacement must not be selected when the original
        // narrow source runs out of support during a melting step.
        claims.insert_claim(claim(1.0, 10.0 * scale)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid-liquid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = Arc::new(
        EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            [
                (0.0, 1.0, 0.0),
                (4.0, 5.0, 0.0),
                (14.0, 5.0, 1.0),
                (19.0, 10.0, 1.0),
            ]
            .map(|(h, t, phi)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: phi,
                bulk_density_kg_m3: 20.0,
            })
            .to_vec(),
        )
        .unwrap(),
    );
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 5.0)
        .unwrap()
        .with_quantity("phi", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    Arc::new(
        PhaseWallConductionSource::new(
            curve,
            ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("phi"),
            },
            &key,
            MaterialPropertySelection::PinnedByProperty(vec![(key.name().to_owned(), selected)]),
            anchor,
        )
        .unwrap(),
    )
}

fn cell(h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 1.0,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.0, 0.0, 0.2],
            volume_m3: 0.05,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 1.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn wall(
    n: &MixtureNetwork1d,
    source: &Arc<PhaseWallConductionSource>,
    h: f64,
) -> NetworkWallOperator<'static> {
    let caloric =
        MaterialWallCaloric::from_phase_curve(source.curve().clone(), cell(h), [0.0, 19.0])
            .unwrap();
    assert!(Arc::ptr_eq(caloric.phase_curve().unwrap(), source.curve()));
    let contact = caloric.initial_contact(0.3, WallDrag::None).unwrap();
    let layer = MaterialWallConduction::from_phase_source(
        source.clone(),
        Length::new(0.1),
        Area::new(0.25),
    )
    .unwrap();
    uniform_operator(n, contact)
        .with_material_heat_storage(
            n.branches()
                .iter()
                .map(|b| vec![Some(caloric.clone()); b.cell_count()])
                .collect(),
            vec![Some(caloric); n.plenums().len()],
            || false,
        )
        .unwrap()
        .with_material_conduction(
            n.branches()
                .iter()
                .map(|b| vec![Some(layer.clone()); b.cell_count()])
                .collect(),
            vec![Some(layer); n.plenums().len()],
            || false,
        )
        .unwrap()
}
fn energy(contact: WallContact) -> f64 {
    let WallThermal::Caloric {
        specific_energy_j_kg,
        ..
    } = contact.thermal
    else {
        panic!("original enthalpy coordinate required")
    };
    specific_energy_j_kg
}
fn combined(k: f64) -> f64 {
    1.0 / (1.0 / 0.3 + 0.1 / (0.25 * k))
}
fn hot_network(g: &MixtureGas, two_species: bool) -> MixtureNetwork1d {
    let fractions = if two_species {
        vec![1.0, 0.0]
    } else {
        vec![1.0]
    };
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&fractions, 0.0, 8.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap()
}

#[test]
fn equal_temperature_latent_inventories_and_source_swaps_change_network_heat() {
    let g = gas(1, false, 0.0);
    let mut outcomes = Vec::new();
    for scale in [1.0, 3.0] {
        let source = source(scale, 1.0, 5.0);
        for h in [5.0, 13.0] {
            let mut n = network(&g, &[0.0], true);
            let mut wall = wall(&n, &source, h);
            let before = n.totals().absolute_energy(&g).unwrap();
            let report = wall
                .advance_for(&mut n, Time::new(0.2), controls(), || false)
                .unwrap();
            let phi = (h - 4.0) / 10.0;
            let conductance = combined(scale * (0.1 + 1.9 * phi));
            // On the actual latent plateau Tw=5 exactly while h loses Q.
            let q = 0.2 * conductance * 3.0 / (1.0 + 0.2 * conductance / 0.5);
            for contact in wall
                .branch_contacts(0)
                .unwrap()
                .iter()
                .chain(wall.plenum_contacts())
            {
                close(energy(*contact), h - q, 2e-11);
                close(
                    source
                        .curve()
                        .temperature_at_specific_enthalpy(energy(*contact))
                        .unwrap(),
                    5.0,
                    0.0,
                );
            }
            close(
                n.totals().absolute_energy(&g).unwrap() - before,
                3.0 * q,
                3e-11,
            );
            close(report.wall_energy_change_j, -3.0 * q, 3e-11);
            close(report.energy_residual_j, 0.0, 3e-11);
            let pressure = g
                .primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .pressure;
            close(pressure, 2.0 + q / 0.5, 3e-11);
            let binding = wall.branch_conduction(0, 0).unwrap();
            assert!(binding.source().is_none());
            let exact = binding
                .phase_source()
                .unwrap()
                .at_specific_enthalpy(h)
                .unwrap();
            assert_eq!(
                exact.properties()[0].answer().receipt.selected,
                source.reference().properties()[0].answer().receipt.selected
            );
            outcomes.push(pressure);
        }
    }
    assert!(outcomes[1] > outcomes[0]);
    assert!(outcomes[2] > outcomes[0]);
    assert!(outcomes[3] > outcomes[1]);
}

#[test]
fn melting_changes_conduction_at_the_same_temperature_without_extra_latent_heat() {
    let g = gas(1, false, 0.0);
    let mut n = hot_network(&g, false);
    let source = source(1.0, 1.0, 3.99);
    let mut wall = wall(&n, &source, 3.99);
    let initial = n.totals().absolute_energy(&g).unwrap() + 2.0 * 3.99;
    let (mut h, mut tg) = (3.99, 8.0);
    let mut coefficients = Vec::new();
    for _ in 0..3 {
        let phi = if h < 4.0 { 0.0 } else { (h - 4.0) / 10.0 };
        let coefficient = combined(0.1 + 1.9 * phi);
        let q = 0.2 * coefficient * (5.0 - tg) / (1.0 + 0.2 * coefficient / 0.5);
        h -= q;
        tg += q / 0.5;
        assert!((4.0..14.0).contains(&h));
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        for contact in wall.branch_contacts(0).unwrap() {
            close(energy(*contact), h, 3e-11);
            let WallThermal::Caloric {
                conductance_w_per_k,
                ..
            } = contact.thermal
            else {
                unreachable!()
            };
            close(conductance_w_per_k, coefficient, 1e-14);
        }
        close(
            g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .temperature,
            tg,
            3e-11,
        );
        close(
            n.totals().absolute_energy(&g).unwrap() + 2.0 * h,
            initial,
            4e-11,
        );
        coefficients.push(coefficient);
    }
    assert!(coefficients[1] > coefficients[0] && coefficients[2] > coefficients[1]);
}

#[test]
fn phase_feedback_refines_toward_independent_continuous_energy_solution() {
    let g = gas(1, false, 0.0);
    let source = source(1.0, 1.0, 5.0);
    // Closed two-body energy gives Tg=(9-h)/0.5. On the plateau,
    // dh/dt=G(h)*(6.5-h)/0.5 and k(h)=0.19*h-0.66. Integrating
    // dt/dh yields these logarithms independently of the discrete exchange.
    let conductivity = |h: f64| 0.1 + 1.9 * (h - 4.0) / 10.0;
    let elapsed = |h: f64| {
        let ratio = ((5.0 - 6.5) / (h - 6.5)).ln();
        0.5 / 0.3 * ratio
            + 0.5 * 0.1 / (0.25 * conductivity(6.5))
                * (ratio - (conductivity(5.0) / conductivity(h)).ln())
    };
    let (mut lower, mut upper) = (5.0, 6.5 - 1e-12);
    for _ in 0..90 {
        let middle = 0.5 * (lower + upper);
        if elapsed(middle) < 0.4 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let exact = 0.5 * (lower + upper);
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut n = hot_network(&g, false);
        let mut wall = wall(&n, &source, 5.0);
        for _ in 0..steps {
            wall.advance_for(&mut n, Time::new(0.4 / steps as f64), controls(), || false)
                .unwrap();
        }
        errors.push((energy(wall.branch_contacts(0).unwrap()[0]) - exact).abs());
    }
    assert!(
        errors[2] < 0.0016,
        "finest phase-feedback error={:e}",
        errors[2]
    );
    assert!(
        errors[1] < 0.56 * errors[0] && errors[2] < 0.56 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn phase_source_endpoint_and_chart_mismatch_refuse_atomically_without_source_switching() {
    let g = gas(1, false, 0.0);
    let source = source(1.0, 0.001, 3.99);
    let mut n = hot_network(&g, false);
    let mut wall = wall(&n, &source, 3.99);
    wall.advance_for(&mut n, Time::new(0.001), controls(), || false)
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::PhaseConduction(_))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);
    // The unsupported candidate does not consume latent inventory or a film.
    wall.advance_for(&mut n, Time::new(0.001), controls(), || false)
        .unwrap();

    let mut knots = source.curve().knots().to_vec();
    knots[2].specific_enthalpy_j_kg += 1.0;
    let foreign = Arc::new(
        EquilibriumEnthalpyPhaseCurve::try_new(source.material().content_hash(), knots).unwrap(),
    );
    let caloric = MaterialWallCaloric::from_phase_curve(foreign, cell(3.99), [0.0, 19.0]).unwrap();
    let layer =
        MaterialWallConduction::from_phase_source(source, Length::new(0.1), Area::new(0.25))
            .unwrap();
    let candidate = uniform_operator(&n, caloric.initial_contact(0.3, WallDrag::None).unwrap())
        .with_material_heat_storage(vec![vec![Some(caloric); 2]], vec![], || false)
        .unwrap();
    assert!(
        candidate
            .with_material_conduction(vec![vec![Some(layer); 2]], vec![], || false)
            .is_err()
    );
}

#[test]
fn late_cancellation_preserves_phase_inventory_and_retry_replays_exactly() {
    let g = gas(1, false, 0.0);
    let source = source(1.0, 1.0, 3.99);
    let mut reference = hot_network(&g, false);
    let mut reference_wall = wall(&reference, &source, 3.99);
    let mut count = 0;
    let expected = reference_wall
        .advance_for(&mut reference, Time::new(0.2), controls(), || {
            count += 1;
            false
        })
        .unwrap();
    let mut n = hot_network(&g, false);
    let mut wall = wall(&n, &source, 3.99);
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let mut polls = 0;
    assert_eq!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || {
            polls += 1;
            polls == count
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), None);
    let actual = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(fields(&n), fields(&reference));
    assert_eq!(wall.branch_contacts(0), reference_wall.branch_contacts(0));
}

#[test]
fn phase_conductivity_reaches_native_reaction_rates_and_closed_energy() {
    use crate::chemistry::IsochoricKineticsSolve;
    use crate::multicomponent::network::reactive::{ReactingNetwork1d, ReactiveStepControl};
    use fs_thermochem::kinetics::{
        ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
    };
    use fs_thermochem::{MolarEnergyQuantityV1, ReactionId, SpeciesId};

    let g = gas(2, false, 0.0);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("phase-wall-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(1.0),
            2.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        g.source().clone(),
        "manufactured-phase-wall-isomers",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let chemistry = IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10,
        relative_tolerance: 2e-7,
        absolute_amount_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 0.02,
        maximum_steps: 30_000,
        ..IsochoricKineticsSolve::default()
    };
    let mut outcomes = Vec::new();
    for scale in [0.1, 10.0] {
        let source = source(scale, 1.0, 5.0);
        let n = hot_network(&g, true);
        let mut wall = wall(&n, &source, 5.0);
        let initial = n.totals().absolute_energy(&g).unwrap() + 2.0 * 5.0;
        let initial_mass = n.totals().mass();
        let mut flow = ReactingNetwork1d::new(n, mechanism.clone(), chemistry, || false).unwrap();
        for _ in 0..3 {
            let report = flow
                .step_with_wall(
                    &mut wall,
                    Time::new(0.02),
                    ReactiveStepControl::default(),
                    controls(),
                    || false,
                )
                .unwrap();
            let wall_energy: f64 = wall
                .branch_contacts(0)
                .unwrap()
                .iter()
                .map(|&contact| energy(contact))
                .sum();
            close(
                flow.totals().absolute_energy(&g).unwrap() + wall_energy,
                initial,
                2e-9,
            );
            close(flow.totals().mass(), initial_mass, 1e-12);
            close(report.absolute_energy_residual_j, 0.0, 2e-9);
            assert!(
                wall_energy > 10.0,
                "gas heat must enter actual latent inventory"
            );
        }
        let state = flow.branches()[0].cell(0).unwrap();
        outcomes.push((
            g.primitive(state, || false).unwrap().pressure,
            state.partial_densities[1],
        ));
    }
    assert!(
        outcomes[1].0 < outcomes[0].0,
        "stronger phase-wall conduction cools gas faster"
    );
    assert!(
        outcomes[1].1 < outcomes[0].1,
        "actual cooling must reach temperature-dependent chemistry"
    );
}
