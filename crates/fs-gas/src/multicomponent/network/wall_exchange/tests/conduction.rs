//! G1 series-resistance heat balance and G3 shared gas/solid material changes.
use super::*;
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::{CircularLaminarWall, MaterialWallConduction};
use crate::multicomponent::wall_exchange::{MaterialWallCaloric, WallCaloricConstraint};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, Dims, QuantitySpec, Temperature};
use std::sync::Arc;

// Manufactured source k(T)=a+bT and Cv(T)=c+dT for one explicit kilogram.
// These arbitrary units/state ranges are analytic fixtures, not metal data.
fn source(
    conductivity: [f64; 2],
    capacity: [f64; 2],
    interval: [f64; 2],
) -> (Arc<PreparedThermalTransport>, MaterialWallCaloric<'static>) {
    let keys = [
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("constant_volume_specific_heat", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let [lo, hi] = interval;
    let mut claims = ClaimSet::new();
    for (key, [a, b]) in keys.iter().zip([[1000.0, 0.0], capacity, conductivity]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(lo, a + b * lo), (hi, a + b * hi)],
                    dims: key.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured wall conduction fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "fixed solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                t,
            )
            .unwrap()
    };
    let (lower, upper) = (point(lo), point(hi));
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let transport = PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [selection; 3],
        },
    )
    .unwrap();
    let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 0,
        card: &card,
        capacity: &keys[1],
        amount: HeatStorageAmount::MassKg(1.0),
        temperature_axis: "T",
        lower: &lower,
        upper: &upper,
        selection,
    }])
    .unwrap();
    (
        Arc::new(transport),
        MaterialWallCaloric::new(
            Arc::new(storage),
            WallCaloricConstraint::IsochoricInternalEnergy,
        ),
    )
}

fn layer(k: [f64; 2], interval: [f64; 2]) -> MaterialWallConduction {
    MaterialWallConduction::new(
        source(k, [1.0, 0.0], interval).0,
        Length::new(0.1),
        Area::new(0.25),
    )
    .unwrap()
}
fn bind<'law>(
    wall: NetworkWallOperator<'law>,
    n: &MixtureNetwork1d,
    layer: &MaterialWallConduction,
) -> NetworkWallOperator<'law> {
    wall.with_material_conduction(
        n.branches()
            .iter()
            .map(|b| vec![Some(layer.clone()); b.cell_count()])
            .collect(),
        vec![Some(layer.clone()); n.plenums().len()],
        || false,
    )
    .unwrap()
}
fn conductance(contact: WallContact) -> f64 {
    match contact.thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Reservoir {
            conductance_w_per_k,
            ..
        } => conductance_w_per_k,
        _ => panic!("temperature-coordinate contact expected"),
    }
}

#[test]
fn conductivity_only_swap_changes_gas_pressure_through_series_heat_balance() {
    let g = gas(1, false, 0.0);
    let mut heat = Vec::new();
    for k in [0.05, 5.0] {
        let mut n = network(&g, &[0.0], true);
        let material = layer([k, 0.0], [1.0, 10.0]);
        let mut wall = bind(uniform_operator(&n, lumped()), &n, &material);
        let before = n.totals();
        let original = fields(&n);
        let report = wall
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        // Independent planar layer + film: R=1/(h*A)+L/(k*A).
        let total_g = 1.0 / (1.0 / 0.3 + 0.1 / (k * 0.25));
        let q = 0.2 * total_g * 3.0 / (1.0 + 0.2 * total_g * (2.0 + 1.0));
        for i in 0..2 {
            let state = n.branches()[0].cell(i).unwrap();
            close(
                g.primitive(state, || false).unwrap().pressure,
                2.0 + 2.0 * q,
                2e-10,
            );
            assert_eq!(state.partial_densities, original.0[0][i].partial_densities);
            close(
                wall_temperature(wall.branch_contacts(0).unwrap()[i]),
                5.0 - q,
                2e-10,
            );
            close(
                conductance(wall.branch_contacts(0).unwrap()[i]),
                total_g,
                1e-15,
            );
        }
        close(
            g.primitive(n.plenums()[0].state(), || false)
                .unwrap()
                .temperature,
            2.0 + 2.0 * q,
            2e-10,
        );
        close(wall_temperature(wall.plenum_contacts()[0]), 5.0 - q, 2e-10);
        close(report.heat_to_gas_j, 3.0 * q, 2e-10);
        close(report.wall_energy_change_j, -3.0 * q, 2e-10);
        close(n.totals().energy - before.energy, 3.0 * q, 2e-10);
        close(report.energy_residual_j, 0.0, 1e-10);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        assert_eq!(n.totals().species_masses, before.species_masses);
        assert_eq!(
            wall.branch_conduction(0, 0).unwrap().source(),
            material.source()
        );
        assert_eq!(wall.plenum_conduction(0).unwrap().area(), Area::new(0.25));
        assert!(wall.plenum_conduction(usize::MAX).is_none());
        heat.push(q);
    }
    assert!(heat[1] > 2.0 * heat[0]);
}

#[test]
fn common_material_capacity_and_conductivity_refresh_without_compounding_resistance() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let (transport, caloric) = source([0.2, 0.1], [1.0, 0.2], [1.0, 10.0]);
    let material =
        MaterialWallConduction::new(transport, Length::new(0.1), Area::new(0.25)).unwrap();
    let contact = caloric.contact(5.0, 0.3, WallDrag::None).unwrap();
    let wall = uniform_operator(&n, contact)
        .with_material_heat_storage(
            vec![vec![Some(caloric.clone()); 2]],
            vec![Some(caloric)],
            || false,
        )
        .unwrap();
    let mut wall = bind(wall, &n, &material);
    let (mut tw, mut tg) = (5.0_f64, 2.0_f64);
    for _ in 0..4 {
        let k = 0.2 + 0.1 * tw;
        let total_g = 1.0 / (1.0 / 0.3 + 0.1 / (k * 0.25));
        let alpha = 0.2 * total_g / (1.0 + 0.2 * total_g / 0.5);
        // Exact original Cwall(T)=1+0.2*T integral in the existing BE exchange.
        let constant = tw + 0.1 * tw * tw + alpha * tg;
        let linear = 1.0 + alpha;
        let next_tw = 2.0 * constant / (linear + (linear * linear + 0.4 * constant).sqrt());
        let q = alpha * (next_tw - tg);
        let report = wall
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        tg += q / 0.5;
        tw = next_tw;
        for contact in wall
            .branch_contacts(0)
            .unwrap()
            .iter()
            .chain(wall.plenum_contacts())
        {
            close(wall_temperature(*contact), tw, 3e-10);
            close(conductance(*contact), total_g, 3e-12);
        }
        close(
            g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .temperature,
            tg,
            3e-10,
        );
        close(report.heat_to_gas_j, 3.0 * q, 3e-10);
        close(report.wall_energy_change_j, -3.0 * q, 3e-10);
        close(report.energy_residual_j, 0.0, 2e-10);
    }
}

#[test]
fn gas_transport_and_solid_conduction_reach_the_same_native_exchange() {
    let g = gas(1, false, 0.0);
    let geometry = CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    };
    for gas_k in [0.02, 0.2] {
        let mut n = network(&g, &[0.2], false);
        let wall = uniform_operator(&n, lumped())
            .with_circular_laminar_transport(
                law(&g, 0.004, 0.0, gas_k, 0.0),
                vec![Some(geometry)],
                || false,
            )
            .unwrap();
        let mut wall = bind(wall, &n, &layer([0.1, 0.0], [1.0, 10.0]));
        let report = wall
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        let diameter = geometry.diameter.value();
        let film = 3.66 * gas_k * (4.0 * 0.2 / diameter) / diameter;
        let total_g = 1.0 / (1.0 / film + 0.1 / (0.1 * 0.25));
        let speed = 0.2 * (-32.0 * 0.004 / (diameter * diameter) * 0.2).exp();
        let heated = 2.0 + (0.2 * 0.2 - speed * speed) / 5.0;
        let q = 0.2 * total_g * (5.0 - heated) / (1.0 + 0.2 * total_g * 3.0);
        close(
            conductance(wall.branch_contacts(0).unwrap()[0]),
            total_g,
            1e-14,
        );
        close(
            g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .temperature,
            heated + q / 0.5,
            2e-10,
        );
        close(report.heat_to_gas_j, 2.0 * q, 2e-10);
        close(report.wall_energy_change_j, -2.0 * q, 2e-10);
        close(report.energy_residual_j, 0.0, 1e-10);
        assert!(report.drag_heating_j > 0.0);
    }
}

#[test]
fn unsupported_wall_endpoint_preserves_every_gas_wall_and_report() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let mut wall = bind(
        uniform_operator(&n, lumped()),
        &n,
        &layer([1.0, 0.0], [4.99, 6.0]),
    );
    wall.advance_for(&mut n, Time::new(0.001), controls(), || false)
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let plenums = wall.plenum_contacts().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::Conduction(_))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.plenum_contacts(), plenums);
    assert_eq!(wall.last_report(), report);
    wall.advance_for(&mut n, Time::new(0.001), controls(), || false)
        .unwrap();
}

#[test]
fn late_cancellation_preserves_accepted_state_and_replay_uses_original_resistances() {
    let g = gas(1, false, 0.0);
    let make = || {
        let mut n = network(&g, &[0.0], true);
        let mut wall = bind(
            uniform_operator(&n, lumped()),
            &n,
            &layer([0.1, 0.2], [1.0, 10.0]),
        );
        wall.advance_for(&mut n, Time::new(0.01), controls(), || false)
            .unwrap();
        (n, wall)
    };
    let (mut expected, mut reference) = make();
    let mut calls = 0;
    let report = reference
        .advance_for(&mut expected, Time::new(0.2), controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    let (mut n, mut wall) = make();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let plenums = wall.plenum_contacts().to_vec();
    let previous = wall.last_report();
    let mut count = 0;
    assert_eq!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || {
            count += 1;
            count == calls
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.plenum_contacts(), plenums);
    assert_eq!(wall.last_report(), previous);
    assert_eq!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap(),
        report
    );
    assert_eq!(fields(&n), fields(&expected));
    assert_eq!(wall.branch_contacts(0), reference.branch_contacts(0));
    assert_eq!(wall.plenum_contacts(), reference.plenum_contacts());
}

#[test]
fn disabled_film_stays_disabled_and_unsupported_wall_models_refuse() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], false);
    let zero = WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Reservoir {
            temperature_k: 5.0,
            conductance_w_per_k: 0.0,
        },
    };
    let material = layer([1.0, 0.0], [1.0, 10.0]);
    let mut wall = bind(uniform_operator(&n, zero), &n, &material);
    let before = fields(&n);
    let report = wall
        .advance_for(&mut n, Time::new(1.0), controls(), || false)
        .unwrap();
    assert_eq!(fields(&n), before);
    assert_eq!(report.heat_to_gas_j, 0.0);
    assert_eq!(report.external_heat_to_gas_j, 0.0);
    assert_eq!(wall.branch_contacts(0).unwrap(), [zero, zero]);
    for thermal in [
        WallThermal::Adiabatic,
        WallThermal::Caloric {
            mass_kg: 1.0,
            specific_energy_j_kg: 5.0,
            conductance_w_per_k: 0.3,
        },
    ] {
        let contact = WallContact {
            drag: WallDrag::None,
            thermal,
        };
        assert!(
            uniform_operator(&n, contact)
                .with_material_conduction(vec![vec![Some(material.clone()); 2]], vec![], || false,)
                .is_err()
        );
    }
    assert!(
        MaterialWallConduction::new(
            Arc::new(material.source().unwrap().clone()),
            Length::new(0.0),
            Area::new(0.25)
        )
        .is_err()
    );
}

#[test]
fn wall_material_changes_reacting_network_pressure_and_composition_with_closed_energy() {
    use crate::chemistry::IsochoricKineticsSolve;
    use crate::multicomponent::network::reactive::{ReactingNetwork1d, ReactiveStepControl};
    use fs_thermochem::kinetics::{
        ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
    };
    use fs_thermochem::{MolarEnergyQuantityV1, ReactionId, SpeciesId};

    let g = gas(2, false, 0.0);
    // Reversible isomers with a temperature-dependent rate. Equal calorics
    // isolate wall-driven pressure and chemistry feedback from heat release.
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("wall-conduction-isomerization").unwrap(),
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
        "manufactured-wall-coupled-isomers",
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
    for k in [0.05, 5.0] {
        let branch = MixtureEuler1d::new(
            g.clone(),
            Length::new(0.5),
            0.4,
            vec![g.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap(); 2],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
        let mut wall = bind(
            uniform_operator(&n, lumped()),
            &n,
            &layer([k, 0.0], [1.0, 10.0]),
        );
        let initial_energy = n.totals().absolute_energy(&g).unwrap() + 2.0 * 5.0;
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
            assert_eq!(report.elapsed, Time::new(0.02));
            assert_eq!(wall.last_report().unwrap().elapsed, report.elapsed);
            let wall_energy: f64 = wall
                .branch_contacts(0)
                .unwrap()
                .iter()
                .map(|&contact| wall_temperature(contact))
                .sum();
            close(
                flow.totals().absolute_energy(&g).unwrap() + wall_energy,
                initial_energy,
                2e-9,
            );
            close(flow.totals().mass(), initial_mass, 1e-12);
            close(report.absolute_energy_residual_j, 0.0, 2e-9);
        }
        let state = flow.branches()[0].cell(0).unwrap();
        assert!(
            state.partial_densities[1] > 0.0,
            "the actual reactor must create product"
        );
        outcomes.push((
            g.primitive(state, || false).unwrap().pressure,
            state.partial_densities[1],
        ));
    }
    assert!(
        outcomes[1].0 > outcomes[0].0,
        "higher wall conductivity must heat the gas faster"
    );
    assert!(
        outcomes[1].1 > outcomes[0].1,
        "actual wall heating must reach the chemical rate law"
    );
}
