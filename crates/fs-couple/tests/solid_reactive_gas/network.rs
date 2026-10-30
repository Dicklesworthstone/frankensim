//! G1/G3/G4: independent finite gas inventories meet one spatial solid owner.
use super::*;
use fs_couple::solid_enthalpy::gas::network::*;

const CAPACITY: f64 = 2.5 * 100000.0 * VOLUME / 300.0;

fn controls() -> SolidGasNetworkSolve {
    SolidGasNetworkSolve {
        maximum_gases: 2,
        coupled: solve(),
    }
}

// Eliminate both constant-cv gas temperatures, then solve the remaining
// two-cell BE heat equations analytically. No native heat/EOS/solver call.
// Order is wall-left, wall-right, hot-gas, cold-gas.
fn reference(previous: [f64; 4], wall_k: f64, gas_k: [f64; 2], dt: f64) -> [f64; 4] {
    let g = [
        0.1 / (0.001 / gas_k[0] + 0.02 / wall_k),
        0.2 / (0.001 / gas_k[1] + 0.02 / wall_k),
    ];
    let a = g.map(|g| dt * g * CAPACITY / (CAPACITY + dt * g));
    let link = 0.5 * dt;
    let diagonal = a.map(|a| 10.0 + link + a);
    let rhs = [
        10.0 * previous[0] + a[0] * previous[2],
        10.0 * previous[1] + a[1] * previous[3],
    ];
    let determinant = diagonal[0] * diagonal[1] - link * link;
    let left = (rhs[0] * diagonal[1] + link * rhs[1]) / determinant;
    let right = (rhs[1] * diagonal[0] + link * rhs[0]) / determinant;
    [
        left,
        right,
        (CAPACITY * previous[2] + dt * g[0] * left) / (CAPACITY + dt * g[0]),
        (CAPACITY * previous[3] + dt * g[1] * right) / (CAPACITY + dt * g[1]),
    ]
}

fn temperatures<S>(candidate: &SolidGasNetworkCandidate<'_, '_, S>) -> [f64; 4] {
    [
        candidate.thermal().temperature_k(CellId(30)).unwrap(),
        candidate.thermal().temperature_k(CellId(10)).unwrap(),
        candidate.gas(20).unwrap().state().state().temperature,
        candidate.gas(10).unwrap().state().state().temperature,
    ]
}

#[test]
fn shared_wall_matches_four_reservoir_balance_and_native_solid_continuation() {
    with_cx(|cx, _| {
        let law = reaction(1.0, 0.0, 1.0);
        let mut changes = Vec::new();
        for wall_k in [0.1, 10.0] {
            let (card, caloric) = wall(wall_k, false);
            let mut thermal = initial(caloric.curve());
            let context = point(300.0);
            let model = solid(&thermal, &card, caloric.curve(), &context, cx);
            let mut hot = gas(&law);
            let mut cold = gas(&law);
            hot.add_heat(CAPACITY * 100.0, solve().gas, || false)
                .unwrap();
            cold.add_heat(-CAPACITY * 10.0, solve().gas, || false)
                .unwrap();
            let initial_energy = thermal.energy_j()
                + hot.mass_kg() * hot.state().specific_internal_energy_j_kg()
                + cold.mass_kg() * cold.state().specific_internal_energy_j_kg();
            let patches = contacts();
            let mut expected = [300.0, 300.0, 400.0, 290.0];
            for interval in 1..=3 {
                let inputs = [
                    SolidGasCompartment {
                        id: 20,
                        gas: &hot,
                        external_heat_j: 0.0,
                        contacts: &patches[..1],
                        reaction: None,
                    },
                    SolidGasCompartment {
                        id: 10,
                        gas: &cold,
                        external_heat_j: 0.0,
                        contacts: &patches[1..],
                        reaction: None,
                    },
                ];
                let next = model
                    .step_with_isochoric_gases(
                        &thermal,
                        &inputs,
                        interval as f64,
                        None,
                        controls(),
                        &[0.0; 4],
                        cx,
                    )
                    .unwrap();
                expected = reference(expected, wall_k, [0.03; 2], 1.0);
                for (actual, expected) in temperatures(&next).into_iter().zip(expected) {
                    near(actual, expected, 2e-7);
                }
                assert_eq!(
                    next.gases().iter().map(|g| g.id).collect::<Vec<_>>(),
                    vec![10, 20]
                );
                assert!(next.report().gases[0].solid_to_gas_heat_j > 0.0);
                assert!(next.report().gases[1].solid_to_gas_heat_j < 0.0);
                assert!(next.report().gases.iter().all(|g| g.endpoint.is_none()));
                assert!(
                    next.report().combined_energy_residual_j.abs()
                        <= solve().total_energy_tolerance_j
                );
                near(next.thermal().time_s(), interval as f64, 0.0);
                let loads = model
                    .geometry()
                    .assemble(&[], VolumeLoad3::None, &[], cx)
                    .unwrap();
                let free = next.state().solve(&loads, &FREE, config(), cx).unwrap();
                let u = free.solution.displacement.displacement_m();
                near(u[3][2], 1e-5 * (expected[0] - 300.0), 1e-10);
                near(u[4][2], -1e-5 * (expected[1] - 300.0), 1e-10);
                if interval == 1 {
                    let reversed = [inputs[1], inputs[0]];
                    let repeat = model
                        .step_with_isochoric_gases(
                            &thermal,
                            &reversed,
                            1.0,
                            None,
                            controls(),
                            &[0.0; 4],
                            cx,
                        )
                        .unwrap();
                    assert_eq!(temperatures(&next), temperatures(&repeat));
                    assert_eq!(next.report().trials, repeat.report().trials);
                }
                let (accepted, mut endpoints, _, _) = next.into_parts();
                hot = endpoints.pop().unwrap().gas;
                cold = endpoints.pop().unwrap().gas;
                thermal = accepted;
                near(
                    thermal.energy_j()
                        + hot.mass_kg() * hot.state().specific_internal_energy_j_kg()
                        + cold.mass_kg() * cold.state().specific_internal_energy_j_kg()
                        - initial_energy,
                    0.0,
                    8e-7,
                );
            }
            changes.push(expected[0]);
        }
        assert!(
            changes[1] - changes[0] > 10.0,
            "only wall material data changed"
        );
    });
}

#[test]
fn reacting_chamber_and_inert_coolant_share_one_wall_endpoint_without_mixing_species() {
    with_cx(|cx, _| {
        let (rate, release, product_k) = (1.0, 20000.0, 4.0);
        let law = super::coupled::implicit(reaction(rate, release, product_k));
        let reacting = gas(&law);
        let inert_law = reaction(1.0, 0.0, 1.0);
        let mut inert = gas(&inert_law);
        inert
            .add_heat(-10.0 * CAPACITY, solve().gas, || false)
            .unwrap();
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [30, 10].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.3,
            distance_b_m: 0.7,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        };
        let electrical =
            fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
        let network = fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        let patches = contacts();
        let inputs = [
            SolidGasCompartment {
                id: 20,
                gas: &reacting,
                external_heat_j: 0.0,
                contacts: &patches[..1],
                reaction: Some((&law, chemical())),
            },
            SolidGasCompartment {
                id: 10,
                gas: &inert,
                external_heat_j: 0.0,
                contacts: &patches[1..],
                reaction: None,
            },
        ];
        let next = model
            .step_with_isochoric_gases_heat_network(
                &thermal,
                &inputs,
                1.0,
                &prepared,
                controls(),
                &[0.0; 4],
                cx,
            )
            .unwrap();
        let [left, right, hot, cold] = temperatures(&next);
        let xi = super::coupled::extent(hot, 0.0, rate, release, 1.0);
        let amount =
            100000.0 * VOLUME / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 300.0);
        let kg = 0.03 * (1.0 + (product_k - 1.0) * xi / amount);
        let gh = 0.1 / (0.001 / kg + 0.02);
        let gc = 0.2 / (0.001 / 0.03 + 0.02);
        // Four independent BE balances, including formation energy once.
        near(
            10.0 * (left - 300.0) + 0.5 * (left - right) + gh * (left - hot),
            0.0,
            2e-7,
        );
        near(
            10.0 * (right - 300.0) + 0.5 * (right - left) + gc * (right - cold),
            0.0,
            2e-7,
        );
        near(
            CAPACITY * (hot - 300.0) - release * xi - gh * (left - hot),
            0.0,
            3e-7,
        );
        near(CAPACITY * (cold - 290.0) - gc * (right - cold), 0.0, 2e-7);
        let chemical = next.report().gases[1].endpoint.as_ref().unwrap();
        near(chemical.reaction.extent_mol, xi, 1e-9);
        assert_eq!(chemical.heat, next.report().gases[1].gas);
        assert!(next.report().gases[0].endpoint.is_none());
        assert!(next.report().transport.is_some());
        assert_eq!(next.gas(10).unwrap().model(), inert.model());
        assert_eq!(next.gas(10).unwrap().mass_kg(), inert.mass_kg());
        assert_eq!(next.gas(20).unwrap().mass_kg(), reacting.mass_kg());
        near(thermal.time_s(), 0.0, 0.0);
    });
}

#[test]
fn network_admission_late_refusal_and_cancellation_preserve_every_inventory() {
    with_cx(|cx, gate| {
        let law = reaction(1.0, 0.0, 1.0);
        let hot = gas(&law);
        let cold = gas(&law);
        let saved = [hot.clone(), cold.clone()];
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let patches = contacts();
        let mut inputs = [
            SolidGasCompartment {
                id: 20,
                gas: &hot,
                external_heat_j: 500.0,
                contacts: &patches[..1],
                reaction: None,
            },
            SolidGasCompartment {
                id: 10,
                gas: &cold,
                external_heat_j: 0.0,
                contacts: &patches[1..],
                reaction: None,
            },
        ];
        let mut budget = controls();
        budget.maximum_gases = 1;
        assert!(
            model
                .step_with_isochoric_gases(&thermal, &inputs, 1.0, None, budget, &[0.0; 4], cx)
                .is_err()
        );
        inputs[1].id = 20;
        assert!(
            model
                .step_with_isochoric_gases(&thermal, &inputs, 1.0, None, controls(), &[0.0; 4], cx)
                .is_err()
        );
        inputs[1].id = 10;
        let mut overlap = patches.clone();
        overlap[0].id += 1000;
        overlap[0].geometry.area_m2 = 2.0;
        overlap[1].id += 1000;
        overlap[1].cell = overlap[0].cell;
        overlap[1].geometry.area_m2 = 2.0;
        inputs[0].contacts = &overlap[..1];
        inputs[1].contacts = &overlap[1..];
        assert!(
            model
                .step_with_isochoric_gases(&thermal, &inputs, 1.0, None, controls(), &[0.0; 4], cx)
                .is_err()
        );
        inputs[0].contacts = &patches[..1];
        inputs[1].contacts = &patches[1..];
        let accepted = model
            .step_with_isochoric_gases(&thermal, &inputs, 1.0, None, controls(), &[0.0; 4], cx)
            .unwrap();
        let mut bounded = controls();
        bounded.coupled.maximum_trials = 1;
        assert!(matches!(
            model.step_with_isochoric_gases(&thermal, &inputs, 1.0, None, bounded, &[0.0; 4], cx),
            Err(SolidGasError::TrialBudget)
        ));
        let replay = model
            .step_with_isochoric_gases(&thermal, &inputs, 1.0, None, controls(), &[0.0; 4], cx)
            .unwrap();
        assert_eq!(temperatures(&accepted), temperatures(&replay));
        gate.request();
        assert!(matches!(
            model.step_with_isochoric_gases(
                &thermal,
                &inputs,
                1.0,
                None,
                controls(),
                &[0.0; 4],
                cx
            ),
            Err(SolidGasError::Cancelled)
        ));
        assert_eq!(hot, saved[0]);
        assert_eq!(cold, saved[1]);
        near(thermal.time_s(), 0.0, 0.0);
        assert!(accepted.thermal().time_s() > thermal.time_s());
    });
}
