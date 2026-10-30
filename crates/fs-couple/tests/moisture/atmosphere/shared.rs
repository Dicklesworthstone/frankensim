//! G1/G3/G4: shared finite vapor, independent balance and native material response.
use super::*;

fn participant<'a>(
    id: u64,
    slab: &'a MoistureSlab,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> MoistureParticipant<'a> {
    MoistureParticipant {
        id,
        slab,
        faces,
        transport: config(),
    }
}

#[test]
fn heterogeneous_specimens_adsorb_and_desorb_through_one_conserved_inventory() {
    with_cx(|cx| {
        let mut gas = atmosphere(0.5, 1.0);
        let mut a = slab(card([0.01; 2], 0.3), 1.0, vec![0.02], cx);
        let mut b = slab(card([0.03; 2], 0.3), 0.7, vec![0.2], cx);
        let faces_a = [Some(sorption(0.2)), Some(sorption(0.2))];
        let faces_b = [Some(sorption(0.3)), Some(sorption(0.3))];
        let saturation = saturation_pressure() * WATER_MOLAR_MASS / (R * TEMPERATURE);
        let carrier = carrier_mass(&gas);
        let initial_water =
            saturation * gas.relative_humidity().unwrap() + a.water_mass_kg() + b.water_mass_kg();
        let dt = 2.0;
        let alpha = [4.0 * dt * 0.01, 4.0 * dt * 0.03 / (0.7 * 0.7)];
        for interval in 0..3 {
            let old_rh = gas.relative_humidity().unwrap();
            let old = [a.moisture()[0], b.moisture()[0]];
            let mass = [a.dry_mass_kg(), b.dry_mass_kg()];
            let slope = [0.2, 0.3];
            // Eliminate both one-cell backward-Euler diffusion equations from
            // S*(r-r0) + sum(m_dry*(M-M0)) = 0, independently of the root solver.
            let weight = [
                mass[0] * alpha[0] / (1.0 + alpha[0]),
                mass[1] * alpha[1] / (1.0 + alpha[1]),
            ];
            let expected_rh = (saturation * old_rh + weight[0] * old[0] + weight[1] * old[1])
                / (saturation + weight[0] * slope[0] + weight[1] * slope[1]);
            let expected: [f64; 2] = core::array::from_fn(|i| {
                (old[i] + alpha[i] * slope[i] * expected_rh) / (1.0 + alpha[i])
            });
            let candidate = gas
                .advance_shared(
                    &[participant(7, &a, &faces_a), participant(3, &b, &faces_b)],
                    dt,
                    2,
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                candidate
                    .solids()
                    .iter()
                    .map(|body| body.id)
                    .collect::<Vec<_>>(),
                vec![3, 7]
            );
            near(
                candidate.atmosphere().relative_humidity().unwrap(),
                expected_rh,
                2e-11,
            );
            for (id, i) in [(7, 0), (3, 1)] {
                let body = candidate
                    .solids()
                    .iter()
                    .find(|body| body.id == id)
                    .unwrap();
                near(body.moisture.moisture()[0], expected[i], 2e-12);
                near(
                    body.transport.boundary_water_kg.iter().sum(),
                    mass[i] * (expected[i] - old[i]),
                    2e-13,
                );
                assert_eq!(body.moisture.time_s(), (interval + 1) as f64 * dt);
                for source in body.transport.surface_sorption.iter().flatten() {
                    near(
                        source.property(SORPTION_PROPERTY).unwrap().value_si(),
                        slope[i] * expected_rh,
                        1e-11,
                    );
                }
            }
            assert!(expected[0] > old[0]);
            assert!(expected[1] < old[1]);
            near(candidate.step().water_balance_residual_kg, 0.0, 1e-13);
            near(candidate.step().total_mass_residual_kg, 0.0, 1e-13);
            near(carrier_mass(candidate.atmosphere()), carrier, 2e-14);
            near(
                candidate.step().vapor_mass_kg[1]
                    + candidate
                        .solids()
                        .iter()
                        .map(|body| body.moisture.water_mass_kg())
                        .sum::<f64>(),
                initial_water,
                5e-13,
            );
            let (next_gas, solids, _) = candidate.into_parts();
            let mut solids = solids.into_iter();
            b = solids.next().unwrap().moisture;
            a = solids.next().unwrap().moisture;
            gas = next_gas;
        }
    });
}

#[test]
fn companion_material_changes_the_same_plate_and_gas_without_enumeration_bias() {
    with_cx(|cx| {
        let gas = atmosphere(0.7, 1.0);
        let observed = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 4], cx);
        let binding = plate(&observed, PlateMaterialModel::Isotropic, cx);
        let faces = [Some(sorption(0.25)), Some(sorption(0.25))];
        let mut outcomes = Vec::new();
        for diffusivity in [0.005, 0.08] {
            // Only the OTHER specimen's source diffusivity changes.
            let companion = slab(card([diffusivity; 2], 0.3), 1.0, vec![0.02; 3], cx);
            let candidate = gas
                .advance_shared(
                    &[
                        participant(10, &observed, &faces),
                        participant(20, &companion, &faces),
                    ],
                    2.0,
                    7,
                    controls(),
                    cx,
                )
                .unwrap();
            let reordered = gas
                .advance_shared(
                    &[
                        participant(20, &companion, &faces),
                        participant(10, &observed, &faces),
                    ],
                    2.0,
                    7,
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(candidate.atmosphere(), reordered.atmosphere());
            for (a, b) in candidate.solids().iter().zip(reordered.solids()) {
                assert_eq!(a.id, b.id);
                assert_eq!(a.moisture.moisture(), b.moisture.moisture());
                assert_eq!(a.transport.boundary_water_kg, b.transport.boundary_water_kg);
            }
            let solid = &candidate.solids()[0];
            let mechanical = binding.resolve(&solid.moisture, cx).unwrap();
            near(
                mechanical.wet_mass_kg(),
                solid.moisture.wet_mass_kg(),
                1e-13,
            );
            for (source, moisture) in mechanical
                .elastic_states()
                .iter()
                .zip(solid.moisture.moisture())
            {
                near(
                    source.property("young_modulus").unwrap().value_si(),
                    1e9 * (1.0 - moisture),
                    1e-6,
                );
            }
            outcomes.push((
                candidate.atmosphere().relative_humidity().unwrap(),
                mechanical.wet_mass_kg(),
                mechanical.elastic_states()[0]
                    .property("young_modulus")
                    .unwrap()
                    .value_si(),
                candidate.atmosphere().gas().state().state().sound_speed,
            ));
        }
        assert!(outcomes[1].0 < outcomes[0].0);
        assert!(outcomes[1].1 < outcomes[0].1);
        assert!(outcomes[1].2 > outcomes[0].2);
        assert_ne!(outcomes[1].3, outcomes[0].3);
        assert_eq!(observed.time_s(), 0.0);
        assert_eq!(observed.moisture(), &[0.1; 4]);
    });
}

#[test]
fn late_specimen_refusal_clock_budget_and_cancellation_publish_no_partial_exchange() {
    with_cx(|cx| {
        let gas = atmosphere(0.8, 1.0);
        let saved_gas = gas.clone();
        let a = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 2], cx);
        let b = slab(card([0.02; 2], 0.3), 0.5, vec![0.03; 2], cx);
        let faces = [Some(sorption(0.25)), Some(sorption(0.25))];
        let mut late = participant(2, &b, &faces);
        late.transport.maximum_content_increment = 1e-10;
        assert!(
            gas.advance_shared(&[participant(1, &a, &faces), late], 2.0, 4, controls(), cx)
                .is_err()
        );
        let valid = gas
            .advance_shared(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                2.0,
                4,
                controls(),
                cx,
            )
            .unwrap();
        assert!(valid.solids()[1].moisture.moisture()[0] - b.moisture()[0] > 1e-10);
        assert!(
            gas.advance_shared(
                &[participant(1, &a, &faces), participant(2, &b, &faces)],
                2.0,
                3,
                controls(),
                cx,
            )
            .is_err()
        );
        assert!(
            gas.advance_shared(
                &[participant(1, &a, &faces), participant(1, &b, &faces)],
                2.0,
                4,
                controls(),
                cx,
            )
            .is_err()
        );
        assert!(
            gas.advance_shared(
                &[
                    participant(1, &a, &faces),
                    participant(2, &valid.solids()[1].moisture, &faces)
                ],
                2.0,
                4,
                controls(),
                cx,
            )
            .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                gas.advance_shared(
                    &[participant(1, &a, &faces), participant(2, &b, &faces)],
                    2.0,
                    4,
                    controls(),
                    &cancelled,
                ),
                Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
            ));
        });
        assert_eq!(gas, saved_gas);
        assert_eq!(a.time_s(), 0.0);
        assert_eq!(b.time_s(), 0.0);
        assert_eq!(a.moisture(), &[0.1; 2]);
        assert_eq!(b.moisture(), &[0.03; 2]);
    });
}
