//! One complete coupled multiplier must survive both parameter-owning wrappers.

use super::*;

#[test]
fn robin_feedback_composes_with_radiation_and_reference_density_controls() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let source = ScalarField::nodal("heat density", 4, SOURCE.to_vec()).unwrap();
        let assigned = |density| {
            HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &fixture.mesh,
                &[ReferenceEnthalpyMaterial {
                    curve: &fixture.curve,
                    reference_density_kg_m3: density,
                }],
                &[0],
                budget(),
            )
            .unwrap()
        };
        let stepper = assigned(RHO);
        let patches = [patch(500.0, 520.0)];
        let accepted = stepper
            .advance_with_ambient_radiation(
                cx,
                fixture.problem(&source),
                None,
                &OLD,
                DT,
                config(),
                &patches,
                radiation_config(),
            )
            .unwrap()
            .conduction;
        // Equal face areas on this regular tetrahedron give each vertex the
        // same surface-mean weight. Keep this law anchor FIXED in every solve.
        let anchor = accepted.temperature.iter().sum::<f64>() / 4.0;
        let coupled = stepper
            .linearize_accepted_parameters_with_ambient_radiation(
                cx,
                fixture.problem(&source),
                None,
                &OLD,
                DT,
                config(),
                &patches,
                accepted,
                8,
            )
            .unwrap()
            .with_robin_feedback(cx, &["exterior"], &[vec![0.3]], &[0.02], 24)
            .unwrap();
        let original = coupled
            .transport()
            .robin_response(cx, &["exterior"], linear())
            .unwrap();
        assert_eq!(original.ports()[0].htc_w_m2_k, HTC);
        assert_eq!(original.ports()[0].reference_k, 300.0);
        let seed = coupled.transport().temperature_pullback(cx, &GOAL).unwrap();
        let gradient = coupled.pullback(cx, &seed, linear()).unwrap();
        let actual_dual = coupled
            .transport()
            .apply_jacobian_transpose(cx, &gradient.storage.transport.adjoint)
            .unwrap();
        assert!(
            actual_dual
                .iter()
                .zip(&seed)
                .all(|(a, b)| (a - b).abs() < 2e-12)
        );

        let objective = |density, ambient| {
            let perturbed = assigned(density);
            let mut mean = anchor;
            for _ in 0..64 {
                let boundary = ThermalBoundaryBuilder::new(&fixture.mesh)
                    .region(
                        "exterior",
                        |_| true,
                        ThermalBc::robin(
                            HTC + 0.02 * (mean - anchor),
                            300.0 + 0.3 * (mean - anchor),
                        )
                        .unwrap(),
                    )
                    .unwrap()
                    .finish()
                    .unwrap();
                let solved = perturbed
                    .advance_with_ambient_radiation(
                        cx,
                        ConductionProblem {
                            boundary: &boundary,
                            ..fixture.problem(&source)
                        },
                        None,
                        &OLD,
                        DT,
                        config(),
                        &[patch(ambient, 520.0)],
                        radiation_config(),
                    )
                    .unwrap()
                    .conduction;
                let next = solved.temperature.iter().sum::<f64>() / 4.0;
                if (next - mean).abs() < 1e-10 {
                    return solved
                        .temperature
                        .iter()
                        .zip(GOAL)
                        .map(|(t, w)| w * (t - 350.0))
                        .sum();
                }
                mean = next;
            }
            panic!("bounded physical Robin/radiation iteration did not converge");
        };
        agree(
            gradient.storage.material_reference_density[0],
            five_point(0.01, |delta| objective(RHO + delta, 500.0)),
            "Robin/radiation reference density",
        );
        agree(
            gradient.reservoir_temperatures[0],
            five_point(0.05, |delta| objective(RHO, 500.0 + delta)),
            "Robin/radiation reservoir temperature",
        );
        assert_eq!(
            coupled
                .into_storage()
                .pullback(cx, &seed, linear())
                .unwrap(),
            gradient.storage,
        );
    });
}
