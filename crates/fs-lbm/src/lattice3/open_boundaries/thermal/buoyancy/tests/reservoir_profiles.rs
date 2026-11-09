use super::*;

#[test]
fn buoyant_profile_parity_domain_refusal_and_cancel_preserve_all_epochs_then_retry_exactly() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let heat = vec![0.0004; 72];
    actual
        .set_reservoir_temperature_profile(0, &vec![0.5; 12])
        .unwrap();
    actual
        .set_reservoir_temperature_profile(1, &vec![0.7; 12])
        .unwrap();
    // Uniform profiles retain exact scalar-boundary behavior through momentum,
    // source half-kicks, conducting finite walls, buoyancy and publication.
    for _ in 0..3 {
        assert_eq!(advance(&mut actual, &heat), advance(&mut reference, &heat));
        assert_eq!(accepted(&actual), accepted(&reference));
    }
    let profiles = [
        (0..12).map(|k| 0.45 + 0.01 * k as f64).collect::<Vec<_>>(),
        (0..12).map(|k| 0.8 + 0.01 * k as f64).collect::<Vec<_>>(),
    ];
    let before_inputs = accepted(&actual);
    for face in 0..2 {
        actual
            .set_reservoir_temperature_profile(face, &profiles[face])
            .unwrap();
        reference
            .set_reservoir_temperature_profile(face, &profiles[face])
            .unwrap();
    }
    assert_eq!(accepted(&actual), before_inputs);
    let mut refused = profiles[0].clone();
    // Within the scalar carrier's [-2,3] range, outside the actual affine
    // density-departure envelope. A late offending entry must leave all inputs.
    refused[11] = 2.0;
    assert!(matches!(
        actual.set_reservoir_temperature_profile(0, &refused),
        Err(BoussinesqOpenError3::BoundaryLaw { face: 0, .. })
    ));
    assert_eq!(
        actual.thermal().reservoir_temperature_profile(0).unwrap(),
        profiles[0]
    );
    assert_eq!(accepted(&actual), before_inputs);

    let energy_before = total_energy(&actual);
    let mut background: Vec<_> = (0..72)
        .map(|cell| {
            [
                BACKGROUND[0] + 1e-6 * (1 + cell % 3) as f64,
                BACKGROUND[1] - 1e-6,
                BACKGROUND[2],
            ]
        })
        .collect();
    let authored_background = background.clone();
    let mut expected_background = background.clone();
    let mut polls = 0;
    let expected = reference
        .advance_with_diffusion_poll(
            Some(&mut expected_background),
            None,
            &heat,
            None,
            None,
            [0.0; 3],
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let before = accepted(&actual);
    let mut calls = 0;
    let cancelled = actual.advance_with_diffusion_poll(
        Some(&mut background),
        None,
        &heat,
        None,
        None,
        [0.0; 3],
        &mut || {
            calls += 1;
            if calls == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_eq!(accepted(&actual), before);
    assert_eq!(background, authored_background);
    for face in 0..2 {
        assert_eq!(
            actual
                .thermal()
                .reservoir_temperature_profile(face)
                .unwrap(),
            profiles[face]
        );
    }
    let retry = actual
        .advance_with_diffusion_poll(
            Some(&mut background),
            None,
            &heat,
            None,
            None,
            [0.0; 3],
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(accepted(&actual), accepted(&reference));
    assert_eq!(background, expected_background);
    let solid_source: f64 = retry
        .heat
        .walls
        .into_iter()
        .flatten()
        .filter_map(|wall| wall.solid)
        .map(|solid| solid.source)
        .sum();
    near(
        total_energy(&actual) - energy_before,
        heat.iter().sum::<f64>()
            + solid_source
            + retry
                .heat
                .reservoirs
                .iter()
                .map(|reservoir| reservoir.heat_in)
                .sum::<f64>(),
        2e-11,
    );
    let before_clear = accepted(&actual);
    actual.clear_reservoir_temperature_profile(0).unwrap();
    assert!(actual.thermal().reservoir_temperature_profile(0).is_none());
    assert_eq!(
        actual.thermal().boundaries().faces[0],
        ThermalOpenFace3::Reservoir { temperature: 0.5 }
    );
    assert_eq!(accepted(&actual), before_clear);

    // Profiles may also be installed before attaching buoyancy. Its constructor
    // must admit every authored value before it creates a new force frame.
    let BoussinesqOpenGrid3 {
        flow,
        mut thermal,
        law,
        background,
        limits,
        ..
    } = channel(BETA);
    thermal
        .set_reservoir_temperature_profile(0, &refused)
        .unwrap();
    let result = BoussinesqOpenGrid3::try_new_with_poll(
        flow.inner,
        thermal,
        law,
        background,
        limits,
        &mut || Ok(()),
    );
    assert!(matches!(
        result,
        Err(BoussinesqOpenError3::BoundaryLaw { face: 0, .. })
    ));
}
