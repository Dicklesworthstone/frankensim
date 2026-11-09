use super::*;

fn check_fields_unchanged<const Q: usize>(before: &Accepted<Q>, after: &Accepted<Q>) {
    assert_eq!(after.flow, before.flow);
    assert_eq!(after.heat, before.heat);
    assert_eq!(after.force, before.force);
    assert_eq!(after.background, before.background);
    assert_eq!(after.steps, before.steps);
    assert_eq!(after.heat_receipt, before.heat_receipt);
    assert_eq!(after.receipt, before.receipt);
}

#[test]
fn transient_wall_power_and_inlet_temperature_drive_original_buoyant_convection() {
    let mut coupled = channel(BETA);
    let heat = vec![0.0004; 72];
    advance(&mut coupled, &heat);
    for interval in 0..6 {
        let before = accepted(&coupled);
        let energy_before = total_energy(&coupled);
        let lower: Vec<_> = (0..18)
            .map(|link| (1 + link % 3) as f64 * if interval % 2 == 0 { 0.001 } else { -0.0004 })
            .collect();
        let upper = vec![-0.0001; 18];
        coupled.set_finite_wall_heat(2, &lower).unwrap();
        coupled.set_finite_wall_heat(3, &upper).unwrap();
        coupled
            .set_reservoir_temperature(0, 0.5 + 0.01 * interval as f64)
            .unwrap();
        let after = accepted(&coupled);
        check_fields_unchanged(&before, &after);
        assert_eq!(total_energy(&coupled).to_bits(), energy_before.to_bits());
        for face in [2, 3] {
            for (old, new) in before.walls[face].iter().zip(&after.walls[face]) {
                let (ThermalWall::FiniteSolid(old), ThermalWall::FiniteSolid(new)) = (old, new)
                else {
                    panic!("finite wall")
                };
                assert_eq!(old.temperature.to_bits(), new.temperature.to_bits());
                assert_eq!(old.heat_capacity.to_bits(), new.heat_capacity.to_bits());
                assert_eq!(old.resistance.to_bits(), new.resistance.to_bits());
            }
        }
        let report = advance(&mut coupled, &heat);
        let solid_source: f64 = [2, 3]
            .into_iter()
            .map(|face| report.heat.walls[face].unwrap().solid.unwrap().source)
            .sum();
        near(solid_source, lower.iter().chain(&upper).sum(), 2e-14);
        let reservoir_heat: f64 = report.heat.reservoirs.iter().map(|r| r.heat_in).sum();
        near(
            total_energy(&coupled) - energy_before,
            heat.iter().sum::<f64>() + solid_source + reservoir_heat,
            2e-11,
        );
        near(report.heat.combined_defect, 0.0, 2e-12);
        assert!(coupled.thermal().wall_conduction(1).is_some());
        for cell in 0..72 {
            let temperature = coupled.thermal().temperature(cell).unwrap();
            for axis in 0..3 {
                near(
                    coupled.flow().forces()[cell][axis],
                    BACKGROUND[axis] - REFERENCE_RHO * BETA * (temperature - 0.5) * GRAVITY[axis],
                    2e-15,
                );
            }
        }
    }
}

#[test]
fn wall_input_admission_trace_refusal_and_final_cancellation_preserve_the_joint_epoch() {
    let make_stateless = || {
        let mut faces = channel_faces();
        faces[2] = ThermalOpenFace3::Wall(ThermalWall::Temperature(0.8));
        faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
        make::<D3q19, 19>(
            CollisionModel3::Bgk { tau: 0.9 },
            [6, 4, 3],
            native_channel(),
            faces,
            0.5,
            BETA,
            0.02,
        )
    };
    let mut actual = make_stateless();
    let mut reference = make_stateless();
    let heat = vec![0.0002; 72];
    advance(&mut actual, &heat);
    advance(&mut reference, &heat);
    let before = accepted(&actual);
    let boundaries = actual.thermal().boundaries();
    for refused in [
        actual.set_wall_boundary(2, ThermalWall::Temperature(1.6)),
        actual.set_internal_wall_boundary(ThermalWall::Temperature(1.6)),
    ] {
        assert!(matches!(
            refused,
            Err(BoussinesqOpenError3::BoundaryLaw {
                source: BoussinesqLawError3::DensityDeparture { .. },
                ..
            })
        ));
        assert_eq!(accepted(&actual), before);
        assert_eq!(actual.thermal().boundaries(), boundaries);
    }

    // Flux data is finite and legal; its actual halfway trace must still fit
    // the Boussinesq density envelope when the original heat step resolves it.
    actual
        .set_wall_boundary(2, ThermalWall::HeatFlux(0.3))
        .unwrap();
    let with_flux = accepted(&actual);
    let flux_boundaries = actual.thermal().boundaries();
    let refused = actual.step_with_poll(None, &heat, None, [0.0; 3], &mut || Ok(()));
    assert!(
        matches!(
            refused,
            Err(BoussinesqOpenError3::BoundaryLaw {
                face: 2,
                source: BoussinesqLawError3::DensityDeparture { .. },
            })
        ),
        "{refused:?}"
    );
    assert_eq!(accepted(&actual), with_flux);
    assert_eq!(actual.thermal().boundaries(), flux_boundaries);

    for coupled in [&mut actual, &mut reference] {
        coupled
            .set_wall_boundary(2, ThermalWall::Temperature(0.9))
            .unwrap();
        coupled.set_reservoir_temperature(0, 0.6).unwrap();
    }
    let before = accepted(&actual);
    let boundaries = actual.thermal().boundaries();
    assert_eq!(before, accepted(&reference));
    let mut polls = 0;
    let expected = reference
        .step_with_poll(None, &heat, None, [0.0; 3], &mut || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let refused = actual.step_with_poll(None, &heat, None, [0.0; 3], &mut || {
        calls += 1;
        if calls == polls {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_eq!(accepted(&actual), before);
    assert_eq!(actual.thermal().boundaries(), boundaries);
    let retried = advance(&mut actual, &heat);
    assert_eq!(retried, expected);
    assert_eq!(accepted(&actual), accepted(&reference));
}
