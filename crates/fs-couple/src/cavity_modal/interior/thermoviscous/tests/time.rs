use super::*;

#[test]
fn original_gas_source_reaches_causal_field_work_without_caloric_shortcuts() {
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let original =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let prepared = original
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        prepared.inner().medium(),
        source_medium(original.gas(), 0.0).unwrap()
    );
    let step = prepared
        .step(&prepared.zero_state(), &[1.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.pressure_work_j > 0.0 && step.work.final_energy_j > 0.0);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.pressure_work_j);
    assert!(
        original
            .prepare_time(
                1.0 / BAND[1],
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let molecular = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    assert!(
        molecular
            .prepare_time(
                1.0 / 48_000.0,
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
}
