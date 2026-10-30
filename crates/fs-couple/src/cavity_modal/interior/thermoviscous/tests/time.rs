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
    let prepared = molecular
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let source = molecular.relaxation().unwrap();
    assert_eq!(
        prepared.frozen_cp_j_kg_k(),
        source.frozen_cv_j_kg_k() + gas.specific_gas_constant
    );
    assert_eq!(
        prepared.caloric_reservoirs().len(),
        source.processes().len()
    );
    for (state, original) in prepared.caloric_reservoirs().iter().zip(source.processes()) {
        assert_eq!(
            state.heat_capacity_j_kg_k.to_bits(),
            original.heat_capacity_j_kg_k.to_bits()
        );
        assert_eq!(
            state.relaxation_time_s.to_bits(),
            original.relaxation_time_s.to_bits()
        );
    }
    // The same material law, including the explicit phasor convention change.
    let frequency = 400.0;
    let mut cp = C64::from_re(prepared.frozen_cp_j_kg_k());
    for state in prepared.caloric_reservoirs() {
        cp = cp
            + C64::new(
                1.0,
                std::f64::consts::TAU * frequency * state.relaxation_time_s,
            )
            .recip()
            .scale(state.heat_capacity_j_kg_k);
    }
    close(
        cp,
        molecular.heat_capacity(frequency, &mut || false).unwrap(),
        cp.abs(),
        4e-15,
    );
    let mut initial = prepared.zero_state();
    // Mesh vertex zero is an isothermal wall coefficient. Its molecular state
    // must be admitted independently of the constrained translational field.
    assert_eq!(initial.temperature_k[0], 0.0);
    initial.molecular_temperature_k[0][0] = 1.0;
    let step = prepared
        .step(&initial, &[0.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.molecular_relaxation_heat_j > 0.0);
    assert!(step.work.final_energy_j < step.work.initial_energy_j);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
}
